# coding: utf-8

"""
Custom tasks to derive JER SF.
"""

import law

from columnflow.tasks.framework.base import Requirements
from columnflow.tasks.histograms import MergeHistograms
from columnflow.config_util import get_datasets_from_process
from columnflow.util import maybe_import, DotDict

from dijet.tasks.base import HistogramsBaseTask

ak = maybe_import("awkward")
hist = maybe_import("hist")
np = maybe_import("numpy")
plt = maybe_import("matplotlib.pyplot")
mplhep = maybe_import("mplhep")
it = maybe_import("itertools")


class AlphaExtrapolation(HistogramsBaseTask):
    """
    Task to perform alpha extrapolation.
    Read in and plot asymmetry histograms.
    Cut of non-gaussian tails and extrapolate sigma_A( alpha->0 ).
    """

    # Add nested sibling directories to output path
    output_collection_cls = law.NestedSiblingFileCollection

    # upstream requirements
    reqs = Requirements(
        MergeHistograms=MergeHistograms,
    )

    def requires(self):
        return {
            d: self.reqs.MergeHistograms.req(
                self,
                dataset=d,
                branch=-1,
                _exclude={"branches"},
                _prefer_cli={"variables"},
            )
            for d in self.datasets
        }

    def create_branch_map(self):
        return [
            DotDict({"process": process})
            for process in sorted(self.processes)
        ]

    def workflow_requires(self):
        reqs = super().workflow_requires()
        reqs["merged_hists"] = self.requires_from_branch()
        return reqs

    def load_histogram(self, dataset, variable):
        histogram = self.input()[dataset]["collection"][0]["hists"].targets[variable].load(formatter="pickle")
        return histogram

    def get_datasets(self):
        # Get a samples from process
        process_sets = get_datasets_from_process(self.config_inst, self.branch_data.process, only_first=False)
        process_names = [item.name for item in process_sets]
        # Get all samples from input belonging to process
        samples = set(self.datasets).intersection(process_names)
        if len(samples) == 0:
            samples = process_names
        return (
            samples,
            self.config_inst.get_dataset(process_sets[0]).is_mc,
        )

    def output(self) -> dict[law.FileSystemTarget]:
        # TODO: Unstable for changes like data_jetmet_X
        #       Make independent like in config datasetname groups
        sample = self.extract_sample()
        target = self.target(f"{sample}", dir=True)
        # declare the main target
        outp = {
            "alphas": target.child("alphas.pickle", type="f"),
            "asym": target.child("asym.pickle", type="f", optional=True),
        }
        return outp

    def get_norm_asymmetries(self, histogram, method):

        # NOTE: Alternative loop over alpha like here:
        #       histogram[hist.loc(method), slice(hist.loc(0), hist.loc(0.3), sum), :, :, :].values()
        # TODO: Loose hist structure here. Methd to keep structure here ?
        #       histogram[hist.loc(method), slice(hist.loc(0), hist.loc(0.3)), :, :, :]
        #       For integral not taking .values()
        values = histogram[hist.loc(method), slice(hist.loc(0), hist.loc(0.3)), :, :, :].values()
        # axis = 0 alpha
        inclusiv = np.apply_along_axis(np.cumsum, 0, values)
        # axis = 3 asymmetry
        integral = inclusiv.sum(axis=3, keepdims=True)

        # Store in dictonary to use in alpha extrapolation
        return {"content": inclusiv, "integral": integral}

    def process_asymmetry(self, hists, asyms):

        inclusive_norm = hists["content"] / hists["integral"]
        means = np.nansum(asyms * inclusive_norm, axis=3, keepdims=True)

        # TODO: np.nanaverage ?
        #       -> only numpy.nanmean but no weights
        stds = np.sqrt(np.average(((asyms - means)**2), weights=inclusive_norm, axis=3))
        stds_err = stds / np.sqrt(np.squeeze(hists["integral"]))

        return {"widths": stds, "errors": stds_err}

    def method_index(self, method):
        indices = {"sm": 1, "fe": 2}
        # TODO: check that method is either sm or fe
        return indices[method]

    def run(self):
        # TODO: Gen level for MC
        #       Correlated fit (in jupyter)

        h_all = []
        datasets, isMC = self.get_datasets()

        for dataset in datasets:
            for variable in self.variables:
                h_in = self.reduce_histogram(
                    self.load_histogram(dataset, variable),
                    self.processes,
                    self.shift,
                )
                # Store all hists in a list to sum over after reading
                h_all.append(h_in)
        h_all = sum(h_all)
        amax = np.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3])

        # TODO: Need own task to store asymmetry before this one
        #       New structure of base histogram task necessary
        asym_sm = self.get_norm_asymmetries(h_all, self.method_index("sm"))
        asym_fe = self.get_norm_asymmetries(h_all, self.method_index("fe"))
        results_asym = {
            "sm": asym_sm["content"] / asym_sm["integral"],
            "fe": asym_fe["content"] / asym_fe["integral"],
            "bins": {
                "pt": h_all.axes["dijets_pt_avg"].edges,
                "eta": h_all.axes["probejet_abseta"].edges,
                "alpha": amax,
            },
        }
        self.output()["asym"].dump(results_asym, formatter="pickle")

        # Get binning
        centers_ptavg = h_all.axes["dijets_pt_avg"].centers
        centers_asym = h_all.axes["dijets_asymmetry"].centers
        centers_eta = h_all.axes["probejet_abseta"].centers

        print(f"Number \u03B7 bins: {len(centers_asym)}")
        print(f"Number pT bins: {len(centers_eta)}")
        print(f"Number A bins: {len(centers_ptavg)}")

        sm = self.process_asymmetry(asym_sm, centers_asym)
        fe = self.process_asymmetry(asym_fe, centers_asym)

        # TODO: Use correlated fit
        def fit_linear(subarray):
            fitting = ~np.isnan(subarray)
            if len(subarray[fitting]) < 2:
                coefficients = [0, 0]
            else:
                coefficients = np.polyfit(amax[fitting], subarray[fitting], 1)
            return coefficients
        # TODO: Fit Messungen abspeichern (chi2, ndf, etc.) for diagnostic
        fit_sm = np.apply_along_axis(fit_linear, axis=0, arr=sm["widths"])
        fit_fe = np.apply_along_axis(fit_linear, axis=0, arr=fe["widths"])

        results_alphas = {}
        results_alphas["sm"] = {
            "alphas": sm["widths"],
            "fits": fit_sm,
        }
        results_alphas["fe"] = {
            "alphas": fe["widths"],
            "fits": fit_fe,
        }
        results_asym["bins"] = {
            "pt": h_all.axes["dijets_pt_avg"].edges,
            "eta": h_all.axes["probejet_abseta"].edges,
            "alpha": amax,
        }
        self.output()["alphas"].dump(results_alphas, formatter="pickle")
