# coding: utf-8

import itertools
import law

# from dijet.tasks.base import HistogramsBaseTask
from columnflow.util import maybe_import, DotDict
from columnflow.tasks.framework.base import Requirements
from columnflow.tasks.framework.remote import RemoteWorkflow

from dijet.tasks.asymmetry import Asymmetry
from dijet.constants import eta
from dijet.plotting.base import PlottingBaseTask
from dijet.plotting.util import eta_bin, pt_bin, alpha_bin, add_text, dot_to_p

hist = maybe_import("hist")
np = maybe_import("numpy")
plt = maybe_import("matplotlib.pyplot")
mplhep = maybe_import("mplhep")


class PlotAsymmetries(
    PlottingBaseTask,
    law.LocalWorkflow,
    RemoteWorkflow,
):
    """
    Task to plot asymmetry distributions. Shows the distribution for
    `--asymmetry-variable` for all given `--samples` and `--levels`.
    One plot is produced for each eta, pt and alpha bin for each method (fe,sm).
    The methods to take into account are given as `--categories`.
    """

    # how to create the branch map
    branching_type = "merged"

    # upstream requirements
    reqs = Requirements(
        RemoteWorkflow.reqs,
        Asymmetry=Asymmetry,
    )

    # TODO: use config
    colors = {
        "data": "black",
        "qcdht": "indianred",
    }

    #
    # methods required by law
    #

    def create_branch_map(self):
        """
        Workflow has one branch for each eta bin (eta).
        TODO: Hard coded for now
              Into Base Task
        """
        input_branches = super().create_branch_map()

        branches = []
        for ib in input_branches:
            branches.extend([
                DotDict.wrap(dict(ib, **{
                    "eta": (eta_lo, eta_hi),
                }))
                for eta_lo, eta_hi in zip(eta[:-1], eta[1:])
            ])

        return branches

    def output(self):
        """
        Organize output as a (nested) dictionary. Output files will be in a single
        directory, which is determined by `store_parts`.
        """
        eta_lo, eta_hi = self.branch_data.eta
        n_eta_lo = dot_to_p(eta_lo)
        n_eta_hi = dot_to_p(eta_hi)
        return {
            "dummy": self.target(f"eta_{n_eta_lo}_{n_eta_hi}/DUMMY"),
            "plots": self.target(f"eta_{n_eta_lo}_{n_eta_hi}", dir=True),
        }

    def requires(self):
        return self.reqs.Asymmetry.req_different_branching(self, branch=-1)

    def workflow_requires(self):
        reqs = super().workflow_requires()
        reqs["key"] = self.requires_from_branch()
        return reqs

    #
    # helper methods for handling task inputs/outputs
    #

    def load_input(self, key: str, sample: str, level: str):
        coll_keys = [
            coll_key
            for coll_key, coll in self.input()["collection"].targets.items()
            if sample in coll[key]
        ]
        if len(coll_keys) != 1:
            raise RuntimeError(
                f"found {len(coll_keys)} input collections corresponding to "
                f"sample '{sample}', expected 1",
            )
        return self.input()["collection"][coll_keys[0]][key][sample][level].load(formatter="pickle")

    #
    # task implementation
    #

    def plot_asymmetry(self, content, error, asym):
        fig, ax = plt.subplots()
        plt.bar(
            asym.flatten(),
            content["mc"].flatten(),
            yerr=error["mc"].flatten(),
            align="center",
            width=np.diff(asym)[0],
            alpha=0.6,
            color=self.colors["mc"],
            edgecolor="none",
            label="MC",
        )
        plt.errorbar(
            asym.flatten(),
            content["da"].flatten(),
            yerr=error["da"].flatten(),
            fmt="o",
            marker="o",
            fillstyle="full",
            color=self.colors["da"],
            label="Data",
        )

        ax.set_xlabel("Asymmetry")
        ax.set_ylabel(r"$\Delta$N/N")
        ax.set_yscale("log")
        return fig, ax

    @staticmethod
    def _plot_shim(x, y, xerr=None, yerr=None, method=None, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig, ax = plt.gcf(), plt.gca()

        method = method or "errorbar"

        method_func = getattr(ax, method, None)
        if method_func is None:
            raise ValueError(f"invalid plot method '{method}'")

        if method == "bar":
            kwargs.update(
                align="center",
                width=2 * xerr,
            )
        else:
            kwargs["xerr"] = xerr

        method_func(
            x.flatten(),
            y.flatten(),
            yerr=yerr,
            **kwargs,
        )

        return fig, ax

    def run(self):
        # load inputs (asymmetries and quantiles)
        raw_inputs = {
            key: {
                sample: {
                    level: self.load_input(key, sample=sample, level=level)
                    for level in self.levels
                    if level == "reco" or sample != "data"
                }
                for sample in self.samples
            }
            for key in ("asym", "quantile")
        }

        # dict storing either variables or their gen-level equivalents
        # for convenient access
        vars_ = {
            level: self._make_var_lookup(level=level)
            for level in self.levels
        }

        # eta binning information from branch
        eta_lo, eta_hi = self.branch_data.eta
        eta_midp = 0.5 * (eta_lo + eta_hi)

        # prepare inputs (apply slicing, clean nans)
        def _prepare_input(histogram):
            # # slice histogram to extact bin
            # histogram = histogram[slicer]
            # map `nan` values to zero
            v = histogram.view()
            v.value = np.nan_to_num(v.value, nan=0.0)
            v.variance = np.nan_to_num(v.variance, nan=0.0)
            # return prepared histogram
            return histogram
        inputs = law.util.map_struct(_prepare_input, raw_inputs)

        # binning information from first histogram object
        # (assume identical binning for all)
        def _iter_flat(d: dict):
            if not isinstance(d, dict):
                yield d
                return
            for k, v in d.items():
                yield from _iter_flat(v)
        ref_object = next(_iter_flat(inputs["asym"]))

        alpha_edges = ref_object.axes[vars_["reco"]["alpha"]].edges
        binning_variable_edges = {
            bv: ref_object.axes[bv_resolved].edges
            for bv, bv_resolved in vars_["reco"]["binning"].items()
            if bv != "probejet_abseta"
        }

        plt.style.use(mplhep.style.CMS)

        def iter_bins(edges, **add_kwargs):
            for i, (e_dn, e_up) in enumerate(zip(edges[:-1], edges[1:])):
                yield {"up": e_up, "dn": e_dn, **add_kwargs}

        # loop through bins
        for m, (ia, alpha_up), *bv_bins in itertools.product(
            self.LOOKUP_CATEGORY_ID,
            enumerate(alpha_edges[1:]),
            *[iter_bins(bv_edges, var_name=bv) for bv, bv_edges in binning_variable_edges.items()],
        ):
            # initialize figure and axes
            fig, ax = plt.subplots()
            mplhep.cms.label(
                lumi=41.48,  # TODO: from self.config_inst.x.luminosity?
                com=13,
                ax=ax,
                llabel="Private Work",
                data=True,
            )

            # selector to get current bin
            bin_selectors = {}
            for level in self.levels:
                bin_selectors[level] = {
                    "category": hist.loc(self.LOOKUP_CATEGORY_ID[m]),
                    vars_[level]["alpha"]: hist.loc(alpha_up - 0.001),
                    vars_[level]["binning"]["probejet_abseta"]: hist.loc(eta_midp),
                }
                for bv_bin in bv_bins:
                    bin_selectors[level][vars_[level]["binning"][bv_bin["var_name"]]] = (
                        hist.loc(0.5 * (bv_bin["up"] + bv_bin["dn"]))
                    )

            # loop through samples/levels
            for sample, level in itertools.product(self.samples, self.levels):
                # only use reco level for data
                if sample == "data" and level != "reco":
                    continue

                # get input histogram
                h_in = inputs["asym"][sample][level]
                h_sliced = h_in[bin_selectors[level]]

                # plot
                plot_kwargs = dict(
                    self.config_inst.x("samples", {})
                    .get(sample, {}).get("plot_kwargs", {}),
                )
                if level == "gen":
                    plot_kwargs.update({
                        "color": "forestgreen",
                        "label": "MC (gen)",
                        # "method": "steps",  # TODO
                    })

                self._plot_shim(
                    h_sliced.axes[vars_[level]["asymmetry"]].centers,
                    h_sliced.values(),
                    xerr=h_sliced.axes[vars_[level]["asymmetry"]].widths / 2,
                    yerr=np.sqrt(h_sliced.variances()),
                    ax=ax,
                    **plot_kwargs,
                )

                # plot quantiles
                hs_in_quantiles = inputs["quantile"][sample][level]
                q_lo = hs_in_quantiles["low"][bin_selectors[level]].value
                q_up = hs_in_quantiles["up"][bin_selectors[level]].value
                ax.axvline(q_lo, 0, 0.67, color=plot_kwargs.get("color"), linestyle="--")
                ax.axvline(q_up, 0, 0.67, color=plot_kwargs.get("color"), linestyle="--")

            # annotations
            add_text(ax, 0.05, 0.95, alpha_bin(alpha_up))
            add_text(ax, 0.05, 0.95, eta_bin(eta_lo, eta_hi), offset=0.05)
            for i, bv_bin in enumerate(bv_bins):
                add_text(ax, 0.05, 0.95, pt_bin(bv_bin["dn"], bv_bin["up"]), offset=0.05 * (2 + i))

            # figure adjustments
            ax.set_xlim(-0.5, 0.5)
            ax.set_ylim(5e-5, 10)
            ax.set_yscale("log")
            ax.legend(loc="upper right")
            ax.set_xlabel(
                self.config_inst.get_variable(vars_["reco"]["asymmetry"]).x_title,
            )
            ax.set_ylabel(r"$\Delta$N/N")

            # compute plot filename
            fname_parts = [
                # base name
                "asym",
                # method
                m,
                # alpha
                f"alpha_lt_{dot_to_p(alpha_up)}",
                # abseta bin
                f"abseta_{dot_to_p(eta_lo)}_{dot_to_p(eta_hi)}",
            ]
            # other bins
            for bv_bin in bv_bins:
                fname_parts.append("_".join([
                    bv_bin["var_name"],
                    dot_to_p(bv_bin["dn"]),
                    dot_to_p(bv_bin["up"]),
                ]))

            # save plot to file
            self.save_plot("__".join(fname_parts), fig)
            plt.close(fig)
