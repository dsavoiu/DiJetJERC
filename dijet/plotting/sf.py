# coding: utf-8
from __future__ import annotations

import law

# from dijet.tasks.base import HistogramsBaseTask
from columnflow.util import maybe_import
from columnflow.tasks.framework.base import Requirements

from dijet.tasks.sf import SF
from dijet.plotting.base import PlottingBaseTask
from dijet.plotting.util import eta_bin, add_text, dot_to_p

hist = maybe_import("hist")
np = maybe_import("numpy")
plt = maybe_import("matplotlib.pyplot")
mplhep = maybe_import("mplhep")


class PlotSFs(
    PlottingBaseTask,
):
    """
    Task to plot all SFs.
    One plot for each eta bin for each method (fe,sm).
    """

    # how to create the branch map
    branching_type = SF.branching_type

    # upstream requirements
    reqs = Requirements(
        PlottingBaseTask.reqs,
        SF=SF,
    )

    ### @law.workflow_property(setter=True, cache=True, empty_value=None)
    ### def binning_info(self):
    ###     """
    ###     Open the SF task outputs to determine binning and set up plot. 
    ###     """
    ###     # check if the merging stats are present
    ###     sf_task = self.reqs.SF.req_different_branching(self)
    ###     sf_output = sf_task.output().collection[0]["sfs"]
    ###     if not sf_output.exists():
    ###         return None

    ###     sfs = sf_output.load(formatter="pickle")
    ###     return {
    ###         ax.name: ax.edges
    ###         for ax in sfs.axes
    ###     }

    ### @law.dynamic_workflow_condition
    ### def workflow_condition(self):
    ###     # the workflow can be constructed as soon as the binning information is known
    ###     return self.binning_info is not None

    ### @workflow_condition.create_branch_map
    ### def create_branch_map(self):
    ###     __import__("IPython").embed()
    ###     return []

    #
    # methods required by law
    #

    def output(self):
        """
        Organize output as a (nested) dictionary. Output files will be in a single
        directory, which is determined by `store_parts`.
        """
        return {
            "dummy": self.target("DUMMY"),
            "plots": self.target("plots", dir=True),
        }

    def requires(self):
        return self.reqs.SF.req_different_branching(self, branch=-1)

    def workflow_requires(self):
        reqs = super().workflow_requires()
        reqs["key"] = self.requires_from_branch()
        return reqs

    #
    # helper methods for handling task inputs/outputs
    #

    def load_input(self, key: str):
        return self.input()[key].load(formatter="pickle")


    #
    # task implementation
    #

    def plot_sfs(self, sfs, pt):
        fig, ax = plt.subplots()

        plt.errorbar(pt, sfs["nom"], yerr=sfs["err"], fmt="o", color="black")

        ax.set_xlabel(r"$p_{T}^{ave}$")
        ax.set_ylabel("SF")
        return fig, ax

    def run(self):
        sfs = self.load_input("sfs")["sfs"]

        sfs.view().value = np.nan_to_num(sfs.view().value, nan=0.0)
        sfs.view().variance = np.nan_to_num(sfs.view().variance, nan=0.0)

        # TODO: don't hardcode axis names
        eta_edges = sfs.axes["probejet_abseta"].edges
        pt_centers = sfs.axes["dijets_pt_avg"].centers

        # Set plotting style
        plt.style.use(mplhep.style.CMS)

        pos_x = 0.05
        pos_y = 0.95
        # TODO: setup loops in branch map?
        for m in self.LOOKUP_CATEGORY_ID:
            for ie, (eta_lo, eta_hi) in enumerate(zip(eta_edges[:-1], eta_edges[1:])):

                input_ = {
                    "sfs": {
                        "nom": sfs[hist.loc(self.LOOKUP_CATEGORY_ID[m]), ie, :].values(),
                        "err": sfs[hist.loc(self.LOOKUP_CATEGORY_ID[m]), ie, :].variances(),
                    },
                    "pt": pt_centers,
                }

                fig, ax = self.plot_sfs(**input_)
                mplhep.cms.label(
                    lumi=41.48,  # TODO: from self.config_inst.x.luminosity?
                    com=13,
                    ax=ax,
                    llabel="Private Work",
                    data=True,
                )

                text_eta_bin = eta_bin(eta_lo, eta_hi)
                add_text(ax, pos_x, pos_y, text_eta_bin)
                print(f"Start with eta {text_eta_bin} for {m} method")

                plt.xlim(49, 2100)
                # plt.ylim(0.8, 20)
                ax.set_xscale("log")
                plt.legend(loc="upper right")

                store_bin_eta = f"eta_{dot_to_p(eta_lo)}_{dot_to_p(eta_hi)}"
                self.save_plot(f"sfs_{m}_{store_bin_eta}", fig)
                plt.close(fig)
