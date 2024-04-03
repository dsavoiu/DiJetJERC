# coding: utf-8

import law

# from dijet.tasks.base import HistogramsBaseTask
from columnflow.util import maybe_import, DotDict
from columnflow.tasks.framework.base import Requirements

from dijet.tasks.alpha import AlphaExtrapolation
from dijet.constants import eta
from dijet.plotting.base import PlottingBaseTask
from dijet.plotting.util import eta_bin, pt_bin, add_text, dot_to_p

hist = maybe_import("hist")
np = maybe_import("numpy")
plt = maybe_import("matplotlib.pyplot")
mplhep = maybe_import("mplhep")


class PlotWidths(PlottingBaseTask):
    """
    Task to plot all alphas.
    One plot for each eta and pt bin for each method (fe,sm).
    """

    def load_widths(self):
        return (
            self.input().collection[0]["widths"].load(formatter="pickle"),
            self.input().collection[1]["widths"].load(formatter="pickle"),
        )

    def load_extrapolation(self):
        return (
            self.input().collection[0]["extrapolation"].load(formatter="pickle"),
            self.input().collection[1]["extrapolation"].load(formatter="pickle"),
        )

    def run(self):
        widths_da, widths_mc = self.load_widths()
        extrapol_da, extrapol_mc = self.load_extrapolation()
