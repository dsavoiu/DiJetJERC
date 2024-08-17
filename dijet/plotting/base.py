# coding: utf-8
"""
Custom base tasks for plotting workflow steps from JER SF measurement.
"""
from __future__ import annotations

import law

from columnflow.tasks.framework.base import Requirements
from columnflow.tasks.framework.remote import RemoteWorkflow

from dijet.tasks.base import HistogramsBaseTask


class PlottingBaseTask(
    HistogramsBaseTask,
    law.LocalWorkflow,
    RemoteWorkflow,
):
    """
    Base task to plot histogram from each step of the JER SF workflow.
    An example implementation of how to handle the inputs in a run method can be
    found in columnflow/tasks/histograms.py
    """

    default_plot_extensions = ("pdf", "png")

    # upstream requirements
    reqs = Requirements(
        RemoteWorkflow.reqs,
    )

    #
    # methods required by law
    #

    def create_branch_map(self):
        return super().create_branch_map()

    #
    # helper methods for handling task inputs/outputs
    #

    def save_plot(self, basename: str, fig: object, extensions: tuple[str] | list[str] | None = None):
        extensions = extensions or self.default_plot_extensions
        for ext in extensions:
            target = self.output()["plots"].child(f"{basename}.{ext}", type="f")
            target.dump(fig, formatter="mpl")
            print(f"saved plot: {target.path}")
