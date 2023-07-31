# coding: utf-8

"""
Selection methods for HHtobbWW.
"""

from operator import and_
from functools import reduce
from collections import defaultdict
from typing import Tuple

from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column
from columnflow.production.util import attach_coffea_behavior

from columnflow.selection import Selector, SelectionResult, selector
from columnflow.selection.stats import increment_event_stats
from columnflow.production.cms.mc_weight import mc_weight
from columnflow.production.categories import category_ids
from columnflow.production.processes import process_ids

from dijet.production.weights import event_weights_to_normalize, large_weights_killer
from dijet.production.dijet_balance import dijet_balance
from dijet.selection.jet_selection import jet_selection
from dijet.selection.cutflow_features import cutflow_features
from dijet.selection.stats import dijet_increment_stats

np = maybe_import("numpy")
ak = maybe_import("awkward")

def masked_sorted_indices(mask: ak.Array, sort_var: ak.Array, ascending: bool = False) -> ak.Array:
    """
    Helper function to obtain the correct indices of an object mask
    """
    indices = ak.argsort(sort_var, axis=-1, ascending=ascending)
    return indices[mask[indices]]

@selector(
    uses={
        category_ids, process_ids, attach_coffea_behavior,
        mc_weight, large_weights_killer,  # not opened per default but always required in Cutflow tasks
        jet_selection, dijet_balance, cutflow_features, dijet_increment_stats,
    },
    produces={
        category_ids, process_ids, attach_coffea_behavior,
        mc_weight, large_weights_killer,
        jet_selection, dijet_balance, cutflow_features, dijet_increment_stats,
    },
    exposed=True,
    check_used_columns=False,
    check_produced_columns=False,
)
def default(
    self: Selector,
    events: ak.Array,
    stats: defaultdict,
    **kwargs,
) -> Tuple[ak.Array, SelectionResult]:

    if self.dataset_inst.is_mc:
        events = self[mc_weight](events, **kwargs)
        events = self[large_weights_killer](events, **kwargs)

    # ensure coffea behavior
    events = self[attach_coffea_behavior](events, **kwargs)

    # prepare the selection results that are updated at every step
    results = SelectionResult()

    # build categories
    events = self[category_ids](events, results=results, **kwargs)

    # create process ids
    events = self[process_ids](events, **kwargs)

    # # TODO Implement selection
    # # lepton selection
    # events, lepton_results = self[lepton_selection](events, stats, **kwargs)
    # results += lepton_results

    # jet selection
    events, results_jet = self[jet_selection](events, **kwargs)
    results += results_jet

    # dijet balance for cutflow variables
    # TODO: Remove later
    events = self[dijet_balance](events, **kwargs)

    # produce relevant columns
    events = self[cutflow_features](events, results.objects, **kwargs)

    # results.main.event contains full selection mask. Sum over all steps.
    # Make sure all nans are present, otherwise next tasks fail
    results.main["event"] = reduce(and_, results.steps.values())
    results.main["event"] = ak.fill_none(results.main["event"], False)
    
    self[dijet_increment_stats](events, results, stats, **kwargs)

    return events, results


# @default.init
# def default_init(self: Selector) -> None:
#     if self.config_inst.x("do_cutflow_features", False):
#         self.uses.add(cutflow_features)
#         self.produces.add(cutflow_features)

#     if not getattr(self, "dataset_inst", None) or self.dataset_inst.is_data:
#         return

#     self.uses.add(event_weights_to_normalize)
#     self.produces.add(event_weights_to_normalize)
