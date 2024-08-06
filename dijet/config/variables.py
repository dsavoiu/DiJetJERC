# coding: utf-8

"""
Definition of variables.
"""

import order as od

from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT, Route

np = maybe_import("numpy")
ak = maybe_import("awkward")


from dijet.constants import pt, eta, alpha


def add_feature_variables(config: od.Config) -> None:
    """
    Adds variables to a *config* that are produced as part of the `features` producer.
    """

    # Event properties
    config.add_variable(
        name="n_jet",
        binning=(12, -0.5, 11.5),
        x_title="Number of jets",
        discrete_x=True,
    )

    # jj features
    config.add_variable(
        name="deltaR_jj",
        binning=(40, 0, 5),
        x_title=r"$\Delta R(j_{1},j_{2})$",
    )


def add_variables(config: od.Config) -> None:
    """
    Adds all variables to a *config* that are present after `ReduceEvents`
    without calling any producer
    """

    # (the "event", "run" and "lumi" variables are required for some cutflow plotting task,
    # and also correspond to the minimal set of columns that coffea's nano scheme requires)
    config.add_variable(
        name="event",
        expression="event",
        binning=(1, 0.0, 1.0e9),
        x_title="Event number",
        discrete_x=True,
    )
    config.add_variable(
        name="run",
        expression="run",
        binning=(1, 100000.0, 500000.0),
        x_title="Run number",
        discrete_x=True,
    )
    config.add_variable(
        name="lumi",
        expression="luminosityBlock",
        binning=(1, 0.0, 5000.0),
        x_title="Luminosity block",
        discrete_x=True,
    )

    #
    # Weights
    #

    # TODO: implement tags in columnflow; meanwhile leave these variables commented out (as they only work for mc)
    config.add_variable(
        name="npvs",
        expression="PV.npvs",
        binning=(51, -.5, 50.5),
        x_title="Number of primary vertices",
        discrete_x=True,
    )

    #
    # Object properties
    #

    config.add_variable(
        name="jets_pt",
        expression="Jet.pt",
        binning=(100, 0, 1000),
        unit="GeV",
        x_title="$p_{T}$ of all jets",
    )

    # Jets (3 pt-leading jets)
    for i in range(3):
        config.add_variable(
            name=f"jet{i+1}_pt",
            expression=f"Jet.pt[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(100, 0., 1000.),
            unit="GeV",
            x_title=r"Jet %i $p_{T}$" % (i + 1),
        )
        config.add_variable(
            name=f"jet{i+1}_eta",
            expression=f"Jet.eta[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(100, -5.0, 5.0),
            x_title=r"Jet %i $\eta$" % (i + 1),
        )
        config.add_variable(
            name=f"jet{i+1}_phi",
            expression=f"Jet.phi[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(40, -3.2, 3.2),
            x_title=r"Jet %i $\phi$" % (i + 1),
        )
        config.add_variable(
            name=f"jet{i+1}_mass",
            expression=f"Jet.mass[:,{i}]",
            null_value=EMPTY_FLOAT,
            binning=(40, 0, 200),
            unit="GeV",
            x_title=r"Jet %i mass" % (i + 1),
        )

    #
    # dijet-related variables
    #

    config.add_variable(
        name="dijets_asymmetry",
        expression="dijets.asymmetry",
        binning=(160, -0.8, 0.8),
        x_title=r"A",
    )

    config.add_variable(
        name="dijets_pt_avg",
        expression="dijets.pt_avg",
        binning=pt,
        x_title=r"$p_{T}^{avg}$",
        unit="GeV",
    )

    # eta bin is always defined by probe jet
    # SM: Both jets in eta bin
    # FE: Probejet in eta bin
    config.add_variable(
        name="probejet_abseta",
        expression=lambda events: abs(events.probe_jet.eta),
        binning=eta,
        x_title=r"$|\eta|$",
        aux={
            "inputs": {"probe_jet.eta"},
        },
    )

    config.add_variable(
        name="dijets_alpha",
        expression="alpha",
        #binning=(100, 0, 1),  # ?
        binning=alpha,
        x_title=r"$\alpha$",
    )

    config.add_variable(
        name="dijets_alpha_fine",
        expression="alpha",
        binning=(100, 0, 1),
        x_title=r"$\alpha$",
    )

    config.add_variable(
        name="dijets_mpf",
        expression="dijets.mpf",
        binning=(100, -1, 1),
        x_title=r"MPF",
    )

    config.add_variable(
        name="dijets_mpfx",
        expression="dijets.mpfx",
        binning=(100, -1, 1),
        x_title=r"MPFx",
    )

    config.add_variable(
        name="dijets_response_probe",
        expression="dijets.response_probe",
        binning=(100, 0, 2),
        x_title=r"response probe jet",
    )
    config.add_variable(
        name="dijets_response_reference",
        expression="dijets.response_reference",
        binning=(100, 0, 2),
        x_title=r"response reference jet",
    )

    config.add_variable(
        name="dijets_mpf_gen",
        expression="dijets.mpf_gen",
        binning=(100, -0.025, 0.025),
        x_title=r"MPF gen",
    )

    config.add_variable(
        name="dijets_mpfx_gen",
        expression="dijets.mpfx_gen",
        binning=(100, -0.025, 0.025),
        x_title=r"MPFx gen",
    )

    config.add_variable(
        name="dijets_pt_avg_gen",
        expression="dijets.pt_avg_gen",
        binning=pt,
        x_title=r"$p_{T}^{avg, gen}$",
        unit="GeV",
    )

def add_uhh2_synch_variables(config: od.Config) -> None:
    """
    Add variables needed for synchronization with UHH2.
    """
    variables = {
        "jet1_pt": {
            "expression_cf": lambda x: Route("Jet.pt[:, 0]").apply(x, None),
            "expression_uhh2": "uhh2.jet1_pt",
            "inputs": {"Jet.*"},
            "scale": "log",
            "binning": np.logspace(np.log10(5), np.log10(500), 101),
            "diff_binning": np.linspace(-200, 200, 101),
            "label": r"Jet 1 $p_{T}$",
        },
        "jet2_pt": {
            "expression_cf": lambda x: Route("Jet.pt[:, 1]").apply(x, None),
            "expression_uhh2": "uhh2.jet2_pt",
            "inputs": {"Jet.*"},
            "scale": "log",
            "binning": np.logspace(np.log10(5), np.log10(500), 101),
            "diff_binning": np.linspace(-200, 200, 101),
            "label": r"Jet 2 $p_{T}$",
        },
        "jet3_pt": {
            "expression_cf": lambda x: Route("Jet.pt[:, 2]").apply(x, None),
            "expression_uhh2": "uhh2.jet3_pt",
            "inputs": {"Jet.*"},
            "scale": "log",
            "binning": np.logspace(np.log10(5), np.log10(500), 101),
            "diff_binning": np.linspace(-200, 200, 101),
            "label": r"Jet 3 $p_{T}$",
        },
        "asymmetry": {
            "expression_cf": "dijets.asymmetry",
            "expression_uhh2": "uhh2.asymmetry",
            "binning": np.linspace(-1, 1, 101),
            "diff_binning": np.linspace(-2, 2, 101),
            "label": r"Asymmetry",
        },
        "asymmetry_signflip": {
            "expression_cf": lambda x: -x.dijets.asymmetry,
            "expression_uhh2": "uhh2.asymmetry",
            "inputs": {"dijets.asymmetry"},
            "binning": np.linspace(-1, 1, 101),
            "diff_binning": np.linspace(-2, 2, 101),
            "label": r"Asymmetry (sign-flipped)",
        },
        "pt_avg": {
            "expression_cf": "dijets.pt_avg",
            "expression_uhh2": "uhh2.pt_ave",
            "scale": "log",
            "binning": np.logspace(np.log10(5), np.log10(500), 101),
            "diff_binning": np.linspace(-2, 2, 101),
            "label": r"Dijet $p_{T}^{avg}$",
        },
        "pt_avg_gen": {
            "expression_cf": "dijets.pt_avg_gen",
            "expression_uhh2": "uhh2.gen_pt_ave",
            "scale": "log",
            "binning": np.logspace(np.log10(5), np.log10(500), 101),
            "diff_binning": np.linspace(-2, 2, 101),
            "label": r"Dijet gen $p_{T}^{avg}$",
        },
        "alpha": {
            "expression_cf": "alpha",
            "expression_uhh2": "uhh2.alpha",
            "binning": np.linspace(-1, 1, 101),
            "diff_binning": np.linspace(-2, 2, 101),
            "label": r"$\alpha$"
        },
        "reference_jet_abseta": {
            "expression_cf": lambda x: abs(x.reference_jet.eta),
            "expression_uhh2": lambda x: abs(x.uhh2.barreljet_eta),
            "inputs": {"reference_jet.*", "uhh.*"},
            "binning": np.linspace(0, 5, 50),
            "diff_binning": np.linspace(-5, 5, 101),
            "label": r"Reference jet $\eta$",
        },
        "reference_jet_eta": {
            "expression_cf": "reference_jet.eta",
            "expression_uhh2": "uhh2.barreljet_eta",
            "binning": np.linspace(-5, 5, 101),
            "diff_binning": np.linspace(-2, 2, 101),
            "label": r"Reference jet $|\eta|$",
        },
        "reference_jet_phi": {
            "expression_cf": "reference_jet.phi",
            "expression_uhh2": "uhh2.barreljet_phi",
            "binning": np.linspace(-3.2, 3.2, 101),
            "diff_binning": np.linspace(-2, 2, 101),
            "label": r"Reference jet $\phi$",
        },
        "reference_jet_pt": {
            "expression_cf": "reference_jet.pt",
            "expression_uhh2": "uhh2.barreljet_pt",
            "scale": "log",
            "binning": np.logspace(np.log10(5), np.log10(500), 101),
            "diff_binning": np.linspace(-100, 100, 101),
            "label": r"Reference jet $p_{T}$",
        },
        "reference_gen_jet_eta": {
            "expression_cf": lambda x: ak.firsts(x["GenJet"][ak.singletons(x["reference_jet"]["genJetIdx"])].eta),
            "expression_uhh2": "uhh2.barrelgenjet_eta",
            "inputs": {"reference_jet.*", "GenJet.*"},
            "binning": np.linspace(-5, 5, 101),
            "diff_binning": np.linspace(-2, 2, 101),
            "label": r"Reference gen jet $\eta$",
        },
        "reference_gen_jet_phi": {
            "expression_cf": lambda x: ak.firsts(x["GenJet"][ak.singletons(x["reference_jet"]["genJetIdx"])].phi),
            "expression_uhh2": "uhh2.barrelgenjet_phi",
            "inputs": {"reference_jet.*", "GenJet.*"},
            "binning": np.linspace(-3.2, 3.2, 101),
            "diff_binning": np.linspace(-2, 2, 101),
            "label": r"Reference gen jet $\phi$",
        },
        "reference_gen_jet_pt": {
            "expression_cf": lambda x: ak.firsts(x["GenJet"][ak.singletons(x["reference_jet"]["genJetIdx"])].pt),
            "expression_uhh2": "uhh2.barrelgenjet_pt",
            "inputs": {"reference_jet.*", "GenJet.*"},
            "scale": "log",
            "binning": np.logspace(np.log10(5), np.log10(500), 101),
            "diff_binning": np.linspace(-100, 100, 101),
            "label": r"Reference gen jet $p_{T}$",
        },
        "probe_jet_abseta": {
            "expression_cf": lambda x: abs(x.probe_jet.eta),
            "expression_uhh2": lambda x: abs(x.uhh2.probejet_eta),
            "inputs": {"probe_jet.*", "uhh.*"},
            "binning": np.linspace(0, 5, 50),
            "diff_binning": np.linspace(-5, 5, 101),
            "label": r"Probe jet $|\eta|$"
        },
        "probe_jet_eta": {
            "expression_cf": "probe_jet.eta",
            "expression_uhh2": "uhh2.probejet_eta",
            "binning": np.linspace(-5, 5, 101),
            "diff_binning": np.linspace(-2, 2, 101),
            "label": r"Probe jet $\eta$"
        },
        "probe_jet_phi": {
            "expression_cf": "probe_jet.phi",
            "expression_uhh2": "uhh2.probejet_phi",
            "binning": np.linspace(-3.2, 3.2, 101),
            "diff_binning": np.linspace(-2, 2, 101),
            "label": r"Probe jet $\phi$"
        },
        "probe_jet_pt": {
            "expression_cf": "probe_jet.pt",
            "expression_uhh2": "uhh2.probejet_pt",
            "scale": "log",
            "binning": np.logspace(np.log10(5), np.log10(500), 101),
            "diff_binning": np.linspace(-100, 100, 101),
            "label": r"Probe jet $p_{T}$"
        },
    }

    variable_groups = {
        "synch": set(),
    }
    for var_basename, var_cfg in variables.items():
        # add UHH2 and CF variables
        for fw in ("uhh2", "cf"):
            var_name = f"s_{var_basename}_{fw}"
            var_binning = var_cfg["binning"]
            if not isinstance(var_binning, tuple):
                var_binning = list(var_binning)
            config.add_variable(
                name=var_name,
                expression=var_cfg[f"expression_{fw}"],
                binning=var_binning,
                x_title=rf"{var_cfg['label']} ({fw.upper()})",
                log_x=(var_cfg.get("scale", "linear") == "log"),
                aux={
                    "inputs": var_cfg.get("inputs", set()),
                },
            )
            variable_groups["synch"].add(var_name)

        # # add CF-UHH2 difference (TODO: figure out how to handle lambdas)
        # config.add_variable(
        #     name=f"s_{var_basename}_diff",
        #     expression=lambda x: var_cfg[f"expression_{fw}"],
        #     binning=var_cfg["diff_binning"],
        #     x_title=rf"$\Delta$({var_cfg['label']}) ({fw.upper()})",
        #     scale=var_cfg.get("scale", "linear"),
        #     aux={
        #         "inputs": var_cfg.get("inputs", set()),
        #     },
        # )

    # add variable groups aux if not exists
    if not config.has_aux("variable_groups"):
        config.x.variable_groups = {}

    # update variable groups aux
    config.x.variable_groups.update(variable_groups)
