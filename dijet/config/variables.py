# coding: utf-8

"""
Definition of variables.
"""

import order as od

from columnflow.util import maybe_import
from columnflow.columnar_util import EMPTY_FLOAT

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
