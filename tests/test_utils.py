from __future__ import annotations

import numpy as np
from dbetto import AttrsDict
from lgdo import Array, Table

from legendsimflow import utils


def test_get_parameter_dict():
    popt = [100, 10, 60, 0.6, 100, 0.2, 60]

    popt_dict = utils._curve_fit_popt_to_dict(popt)
    assert len(popt_dict) == 7


def test_hash_string_int():
    int_hash = utils.string_to_remage_seed("blah blah legend simflow")
    assert isinstance(int_hash, int)
    assert int_hash >= 0


def test_add_field_string():
    tab = Table(size=10)

    utils.add_field_string("string", tab, "luigi was here")

    assert isinstance(tab.string, Array)
    assert np.all(tab.string.view_as("np").astype("str") == "luigi was here")


def test_get_dataflow_config(test_l200data):
    config = utils.lookup_dataflow_config(test_l200data / "v2.1.5")

    assert isinstance(config, AttrsDict)
    assert "paths" in config
    assert "$_" not in config.paths.tier
