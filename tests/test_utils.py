from __future__ import annotations

from legendsimflow import utils


def test_get_parameter_dict():
    popt = [100, 10, 60, 0.6, 100, 0.2, 60]

    popt_dict = utils._curve_fit_popt_to_dict(popt)
    assert len(popt_dict) == 7


def test_hash_string_int():
    int_hash = utils.string_to_int("blah blah legend simflow")
    assert isinstance(int_hash, int)
    assert int_hash >= 0
