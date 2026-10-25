from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from legendsimflow import opt_models

# one string of the default array, at (213.5, 0) mm
ONE_STRING = str(
    Path(__file__).parent
    / "dummyprod/inputs/simprod/config/geom/l1000dsg01-geom-config.yaml"
)


@pytest.fixture(scope="module")
def full_array(tmp_path_factory):
    path = tmp_path_factory.mktemp("geom") / "geom-config.yaml"
    path.write_text("executable: legend-pygeom-l1000\n")
    return str(path)


def test_l1000_strings(full_array):
    ids, centres, rots = opt_models._l1000_strings(full_array)
    assert len(ids) == 42
    assert ids[0] == 1
    assert centres[0] == pytest.approx([213.5, 0])
    assert rots[0] == pytest.approx(-np.pi)

    ids, centres, rots = opt_models._l1000_strings(ONE_STRING)
    assert ids.tolist() == [1]
    assert centres[0] == pytest.approx([213.5, 0])
    assert rots[0] == pytest.approx(-np.pi)


def test_nearest_string(full_array):
    _, centres, _ = opt_models._l1000_strings(full_array)
    x = np.array([centres[0, 0] + 30, centres[7, 0]])
    y = np.array([centres[0, 1], centres[7, 1] - 10])
    assert opt_models._nearest_string(x, y, centres).tolist() == [0, 7]


def test_parse_sipm():
    assert opt_models._parse_sipm("S0103B") == (1, 3, "B")
    with pytest.raises(ValueError, match="not a LEGEND-1000 SiPM"):
        opt_models._parse_sipm("S002")


def test_pen_plates(full_array):
    # the plates of string 1, of string 8 and a step between strings 1 and 2
    # but closer to 1
    _, centres, _ = opt_models._l1000_strings(full_array)
    x = np.array(
        [centres[0, 0], centres[7, 0], 0.6 * centres[0, 0] + 0.4 * centres[1, 0]]
    )
    y = np.array(
        [centres[0, 1], centres[7, 1], 0.6 * centres[0, 1] + 0.4 * centres[1, 1]]
    )
    z = np.zeros(3)

    def pen(sipm, geom_config=full_array):
        return opt_models.pen_plates(
            x, y, z, sipm=sipm, prob=0.02, geom_config=geom_config
        ).tolist()

    assert pen("S0101T") == [0.02, 0, 0.02]
    assert pen("S0103B") == [0.02, 0, 0.02]
    assert pen("S0802T") == [0, 0.02, 0]
    assert pen("all") == pytest.approx([0.12] * 3)
    # with a single string, all plates belong to it
    assert pen("S0101T", ONE_STRING) == [0.02] * 3


def test_fibers(full_array):
    # fibers of modules 1, 2 and 3 of string 1, in the middle and at the top
    # end of 1 m long fibers. Modules 2 and 3 from the geometry, module 1 at
    # 240 deg
    x = np.array([213.5 - 52.75 / 2, 251.05, 226.81, 251.05])
    y = np.array([-52.75 * np.sqrt(3) / 2, -37.05, 51.04, -37.05])
    z_top = 349.092
    z = np.array([z_top - 500] * 3 + [z_top])
    pars = {"amplitudes": [0.2, 0.3, 0.5], "absorption_lengths_in_m": [0.05, 0.2, 4]}

    def end_prob(dist_in_mm):
        return sum(
            a * np.exp(-dist_in_mm / (1000 * length))
            for a, length in zip(*pars.values(), strict=True)
        )

    mid, near, far = end_prob(500), end_prob(0), end_prob(1000)
    assert near == pytest.approx(1, rel=1e-3)

    for sipm, expected in [
        ("S0101T", [mid, 0, 0, 0]),
        ("S0101B", [mid, 0, 0, 0]),
        ("S0102T", [0, mid, 0, near]),
        ("S0102B", [0, mid, 0, far]),
        ("S0103B", [0, 0, mid, 0]),
        ("S0201T", [0, 0, 0, 0]),
        ("all", [2 * mid] * 3 + [near + far]),
    ]:
        prob = opt_models.fibers(
            x, y, z, sipm=sipm, **pars, fiber_length_in_m=1, geom_config=full_array
        )
        assert prob == pytest.approx(expected)
        # the strings around string 1 do not change the module of a step
        prob = opt_models.fibers(
            x, y, z, sipm=sipm, **pars, fiber_length_in_m=1, geom_config=ONE_STRING
        )
        assert prob == pytest.approx(expected)
