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
