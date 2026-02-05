from __future__ import annotations

import awkward as ak
import pyg4ometry
import pytest
import reboost
from lgdo import lh5

from legendsimflow import reboost as rutils


def test_remage_hit_range(legend_testdata):
    f_stp = legend_testdata["remage/th228-full-optional-v0_13.lh5"]

    tcm = lh5.read_as("tcm", f_stp, library="ak")

    for det, uid in zip(
        ["det1", "det2", "scint1", "scint2", "optdet1", "optdet2"],
        [11, 12, 1, 2, 101, 102],
        strict=True,
    ):
        n_rows = lh5.read_n_rows(f"stp/{det}", f_stp)

        tcm = lh5.read_as("tcm", f_stp, "ak")

        assert rutils.get_remage_hit_range(tcm, det, uid, [0, len(tcm) - 1]) == (
            0,
            n_rows,
        )

        # divide into groups
        groups = [[0, 10], [11, 40], [41, 101], [102, len(tcm) - 1]]
        n = 0

        for group in groups:
            nen = rutils.get_remage_hit_range(tcm, det, uid, group)[1]

            if nen is not None:
                n += nen

        assert n == n_rows


def test_psd_stuff(legend_testdata):
    dt_map = {}
    for angle in ("000", "045"):
        dt_map[angle] = reboost.hpge.utils.get_hpge_rz_field(
            legend_testdata["lh5/V00048A-drift-time-maps-xtal-axes.lh5"],
            "V00048A",
            f"drift_time_{angle}_deg",
            bounds_error=False,
        )

    xloc = [
        [0.162, 0.162, 0.162, 0.162, 0.162, 0.162, 0.162, 0.162, 0.162],
        [0.183, 0.183],
        [0.197, 0.197, 0.197, 0.197, 0.197, 0.197, 0.197, 0.197],
        [0.201, 0.201, 0.201, 0.201, 0.201, 0.201, 0.201, 0.201],
        [0.213, 0.213, 0.213, 0.213, 0.213],
        [0.176, 0.176, 0.176, 0.176, 0.176, 0.176, 0.177],
        [0.193, 0.193, 0.193, 0.193, 0.193, 0.193],
        [0.201, 0.201, 0.201],
        [0.157, 0.157, 0.157, 0.157, 0.157],
        [0.164, 0.164, 0.164],
    ]

    yloc = [
        [0.107, 0.107, 0.107, 0.107, 0.107, 0.107, 0.107, 0.107, 0.107],
        [0.0906, 0.0906],
        [0.122, 0.122, 0.122, 0.122, 0.122, 0.122, 0.122, 0.122],
        [0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11],
        [0.131, 0.131, 0.131, 0.131, 0.131],
        [0.0946, 0.0946, 0.0946, 0.0947, 0.0947, 0.0947, 0.0947],
        [0.157, 0.157, 0.157, 0.157, 0.157, 0.157],
        [0.145, 0.145, 0.145],
        [0.124, 0.124, 0.124, 0.124, 0.124],
        [0.144, 0.144, 0.144],
    ]

    zloc = [
        [0.535, 0.535, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536, 0.536],
        [0.509, 0.509],
        [0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538, 0.538],
        [0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533, 0.533],
        [0.52, 0.52, 0.52, 0.52, 0.52],
        [0.506, 0.506, 0.506, 0.507, 0.507, 0.507, 0.507],
        [0.557, 0.557, 0.557, 0.557, 0.557, 0.557],
        [0.516, 0.516, 0.516],
        [0.519, 0.519, 0.519, 0.519, 0.519],
        [0.54, 0.54, 0.54],
    ]

    edep = ak.Array(
        [
            [130, 97.2, 233, 179, 133, 85.9, 173, 129, 331],
            [342, 125],
            [42.7, 75.9, 92.6, 261, 193, 49.6, 118, 249],
            [179, 35.6, 72.4, 174, 167, 241, 134, 126],
            [173, 37.3, 184, 156, 259],
            [13.9, 29.9, 91.1, 107, 106, 193, 70.8],
            [218, 13.2, 111, 111, 143, 118],
            [84.3, 56.9, 314],
            [63.7, 17.8, 92.4, 213, 186],
            [51, 95.3, 312],
        ]
    )

    chunk = ak.Array({"xloc": xloc, "yloc": yloc, "zloc": zloc})
    det_loc = pyg4ometry.gdml.Defines.Position(
        "Position",
        183.415,
        125.070,
        490.044,
        unit="mm",
    )

    dt = rutils.hpge_corrected_drift_time(
        chunk,
        dt_map,
        det_loc,
    )

    assert ak.all((dt > 0) & (dt < 3000))

    pars = {
        "amax": 852.0,
        "mu": 0.0999,
        "sigma": 54.5,
        "tail_fraction": 0.223,
        "tau": 507.0,
        "high_tail_fraction": 0.00611,
        "high_tau": 593.0,
    }
    amax = rutils.hpge_max_current(edep, dt, pars)

    assert ak.all((amax > 0) & (amax < 3000))


def test_cluster_photoelectrons_does_not_cross_subarrays():
    """Test that clustering does not merge elements across subarray boundaries."""
    times = ak.Array([[[0.0, 0.6], [0.7, 0.9]]])
    amps = ak.Array([[[1.0, 2.0], [3.0, 4.0]]])

    t_out, a_out = rutils.cluster_photoelectrons(times, amps, thr=1.0)

    assert ak.to_list(t_out) == [[[0.0], [0.7]]]
    assert ak.to_list(a_out) == [[[3.0], [7.0]]]


def test_cluster_photoelectrons_enforces_max_span():
    """Test that clusters respect the maximum time span threshold."""
    times = ak.Array([[0.0, 0.6, 1.1, 1.4, 2.3]])
    amps = ak.Array([[1.0, 2.0, 3.0, 4.0, 5.0]])

    t_out, a_out = rutils.cluster_photoelectrons(times, amps, thr=1.0)

    assert ak.to_list(t_out) == [[0.0, 1.1, 2.3]]
    assert ak.to_list(a_out) == [[3.0, 7.0, 5.0]]


def test_cluster_photoelectrons_empty_and_boundary():
    """Test clustering with empty arrays and exact boundary conditions."""
    times = ak.Array([[], [0.0, 1.0, 1.0001]])
    amps = ak.Array([[], [1.0, 2.0, 3.0]])

    t_out, a_out = rutils.cluster_photoelectrons(times, amps, thr=1.0)

    # [0.0, 1.0] spans exactly 1.0 -> same cluster; 1.0001 starts new
    assert ak.to_list(t_out) == [[], [0.0, 1.0001]]
    assert ak.to_list(a_out) == [[], [3.0, 3.0]]


def test_cluster_photoelectrons_mismatched_shapes():
    """Test that mismatched array shapes raise ValueError."""
    # Different nesting depths
    times_1d = ak.Array([0.0, 1.0, 2.0])
    amps_2d = ak.Array([[1.0, 2.0, 3.0]])

    with pytest.raises(ValueError, match="nesting depth"):
        rutils.cluster_photoelectrons(times_1d, amps_2d, thr=1.0)

    # Same nesting but different list lengths
    times = ak.Array([[0.0, 1.0], [2.0]])
    amps = ak.Array([[1.0], [2.0, 3.0]])

    with pytest.raises(ValueError, match="mismatched list lengths"):
        rutils.cluster_photoelectrons(times, amps, thr=1.0)
