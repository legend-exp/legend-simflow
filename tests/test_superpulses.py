from __future__ import annotations

from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import pytest
from lgdo import Array, Scalar, Struct, lh5

from legendsimflow.superpulses import (
    Slice,
    Superpulse,
    _get_nested_field,
    compute_chi2,
    get_wfs_for_slice,
    lookup_wfs_indices,
    plot_chi2_cut,
    plot_wfs_and_superpulse,
    read_superpulses,
    write_superpulses,
)

dsp = {
    "processors": {
        "wf_mean , wf_std, wf_slope, wf_intercept": {
            "description": "finds mean and rms of whole waveform",
            "function": "linear_slope_fit",
            "module": "dspeed.processors",
            "args": [
                "waveform",
                "wf_mean",
                "wf_std",
                "wf_slope",
                "wf_intercept",
            ],
            "unit": [
                "ADC",
                "ADC",
                "ADC",
                "ADC",
            ],
        },
        "wf_slope_diff , wf_slope_rms": {
            "description": (
                "finds mean and rms relative to linear fit of the waveform"
            ),
            "function": "linear_slope_diff",
            "module": "dspeed.processors",
            "args": [
                "waveform",
                "wf_slope",
                "wf_intercept",
                "wf_slope_diff",
                "wf_slope_rms",
            ],
            "unit": [
                "ADC",
                "ADC",
            ],
        },
        "tp_min, tp_max, wf_min, wf_max": {
            "description": (
                "find max and min of windowed waveform with corresponding time points"
            ),
            "function": "min_max",
            "module": "dspeed.processors",
            "args": [
                "waveform(len(waveform),'f')",
                "tp_min",
                "tp_max",
                "wf_min",
                "wf_max",
            ],
            "unit": [
                "ns",
                "ns",
                "ADC",
                "ADC",
            ],
        },
        "bl_mean_win , bl_std_win, bl_slope_win, bl_intercept_win": {
            "description": (
                "finds mean and rms of windowed waveform baseline (first 3us) "
                "as well as linear fit to this"
            ),
            "function": "linear_slope_fit",
            "module": "dspeed.processors",
            "args": [
                "waveform(len(waveform),'H')[0: 187]",
                "bl_mean_win",
                "bl_std_win",
                "bl_slope_win",
                "bl_intercept_win",
            ],
            "unit": [
                "ADC",
                "ADC",
                "ADC",
                "ADC",
            ],
        },
    }
}

# ===========================================================================
# Slice
# ===========================================================================


def test_slice_properties():
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))
    assert sl.energy_center == pytest.approx(1750.0)
    assert sl.drift_time_center == pytest.approx(1000.0)


def test_slice_str():
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))
    s = str(sl)
    assert "1500" in s
    assert "2000" in s
    assert "900" in s
    assert "1100" in s


def test_slice_is_hashable():
    sl1 = Slice((1500.0, 2000.0), (900.0, 1100.0))
    sl2 = Slice((1500.0, 2000.0), (900.0, 1100.0))
    d = {sl1: "value"}
    # same content -> same hash -> key lookup works
    assert d[sl2] == "value"


def test_slice_is_frozen():
    sl = Slice((1500.0, 2000.0), (900.0, 1100.0))
    with pytest.raises((AttributeError, TypeError)):
        sl.energy_range = (0.0, 1.0)


def test_read_nested_field():
    data = ak.Array({"geds": {"energy": [1, 2, 3]}, "spms_pe": [10, 2, 3]})

    assert ak.all(_get_nested_field(data, "geds/energy") == ak.Array([1, 2, 3]))
    assert ak.all(_get_nested_field(data, "spms_pe") == ak.Array([10, 2, 3]))


def test_lookup_wfs_indices(legend_testdata):
    ref_path = legend_testdata.get_path("lh5/")
    path = ref_path / Path("l200-p16-r008-ssc-20251006T205904Z-tier_evt.lh5")
    files = [str(path)]

    slices = [
        Slice(
            energy_range=(300.0, 2000.0),
            drift_time_range=(float(dt_start), float(dt_start + 50)),
        )
        for dt_start in range(900, 1900, 50)
    ]

    file_idx = lookup_wfs_indices(
        slices,
        evt_files=files,
        n_target=10,
        detector="V07302A",
        t0_field="spms/first_t0",
        end_time_field="geds/psd/drift_time",
    )

    assert isinstance(file_idx, list)
    assert len(file_idx) == len(slices)

    for fi in file_idx:
        assert len(fi.file_idx) == len(fi.hit_idx) == fi.n_sel


def test_get_wfs_for_slice(legend_testdata):
    ref_path = legend_testdata.get_path("lh5/prod-ref-l200/")
    path = ref_path / Path("generated/tier/raw/phy/p03/r001")

    files = [str(p) for p in path.glob("*.lh5")]

    wf_data = get_wfs_for_slice(
        files,
        lh5_group="ch1084804/raw",
        hit_indices=[0, 1],
        file_indices=[0, 0],
        dsp_config=dsp,
        charge_output="waveform",
        current_output="waveform",
        energy_output="wf_max",
        align=None,
    )

    assert isinstance(wf_data.charge_times, np.ndarray)
    assert isinstance(wf_data.current_times, np.ndarray)
    assert isinstance(wf_data.charge_wfs, np.ndarray)
    assert isinstance(wf_data.current_wfs, np.ndarray)

    assert wf_data.charge_wfs.shape == (2, len(wf_data.charge_times))
    assert wf_data.current_wfs.shape == (2, len(wf_data.current_times))


# ===========================================================================
# Superpulse construction and validation
# ===========================================================================


@pytest.fixture(scope="session")
def test_make_superpulse():
    n_charge = 100
    n_current = 80
    energy_range = (1500.0, 2000.0)
    drift_time_range = (900.0, 1100.0)
    detector = "V03422A"
    n_preliminary = 50
    n_final = 40

    """Return a minimal but valid Superpulse for reuse across tests."""
    rng = np.random.default_rng(0)
    return Superpulse(
        charge_wf=rng.random(n_charge),
        current_wf=rng.random(n_current),
        charge_time_axis=np.linspace(-500.0, 500.0, n_charge),
        current_time_axis=np.linspace(-500.0, 500.0, n_current),
        slice=Slice(energy_range, drift_time_range),
        detector=detector,
        n_events_preliminary=n_preliminary,
        n_events_final=n_final,
    )


def test_superpulse_init_valid(test_make_superpulse):
    assert test_make_superpulse.detector == "V03422A"
    assert test_make_superpulse.n_events_preliminary == 50
    assert test_make_superpulse.n_events_final == 40


def test_superpulse_mismatched_charge_length():
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="charge_wf"):
        Superpulse(
            charge_wf=rng.random(100),
            current_wf=rng.random(80),
            charge_time_axis=np.linspace(-500, 500, 99),  # wrong length
            current_time_axis=np.linspace(-500, 500, 80),
            slice=Slice((1500.0, 2000.0), (900.0, 1100.0)),
            detector="V03422A",
            n_events_preliminary=50,
            n_events_final=40,
        )


def test_superpulse_mismatched_current_length():
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="current_wf"):
        Superpulse(
            charge_wf=rng.random(100),
            current_wf=rng.random(80),
            charge_time_axis=np.linspace(-500, 500, 100),
            current_time_axis=np.linspace(-500, 500, 50),  # wrong length
            slice=Slice((1500.0, 2000.0), (900.0, 1100.0)),
            detector="V03422A",
            n_events_preliminary=50,
            n_events_final=40,
        )


def test_superpulse_final_exceeds_preliminary():
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="n_events_final"):
        Superpulse(
            charge_wf=rng.random(100),
            current_wf=rng.random(80),
            charge_time_axis=np.linspace(-500, 500, 100),
            current_time_axis=np.linspace(-500, 500, 80),
            slice=Slice((1500.0, 2000.0), (900.0, 1100.0)),
            detector="V03422A",
            n_events_preliminary=10,
            n_events_final=20,  # more than preliminary
        )


def test_superpulse_repr(test_make_superpulse):
    r = repr(test_make_superpulse)
    assert "V03422A" in r
    assert "40/50" in r


# ===========================================================================
# Superpulse.to_lgdo
# ===========================================================================


def test_to_lgdo_returns_struct(test_make_superpulse):
    result = test_make_superpulse.to_lgdo()
    assert isinstance(result, Struct)


def test_to_lgdo_fields_present(test_make_superpulse):
    result = test_make_superpulse.to_lgdo()
    expected_keys = {
        "charge_wf",
        "current_wf",
        "charge_time_axis",
        "current_time_axis",
        "drift_time_center",
        "drift_time_lo",
        "drift_time_hi",
        "energy_lo",
        "energy_hi",
        "detector",
        "n_events_preliminary",
        "n_events_final",
    }
    assert expected_keys == set(result.keys())


def test_to_lgdo_array_types(test_make_superpulse):
    result = test_make_superpulse.to_lgdo()
    assert isinstance(result["charge_wf"], Array)
    assert isinstance(result["current_wf"], Array)
    assert isinstance(result["charge_time_axis"], Array)
    assert isinstance(result["current_time_axis"], Array)


def test_to_lgdo_scalar_types(test_make_superpulse):
    result = test_make_superpulse.to_lgdo()
    for key in (
        "drift_time_center",
        "drift_time_lo",
        "drift_time_hi",
        "energy_lo",
        "energy_hi",
        "detector",
        "n_events_preliminary",
        "n_events_final",
    ):
        assert isinstance(result[key], Scalar), f"{key} should be Scalar"


def test_to_lgdo_time_axis_units(test_make_superpulse):
    result = test_make_superpulse.to_lgdo()
    assert result["charge_time_axis"].attrs.get("units") == "ns"
    assert result["current_time_axis"].attrs.get("units") == "ns"


def test_to_lgdo_waveform_values(test_make_superpulse):
    result = test_make_superpulse.to_lgdo()
    np.testing.assert_array_equal(
        result["charge_wf"].view_as("np"), test_make_superpulse.charge_wf
    )
    np.testing.assert_array_equal(
        result["current_wf"].view_as("np"), test_make_superpulse.current_wf
    )


# ===========================================================================
# chi2
# ===========================================================================


def test_compute_chi2_returns_array():
    n_events = 20
    bl_std = np.ones(n_events) * 2.0
    cuspEmax = np.ones(n_events) * 1000.0
    wfs = np.random.default_rng(0).random((n_events, 100))

    sp = Superpulse(
        charge_wf=np.random.default_rng(0).random(100),
        current_wf=np.random.default_rng(0).random(80),
        charge_time_axis=np.linspace(-500.0, 500.0, 100),
        current_time_axis=np.linspace(-500.0, 500.0, 80),
        slice=Slice((1500.0, 2000.0), (900.0, 1100.0)),
        detector="V03422A",
        n_events_preliminary=n_events,
        n_events_final=n_events,
    )
    chi2 = compute_chi2(wfs, sp, bl_std=bl_std, cuspEmax=cuspEmax)
    assert isinstance(chi2, np.ndarray)
    assert chi2.shape == (n_events,)


# ===========================================================================
# write_superpulses
# ===========================================================================


def test_write_superpulses_creates_file(test_make_superpulse, tmp_path):
    output_path = str(tmp_path / "test_superpulses.lh5")
    superpulses = {test_make_superpulse.slice: test_make_superpulse}
    write_superpulses(superpulses, output_path, detector="V03422A")
    assert (tmp_path / "test_superpulses.lh5").exists()


def test_write_superpulses_lh5_structure(test_make_superpulse, tmp_path):
    output_path = str(tmp_path / "test_superpulses.lh5")
    write_superpulses(
        {test_make_superpulse.slice: test_make_superpulse},
        output_path,
        detector="V03422A",
    )

    # The expected group path is V03422A/dt_900_1100_ns
    result = lh5.read("V03422A/dt_900_1100_ns", output_path)
    assert isinstance(result, Struct)

    assert "charge_wf" in result
    assert "current_wf" in result
    assert "drift_time_center" in result
    assert "energy_lo" in result
    assert "energy_hi" in result
    assert "n_events_preliminary" in result
    assert "n_events_final" in result


def test_read_superpulses(test_make_superpulse, tmp_path):
    output_path = str(tmp_path / "test_superpulses.lh5")
    write_superpulses(
        {test_make_superpulse.slice: test_make_superpulse},
        output_path,
        detector="V03422A",
    )

    # Read back the file and check contents
    result = read_superpulses(output_path, detector="V03422A")
    assert isinstance(result, dict)
    assert all(isinstance(v, Superpulse) for v in result.values())


def test_plot(test_make_superpulse):
    # Just test that the plotting functions run without error on a valid superpulse
    # Detailed tests of the plot contents would require image comparison which is beyond scope here

    fig, _ = plot_wfs_and_superpulse(
        test_make_superpulse.charge_time_axis,
        test_make_superpulse.current_time_axis,
        np.random.default_rng(0).random(
            (
                test_make_superpulse.n_events_preliminary,
                len(test_make_superpulse.charge_time_axis),
            )
        ),
        np.random.default_rng(0).random(
            (
                test_make_superpulse.n_events_preliminary,
                len(test_make_superpulse.current_time_axis),
            )
        ),
        test_make_superpulse,
    )
    plt.close(fig)

    chi2_values = (
        np.random.default_rng(0).random(test_make_superpulse.n_events_preliminary) * 5
    )  # some above and some below threshold

    fig, _ = plot_chi2_cut(
        chi2_values,
        5,
        test_make_superpulse.charge_time_axis,
        np.random.default_rng(0).random(
            (
                test_make_superpulse.n_events_preliminary,
                len(test_make_superpulse.charge_time_axis),
            )
        ),
        test_make_superpulse,
        "charge",
    )
    plt.close(fig)
