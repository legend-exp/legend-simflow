from __future__ import annotations

import copy
import shutil
from pathlib import Path

import awkward as ak
import legenddataflowscripts
import numpy as np
import pytest
import yaml
from dbetto import AttrsDict
from legendmeta import LegendMetadata
from legendtestdata import LegendTestData
from lgdo import Array, Table, WaveformTable, lh5
from pygeoml200 import core
from reboost.hpge.psd import _current_pulse_model as currmod
from scipy.stats import norm

from legendsimflow.utils import apply_path_defaults

l200data = Path(__file__).parent / "l200data" / "v3.0.0"
dummyprod = Path(__file__).parent / "dummyprod"
config_filename = dummyprod / "simflow-config.yaml"


@pytest.fixture(scope="session")
def legend_testdata():
    ldata = LegendTestData()
    ldata.checkout("68d8b49")
    return ldata


@pytest.fixture(scope="session", autouse=True)
def dummyprod_optmap(legend_testdata):
    """Copy the real optical map from legend_testdata into dummyprod.

    All configs referencing ``$_/inputs/simprod/l200cfg01-optmap-dummy.lh5``
    find a valid LH5 file.  The file is gitignored; this fixture is the sole
    source of truth for test runs.
    """
    src = Path(legend_testdata.get_path("remage/l200cfg01-optmap-dummy.lh5"))
    dst = dummyprod / "inputs/simprod/l200cfg01-optmap-dummy.lh5"
    shutil.copy2(src, dst)
    yield
    dst.unlink(missing_ok=True)


@pytest.fixture(scope="session")
def test_generate_gdml(config):
    geom_config = config.metadata.simprod.config.geom["legend-geom-config"]

    return core.construct(
        use_detailed_fiber_model=False, config=geom_config, public_geometry=True
    )


def make_config():
    with config_filename.open() as f:
        config = yaml.safe_load(f)

    legenddataflowscripts.subst_vars(config, var_values={"_": dummyprod})
    assert config is not None

    # convert all strings in the "paths" block to pathlib.Path
    def _make_path(d):
        for k, v in d.items():
            if isinstance(v, str):
                d[k] = Path(v)
            else:
                d[k] = _make_path(v)
        return d

    config["paths"] = _make_path(config["paths"])
    apply_path_defaults(config["paths"])

    metadata = LegendMetadata(dummyprod / "inputs")

    config["metadata"] = metadata
    config["_proctime"] = "now"

    return AttrsDict(config)


@pytest.fixture(scope="session")
def config():
    return make_config()


@pytest.fixture
def fresh_config():
    return make_config()


class mock_workflow_class:
    def __init__(self):
        self.overwrite_configfiles = [config_filename]


@pytest.fixture(scope="module")
def mock_workflow():
    return mock_workflow_class()


@pytest.fixture(scope="session")
def test_l200data():
    return Path(__file__).parent / "l200data"


@pytest.fixture(scope="session")
def test_make_ssc_data():
    # first an evt file, lets make 10000 events
    size = 100000

    rng = np.random.default_rng(seed=42)

    # basic info
    detector_name = ak.unflatten(np.full(size, b"V03422A"), np.ones(size, dtype=int))
    hit_idx = ak.unflatten(np.arange(size), np.ones(size, dtype=int))

    energy_sum = rng.uniform(5, 100, size=size)
    t0 = rng.uniform(0, 100, size=size)

    # coincident and trigger
    is_forced = np.full(size, False, dtype=bool)
    pulser = np.full(size, False, dtype=bool)

    # muon
    muon = np.full(size, False, dtype=bool)
    muon_offline = np.full(size, False, dtype=bool)

    # geds info
    is_bb_like = np.full(size, True, dtype=bool)
    multiplicity = np.full(size, 1, dtype=int)

    is_good_channel = ak.unflatten(
        np.full(size, True, dtype=bool), np.ones(size, dtype=int)
    )
    energy = ak.unflatten(rng.uniform(-25, 5000, size=size), np.ones(size, dtype=int))
    aoe = ak.unflatten(rng.uniform(-5, 5, size=size), np.ones(size, dtype=int))

    end_time = ak.unflatten(rng.uniform(100, 3000, size=size), np.ones(size, dtype=int))

    evts = ak.Array(
        {
            "trigger": {"is_forced": is_forced},
            "coincident": {"puls": pulser, "muon": muon, "muon_offline": muon_offline},
            "spms": {"energy_sum": energy_sum, "first_t0": t0},
            "geds": {
                "hit_idx": hit_idx,
                "quality": {
                    "is_good_channel": is_good_channel,
                    "is_bb_like": is_bb_like,
                },
                "energy": energy,
                "psd": {
                    "low_aoe": {"time": end_time, "value": aoe},
                    "high_aoe": {"value": aoe},
                },
                "multiplicity": multiplicity,
                "detector_name": detector_name,
            },
        }
    )

    tab = Table(evts)

    Path(l200data / "generated" / "tier" / "pet" / "ssc" / "p16" / "r008").mkdir(
        parents=True,
        exist_ok=True,
    )

    lh5.write(
        tab,
        "evt",
        str(
            l200data
            / "generated"
            / "tier"
            / "pet"
            / "ssc"
            / "p16"
            / "r008"
            / "l200-p16-r008-ssc-20230322T170202Z-tier_evt.lh5"
        ),
        wo_mode="of",
    )

    # now also make a raw file with the same number of entries, we can just make empty waveforms

    rawid = 1108804

    mu = 51 * 10**3
    sigma = 100

    # 1D grid
    t = np.linspace(0, 99968, 6249)

    # Gaussian PDF
    wf_full = np.cumsum(norm.pdf(t, loc=mu, scale=sigma))
    wf = copy.deepcopy(wf_full)

    for i in range(1, len(wf_full)):
        wf_full[i] = wf[i] - wf[i - 1] * (np.exp(16 / 450000)) + wf_full[i - 1]

    wf_win = wf_full[2625:4025]
    wf_presum = wf_full[::8][:-1] * 8

    energy = ak.flatten(energy).to_numpy()

    noise_win = rng.normal(0, 0.5, size=(size, len(wf_win)))
    noise_ps = rng.normal(0, 0.5, size=(size, len(wf_presum)))

    wf_win_tabl = np.vstack([wf_win] * size) * energy[:, np.newaxis] + noise_win
    wf_presum_tabl = np.vstack([wf_presum] * size) * energy[:, np.newaxis] + noise_ps

    wfs_win = WaveformTable(
        t0=42000.0, dt=16.0, t0_units="ns", dt_units="ns", values=wf_win_tabl
    )

    wfs_presum = WaveformTable(
        t0=0.0, t0_units="ns", dt=128.0, dt_units="ns", values=wf_presum_tabl
    )
    out = Table(
        {"raw": {"waveform_presummed": wfs_presum, "waveform_windowed": wfs_win}}
    )

    raw_ssc_dir = l200data / "generated" / "tier" / "raw" / "ssc" / "p16" / "r008"
    raw_ssc_dir.mkdir(parents=True, exist_ok=True)

    raw_ssc_file = raw_ssc_dir / "l200-p16-r008-ssc-20230322T170202Z-tier_raw.lh5"
    lh5.write(
        out,
        f"ch{rawid}",
        str(raw_ssc_file),
        wo_mode="of",
    )

    # hit-tier inputs needed by extract_hpge_current_pulse_model when l200data
    # is configured (noise-waveform selection uses cuspEmax_cal < threshold).
    aoe_classifier = rng.normal(0, 0.3, size=size).astype(np.float32)
    dt_eff = rng.uniform(100, 3000, size=size).astype(np.float32)
    # keep most entries well above threshold and only a small low-energy subset
    # (< 5 keV) so get_noise_maxima_and_sample does not scan all events.
    cusp_emax_cal = np.full(size, 100, dtype=np.float32)
    cusp_emax_cal[:200] = 1

    hit_tab = Table(
        {
            "cuspEmax_ctc_cal": Array(energy.astype(np.float32)),
            "AoE_Classifier": Array(aoe_classifier),
            "dt_eff": Array(dt_eff),
            "cuspEmax_cal": Array(cusp_emax_cal),
        }
    )

    hit_ssc_dir = l200data / "generated" / "tier" / "hit" / "ssc" / "p16" / "r008"
    hit_ssc_dir.mkdir(parents=True, exist_ok=True)

    lh5.write(
        hit_tab,
        f"ch{rawid}/hit",
        str(hit_ssc_dir / "l200-p16-r008-ssc-20230322T170202Z-tier_hit.lh5"),
        wo_mode="of",
    )
    return dummyprod


@pytest.fixture(scope="session")
def make_cal_data():
    """Generate synthetic calibration raw + hit LH5 files for p16/r007.

    The raw waveforms are scaled by the same energy array that is written into
    the hit tier, so the two tiers are internally consistent.  A subset of
    events is placed in the 1593 ± 5 keV fit window required by
    ``lookup_currmod_fit_data``.
    """
    size = 100000
    rawid = 1108804

    rng = np.random.default_rng(seed=142)

    # cal hit energies: fresh array with events in the 1593 keV fit window
    cal_energy = rng.uniform(25, 5000, size=size).astype(np.float32)
    cal_energy[:20] = np.linspace(1589, 1597, 20).astype(np.float32)

    aoe_classifier = rng.normal(0, 0.3, size=size).astype(np.float32)
    dt_eff = rng.uniform(100, 3000, size=size).astype(np.float32)

    hit_tab = Table(
        {
            "cuspEmax_ctc_cal": Array(cal_energy),
            "AoE_Classifier": Array(aoe_classifier),
            "dt_eff": Array(dt_eff),
            "cuspEmax_cal": Array(cal_energy),
        }
    )

    # raw waveforms based on the cal hit energies
    mu = 51 * 10**3
    sigma = 50
    t = np.linspace(0, 99968, 6249)

    curr = currmod(
        t,
        amax=1,
        mu=mu,
        sigma=sigma,
        tail_fraction=0.55,
        tau=150.0,
        high_tail_fraction=0.2,
        high_tau=80.0,
    )
    curr /= np.sum(curr) * 16
    wf_full = np.cumsum(curr)

    wf_win = wf_full[2625:4025]
    wf_presum = wf_full[::8][:-1] * 8

    noise_win = rng.normal(0, 0.5, size=(size, len(wf_win)))
    noise_ps = rng.normal(0, 0.5, size=(size, len(wf_presum)))

    wf_win_tabl = np.vstack([wf_win] * size) * cal_energy[:, np.newaxis] + noise_win
    wf_presum_tabl = (
        np.vstack([wf_presum] * size) * cal_energy[:, np.newaxis] + noise_ps
    )

    wfs_win = WaveformTable(
        t0=42000.0, dt=16.0, t0_units="ns", dt_units="ns", values=wf_win_tabl
    )
    wfs_presum = WaveformTable(
        t0=0.0, t0_units="ns", dt=128.0, dt_units="ns", values=wf_presum_tabl
    )
    raw_out = Table(
        {"raw": {"waveform_presummed": wfs_presum, "waveform_windowed": wfs_win}}
    )

    hit_cal_dir = l200data / "generated" / "tier" / "hit" / "cal" / "p16" / "r007"
    raw_cal_dir = l200data / "generated" / "tier" / "raw" / "cal" / "p16" / "r007"
    hit_cal_dir.mkdir(parents=True, exist_ok=True)
    raw_cal_dir.mkdir(parents=True, exist_ok=True)

    lh5.write(
        hit_tab,
        f"ch{rawid}/hit",
        str(hit_cal_dir / "l200-p16-r007-cal-20230322T170202Z-tier_hit.lh5"),
        wo_mode="of",
    )
    lh5.write(
        raw_out,
        f"ch{rawid}",
        str(raw_cal_dir / "l200-p16-r007-cal-20230322T170202Z-tier_raw.lh5"),
        wo_mode="of",
    )
    return l200data
