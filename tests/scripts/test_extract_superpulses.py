from __future__ import annotations

from pathlib import Path

import awkward as ak
import numpy as np
from legendmeta import LegendMetadata
from lgdo import Table, WaveformTable, lh5
from scipy.stats import norm

from legendsimflow.scripts.build_superpulses_from_data import main
from legendsimflow.superpulses import lookup_superpulse_inputs

l200data = Path(__file__).parent.parent / "l200data" / "v3.0.0"
dummyprod = Path(__file__).parent.parent / "dummyprod"


def test_make_data():
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
    energy = ak.unflatten(rng.uniform(25, 5000, size=size), np.ones(size, dtype=int))
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

    wf_win = wf_full[2625:4025]
    wf_presum = wf_full[::8][:-1]

    energy = ak.flatten(energy).to_numpy()

    noise_win = rng.normal(-2, 2, size=(size, len(wf_win)))
    noise_ps = rng.normal(-2, 2, size=(size, len(wf_presum)))

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

    Path(l200data / "generated" / "tier" / "raw" / "ssc" / "p16" / "r008").mkdir(
        parents=True, exist_ok=True
    )
    lh5.write(
        out,
        f"ch{rawid}",
        str(
            l200data
            / "generated"
            / "tier"
            / "raw"
            / "ssc"
            / "p16"
            / "r008"
            / "l200-p16-r008-ssc-20230322T170202Z-tier_raw.lh5"
        ),
        wo_mode="of",
    )


def test_lookup_inputs():
    meta = LegendMetadata(dummyprod / "inputs")
    raw_files, evt_files, dsp_config, tab_map = lookup_superpulse_inputs(
        l200data, meta, "l200-p16-r008-ssc", "V03422A", evt_tier_name="pet"
    )

    assert len(raw_files) == 1
    assert len(evt_files) == 1

    assert raw_files[0].exists()
    assert evt_files[0].exists()
    assert dsp_config.exists()

    assert raw_files[0].name == "l200-p16-r008-ssc-20230322T170202Z-tier_raw.lh5"
    assert evt_files[0].name == "l200-p16-r008-ssc-20230322T170202Z-tier_evt.lh5"

    assert tab_map["V03422A"] == 1108804


def test_cli(tmp_path):
    main(
        [
            "--l200data",
            str(l200data),
            "--meta",
            str(dummyprod / "inputs"),
            "--runid",
            "l200-p16-r008-ssc",
            "--detector",
            "V03422A",
            "--outdir",
            str(tmp_path / "outputs"),
        ]
    )
