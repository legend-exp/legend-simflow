(analysis)=

# Analysing the output

How to read the output of a production (see [](output.md)) and compare it with
the LEGEND-200 data, in Python. The paths are those of a production in
`/path/to/prod` of the `l200cfg09` configuration, with the default layout.

## The production configuration and the metadata

{func}`~legendsimflow.utils.load_simflow_config` reads the configuration file of
a production, with the defaults filled. Its `metadata` field is the LEGEND
metadata the production ran with (the copy at `paths.metadata`), as a
{class}`~legendmeta.LegendMetadata` object. This is a
[dbetto](https://dbetto.readthedocs.io) text database: a folder of YAML and JSON
files read as one nested dictionary, with directories, files and keys accessed
in the same way. Read all the metadata through it: the Simflow settings below
`simprod/config` (see [](meta.md)), the detector database, the channel maps and
the run information.

```python
from legendsimflow import load_simflow_config

config = load_simflow_config("/path/to/prod")
meta = config.metadata

meta.simprod.config.tier.evt[config.experiment].settings  # settings of evt
meta.simprod.config.tier.stp[config.experiment].simconfig.birds_nest_K40
meta.datasets.runinfo.p19.r006.phy.livetime_in_s  # s
meta.channelmap("20260123T200418Z").V00048A.daq.rawid
meta.hardware.detectors.germanium.diodes.V00048A.type  # "icpc"
```

A settings file holds only the keys that differ from the defaults, listed in
{mod}`legendsimflow.defaults`. The files written by the `par` step are not part
of the metadata (see {ref}`analysis-pars`).

## Reading LH5 files

```python
import lh5

f = "/path/to/prod/generated/tier/cvt/l200cfg09-birds_nest_K40-tier_cvt.lh5"
lh5.show(f)  # the content of the file
evt = lh5.read("evt", f, field_mask=["geds/energy", "geds/rawid"])
evt = evt.view_as("ak")  # as an awkward array
```

`field_mask` reads only the listed fields, which saves memory: the SiPM pulses
and the surface steps are large. `lh5.read` also takes a list of files, and
{class}`~lh5.io.iterator.LH5Iterator` reads large files in chunks. In Julia,
[LegendHDF5IO.jl](https://github.com/legend-exp/LegendHDF5IO.jl) reads the same
files with `lh5open`.

## Events

The event table of a simulation is in its `cvt` file, one row per event (see
{ref}`evt-tier`). Once viewed as an awkward array, it can be sliced and masked
like a NumPy array:

```python
evt = lh5.read_as("evt", f, "ak")
energy = evt.geds.energy[evt.geds.multiplicity == 1]  # keV, one list per event
```

The detectors are identified by their identifier in the simulated geometry
(`rawid`). `detector_uids` at the root of the file maps the detector names to
them:

```python
uids = {name: s.value for name, s in lh5.read("detector_uids", f).items()}
```

## Detector tables

The `hit` and `opt` files of a job hold one table per detector (see
{ref}`detector-tables-and-events`). To build events from them on the fly, use
the TCM of the file and {func}`reboost.read_hit_field_by_tcm`, which returns a
field of the detector tables with one list per event (see also the _remage_
[output documentation](https://remage.readthedocs.io/en/stable/manual/output.html)):

```python
import reboost

hit = "/path/to/prod/generated/tier/hit/birds_nest_K40/l200cfg09-birds_nest_K40-job_0000-tier_hit.lh5"
tcm = lh5.read_as("tcm", hit, "ak")
uids = reboost.get_remage_detector_uids(hit, lh5_table="hit")
energy = reboost.read_hit_field_by_tcm(tcm, hit, "energy", uids)
aoe = reboost.read_hit_field_by_tcm(tcm, hit, "psd/single_temp/aoe", uids)
```

(analysis-pars)=

## Parameters

The `par` step writes YAML files below `paths.pars` (see {ref}`par-output-files`
and {ref}`par-detinfo`). [dbetto](https://dbetto.readthedocs.io) reads the
folder as one database:

```python
from dbetto import TextDB

pars = TextDB("/path/to/prod/generated/pars")
pars.detinfo.usability["l200-p19-r006-phy"].V00048A  # "on"
pars["hpge/eresmod/l200-p19-r006-phy-model"].V00048A  # the energy resolution
pars["simstat/partitions_birds_nest_K40"]  # the run partitions
```

(analysis-pulses)=

## Pulse shapes

With the PSL-based PSD (see {ref}`hit-tier-settings`), the `par` step writes the
average current and charge pulses of the calibration data, the superpulses (see
{ref}`hpge-superpulses-extraction`), and the realistic pulse-shape libraries
(see {ref}`hpge-realistic-psl-extraction`).
{func}`~legendsimflow.superpulses.read_superpulses` reads the superpulses of a
detector, one per drift-time slice.
{func}`reboost.hpge.load_hpge_pulse_shape_libraries` reads the library of a
detector: one current pulse per point of an $(r, z)$ grid in the crystal, per
azimuth with respect to the crystal axes.

```python
import numpy as np
from legendsimflow.superpulses import read_superpulses

pars = "/path/to/prod/generated/pars/hpge"
sp_file = f"{pars}/superpulses/l200-p19-r006-phy-V00050A-superpulses.lh5"
psl = f"{pars}/psl/realistic/l200-p19-r006-phy-hpge-pulse-shape-lib.lh5"

superpulses = read_superpulses(sp_file, "V00050A")  # slice -> superpulse
sl, sp = next(iter(superpulses.items()))  # the longest drift times
sp.current_wf, sp.current_time_axis  # time in ns

lib = reboost.hpge.load_hpge_pulse_shape_libraries(psl, "V00050A", angles=[0])[0]
lib.waveforms, lib.t  # shape (n_r, n_z, n_samples), time in ns
drift_time = lh5.read("V00050A/drift_time_000_deg", psl).view_as("np")  # ns
```

The superpulse file is `V00050A-superpulses.lh5` when it is built from all runs
together (see {ref}`superpulses-settings-meta`). The library file holds all the
detectors of the run and can be very large: read one detector and only the
azimuths you need.

Both time axes have their zero at the current maximum. The mean of the library
pulses with a drift time in a slice is the simulated superpulse that the
electronics fit compares with the data (see {ref}`hpge-elecmod-extraction`).
Each grid point counts once, as in the fit:

```python
lo, hi = sl.drift_time_range  # ns
sim = np.nanmean(lib.waveforms[(drift_time >= lo) & (drift_time < hi)], axis=0)
```

The two pulses have different normalisations. Scale both to a maximum of 1
before comparing them.

(selections)=

## Selecting events

The functions of {mod}`legendsimflow.selections` apply the selections of the
LEGEND-200 analysis (multiplicity, LAr veto, A/E) to the event tables of the
simulations and of the data production alike, so that the two can be compared.
They are the functions the pdf step uses (see [](pdf.md)).

```python
import glob

from legendsimflow.selections import (
    aoe_flags,
    attach_modelable_flag,
    has_hpge_multiplicity,
    is_clean_event,
    is_spms_coincident,
    one_hit,
)

fields = [
    "trigger",
    "coincident",
    "spms/multiplicity",
    "spms/energy_sum",
    "geds/energy",
    "geds/rawid",
    "geds/quality",
    "geds/psd",
]
sevt = lh5.read_as("evt", f, "ak", field_mask=fields)
dfiles = glob.glob("/path/to/l200data/generated/tier/pet/phy/p19/r006/*-tier_pet.lh5")
devt = lh5.read_as("evt", dfiles, "ak", field_mask=fields)

devt = devt[is_clean_event(devt)]  # data only: no forced, pulser, muon events
devt = attach_modelable_flag(devt, config)  # which detectors the simulation models
hit_settings = config.metadata.simprod.config.tier.hit[config.experiment].settings
model = "pulse_lib" if hit_settings.get("simulate_psd_with_psl") else "single_temp"


def ss_energies(evt):
    """Energies of the single-site events that pass the LAr veto, in keV."""
    evt = evt[has_hpge_multiplicity(evt, 1)]  # one hit, in a detector that is on
    vetoed = is_spms_coincident(evt)  # LAr veto, thresholds of the evt settings
    valid, ss, high = aoe_flags(evt, model)  # high is None with the single template
    cut = ~vetoed & valid & ss
    if high is not None:
        cut = cut & ~high
    return one_hit(evt).energy[cut]


e_sim, e_data = ss_energies(sevt), ss_energies(devt)
```

- `has_hpge_multiplicity(evt, 2, 4)` selects a range of multiplicities. The
  energy threshold is that of the production (`energy_thr_kev` argument); the
  data holds only hits above 25 keV, so a lower one cannot be applied to it.
- The data production writes its event tables in the `pet` tier, or in `evt` in
  older productions.
- The LAr veto is recomputed on both tables with the same formula; the flag of
  the data production is not used.
- An event without a valid A/E (no usable PSD in the run, or a detector the
  simulation does not model) counts as background for an A/E cut.
- {func}`~legendsimflow.selections.pick_detector` keeps the events with a hit in
  one detector, with the fields of that hit only.

## The pdfs

The pdf files hold these selections already, as non-overlapping categories (see
{ref}`pdf-tier`). {func}`legendsimflow.pdf.read_pdf` sums the categories below
one or more paths, with `*` matching any category of a level. It reads the
simulated and the data files in the same way:

```python
from legendsimflow.pdf import read_pdf

sim = "/path/to/prod/generated/tier/pdf/l200cfg09-birds_nest_K40-tier_pdf.lh5"
data = "/path/to/prod/generated/tier/pdf/l200data/l200cfg09-l200data-tier_pdf.lh5"

h_sim = read_pdf("mul1/not_spms_coincident/is_ss_aoe", sim)  # as above
h_data = read_pdf("mul1/not_spms_coincident/is_ss_aoe", data)
read_pdf(["mul1/*/is_ss_aoe", "mul1/*/is_low_aoe"], sim, group="icpc")
```

The simulated counts divided by `nr_sim_events` give the probability per
simulated primary; the data file stores the livetime in `livetime_in_s`:

```python
n_sim = lh5.read("nr_sim_events", sim).value
livetime = lh5.read("livetime_in_s", data).value  # s
```

{func}`legendsimflow.pdf.merge_pdf_files` sums pdf files with the same
histograms, for example the data files of a subset of runs.
