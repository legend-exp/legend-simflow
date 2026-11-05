(selections)=

# Event selections on simulations and data

The pdf tier histograms the simulated events after the selections of the
LEGEND-200 analysis: multiplicity, LAr veto and A/E (see {ref}`pdf-tier`). The
functions of {mod}`legendsimflow.selections` apply the same selections to the
event tables of the simulations (`cvt` tier) and of the data (`evt` or `pet`
tier of the data production). The call is the same on both tables: the functions
read the fields of the table they are given (see the field table in
{ref}`pdf-tier-data`). This page goes through the steps on one simulation and
one run of data.

## The production

Everything comes from the production: the A/E model simulated in the `hit` tier,
the list of runs and the detectors the simulation models.
{func}`legendsimflow.utils.load_simflow_config` reads its `simflow-config.yaml`:

```python
import awkward as ak
import lh5

from legendsimflow import load_simflow_config, patterns, utils
from legendsimflow.metadata import get_tier_settings
from legendsimflow.selections import (
    aoe_flags,
    aoe_model,
    attach_modelable_flag,
    has_hpge_multiplicity,
    is_clean_event,
    is_spms_coincident,
    one_hit,
)

config = load_simflow_config("/path/to/the/production")
model = aoe_model(
    get_tier_settings(config, "hit")
)  # "pulse_lib", "single_temp" or None
```

## Loading the events

The simulated events of a `simid` are in one `cvt` file. The data events of a
run are in the `evt` (or `pet`) files of the data production at
`paths.l200data`. Read only the fields the selections need, the SiPM pulses and
the surface steps of the simulation are large:

```python
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

sevt = lh5.read(
    "evt",
    patterns.output_tier_cvt_filename(config, simid="l200a-hpge-bulk-k40"),
    field_mask=fields,
).view_as("ak")

l200data = config.paths.l200data
files = utils.lookup_evt_files(
    l200data, "l200-p16-r000-phy", utils.get_evt_tier_name(l200data)
)
devt = lh5.read("evt", files, field_mask=fields).view_as("ak")
```

Several runs of data can be concatenated with {func}`ak.concatenate`: the events
keep their run in `trigger.period` and `trigger.run`, which the A/E step below
uses.

## Data cleaning

The data has forced triggers, pulser and muon events and hits of bad quality
that the simulation does not have.
{func}`~legendsimflow.selections.is_clean_event` removes them, on the data only:

```python
devt = devt[is_clean_event(devt)]
```

After this line the two tables go through the same steps.

## Multiplicity

A hit is an HPGe signal above 25 keV, the threshold of the data production and
of the `evt` tier (`geds_energy_thr_kev`).
{func}`~legendsimflow.selections.has_hpge_multiplicity` selects the events by
the number of hits, and requires every hit to be in an ON detector (an event
with an AC detector above threshold is never selected):

```python
sevt = sevt[has_hpge_multiplicity(sevt, 1)]
devt = devt[has_hpge_multiplicity(devt, 1)]
```

A range of multiplicities is `has_hpge_multiplicity(evt, 2, 4)`, both included.
From here on each event has one hit, and
{func}`~legendsimflow.selections.one_hit` returns its HPGe record:

```python
energy = one_hit(sevt).energy  # keV
```

Every function that looks for the hits takes the threshold as the
`energy_thr_kev` argument, for a production with another `geds_energy_thr_kev`.
The data production stores only the hits above 25 keV: a lower threshold cannot
be applied to the data, a higher one applies to both tables in the same way.

## LAr veto

{func}`~legendsimflow.selections.is_spms_coincident` flags the events in which
the SiPMs see enough light: at least 4 channels with a pulse, or at least 4
photoelectrons in total, the `lar_veto_multiplicity_thr` and
`lar_veto_energy_sum_pe_thr` settings of the `evt` tier, to pass when they
differ. The same formula runs on both tables; the flag of the data production is
not used, since its definition can differ:

```python
s_vetoed = is_spms_coincident(sevt)
d_vetoed = is_spms_coincident(devt)
```

## A/E

A valid A/E needs a usable PSD in the run (`psd.is_good`, the same information
in both tables) and an A/E the simulation models. The simulation stores this per
hit (`has_aoe`). The data does not know which detectors the simulation models,
and the set changes with the run:
{func}`~legendsimflow.selections.attach_modelable_flag` adds the flag to the
data events, from the modelable-HPGe cache of the `par` step (see
{ref}`pdf-tier-data`) and the run of each event:

```python
devt = attach_modelable_flag(devt, config)
```

{func}`~legendsimflow.selections.aoe_flags` then returns, for the hit of each
event, whether the A/E is valid, above the low-side cut (single-site) and above
the high-side cut. The last is `None` with the single template, which has no
high-side cut:

```python
s_valid, s_ss, s_high = aoe_flags(sevt, model)
d_valid, d_ss, d_high = aoe_flags(devt, model)
```

An event in a detector that is not modeled, or whose PSD is not usable, is
invalid on both sides and treated as background by an A/E cut.

## One detector

{func}`~legendsimflow.selections.pick_detector` keeps the events with a hit in
one detector and, for those, the fields of that hit alone: one energy per event,
one list of surface steps per event. The hits are taken as stored, with no
energy threshold, and the other tables of the events stay as they are. The table
argument picks a SiPM channel instead:

```python
v = pick_detector(sevt, 1104000)
v.geds.energy  # keV, one value per event
v.geds.surface.distance_to_nplus  # mm, one list of steps per event
pick_detector(sevt, 1057600, "spms").spms.energy  # photoelectrons of one channel
```

## The spectrum after all cuts

The flags are boolean arrays, one per event, so any selection is their product.
The spectrum of the single-site events that pass the LAr veto and the high-side
cut:

```python
import hist


def spectrum(evt, vetoed, valid, ss, high):
    energy = one_hit(evt).energy
    h = hist.Hist.new.Reg(6000, 0, 6000, name="energy").Double()
    return h.fill(energy[~vetoed & valid & ss & ~high])


h_sim = spectrum(sevt, s_vetoed, s_valid, s_ss, s_high)
h_data = spectrum(devt, d_vetoed, d_valid, d_ss, d_high)
```

The pdf tier stores the same selection as the category
`mul1/not_spms_coincident/is_ss_aoe`, which {func}`legendsimflow.pdf.read_pdf`
reads, and the data histograms of the pdf tier are built with the functions of
this page (see {ref}`pdf-tier-data`).
