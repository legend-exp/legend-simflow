# Pipeline steps

The Simflow produces a list of simulations, each labelled by a `simid`. A
simulation describes one physical process in one part of the experiment, for
example the decays of a radioactive isotope in the HPGe detector bulk. Every
simulation goes through the same chain of steps:

```{raw} html
:file: ../img/pipeline.svg
```

- [`vtx`](vtx.md) and [`stp`](stp.md) simulate the particle interactions in the
  experiment. `stp` runs _remage_ and records the energy depositions (steps) in
  the detectors. `vtx` runs before it only when the primary vertices come from
  an external generator.
- [`opt`](opt.md) and [`hit`](hit.md) turn the energy depositions into the
  observables of the SiPMs and of the HPGe detectors. They split the events
  among the LEGEND data taking runs (see {ref}`run-partitioning`) and apply the
  detector parameters of each run.
- [`par`](par.md) extracts those detector parameters, per run and per detector,
  mostly from the LEGEND data. All simulations share them.
- [`evt`](evt.md) combines the SiPM and HPGe observables into events, `cvt`
  concatenates the event files of a simulation, and [`pdf`](pdf.md) fills the
  histograms (the pdfs). The metadata that steers each step is described in
  [](meta.md), and the output files in [](output.md).

Each simulation is split into independent jobs with the same number of
primaries. The steps drawn as stacked cards write one file per job, so their
jobs can run in parallel. `cvt` joins them into a single file per simulation.

(detector-tables-and-events)=

## Detector tables and events

```{raw} html
:file: ../img/tcm.svg
```

Up to the `opt` and `hit` steps, the data are stored per detector, not per
event. _remage_ writes one table per detector in the `stp` file, with one row
per hit (the energy depositions of one event in one detector). `opt` and `hit`
keep this layout: they read one detector table at a time, in chunks of rows,
apply the model and parameters of that detector, and write one output table per
detector. A detector model does not need the hits of the other detectors, so
these steps never build events.

To keep track of which hits belong together, each `stp`, `opt` and `hit` file
holds a time-coincidence map (TCM, see {ref}`tcm-output`). For every event, the
TCM lists the detector tables and rows of its hits.

The `evt` step builds the events. It combines the TCMs of the `stp` and `opt`
files into one TCM that points to both the HPGe and the SiPM tables, and loops
over the events, in chunks. For each event, it uses the TCM to read the
observables of its hits from the detector tables, and it writes one row per
event. `cvt` keeps this layout, and `pdf` fills its histograms from these
events. Quantities that combine detectors, like the HPGe multiplicity or the LAr
veto, are therefore computed only from the `evt` step on.

(run-partitioning)=

## Run partitioning

```{raw} html
:file: ../img/runpart.svg
```

The detector parameters change during the data taking. From one LEGEND run to
the next, detectors are switched on or off, and their energy resolution and
pulse-shape response change. A simulation is therefore compared with a list of
data taking runs, and each run must be represented in the simulated events in
proportion to its livetime.

The Simflow splits the simulated events of each simulation into partitions, one
per run. It counts the events in the TCM of all jobs and gives each run a share
proportional to its livetime fraction, taken from
[legend-datasets](https://github.com/legend-exp/legend-datasets). Each partition
is a range of consecutive events. A partition can extend over two or more jobs,
and a job can hold two or more partitions: the first job might hold all the
events of the first run and the first events of the second run.

The `par` step computes the partitions once per simulation. The data of all the
later tiers are partitioned in this way, the hit-oriented and the event-oriented
ones alike. In the detector tables of `opt` and `hit`, the rows of the events of
a partition are processed with the parameters of its run: the usability of each
detector, and the energy resolution and pulse-shape parameters of the HPGe
detectors. In the event table of `evt`, the events of a partition are built with
the SiPM channels that are on in that run and, when random coincidences are
added, with coincidences drawn from the LEGEND data of that run. The partitions
stay in the same tables, split at the same events in every tier, and each row
records the period and run it belongs to.

Every simulated event belongs to exactly one partition, so no event is lost or
used twice. The simulated events can then be compared directly with the data of
the same runs.

## Steps

```{toctree}
:maxdepth: 1

vtx
stp
par
opt
hit
evt
pdf
```
