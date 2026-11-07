# `evt` and `cvt` — event building

The `evt` step combines the HPGe and SiPM observables of the `hit` and `opt`
steps into events, and the `cvt` step concatenates the event files of all the
jobs of a simulation into one. From here on the data are organized by event (see
{ref}`detector-tables-and-events`). The settings are described in
{ref}`evt-tier-settings-meta` and the output fields in {ref}`evt-tier`.

## Events

An event is a row of the TCM: the hits of one simulated event that are close in
time (see {ref}`tcm-output`). The `evt` step combines the TCMs of the `stp` and
`opt` files into one, which points to the HPGe tables of the `hit` file and to
the SiPM tables of the `opt` file, and writes one row per TCM row. An event is
therefore not a primary: a primary whose energy depositions are far apart in
time, as in a decay chain, gives several events, and a primary without energy in
any detector gives none.

Each event records the identifier, period and run of its first hit, and its time
relative to the start of the simulated primary.

## HPGe detectors

The event keeps the HPGe hits above the energy threshold of the data analysis
(for LEGEND-200, 25 keV) in detectors that are not off in the run. Hits in
detectors that are on and hits in detectors in anti-coincidence-only mode (AC)
are both kept, and each hit records whether its detector is on. From these, the
event gets:

- the energy and the PSD observables of each kept hit, as computed by `hit`;
- the multiplicity, the number of kept hits, AC detectors included;
- the energy sum, over the kept hits of the detectors that are on only.

An event without kept HPGe hits stays in the output, with multiplicity zero.

With the surface depositions of the `hit` step stored (see
{ref}`hit-tier-settings`), each kept hit also carries its energy depositions
near the n+ surface, as deposited, with their distance to the surface. They make
it possible to recompute the energy of the hit for another model of the surface
layer, for example another FCCD or dead layer fraction, with
{func}`reboost.math.piecewise_linear_activeness`.

:::{warning}

The hits are selected with the energy threshold computed with the nominal
surface model. A hit below the threshold for the nominal model can pass it for a
thinner FCCD, but it is missing from the event. Energies recomputed for another
model are therefore biased near the threshold, and so is the multiplicity of the
events.

:::

## SiPMs

The event lists the photoelectron pulses of every SiPM channel that is not off
in the run, above the SiPM threshold of the analysis, with their amplitudes in
photoelectrons and their times relative to the event time. From these, the event
gets the SiPM multiplicity, the number of channels with at least one pulse, and
the sum of the amplitudes in photoelectrons.

The event is coincident with the liquid argon (the LAr veto of the data
analysis) when the SiPM multiplicity or the amplitude sum reaches a threshold;
LEGEND-200 uses 4 channels or 4 photoelectrons, LEGEND-1000 1 and 1. The event
is coincident with the HPGe array when it has at least one kept HPGe hit.

## Random coincidences

The SiPMs of the real experiment see light that has nothing to do with the
event: argon-39 decays, dark counts, other backgrounds. A simulated event alone
does not contain it, so the LAr veto would be too lenient. The `evt` step can
add it from the data: it overlays on each simulated event the SiPM pulses of a
random-coincidence event taken from the LEGEND-200 data of the run.

The random-coincidence events come either from the forced-trigger and pulser
events of the physics data, and from the time before the trigger of the
HPGe-triggered events, cut into windows of 6 µs, or from a dedicated
noise-trigger data stream. They are used in the order they are stored, each once
until they run out, and are added to the multiplicity, the amplitude sum and the
LAr veto of the event. Their pulses are also stored separately from the
simulated ones.

:::{note}

The random-coincidence events are taken from the start of the data again for
every job and every run partition. The jobs of a simulation therefore share the
same random-coincidence events, which matters only when a job holds more
simulated events than there are random-coincidence events in the data.

:::

## Run partitioning

The events of each run partition (see {ref}`run-partitioning`) are built with
the SiPM channels that are on in that run, and with random coincidences from
that run (or from the run configured for it). The partitions are written one
after the other, in the order of the runs, into the same table.

## Concatenation

The `cvt` step concatenates the event tables of all the jobs of a simulation, in
the order of the jobs, into one file. It does not select, reorder or reweight
events. The number of simulated primaries of the jobs is summed, so that the
file can be normalised (see [](pdf.md)). The content of the file is described in
{ref}`cvt-tier`.

:::{note}

The fields that point back to rows of the `hit` and `opt` files of a job keep
their per-job values, so they are not valid in the concatenated file.

:::
