# `pdf` — PDF generation

The `pdf` step turns the events of a simulation into histograms of the HPGe
energy: the probability density functions (pdfs) of the background model. It can
also histogram the LEGEND-200 data in the same format, with the same event
selection, so that simulations and data can be compared directly. The settings
are described in {ref}`pdf-tier-settings` and the structure of the output in
{ref}`pdf-tier`; [](selections.md) shows how to apply the same selections by
hand.

## Normalisation

The histograms hold raw counts, one file per simulation, summed over all the
runs of its run list (the events are already shared among the runs in proportion
to their livetime, see {ref}`run-partitioning`). The file also stores the number
of simulated primaries. Dividing the counts by it gives the probability per
primary; multiplying by the number of decays of the source in the livetime of
the runs gives the expected counts.

## Event selection

Only events with HPGe hits above the energy threshold enter the histograms:

- events with exactly one hit (multiplicity 1) fill the energy histograms;
- events with exactly two hits (multiplicity 2) fill a two-dimensional histogram
  of the lower against the higher energy;
- events with a hit in a detector in anti-coincidence-only mode are left out.

The multiplicity-1 events are further split into categories, one histogram for
each combination:

- **LAr veto**: coincident with the liquid argon or not, as defined in
  [](evt.md);
- **A/E**: high A/E (only with the pulse-shape-library PSD), single-site, low
  A/E, or no valid A/E (the detector has no PSD model or no valid PSD in the
  run). Each event is in exactly one category;
- **energy in the FCCD layer** (optional): whether the event deposits more than
  a threshold in the layer between the n+ surface and the FCCD, before the
  dead-layer model. This needs the surface energy depositions of the `hit` step.

A category level is absent when its information is: the LAr veto without SiPMs,
the A/E without PSD simulation.

## Histograms

The energy histograms span 0 to 6000 keV in 1 keV bins; energies outside this
range are not stored. Each category has one histogram per detector group: all
detectors, and the groups defined in the settings by patterns on the detector
names. The categories select whole events, the groups select the detector of the
hit. The two-dimensional multiplicity-2 histogram has no categories and no
groups.

## LEGEND-200 data

When the LEGEND-200 data production is available, the step also histograms the
data of the runs in the run lists of all the simulations, one file per run and
one merged file. It applies the same selection functions as for the simulations,
with the settings of the simulations, plus the data-quality cuts of the physics
analysis: no forced triggers, no pulser or muon coincidences, and events that
pass the quality cuts. The LAr veto is recomputed from the SiPM multiplicity and
amplitude sum of the data with the thresholds of the simulation, rather than
taken from the data production. The data files store the livetime of the runs
instead of the number of primaries.

:::{note}

The merged data file covers the runs of all the simulations, while the pdf of a
simulation covers only the runs of its own run list. For a simulation with a
reduced run list, compare it with the data files of its runs.

:::
