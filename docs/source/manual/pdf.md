# `pdf` — PDF generation

The `pdf` step turns the events of a simulation into histograms of the HPGe
energy: the probability density functions (pdfs) of the background model. It can
also histogram the LEGEND-200 data in the same format, with the same event
selection, so that simulations and data can be compared directly. The settings
are described in {ref}`pdf-tier-settings` and the structure of the output in
{ref}`pdf-tier`; [](analysis.md) shows how to apply the same selections by hand.

## Normalisation

The events of a simulation are shared among the runs of its run list in
proportion to their livetime (see {ref}`run-partitioning`). The step fills one
file per run with the events of that run, and sums these files into one file per
simulation. The histograms hold raw counts. Each file also stores the number of
simulated primaries it covers. The partitioning counts only the events with
energy in a detector, so the primaries of the simulation are shared among the
runs in proportion to their events, that is to their livetime (see
{func}`legendsimflow.partitioning.primaries_per_run`). The file of the
simulation holds the total. Dividing the counts by it gives the probability per
primary; multiplying by the number of decays of the source in the livetime of
the runs gives the expected counts.

One job per simulation reads its `evt` files once, in the order of the jobs, and
writes the file of each run as soon as all its events are read. The files of a
simulation can also be built by hand:

```console
> pixi run tier-pdf --simid birds_nest_K40 \
    --evt-files generated/tier/evt/birds_nest_K40/*.lh5 \
    --simstat-part-file generated/pars/simstat/partitions_birds_nest_K40.yaml \
    --pdf-file pdf.lh5 --run-dir pdf-runs --simflow-config simflow-config.yaml
```

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
instead of the number of primaries (see {ref}`pdf-tier-data`). One run can also
be histogrammed by hand, into the file the workflow writes:

```console
> pixi run tier-pdf-l200data --runid l200-p16-r000-phy --simflow-config simflow-config.yaml
```

:::{note}

The merged data file covers the runs of all the simulations, while the pdf of a
simulation covers only the runs of its own run list. For a simulation with a
reduced run list, compare it with the data files of its runs. The files of one
run, simulated and data, can always be compared directly.

:::
