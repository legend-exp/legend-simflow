# `par` — parameter extraction

The `par` step builds the detector response model that the `hit` step applies to
the simulated energy depositions: for every run in the run list and every HPGe
detector, the energy resolution, the resolution and the cuts of the A/E
pulse-shape discrimination (PSD), and the pulse shapes from which the A/E of a
simulated event is computed. It also records which detectors are on and which
can get a pulse-shape model, and computes the run partitions (see
{ref}`run-partitioning`). Everything is computed once per production and shared
by all simulations. The output files are listed in {ref}`par-output-files` and
the rules in the [rule reference](../api/snakemake_rules.md).

All the parameters on this page concern the HPGe detectors. The SiPMs need no
run-dependent parameters besides their usability, which the step caches together
with that of the HPGe detectors in the detector-info files (see
{ref}`par-detinfo`).

(par-collection-model)=

## Where the parameters come from

The figure below shows the parameters, where they come from and how they reach
the `hit` step. Dashed arrows are inputs from the data production that the
metadata can replace, as explained below.

```{raw} html
:file: ../img/pardag.svg
```

Most parameters are measured ones. The Simflow reads them from the LEGEND data
production of the same runs (see [](setup.md)), where the calibration analysis
has fitted them to the calibration data. Each parameter can also be set in the
Simflow metadata, per detector or with a `default` for all detectors (see the
metadata directories in [](meta.md)):

- without metadata, every detector with a calibration entry in the data
  production gets its measured parameters;
- with metadata but no `default`, the measured parameters are collected as
  above, and the listed detectors are then overridden;
- with a `default`, the metadata replaces the data production entirely: every
  HPGe detector in the channel map gets the `default` or its own entry, and the
  data production is not needed.

The choice is made per parameter, so the energy resolution can come from the
metadata while the PSD cuts come from the data. Two parameters, the
current-pulse template and the electronics response, are fitted by the Simflow
itself rather than collected: for them, the metadata is used only when it gives
a `default`. The `hit` step checks that every detector it simulates has the
parameters it needs, and falls back to defaults or stops (see
{ref}`build-tier-hit-hpge`).

(hpge-modeling-criteria)=

## Which detectors get a pulse-shape model

The energy resolution, the A/E resolution and the PSD cuts are collected for
every detector that has them. The pulse-shape simulation, that is the
current-pulse template, the drift-time maps and the pulse-shape libraries, is
instead done only for the detectors that can be modelled. An HPGe detector can
be modelled in a run when:

- it is on in that run;
- it is not excluded by hand for that run in the skip list of the metadata (see
  {ref}`skip-metadata-dir`; the reason given there is documentation only, the
  workflow does not read it);
- it has an operational voltage;
- its crystal metadata gives an impurity profile;
- it is operated at least 100 V above its depletion voltage, taken from the
  detector metadata. Each experiment can change this margin (see
  {ref}`modeling-settings-meta`).

When the production uses impurity profiles tuned on data (see
[](impurity-tuning.md)), the last condition is replaced: the detector must be
listed in the file of tuned profiles, whose depletion voltage then replaces the
one of the metadata in the simulation. The margin above depletion is not
checked. The tuning production itself selects the detectors in its own way,
described on the same page.

The outcome is cached in the detector-info files ({ref}`par-detinfo`).

```{toctree}
:hidden:

impurity-tuning
```

(hpge-eresmod-extraction)=

## Energy resolution

The energy resolution of a detector follows the model of the data calibration,

$$
\text{FWHM}(E) = \sqrt{a + b E},
$$

with $E$ in keV. The parameters $a$ and $b$ are taken, per run and per detector,
from the energy calibration of the data production (the fit of the calibration
peaks after the charge-trapping correction). The `hit` step smears the energy
deposited in a detector with a Gaussian of this width. Overrides go in
{ref}`eresmod-metadata-dir`.

(hpge-aoeresmod-extraction)=

## A/E resolution

The A/E resolution is the width of the A/E distribution of single-site events as
a function of the energy,

$$
\sigma_{A/E}(E) = \sqrt{a + (b/E)^c},
$$

with $E$ in keV, taken per run and per detector from the A/E calibration of the
data production. The `hit` step expresses the simulated A/E in units of it, as
the data analysis does. Overrides go in {ref}`aoeresmod-metadata-dir`.

(hpge-psdcuts-extraction)=

## PSD cuts

The PSD cuts are the two thresholds of the A/E classifier of the data analysis,
in units of the A/E resolution: the low cut, which rejects multi-site events,
and the high cut, which rejects events on the detector surface. They are taken
per run and per detector from the A/E calibration of the data production. The
`hit` step applies them to the simulated classifier, so that the simulated
events are selected like the data. Overrides go in {ref}`psdcuts-metadata-dir`.

(hpge-currmod-extraction)=

## Current-pulse template

The single-template PSD simulation (the default, see {ref}`hit-tier-settings`)
represents the current pulse of a single-site event in a detector with one
analytic template. The Simflow fits it, per run and per detector, to one
measured pulse: among the single-site events of the calibration data in a 10 keV
window around the 1593 keV double-escape peak of Tl-208, the one with the drift
time closest to the median. For physics runs the calibration data are those of
the reference calibration run. The smoothed current waveform of that event is
fitted with a Gaussian core and an exponential tail on each side, seven
parameters in all. The step also records the mean A/E that the template gives at
this energy and the A/E resolution due to the electronic noise, measured on
baseline waveforms of the same data. The `hit` step uses the template and the
drift-time map of the detector to compute the A/E of each simulated event.

A metadata entry replaces the fit only when the metadata gives a `default` (see
{ref}`currmod-metadata-dir`).

(hpge-aoemeanmod-extraction)=

## A/E mean energy dependence

In the data, the A/E classifier is normalised so that its mean does not depend
on the energy. The mean of the simulated A/E does depend on it. The Simflow can
correct for this with a second pass (see {ref}`hit-tier-settings`): a first
`hit` pass is built for a calibration-source simulation (by default the Tl-208
decays in the calibration sources), the mean simulated A/E is fitted as a
function of the energy in Compton-continuum bands, as the data calibration does,
and the final `hit` pass subtracts it from the simulated A/E. The fit is done
once per detector and holds for all runs. When the same data exist in the data
production, the same fit on data is stored next to it for validation. The
settings are described in {ref}`aoemeanmod-metadata-dir`.

(hpge-psl-overview)=

## Pulse-shape simulation

The PSD simulation needs to know how the current pulse of an event depends on
where in the crystal the energy is deposited. The Simflow gets this from a
simulation of the charge collection in each detector with
[`SolidStateDetectors.jl`](https://juliaphysics.github.io/SolidStateDetectors.jl/stable/)
(SSD), at every point of a grid in the crystal, in two forms:

- a drift-time map, which the single-template PSD uses together with the
  current-pulse template;
- a pulse-shape library, which gives the full current pulse of every point after
  the readout electronics. It is built only when the PSL-based PSD is enabled
  (see {ref}`hit-tier-settings`), through four stages: the ideal pulse shapes
  are simulated, average pulses (superpulses) are built from the calibration
  data, the response of the readout electronics is fitted by comparing the two,
  and the ideal pulses are convolved with the fitted response.

Both simulations use the crystal geometry and the impurity profile from the
LEGEND metadata, and share the same settings (see {ref}`ssd-settings-meta`).
Before the field is solved, the impurity profile is scaled so that the simulated
depletion voltage reproduces the measured one, or the one of a tuned impurity
model when the production uses one; the job stops if the two still disagree. The
scale factor and the depletion voltages are recorded in the detector-info files
({ref}`par-detinfo`).

(hpge-dtmap-extraction)=

## Drift-time maps

For every modelled detector and operating voltage, SSD solves the electric field
and simulates a point-like 2039 keV deposition at every point of an $(r, z)$
grid, for two azimuths (0° and 45° to the crystal axes). The map stores, for
each point, the time at which the hole current is largest, the quantity the data
analysis measures as the time of the current maximum. The maps of all the
detectors of a run are merged into one file per run.

(hpge-ideal-psl-extraction)=

## Ideal pulse shapes

The same simulation also records the full charge signal at every grid point: the
pulse that the detector would give before the readout electronics. These ideal
pulse shapes are the input of the electronics fit and of the realistic library.

(hpge-superpulses-extraction)=

## Data superpulses

To measure the response of the readout electronics, the Simflow builds average
current and charge pulses from the calibration data of each detector, binned in
drift time. The events come from the reference calibration runs of the run list
(each counted once): single-detector events between 1500 and 2000 keV that pass
the data quality and single-site A/E selection. The drift time of an event is
the time of its current maximum minus its start time. In each drift-time slice,
up to 100 waveforms are averaged, the waveforms that differ too much from this
average are discarded with a $\chi^2$ cut on the charge pulses, and the mean of
the remaining ones is the superpulse. By default one set of superpulses per
detector is built from all runs together; they can also be built per run (see
{ref}`superpulses-settings-meta`).

(hpge-elecmod-extraction)=

## Electronics response

The readout electronics are described by the convolution of a Gaussian
(digitizer bandwidth $\sigma$) and a causal exponential decay (preamplifier
response, decay constant $\tau$),

$$
h = h_\text{digi} * h_\text{preamp}, \qquad
h_\text{digi}(t) \propto \exp\!\left[-\frac{1}{2}\left(\frac{t}{\sigma}\right)^2\right], \qquad
h_\text{preamp}(t) \propto e^{-t/\tau}\,\Theta(t),
$$

with $\sigma$ and $\tau$ in ns, $\Theta$ the Heaviside step, and each factor
normalised to unit sum. Convolving an ideal charge pulse with $h$ and
differentiating gives the simulated current pulse.

The two parameters are fitted per run and per detector: the ideal pulses whose
drift time falls in a slice are convolved, differentiated and smoothed like the
current in the data analysis, then compared with the superpulse of that slice
within ±500 ns of the current maximum. The cost is the RMS difference, weighted
by the data amplitude squared, averaged over at most five slices with drift
times between 600 and 3000 ns, and it is minimised with Minuit (the simplex
algorithm by default), with $\sigma$ between 1 and 50 ns and $\tau$ between 1
and 200 ns. Since the superpulses are the same for all runs by default, the fit
of a detector differs between runs only through the ideal pulses of the run's
voltage. The settings are in {ref}`elecmod-fit-settings-meta`. A metadata entry
replaces the fit only when the metadata gives a `default` (see
{ref}`elecmod-metadata-dir`).

(hpge-realistic-psl-extraction)=

## Realistic pulse shapes

The realistic pulse-shape library of a detector in a run is built from the ideal
pulse shapes at the operating voltage of the run and the fitted electronics
response:

1. the ideal charge pulses are convolved with the response and differentiated
   into current pulses;
2. the current pulses are smoothed with the same moving-window average as the
   current used for A/E in the data analysis (taken from the data production
   when available, see {ref}`elecmod-fit-settings-meta`);
3. the pulses are aligned at their maximum, and the drift time of each grid
   point is the time of the maximum, as in the data;
4. the amplitudes are normalised by the mode of the A/E distribution over the
   crystal.

The libraries of all the detectors of a run are merged into one file per run,
which the `hit` step reads to assign to each simulated energy deposition the
pulse of its position and compute its A/E.
