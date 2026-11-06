(impurity-tuning)=

# HPGe impurity tuning

The electric field in an HPGe crystal, and with it the drift times and the pulse
shapes of the detector, depend on the net impurity density in the crystal. The
impurity profile along the crystal is known only approximately: the crystal
metadata gives a model fitted to a few measurements on the boule by the crystal
vendor, and the detector metadata gives the measured depletion voltage.

The impurity tuning determines the best-fit impurity profile of each detector by
matching two things to LEGEND-200 data: the shape of the current pulses, and the
distribution of the drift times. Both are simulated for a grid of impurity
profiles, and the profile that reproduces the data best is kept. The profiles of
the grid differ in two parameters:

- the slope, a factor that scales the variation of the impurity density along
  the crystal (0 keeps the profile of the metadata; the density at the seed end
  of the boule is kept fixed);
- the depletion voltage, reached by rescaling the whole impurity density.

```{image} ../img/impurity-model-light.svg
:alt: Impurity density along the boule of crystal V04545 and the effect of the two tuning parameters
:class: only-light
```

```{image} ../img/impurity-model-dark.svg
:alt: Impurity density along the boule of crystal V04545 and the effect of the two tuning parameters
:class: only-dark
```

The figure shows the impurity model of one crystal, fitted to the measurements
of the vendor, and how the two parameters change it.

The tuned values are then used by the pulse-shape simulation of later
productions (see {ref}`hpge-psl-overview`).

## A tuning-only production

:::{warning}

The tuning and a standard production cannot run in the same production. With the
tuning switched on, the file of tuned profiles is the only output: no `hit`,
`evt` or `pdf` files are built, and the detectors are selected with the tuning
criteria below. Run the tuning first, then a separate standard production that
reads the tuned profiles.

:::

The tuning runs as a separate production, switched on in the `hit` settings (see
{ref}`hit-tier-settings` and {ref}`impurityscan-settings-meta`). Each simulation
of this production is compared with the data of the single run in its run list.
The production needs the LEGEND-200 data production, and the pulse-shape
simulation of the `hit` step to be enabled. Each detector must be operated at
the same voltage in all the runs, because the superpulses and the scan below
combine the runs.

A detector is tuned when it is on in the run, is not excluded in the skip list,
has an operational voltage and an impurity profile in its crystal metadata, and
its pulse-shape discrimination is marked as valid in the data. Unlike in a
standard production (see {ref}`hpge-modeling-criteria`), the margin above the
depletion voltage is not checked, since the depletion voltage is one of the
fitted parameters.

Each tuned detector goes through the chain in the figure below; the following
sections describe its steps.

```{raw} html
:file: ../img/imptune.svg
```

## Pulse-shape scan

For each detector, the ideal pulse-shape library (see
{ref}`hpge-ideal-psl-extraction`) is simulated at the operating voltage for
every point of a grid in the two parameters. By default the slope goes from −0.9
to 3.0 in steps of 0.2, and the depletion voltage from 1000 V to 20 V below the
operating voltage in steps of 20 V: 1000 simulations per detector (see
{ref}`psl-scan-settings-meta`).

## Electronics response at every grid point

The data superpulses are built from the calibration data as in a standard
production (see {ref}`hpge-superpulses-extraction`). At every grid point, the
response of the readout electronics is fitted to them (see
{ref}`hpge-elecmod-extraction`). The RMS difference left by the fit measures how
well the pulse shapes of that grid point reproduce the measured ones. The
response of the grid point with the smallest RMS is used to compute the
simulated drift times.

## Drift times

For every event of the tuning simulations that deposits more than 1500 keV in
the active volume of the detector, the drift time is computed at every grid
point, from the pulse-shape library of that point convolved with the best-fit
electronics response.

In the data, the events of the run between 1500 and 2500 keV are selected:
single-site-like events that are not forced triggers, not in coincidence with a
muon, and with more than 10 photoelectrons in the SiPMs. The drift time of an
event is the time of the current maximum in the HPGe detector minus the time of
the event given by the argon scintillation light in the SiPMs.

The two distributions are compared through two observables: the position of the
first peak of the smoothed distribution, and its 90th percentile.

## Fit

At every grid point, the cost is the sum of a drift-time term and a pulse-shape
term,

$$
\chi^2 = \frac{(\Delta t_\text{peak})^2 + (\Delta t_{90})^2}{(50\,\text{ns})^2}
+ \frac{\text{RMS}^2}{0.5^2},
$$

where $\Delta t_\text{peak}$ and $\Delta t_{90}$ are the differences between
data and simulation in the two drift-time observables, and RMS is the residual
of the electronics fit. The two denominators are the systematic uncertainties of
the comparison: a closer match than about 50 ns in the drift times, or than an
RMS of 0.5 in the pulse shapes, is not meaningful given the limits of the
simulation and of the data. They set the scale of each term, so that the two
terms weigh in the sum according to how well each can be trusted. The grid point
with the smallest $\chi^2$ gives the tuned slope and depletion voltage of the
detector; there is no interpolation between grid points. The uncertainties can
be changed (see {ref}`impurityscan-settings-meta`).

The tuned values of all the detectors are written to one file. A later
production reads it (see the `hpge_impurity_model` path in [](setup.md)): its
drift-time maps and pulse-shape libraries use the tuned profiles, and its
detectors are selected by their presence in the file (see
{ref}`hpge-modeling-criteria`).

:::{note}

The data files of the run are found when the tuning runs and are not tracked as
inputs of the workflow. After an update of the data production, rerun the fit
explicitly.

:::
