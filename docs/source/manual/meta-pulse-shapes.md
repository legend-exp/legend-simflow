# Pulse-shape simulation settings

These settings control the simulation of the HPGe pulse shapes in the `par` step
(see {ref}`hpge-psl-overview`) and the impurity tuning (see
[](impurity-tuning.md)). Each is a settings file in
`simprod/config/pars/{experiment}/geds/` (see [](meta.md)). The defaults suit
most productions.

(ssd-settings-meta)=

## Field simulation

`ssd/settings.yaml` sets the grid of the
[SolidStateDetectors.jl](https://juliaphysics.github.io/SolidStateDetectors.jl/stable/)
simulation, shared by the drift-time maps, the pulse-shape libraries and the
impurity scan:

`grid_size_in_mm` (default 0.5)

: the spacing of the $(r, z)$ grid in mm. The computing time grows with the
square of `1/grid_size_in_mm`; a coarse grid such as 10 mm is enough for tests.

`ssd_refinement_limits` (default `[0.2, 0.1, 0.05, 0.02]`)

: the thresholds of the successive refinements of the field grid. Too coarse a
refinement can leave a detector undepleted.

`padding` (default 3)

: grid points added around the crystal.

(superpulses-settings-meta)=

## Superpulses

`superpulses/settings.yaml` controls the average pulses built from the
calibration data (see {ref}`hpge-superpulses-extraction`). They are needed for
the electronics fit of the pulse-shape-library PSD and for the impurity tuning.

`build_per_runid` (default `false`)

: build the superpulses of each run separately, rather than one set per detector
from all runs together. Not allowed in the impurity tuning.

`drift_time_slices` (default `"1000:200:2000"`)

: the drift-time slices, as `start:step:stop` in integer ns.

`tier_name` (default `auto`)

: the data tier the events are selected from: `evt`, `pet`, `hit` or `pht`;
`auto` takes `pet` if the data production has it, else `evt`.

`t0_field`, `end_time_field` (default `auto`, `geds/psd/low_aoe/time`)

: the fields of the start and end time of an event, whose difference is its
drift time. `auto` takes the SiPM event time if present, else the HPGe start
time; it works only with the `evt` and `pet` tiers.

Advanced: `target_wfs` (100) and `min_number_wfs` (10), the waveforms averaged
per slice and the fewest for a slice to be kept; `chi2_threshold` (3), the
outlier cut; `max_files` (all), a limit on the data files per run for tests;
`charge_output`, `curr_output`, `energy_output`, the names of the processed
waveforms and energy in the data production.

(elecmod-fit-settings-meta)=

## Electronics fit

`elecmod/settings.yaml` controls the fit of the electronics response to the
superpulses (see {ref}`hpge-elecmod-extraction`), both in the
pulse-shape-library PSD and in the impurity scan. When the
{ref}`elecmod-metadata-dir` directory gives a `default`, the PSD fit does not
run.

`curr_mwa` (default: from the data production)

: the moving-window average applied to the simulated currents, as
`{length, num_mw, mw_type}`, with `length` in ns. It is taken from the
processing of the current used for A/E in the data production when that is
available, and the setting is then ignored. Without the data production, it is
taken from this setting, or else from built-in values (48 ns, 3 windows).

`sigma_start`, `tau_start` (default 10 ns, 50 ns)

: the start values of the fit.

`sigma_limits`, `tau_limits` (default `[1, 50]` ns, `[1, 200]` ns)

: the allowed ranges of the two parameters.

`dt_range_tuning`, `max_num_superpulses` (default `[600, 3000]` ns, 5)

: the fit uses the slices with drift times in this range, at most this many,
those with the longest drift times.

`comparison_window` (default `[-500, 500]` ns)

: the time window, around the current maximum, in which simulated and measured
pulses are compared.

Advanced: `weight_power` (2), the power of the data amplitude that weights the
residuals; `minimiser_mode` (`simplex`; also `migrad` or `both`) and `max_calls`
(1000), the minimiser; `angle` (`"000"`), the crystal axis of the ideal pulses
used; `plot_window`, the range of the control plots. The impurity scan uses its
own defaults for `tau_start` (30 ns), `tau_limits` (`[2, 100]` ns) and
`max_calls` (100), and starts each fit from the result of the previous grid
point.

(psl-scan-settings-meta)=

## Impurity scan grid

`ssd/scan_settings.yaml` sets the grid of impurity profiles simulated by the
impurity tuning (see [](impurity-tuning.md)). The file must exist in a tuning
production; missing keys keep their defaults. Each key is a range written
`"start:step:stop"`, quoted.

```{code-block} yaml
:caption: pars/{experiment}/geds/ssd/scan_settings.yaml

depv_shift: "-1000:20:-20"
slope: "-0.9:0.2:3.0"
```

`depv_shift`

: the depletion voltages, in V, relative to the operational voltage of the
detector.

`slope`

: the slope factor $s$: the variation of the impurity density along the crystal
is multiplied by $1 + s$, so 0 keeps the profile of the metadata and $-1$ makes
it flat.

The defaults give $50 \times 20 = 1000$ simulations per detector.

(impurityscan-settings-meta)=

## Impurity fit

`impurityscan/settings.yaml` controls the comparison with the data in the
impurity tuning (see [](impurity-tuning.md)). The file must exist in a tuning
production.

```{code-block} yaml
:caption: pars/{experiment}/geds/impurityscan/settings.yaml (the defaults)

ssd:
  grid_size_in_mm: null
drift_time_scan:
  min_energy_in_keV: 1500
  max_events: null
impurity_fit:
  energy_range_in_keV: [1500, 2500]
  drift_time_uncertainty_in_ns: 50
  waveform_uncertainty: 0.5
  drift_time_observables:
    percentile: 90
    smoothing_in_ns: 50
    min_peak_fraction: 0.25
```

`ssd.grid_size_in_mm`

: the grid spacing of the scan; `null` uses the one of the field simulation.

`drift_time_scan.min_energy_in_keV`, `drift_time_scan.max_events`

: the simulated events whose drift times are computed: above this energy in the
active volume, at most this many per detector and simulation (`null`: all). The
energy must not be above the lower edge of `impurity_fit.energy_range_in_keV`.

`impurity_fit.energy_range_in_keV`

: the energy window of the events compared between data and simulation.

`impurity_fit.drift_time_uncertainty_in_ns`, `impurity_fit.waveform_uncertainty`

: the systematic uncertainties that scale the drift-time and the pulse-shape
terms of the $\chi^2$.

`impurity_fit.drift_time_observables`

: how the two drift-time observables are taken: the `percentile` of the
distribution, and the first peak of the distribution smoothed with a Gaussian of
`smoothing_in_ns`, the first peak higher than `min_peak_fraction` of the maximum
(see {func}`legendsimflow.drift_time.drift_time_observables`).
