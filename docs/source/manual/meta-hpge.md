# HPGe parameters

The `par` step collects the parameters of the HPGe detectors per run and per
detector, mostly from the LEGEND-200 data production (see [](par.md)). The
metadata in `simprod/config/pars/{experiment}/geds/` can replace them, and
selects which detectors get a pulse-shape model.

## Parameter directories

Each parameter has a validity-based directory (see {ref}`validity-dirs-meta`)
whose data files map detector names to entries:

```{code-block} yaml
:caption: pars/{experiment}/geds/eresmod/l200-p03-r%-T%-all-eresmod.yaml

default:
  expression: FWHMLinear
  parameters: {a: 0.5, b: 0.001}
V02160A:
  expression: FWHMLinear
  parameters: {a: 0.3, b: 0.0009}
```

For each run, the directory decides where the parameters come from:

- no directory, or no file valid for the run: every detector with an entry in
  the data production gets its measured parameters;
- a file without `default`: the measured parameters are collected as above, and
  the listed detectors get their entry instead;
- a file with `default`: the data production is not used: every HPGe detector of
  the channel map gets its own entry or the `default`.

An entry replaces the measured one whole. Without the data production, each
directory needs a `default`. The current-pulse model and the electronics
response are fitted by the Simflow rather than collected, so their directories
are used only when they give a `default`.

(eresmod-metadata-dir)=

### Energy resolution

`eresmod/`. An entry has `expression: FWHMLinear` and `parameters: {a, b}`, for

$$
\text{FWHM}(E) = \sqrt{a + b E},
$$

with $E$ and the FWHM in keV ($a$ in keV², $b$ in keV). No other expression is
supported here. `uncertainties` may be given and are not used.

(aoeresmod-metadata-dir)=

### A/E resolution

`aoeresmod/`. An entry has `expression: SigmaFit` and `parameters: {a, b, c}`,
for the width of the A/E distribution of single-site events,

$$
\sigma_{A/E}(E) = \sqrt{a + (b/E)^c},
$$

with $E$ in keV and $\sigma_{A/E}$ in units of A/E.

(psdcuts-metadata-dir)=

### PSD cuts

`psdcuts/`. An entry has `aoe: {low_side, high_side}`, the low and the high cut
on the A/E classifier, in units of the A/E resolution. Missing values are taken
from `psdcuts_default` of the `hit` settings.

(currmod-metadata-dir)=

### Current-pulse model

`currmod/`. Used only for the single-template PSD, for the detectors that get a
pulse-shape model, and only when the file gives a `default`: otherwise the model
is fitted on the data. An entry has:

- `current_pulse_pars`: the parameters of the template, `amax`, `mu`, `sigma`,
  `tail_fraction`, `tau`, and optionally `high_tail_fraction` and `high_tau`
  (default 0); times in ns, fractions dimensionless. The function is
  `reboost.hpge.psd._current_pulse_model`;
- `mean_aoe`: the A/E that the template gives at 1593 keV;
- `current_reso`: the electronic noise on the maximum of the current. The A/E is
  smeared with `current_reso / mean_aoe`.

(elecmod-metadata-dir)=

### Electronics response

`elecmod/`. Used only for the pulse-shape-library PSD of a standard production,
and only when the file gives a `default`: otherwise the response is fitted to
the superpulses of the data. An entry has `sigma` and `tau` in ns, the digitizer
bandwidth and the preamplifier decay constant (see
{ref}`hpge-elecmod-extraction`). With a `default`, no superpulses are built for
the run. The impurity tuning always fits the response and does not read this
directory.

(aoemeanmod-metadata-dir)=

## A/E mean energy dependence

There is no directory for this parameter. Without the two-pass correction (see
{ref}`hpge-aoemeanmod-extraction`), every detector uses `aoemeanmod_default` of
the `hit` settings. With it, the correction is fitted per detector on the
simulations whose simid matches `simid_regex` in
`pars/{experiment}/geds/aoemeanmod/settings.yaml` (default
`sis*_z*_slot*_Pb212_to_Pb208`), and the simulation list of the production must
contain one.

(modeling-settings-meta)=

## Modelled detectors

`modeling/settings.yaml` sets `min_voltage_above_depletion_in_V` (default 100):
how far above its depletion voltage a detector must be operated to get a
pulse-shape model (see {ref}`hpge-modeling-criteria`). It is not used when the
production uses tuned impurity profiles, nor in the impurity tuning itself.

(skip-metadata-dir)=

## Skip list

`skip/` is a validity-based directory that excludes detectors from the
pulse-shape model for given runs, whatever their status in the channel map. Its
files map detector names to a reason:

```{code-block} yaml
:caption: pars/{experiment}/geds/skip/l200-p03-r005-T%-all-skip.yaml

V02160A: "simulation does not converge"
```

The reason is documentation only; the workflow does not read it. A skipped
detector gets no PSD observables and never stops the `hit` step (see
{ref}`build-tier-hit-hpge`). The skip list also applies to the impurity tuning.
After editing it, rebuild the cache of the modelled detectors (see [](meta.md)).
