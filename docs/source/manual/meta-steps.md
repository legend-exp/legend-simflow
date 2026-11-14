# Step settings

Each processing step reads a settings file,
`simprod/config/tier/{step}/{experiment}/settings.yaml`, that applies to all the
simulations of an experimental configuration. What the steps do is described in
[](pipeline.md); this page lists what each setting changes. Keys marked required
have no default.

All step settings have a required `buffer_len`: how much data is read at once,
e.g. `"50*MB"`. It sets the memory used by the step and does not change the
output.

(opt-tier-settings)=

## `opt`

See [](opt.md) for the step and {ref}`opt-tier` for its output.

```{code-block} yaml
:caption: tier/opt/{experiment}/settings.yaml

scintillator_volume_name: liquid_argon
optmap_per_sipm: true
optmap_scaling_factor: 0.3
photoelectron_resolution_sigma: 0.3
time_resolution_in_ns: 16
max_pes_per_hit_per_sipm: 5
max_pes_per_hit_combined: 100
buffer_len: "10*MB"
```

`light_source` (default `optmap`)

: where the photoelectrons come from: `optmap`, sampled from the optical maps or
light models of the light sources, or `tracked_photons`, the photons tracked up
to the SiPMs in the `stp` step (see [](opt.md)). With tracked photons, the SiPM
efficiency must be applied in the geometry, and the settings that act on the
maps below are ignored.

`scintillator_volume_name` (required without `scintillators`)

: the name of the liquid argon table in the `stp` files, e.g. `liquid_argon`. It
defines the only light source, `lar`.

`scintillators` (optional)

: the light sources, as `<name>: {table, material, ...}`: `table` is the name of
the table of the source in the `stp` files, `material` its scintillation model
(`lar`, `pen` or `fiber`). A source can set its own `optmap_scaling_factor`, and
a `model` that replaces its optical map with a light model (a Python function of
the position, see {mod}`legendsimflow.opt_models` and
{func}`legendsimflow.metadata.get_opt_scintillators`). The optical map of a
source `<name>` is read from `paths.optical_maps.<name>` of the configuration
file.

`optmap_per_sipm` (required)

: sample the photoelectrons of each SiPM with its own map (`true`) or of all
SiPMs together with one map (`false`).

`optmap_scaling_factor` (required with optical maps, unless every source sets
its own)

: the SiPM efficiency, multiplied with the detection probabilities of the maps,
which are usually simulated for an ideal SiPM. A number applies to all channels;
a mapping `<sipm>: value` sets one per channel, with a `default` for the others.

`photoelectron_resolution_sigma` (required)

: the relative resolution of a single-photoelectron amplitude.

`time_resolution_in_ns` (required)

: photoelectrons closer than this are merged into one.

`max_pes_per_hit_per_sipm`, `max_pes_per_hit_combined` (required)

: the largest number of photoelectrons per event window, for each SiPM (per-SiPM
maps or tracked photons) or for all of them (one map). Windows that reach it are
flagged as saturated; zero or less removes the limit.

`store_expected_pes` (default `false`)

: also store the expected number of photoelectrons of each window, before the
Poisson sampling.

(hit-tier-settings)=

## `hit`

See [](hit.md) for the step and {ref}`hit-tier-output` for its output.

```{code-block} yaml
:caption: tier/hit/{experiment}/settings.yaml

dead_layer_fraction: 0.5
buffer_len: "500*MB"

eresmod_default:
  expression: FWHMLinear
  parameters: {a: 0.5, b: 0.001}
aoeresmod_default:
  expression: SigmaFit
  parameters: {a: 0.0001, b: 0, c: 1}
psdcuts_default:
  aoe: {low_side: -1.5, high_side: 3}
aoemeanmod_default:
  single_template: {expression: a * x + b, pars: {a: 0, b: 1}}
  psl: {expression: a * x + b, pars: {a: 0, b: 1}}
```

`dead_layer_fraction` (required)

: the fraction of the full charge collection depth (FCCD) that is fully dead;
the charge collection rises linearly from there to the FCCD (see [](hit.md)).

`store_surface_steps` (default `false`)

: also store the energy depositions near the n+ surface, step by step, to
recompute the energy later with another surface model.

`surface_depth_in_mm` (default 3)

: the depth of the surface layer stored with `store_surface_steps`. It must be
larger than the largest FCCD that will be studied.

`simulate_psd` (default `true`)

: simulate the A/E with the single-template method.

`simulate_psd_with_psl` (default `false`)

: simulate the A/E with the pulse-shape library.

`two_pass_aoe_correction` (default `false`)

: correct the energy dependence of the simulated A/E with a first `hit` pass
(see {ref}`hpge-aoemeanmod-extraction`).

`tune_hpge_impurities_on_data` (default `false`)

: make the production a tuning-only production of the HPGe impurity profiles
(see [](impurity-tuning.md)). It needs `simulate_psd` and the LEGEND-200 data
production.

`eresmod_default`, `aoeresmod_default`, `psdcuts_default`, `aoemeanmod_default`
(required)

: the fallback energy resolution, A/E resolution, PSD cuts and A/E mean energy
dependence, used for the detectors that have none (see
{ref}`build-tier-hit-hpge`). The formats are those of the HPGe parameters (see
[](meta-hpge.md)); `aoemeanmod_default` has one function of the energy `x` in
keV for each PSD method.

(evt-tier-settings-meta)=

## `evt`

See [](evt.md) for the step and {ref}`evt-tier` for its output.

```{code-block} yaml
:caption: tier/evt/{experiment}/settings.yaml

geds_energy_thr_kev: 25
spms_energy_thr_pe: 0
lar_veto_multiplicity_thr: 4
lar_veto_energy_sum_pe_thr: 4
buffer_len: "50*MB"
```

`geds_energy_thr_kev` (required, unless `skip_hit`)

: the HPGe energy threshold in keV; hits at or below it are left out.

`spms_energy_thr_pe` (required)

: the SiPM threshold in photoelectrons; pulses at or below it are left out.

`lar_veto_multiplicity_thr`, `lar_veto_energy_sum_pe_thr` (required)

: the event is coincident with the liquid argon when the number of SiPM channels
with light, or the sum of their amplitudes in photoelectrons, reaches these
values.

`skip_opt`, `skip_hit` (default `false`)

: build events without SiPM or without HPGe data; the skipped step does not run.
They cannot both be `true`.

`add_random_coincidences` (default `false`)

: add random-coincidence SiPM pulses from the LEGEND-200 data (see [](evt.md)).
It needs the LEGEND-200 data production.

`random_coincidence_mode` (default `forced_trigger`)

: where the random coincidences come from: `forced_trigger`, the forced-trigger,
pulser and pre-trigger data of the physics runs, or `noise_trigger`, a dedicated
noise-trigger stream.

`random_coincidence_runid` (default: the run of the events)

: the run to take the random coincidences from, when it is not the run being
simulated: one run for all, or a mapping from each simulated run to its source
run, with an entry for every run.

`max_events_per_chunk` (default 200000)

: the largest number of events processed at once; it limits the memory and does
not change the output.

(cvt-tier-settings)=

## `cvt`

`buffer_len` is the only setting (see [](evt.md) for the step).

(pdf-tier-settings)=

## `pdf`

See [](pdf.md) for the step and {ref}`pdf-tier` for its output.

```{code-block} yaml
:caption: tier/pdf/{experiment}/settings.yaml

buffer_len: "50*MB"
detector_groups:
  bege: "B.*"
  icpc: "V.*"
```

`detector_groups` (optional)

: the detector groups that get their own histograms, as `<name>: <regex>`,
matched against the whole detector name. The group `all` is always added.

`fccd_tag_energy_thr_kev` (optional)

: also split the events by whether they deposit more than this energy, in keV,
in the FCCD layer (see [](pdf.md)). It needs `store_surface_steps` in the `hit`
settings.

`build_l200data_pdfs` (default `true`)

: also histogram the LEGEND-200 data. Without the data production, a warning is
issued and nothing is built.

`morph_dead_layer` (default `false`)

: also build the multiplicity-1 PDFs on a grid of dead-layer models (see
[](pdf.md)). It needs `store_surface_steps` in the `hit` settings.

`morph_grid` (default `fccd_sigma: [-5, 5, 0.5]`, `dlf: [0, 1, 0.1]`)

: the `[first, last, step]` values of the two axes of the grid. `fccd_sigma` is
the shift of the FCCD of each detector, in standard deviations of its FCCD
uncertainty. `dlf` is the dead layer fraction.

`morph_min_fccd_in_mm` (default `0.5`)

: the smallest FCCD of the grid, in mm. The shifted FCCD is limited to the range
from this value to `surface_depth_in_mm` of the `hit` settings.
