# `hit` — hit tier building

The `hit` step simulates the response of the HPGe detectors: it turns the energy
depositions of the `stp` step into the energy and the pulse-shape observables
that the data analysis measures. Like `opt`, it works one detector table at a
time and builds no events (see {ref}`detector-tables-and-events`). It processes
the events of each run partition with the parameters of that run (see
{ref}`run-partitioning`), which the `par` step provides (see [](par.md)). The
settings are described in {ref}`hit-tier-settings` and the output fields in
{ref}`hit-tier-output`.

## Active volume

Near the n+ surface of a detector the charge is not fully collected. For every
energy deposition (step), the step computes the distance to the n+ surface and
weights the deposited energy with the charge collection efficiency at that
depth: zero in the dead layer, rising linearly to one at the full charge
collection depth (FCCD) of the detector, taken from the detector metadata. The
dead layer is a fixed fraction of the FCCD, set in the settings. The sum of the
weighted energies of a hit is its true energy.

Optionally, the energies deposited within a few millimetres of the n+ surface
are also stored step by step, with their distance to the surface, so that the
energy can be recomputed later for another model of the surface layer.

## Energy resolution

The true energy is smeared with a Gaussian of standard deviation FWHM(E)/2.355,
with the energy resolution of the detector in the run.

## Pulse-shape discrimination

The A/E of a simulated event, the maximum A of the current pulse divided by the
energy E, is computed by building the current pulse of the event as the sum of
the pulses of its energy depositions, each weighted by its active energy. Two
methods give the pulse of a single deposition (see {ref}`hpge-psl-overview`):

- **single template** (the default): the current-pulse template of the detector,
  shifted in time by the drift time of the deposition. The drift time is read
  from the drift-time map at the position of the deposition. The map is
  simulated along two crystal axes, 45° apart, and the drift time in between is
  interpolated with the four-fold symmetry of the germanium crystal;
- **pulse-shape library**: the realistic current pulse of the deposition, taken
  from the library at its position, which includes the drift time.

The two methods can run together, each with its own output fields. In both, the
maximum of the current is smeared with the electronic noise measured for the
current-pulse template, and then:

- the A/E is corrected for its energy dependence: its mean at the energy of the
  event is subtracted and one is added, so that the corrected A/E is centred at
  one (see {ref}`hpge-aoemeanmod-extraction`);
- the A/E classifier is the distance of the corrected A/E from one, in units of
  the A/E resolution of the detector in the run, as in the data analysis;
- an event is single-site when the classifier is above the low PSD cut. The
  pulse-shape-library method also applies the high cut, which rejects events
  near the p+ contact, and flags the events that pass both cuts.

The drift time at the maximum of the current is stored too. A detector without a
drift-time map or a current-pulse template, for example one that cannot be
modelled (see {ref}`hpge-modeling-criteria`), gets no PSD observables.

## Bookkeeping

Each output row records the period and the run of its partition, the usability
of the detector and of its PSD in that run, and whether the crystal metadata of
the detector is complete enough for the simulation. The detectors that are off
in a run are processed as well; the `evt` step uses the usability to treat them.
The step also builds a TCM of its output.

(build-tier-hit-hpge)=

## Hard errors and fallbacks

The step checks that every simulated detector has the parameters it needs from
the `par` step. Depending on the detector, a missing parameter is either a hard
error, which stops the step, or replaced by a fallback from
{ref}`hit-tier-settings`, with a warning:

| Parameter             | Hard error                                                                                                    | Fallback (+ warning) | Fallback key         |
| --------------------- | ------------------------------------------------------------------------------------------------------------- | -------------------- | -------------------- |
| Energy resolution     | ON detector missing entry                                                                                     | all other detectors  | `eresmod_default`    |
| A/E resolution        | ON detector with `psd_usability ≠ "missing"` and PSD modelable (drift-time map + current model) missing entry | all other cases      | `aoeresmod_default`  |
| PSD cuts              | ON detector with `psd_usability ≠ "missing"` and PSD modelable (drift-time map + current model) missing entry | all other cases      | `psdcuts_default`    |
| A/E energy dependence | never                                                                                                         | all detectors        | `aoemeanmod_default` |

The energy resolution is checked for every detector, whether PSD is simulated or
not. The A/E resolution and the PSD cuts can only stop the step for a detector
whose PSD is modelled, so without PSD simulation they always fall back. A
fallback is never needed for a parameter whose metadata gives a `default` (see
{ref}`par-collection-model`), because every detector then has an entry.

:::{note}

An ON detector with `psd_usability = "missing"` explicitly signals that PSD data
are unavailable for that detector (e.g. a known hardware issue). Similarly, a
detector missing a drift-time map or current-pulse model cannot have its PSD
response simulated at all: its PSD output fields are NaN, with a warning when
the detector is ON with valid PSD in the data. Both situations fall back to the
default A/E resolution and PSD cuts rather than raising a hard error. A detector
listed in the {ref}`skip-metadata-dir` directory is treated in the same way: it
has no drift-time map or current-pulse model, so it gets NaN PSD fields and the
fallback A/E resolution and PSD cuts, with no hard error.

:::
