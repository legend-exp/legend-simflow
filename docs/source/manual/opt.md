# `opt` — optical hit building

The `opt` step simulates the light detected by the SiPMs. Energy deposited in a
scintillator (the liquid argon, and optionally other light sources such as PEN
plates or scintillating fibers) produces scintillation photons, some of which
reach the SiPMs and give photoelectrons. The step turns the energy depositions
of the `stp` step into these photoelectrons, with their times and amplitudes,
and applies the response of the SiPMs. Like `hit`, it works one detector table
at a time and builds no events (see {ref}`detector-tables-and-events`): its
output holds one table per SiPM, or one table for all the SiPMs together. The
settings are described in {ref}`opt-tier-settings` and the output fields in
{ref}`opt-tier`.

## From energy to photoelectrons

The number of photons emitted in each energy deposition (step) follows from the
deposited energy, the light yield of the scintillator and its quenching, which
depends on the particle type. The probability that a photon emitted at the
position of the step gives a photoelectron in a SiPM is read from an optical
map: a table of detection probabilities on a grid in the volume of the light
source, simulated beforehand by tracking optical photons through the geometry.
Maps are usually simulated for an ideal SiPM, so the probabilities are
multiplied by the efficiency of the SiPMs, one value for all channels or one per
channel. The number of photoelectrons of each step is drawn from a Poisson
distribution with this mean, and their times follow the scintillation time
profile of the material, for liquid argon its fast and slow components.

The photoelectrons can be sampled for each SiPM, with one map per channel, or
for all the SiPMs together with a single map, which is faster. A light source
without an optical map can use a light model instead: a function of the position
that gives the detection probability directly, for example a uniform one for the
PEN plates.

The light of all the sources in the same event window is pooled, so that liquid
argon and PEN light arriving together forms one signal. An event window groups
the hits of one simulated event closer than 10 µs in time, as listed by the TCM
of the `stp` file. Each event window with energy deposited in a light source
gives one row of the output.

Some simulations place hardware in the argon that the optical map does not
contain, such as a calibration source and its absorber. For them, the region
around the hardware can be replaced by a map simulated with the hardware in
place.

## SiPM response

The photoelectrons are then turned into the signals the data analysis sees:

- each photoelectron gets an amplitude drawn from a Gaussian around one, with
  the single-photoelectron resolution of the SiPMs;
- photoelectrons closer in time than the time resolution of the system are
  merged into one, at the time of the first and with the summed amplitude.

The number of photoelectrons per event window is capped, per SiPM or for all
SiPMs together, to keep the memory and the processing time bounded. Windows that
reach the cap are flagged as saturated, so that their amplitude can be taken as
a lower limit.

## Tracked photons

Instead of an optical map, the `stp` step can track the optical photons up to
the SiPMs, applying the SiPM efficiency in the geometry. The `opt` step then
takes every photon recorded in a SiPM as a photoelectron and applies only the
SiPM response above. As in `hit`, each hit of a SiPM gives one output row, so
the light of any origin is kept, not only that of the scintillators. The output
is always per SiPM. This is slower to simulate but needs no optical map.

## Run partitioning

The events of each partition (see {ref}`run-partitioning`) record the period and
the run they belong to, and the usability of each SiPM in that run. A SiPM that
is off in a run is kept in the `opt` output; the `evt` step leaves it out when
it builds the events of that run.
