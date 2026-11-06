# `stp` — particle simulation

The `stp` step simulates the passage of particles through the experiment with
[_remage_](https://remage.readthedocs.io), a simulation program based on Geant4.
It records the energy depositions (steps) in the detectors.

Each simulation has its own geometry. The geometry generator of the experiment
(for example [legend-pygeom-l200](https://legend-pygeom-l200.readthedocs.io) for
LEGEND-200) builds it from the LEGEND metadata, with the changes requested by
that simulation (see {ref}`geom-config.yaml`). It writes the geometry to a GDML
file (a text format for Geant4 geometries). The `vtx`, `opt` and `hit` steps
read the same GDML file, so the whole chain uses the same geometry.

The Simflow then writes one _remage_ macro per simulation, which sets the
primary particles and where they start (see {ref}`simconfig.yaml`). All jobs of
a simulation run this macro. Each job gets its number of primaries and a random
seed for _remage_. The seed is a hash of the path of the output file of the job,
so every job has a different seed, and a job that runs again with the same
inputs simulates the same events.

Each job runs _remage_ in a single thread: the Simflow gets its parallelism from
running many jobs at the same time. A job writes one `stp` file, with one table
per detector and a TCM (see {ref}`detector-tables-and-events`).
