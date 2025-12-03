# legend-simflow

End-to-end Snakemake workflow to run Monte Carlo simulations of signal and
background signatures in the LEGEND experiment and produce probability-density
functions (pdfs). Configuration metadata (e.g. rules for generating simulation
macros or post-processing settings) is stored at
[legend-simflow-config](https://github.com/legend-exp/legend-simflow-config).

## Key concepts

- Simulations are labeled by an unique identifier (e.g. `hpge-bulk-2vbb`), often
  referred as `simid` (simID). The identifiers are defined in
  [legend-simflow-config](https://github.com/legend-exp/legend-simflow-config)
  through `simconfig.yaml` files in tier directories `stp` and `vtx`.
- A simulation (`simid`) can consist of several jobs (simulation jobs with
  different random seeds run in parallel). Each job is assigned its own `jobid`
  (integer number).
- Each simulation is defined by a template macro (also stored as metadata) and
  by a set of rules (in `simconfig.yaml`) needed to generate the actual macros
  (template variable substitutions, number of primaries, number of jobs, etc).
- The production is organized in tiers. The state of a simulation in a certain
  tier is labeled as `<tier>.<simid>`. Snakemake understands this syntax.
- The generated pdfs refer to a user-defined selection of LEGEND data taking
  runs. Such a list of runs is specified through the configuration file.
- The production can be restricted to a subset of simulations by passing a list
  of identifiers to Snakemake.

### Workflow steps (tiers)

1. Tier `vtx` building: run simulations that generate Monte Carlo event vertices
   needed to some simulations in the next tier. Simulations that do not need a
   special event vertices will directly start from tier `raw`.

1. Tier `stp` building: run full event simulations. Simulation macro commands
   are generated according to rules defined in the metadata.

1. Run the first (hit-oriented) step of simulation post-processing. Here, a
   "hit" represents a collection of Geant4 "step" in a single detector
   (sensitive volume). In these tiers, hit-wise operations like optical map or
   HPGe detector models application are typically performed. The "run
   partitioning" is also performed at this stage (see below). Two tiers belong
   to this group:
   - Tier `opt`: convolution of the optical models (i.e. optical "maps")
   - Tier `hit`: convolution of HPGe detector models (energy and pulse shape)

1. Tier `evt` building: multiple operations are performed in order to build
   actual events.

1. Tier `pdf` building: summarize `evt`-tier output into histograms (the pdfs).

### Run partitioning

"Run partitioning" refers to incorporating information about the experiment's
data taking runs for which the user wants to build pdfs:

- Partition the simulated event statistics into fractions corresponding to the
  actual total livetime fraction spanned by each selected run. This information
  is extracted from
  [`legend-metadata/datasets/runinfo.yaml`](https://github.com/legend-exp/legend-datasets/blob/main/runinfo.yaml)
- For each partition, apply HPGe models such as energy resolution or
  pulse-shape.
- ...apply optical models (detection probability lookup tables) for the
  scintillators.
- ...apply detector status flags (available in
  [`legend-metadata/datasets/statuses`](https://github.com/legend-exp/legend-datasets/blob/main/statuses))

## Next steps

```{toctree}
:maxdepth: 1

Package API reference <api/modules>
Snakemake rules <api/snakemake_rules>
```

```{toctree}
:maxdepth: 1

manual/index
```

```{toctree}
:maxdepth: 1
:caption: Related projects

remage <https://remage.readthedocs.io/>
```

```{toctree}
:maxdepth: 1
:caption: Development

Source Code <https://github.com/legend-exp/legend-simflow>
```
