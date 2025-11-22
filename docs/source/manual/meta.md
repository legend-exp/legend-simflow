# Data names and paths

## File naming conventions

The file naming convention for the `stp` and `hit` tiers is:

```
{experiment}-{simid}-job_{jobid}-tier_{tier}.{extension}
```

where each label (in curly brackets `{}`) is alphanumeric (including
underscores: `_`). Do not use dashes (`-`) or other characters.

- `experiment` — name or label for the experimental configuration being
  simulated.

- `simid` — stands for "simulation identifier", i.e. a string to uniquely label
  a simulation. Typically specifies the physical process being simulated and the
  experiment's components involved.

- `jobid` — stands for "(simulation) job identifier". It's a zero-padded integer
  that labels independent jobs across which the simulation is split.

- `tier` — the three-character label of the tier. At the moment the simflow
  supports `vtx`, `stp` and `hit` tiers.

- `extension` — file extension. `lh5` for LEGEND HDF5 files, `gdml` for GDML
  geometry files, `yaml` for plain-text YAML configuration files, `log` for log
  files.

# Metadata

The workflow configuration metadata is stored in
[legend-simflow-config](https://github.com/legend-exp/legend-simflow-config).
Documentation about which metadata is stored there (e.g. which LEGEND
experimental configuration are supported) can be found in the `README.md`.

In this section, the specification of the metadata format is documented.

## `tier` — static tier configuration

Metadata is organized in this directory by tier (first level) and experimental
configuration (second level).

### `stp` tier

The simulation jobs are configured centrally with a `simconfig.yaml` file.
