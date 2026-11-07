# Metadata

The metadata steers a production: which simulations to run, how each step
processes them, and which detector parameters to use. It lives in the
[legend-simflow-config](https://github.com/legend-exp/legend-simflow-config)
repository, inside the LEGEND metadata (see [](setup.md)); its `README.md` lists
the experimental configurations it supports. This page describes how the
metadata is organized, and the pages below document what each file sets.

## Layout

All the files are below `simprod/config/`, with one folder per experimental
configuration (`{experiment}`, e.g. `l200cfg01`):

```
simprod/config/
├── simlists.yaml                     # named lists of simulations
├── tier/
│   ├── vtx/{experiment}/simconfig.yaml   # vertex generators
│   ├── stp/{experiment}/
│   │   ├── simconfig.yaml            # the simulations
│   │   ├── generators.yaml           # primary particles
│   │   ├── confinement.yaml          # where primaries start
│   │   └── *.mac                     # macro templates
│   ├── hit/{experiment}/simconfig.yaml   # run lists per simulation
│   └── {opt,hit,evt,cvt,pdf}/{experiment}/settings.yaml   # step settings
├── geom/{experiment}-geom-config.yaml    # geometry
└── pars/{experiment}/geds/           # HPGe parameters
    ├── {eresmod,aoeresmod,psdcuts,currmod,elecmod,skip}/  # validity-based
    └── {ssd,superpulses,elecmod,modeling,impurityscan,aoemeanmod}/settings.yaml
```

## Settings files

A `settings.yaml` file sets the options of a step or of a parameter extraction.
A missing key keeps its default, and a missing file is the same as an empty one.
Nested blocks are merged key by key; a list replaces the default whole. Keys
without a default are required, and the workflow stops with an error when it
needs one that is missing. The defaults are listed in
{mod}`legendsimflow.defaults`.

(validity-dirs-meta)=

## Validity-based directories

Some HPGe parameters change with time, so they are stored in validity-based
directories: data files with a `validity.yaml` that says which files apply from
which date.

```{code-block} yaml
:caption: validity.yaml

- valid_from: 20220102T000000Z
  apply:
    - l200-p02-r%-T%-all-eresmod.yaml
- valid_from: 20230312T000000Z
  apply:
    - l200-p03-r%-T%-all-eresmod.yaml
```

For each run, the Simflow uses the files that apply at the start time of the
run, for the data type of the run (`phy`, `cal`, ...). The file names are free;
the ones above follow a convention. A data file maps detector names to entries,
and can give a `default` entry for all the other detectors. A missing directory,
or no file valid for the run, means that the metadata gives nothing for that
run.

:::{note}

Some metadata, such as the detector status and the skip list, is read once and
cached by the workflow. After editing it, rebuild the cache with
`--forcerun cache_modelable_hpges cache_detector_usabilities` in a separate,
earlier call.

:::

## Contents

```{toctree}
:maxdepth: 1

meta-simulations
meta-steps
meta-hpge
meta-pulse-shapes
```
