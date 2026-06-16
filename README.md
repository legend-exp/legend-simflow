# legend-simflow

<img src=".github/logo.jpg" alt="legend-simflow logo" align="left" height="170">

[![PyPI](https://img.shields.io/pypi/v/legend-simflow?logo=pypi)](https://pypi.org/project/legend-simflow/)
![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/legend-exp/legend-simflow?logo=git)
[![GitHub Workflow Status](https://img.shields.io/github/checks-status/legend-exp/legend-simflow/main?label=main%20branch&logo=github)](https://github.com/legend-exp/legend-simflow/actions)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Codecov](https://img.shields.io/codecov/c/github/legend-exp/legend-simflow?logo=codecov)](https://app.codecov.io/gh/legend-exp/legend-simflow)
[![Read the Docs](https://img.shields.io/readthedocs/legend-simflow?logo=readthedocs)](https://legend-simflow.readthedocs.io)
![GitHub issues](https://img.shields.io/github/issues/legend-exp/legend-simflow?logo=github)
![GitHub pull requests](https://img.shields.io/github/issues-pr/legend-exp/legend-simflow?logo=github)
![License](https://img.shields.io/github/license/legend-exp/legend-simflow)

End-to-end Snakemake workflow to run Monte Carlo simulations of signal and
background signatures in the LEGEND experiment and produce probability-density
functions (pdfs). Configuration metadata (e.g. rules for generating simulation
macros or post-processing settings) is stored at
[legend-simflow-config](https://github.com/legend-exp/legend-simflow-config).

## Features

- Tier-based Snakemake workflow taking Geant4
  ([_remage_](https://remage.readthedocs.io)) Monte Carlo events all the way to
  analysis-ready pdfs.
- Simulated statistics weighted by the livetime of a user-selected list of data
  taking runs (run partitioning).
- Fully metadata-driven: a production is configured by editing a single YAML
  file, no code required. Runs locally or on HPC sites through ready-made
  Snakemake profiles.

### Detector response models

Physics and detector models tuned to real LEGEND-200 data and applied during
post-processing:

- **HPGe energy**: per-detector energy scale and measured energy resolution
  FWHM(E) used to smear the simulated energy.
- **HPGe active volume**: dead-layer / active-volume model from detector
  geometry and metadata.
- **HPGe pulse shape and PSD**: drift-time maps and ideal pulse-shape libraries
  computed with
  [_SolidStateDetectors.jl_](https://juliaphysics.github.io/SolidStateDetectors.jl),
  convolved with an electronics-response model fitted to data superpulses; A/E
  response from current-signal templates and measured A/E resolution, with
  data-driven PSD survival cuts.
- **Liquid-argon scintillation and SiPMs**: scintillation photon generation,
  photoelectron detection sampled from optical maps, per-photoelectron amplitude
  resolution, and time clustering reproducing the SiPM time response.
- **Detector status**: per-run usability and PSD-usability flags.
- **Event building**: time-coincidence maps (TCM) across detectors to group hits
  into physics events.

## Documentation

Full documentation is hosted at
[legend-simflow.readthedocs.io](https://legend-simflow.readthedocs.io). Good
entry points:

- [Overview and key concepts](https://legend-simflow.readthedocs.io/en/latest/)
- [Installation and configuration](https://legend-simflow.readthedocs.io/en/latest/manual/setup.html)
- [Running a production](https://legend-simflow.readthedocs.io/en/latest/manual/prod.html)
- [Configuration metadata and file naming](https://legend-simflow.readthedocs.io/en/latest/manual/meta.html)
