# Validation

A production must be checked before its output is used. The Simflow helps with
two products: control plots, made by the steps along the way, and a report of
the warnings logged during the run. Both should be inspected after every
production.

## Control plots

The steps write control plots next to their output, in `plots/` folders:

| Step                | Plots                                                                                                                                                                                                                                                          |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `stp`               | the positions of the primary vertices of the first job, a rendering of the geometry, and the simulated against the measured mass of each HPGe detector                                                                                                         |
| `par`               | the drift-time maps, the data superpulses, the realistic pulse shapes, and the fits of the current-pulse template, the electronics response, the A/E mean energy dependence and the impurity profiles, for the parameters a production builds (see [](par.md)) |
| `opt`, `hit`, `cvt` | the distributions of the observables of each simulation                                                                                                                                                                                                        |

They show, for example, whether the primaries start in the intended volumes,
whether the geometry and the detector masses are right, and whether the fits of
the `par` step converged. The `stp` plots are made with priority, so they can be
checked early, while the production runs. To collect all the plots of a
production in one archive, see [](prod.md).

## Warnings report

At the end of a successful run, the Simflow collects every line of the log files
of the run that contains a warning into one file,
`paths.log/<time>/warnings-report.log`, grouped by log file and with the line
numbers. `paths.log/latest` points to the folder of the last run, and the number
of warnings is printed at the end of the run.

Many warnings point to a choice the Simflow made in place of the user, for
example a detector that got a fallback parameter (see
{ref}`build-tier-hit-hpge`), a detector missing from the metadata of a run, or a
detector without a PSD model. Read the report and understand every warning
before using the output.
