# ruff: noqa: I002

# Copyright (C) 2026 Luigi Pertoldi <gipert@pm.me>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Simulate mono-energetic electrons in the bulk of all HPGe detectors.

The electron-gun simulations feed the A/E mean energy-dependence extraction
(see :mod:`.extract_hpge_aoemean_energy_dependence`). This script builds the
production geometry, renders the remage macro from the ``aoemeanmod`` macro
template and runs remage (through its Python API) once per electron energy,
then plots the primary vertices as a check of the confinement.
"""

import argparse
import sys
from pathlib import Path

import legenddataflowscripts as ldfs
import legenddataflowscripts.utils  # ensures ldfs.utils is loaded
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import commands, geometry, nersc, patterns, utils
from legendsimflow.metadata import electron_gun_primaries
from legendsimflow.plot import plot_primary_vertices
from legendsimflow.scripts import log_script_invocation


@snakemake_compatible(
    mapping={
        "geom_config": "input.geom_config",
        "geom_file": "output.geom",
        "macro_file": "output.macro",
        "stp_files": "output.stp_files",
        "plot_file": "output.plot_file",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Simulate mono-energetic electrons in the bulk of all HPGe detectors."
    )
    parser.add_argument(
        "--geom-config", required=True, help="geometry configuration YAML file"
    )
    parser.add_argument("--geom-file", required=True, help="output GDML geometry file")
    parser.add_argument("--macro-file", required=True, help="output remage macro file")
    parser.add_argument(
        "--stp-files",
        nargs="+",
        required=True,
        help=(
            "output stp tier files, one per electron energy; the energy is read "
            "from the file name"
        ),
    )
    parser.add_argument(
        "--plot-file", required=True, help="output primary vertices plot file"
    )
    parser.add_argument("--log-file", default=None, help="log file")
    parser.add_argument(
        "--simflow-config",
        "--config",
        dest="simflow_config",
        required=True,
        help="simflow config YAML path",
    )
    args = parser.parse_args()

    # imported here so that the module can be imported (and the CLI inspected)
    # in environments without remage, e.g. to collect the tests
    import remage  # noqa: PLC0415

    config = utils.init_simflow_context(args.simflow_config, workflow=None).config
    metadata = config.metadata

    log = ldfs.utils.build_log(metadata.simprod.config.logging, args.log_file)
    log_script_invocation(log, "simulate-electron-gun", parser, args)

    geom_config = nersc.dvs_ro(config, args.geom_config)
    geom_file = Path(args.geom_file)
    macro_file = Path(args.macro_file)
    stp_files = [Path(f) for f in args.stp_files]

    # same logic as the stp tier (see commands.remage_run)
    n_events = electron_gun_primaries(config)
    if config.benchmark.get("enabled", False):
        n_events = config.benchmark.n_primaries["stp"]

    # Snakemake redirects the standard streams to logger proxies without
    # isatty(), which the remage logging setup queries
    for stream in (sys.stdout, sys.stderr):
        if not hasattr(stream, "isatty"):
            stream.isatty = lambda: False

    log.info("building the geometry")
    geometry.build_gdml(config, geom_config, geom_file)

    log.info("rendering the remage macro")
    commands.make_electron_gun_macro(config, geom_file, macro_file)

    for stp_file in stp_files:
        energy = patterns.electron_gun_energy_from_path(stp_file)
        log.info("simulating %d electrons of %d keV", n_events, energy)

        output, move2cfs = nersc.make_on_scratch(config, stp_file)
        remage.remage_run(
            [str(macro_file)],
            gdml_files=[str(geom_file)],
            output=str(output),
            macro_substitutions={
                "N_EVENTS": str(n_events),
                "SEED": str(
                    utils.string_to_remage_seed(
                        str(stp_file), seed=config.get("simflow_rng_seed", 0)
                    )
                ),
                commands.ELECTRON_GUN_ENERGY_ALIAS: str(energy),
            },
            threads=1,
            overwrite_output=True,
            merge_output_files=True,
            log_level="detail",
        )
        move2cfs()

    log.info("plotting the primary vertices")
    Path(args.plot_file).parent.mkdir(parents=True, exist_ok=True)
    plot_primary_vertices([nersc.dvs_ro(config, f) for f in stp_files], args.plot_file)


if __name__ == "__main__":
    main()
