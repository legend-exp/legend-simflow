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
production geometry and runs remage once per electron energy. remage writes one
file per run, so the energies end up in separate files, named after the energy
they were simulated at.

The energies to simulate come from the command line only: either from the names
of the output files (``--stp-files``, as the Snakemake rule passes them) or from
``--energies``, in which case the script names the files itself in
``--output-dir``.
"""

import argparse
import sys
from pathlib import Path

import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils  # ensures ldfs.utils is loaded
import pygeomtools
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import commands, geometry, nersc, patterns, utils
from legendsimflow.metadata import electron_gun_macro_template, electron_gun_primaries
from legendsimflow.scripts import log_script_invocation


@snakemake_compatible(
    mapping={
        "geom_config": "input.geom_config",
        "stp_files": "output.stp_files",
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
    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument(
        "--stp-files",
        nargs="+",
        help=(
            "output stp tier files, one per electron energy; the energy is read "
            "from the file name"
        ),
    )
    output.add_argument(
        "--output-dir",
        help="directory for the output stp tier files, which are named by energy",
    )
    parser.add_argument(
        "--energies",
        nargs="+",
        type=int,
        metavar="KEV",
        help="electron kinetic energies to simulate, in keV (only with --output-dir)",
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

    if args.stp_files is not None and args.energies is not None:
        parser.error("--energies can only be used together with --output-dir")
    if args.output_dir is not None and args.energies is None:
        parser.error("--energies is required together with --output-dir")

    # imported here so that the module can be imported (and the CLI inspected)
    # in environments without remage, e.g. to collect the tests
    import remage  # noqa: PLC0415

    config = utils.init_simflow_context(args.simflow_config, workflow=None).config

    log = ldfs.utils.build_log(config.metadata.simprod.config.logging, args.log_file)
    log_script_invocation(log, "simulate-electron-gun", parser, args)

    # HACK: remage decides whether to colour its output with
    # sys.stderr.isatty(), but build_log() above (and Snakemake, when the script
    # runs as a rule) replaces the stream with a proxy that has no isatty().
    # Drop this once remage guards that call.
    for stream in (sys.stdout, sys.stderr):
        if not hasattr(stream, "isatty"):
            stream.isatty = lambda: False

    if args.stp_files is not None:
        stp_files = [Path(f) for f in args.stp_files]
    else:
        out_dir = Path(args.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        stp_files = [
            out_dir / patterns.output_electron_gun_stp_filename(config, energy=e).name
            for e in args.energies
        ]

    # the GDML is written for the extraction step, which reads the detector
    # metadata back from it (it differs from the plain legend-metadata one:
    # the geometry generator reworks it while building)
    log.info("building the geometry")
    geom_config_file = nersc.dvs_ro(config, args.geom_config)
    geom_config = dbetto.utils.load_dict(geom_config_file)
    # expand $_ to the directory holding the template, as dbetto does
    dbetto.Props.subst_vars(
        geom_config, var_values={"_": Path(geom_config_file).parent.resolve()}
    )
    registry = geometry.construct_geometry(config, geom_config)

    geom_file = patterns.output_electron_gun_geom_filename(config)
    if args.output_dir is not None:
        geom_file = Path(args.output_dir) / geom_file.name
    geom_file.parent.mkdir(parents=True, exist_ok=True)
    pygeomtools.write_pygeom(registry, geom_file)

    # the germanium volumes are the only macro commands the Simflow fills in:
    # they are read from the geometry, everything else is in the template
    volumes = sorted(pygeomtools.detectors.get_all_sensvols(registry, "germanium"))
    if not volumes:
        msg = f"no germanium sensitive volumes in the {config.experiment} geometry"
        raise RuntimeError(msg)

    confinement = ["/RMG/Generator/Confine Volume"]
    confinement += [
        f"/RMG/Generator/Confinement/Physical/AddVolume {v}" for v in volumes
    ]

    macro_text = commands.render_macro_template(
        electron_gun_macro_template(config), {"CONFINEMENT": "\n".join(confinement)}
    )
    n_events = electron_gun_primaries(config)

    for stp_file in stp_files:
        energy = patterns.electron_gun_energy_from_path(stp_file)
        output, move2cfs = nersc.make_on_scratch(config, stp_file)

        log.info("simulating %d electrons of %d keV", n_events, energy)
        remage.remage_run(
            macro_text,
            gdml_files=str(geom_file),
            output=str(output),
            macro_substitutions={
                "N_EVENTS": str(n_events),
                # the template asks for the electron energy by this alias
                "ENERGY_KEV": str(energy),
            },
            threads=1,
            overwrite_output=True,
            merge_output_files=True,
            log_level="detail",
        )

        move2cfs()


if __name__ == "__main__":
    main()
