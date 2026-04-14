# ruff: noqa: I002

# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>,
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

import argparse

import dbetto
import legenddataflowscripts as ldfs
import legenddataflowscripts.utils
from snakemake_argparse_bridge import snakemake_compatible

from legendsimflow import hpge_pars, utils
from legendsimflow import metadata as mutils
from legendsimflow.scripts import log_script_invocation


@snakemake_compatible(
    mapping={
        "runid": "wildcards.runid",
        "eresmod_file": "output.eresmod_file",
        "aoeresmod_file": "output.aoeresmod_file",
        "psdcuts_file": "output.psdcuts_file",
        "log_file": "log[0]",
        "simflow_config": "config",
    }
)
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract HPGe observable models for a LEGEND run."
    )
    parser.add_argument(
        "--runid",
        required=True,
        help="LEGEND run identifier (e.g. l200-p03-r000-phy)",
    )
    parser.add_argument(
        "--eresmod-file",
        required=True,
        help="output YAML file for energy resolution model parameters",
    )
    parser.add_argument(
        "--aoeresmod-file",
        required=True,
        help="output YAML file for A/E resolution model parameters",
    )
    parser.add_argument(
        "--psdcuts-file",
        required=True,
        help="output YAML file for PSD cut values",
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

    config = utils.init_simflow_context(args.simflow_config, workflow=None).config
    metadata = config.metadata

    log = ldfs.utils.build_log(metadata.simprod.config.logging, args.log_file)
    log_script_invocation(log, "extract-hpge-obs-models", parser, args)

    runid = args.runid
    l200data = config.paths.get("l200data", None)

    # Collection step: gather what is available; completeness checked in hit.py.
    # For each observable, behaviour depends on simprod/config/pars/geds/<obs>/:
    #  - absent: use l200data only
    #  - present, no "default": l200data base + per-detector overrides
    #  - present, with "default": expand all channelmap geds detectors;
    #    l200data not consulted
    raw_eresmod = mutils.simpars(metadata, "geds.eresmod", runid, default=None)
    raw_aoeresmod = mutils.simpars(metadata, "geds.aoeresmod", runid, default=None)
    raw_psdcuts = mutils.simpars(metadata, "geds.psdcuts", runid, default=None)

    eresmod_default = (
        raw_eresmod.get("default", None) if raw_eresmod is not None else None
    )
    aoeresmod_default = (
        raw_aoeresmod.get("default", None) if raw_aoeresmod is not None else None
    )
    psdcuts_default = (
        raw_psdcuts.get("default", None) if raw_psdcuts is not None else None
    )

    # resolve the channelmap once if any observable needs to expand all geds detectors
    chmap = None
    if (
        eresmod_default is not None
        or aoeresmod_default is not None
        or psdcuts_default is not None
    ):
        tstamp = mutils.runinfo(metadata, runid).start_key
        chmap = metadata.channelmap(tstamp, skip_version_check=True)

    # --- energy resolution model ---
    if eresmod_default is not None:
        # metadata with default: expand all channelmap geds detectors
        msg = f"using eresmod defaults from simprod/config/pars for {runid}"
        log.info(msg)
        eresmod_dict = {
            name: raw_eresmod.get(name, eresmod_default).to_dict()
            for name, info in chmap.items()
            if getattr(info, "system", None) == "geds"
        }
        dbetto.utils.write_dict(eresmod_dict, args.eresmod_file)
    else:
        # l200data base (always), then apply per-detector metadata overrides if any
        if l200data is None:
            msg = (
                "l200data is not configured and no eresmod metadata with a 'default' "
                "key was found — cannot extract energy resolution model"
            )
            raise RuntimeError(msg)
        msg = f"extracting eresmod from l200data for {runid}"
        log.info(msg)
        hit_tier_name = utils.get_hit_tier_name(l200data)
        pars_db = utils.init_generated_pars_db(l200data, tier=hit_tier_name, lazy=True)

        eres_pars_dict = hpge_pars.lookup_energy_res_metadata(
            l200data,
            metadata,
            runid,
            hit_tier_name=hit_tier_name,
            pars_db=pars_db,
        )

        out_dict = dbetto.AttrsDict({})
        fields = ["expression", "parameters", "uncertainties"]
        for hpge, meta in eres_pars_dict.items():
            out_dict[hpge] = {f: meta[f] for f in fields}

        if raw_eresmod is not None:
            msg = "applying per-detector eresmod overrides from simprod/config/pars"
            log.info(msg)
            for hpge, meta in raw_eresmod.items():
                out_dict[hpge] = meta.to_dict()

        dbetto.utils.write_dict(out_dict.to_dict(), args.eresmod_file)

    # --- A/E resolution model ---
    if aoeresmod_default is not None:
        msg = f"using aoeresmod defaults from simprod/config/pars for {runid}"
        log.info(msg)
        aoeresmod_dict = {
            name: raw_aoeresmod.get(name, aoeresmod_default).to_dict()
            for name, info in chmap.items()
            if getattr(info, "system", None) == "geds"
        }
        dbetto.utils.write_dict(aoeresmod_dict, args.aoeresmod_file)
    else:
        if l200data is None:
            msg = (
                "l200data is not configured and no aoeresmod metadata with a 'default' "
                "key was found — cannot extract A/E resolution model"
            )
            raise RuntimeError(msg)
        msg = f"extracting aoeresmod from l200data for {runid}"
        log.info(msg)
        hit_tier_name = utils.get_hit_tier_name(l200data)
        pars_db = utils.init_generated_pars_db(l200data, tier=hit_tier_name, lazy=True)

        aoeres_pars_dict = hpge_pars.lookup_aoe_res_metadata(
            l200data,
            metadata,
            runid,
            hit_tier_name=hit_tier_name,
            pars_db=pars_db,
        )

        out_dict = dbetto.AttrsDict({})
        for hpge, meta in aoeres_pars_dict.items():
            out_dict[hpge] = {
                "expression": meta["expression"],
                "parameters": meta["pars"],
                "errs": meta["errs"],
            }

        if raw_aoeresmod is not None:
            msg = "applying per-detector aoeresmod overrides from simprod/config/pars"
            log.info(msg)
            for hpge, meta in raw_aoeresmod.items():
                out_dict[hpge] = meta.to_dict()

        dbetto.utils.write_dict(out_dict.to_dict(), args.aoeresmod_file)

    # --- PSD cut values ---
    if psdcuts_default is not None:
        msg = f"using psdcuts defaults from simprod/config/pars for {runid}"
        log.info(msg)
        psdcuts_dict = {
            name: raw_psdcuts.get(name, psdcuts_default).to_dict()
            for name, info in chmap.items()
            if getattr(info, "system", None) == "geds"
        }
        dbetto.utils.write_dict(psdcuts_dict, args.psdcuts_file)
    else:
        if l200data is None:
            msg = (
                "l200data is not configured and no psdcuts metadata with a 'default' "
                "key was found — cannot extract PSD cut values"
            )
            raise RuntimeError(msg)
        msg = f"extracting PSD cut values from l200data for {runid}"
        log.info(msg)
        hit_tier_name = utils.get_hit_tier_name(l200data)
        pars_db = utils.init_generated_pars_db(l200data, tier=hit_tier_name, lazy=True)

        aoecuts_pars_dict = hpge_pars.lookup_psd_cut_values(
            l200data,
            metadata,
            runid,
            hit_tier_name=hit_tier_name,
            pars_db=pars_db,
        )

        out_dict = dbetto.AttrsDict(aoecuts_pars_dict)

        if raw_psdcuts is not None:
            msg = "applying per-detector psdcuts overrides from simprod/config/pars"
            log.info(msg)
            for hpge, meta in raw_psdcuts.items():
                out_dict[hpge] = meta.to_dict()

        dbetto.utils.write_dict(out_dict.to_dict(), args.psdcuts_file)


if __name__ == "__main__":
    main()
