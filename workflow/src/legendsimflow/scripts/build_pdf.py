# Copyright (C) 2023 Luigi Pertoldi <gipert@pm.me>
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

# ruff: noqa: F821, T201, S307

from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path

import awkward as ak
import numpy as np
import ROOT
import uproot
from legendmeta import LegendMetadata

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s: %(message)s")
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)


def process_mage_id(mage_ids: np.ndarray | list[int]):
    """
    Function to extract the channel mapping from MaGe IDs
    Parameters
    ----------
        mage_ids - a list (or numpy array) of the integer mage ids
    Returns
    ----------
        a dictionary of the form
        "channel" : {
                    mage_id : rawid,
                    ...
                    },
        "name"    : {
                    mage_id : detector name (eg V021...)
                    },
        "string"  : {
                    mage_id :string,
                    },
        "position": {
                    mage_id : position
                    },
        "type":     {
                    mage_id : type
                    }

    """
    mage_names = {"name": {}, "channel": {}, "position": {}, "string": {}, "type": {}}
    for _mage_id in mage_ids:
        m_id = str(_mage_id)
        is_ged = bool(int(m_id[0]))
        if not is_ged:
            # This should never be the case
            continue

        string = int(m_id[3:5])
        pos = int(m_id[5:7])

        for _name, _meta_dict in chmap.items():
            if _meta_dict["system"] == "geds":
                location = _meta_dict["location"]
                if location["string"] == string and location["position"] == pos:
                    mage_names["channel"][_mage_id] = f"ch{_meta_dict['daq']['rawid']}"
                    mage_names["name"][_mage_id] = _name
                    mage_names["string"][_mage_id] = int(string)
                    mage_names["position"][_mage_id] = int(pos)
                    mage_names["type"][_mage_id] = _meta_dict["type"]

    return mage_names


def get_run(text):
    pattern = re.compile(r"r\d\d\d")
    return re.findall(pattern, text)


def get_period(text):
    pattern = re.compile(r"p\d\d")
    return re.findall(pattern, text)


def get_m2_categories(channel_array, channel_to_string, channel_to_position):
    """
    Get the categories for the m2 data based on 3 categories defined in the config
    Categories:
        1 - Same string vertical neighbour
        2 - Same string not vertical neighbor
        3 - Different string
    Parameters:
        channel_array: 2D numpy array of channels
        channel2string: vectorised numpy function to convert channel into string
        channel2position: vectorised numpy function to convert channel into position
    Returns:
        categories: list of categories per event
    """

    channel_array = np.vstack(channel_array)
    channel_one = channel_array[:, 0].T
    channel_two = channel_array[:, 1].T

    ## convert to the list of strings
    string_one = channel_to_string(channel_one)
    string_two = channel_to_string(channel_two)

    same_string = string_one == string_two
    position_one = channel_to_position(channel_one)
    position_two = channel_to_position(channel_two)
    neighbour = np.abs(position_one - position_two) == 1

    is_cat_one = (same_string) & (neighbour)
    is_cat_two = (same_string) & (~neighbour)
    is_cat_three = ~same_string
    category = 1 * is_cat_one + 2 * is_cat_two + 3 * is_cat_three
    return np.array(category)


def get_string_row_diff(channel_array, channel2string, channel2position):
    """
    Get the categories for the m2 data based on 3 categories (should be in the cfg)
    1) Same string vertical neighbour
    2) Same string not vertical neighbor
    3) Different string
    Parameters:
        channel_array: 2D numpy array of channels
        channel2string: vectorised numpy function to convert channel into string
        chnanel2position: vectorised numpy function to convert channel into position
    Returns:
        categories: list of categories per event
    """

    channel_array = np.vstack(channel_array)
    channel_one = channel_array[:, 0].T
    channel_two = channel_array[:, 1].T

    ## convert to the list of strings
    string_one = channel2string(channel_one)
    string_two = channel2string(channel_two)
    string_diff_1 = (string_one - string_two) % 11
    string_diff_2 = (-string_one + string_two) % 11
    string_diff = np.array([min(a, b) for a, b in zip(string_diff_1, string_diff_2)])

    position_one = channel2position(channel_one)
    position_two = channel2position(channel_two)

    floor_diff = np.abs(position_one - position_two)

    return np.array(string_diff), np.array(floor_diff)


def get_vectorised_converter(mapping):
    """Create a vectorized function converting channel to some other quantity based on a dict
    Parameters:
        - mapping: a python dictionary of the mapping
    Return:
        - a numpy vectorised function of this mapping
    """

    def channel_to_other(mage_id):
        """Extract which string a given channel is in"""

        return mapping[mage_id]

    return np.vectorize(channel_to_other)


parser = argparse.ArgumentParser(
    prog="build_pdf", description="build LEGEND pdf files from evt tier files"
)
parser.add_argument(
    "--raw-files",
    "-r",
    default=None,
    nargs="+",
    help="path to raw simulation files for number-of-primaries determination",
)
parser.add_argument("--config", "-c", required=True, help="configuration file")
parser.add_argument("--output", "-o", required=True, help="output file name")
parser.add_argument(
    "--metadata", "-m", required=False, help="path to legend-metadata", default=None
)
parser.add_argument("input_files", nargs="+", help="evt tier files")

args = parser.parse_args()
remove_ac_hits = True

if not isinstance(args.input_files, list):
    args.input_files = [args.input_files]

with Path(args.config).open() as f:
    rconfig = json.load(f)

meta = LegendMetadata() if args.metadata is None else LegendMetadata(args.metadata)

chmap = meta.channelmap(rconfig["timestamp"])
globs = {"ak": ak, "np": np}


chmap_mage = None
channel_to_string = None
channel_to_position = None

geds_mapping = {
    f"ch{_dict['daq']['rawid']}": _name
    for _name, _dict in chmap.items()
    if chmap[_name]["system"] == "geds"
}
geds_strings = {
    f"ch{_dict['daq']['rawid']}": _dict["location"]["string"]
    for _name, _dict in chmap.items()
    if chmap[_name]["system"] == "geds"
}
geds_positions = {
    f"ch{_dict['daq']['rawid']}": _dict["location"]["position"]
    for _name, _dict in chmap.items()
    if chmap[_name]["system"] == "geds"
}
geds_types = {
    f"ch{_dict['daq']['rawid']}": _dict["type"]
    for _name, _dict in chmap.items()
    if chmap[_name]["system"] == "geds"
}

strings = np.sort([item[1] for item in geds_strings.items()])

n_primaries_total = 0

logger.info("... computing number of simulated primaries from raw files")
if args.raw_files:
    for file in args.raw_files:
        with uproot.open(f"{file}:fTree", object_cache=None) as fTree:
            n_primaries_total += fTree["fNEvents"].array(entry_stop=1)[0]
msg = f"... nprimaries {n_primaries_total}"
logger.info(msg)

# So there are many input files fed into one pdf file
# set up the hists to fill as we go along
# Creat a hist for all dets (even AC ones)

logger.info("... initializing histograms")
hists = {
    _cut_name: {
        _rawid: ROOT.TH1F(
            f"{_cut_name}_{_rawid}",
            f"{_name} energy deposits",
            rconfig["hist"]["nbins"],
            rconfig["hist"]["emin"],
            rconfig["hist"]["emax"],
        )
        for _rawid, _name in sorted(geds_mapping.items())
    }
    for _cut_name in rconfig["cuts"]
    if rconfig["cuts"][_cut_name]["is_sum"] is False
    and rconfig["cuts"][_cut_name]["is_2d"] is False
}

runs = meta.dataprod.runinfo
run_hists = {}
for _cut_name in rconfig["cuts"]:
    if not rconfig["cuts"][_cut_name]["is_sum"]:
        run_hists[_cut_name] = {}
        for _period, _run_list in runs.items():
            for run in _run_list:
                hist_name = f"{_cut_name}_{_period}_{run}"
                hist_title = f"{_period} {run} energy deposits"
                nbins = rconfig["hist"]["nbins"]
                emin = rconfig["hist"]["emin"]
                emax = rconfig["hist"]["emax"]
                run_hists[_cut_name][f"{_period}_{run}"] = ROOT.TH1F(
                    hist_name, hist_title, nbins, emin, emax
                )


# When we want to start summing the energy of events we have to treat them differently

sum_hists = {}
hists_2d = {}

## categories for m2
string_diff = np.arange(7)
floor_diff = np.arange(8)
names_m2 = [f"sd_{item1}" for item1 in string_diff]
names_m2.extend(["all", "cat_1", "cat_2", "cat_3"])
names_m2.extend(["e1_icpc", "e1_bege", "e1_ppc", "e1_coax"])

for _cut_name in rconfig["cuts"]:
    if rconfig["cuts"][_cut_name]["is_sum"] is True:
        sum_hists[_cut_name] = {}
        sum_hists[_cut_name]["all"] = ROOT.TH1F(
            f"{_cut_name}_all_summed",
            "summed energy deposits",
            rconfig["hist"]["nbins"],
            rconfig["hist"]["emin"],
            rconfig["hist"]["emax"],
        )
    if rconfig["cuts"][_cut_name]["is_2d"] is True:
        hists_2d[_cut_name] = {}

        if (
            "lar_energy" not in rconfig["cuts"][_cut_name]
            or rconfig["cuts"][_cut_name]["lar_energy"] is False
        ):
            for cat in names_m2:
                hists_2d[_cut_name][cat] = ROOT.TH2F(
                    f"{_cut_name}_{cat}_2d",
                    "energy deposits",
                    rconfig["hist"]["nbins"],
                    rconfig["hist"]["emin"],
                    rconfig["hist"]["emax"],
                    rconfig["hist"]["nbins"],
                    rconfig["hist"]["emin"],
                    rconfig["hist"]["emax"],
                )
        else:
            hists_2d[_cut_name]["all"] = ROOT.TH2F(
                f"{_cut_name}_all_2d",
                "energy deposits",
                rconfig["hist"]["nbins"],
                rconfig["hist"]["emin"],
                rconfig["hist"]["emax"],
                rconfig["hist"]["nbins"],
                rconfig["hist"]["emin"],
                rconfig["hist"]["emax"],
            )

for file_name in args.input_files:
    msg = f" >>> loading file {file_name}"
    logger.info(msg)

    ## get the run and period
    file_end = file_name.split("/")[-1]

    run = get_run(file_end)
    period = get_period(file_end)

    if len(period) > 0:
        period = period[0]
        run = run[0]
    else:
        msg = "filename doesn't contain run / period"
        raise ValueError(msg)

    # now open the file

    with uproot.open(f"{file_name}:simTree", object_cache=None) as pytree:
        if pytree.num_entries == 0:
            msg = f"ERROR: MPP evt file {file_name} has 0 events in simTree"
            raise RuntimeError(msg)

        n_primaries_total += pytree["mage_n_events"].array()[0]

        for array in pytree.iterate(step_size="100 MB"):

            array_copy = ak.copy(array)

            rng = np.random.default_rng()
            array_copy["npe_tot_poisson"] = rng.poisson(array_copy.npe_tot)

            # compute some channel mappings
            mage_ids = np.unique(ak.flatten(array_copy["mage_id"]).to_numpy())

            # get channel mappings
            chmap_mage = process_mage_id(mage_ids)
            channel_to_string = get_vectorised_converter(chmap_mage["string"])
            channel_to_position = get_vectorised_converter(chmap_mage["position"])
            channel_to_type = get_vectorised_converter(chmap_mage["type"])

            # remove below threshold hits

            array_copy["energy"] = array_copy["energy"][
                eval(f"energy > {rconfig['energy_threshold']}", globs, array_copy)
            ]
            array_copy["mage_id"] = array_copy["mage_id"][
                eval(f"energy > {rconfig['energy_threshold']}", globs, array_copy)
            ]
            array_copy["is_off"] = array_copy["is_off"][
                eval(f"energy > {rconfig['energy_threshold']}", globs, array_copy)
            ]
            array_copy["is_ac"] = array_copy["is_ac"][
                eval(f"energy > {rconfig['energy_threshold']}", globs, array_copy)
            ]

            # remove hits in off detectors
            array_copy["energy"] = array_copy["energy"][~array_copy["is_off"]]
            array_copy["mage_id"] = array_copy["mage_id"][~array_copy["is_off"]]
            array_copy["is_ac"] = array_copy["is_ac"][~array_copy["is_off"]]

            # remove records without any hit
            array_copy = array_copy[ak.num(array_copy.energy, axis=-1) > 0]

            # compute multiplicity, the length of the
            array_copy["mul"] = ak.num(array_copy["energy"], axis=-1)
            array_copy["mul_is_good"] = ak.sum(~array_copy["is_ac"], axis=-1)

            # remove events with hits in AC channels
            if remove_ac_hits:
                array_copy = array_copy[ak.all(~array_copy["is_ac"], axis=-1)]

            # loop over the cuts
            for _cut_name, _cut_dict in rconfig["cuts"].items():
                msg = f"... processing cut {_cut_name}"
                logger.debug(msg)
                _cut_string = _cut_dict["cut_string"]

                # if the cut string is empty return a copy, else query the array
                if _cut_string == "":
                    array_cut = ak.copy(array_copy)
                else:
                    array_cut = array_copy[eval(_cut_string, globs, array_copy)]

                if len(array_cut) == 0:
                    continue

                # if the cut is not sum or 2d false then flatten (by channel)
                if _cut_dict["is_sum"] is False and _cut_dict["is_2d"] is False:
                    for __mage_id in mage_ids:
                        _rawid = chmap_mage["channel"][__mage_id]

                        # extract the energy array (flattening)
                        _energy_array = (
                            ak.flatten(
                                array_cut["energy"][
                                    (array_cut["mage_id"] == __mage_id)
                                ],
                                axis=-1,
                            ).to_numpy()
                            * 1000
                        )

                        if len(_energy_array) == 0:
                            continue
                        hists[_cut_name][_rawid].FillN(
                            len(_energy_array),
                            _energy_array,
                            np.ones(len(_energy_array)),
                        )

                        # fill also time dependent hists
                        run_hists[_cut_name][f"{period}_{run}"].FillN(
                            len(_energy_array),
                            _energy_array,
                            np.ones(len(_energy_array)),
                        )

                elif _cut_dict["is_2d"] is True and (
                    "lar_energy" not in _cut_dict or _cut_dict["lar_energy"] is False
                ):

                    _energy_1_array = (
                        ak.max(array_cut["energy"], axis=-1).to_numpy() * 1000
                    )
                    _energy_2_array = (
                        ak.min(array_cut["energy"], axis=-1).to_numpy() * 1000
                    )

                    # sort it
                    _mult_channel_array = array_cut["mage_id"][
                        ak.argsort(array_cut["energy"], axis=-1)
                    ].to_numpy()

                    # loop over categories
                    for name in names_m2:
                        if name != "all" and "e1" not in name:
                            categories = get_m2_categories(
                                _mult_channel_array,
                                channel_to_string,
                                channel_to_position,
                            )
                            string_diff, floor_diff = get_string_row_diff(
                                _mult_channel_array,
                                channel_to_string,
                                channel_to_position,
                            )

                            if "cat" in name:
                                cat = int(name.split("_")[1])
                                _energy_1_array_tmp = np.array(_energy_1_array)[
                                    np.where(categories == cat)[0]
                                ]
                                _energy_2_array_tmp = np.array(_energy_2_array)[
                                    np.where(categories == cat)[0]
                                ]

                            elif "sd" in name:
                                sd = int(name.split("_")[1])

                                ids = np.where(string_diff == sd)[0]
                                _energy_1_array_tmp = np.array(_energy_1_array)[ids]
                                _energy_2_array_tmp = np.array(_energy_2_array)[ids]

                        # select by type
                        elif "e1" in name:
                            e1_mage_id = _mult_channel_array[:, 1]
                            types = channel_to_type(e1_mage_id)
                            ids = np.where(types == name.split("_")[1])
                            _energy_1_array_tmp = np.array(_energy_1_array)[ids]
                            _energy_2_array_tmp = np.array(_energy_2_array)[ids]
                        else:
                            _energy_1_array_tmp = np.array(_energy_1_array)
                            _energy_2_array_tmp = np.array(_energy_2_array)

                        if len(_energy_1_array_tmp) == 0:
                            continue
                        hists_2d[_cut_name][name].FillN(
                            len(_energy_1_array_tmp),
                            _energy_2_array_tmp,
                            _energy_1_array_tmp,
                            np.ones(len(_energy_1_array_tmp)),
                        )

                # th2 of lar energy
                elif _cut_dict["is_2d"] is True and _cut_dict["lar_energy"] is True:
                    _energy_array = (
                        ak.sum(array_cut["energy"], axis=-1).to_numpy() * 1000
                    )
                    _lar_energy_array = array_cut["edep_lar_MeV"].to_numpy() * 1000

                    if len(_energy_array) == 0:
                        continue
                    hists_2d[_cut_name]["all"].FillN(
                        len(_energy_array),
                        _energy_array,
                        _lar_energy_array,
                        np.ones(len(_energy_array)),
                    )
                else:
                    _summed_energy_array = (
                        ak.sum(array_cut["energy"], axis=-1).to_numpy() * 1000
                    )

                    if len(_summed_energy_array) == 0:
                        continue

                    sum_hists[_cut_name]["all"].FillN(
                        len(_summed_energy_array),
                        _summed_energy_array,
                        np.ones(len(_summed_energy_array)),
                    )

# The individual channels have been filled
# now add them together to make the grouped hists
# We don't need to worry about the AC dets as they will have zero entries
logger.info("... making grouped pdfs")
for _cut_name, _ in hists.item():
    hists[_cut_name]["all"] = ROOT.TH1F(
        f"{_cut_name}_all",
        "All energy deposits",
        rconfig["hist"]["nbins"],
        rconfig["hist"]["emin"],
        rconfig["hist"]["emax"],
    )
    for _type in ["bege", "coax", "icpc", "ppc"]:
        hists[_cut_name][_type] = ROOT.TH1F(
            f"{_cut_name}_{_type}",
            f"All {_type} energy deposits",
            rconfig["hist"]["nbins"],
            rconfig["hist"]["emin"],
            rconfig["hist"]["emax"],
        )
    for _rawid, _name in geds_mapping.items():
        hists[_cut_name][chmap[geds_mapping[_rawid]]["type"]].Add(
            hists[_cut_name][_rawid]
        )
        hists[_cut_name]["all"].Add(hists[_cut_name][_rawid])


# write the hists to file (but only if they have none zero entries)
# Changes the names to drop type_ etc
msg = f"... writing to file {args.output}"
logger.info(msg)
out_file = uproot.recreate(args.output)
for _cut_name, _hist_dict in hists.items():
    dir = out_file.mkdir(_cut_name)
    for key, item in _hist_dict.items():
        dir[key] = item

## fill run based histos
for _cut_name, _hist_dict in run_hists.items():
    for key, item in _hist_dict.items():
        out_file[_cut_name + "/" + key] = item

# All other hists
for dict in [sum_hists, hists_2d]:
    for _cut_name, _sub_dirs in dict.items():
        dir = out_file.mkdir(_cut_name)

        for _sub_dir, _hist in _sub_dirs.items():
            dir[_sub_dir] = _hist

msg = f"... nprimaries {n_primaries_total}"
logger.info(msg)
out_file["number_of_primaries"] = str(int(n_primaries_total))
out_file.close()
