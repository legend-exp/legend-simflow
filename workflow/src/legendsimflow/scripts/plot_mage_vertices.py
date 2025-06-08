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

from __future__ import annotations

import argparse

import ROOT

parser = argparse.ArgumentParser(
    prog="plot_mage_vertices", description="plot MaGe primary event vertices"
)
parser.add_argument("--output", "-o", help="output file name")
parser.add_argument(
    "--batch", "-b", action="store_true", help="run in graphics batch mode"
)
parser.add_argument("input_files", nargs="+", help="MaGe output files")

args = parser.parse_args()

if not isinstance(args.input_files, list):
    args.input_files = [args.input_files]

if args.batch:
    ROOT.gROOT.SetBatch()

# do not display any of the standard histogram decorations
ROOT.gStyle.SetOptStat(False)
ROOT.gStyle.SetOptFit(0)

tree = ROOT.TChain("fTree")
for f in args.input_files:
    tree.Add(f)

c = ROOT.TCanvas("c", "MaGe event vertices", 1000, 1000)
c.Divide(2, 2)

exprs = [
    "eventPrimaries.fSteps.fY:eventPrimaries.fSteps.fX",
    "eventPrimaries.fSteps.fZ:eventPrimaries.fSteps.fY",
    "eventPrimaries.fSteps.fZ:eventPrimaries.fSteps.fX",
    "eventPrimaries.fSteps.fZ:eventPrimaries.fSteps.fY:eventPrimaries.fSteps.fX",
]

for i, expr in enumerate(exprs):
    p = c.cd(i + 1)
    tree.Draw(expr)
    p.SetGrid()

if args.output:
    c.SaveAs(args.output)
