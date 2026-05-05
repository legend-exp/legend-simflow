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

import tarfile
from pathlib import Path


def create_plots_tarball(generated_dir: Path, output: Path, prefix: str) -> None:
    """Archive all plots/ directories under generated_dir into a .tar.xz.

    Parameters
    ----------
    generated_dir
        The ``generated/`` directory of the production cycle.
    output
        Path to write the ``.tar.xz`` tarball.
    prefix
        Prefix directory name inside the archive (e.g. ``prod-v1-plots``).
    """
    plots_dirs = sorted(
        p for p in generated_dir.rglob("*") if p.is_dir() and p.name == "plots"
    )

    output.parent.mkdir(parents=True, exist_ok=True)

    with tarfile.open(output, "w:xz") as tar:
        for plots_dir in plots_dirs:
            rel = plots_dir.relative_to(generated_dir)
            tar.add(plots_dir, arcname=f"{prefix}/{rel}")


def create_pdfs_tarball(pdf_dir: Path, output: Path, prefix: str) -> None:
    """Archive all pdf tier LH5 files into a .tar.xz.

    Parameters
    ----------
    pdf_dir
        The ``tier/pdf/`` directory of the production cycle.
    output
        Path to write the ``.tar.xz`` tarball.
    prefix
        Prefix directory name inside the archive (e.g. ``prod-v1-pdfs``).
    """
    output.parent.mkdir(parents=True, exist_ok=True)

    with tarfile.open(output, "w:xz") as tar:
        tar.add(pdf_dir, arcname=prefix)
