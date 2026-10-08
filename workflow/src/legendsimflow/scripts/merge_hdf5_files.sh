#!/usr/bin/env bash

# Copyright (C) 2025 Luigi Pertoldi <gipert@pm.me>
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

# Merge HDF5 files by copying the top-level objects of each input file into
# the output file. The first input is copied as a whole, so its root-group
# attributes are kept. Two inputs must not contain a top-level object with the
# same name. With no inputs, an empty HDF5 file is written.
#
# Requires the HDF5 command-line tools (h5ls, h5copy). Python with h5py is
# needed only to write the empty file.

set -euo pipefail

if [ "$#" -lt 1 ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
  echo "usage: $(basename "$0") OUTPUT [INPUT...]" >&2
  exit 1
fi

out="$1"
shift

if [ "$#" -eq 0 ]; then
  python -c "import h5py, sys; h5py.File(sys.argv[1], 'w').close()" "$out"
  exit 0
fi

cp "$1" "$out"
shift

for f in "$@"; do
  # h5ls prints one line per top-level object, "NAME   TYPE", with spaces in
  # NAME escaped as "\ ". Keep the NAME column and let `read` (without -r)
  # remove the escapes.
  # shellcheck disable=SC2162
  h5ls "$f" | sed -E 's/^((\\.|[^ \\])+) +.*$/\1/' | while IFS= read o; do
    if h5ls "$out/$o" > /dev/null 2>&1; then
      echo "error: object '/$o' of $f already exists in $out" >&2
      exit 1
    fi
    h5copy -i "$f" -o "$out" -s "/$o" -d "/$o"
  done
done
