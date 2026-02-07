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
import Pkg

if !ispath(joinpath(first(DEPOT_PATH), "registries", "LegendJuliaRegistry"))
    @info("Installing Legend Julia package registry")
    Pkg.Registry.add("General")
    Pkg.Registry.add(url = "https://github.com/legend-exp/LegendJuliaRegistry.git")
end

Pkg.instantiate()
Pkg.precompile()
