"""Plot the HPGe impurity model and the effect of the two tuning parameters."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# crystal V04545 of legend-metadata: measurements by the vendor and the
# parabolic-exponential boule model fitted to them, densities in 1e9 e/cm^3
MEAS_Z_MM = [0, 30, 50, 70, 100, 110]
MEAS_RHO = [6, 8, 10, 12, 20, 30]
PARS = {
    "a": 5.979,
    "b": 0.06043,
    "c": 0.0003386,
    "n": 0.0005262,
    "l": 14.65,
    "m": 9.407,
}


def impurity(z, pars, slope=0.0, scale=1.0):
    """Impurity density along the boule after the slope and scale tuning."""
    a, b, c, n = pars["a"], pars["b"], pars["c"], pars["n"]
    z0, zs = pars["l"], pars["m"]
    rho = a + b * z + c * z**2 + n * np.exp((z - z0) / zs)
    rho0 = a + n * np.exp(-z0 / zs)
    return scale * (rho0 + (1 + slope) * (rho - rho0))


def plot(out: Path, dark: bool) -> None:
    fg, muted = ("#cfd0d0", "#9ca0a5") if dark else ("#1f2328", "#5a5c63")
    colors = (
        ["#4493f8", "#f0883e", "#3fb950"] if dark else ["#0969da", "#bc4c00", "#1a7f37"]
    )
    mpl.rcParams.update(
        {
            "text.color": fg,
            "axes.labelcolor": fg,
            "axes.edgecolor": muted,
            "xtick.color": muted,
            "ytick.color": muted,
            "font.size": 10,
            "svg.fonttype": "none",
        }
    )
    z = np.linspace(0, 112, 300)
    fig, axs = plt.subplots(1, 2, figsize=(8, 3.2), sharey=True, layout="constrained")
    fig.patch.set_alpha(0)

    variants = [
        (
            axs[0],
            "slope",
            [
                (0.0, colors[0], "model, $s = 0$"),
                (-0.5, colors[1], "$s = -0.5$"),
                (1.0, colors[2], "$s = +1$"),
            ],
            lambda v: {"slope": v},
        ),
        (
            axs[1],
            "depletion voltage",
            [
                (1.0, colors[0], "model"),
                (0.8, colors[1], "density $\\times$ 0.8: lower"),
                (1.2, colors[2], "density $\\times$ 1.2: higher"),
            ],
            lambda v: {"scale": v},
        ),
    ]
    for ax, title, vals, kw in variants:
        ax.patch.set_alpha(0)
        for i, (v, c, label) in enumerate(vals):
            ls = "-" if i == 0 else "--"
            ax.plot(z, impurity(z, PARS, **kw(v)), color=c, ls=ls, lw=1.6, label=label)
        ax.plot(MEAS_Z_MM, MEAS_RHO, "o", ms=4.5, color=fg, label="vendor measurements")
        ax.set_title(title, color=fg, fontsize=10.5)
        ax.set_xlabel("distance from the seed end [mm]")
        ax.legend(frameon=False, fontsize=8.5, loc="upper left")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axs[0].set_ylabel("impurity density [$10^9$ e/cm$^3$]")
    fig.savefig(out)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out_dir", type=Path)
    args = parser.parse_args()
    for theme in ("light", "dark"):
        plot(args.out_dir / f"impurity-model-{theme}.svg", dark=theme == "dark")


if __name__ == "__main__":
    main()
