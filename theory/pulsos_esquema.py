"""Schematic of the repeated-lensing pulse train, for theory.pdf.

Black and white, no numbers on either axis. The shape is generic for a
point-mass lens in the diffractive regime (w < 1) crossed by a source on a
near-edge-on circular orbit: one amplification pulse per outer period, each
flanked by the diffraction ringing that geometric optics does not produce.
Case B's figure of the same quantity carries this project's actual numbers
(cases/case_B_monochromatic/caseB_repeated_pulses.png); this one deliberately
does not, since theory.pdf states the phenomenon and not the instance.

Writes pulsos_esquema.pdf next to this file. Run: python3 theory/pulsos_esquema.py
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))

from gwlens import geometry as geo, system, units
from gwlens import waveoptics as wo


def main():
    n_periods = 3
    t = np.linspace(0.0, n_periods * system.P_OUT_S, 6001, endpoint=False)
    y, lensed, _ = geo.impact_parameter_of_time(
        t, system.A_OUT_M / units.PC_SI, system.E_OUT, system.P_OUT_S,
        np.deg2rad(system.I_OUT_DEG), system.OMEGA_OUT, system.LITTLE_OMEGA_OUT,
        system.T_PERI_OUT, system.M_LENS_MSUN, system.D_L_PC)
    w = wo.w_of_frequency(system.F_B_HZ, system.M_LENS_MSUN)

    mag2 = np.ones_like(t)
    near = lensed & (y < 60.0)
    mag2[near] = np.abs([complex(wo.F_hybrid(w, yy)) for yy in y[near]]) ** 2

    tp = t / system.P_OUT_S

    fig, ax = plt.subplots(figsize=(6.4, 2.5))
    ax.plot(tp, mag2, color="black", lw=0.9)
    ax.axhline(1.0, color="black", lw=0.6, ls=(0, (4, 3)))

    # The period, marked between two consecutive peaks rather than stated.
    peaks = [0.25 + k for k in range(n_periods)]
    ax.annotate("", xy=(peaks[0], 1.46), xytext=(peaks[1], 1.46),
                arrowprops=dict(arrowstyle="<->", color="black", lw=0.7))
    ax.text(0.5 * (peaks[0] + peaks[1]), 1.49, r"$P_\mathrm{ext}$",
            ha="center", va="bottom", fontsize=9)

    ax.text(0.02, 0.965, r"$|F|^{2}=1$: sin lente", fontsize=8,
            va="top", ha="left")
    ax.annotate("anillado de difracción", xy=(peaks[1] + 0.16, 1.05),
                xytext=(peaks[1] + 0.34, 1.30), fontsize=8,
                arrowprops=dict(arrowstyle="->", color="black", lw=0.6))

    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$|F(w,y(t))|^{2}$")
    ax.set_xlim(-0.06, n_periods)
    ax.set_ylim(0.86, 1.62)
    # Schematic: the shape is the content, the values are not. No tick
    # labels on either axis, and only the two spines that carry meaning.
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.8)

    fig.tight_layout()
    out = HERE / "pulsos_esquema.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
