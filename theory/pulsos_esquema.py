"""Schematic of the repeated-lensing pulse train, for theory.pdf.

Black and white, no numbers on either axis. The shape is generic for a
point-mass lens in the diffractive regime (w < 1) crossed by a source on a
near-edge-on circular orbit: one amplification pulse per outer period, each
flanked by the diffraction ringing that geometric optics does not produce.
Case B's figure of the same quantity carries this project's actual numbers
(cases/case_B_monochromatic/caseB_repeated_pulses.png); this one deliberately
does not, since theory.pdf states the phenomenon and not the instance.

Two panels: two periods on the left, so that "one pulse per period" is a
statement the reader can see rather than take, and one pulse enlarged on the
right, where the ringing is resolved. No unlensed reference line -- with
|F|^2 on the axis the unlensed case is the value 1, and drawing it as a
dashed rule through the middle of the figure added a horizontal line and no
information.

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


def amplification(t):
    """|F(w, y(t))|^2 on the given times."""
    y, lensed, _ = geo.impact_parameter_of_time(
        t, system.A_OUT_M / units.PC_SI, system.E_OUT, system.P_OUT_S,
        np.deg2rad(system.I_OUT_DEG), system.OMEGA_OUT, system.LITTLE_OMEGA_OUT,
        system.T_PERI_OUT, system.M_LENS_MSUN, system.D_L_PC)
    w = wo.w_of_frequency(system.F_B_HZ, system.M_LENS_MSUN)
    mag2 = np.ones_like(t)
    near = lensed & (y < 60.0)
    mag2[near] = np.abs([complex(wo.F_hybrid(w, yy)) for yy in y[near]]) ** 2
    return mag2


def esqueleto(ax):
    """Schematic axes: the shape is the content, the values are not."""
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.8)


def main():
    P = system.P_OUT_S
    n_periods = 2

    t = np.linspace(0.0, n_periods * P, 8001, endpoint=False)
    mag2 = amplification(t)
    tp = t / P

    # one pulse, enlarged: conjunction sits at a quarter period
    half = 0.075 * P
    t_z = np.linspace(0.25 * P - half, 0.25 * P + half, 4001)
    mag2_z = amplification(t_z)
    tz = (t_z - 0.25 * P) / P

    fig, (ax, axz) = plt.subplots(
        1, 2, figsize=(6.6, 2.4), gridspec_kw={"width_ratios": [1.7, 1.0]})

    ax.plot(tp, mag2, color="black", lw=0.9)
    peaks = [0.25, 1.25]
    ax.annotate("", xy=(peaks[0], 1.47), xytext=(peaks[1], 1.47),
                arrowprops=dict(arrowstyle="<->", color="black", lw=0.7))
    ax.text(0.5 * sum(peaks), 1.50, r"$P_\mathrm{ext}$",
            ha="center", va="bottom", fontsize=9)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$|F(w,y(t))|^{2}$")
    ax.set_xlim(-0.05, 1.55)
    ax.set_ylim(0.89, 1.56)
    esqueleto(ax)

    axz.plot(tz, mag2_z, color="black", lw=0.9)
    axz.annotate("anillado de\ndifracción", xy=(0.030, 1.055),
                 xytext=(0.040, 1.26), fontsize=8, ha="left",
                 arrowprops=dict(arrowstyle="->", color="black", lw=0.6))
    axz.set_xlabel(r"$t$ (ampliado)")
    axz.set_xlim(tz[0], tz[-1])
    axz.set_ylim(0.89, 1.56)
    esqueleto(axz)

    fig.tight_layout()
    out = HERE / "pulsos_esquema.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
