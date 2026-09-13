# The figures, one by one

What each figure is for, what it shows, and what it deliberately does not.
Eight figures: four per case, plus one schematic in `theory/`. The figures
themselves are in Spanish, because they go into `report.pdf`; this page is in
English like the rest of the repo's own documentation.

Two of them answer the same question at different depths, and the difference
is worth stating once here because it is the most common confusion:

- **`caseB_repeated_pulses.png` is the lens alone.** It plots `|F|²`, the
  amplification factor — a property of the lens and the geometry, with no
  waveform in it at all.
- **`caseB_detector_view.png` is the measured signal.** It plots the strain,
  `|F|` multiplied into an actual waveform. The pulses in it are the pulses of
  the first figure, riding on a carrier.

---

## Case A — the chirp, static lens

The inner binary sweeps 10 → 125.6 Hz in 15.2 s, against a 4-day outer
period. The lens does not move: the whole signal arrives at one fixed impact
parameter, `y_A = 1.589`. `w` runs from 62 to 778, so this is geometric
optics: two images, one delayed by `ΔT = 3.43 s`.

### `caseA_F_of_f.png` — the amplification across the band

Three panels. The first two are `|F(f)|` over 3 Hz windows at **opposite ends
of the band**, at 10 Hz (`w = 62`) and at `f_isco` (`w = 778`); the third is
`arg F(f)` over the first window.

Two windows rather than one because the claim being made is that the envelope
width and the fringe period are *the same everywhere in the band* — that is
the geometric-optics statement, and one window cannot show it. The dashed
lines are the exact envelope `√μ₊ ± √|μ₋|`, and the fringe period is
`1/ΔT = 0.29 Hz`, independent of `f`.

The previous version had a full-band panel instead, which came out as a
featureless blue rectangle: it asserted that the envelope is constant rather
than showing it. The curves were also drawn from the waveform's own FFT grid,
7 points per fringe, and looked polygonal. `F(f)` is a closed form, so these
panels now use their own dense grid.

### `caseA_strain_time.png` — what the lens does to the waveform

Four panels: unlensed strain, lensed strain **on the same scale**, a zoom
centred on the merger, and the **difference**.

The first two are separate panels rather than one overlay. At this time
resolution the two curves are indistinguishable — which is the point, and is
worth seeing — and drawing one on top of the other simply hid the lower one.

The fourth panel is where the lensing becomes visible, and it is the answer to
"does subtracting the two make sense?". It does, and for a specific reason. In
geometric optics `F = √μ₊ − i√|μ₋| e^{iwΔT}`, so

    h_lensed − h_unlensed = (√μ₊ − 1)·h(t) + √|μ₋|·h(t − ΔT)

The first term is 3% of the signal and rides under the chirp. The second is a
**full copy** of the signal at 24%, delayed by `ΔT`. So the subtraction
isolates the second image and shows it replaying the entire inspiral, merger
and ringdown against a flat zero, `ΔT` late. That is the causality signature
of lensing, as directly as it can be drawn.

### `caseA_diffraction_pattern.png` — the pattern in the source plane

`|F(w,y)|²` on the sky, at the two ends of the band. Concentric rings around
the lens, with the source marked at `y_A`. Both panels are windowed so the
fringes can be counted rather than merely be present: `|y| < 1.6` on the left,
and `|Δy| < 0.02` on the right, where `w = 778` makes the fringe period
~0.003 and a wider window packed a hundred rings into 500 pixels. Nothing is
lost by cropping — the pattern is axisymmetric, and the outer rings are more
of the same.

### `caseA_detector_envelope.png` — what a detector would record

The strain envelope, log scale, lensed against unlensed. Two peaks: the merger,
and 3.43 s later the echo, a factor ~4 lower (`√|μ₋| = 0.24`). No noise, no
antenna pattern, no PSD — an idealized view, and labelled as one.

---

## Case B — quasi-monochromatic source, moving lens

The same binary, earlier: `f_B = 0.05 Hz`, 240 days from merger. Now the lens
*moves*, sweeping `y(t)` through the pattern once per 4-day outer period.
`w_B = 0.31`, so this is diffraction, not geometric optics: no separate images,
and `F` is the exact wave-optics expression throughout.

### `caseB_repeated_pulses.png` — the lens alone

`|F(w_B, y(t))|²` over the 12-day observation: **three pulses, one per outer
period**, each peaking at 1.384 and flanked by diffraction ringing. The lower
panel zooms one pulse so the ringing is resolved — that ringing is the
signature of wave optics, and geometric optics does not produce it.

Each pulse lasts about **1.85 h** (FWHM), which is 1.9% of the period and 333
carrier cycles. This is the "lensing time": long compared to the lens response
time (0.985 s) by a factor of ~6700, which is why the quasi-static treatment
of a moving lens is legitimate.

### `caseB_detector_view.png` — the measured signal

The same three pulses, now in the strain: the envelope over 12 days on top,
and below a zoom on one pulse peak resolving individual 20 s carrier cycles.

The zoom shows an **amplitude** difference — crests 17.7% taller with the
lens. It does not show a dephasing, and a previous version of this figure
claimed it did. Both curves are evaluated at the same retarded time, so the
90-cycle Roemer delay is common to them and cancels; all that remains between
them is `F`, whose phase at the peak is `arg F = 1.01°`, a 0.056 s shift on a
20 s carrier — under a tenth of a pixel.

### `caseB_doppler_vs_lensing.png` — the orbit's two imprints

The point of this figure is that the outer orbit does **two** things to the
waveform, and the one this project is about is by far the smaller. In GW
cycles: Roemer/Doppler writes 90.5, lensing writes 0.017 — a ratio of 5206.

One curve per panel, on a shared time axis, each at its own scale. A ratio of
5206 cannot be drawn on a common axis: either the small curve is a flat line
at zero or the large one is off-screen. The previous version tried both at
once plus an inset, and the inset covered the sinusoid it was an inset of.

Conjunction is marked on both panels, because it answers the question the
figure otherwise raises: **why do the lensing pulses land on the Roemer
maximum?** Because the Roemer delay is extremal exactly where the
line-of-sight velocity vanishes, and that is the same instant the source
passes behind the lens. It is not a coincidence, and it is also the reason
evaluating `F` at a fixed `w_B` is safe here.

### `caseB_pattern_only.png` — the diffraction pattern

The fixed `|F(w_B,y)|²` pattern: concentric rings with the Einstein-ring peak
at the centre. The source's orbit track used to be drawn on top of it; that is
gone. The track is a near-straight line crossing rings whose spacing is the
whole content of the figure, and it hid exactly what it was meant to locate.
Where the source sits on the pattern is better said in time — that is
`caseB_repeated_pulses.png` — and the animated version in `report/report.html`
puts the two together properly, with a moving marker over this background.

---

## `theory/pulsos_esquema.pdf` — the schematic

A black-and-white, number-free version of the repeated-pulse train, for
`theory.pdf`. Same quantity as `caseB_repeated_pulses.png`, but the shape is
the content: one pulse per outer period, diffraction ringing on each flank,
baseline at `|F|² = 1`. `theory.pdf` states the phenomenon; the instance, with
this project's numbers, is the case figure.

Generated by `theory/pulsos_esquema.py`. `graphicx` is already loaded in
`theory.tex`, so including it is one line.

---

## Regenerating

Case B: `python3 cases/case_B_monochromatic/run.py` — runs anywhere, ~10 s.

Case A: `python3 cases/case_A_chirp/run.py` — **must be run where `pycbc` is
importable** (WSL2, per README), or it silently falls back to the
inspiral-only TaylorF2 waveform and the committed figures lose the merger and
ringdown entirely. `waveform_source` in `provenance/numbers.json` records
which path ran.
