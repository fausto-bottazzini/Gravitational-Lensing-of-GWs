# The figures, one by one

What each figure is for, what it shows, and what it deliberately does not.
Nine figures: four per case, plus one schematic in `theory/`. The figures
themselves are in Spanish, because they go into `report.pdf`; this page is in
English like the rest of the repo's own documentation.

Panel titles in the figures are short labels, on purpose. A title names the
panel; it is not a caption. Everything a reader needs beyond the label is
here.

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

### `caseA_F_of_f.png` — where the fringes come from

Two panels: `F` in the complex plane, and `arg F(f)`.

Across the whole Case A band this is geometric optics, so
`F = √μ₊ − i√|μ₋| e^{iwΔT}`: sweeping `f` walks the point around a **circle**
of radius `√|μ₋| = 0.2396` centred on `√μ₊ = 1.0283`, one turn per fringe.
Three facts come off that one picture:

- `|F|` lies between the nearest and farthest points of the circle from the
  origin, `0.789` and `1.268`;
- `arg F` never winds, because the circle does not enclose the origin. The two
  tangents from the origin (dashed) bound it by
  `arcsin(√|μ₋|/√μ₊) = 0.2352 rad`, which is exactly what it reaches;
- the fringe period is `1/ΔT = 0.2912 Hz` anywhere in the band, since `wΔT` is
  linear in `f` and `ΔT` is fixed — the same 3.434 s as the echo.

It replaces three panels of `|F(f)|` and `arg F(f)` in different windows,
which measurement showed to be the same curve three times: the fringe period
is 0.29120 Hz at 10 Hz, at 60 and at 123, and the `|F|` and `arg F` bounds
agree to four decimals in all three. The circle is drawn over exactly one turn
— the 3 Hz window holds ten fringes, so drawing the whole locus laid ten
circles on top of each other and the chords between points nearly a turn apart
cut straight across it. Both panels use a dense grid of their own rather than
the waveform's FFT grid, whose 0.031 Hz spacing is nine samples per fringe.

### `caseA_strain_time.png` — what the lens does to the waveform

Five panels: unlensed strain, lensed strain **on the same scale**, a zoom
centred on the merger, the **measured ratio** of the two envelopes, and that
same ratio **in closed form** directly below it.

The first two are separate panels rather than one overlay. At this time
resolution the two curves are indistinguishable — which is the point, and is
worth seeing — and drawing one on top of the other simply hid the lower one.

The fourth panel is where the lensing becomes visible. It is the **ratio**,
not the difference, and the distinction matters. Dividing the lensed analytic
signal by the unlensed one leaves `|F|`, and because the chirp sweeps
frequency monotonically, plotting that against time sweeps out F's
interference fringes. The curve oscillates between `√μ₊ ∓ √|μ₋|` (drawn as
dashed lines) from end to end, and the fringes crowd together towards the
merger because `df/dt` does: **1.0 fringes per second at 10 Hz against 93 per
second at 36 Hz**. The fringe period is constant in *frequency* (0.29 Hz,
`caseA_F_of_f.png`); this is the same curve seen through `f(t)`.

A previous version plotted the difference instead. That is a different
statement — `h_lensed − h_unlensed = (√μ₊−1)h(t) + √|μ₋|h(t−ΔT)`, so
subtracting isolates the second image as a delayed copy — which is about
arrival times, not about interference. The second image is already visible in
the second panel and is a whole second track in `caseA_spectrogram.png`.

**The fifth panel is the same quantity evaluated in closed form**, `|F(f(t))|`
with `f(t)` the leading-order chirp. The obvious question about the fourth
panel is whether the measured ratio simply *is* `|F|`; the answer is that it
is, in the middle, and the two shaded regions are where it cannot be:

- **Before `t − t_merger = −11.75 s`** there is only one image. The second is
  the first delayed by `ΔT = 3.43 s`, the signal starts 15.18 s before the
  merger, and until the copy arrives there is nothing to interfere with. The
  measured ratio sits flat at `√μ₊`: between **1.0041 and 1.0342** there,
  against `√μ₊ = 1.0283`, and it opens to the full **0.7895…1.2640** band from
  −11.5 s on. The closed form knows nothing about this — `|F(f)|` is a
  frequency-domain statement about the whole signal — so it oscillates there
  anyway, and the two panels disagree.
- **After `t − t_merger = −3.22 s`** the envelope can no longer follow the
  fringes. `env_l` and `env_u` are analytic-signal magnitudes, and an analytic
  envelope only tracks modulation *slower* than its own carrier. The figure of
  merit is `ΔT·(df/dt)/f`: 0.10 at 13 s before the merger, 0.43 at 3 s, 0.86 at
  1.5 s, 2.17 at 0.6 s. Past ~0.25 the drawn fringes shrink toward `√μ₊`, and
  that shrinking is a property of the Hilbert transform, not of the lens. The
  closed form keeps filling the band, which is what really happens.

An earlier version simply truncated the fourth panel at −3.22 s. Showing the
boundary is better than hiding it, and the closed-form panel underneath is
what makes it legible.

Neither panel runs past the merger, and that is not a choice: after it there is
no denominator. The unlensed signal dies in about 20 ms — `env_u` is 1.6e-5 of
its peak 50 ms later — and the quotient runs to thousands. That is no longer
the ratio of two signals, it is the second image divided by numerical noise.

Both are drawn as a **min–max band per pixel column**, not as sampled curves:
near the merger one column spans dozens of fringes and point sampling there
would be pure alias. For the closed form the per-column extremes are exact
(the fringe phase `w·ΔT̂ = 2πfΔT` is linear in `f`, so the extremes of the sine
over a column follow in closed form); for the measured one they are the
extremes of the samples that fall in the column.

### `caseA_diffraction_pattern.png` — the pattern in the source plane

`|F(w,y)|²` on the sky, at the two ends of the band. Concentric rings around
the lens, with the source marked at `y_A`. Both panels are windowed so the
fringes can be counted rather than merely be present: `|y| < 2.0` on the left,
and `|Δy| < 0.02` on the right, where `w = 778` makes the fringe period
~0.003 and a wider window packed a hundred rings into 500 pixels. Nothing is
lost by cropping — the pattern is axisymmetric, and the outer rings are more
of the same.

### `caseA_spectrogram.png` — the chirp in time-frequency

The chirp drawn the way a chirp is normally drawn: time across, frequency up
(log), colour for amplitude. Unlensed above, lensed below, one colour scale.

The unlensed panel is one track sweeping up to the merger. The lensed panel
has **two** — the same sweep repeated `ΔT = 3.43 s` later at `√|μ₋| = 0.24`
of the amplitude. Two images become two tracks. That is the whole of
geometric-optics lensing in one picture, and it is the figure to show first.

The stipple along the lensed track at low frequency is the interference
between the two images — the same fringes `caseA_strain_time.png`'s ratio
panel measures, here seen as a texture rather than a curve.

Window: 0.5 s (`nperseg=1024` at `fs=2048`), which is five cycles at the
10 Hz band start — the shortest that still localises the low end, and short
enough that the 0.29 Hz fringes do not resolve into the picture. Hann rather
than scipy's default Tukey(0.25): its sidelobes are ~30 dB lower, and the
Tukey ones were drawing a fan of spurious arcs above the real track.

It replaces an amplitude-vs-time plot of the same two waveforms, which added
nothing to the strain figure. No noise, no antenna pattern, no PSD.

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
and below a ±300 s window at one pulse's peak, resolving individual 20 s
carrier cycles. The lensed crests are **17.7% taller**. That is the panel.

It cannot also show the pulse's shape, and the arithmetic is why. The pulse
is 6653 s wide (FWHM) against a 20 s carrier: **333 cycles per pulse**. Across
the figure's ~1360 usable pixels, ±300 s gives 45 px per cycle and `|F|` moves
0.4% of its peak height — clean cycles, flat amplification. Widening until
`|F|` visibly rises and falls (±2500 s, 29%) drops the cycles to 5 px and the
carrier fills in as a solid block; drawing envelopes over it then just adds
four horizontal rules. Both were tried. Neither works, because no window does
both — the pulse shape belongs to `caseB_repeated_pulses.png`.

What separates the two curves is **amplitude**. It is not a dephasing, and a
previous version of this figure claimed it was. Both are evaluated at the same
retarded time, so the 90-cycle Roemer delay is common to them and cancels; all
that remains between them is `F`, whose phase at the peak is `arg F = 1.01°`,
a 0.056 s shift on a 20 s carrier — under a tenth of a pixel.

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

### `caseB_pattern.png` — the diffraction pattern

The fixed `|F(w_B,y)|²` pattern: concentric rings with the Einstein-ring peak
at the centre, axes, colour bar, and a single marker at the closest approach
the source makes, `y = 1.59`.

The source's whole orbit track used to be drawn over it. That is gone: the
track is a near-straight line crossing rings whose spacing is the entire
content of the figure, and it hid exactly what it was meant to locate. The
marker is the one thing the track carried that the pattern does not say by
itself. Where the source sits on the pattern over time is better said in time,
which is `caseB_repeated_pulses.png`.

`caseB_pattern_only.png` is the same pattern with no axes, labels or colour
bar, on a square canvas. It is not a figure for the report: it is the
background that `report/report.html` animates a moving marker over.

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
