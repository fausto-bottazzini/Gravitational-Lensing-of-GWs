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
merger because `df/dt` does: **0.85 fringes per second at 10 Hz against 93
per second at 36 Hz**. The fringe period is constant in *frequency* (0.29 Hz,
`caseA_F_of_f.png`); this is the same curve seen through `f(t)`.

A previous version plotted the difference instead. That is a different
statement — `h_lensed − h_unlensed = (√μ₊−1)h(t) + √|μ₋|h(t−ΔT)`, so
subtracting isolates the second image as a delayed copy — which is about
arrival times, not about interference. The second image is already visible in
the second panel and is a whole second track in `caseA_spectrogram.png`.

**The fifth panel is the same quantity evaluated in closed form**, `|F(f(t))|`
with `f(t)` the leading-order chirp. The obvious question about the fourth
panel is whether the measured ratio simply *is* `|F|`. It is, in the middle,
and where it is not, one quantity explains both ends at once.

In geometric optics the lensed signal is two copies, `h_l(t) = √μ₊·h(t) +
√|μ₋|·h(t−ΔT)` up to the Morse phase, and the analytic signal is linear, so
the ratio of envelopes is *exactly*

    env_l/env_u = |√μ₊ + √|μ₋|·q(t)·e^{iΔφ(t)}|,   q(t) = A(t−ΔT)/A(t)

The modulation **depth** is `√|μ₋|·q(t)`, not `√|μ₋|`. The closed form is the
`q = 1` limit: in the frequency domain both images contribute at the same `f`
with the same `|H(f)|`, which in the time domain means the chirp amplitude must
barely change over `ΔT`.

`q → 0` at both ends, the same mechanism seen twice:

- **Before `t − t_merger = −11.75 s`** the delayed copy does not exist yet. The
  signal starts 15.18 s before the merger and the copy enters `ΔT = 3.43 s`
  later, so there is nothing to interfere with and the ratio sits flat at
  `√μ₊`: between **1.0041 and 1.0342** there, against `√μ₊ = 1.0283`, opening to
  the full **0.7895…1.2640** band from −11.5 s on.
- **Towards the merger** `A(t)` runs away while `A(t−ΔT)` is still inspiral, so
  `q` falls — **0.91 at −7 s, 0.83 at −2.7 s, 0.72 at −1.1 s, 0.41 at −0.07 s**
  — and the band closes back towards `√μ₊ = 1.028`, which is very nearly 1.

Both ends read as "the lens does nothing" for the same reason: at that instant
there is effectively only one image. The closed form knows nothing about either
— `|F(f)|` is a frequency-domain statement about the whole signal — so it keeps
filling the full band, and the two panels part company.

This is measured rather than asserted: `ratio_depth_model_max_relerr` compares
the predicted half-depth `√|μ₋|·q(t)` against the measured extremes over windows
that each hold several fringes, and it agrees to **1.3 %** from −8 s to 60 ms
before the merger.

So the fourth panel carries `√μ₊ ± √|μ₋|·q(t)` drawn over the measured band,
which fills it from end to end. That is the whole explanation and it needs no
words on the figure.

Two earlier versions got this wrong in ways worth recording. The first
truncated the panel at −3.22 s. The second ran it to the merger but shaded the
two ends as regions where the measurement "could not be trusted", and blamed
the closing-up on the Hilbert transform being unable to follow fringes faster
than its own carrier — which reads as a numerical problem, and is not what
happens. The envelope follows the fringes fine; there is simply less modulation
left to follow.

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

Each pulse lasts about **1.85 h** (FWHM of `|F|`; the panel plots `|F|²`,
whose FWHM is a little narrower at 1.81 h), which is 1.9% of the period and 333
carrier cycles. This is the "lensing time": long compared to the lens response
time (0.985 s) by a factor of ~6700, which is why the quasi-static treatment
of a moving lens is legitimate.

### `caseB_detector_view.png` — the measured signal

The same three pulses, now in the strain: the envelope over 12 days on top,
and below a ±50 s window at one pulse's peak — **five carrier cycles**. The
lensed crests are **17.7% taller** (`max_abs_F = 1.1766`). That is the panel.

**The baseline in the top panel drifts upward**, and that is the source, not
the lens. "Quasi-monochromatic" is not monochromatic: over the 12 days the
binary is still inspiralling, from 0.05000 to 0.05097 Hz (+1.94%), and the
amplitude goes as `f^(2/3)`, so it rises **1.287%**. That matches
`(t_c/(t_c−T))^(1/4) = 1.01287` to five digits. It is the only place in Case B
where the "quasi" is visible.

The window was ±300 s (30 cycles) until it was noticed that the figure looked
like it was *understating* the amplification: at that density the two curves
cross sixty times and read as one sinusoid drawn twice. Nothing is given up by
narrowing, because `|F|` was flat either way — it moves 0.064% of its peak
across ±300 s and 0.002% across ±50 s. Going the other way, to where `|F|`
really does move (±2500 s, 4.4%), drops the cycles to 5 px and the carrier
fills in as a solid block; drawing envelopes over that just adds four
horizontal rules. The pulse is 6653 s wide (FWHM) against a 20 s carrier —
**333 cycles per pulse** — so no window shows both, and the pulse shape belongs
to `caseB_repeated_pulses.png`.

What separates the two curves is **amplitude**, and the panel cannot show a
dephasing — not because there is none to see, but because there could not be
one at this scale even if the effect were larger. Both curves are evaluated at
the same retarded time, so the 90-cycle Roemer delay is common to them and
cancels; what is left is `F`, whose phase runs over `lens_phase_ptp_cycles =
0.0174` cycles across the *whole* observation, sub-pixel here by construction.
Crests that line up in this panel are therefore not evidence of anything. The
evidence that the lens writes amplitude and not phase is that number against
Roemer's 90.5, and its figure is `caseB_doppler_vs_lensing.png`. An earlier
title claimed the dephasing was visible cycle by cycle; it is not.

### `caseB_doppler_vs_lensing.png` — the orbit's two imprints

The point of this figure is that the outer orbit does **two** things to the
waveform, and the one this project is about is by far the smaller. In GW
cycles: Roemer/Doppler writes 90.5, lensing writes 0.017 — a ratio of 5206.

One mechanism per panel, on a shared time axis, each at its own scale. A ratio
of 5206 cannot be drawn on a common axis: either the small curve is a flat line
at zero or the large one is off-screen. An earlier version tried both at once
plus an inset, and the inset covered the sinusoid it was an inset of.

**Delay and shift are not the same thing, and they sit a quarter period apart.**
The top panel's y-label used to read "Roemer / Doppler", which conflates them
and invites a correct objection: the Doppler *shift* is largest at quadrature,
with the source running along the line of sight, while the lensing happens at
conjunction — so marking conjunction looked like a claim that the two peak
together. They do not. What peaks at conjunction is the *delay*, and the two
are a derivative apart, so where one is extremal the other vanishes. Measured,
at conjunction: the Roemer delay is **+904.9 s**, its extremum over the orbit
(±904.9), and `|v_los|/c` is **1.3×10⁻⁶** against an orbital maximum of
**1.645×10⁻²**. A quarter period away it is exactly the other way round. The
panel now draws both, on twin axes, and the crossing pattern says it without a
caption. It is also why evaluating `F` at a fixed `w_B` is safe: the frequency
shift is at its minimum precisely where `F` matters most.

**What can be read off the lower panel, and what cannot.** Its dense ringing is
real — a grid 38× finer reproduces it, 0.01738 cycles peak to peak against
0.01733 on the old grid — but not all of it is above the systematics, and the
panel now says so. Evaluating `F` at fixed `w_B` rather than at the true
instantaneous observed `w` moves `arg F` by up to
`arg_F_cycles_err_from_fixed_w_max = 5.6×10⁻⁴` cycles, drawn as the grey band.
Against that:

- the two large excursions (**0.0121 cycles**, at `y ≈ 3`, ±1.28 h either side
  of conjunction) are **21× the systematic** and are a result;
- the low ringing far from the pulses sits *inside* the band and is not. Out
  there the lensing has already gone to zero while the error has not — which is
  precisely why it never mattered that the approximations loosen away from
  alignment.

The lower panel gets its own grid for the same reason `caseA_F_of_f.png` does:
`F(w,y)` is a closed form, so the drawing grid need not inherit the waveform's
sampling. The shared 576 s grid resolved the ringing but swallowed one
oscillation (30 zero crossings against 32).

**Why the phase's shape differs from the amplitude's.** `caseB_repeated_pulses.png`
shows a single upward peak at conjunction; this panel shows two deep troughs
*flanking* conjunction with a local maximum between them. They are the same
diffraction ringing seen in the two components of one complex `F` — over
`y ∈ [1.55, 30]` both cross their asymptote **45 times** — but in quadrature
with each other. At closest approach `|F|² = 1.384`, its maximum, while
`arg F/2π = +0.0028`, essentially zero; `arg F`'s extremum (−0.0121) falls at
`y ≈ 3`, which the source crosses twice, once on each side. That is the general
rule: where `|F|` is stationary, `F` is nearly real, so its phase passes through
zero there.

### `caseB_one_period.html` — one outer turn, at any magnification

The only figure here that is not an image, and for a reason that is arithmetic
rather than taste: **one outer turn holds more than 17 000 carrier cycles**. No
fixed picture shows the turn and the wave at once — seeing both means being able
to zoom four orders of magnitude — so this one is a page. Wheel to zoom on the
pointer, drag to pan, click the overview to jump, or use the buttons, which go
straight to each stage:

| stage | `|F|` | Roemer delay | `f_obs/f_em − 1` |
|---|---|---|---|
| carrier at the pulse (conjunction) | **1.1762** (max) | **+905 s** (max) | −0.02 % (≈0) |
| carrier at quadrature | 1.0000 | 4 s (≈0) | **−1.62 %** (max) |
| carrier half a turn on (source in front) | 1.0000 | **−905 s** (min) | −0.03 % (≈0) |

That table is the whole of §"the orbit does two things" made clickable: lensing
and Doppler peak a quarter period apart, and what coincides with the pulse is
the *delay*, not the shift.

**The signal is not sampled into the page, it is generated in it.** What varies
slowly — `F(t)` and the observed-to-emission time map — is tabulated every 69 s;
what varies fast has a closed form,

    f(t) = f_c (1 − t/τ_c)^(−3/8),   φ(t) = (16π/5) f_c τ_c [1 − (1 − t/τ_c)^(5/8)]

with `A ∝ f^(2/3)`. Sampling `h(t)` instead would be 200 000 numbers per turn.
That the closed form agrees with the phase the repository computes is measured,
not assumed: `one_period_closed_form_phase_err_cycles` = **5.8×10⁻¹⁰ cycles**
over the whole turn. The first version of it was wrong by a factor of two in the
prefactor and that check is what caught it — it read 4407 cycles.

**The vertical scale is fixed for the whole turn, not per window.** An earlier
version auto-scaled each window to its own maximum, so panning off the pulse
made the wave grow — the axis meant something different at every position,
which destroys the one comparison the page exists to allow.

**A fixed-period comb appears when few cycles are in view**, anchored to the
window centre. It is the clock against which the Doppler *drag* is visible: at
conjunction and half a turn on the shift vanishes and the wave stays locked to
the comb (−0.001 and +0.015 cycles across the window); at quadrature it walks
off it by **−0.089 cycles**, which is the −1.62 % shift plus the 0.16 % of
intrinsic drift accumulated over the quarter turn. Without a reference the
stretching is invisible.

It is worth saying what the comb is *not* for. Comparing the position of a
crest *between* stages reads nothing: two windows a quarter turn apart are
separated by thousands of accumulated cycles, of which Roemer contributes about
45, so what is left mod 2π is dominated by the carrier. The drag is legible only
*within* a window, against the comb.

Drawing follows one rule: whether more than one carrier cycle falls in a pixel
column. If it does, the curve would fill the column and the only honest thing is
the band, which is then drawn **exactly** as ±|F|·A rather than sampled; if it
does not, the wave itself is drawn. The page is its own source and its own
output — `run.py` replaces the contents of its `<script id="datos">` block in
place, the same arrangement as `report/report.html`, and for the same reason:
the page is meant to open as a `file://` URL, where `fetch()` of a local JSON is
blocked by CORS.

It replaces a three-panel PNG of the same turn. That figure was correct and
unusable: the pulse it exists to show occupies 1.9 % of the axis.

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
