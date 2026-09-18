# Case B: quasi-monochromatic source, moving lens

Produced by `run.py` (`python3 cases/case_B_monochromatic/run.py`, ~7 s);
every number below is in `provenance/numbers.json`. Same triple as
[Case A](../case_A_chirp/RESULTS.md), same inner binary, earlier in its
inspiral: f=0.05 Hz instead of the chirp's 10-125.6 Hz.

## Regime

Over the **3-outer-period (12 d)** observation used here, the GW frequency
drifts by only **1.9%** (0.0500 -> 0.0510 Hz) while there are **60.2** outer
periods left before merger — quasi-monochromatic, confirmed in
`tests/test_system.py::check_quasi_monochromatic_regime`. Unlike Case A, the
lens genuinely **moves**: y(t) sweeps through the diffraction pattern once
per 4-day outer period.

Three periods, not more, and the reason is the waveform rather than the
lens. Three is the fewest that shows the lensing is *repeated* rather than
a one-off. Beyond that the baseline only costs: Case B's waveform is the
leading-order quadrupole one (`chirp.restricted_pn_amplitude_td` for the
envelope, `chirp.freq_of_time` for the phase), and the phase a 0PN model
omits accumulates with the window while the phase the lens writes does not.
At f_B the source is at v/c = 0.030, so the 1PN correction is only 5.8e-3
of the leading term — but over 52 000 carrier cycles even that is large
compared with `arg F` = 0.018 rad. **This is why Case B is an amplitude
result.** |F| is untouched by phase truncation, so every number and figure
below that concerns amplitude stands on its own; the lensing *phase* is
reported (it is what makes the Doppler comparison below meaningful) but is
not claimed to be recoverable at this waveform accuracy.

w_B = **0.309** throughout (fixed — only y(t) varies) — solidly in the
diffraction-dominated regime (`F_hybrid` never switches away from the exact
closed form here; `w_geo_threshold=30`).

Holding w fixed is an approximation: neither the 1.9% intrinsic drift nor the
Doppler shift is fed back into F. D'Orazio & Loeb (2020) make the same
omission and give the reason — lensing only happens where the line-of-sight
velocity crosses zero, which for this near-edge-on circular orbit is exactly
the orbital phase of the lensing pulse, so the Doppler shift is at its
minimum precisely where F matters most. That is an argument, not a number,
so here are the numbers, and there are two: evaluating F at the true
instantaneous *observed* frequency instead of at w_B changes |F| by
`abs_F_relerr_from_fixed_w_at_pulse` = **1.0e-4** at the pulse (the argument
holds), but by up to `abs_F_relerr_from_fixed_w_max` = **3.5e-3** at the
worst point of the lensed half — far from the pulse, where |F| sits on the
steep flank of a diffraction fringe, on a value close to 1 that carries no
weight in any result quoted here. Quoting only the first would have been the
overclaim the argument invites.

## Repeated lensing

Exactly **half** the orbit is lensed at all (`fraction_of_time_lensed
=0.4994`): the source spends the other half in front of the lens
(`D_LS<0`, `geometry.impact_parameter_of_time`) where there is no lensing
geometry and F=1 by construction, not approximation; the impact parameter
formally diverges (`y->infinity`) exactly at each front/back crossing,
where `D_LS->0`.

The minimum impact parameter reached, **y=1.589**, coincides with Case A's
y_A to 4 significant figures — not a coincidence: both are the same
geometric minimum of the same (circular, e_out=0) outer orbit, found by an
`argmin` search over `geometry.impact_parameter_of_time` in each case
script's `main()` (same underlying function both times, so this is a
statement about grid resolution -- 600 samples per outer period here vs.
2000 in Case A -- not an independent confirmation of the geometry itself; a
previous version of this paragraph overstated it as "independently-written"
search code, corrected after an independent review). Kept here as an
internal consistency check between the two cases' scripts.

`caseB_repeated_pulses.png` shows the resulting **periodic amplification
pulses**, one per 4-day outer period, `|F|²` peaking at **1.384**
(`max_abs_F²`, from `max_abs_F=1.177`) with visible diffraction ringing on
each side of the main peak — the same ring structure as Case A's Figure 3,
just swept through in time instead of shown as a static map. Averaged over
the whole 12-day observation, `<|F|²>=1.006` and the RMS amplitude
amplification is **1.003**: close to unity, because the lensed half and the
completely-unlensed half (F=1) roughly average out — most of the
observation is spent either unlensed or only mildly amplified, with sharp
but brief peaks.

## The outer orbit's *other* imprint — and it is the larger one

Sweeping y(t) through the diffraction pattern is only half of what the
outer orbit does to the waveform, and it is the smaller half.

The same orbit also carries the source along the line of sight, at
`beta_los_max` = **0.0165**, with a light-crossing time `a_out/c`. The
resulting Roemer (light-travel-time) modulation of the arrival time is
`roemer_delay_ptp_s` = **1810 s peak-to-peak**, which at f_B is
`roemer_phase_ptp_cycles` = **90.5 GW cycles** of phase written on the
waveform. Set against that:

| | |
|---|---|
| Roemer/Doppler phase, peak-to-peak | **90.5 cycles** (`roemer_phase_ptp_cycles`) |
| lensing phase `arg F/2π`, peak-to-peak | **0.0174 cycles** (`lens_phase_ptp_cycles`) |
| ratio | **5206×** (`roemer_over_lens_phase_ratio`) |
| Roemer delay vs. image delay `ΔT` at the pulse | 1810 s vs **3.43 s** → **527×** (`roemer_over_image_delay_ratio`) |

`caseB_doppler_vs_lensing.png` is that comparison as one figure, one
mechanism per panel and each at its own scale, since a ratio of 5206 cannot
be drawn on a common axis. The top panel carries the Roemer **delay** in GW
cycles *and* the Doppler **shift** on a twin axis, because they are not the
same thing and sit a quarter period apart — see below. The bottom panel is
the lensing phase on its own dense grid, with the fixed-`w` systematic drawn
as a band so that what can be read off it is separated from what cannot.
This is not a defect of the lensing calculation — it is
what a hierarchical triple actually looks like, and it is why
[D'Orazio & Loeb (2020)](https://arxiv.org/abs/1910.02966), the paper
this case follows for repeated lensing, models the Doppler boost alongside
it (their Appendix A) rather than either one alone.

The Doppler swing in observed frequency, `doppler_freq_swing_Hz` =
**0.0025 Hz**, is itself larger than the entire intrinsic chirp drift across
the 12-day observation (`intrinsic_freq_drift_Hz` = **0.00097 Hz**).

**Delay and shift are not the same thing, and they are a quarter period
apart.** They are the same quantity separated by a derivative — the delay
goes as the separation along the line of sight, the shift as its velocity —
so where one is extremal the other vanishes. Measured, at conjunction: the
Roemer delay is at its extremum, **+905 s of ±905**, while `|v_los|/c` is
**1.3e-6** against an orbital maximum of **1.6e-2**. A quarter period away it
is exactly the other way round, and that is where `D_LS = 0` and there is no
lensing at all. So the lensing pulse coincides with the *delay*, not with the
shift. Calling the pair "the Roemer-Doppler effect", as though it were one
thing with one maximum, puts it at the wrong orbital phase; an earlier version
of the figure's y-label did exactly that.

That quadrature is also why evaluating `F` at the fixed nominal `w_B` is
legitimate: the frequency shift is smallest precisely where `F` matters most.
That is an argument, so here are the numbers, and there are four —
`abs_F_relerr_from_fixed_w_at_pulse` = **1.1e-4** and
`abs_F_relerr_from_fixed_w_max` = **3.5e-3** for the modulus,
`arg_F_cycles_err_from_fixed_w_at_pulse` = **2.5e-5** cycles and
`arg_F_cycles_err_from_fixed_w_max` = **5.6e-4** cycles for the phase. Both
worst cases fall far from the pulse, on the steep flank of a fringe where
`|F|` is near 1 and nothing quoted here is read off. The phase bound is new:
`|F|` had one, `arg F` did not, and it is what decides how much of the lensing
panel's ringing is a result — the two large excursions are 21 times it, the
far ringing sits inside it.

**How it is applied, and what it does not change.** The waveform is
evaluated at the retarded emission time solved from
`t_obs = t_em + z_src(t_em)/c` (`src/gwlens/doppler.py::emission_time`,
converging in `emission_time_iterations` = 8 fixed-point steps). That single
substitution *is* the complete first-order treatment — differentiating the
relation gives `dt_em/dt_obs = 1/(1+beta_los)` identically, so no separate
Doppler factor is multiplied in anywhere, which is checked to 5.5e-11 in
`tests/test_doppler.py::check_retarded_time_is_the_whole_first_order_doppler`.

Because it is pure *timing*, it leaves the strain **envelope** alone: the
envelope evolves on the 240-day time-to-merger scale and the delay only
reshuffles it by ~900 s, for a fractional change of
`doppler_envelope_max_fractional_change` = **1.1e-5**. Every amplitude
result in this file — the pulse heights, `<|F|²>`, the RMS amplification,
every diffraction figure — is therefore numerically unchanged by adding
Doppler, and demonstrably so rather than by assumption
(`tests/test_doppler.py::check_envelope_unaffected_by_roemer`). What *is*
changed is the waveform's phase, by
`doppler_phase_ptp_cycles_vs_static` = **91.9 cycles** measured directly
against the same waveform built without the substitution.

**What is deliberately not applied.** Everything beyond first order is,
for a circular orbit, exactly constant: the gravitational redshift in the
lens's potential plus the transverse (second-order) Doppler shift combine
into the Schwarzschild circular-geodesic factor
`(1 - 3GM_L/(a c²))^(-1/2)`, giving
`orbital_redshift_factor_minus_1` = **4.08e-4**. A constant redshift is
exactly degenerate with a rescaling of the chirp mass and of t_c, so it is
a bias on any inferred M_chirp, not a feature of any waveform — reported,
not applied (`src/gwlens/doppler.py::orbital_redshift_factor`).

## The extended pattern, and how a detector would see it

`caseB_pattern_only.png` is the fixed |F(w_B,y)|² diffraction pattern:
concentric rings with the Einstein-ring peak at the centre, the same
structure Case A shows at its own frequency. An earlier version of this
figure drew the source's sky-projected orbit track on top of it. That is
gone: the track is a near-straight line crossing rings whose spacing is
the entire content of the figure, and drawing it over them hid exactly
what it was supposed to locate. Where the source sits on the pattern is
better said in time, which is `caseB_repeated_pulses.png`; the animated
version in `report/report.html` puts the two together properly, with a
moving marker over this same background.

`caseB_detector_view.png` is the idealized-detector view asked for: top
panel, the strain envelope over the full 12-day observation, showing three
clean amplification pulses; bottom panel, **five** cycles of the 20 s carrier
at one pulse's peak, on its own fine grid (the 1800-point envelope grid, 576 s
spacing, undersamples a 20 s carrier by a factor of 30).

The baseline in the top panel **drifts upward by 1.29%** across the
observation, and that is the source rather than the lens: quasi-monochromatic
is not monochromatic. The binary is still inspiralling, `f_start_Hz`=0.05000 to
`f_end_Hz`=0.05097 (`fractional_freq_drift` = 1.9%), and the amplitude goes as
`f^(2/3)`, which predicts `(t_c/(t_c-T))^(1/4)` = 1.0129 — matching the drift
to five digits. It is the only place in this case where the "quasi" is visible.

What the bottom panel shows is the **amplitude** difference: crests 17.7%
taller with the lens than without (`max_abs_F` = 1.1766). The window was 30
cycles until it was noticed that the figure looked like it was *understating*
the amplification — at that density the two curves cross sixty times and read
as one sinusoid drawn twice. Nothing was given up by narrowing, because `|F|`
was flat either way: it moves 0.064% of its peak across ±300 s and 0.002%
across ±50 s.

It does *not* show a dephasing — and, to be exact about it, it **could not**,
even if the effect were larger. Both curves are evaluated at the same retarded
time, so the 90-cycle Roemer delay is common to them and cancels; what is left
is `F`, whose phase runs over `lens_phase_ptp_cycles` = 0.0174 cycles across
the *whole* observation, sub-pixel at any window that resolves cycles. Crests
that line up here are therefore not evidence of anything. The evidence that
the lens writes amplitude and not phase is that number against Roemer's 90.5.
Two earlier versions of this file claimed the panel demonstrated it.

## One outer turn, at any magnification (`caseB_one_period.html`)

The only figure in this case that is not an image, and for an arithmetic
reason: one outer turn holds **17280 carrier cycles**
(`one_period_carrier_cycles`), so no fixed picture shows the turn and the wave
at once. It is a page — wheel to zoom, drag to pan, a button per stage — and it
makes the quadrature above something to click rather than read:

| stage | `|F|` | Roemer delay | `f_obs/f_em − 1` | drag vs a fixed-period comb |
|---|---|---|---|---|
| at the pulse (conjunction) | **1.1762** | **+905 s** | −0.02 % | −0.001 cycles |
| at quadrature | 1.0000 | 4 s | **−1.62 %** | **−0.089 cycles** |
| half a turn on (source in front) | 1.0000 | **−905 s** | −0.03 % | +0.015 cycles |

The last column answers the question the page otherwise invites. A comb of
fixed period, anchored to the window centre, measures `dτ/dt` and not `τ`: a
constant delay shifts the wave rigidly and a rigid shift is invisible against
a comb that re-anchors with it, so at conjunction — where the delay is at its
*largest* — the wave stays locked, and at quadrature, where the delay passes
through zero but is changing fastest, it walks off. Comparing a crest's
position *between* stages reads nothing at all: two windows a quarter turn
apart are separated by thousands of accumulated cycles, of which Roemer
contributes about 45.

The signal is generated in the page rather than sampled into it — sampling
`h(t)` at 11 points per cycle would be 200000 numbers per turn — so `F(t)` and
the observed-to-emission time map are tabulated every 69 s and the carrier
comes from closed forms. That those agree with the phase this repository
computes is measured:
`one_period_closed_form_phase_err_cycles` = **5.8e-10** cycles over the whole
turn. The first version had the prefactor wrong by a factor of two and that
number is what caught it, reading 4407.
