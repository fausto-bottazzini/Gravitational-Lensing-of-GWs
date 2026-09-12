# Case B: quasi-monochromatic source, moving lens

Produced by `run.py` (`python3 cases/case_B_monochromatic/run.py`, ~7 s);
every number below is in `provenance/numbers.json`. Same triple as
[Case A](../case_A_chirp/RESULTS.md), same inner binary, earlier in its
inspiral: f=0.05 Hz instead of the chirp's 10-125.6 Hz.

## Regime

Over the **6-outer-period (24 d)** observation used here, the GW frequency
drifts by only **4.0%** (0.0500 -> 0.0520 Hz) while there are **60.2** outer
periods left before merger — quasi-monochromatic, confirmed in
`tests/test_system.py::check_quasi_monochromatic_regime`. Unlike Case A, the
lens genuinely **moves**: y(t) sweeps through the diffraction pattern once
per 4-day outer period.

w_B = **0.309** throughout (fixed — only y(t) varies) — solidly in the
diffraction-dominated regime (`F_hybrid` never switches away from the exact
closed form here; `w_geo_threshold=30`).

Holding w fixed is an approximation: neither the 4% intrinsic drift nor the
Doppler shift is fed back into F. D'Orazio & Loeb (2020) make the same
omission and give the reason — lensing only happens where the line-of-sight
velocity crosses zero, which for this near-edge-on circular orbit is exactly
the orbital phase of the lensing pulse, so the Doppler shift is at its
minimum precisely where F matters most. That is an argument, not a number,
so here are the numbers, and there are two: evaluating F at the true
instantaneous *observed* frequency instead of at w_B changes |F| by
`abs_F_relerr_from_fixed_w_at_pulse` = **1.0e-4** at the pulse (the argument
holds), but by up to `abs_F_relerr_from_fixed_w_max` = **6.4e-3** at the
worst point of the lensed half — far from the pulse, where |F| sits on the
steep flank of a diffraction fringe, on a value close to 1 that carries no
weight in any result quoted here. Quoting only the first would have been the
overclaim the argument invites.

## Repeated lensing

Exactly **half** the orbit is lensed at all (`fraction_of_time_lensed
=0.4997`): the source spends the other half in front of the lens
(`D_LS<0`, `geometry.impact_parameter_of_time`) where there is no lensing
geometry and F=1 by construction, not approximation; the impact parameter
formally diverges (`y->infinity`) exactly at each front/back crossing,
where `D_LS->0`. `caseB_pattern_with_orbit.png` (below) shows the lensed
and unlensed halves of the orbit directly, distinguished by color.

The minimum impact parameter reached, **y=1.589**, coincides with Case A's
y_A to 4 significant figures — not a coincidence: both are the same
geometric minimum of the same (circular, e_out=0) outer orbit, found by an
`argmin` search over `geometry.impact_parameter_of_time` in each case
script's `main()` (same underlying function both times, so this is a
statement about grid resolution -- 2000 samples/period here vs. 3600 in
Case A -- not an independent confirmation of the geometry itself; a
previous version of this paragraph overstated it as "independently-written"
search code, corrected after an independent review). Kept here as an
internal consistency check between the two cases' scripts.

`caseB_repeated_pulses.png` shows the resulting **periodic amplification
pulses**, one per 4-day outer period, `|F|²` peaking at **1.385**
(`max_abs_F²`, from `max_abs_F=1.177`) with visible diffraction ringing on
each side of the main peak — the same ring structure as Case A's Figure 3,
just swept through in time instead of shown as a static map. Averaged over
the whole 24-day observation, `<|F|²>=1.006` and the RMS amplitude
amplification is **1.003**: close to unity, because the lensed half and the
completely-unlensed half (F=1) roughly average out — most of the
observation is spent either unlensed or only mildly amplified, with sharp
but brief peaks.

## The outer orbit's *other* imprint — and it is the larger one

This section corrects the headline this case used to carry. Until
2026-09-07 this script modelled only what the outer orbit does *optically*
— sweep y(t) through the diffraction pattern — and described that as the
outer orbit's signature on the waveform. It is not; it is the smaller half
of the answer.

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

`caseB_doppler_vs_lensing.png` is that comparison as one figure: the
Roemer sinusoid at full scale, the lensing phase correctly flat on zero
beside it, and an inset at 3392× so the lensing phase is still shown
resolved (its diffraction ringing around each pulse) rather than only shown
to be invisible. This is not a defect of the lensing calculation — it is
what a hierarchical triple actually looks like, and it is why
[D'Orazio & Loeb (2020)](https://arxiv.org/abs/1910.02966), the paper
this case follows for repeated lensing, models the Doppler boost alongside
it (their Appendix A) rather than either one alone.

The Doppler swing in observed frequency, `doppler_freq_swing_Hz` =
**0.0035 Hz**, is itself larger than the entire intrinsic chirp drift across
the 24-day observation (`intrinsic_freq_drift_Hz` = **0.0020 Hz**).

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
`doppler_envelope_max_fractional_change` = **1.2e-5**. Every amplitude
result in this file — the pulse heights, `<|F|²>`, the RMS amplification,
every diffraction figure — is therefore numerically unchanged by adding
Doppler, and demonstrably so rather than by assumption
(`tests/test_doppler.py::check_envelope_unaffected_by_roemer`). What *is*
changed is the waveform's phase, by
`doppler_phase_ptp_cycles_vs_static` = **93.8 cycles** measured directly
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

`caseB_pattern_with_orbit.png` puts Case A's ring pattern and Case B's orbit
together explicitly: the fixed |F(w_B,y)|² diffraction pattern (concentric
rings, central Einstein-ring spike) with the source's actual sky-projected
track overlaid, split into its lensed (cyan) and unlensed (dashed green)
halves — this is the "extended, in-space" pattern and the trajectory a real
triple would trace through it, in one figure.

`caseB_detector_view.png` is the idealized-detector view asked for: top
panel, the strain envelope over the full 24-day observation, showing six
clean amplification pulses riding on the slowly-drifting quasi-monochromatic
carrier; bottom panel, a dedicated fine time grid (dt=0.1 s, since the
coarse 3600-point grid used for the envelope — 576 s spacing — wildly
undersamples the 20 s carrier period; using it directly for a "raw
waveform" panel was an aliasing bug caught by eye and fixed, see
wiki/log.md) zoomed on one pulse peak, where the lensed and unlensed
waveforms visibly dephase over just a few carrier cycles.
