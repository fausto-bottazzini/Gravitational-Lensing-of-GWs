# Case A: the chirp, static lens

Everything below is produced by `run.py` (`python3 cases/case_A_chirp/run.py`,
~5 s); every number quoted here is in `provenance/numbers.json`, written by
that same run, not hand-typed. System parameters: `src/gwlens/system.py`
(m1=20, m2=15 Msun; lens M_L=5e4 Msun at D_L=5000 pc; outer orbit a_out=1.82
AU, P_out=4 d, i_out=87 deg).

## Regime

The inner binary inspirals from f=10 Hz to f_isco=125.6 Hz in
**15.18 s** — **4.4e-5** of the 4-day outer period (confirmed in
`tests/test_system.py::check_static_lens_regime`). The lens is frozen for
the entire chirp: a single impact parameter,

    y_A = 1.589

taken (see `run.py`'s comments) at the outer orbit's phase of closest
projected approach — an explicit, stated choice of merger phase, not a fit;
a different phase would move y_A elsewhere in Figure 3's pattern without
changing the story.

## Wave-optics parameter across the band

    w(f=10 Hz)      = 61.9
    w(f=f_isco)      = 777.6

Both are well above the `w_geo_threshold=30` where `F_hybrid` switches from
the exact closed form to the (independently validated, `tests/test_
waveoptics.py::check_geometric_optics_limit`) geometric-optics asymptotic
form — i.e. every F(f) value plotted below is the fast geometric-optics
formula, not a slow mpmath evaluation, and is accurate to <4% by that test.

This means Case A does **not** show a clean, few-fringe diffraction pattern
(that is Case B's regime). It shows the two-image **interference** term
`exp(i*w*DeltaT)` oscillating very rapidly in frequency —
`caseA_F_of_f.png` looks like a filled band over the full 10-125 Hz sweep and
only resolves into distinct fringes when zoomed to a ~3 Hz window (bottom
panel of that figure). This is still genuinely wave optics (interference,
not the smooth `w->infinity` envelope), just deep enough into it that the
fringes are individually unresolvable at the full-band frequency scale — the
predicted, checked, transition toward geometric optics as merger approaches.

|F(f)| itself stays modest: it oscillates between about 0.79 and 1.27
throughout (`abs_F_at_f_start=1.236`, `abs_F_at_f_isco=1.137` are two
individual samples of that oscillation, not its envelope).

## What lensing does to the waveform

- **Peak-sample amplification: 0.487** (the single largest |h(t)| sample is
  *smaller* lensed than unlensed).
- **RMS amplification: 0.960** (power averaged over the whole chirp is
  close to, and slightly below, unlensed).

These two numbers looking different from each other, and both being close to
or below 1 despite |F|>1 through much of the band, is the point, not an
error: F(f) modulates **phase** as well as amplitude (`caseA_F_of_f.png`,
bottom sub-panel: arg F sweeps through a full 2*pi cycle roughly once per
fringe). The lensed waveform is a dephased, not simply rescaled, copy of the
unlensed one, so a single peak sample is not a robust amplification
statistic here — a genuinely wave-optics (not geometric-optics) statement:
in the geometric-optics limit lensing IS just a magnification (a real
scalar), but Case A sits close enough to that limit that this distinction
is subtle rather than dramatic. `caseA_strain_time.png` shows both curves;
`caseA_detector_envelope.png` shows the analytic-signal envelope of each
(Hilbert transform) as the idealized single-detector view asked for — note
this envelope construction assumes a slowly-varying carrier and becomes less
reliable in the last ~0.1 s where the instantaneous frequency itself is
changing fastest, visible as the sharp upturn right at merger.

## The extended diffraction pattern (`caseA_diffraction_pattern.png`)

The point lens is axisymmetric, so the amplification pattern on the source
plane is a set of concentric rings around the lens, with a central bright
spike at y=0 (the Einstein-ring diffraction spike — finite in wave optics
where the geometric-optics magnification formally diverges, see
`tests/test_waveoptics.py::check_einstein_ring_regularization`; this is the
GW analogue of the Arago/Poisson spot in optical diffraction, Deguchi &
Watson 1986). Left panel: the full pattern at w=61.9, source marked at y_A.
Right panel: zoomed to |Δy|<0.15 around the source at w=777.6 — here the
rings are so closely spaced (radial period ~2π/(w·y) ≈ 0.005 in y) that,
locally, they look like near-straight parallel fringes rather than curved
rings; both are the same physical pattern, just resolved at very different
scales, which is itself the wave-to-geometric-optics story of this whole
case.
