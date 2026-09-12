# Case A: the chirp, static lens

Everything below is produced by `run.py` (`python3 cases/case_A_chirp/run.py`,
~5 s); every number quoted here is in `provenance/numbers.json`, written by
that same run, not hand-typed. System parameters: `src/gwlens/system.py`
(m1=20, m2=15 Msun; lens M_L=5e4 Msun at D_L=5000 pc; outer orbit a_out=1.82
AU, P_out=4 d, i_out=87 deg).

**Waveform**: `src/gwlens/imr_waveform.py` tries `pycbc`'s `IMRPhenomD`
(Khan et al. 2016) first -- a published, NR-calibrated approximant with a
genuine inspiral-merger-**ringdown**, exactly what LIGO/Virgo search and
parameter-estimation pipelines use -- and falls back to the hand-built
`src/gwlens/taylorf2.py` (restricted 2PN, inspiral-only, cut off at f_isco)
if `pycbc` is not importable. `pycbc` is not pip-installable on native
Windows in any reasonable time (a first attempt hung indefinitely resolving
`lalsuite`); this repository's committed numbers/figures were generated in
WSL2 Ubuntu, where `pip install pycbc` works cleanly (see **README.md** for
the exact reproduction path on each platform). `numbers.json`'s
`waveform_source` field always records which path actually ran for a given
run -- here, `IMRPhenomD (pycbc/LALSimulation, Khan et al. 2016)`.

Both waveform generators are frequency-domain, on the same grid lensing is
applied on (`h_lensed(f) = F(f,y_A) h_unlensed(f)`, one `irfft` each for the
time-domain figures) -- `pycbc`'s FD waveform is not delivered pre-aligned
to a convenient positive merger time, which was diagnosed (the same
symptom as an earlier `taylorf2.py` bug: the reconstructed merger time sat
at the very edge of the FFT window, most of the inspiral wrapped to the
"wrong" end) and fixed with an explicit frequency-domain time shift,
checked to land the merger within 0.01% of the requested target time. See
`wiki/log.md` for both this and the earlier `taylorf2.py` Fourier-sign bug
(still relevant: it is what the fallback path uses).

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

## The source is frozen too — a separate claim, checked separately

"The lens is frozen" is about the *geometry*: y_A does not move. It says
nothing about the source's **motion**, and in [Case B](../case_B_monochromatic/RESULTS.md)
that motion — the same outer orbit, seen over days instead of seconds —
turns out to imprint 90 GW cycles of Roemer/Doppler phase, dwarfing
everything the lens does. So the claim that Case A can ignore it is not
inherited from the static-lens argument and is not assumed; it is computed
(`src/gwlens/doppler.py`, and
`tests/test_doppler.py::check_case_A_roemer_is_negligible`):

- The merger is placed at the outer orbit's phase of closest approach
  (t/P_out = 0.25 — chosen for the lensing, not for this), which is exactly
  where the line-of-sight velocity crosses zero. The mean over the chirp is
  `mean_beta_los_over_chirp` = **2.3e-6**.
- Over a 15 s window out of a 4-day orbit, a constant delay and a constant
  slope in the delay are anyway **exactly degenerate** with t_c and with the
  chirp mass — unmeasurable in a single event. Only the curvature (the
  line-of-sight acceleration) is non-degenerate, and that residual is
  `roemer_residual_ptp_s` = **8.6e-6 s**, i.e.
  `roemer_residual_phase_at_fisco_rad` = **0.0068 rad** at the highest
  in-band frequency.
- It is therefore reported, not applied. Applying it would mean resampling
  h(t) onto a retarded time grid, whose interpolation error would itself
  exceed the 0.0068 rad being modelled.
- The one kinematic effect that is *not* small is also not a waveform
  feature: the gravitational redshift in the lens's potential plus the
  transverse Doppler shift combine, for a circular geodesic, into
  `(1-3GM_L/(a c²))^(-1/2)`, giving
  `orbital_redshift_factor_minus_1` = **4.08e-4**. Constant ⇒ exactly
  degenerate with the chirp mass ⇒ a bias on any inferred M_chirp, not
  something visible in any figure here.

## Wave-optics parameter across the band

    w(f=10 Hz)      = 61.9
    w(f=f_isco)      = 777.6

Both are well above the `w_geo_threshold=30` where `F_hybrid` switches from
the exact closed form to the (independently validated, `tests/test_
waveoptics.py::check_geometric_optics_limit`) geometric-optics asymptotic
form — i.e. every F(f) value plotted below is the fast geometric-optics
formula, not a slow mpmath evaluation.

How accurate that is, stated at the y this case actually uses rather than as
a single worst-case number: at y_A=1.589 the exact and geometric-optics
evaluators agree to `geo_vs_exact_relerr_at_f_start` = **1.9e-4** and
`geo_vs_exact_relerr_at_f_isco` = **1.5e-5**, the error falling with w as an
asymptotic expansion should. The `check_hybrid_matches_at_threshold` figure
of 3.1% is the *worst case over y* at the threshold itself, and it is set
entirely by y=0.3 — where the two images sit close together and highly
magnified, so stationary phase is at its weakest, a y nothing here
evaluates. At y_A that same threshold discontinuity is
`geo_vs_exact_relerr_at_threshold_w30` = **3.1e-4**, and the threshold is
never crossed in band anyway.

This means Case A does **not** show a clean, few-fringe diffraction pattern
(that is Case B's regime), and -- worth being precise about, since an
earlier version of this section overclaimed it -- every in-band `F` value
here comes from `F_geometric_optics` (`F_hybrid` never drops below
`w_geo_threshold=30`; confirmed independently by reconstructing
`h_lensed` from nothing but two delayed, rescaled copies of the unlensed
waveform and finding `max|h_lensed - h_two_image|/max|h_lensed| = 1.1e-13`,
machine precision). What's shown is the two-image **interference** term
`exp(i*w*DeltaT)` oscillating very rapidly in frequency, NOT the exact
diffraction integral contributing anything measurable beyond that limit.
Because y_A (and so mu_+, mu_-) is fixed for the whole chirp, the *envelope*
of `|F(f)|` is exactly constant across the entire band --
`sqrt(mu_+)-sqrt(|mu_-|)` to `sqrt(mu_+)+sqrt(|mu_-|)`, shown as a shaded
band in `caseA_F_of_f.png`'s top panel -- and the fringes inside that band
have a fixed period in frequency, `1/image_time_delay_seconds=0.29 Hz`,
the same at every f (Delta_T is fixed; the oscillation phase 2*pi*f*Delta_T
is exactly linear in f). The middle panel zooms a 3 Hz window at the start
of the band purely for concreteness -- an equal-width window anywhere else
in the band looks the same. Case B, at `w_B=0.309`, is where the exact
diffraction integral actually matters (`F_hybrid` uses `F_point_lens`
throughout there).

|F(f)| itself stays modest: it oscillates between about 0.79 and 1.27
throughout (`abs_F_at_f_start=1.236`, `abs_F_at_f_isco=1.137` are two
individual samples of that oscillation, not its envelope).

## What lensing does to the waveform

- **Peak-sample amplification: 1.065** (the single largest |h(t)| sample is
  larger lensed than unlensed).
- **RMS amplification: 1.0555** (power averaged over the whole waveform —
  by Parseval's theorem this depends only on |F(f)|, so it is identical
  whichever Fourier-phase reference convention is used; see below).

  This one is **not a free number**, and checking it against its closed form
  is the only end-to-end validation of this entire pipeline — waveform
  generator, F evaluator, Fourier convention, reference phase, window —
  against an independent textbook result. By Parseval, the squared RMS ratio
  is the |H(f)|²-weighted mean of |F(f)|² across the band; in the
  geometric-optics regime (which the whole in-band F is, above),
  |F|² = μ₊ + |μ₋| + 2√(μ₊|μ₋|)·sin(2πf·ΔT), and the oscillating term
  averages away over the ~400 fringes in the band, leaving exactly
  μ₊ + |μ₋| — the Paczynski (1986) **total magnification** A(y_A), which
  `src/gwlens/waveoptics.py::total_magnification_paczynski` implements
  independently and `tests/test_waveoptics.py::check_paczynski_magnification`
  separately verifies against the image sum. So

      rms_strain_amplification  =  sqrt(A(y_A))
      1.0554959                 vs 1.0558485        (relative 3.3e-4)

  (`total_magnification_paczynski`, `rms_amplification_predicted`,
  `rms_amplification_relerr` in `numbers.json`), with nothing fitted and no
  free parameter. The residual is the finite-fringe-count sampling of that
  average — running the same script on the TaylorF2 fallback, whose band is
  exactly [10 Hz, f_isco] with no ringdown power above it, brings the
  agreement to 9.6e-9.

F(f) modulates **phase** as well as amplitude, so the lensed waveform is a
dephased, not simply rescaled, copy of the unlensed one
(`caseA_F_of_f.png`, bottom panel: `arg F` oscillates in a *bounded* range
(~±0.23 rad here), in phase with each `|F|` fringe above it, NOT a full
2*pi wind per fringe -- `F = sqrt(mu_+) + i*sqrt(|mu_-|)*exp(-iw*DeltaT)`
traces a circle of radius `sqrt(|mu_-|)` centered on `sqrt(mu_+)`
(`F_geo`, `waveoptics.py`) as `w*DeltaT` sweeps `2*pi` per fringe; since
`|mu_-|<mu_+` here that circle does not enclose the origin, so the phase
cannot wind all the way around).

**A note on the reference phase, causality, and a genuine echo.** `F(w,y)`
is referenced to the strong (minimum-time) image's own arrival: it carries
zero phase there by construction (`src/gwlens/waveoptics.py`,
"Reference phase"), so the combined lensed merger peaks
at *exactly* the unlensed merger time
(`peak_time_lensed_minus_unlensed_s=0.0`). Case A's lensed signal is the
coherent sum of that strong image and a weak (saddle-point) image,
separated by a **fixed** delay `ΔT(y_A)=image_time_delay_seconds=3.434 s`
(`src/gwlens/waveoptics.py::time_delay_difference`); both are built only
from the source's own past -- always causal. The independent,
unambiguous confirmation: a second, genuinely separate, `echo_to_
unlensed_peak_ratio=0.242`-as-faint (matching the weak image's own
magnification, `sqrt_mu_minus_geometric=0.240`, to 1.1%) copy of the
*entire* merger/ringdown waveform appears at `echo_peak_time_s`, exactly
`ΔT` after the merger to within `echo_peak_time_minus_prediction_s=
-4.0e-5 s` -- four hundred-thousandths of a second, not the ~0.6 s an
earlier, unnormalized version of this analysis found (and mis-attributed
to interference; see `wiki/log.md` for that whole episode, including how
it was caught). At the echo instant the *unlensed* envelope
(`echo_time_env_unlensed=1.6e-22`) is at its FFT-array noise floor -- the
unlensed signal has long since ended -- while the lensed envelope
(`echo_peak_env_lensed=9.68e-21`) is ~60x above it: a real, separate,
later-arriving signal, present only because of the lens, and quantitatively
right (the `sqrt(|mu_-|)` agreement above) to better than 1.1%.

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
