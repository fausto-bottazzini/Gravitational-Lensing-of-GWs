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

## Wave-optics parameter across the band

    w(f=10 Hz)      = 61.9
    w(f=f_isco)      = 777.6

Both are well above the `w_geo_threshold=30` where `F_hybrid` switches from
the exact closed form to the (independently validated, `tests/test_
waveoptics.py::check_geometric_optics_limit`) geometric-optics asymptotic
form — i.e. every F(f) value plotted below is the fast geometric-optics
formula, not a slow mpmath evaluation, and is accurate to <4% by that test.

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
- **RMS amplification: 1.056** (power averaged over the whole waveform —
  by Parseval's theorem this depends only on |F(f)|, so it is identical
  whichever Fourier-phase reference convention is used; see below).

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
"REFERENCE-PHASE NORMALIZATION FIX"), so the combined lensed merger peaks
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
