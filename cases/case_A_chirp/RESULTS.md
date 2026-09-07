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

- **Peak-sample amplification: 1.091** (the single largest |h(t)| sample is
  larger lensed than unlensed).
- **RMS amplification: 1.056** (power averaged over the whole waveform —
  by Parseval's theorem this depends only on |F(f)|, so it is identical
  whichever Fourier-phase convention is used; see the causality note below).

F(f) modulates **phase** as well as amplitude (`caseA_F_of_f.png`, bottom
sub-panel: arg F sweeps through a full 2*pi cycle roughly once per fringe),
so the lensed waveform is a dephased, not simply rescaled, copy of the
unlensed one — a genuinely wave-optics (not geometric-optics) statement: in
the geometric-optics limit lensing IS just a real magnification.

**A note on causality, and a genuine echo.** An earlier version of this
analysis had the project's Fourier sign convention backwards (documented in
`wiki/log.md` and `wiki/conventions.md`); fixed, `caseA_strain_time.png`'s
merger/ringdown zoom now shows something that can look alarming at first
glance: the combined lensed merger peak sits **0.60 s *before*** the
unlensed one (`peak_time_lensed_minus_unlensed_s=-0.6025`). This is *not*
acausal. Case A's lensed signal is the coherent sum of two images — a
strong (minimum-time) and a weak (saddle-point) image — separated by a
**fixed** delay `ΔT(y_A)=image_time_delay_seconds=3.434 s`
(`src/gwlens/waveoptics.py::time_delay_difference`); both are built only
from the source's own past. Near merger the GW frequency sweeps fast enough
that F(f)'s phase varies rapidly across the band, so the *interference*
between the two images can shift where their *sum*'s envelope peaks, by an
amount unrelated to ΔT and with no information arriving early. The
independent, unambiguous confirmation is `caseA_strain_time.png`'s third
panel and `caseA_detector_envelope.png` (log-scale): a second, genuinely
separate, ~2-orders-of-magnitude-fainter copy of the *entire* merger/
ringdown waveform appears at `echo_peak_time_s` — this project's proxy for
the weak image alone, since the two lensed images share no waveform
features until they've each individually finished. Naively predicting this
echo's arrival as `t_peak_unlensed + ΔT` misses by ~0.6 s
(`echo_peak_time_minus_prediction_s=-0.603`) — the *same* 0.6 s as the
near-merger interference shift above, to three digits, because ΔT applies
from the strong image's own (interference-free) arrival, and the combined
near-merger peak used as its proxy is shifted by that same interference.
Rather than assume either reference point, `run.py` searches the lensed
envelope directly for this echo's actual peak and logs both the observed
time and its offset from the naive prediction. At that instant the
*unlensed* envelope (`echo_time_env_unlensed=1.5e-22`) is at its FFT-array
noise floor — the unlensed signal has long since ended — while the lensed
envelope (`echo_peak_env_lensed=9.65e-21`) is ~60x above it: a real,
separate, later-arriving signal, present only because of the lens.

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
