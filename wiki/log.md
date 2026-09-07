# Log

Append-only. What was decided, ingested, corrected, and when.

- **2026-09-07** — Project started. Read `GW-AI-course/final-project.html` (the
  assignment) and `GW-AI-course/day5/exercise/{README.md,b/,c/}` for the
  provenance pattern (`RESULTS.md` + `provenance/claims.yaml` +
  `provenance/numbers.json`) this repo now follows. No solved content copied,
  only the scaffold shape.
- **2026-09-07** — User decisions locked in (see chat, not reproduced here):
  repo name `gw-lensing-hierarchical-triple`, public, pushed by the user at the
  end (no `gh` CLI available in this environment); lens mass in the IMBH range
  so wave-optics diffraction is visible in *both* cases rather than picking the
  astrophysically-standard SMBH mass (which would put the chirp case deep in
  the geometric-optics limit); the chirp case and the monochromatic case are
  the *same* hierarchical triple at two epochs, not two unrelated examples.
  Third deliverable added on request: a theory PDF written like a textbook
  chapter (`theory/theory.tex` → `theory.pdf`), separate from the results
  HTML+PDF.
- **2026-09-07** — `wiki/conventions.md` seeded: units, Fourier sign
  convention (`h(f)=∫h(t)e^{-2πift}dt`), lensing notation (`w`, `y`, `F(f,y)`),
  orbit element convention, scalar-field-reduction statement to be derived in
  `theory.tex` Sec. 2.
- **2026-09-07** — `src/gwlens/waveoptics.py` built and self-checked:
  - Paczynski (1986) total-magnification formula vs. the image-sum from
    `image_positions`+`magnification` agree to float precision (all `y`
    tested) — confirms the point-lens Jacobian formula `mu=x^4/(x^4-1)`.
  - `F_point_lens` (closed-form, mpmath) -> 1 as `w -> 0`, as required physically.
  - **Caught and fixed a bug**: first version of `F_geometric_optics` matched
    `F_point_lens` in `|F|` at large `w` but NOT in phase (off by
    `w*phi_+(y)`, growing with `w`). Diagnosed by comparing the *ratio*
    `F_point_lens/F_geometric_optics` and recognizing `w*phi_+(y)` in the
    residual phase; traced to `F_point_lens` (the standard closed form)
    carrying the minimum image's own, un-subtracted Fermat-potential phase,
    which the naive stationary-phase writeup silently drops. Fixed by
    multiplying `F_geometric_optics` by `exp(i*w*phi_+)`. After the fix,
    relative error between the two independent methods is 3e-4 at `w=10` and
    3e-5 at `w=1000`, shrinking with `w` as expected for an asymptotic
    expansion — recorded as the check in
    `tests/test_waveoptics.py::test_geometric_optics_limit`.
- **2026-09-07** — `src/gwlens/chirp.py` (leading-order chirp): inverting
  `time_to_merger(f)` for `f(t)` by hand produced a wrong exponent on both
  `(5/256)` and `Mchirp` (a bracket-power slip); caught immediately because
  `freq_of_time(t=0)` did not reproduce the input `f0`. Fixed by writing the
  inversion out as separate factors instead of one clever bracket. Then
  `tests/test_chirp.py`'s own first draft had two more errors (this time in
  the *tests*, not the code): expected `f_isco(60 Msun)` in 100–300 Hz
  (checked by hand: correct value is ~73 Hz, matching the standard
  `f_isco ~ 4400 Hz * (Msun/M)` rule of thumb), and expected
  `Mchirp/Mtotal=2^(-1/5)` at equal mass instead of the algebraically correct
  `2^(-6/5)`. Both test bugs fixed by re-deriving the expected numbers by
  hand rather than guessing; `src/gwlens/chirp.py` itself needed no further
  changes. All 4 checks in `tests/test_chirp.py` now pass.
- **2026-09-07** — `src/gwlens/system.py`: pinned the one triple used
  everywhere. m1=20, m2=15 Msun (Mchirp=15.05, f_isco=125.6 Hz); lens
  M_L=5e4 Msun ("massive IMBH", between the classic IMBH and SMBH ranges --
  needed to keep w=O(1) at the Case B frequency while staying a physically
  legitimate point-mass lens; logged here rather than silently drifting past
  the "10^2-10^4" range floated when the user picked the IMBH option) at
  D_L=5000 pc; outer orbit a_out=1.82 AU from Kepler's third law at
  P_out=4 d, i_out=87 deg (near edge-on, needed for the front/back "repeated
  lensing" split in `geometry.impact_parameter_of_time` to actually bite),
  e_out=0. Case A starts the chirp at f=10 Hz (merger takes 15 s, 4.4e-5 of
  P_out -- static lens, confirmed). Case B observes at fixed f=0.05 Hz for 6
  outer periods (24 d) out of 240 d to merger (4% frequency drift -- quasi
  -monochromatic, confirmed). w_B=0.31 (clear diffraction regime); w at
  Case A's start is 62, rising to 778 at f_isco -- NOT O(1) as originally
  hoped when the IMBH-mass option was picked, but not a problem: at those w
  the exact F(f) still shows real, resolved interference fringes in
  frequency (consistent with the validated geometric-optics limit), just
  past the point where they average into a smooth envelope -- reframed Case
  A's story from "diffraction pattern" to "dense interference fringes
  approaching the geometric-optics limit as merger approaches," arguably the
  more interesting story for a chirp anyway. All 4 checks in
  `tests/test_system.py` pass (one iteration: `quasi_monochromatic_regime`
  first failed testing frequency drift over `tau_B/2`, close to formal
  coalescence where f always diverges regardless of parameters; fixed to
  test drift over the actual 6-period Case B observing window instead).
- **2026-09-07** — User asked (mid-session) to use a waveform from an actual
  search-template family rather than the hand-built leading-order chirp.
  Tried `pip install pycbc` first (the course's own tool for days 2/5); it
  hung for a very long time resolving/building the lalsuite dependency
  chain with no output, consistent with lalsuite's historically poor native
  -Windows support, and was killed rather than waited out indefinitely.
  Implemented `src/gwlens/taylorf2.py` instead: the standard frequency
  -domain, stationary-phase-approximation (SPA) PN inspiral waveform, the
  actual functional form LIGO/Virgo low-mass CBC template banks use for the
  inspiral part (LALSimulation/PyCBC's "TaylorF2" approximant), truncated
  at 2PN (the phase coefficients below that order are ones I am confident
  reproducing correctly from memory without a verified reference to check
  numbers against; 2.5-3.5PN would add several more intricate rational/log
  terms and real risk of a silent error, so left out and stated as a scope
  choice, not hidden).
  **A real, serious bug, caught and fixed**: the first version used the
  literature's usual sign for the SPA phase, which assumes the OPPOSITE
  Fourier convention from this project's own (`wiki/conventions.md`).
  Feeding it straight to `numpy.fft.irfft` reconstructed a "chirp" whose
  peak drifted with the amount of zero-padding (17s at 2x padding, 245s at
  16x padding) instead of sitting still at `t_c` — the unmistakable
  signature of a sign error aliasing power to the wrong end of the
  (periodic) IFFT window. Re-derived the correct sign directly from the
  stationary-phase condition applied to THIS project's inverse-transform
  convention (`h(t)=int h(f) e^{+2pi i f t} df`, same sign as `numpy.
  ifft`), confirmed empirically by flipping the overall sign of `Psi(f)`
  and checking the reconstructed peak lands at `t_c` (now to <0.2%) and,
  critically, STOPS moving as padding increases (checked at 2x, 4x, 8x
  padding: peak time now identical to machine precision across all three —
  `tests/test_taylorf2.py::check_ifft_reconstructs_chirp_at_tc`). Also
  caught by this: the earlier `leading_order_group_delay` check (matching
  the LEADING-order term's derivative against `chirp.time_to_merger`) had
  passed throughout, on BOTH sign conventions after appropriately negating
  the comparison target -- it only tests internal algebraic self
  -consistency of one formula, not which overall sign the physical IFFT
  needs, so it could not have caught this on its own. Recorded here so the
  next reader trusts `check_ifft_reconstructs_chirp_at_tc`, not just the
  group-delay check, as the test that actually validates the sign.
  `cases/case_A_chirp/run.py` rewritten to build the waveform natively in
  the frequency domain with `taylorf2.htilde` (banded to
  `[F_A_START_HZ, f_isco]`, half-cosine-tapered edges) and apply lensing as
  a single `F(f)` multiplication on that same grid, replacing the previous
  time-domain-first-then-FFT construction — simpler, and now uses a real
  search-template phase instead of the leading-order-only one.
- **2026-09-07** — Independent fresh-agent review run (`checks/independent_review/REVIEW.md`):
  `reproduce.sh` (checks + cases + theory + report) reproduced everything with zero
  intervention, all 26 checks passed, and all `RESULTS.md`/report/theory numbers traced
  correctly to `provenance/numbers.json`. Verdict: the project holds up — the physics
  checked out against the cited literature forms with no errors found, provenance
  discipline was followed consistently, and scope/limitations are stated honestly; the
  only issues found were minor stale-comment/docstring inconsistencies
  (`waveoptics.py::F_hybrid`'s docstring still says `w_geo_threshold=20` though the code
  default is 30; `test_system.py`'s comment says `T_OBS_B_S = 8 outer periods` though the
  code and every result say 6) plus two nitpick-level rounding/prose mismatches, none of
  which affect any reported number or figure.
- **2026-09-07** — User pushed back on the TaylorF2-only waveform (fair: no merger, no
  ringdown, formally diverges and has to be cut off by hand at f_isco) and asked to try
  `pycbc` again via WSL2, since this machine has it installed and the professor uses Linux
  anyway. `pip install pycbc` inside a fresh WSL2 Ubuntu venv worked cleanly (~2 min,
  prebuilt wheels for numpy 2.3.5/scipy 1.16.3/lalsuite 7.26.15/pycbc 2.11.0, Python
  3.14) — confirms the native-Windows failure earlier was Windows-specific, not
  pycbc/lalsuite being broken in general. `get_fd_waveform(approximant="IMRPhenomD", ...)`
  generates a real inspiral-merger-ringdown waveform in ~1 s.
  **A second alignment bug, same symptom as the earlier TaylorF2 sign bug**: pycbc's FD
  waveform is not delivered pre-aligned to a convenient positive merger time — an
  unshifted `irfft` put the peak right at the edge of the FFT window (`t=32.59s` for a
  32.6s window), with most of the inspiral wrapped to the "wrong" end. Fixed with an
  explicit frequency-domain time shift, `H *= exp(-i*2*pi*f*(t_target - t_measured))`
  (the sign of the exponent found empirically by testing both and checking which one
  actually moved the peak to `t_target` — see `src/gwlens/imr_waveform.py`); confirmed
  the shifted peak lands within 0.003s of the 27.71s target.
  `src/gwlens/imr_waveform.py` added: tries `pycbc`/IMRPhenomD (Khan et al. 2016) first,
  falls back to `taylorf2.py` (inspiral-only) if `pycbc` is not importable — so
  `cases/case_A_chirp/run.py` runs unmodified on any platform, just with a better
  waveform where `pycbc` is available. This repo's own committed Case A figures/numbers
  were generated inside WSL2 (`.venv-wsl/`, gitignored, not part of the reproducibility
  contract for platforms without WSL) — see README.md for exactly which command runs
  where. Payoff: with a real merger, the wave-optics dephasing already documented for the
  inspiral is now visible by eye in the merger itself — the lensed merger transient comes
  out visibly smeared (broader, lower peak) compared to the unlensed one's sharp spike,
  because the broadband merger samples many F(f) interference fringes at once where the
  narrowband inspiral only sampled a slowly-varying few. Not derived in detail, reported
  as an observed, physically-plausible consequence of the already-validated F(f), not an
  independently re-derived result on its own.
- **2026-09-07** — User asked for `report.html` to be a lighter, visual "presentation" (not
  dense), with interactive visualizations: the outer orbit's motion, the extended
  interference/diffraction pattern in space, and an idealized detector view, synced -- mainly
  for Case B, something similar for Case A. Built this as: `cases/case_B_monochromatic/run.py`
  now also exports `caseB_animation_data.json` (orbit track y1(t)/y2(t)/lensed over ~1.5
  periods, |F(t)|^2 and envelope over the full 6) and a clean track-free
  `caseB_pattern_only.png`; `case_A_chirp/run.py` exports `caseA_animation_data.json`
  (windowed waveform + instantaneous-frequency mapping, zoomed F(f) fringes). `report/
  report_template.html` + `report/build_report.py` assemble the final `report.html` by
  INLINING both JSON payloads as `<script>` data rather than `fetch()`-ing them: a page opened
  directly as a `file://` URL (the expected way to view this) cannot `fetch()` a local JSON
  file (CORS blocks it), so inlining is required for the page to work with no local web server
  -- added `run_report` in `reproduce.sh` runs `build_report.py` alongside the LaTeX build.
  Vanilla-JS canvas plots + a CSS-positioned marker over the pattern PNG, no chart library
  (stays inside the project's own dependency footprint).
  **Could not visually test in an actual browser** (`claude-in-chrome` reported no extension
  connected in this environment) -- verified what could be verified without one instead:
  `node --check` on the extracted script (syntax), and running the actual page script under
  Node with a minimal DOM/canvas-context mock (`document.getElementById`, `classList.toggle`,
  a no-op 2D context) driving every control's event handler across the full scrub range
  (0-1000 in steps of 5) with no exceptions thrown. That pass caught one real bug: `y` (and so
  the marker's plotted position) formally diverges at each `D_LS=0` crossing (the same
  divergence documented in `theory.tex` Sec. 4.2), which sent the marker's CSS `left`/`top` to
  values like `191.93%` -- clipped invisibly by the stage's `overflow:hidden` rather than
  crashing, but a marker silently vanishing reads as a bug to a viewer. Fixed by clamping the
  marker's position to the visible stage edge (`[1.5%, 98.5%]`) instead, re-verified the clamp
  holds across the full scrub range. This DOM-mock pass is a real but partial substitute for
  actually looking at the rendered page (no verification of layout, CSS, or that things look
  right) -- flagged here rather than implied to be equivalent to visual QA.
- **2026-09-07** — Second independent fresh-agent review run, after the pycbc integration, the
  Spanish translation and the interactive `report.html`
  (`checks/independent_review/REVIEW_v2.md`). This one *did* render `report.html` for real
  (headless Chrome `--screenshot`, plus four temporary copies with an injected scrub-setter to
  drive both animations). Verdict: the reproducibility machinery still holds up --
  `reproduce.sh` (checks + cases + theory + report) ran end to end with no intervention, 26/26
  checks passed, `report.html` rebuilt byte-identical, Case B reproduced bit-for-bit, and every
  number in `RESULTS.md`/`report`/`theory.tex` Table 4.1 traced to its `numbers.json` (the first
  review's two nitpicks are fixed) -- but it found one **critical physics bug the first review
  missed**: `F` is written in Takahashi & Nakamura's Fourier convention (delay = `+i*w*DeltaT`)
  and then multiplied into a waveform that is inverted with `numpy.fft.irfft`, whose convention
  is the opposite, so Case A's saddle-point image arrives 3.43 s *before* the merger instead of
  after it, and `wiki/conventions.md`'s Fourier paragraph states the wrong pairing as its
  justification -- the same class of error the log records catching twice in `taylorf2.py` and
  `imr_waveform.py`, missed in the one file the project is actually about, because every
  `tests/test_waveoptics.py` check compares evaluators of `F` to each other (all sharing the
  convention) and none checks the time-domain response. `|F|` is unaffected, so all of Case B,
  every diffraction figure and `rms_strain_amplification` stand; `peak_strain_amplification`
  (0.962) and the "smeared lensed merger" story do not, and must be recomputed in WSL. Also
  found: the Case A animation is decimated to 84 Hz against a 125.6 Hz signal (the merger panel
  is aliased, and its apparent peak ratio 0.60 contradicts the 0.962 printed beside it, while
  the already-exported Hilbert envelopes that would have fixed it go unused); the Case B marker
  is clamped to the frame edge for 78% of the orbit and is plotted at a fixed median `theta_E`
  rather than `theta_E(t)` (19% too far out at peak lensing); Ch. 5's retro-lensing argument
  forward-references an estimate that exists in neither section; two hard-coded `Capítulo~5`
  refs now point at the wrong chapter; §2.1's `lambda_GW << R_E` criterion is violated by
  Case B's own `w_B=0.31` design; `reproduce.sh` never actually uses the venv it creates on
  Windows (no `python3` in a Windows venv, so pip installs into the global interpreter); and
  `wiki/todo.md`/`index.md` still describe pycbc as abandoned and the animation as dropped.
  Full list, with file/line for each, in `REVIEW_v2.md`.
- **2026-09-07** — Fixed the critical `REVIEW_v2.md` sign-convention bug and everything it
  touched, code-side (exposition/results files deliberately left for a later pass, per
  instruction). `wiki/conventions.md`'s Fourier-convention paragraph had it backwards: in this
  project's stated convention (`h(f)=∫h(t)e^{-2πift}dt`, matching `numpy.fft.rfft`/`irfft`), a
  LATER-arriving copy of a signal gets phase `e^{-i2πfΔt}`, not `+i2πfΔt` as the doc claimed --
  standard delay theorem, checked directly against `numpy.fft.irfft(numpy.fft.rfft(pulse) *
  numpy.exp(-2j*np.pi*f*dt))` on a test pulse. `waveoptics.py`'s four `F(w,y)` evaluators
  (`F_bruteforce_2d`, `F_radial_1d`, `F_point_lens`, `F_geometric_optics`) had all been derived
  self-consistently from the wrong pairing, so they agreed with each other perfectly while
  disagreeing with `numpy`'s actual FFT convention -- fixed by conjugating each evaluator's
  output (`|F|` is conjugation-invariant, so this changes nothing about Case B, the diffraction
  figures, or `rms_strain_amplification`; it only changes phase-sensitive time-domain
  reconstructions, i.e. Case A's waveform shape). Added
  `tests/test_waveoptics.py::check_causality_of_lensed_pulse`, first in the CHECKS list: builds
  a short Gaussian pulse, lenses it exactly as the case script does (FFT, multiply by `F`,
  IFFT), and checks the weak (saddle-point) image lands AFTER the main pulse, not before --
  this is the only kind of check that can catch a *global* sign-convention bug, since every
  other check in the file compares `F` evaluators to each other and both share whichever
  convention is (mis)implemented. Re-ran the full case A pipeline in WSL (pycbc/IMRPhenomD):
  `peak_strain_amplification` goes from 0.962 (wrong-sign artifact) to 1.091;
  `peak_time_lensed_minus_unlensed_s=-0.60s` (the combined near-merger peak now looks slightly
  *earlier*, not later -- this is real, but it's an interference effect between two
  always-past images, not acausality; see below). Also caught and fixed while re-deriving these
  numbers: the sample rate `fs=4*f_isco=502.5 Hz` (chosen back when only an inspiral-only
  waveform existed) undersamples real IMRPhenomD ringdown content, which the review's own
  screenshots showed carrying >0.1% power out to ~640 Hz -- raised to `fs=2048 Hz`.
  User asked, on seeing the -0.60s shift: "tiene sentido eso?" (does that make sense?). It does:
  Case A's lensed signal is the coherent SUM of a strong (minimum-time) image and a weak
  (saddle-point) image separated by a FIXED delay `ΔT(y)=image_time_delay_seconds=3.434s`
  (both images individually causal, both built only from the pulse's own past). Near merger the
  frequency sweeps fast enough that `F(f)`'s phase varies rapidly across the band, so the two
  images' interference pattern can shift where the *combined* envelope peaks, by an amount
  unrelated to `ΔT`; it does not mean any information arrived early. Independent confirmation:
  a real, physically separate SECOND copy of the whole merger *does* appear later, as
  `image_time_delay_seconds` predicts -- a distinct echo of the full IMR waveform, ~2 orders of
  magnitude fainter, visible in the lensed strain and absent (down to a ~1e-22 FFT-artifact
  floor) in the unlensed strain at the same instant. Where exactly it lands was itself an
  instructive check: naively adding `ΔT` to the UNLENSED merger time misses it by ~0.6s --
  because the fixed geometric delay applies from the strong image's own (interference-free)
  arrival, and the combined near-merger peak (used as its proxy) is itself shifted by
  interference, by almost exactly that same ~0.6s (`echo_peak_time_minus_prediction_s=-0.603s`,
  matching `peak_time_lensed_minus_unlensed_s=-0.6025s` to three digits -- not a coincidence,
  the same interference effect explains both numbers). Rather than assume a reference point,
  `case_A_chirp/run.py` now searches the lensed envelope for the actual local maximum near the
  predicted delay and logs both the observed time and its offset from the naive prediction
  (`echo_peak_time_s`, `echo_peak_time_minus_prediction_s`, `echo_peak_env_lensed`,
  `echo_time_env_unlensed`) -- more honest, and the discrepancy is itself a small physics
  lesson worth keeping visible rather than silently fixing away. `caseA_strain_time.png` gained
  a third panel zoomed on this echo (a fainter but recognisable replay of the whole chirp/
  merger/ringdown shape, exactly as the geometric-optics superposition picture predicts) and
  `caseA_detector_envelope.png` switched to a log-scale y-axis spanning both peaks, since the
  echo is ~2 orders of magnitude below merger and invisible on a linear scale. Regenerated
  `tests/CHECKS_imr_waveform.json` via WSL so the committed file shows
  `pycbc_merger_time_alignment` genuinely passing rather than `SKIPPED` (it had last been
  committed from a native-Windows run, where pycbc isn't installed) -- and re-learned, again,
  that running `reproduce.sh checks` on native Windows afterwards silently overwrites this one
  file back to `SKIPPED` (same class of issue as the earlier `case_A` overwrite, see above);
  the fix there is procedural, not code: always regenerate this specific file from WSL *last*,
  right before committing. Full 6-file test suite re-run clean on both native Windows
  (TaylorF2 fallback path) and WSL (pycbc path).
- **2026-09-07** — User pass over the presentation layer (not the physics): reordered
  `README.md` (physics before logistics), trimmed `wiki/conventions.md`'s Fourier section to
  state-the-convention-and-move-on (full narrative stays in this log), fixed a real duplication
  bug in that same page (`y` defined twice under different-looking names with no note they're
  the same quantity; `psi(x)=ln|x|`'s "defined up to a constant" caveat missing entirely, which
  is exactly the fact behind the reference-phase bug two entries up). Removed
  `checks/independent_review/` from the repo (kept locally, `.gitignore`d): three working
  documents from getting here, not part of the deliverable.
  Added `bibliography/` with the PDF of every arXiv-available cited paper, and switched both
  `theory.tex`/`report.tex` from author-year to numbered, hyperlinked citations
  (`natbib`'s `numbers` option; each downloaded paper's `refs.bib` entry gets a `url` pointing
  at its local PDF, rendered as a link by `unsrtnat`). This surfaced two real bugs: both
  documents were printing an auto-generated bibliography heading (natbib's own, English,
  `\bibname`/`\refname`) directly under a manually-added one, wasting a whole near-blank page
  in `theory.tex` (`\chapter*` always starts fresh); and `report.tex`'s figures were `[h]`-only
  with no `t`/`b`/`p` fallback, so the ones that didn't fit backed up into a float queue that
  only drained after the bibliography, dumping six pages of pure floating figures at the very
  end instead of each figure sitting near the text that refers to it (fixed with the `float`
  package's `[H]`).
  While fixing the citation links, re-grepped `theory.tex`'s own `\chapter{...}` commands to
  double check chapter numbers before touching any `Cap\'itulo~N` cross-reference, and got a
  different count than earlier in this session: one chapter (`M\'as all\'a de este proyecto`,
  labelled `ch:beyond`) is written as `\chapter[short]{long title}`, which a plain
  `grep "chapter{"` silently misses. That means the *first* count (used earlier to reject
  `REVIEW_v2.md` Finding 6, "these `Cap\'itulo~5` refs should be `Cap\'itulo~6`") was wrong and
  the review was right: the binaria-interna/waveform chapter genuinely is Chapter 6, not 5.
  Fixed the two hardcoded refs, and — so this can't silently happen again — replaced every
  hardcoded chapter number anywhere in `theory.tex` with `\label`/`\ref` (the two chapters
  that had no label yet, "Planteo y alcance" and "La binaria interna como fuente", now do).
  Lesson logged plainly: a grep-based structural check is only as good as the pattern, and a
  literal-substring pattern missed an optional-argument macro call outright, silently, with no
  error -- worth remembering next time a "let me just grep for X" verification feels
  sufficient.
  Separately, fixed several figures after actually looking at the committed PNGs rather than
  just the code: legend boxes with no headroom covering real data in three panels (two in
  `caseB_detector_view.png`, one in the now-removed `caseB_y_of_t.png`); a diffraction-pattern
  colorbar labelling its unit-value tick "$10^0$" instead of "$1$"; `caseB_y_of_t.png` itself
  dropped entirely (the two formal divergences dominate the plot and it added nothing beyond
  what `fraction_of_time_lensed` and `caseB_pattern_with_orbit.png` already show);
  `caseB_detector_view.png`'s zoom-panel title claimed to show "the peak of one lensing pulse"
  when a pulse's actual shape (sinc-like diffraction ringing) lives on a ~day timescale, far
  wider than that panel's +-200s window -- retitled to say what it actually shows (carrier-cycle
  dephasing at the peak instant, not the pulse's own shape). Most substantively,
  `caseA_F_of_f.png`'s full-band panel plotted the raw, ~400-fringe |F(f)| curve at a scale
  where it renders as a dense, uninformative-looking blob; since F is in the geometric-optics
  regime for the whole chirp (already established), that envelope is not merely hard to
  resolve visually, it is EXACTLY two constant lines, `sqrt(mu_+) +- sqrt(|mu_-|)`, for every
  f in the band -- redrawn to show that constant band directly. Also derived and logged that
  the fringe period itself is exactly constant across the band too
  (`fringe_period_hz = 1/image_time_delay_seconds`, since `Delta_T` is fixed for a static lens
  and the oscillation phase `2*pi*f*Delta_T` is exactly linear in `f`), so the existing zoom
  inset -- drawn at the start of the band for concreteness -- is representative of any
  equal-width window anywhere in it, not something special about the start.
  Finally, `theory.tex`'s "checkbox" asides describing bugs found and fixed during development
  (the sign-convention fix, the reference-phase fix, the TaylorF2/pycbc alignment fixes) were
  rewritten as plain formal derivation/verification notes -- theory.tex is meant to read as a
  textbook chapter, not a debugging diary; the full development narrative already belongs here,
  in this log, and stays here.
