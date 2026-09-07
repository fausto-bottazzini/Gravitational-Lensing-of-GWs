# Independent review

Fresh session, no prior context, no communication with whoever built this repo. Ran on
2026-09-07 on a clean-ish Windows checkout (Git Bash + PowerShell available; MiKTeX present).

## Verdict

This project holds up well. `./reproduce.sh` (checks + cases), `./reproduce.sh theory`,
and `./reproduce.sh report` all ran to completion on the first try with zero manual
intervention, all 26 correctness checks across five `tests/test_*.py` files passed, and
every number I spot-checked in `cases/*/RESULTS.md`, `report/report.html`,
`report/report.pdf`, and `theory/theory.pdf` traced back correctly to the corresponding
`provenance/numbers.json` (or, for the theory PDF's paraxial-validity table,
`theory/paraxial_validity_numbers.json`) — I found no case where a quoted number
disagreed with its source file. The physics itself is standard and, as far as I can
check it, correctly implemented: the point-lens amplification factor, its geometric-optics
and low-w limits, the Paczynski magnification cross-check, the Kepler-orbit geometry, the
TaylorF2 PN phase coefficients, and the ISCO frequency formula all match the literature
forms they cite. The project is unusually honest about its own scope and limitations
(no absolute strain calibration, 2PN truncation stated and justified, the paraxial
approximation's ~60s/crossing breakdown window quantified rather than hand-waved past,
`pycbc` install abandoned and logged rather than silently dropped). The issues I found are
all minor documentation/prose inconsistencies (stale docstring numbers, one rounding slip,
one imprecise restatement) — none of them affect any plotted figure, reported number, or
physics conclusion. I would pass this on reproducibility and scientific rigor grounds; a
strict grader might dock a fraction of a point for the stale-docstring class of issue
since it is exactly the kind of thing "checked" documentation is supposed to prevent.

## Reproduction log

Environment: Windows 10, Python 3.11.9, Git Bash (`bash` tool) for `reproduce.sh`,
MiKTeX (`pdflatex`, `bibtex` both on PATH). Repo already had a `.venv/` and populated
`cases/*/*.png` + `provenance/numbers.json` present before I started (from a prior run by
the project's own author, not by me) — `checks/independent_review/` itself was empty.

1. Read `README.md` top to bottom — clear, one command to reproduce everything, repo map
   accurate.
2. `bash reproduce.sh` (`all` = checks + cases). Completed in well under a minute.
   - venv already existed so it reused it; `pip install -r requirements.txt` succeeded
     quietly (numpy 1.26.4, scipy 1.12.0, matplotlib 3.8.3, mpmath 1.3.0).
   - All 5 `tests/test_*.py` printed `ALL CHECKS PASSED` (26 individual `[PASS]` lines,
     0 `[FAIL]`). No warnings, no exceptions.
   - `cases/case_A_chirp/run.py` and `cases/case_B_monochromatic/run.py` both ran,
     printed their `numbers.json` dict to stdout, and wrote all 8 figures (4 each).
   - No workaround needed anywhere in this step.
3. `bash reproduce.sh theory`: `pdflatex` → `bibtex` → `pdflatex` ×2 on `theory.tex`.
   Completed with only routine warnings (hyperref "token not allowed in PDF string",
   underfull hbox, a duplicate-destination warning from `\tableofcontents`) — no errors.
   Produced `theory/theory.pdf`, 18 pages, matching the version already in the repo.
4. `bash reproduce.sh report`: same pattern on `report.tex`, needs the case figures to
   already exist (they did, from step 2). Completed cleanly, 6 pages, `report/report.pdf`.
5. Compared every number in both `cases/*/RESULTS.md` against the corresponding
   `provenance/numbers.json` by hand (see Findings) — all matched to the stated
   precision except one rounding slip (Finding 4, minor).
6. Read `report/report.html`/`report/report.pdf` and grepped the key numbers
   (61.9, 777.6, 15.18, 1.038, 1.056, 0.309, 1.385, 1.589) out of `report.html` —
   all present and matching `numbers.json`.
7. Read `theory/theory.pdf` in full (18 pages) and checked its Table 4.1 (paraxial
   validity numbers) against `theory/paraxial_validity_numbers.json` — all six values
   matched to the table's stated precision.

No step required any intervention beyond what `README.md` and `reproduce.sh` already
document. This is a genuinely reproducible repo.

## Findings

### Physics

- **None found.** Spot-checked: `F_point_lens` (Deguchi & Watson 1986 / Takahashi &
  Nakamura 2003 closed form, `src/gwlens/waveoptics.py:183`) against the formula quoted
  in `theory.tex` Eq. (3.4) — identical. `magnification(x) = x^4/(x^4-1)`
  (`waveoptics.py:49`) against the standard point-lens image magnification
  `[1-(θ_E/θ)^4]^-1` — algebraically identical, and it's cross-checked at runtime against
  the independent Paczynski (1986) total-magnification formula
  (`check_paczynski_magnification`, relative error 2.13e-15). `isco_frequency`
  (`chirp.py:27`) reduces to the standard `f_isco ≈ 4400 Hz × (M_sun/M)` rule of thumb
  (verified by hand: 4400.2 Hz at M=1 Msun) and the test confirms 73.29 Hz vs. 73.33 Hz
  expected at 60 Msun. The TaylorF2 2PN phase coefficients
  (`taylorf2.py:59-61`, `φ2=3715/756+55η/9`, `φ3=-16π`,
  `φ4=15293365/508032+27145η/504+3085η²/72`) match the standard literature values I know
  (Buonanno, Iyer, Ochsner, Pan & Sathyaprakash 2009, as cited). Constants in `units.py`
  (G, c, GMsun, pc, AU) match CODATA 2018 / IAU 2015 / IAU 2012 reference values to the
  precision quoted. The paraxial-validity numbers in `theory.tex` Table 4.1 (Chapter 4.3)
  all trace correctly to `theory/paraxial_validity_numbers.json`.

- **[minor] Stale docstring number in `src/gwlens/waveoptics.py::F_hybrid`
  (lines 206-221).** The function signature default is `w_geo_threshold=30.0`, and every
  other place in the repo (tests, RESULTS.md, theory.tex, report) consistently says 30.
  But the docstring prose itself says: *"`w_geo_threshold=20` sits comfortably inside
  that validated agreement."* This is a leftover from an earlier version of the code
  (the threshold was apparently moved from 20 to 30 without updating this one sentence).
  It has zero effect on any actual result — the executable default and every consumer
  agree on 30 — but it is exactly the kind of inconsistency that erodes trust in
  "this code is self-documenting" claims, and a careful reader/grader diffing the
  docstring against the signature would flag it immediately.

- **[minor] Stale comment in `tests/test_system.py::check_quasi_monochromatic_regime`
  (line 45-46).** The comment reads *"not over the full remaining time to merger (an
  earlier version of this check used tau_B/2..."* and, critically, says
  *"(system.T_OBS_B_S = 8 outer periods)"*. But `src/gwlens/system.py:51` defines
  `T_OBS_B_S = 6.0 * P_OUT_S` with the comment `# Case B observing baseline: 6 outer
  periods`, and every RESULTS.md/report/theory.tex statement says 6 outer periods (24
  days), which is also what actually printed when I ran it
  (`n_outer_periods_observed: 6.0`). The "8" in the test comment is simply wrong —
  another leftover from an earlier parameter choice, never touching the actual
  numbers used or reported, but a genuine internal inconsistency a grader diffing
  comments against code would catch.

### Internal consistency (numbers across documents)

- **[nitpick] Case B `RESULTS.md` rounds `f_end_Hz` incorrectly.**
  `cases/case_B_monochromatic/RESULTS.md` states the frequency drifts "0.0500 -> 0.0521
  Hz", but `provenance/numbers.json`'s `f_end_Hz = 0.05200922081649974` rounds to
  **0.0520**, not 0.0521 (verified: `f'{0.05200922081649974:.4f}'` → `'0.0520'`). A
  one-digit rounding slip in the prose, not in the underlying number — the 4.0% drift
  figure quoted alongside it is correct (`fractional_freq_drift = 0.0402`, rounds to
  4.0%). No figure or downstream conclusion depends on the wrong digit.

- **[nitpick] `theory.tex` §3.3 quotes "3.5% agreement at w=30" for the hybrid-evaluator
  threshold check, but the test it cites
  (`tests/test_waveoptics.py::check_hybrid_matches_at_threshold`) actually prints
  `3.10e-2` (3.1%) when run.** Both are under the test's own 4% tolerance and the
  qualitative claim ("agrees well at the switch point") is correct either way, but the
  theory document's number does not exactly match what the check it cites currently
  produces — likely a rounding-up-for-safety choice or a stale number from when the
  check's `y` sample set was different, rather than a wrong claim.

- Everything else I checked agreed exactly: `y_A=1.589` (Case A) appears identically in
  `RESULTS.md`, `report.html`/`report.pdf`, and `numbers.json`; Case B's independently
  found `min_y_during_observation=1.58900...` agreeing with Case A's `y_A=1.58895...` to
  4 significant figures is correctly described as agreement "to 4 significant figures,"
  not overstated as exact; `max_abs_F²=1.385` in the figure caption and RESULTS.md
  correctly equals `max_abs_F² = 1.176655857549728² = 1.3845...`; the RMS-vs-sqrt(mean
  |F|²) values in Case B differ at the 5th significant figure (1.0030 vs 1.0030 to 4
  figures, 1.00304 vs 1.00305 at 5) — consistent with the two being different statistics
  (RMS of a time-domain sample vs. sqrt of a frequency/orbit-averaged mean-square), not
  an error.

### Unsupported claims

- **None found.** Every quoted number in both `RESULTS.md` files traces to
  `provenance/numbers.json`, written by the same script, and every qualitative claim I
  checked (repeated-lensing origin, Einstein-ring diffraction-spike regularization,
  geometric-optics limit validity range, paraxial-approximation breakdown window) has an
  explicit pointer to a `tests/` check or a cited reference, matching the discipline
  `wiki/conventions.md` states for itself ("nothing in RESULTS.md may state a number
  that is not also in some numbers.json").

### Scope honesty

- Good. The report's "Scope and limitations" section (`report/report.pdf` p.3) and
  `theory.tex` §4.3 both explicitly flag: point-mass lens only, TaylorF2 truncated at
  2PN, no absolute strain calibration, no noise/antenna-pattern/matched-filtering (stated
  as illustrative), and — most substantively — that the paraxial/thin-lens approximation
  is quantifiably pushed near each Case B front/back crossing, with an actual computed
  number (~60s window, 1.7e-4 of the outer period) rather than a vague caveat. This is
  the strongest part of the project's rigor: it doesn't just assert an approximation is
  fine, it computes how close to the edge it is.
- `wiki/log.md` documents three real bugs caught and fixed during development (a
  frequency-domain-vs-time-domain amplitude scaling mixup, a Fourier-sign error in the
  TaylorF2 phase, a geometric-optics phase-normalization omission) with enough technical
  detail to verify each account against the current code — I checked the sign-convention
  fix by reading `taylorf2.py::spa_phase`'s docstring and it is consistent with
  `wiki/conventions.md`'s stated Fourier convention.

### Anything a grader would flag

- The two stale-comment/docstring issues above (Findings under "Physics") are the kind
  of thing that would cost a small amount of credit under a strict "does the
  documentation actually match the code" rubric, even though neither affects any
  reported result.
- I could not find anything resembling an inflated claim, a missing citation, or a
  result presented without a check — the provenance discipline (`claims.yaml` +
  `numbers.json` + pointer to a `tests/` function for every stated number) is followed
  consistently in both cases I read in full.

## What I did NOT have time/ability to check

- I did not re-derive the wave equation → Helmholtz reduction (`theory.tex` Chapter 2)
  from first principles myself; I only checked that the stated result (identical scalar
  amplification factor for both GW polarizations) is consistent with the citations given
  (Isaacson 1968; Peters 1974; Takahashi & Nakamura 2003) and is used consistently
  downstream. A reviewer with more GR background should re-derive Eq. (2.1)→(2.2) by
  hand.
- I did not independently re-run the `F_bruteforce_2d` / `F_radial_1d` /
  `F_point_lens` three-way comparison at parameter values outside what `tests/` already
  covers — I trust the test's own tolerances and cited convergence behavior
  (`O(eta)` scaling verified in the test output: ratios 1.98, 1.99) rather than
  re-deriving the Bessel-reduction integral myself.
- I did not check the TaylorF2 2PN coefficients against a second independent source
  beyond my own recollection of the literature form cited (Buonanno et al. 2009); I did
  not have network access to pull the paper and diff coefficients character-by-character.
- I did not evaluate whether the choice of `M_L=5e4 Msun`, `i_out=87°`, `D_L=5 kpc`, etc.
  is astrophysically well-motivated beyond what the project itself claims (it explicitly
  disclaims astrophysical typicality and states the numbers were chosen so `w` is O(1) at
  both epochs) — I only checked that the stated regime consequences (static lens, 4.4e-5
  duty cycle; quasi-monochromatic, 4% drift over 6 periods) follow correctly from those
  choices, which they do.
- I did not open every `.png` figure pixel-by-pixel against its caption's numeric claims
  (e.g., I did not measure "0.79 to 1.27" off `caseA_F_of_f.png` myself) — I visually
  inspected the rendered report figures and they are qualitatively consistent with what
  the text describes, but did not extract pixel data to verify the exact bounds.
- I have no way to independently verify authorship/timeline claims in `wiki/log.md`
  (e.g., that the described bugs were actually caught in the order and manner described)
  — I only checked that the current code state is consistent with the log's account of
  what was fixed and why.
- LaTeX build: I did not check the PDF's bibliography formatting/citation completeness
  beyond visual inspection of the rendered `References`/`Bibliography` pages.
