# Todo / open items

- [x] Pin exact `(m1, m2, M_L, a_out, e_out, i_out)` numbers — `src/gwlens/system.py`.
- [x] `theory.tex`: full derivation, `theory/theory.pdf` (23 pages).
- [x] `src/gwlens`: wave optics core, orbits, chirp waveform, TaylorF2, one
  pinned system. No separate "detector toy model" module ended up needed --
  the idealized detector view (Hilbert-transform envelope, no noise/antenna
  pattern) is a few lines in each case script; a whole module for it would
  have been overkill for what was actually shown.
- [x] Three independent evaluators of the point-lens `F(w,y)` agree
  (brute-force 2D quadrature, 1D Bessel-reduced radial integral, closed-form
  confluent hypergeometric via `mpmath`) — `tests/test_waveoptics.py`.
- [x] Case A figures — `cases/case_A_chirp/RESULTS.md`.
- [x] Case B figures — `cases/case_B_monochromatic/RESULTS.md`.
- [x] Three fresh-agent independent reviews, run with no prior context each
  time, findings fixed and logged in `wiki/log.md`. The review reports
  themselves (`checks/independent_review/`) were working documents for
  this development process, not part of the deliverable, and are kept
  locally rather than in the repo.
- [x] `report/` HTML + PDF assembled from the case RESULTS.md + figures.
- [x] `README.md` at repo root with a `reproduce` path a stranger can run.
- [x] `requirements.txt` pinned; `reproduce.sh` tested from a clean venv --
  found and fixed two real bugs this way: `pip install --upgrade pip` fails
  on Windows when pip is running, and (more serious, caught by an
  independent review) a bare `python3` call after `source
  .venv/Scripts/activate` on Windows silently ran the SYSTEM Python instead
  of the venv's (Windows venvs only ever create `python.exe`, never
  `python3.exe`/`python3`), so every `pip install`/script invocation was
  installing into and running from outside the venv the script had just
  created, defeating the isolation `reproduce.sh` claimed to give. Fixed by
  locating the venv's own interpreter explicitly instead of relying on
  `PATH` after activation; see `wiki/log.md`.

## What's left, if there's time later (none of this blocks the deliverable)

- [ ] Higher PN order (2.5-3.5PN) for `taylorf2.py`, IF a verified reference
  implementation becomes available to check the extra coefficients against
  (deliberately not attempted from memory alone -- see `wiki/log.md`).
- [ ] A real non-paraxial treatment of the ~60s window per Case B crossing
  (`theory.tex` Sec. 4.3) -- quantified as not worth the complexity given
  what it would and wouldn't change, but a genuinely open problem if
  someone wanted to push further.
- [x] `pycbc`/LALSimulation as the real inspiral-merger-ringdown waveform
  (IMRPhenomD, Khan et al. 2016) for Case A -- native Windows install still
  hangs indefinitely resolving `lalsuite`, but pycbc installs cleanly under
  WSL2, which is what the professor's own machine runs anyway; see
  `wiki/log.md` for how this was found and fixed.
  `taylorf2.py`'s 2PN inspiral-only waveform remains the automatic fallback
  when pycbc isn't importable, so `reproduce.sh` still runs end to end with
  no manual steps on any platform -- only the WSL path reproduces the exact
  committed Case A merger/ringdown/echo figures. WSL setup documented in
  `README.md`.
- [x] An interactive visualization of Case B's (and Case A's) orbit,
  interference/diffraction pattern, and idealized detector view, synced to
  a scrub control -- built as `report/report.html` (vanilla JS + canvas,
  data inlined from `caseX_animation_data.json`, no video file or chart
  library in the repo); see `wiki/log.md`.
