# Todo / open items

- [ ] Pin exact `(m1, m2, M_L, a_out, e_out, i_out)` numbers — done in
  `cases/*/run.py` directly (single source of truth), not duplicated here.
- [ ] `theory.tex`: full derivation, Sections 1-6 (see `theory/OUTLINE.md`).
- [ ] `src/gwlens`: wave optics core, orbits, waveform, detector toy model.
- [ ] Three independent evaluators of the point-lens `F(w,y)` must agree
  (brute-force 2D quadrature, 1D Bessel-reduced radial integral, closed-form
  confluent hypergeometric via `mpmath`) — this is the main correctness check
  for the whole project and gates everything downstream.
- [x] Case A figures: lens-plane diffraction/interference pattern (extended,
  2D), F(f) vs f across the chirp, lensed vs unlensed strain, simple detector
  view. Done -- see cases/case_A_chirp/RESULTS.md.
- [x] Case B figures: F(w_B,y(t)) with orbital modulation (repeated lensing
  pulses), the fixed pattern with the orbit track overlaid, idealized
  detector view. Done -- see cases/case_B_monochromatic/RESULTS.md. (A
  literal "movie" was dropped in favor of the static pattern+track figure,
  which shows the same information without needing a video file in the
  repo -- logged as a scope choice, not an oversight.)
- [ ] `checks/independent_review`: fresh subagent, no prior context, told to
  reproduce and critique — run near the end, log its verdict here whether it
  passes or not.
- [ ] `report/` HTML + PDF assembled from the case RESULTS.md + figures.
- [ ] `README.md` at repo root with a `reproduce` path a stranger can run.
- [ ] requirements.txt pinned; test reproduce from a clean venv if time
  allows.
- [ ] Final commit; leave push-to-GitHub as the one step for the user.
