# Todo / open items

- [ ] Pin exact `(m1, m2, M_L, a_out, e_out, i_out)` numbers — done in
  `cases/*/run.py` directly (single source of truth), not duplicated here.
- [ ] `theory.tex`: full derivation, Sections 1-6 (see `theory/OUTLINE.md`).
- [ ] `src/gwlens`: wave optics core, orbits, waveform, detector toy model.
- [ ] Three independent evaluators of the point-lens `F(w,y)` must agree
  (brute-force 2D quadrature, 1D Bessel-reduced radial integral, closed-form
  confluent hypergeometric via `mpmath`) — this is the main correctness check
  for the whole project and gates everything downstream.
- [ ] Case A figures: lens-plane diffraction/interference pattern (extended,
  2D), F(f) vs f across the chirp, lensed vs unlensed strain, simple detector
  view.
- [ ] Case B figures: same but F(f,t) with orbital modulation, a
  time-sequence/movie of the moving diffraction pattern, spectrogram showing
  sidebands at the outer orbital frequency.
- [ ] `checks/independent_review`: fresh subagent, no prior context, told to
  reproduce and critique — run near the end, log its verdict here whether it
  passes or not.
- [ ] `report/` HTML + PDF assembled from the case RESULTS.md + figures.
- [ ] `README.md` at repo root with a `reproduce` path a stranger can run.
- [ ] requirements.txt pinned; test reproduce from a clean venv if time
  allows.
- [ ] Final commit; leave push-to-GitHub as the one step for the user.
