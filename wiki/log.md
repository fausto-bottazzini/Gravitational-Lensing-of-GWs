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
