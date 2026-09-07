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
