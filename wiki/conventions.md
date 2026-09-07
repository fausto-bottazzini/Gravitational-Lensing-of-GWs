# Conventions

The notation contract for this project. If code, `theory.tex` or `report` ever
disagree with this page, this page wins until corrected here first. Settled 2026-09-07.

## Physical scenario

One hierarchical triple, followed across its own evolution, not two unrelated
examples:

- **Inner binary** (the GW *source*): two compact objects, masses `m1, m2`,
  total mass `Mb = m1+m2`, that inspiral and eventually merge.
- **Outer body** (the GW *lens*): a third, more massive point mass `M_L`
  (IMBH range; `system.py` uses `M_LENS_MSUN = 5e4` — chosen in `wiki/log.md`
  so that the wave-optics parameter `w` is O(1) at both frequencies used
  below, not to be astrophysically "typical").
- Hierarchy: `a_in << a_out` (inner semi-major axis vs. the outer orbit of the
  binary's center of mass around `M_L`). This is standard triple-stability
  hierarchy and is *also* the condition that lets us treat the binary as a
  single point source for lensing (both members deflected identically — no
  differential/microlensing across the binary itself).
- Two snapshots of the SAME system:
  - **Case A ("chirp")**: late inspiral, `Mb` merges in the LIGO band. Merger
    duration (seconds) `<<` outer orbital period `P_out` (years+) ⇒ lens
    position frozen during the whole observed chirp ⇒ **static lens**, single
    impact parameter `y`.
  - **Case B ("quasi-monochromatic")**: early inspiral, `Mb` wide, GW frequency
    in the LISA band and nearly constant over the observation. Observation
    baseline `T_obs` is now comparable to or longer than `P_out` ⇒ the outer
    orbit **is** visible: `y(t)` sweeps an ellipse and periodically modulates
    the lensing pattern.
- We work in the regime where cosmological redshift is irrelevant (the triple
  is Galactic-scale), so `D_L`, `D_S`, `D_LS = D_S - D_L` are ordinary Euclidean
  distances and the lens redshift `z_L = 0` (drops out of the usual
  `M_Lz = (1+z_L) M_L` factor in the lensing literature).

## Units

- Derivations: geometrized units `G = c = 1` throughout `theory.tex`, restored
  explicitly whenever a formula is first written in a form meant to be
  evaluated numerically.
- Code and numbers: SI-adjacent astro units — solar masses (`Msun`), seconds,
  Hz, Mpc/pc for distances. Constants pinned once in `src/gwlens/units.py`
  from CODATA 2018 (`G`, `c`) and the IAU 2015 nominal solar mass parameter
  (`GMsun`); no `astropy` dependency (kept out to keep the environment small —
  logged as a choice, not an oversight).
- `Msun` in geometrized *length* units: `GMsun/c^2 = 1.4766...  km`; in
  geometrized *time* units: `GMsun/c^3 = 4.9255...  us`. Both derived in code
  from the pinned `G, c, GMsun`, never hand-typed twice.

## Fourier convention

`h(f) = ∫ h(t) e^{-2*pi*i*f*t} dt`, `h(t) = ∫ h(f) e^{+2*pi*i*f*t} df` — same
sign as `numpy.fft.rfft`/`irfft`, which is why the case scripts can multiply
by `F(f)` and inverse-FFT directly with no extra bookkeeping.

**A later-arriving image gets a `-i*2*pi*f*Delta_t` phase in THIS
convention** (standard delay theorem: `h(t-Delta_t)` transforms to
`e^{-2*pi*i*f*Delta_t} H(f)` when the forward transform carries
`e^{-2*pi*i*f*t}`). Most of the cited lensing literature (Takahashi &
Nakamura 2003 included) uses the *opposite* sign convention: a formula
copied from a paper needs its overall `i` flipped (conjugated) before it
lands in this repo, checked with an actual time-domain causality test
(`tests/test_waveoptics.py::check_causality_of_lensed_pulse`) — agreement
between two evaluators is not enough to catch a flip, since `|F(w,y)|` is
conjugation-invariant and would look identical either way. `wiki/log.md`
has the full account of the two times this repo got it wrong regardless.

## Reference-phase convention

`F(w,y)` as computed here is referenced to the strong (minimum-time)
image's own arrival: it carries zero phase there by construction, so only
the physical, convention-independent delay `Delta_T(y)` between images
shows up. This matters because `psi(x)=ln|x|` (below) fixes the deflection
potential only up to an arbitrary additive constant, which otherwise
leaks into `F`'s ABSOLUTE phase — never into `|F|` or `Delta_T`, but very
much into any absolute lensed-vs-unlensed timing read off a time-domain
reconstruction (Case A's merger peak, its echo) if this reference isn't
applied. `wiki/log.md` has the account of the one time this repo got that
wrong and mistook the resulting rigid time-shift for interference.

## Lensing / wave-optics notation

| symbol | meaning |
|---|---|
| `D_L, D_S, D_LS` | observer-lens, observer-source, lens-source distances |
| `xi_0` | length scale on the lens plane (fixed to the Einstein radius `R_E`) |
| `eta_0 = xi_0 D_S/D_L` | corresponding scale on the source plane |
| `x = xi/xi_0` | dimensionless lens-plane position (2-vector) |
| `y = eta/eta_0 = beta/theta_E` | dimensionless source position / impact parameter (2-vector; `y=|y|` for a point lens by symmetry). Two equivalent forms of the same quantity: `eta` is the physical source-plane offset, `eta_0` its scale above; `beta` is the same offset as a true angular position on the sky and `theta_E` the angular Einstein radius (`theta_E = xi_0/D_L`) -- substituting `eta=beta*D_S`, `eta_0=xi_0*D_S/D_L` shows these are identical, not two different `y`s. |
| `psi(x)` | (dimensionless) deflection potential, `psi(x) = ln|x|` for a point lens -- defined only up to an arbitrary additive constant (any solution of the underlying Poisson equation is; see "Reference-phase convention" above for where that constant matters and where it cancels) |
| `t_d(x,y)` | time delay (geometric + Shapiro) of the ray through `x`, relative to no lens |
| `F(f,y)` | (dimensionless, complex) amplification factor: `h_lensed(f) = F(f,y) * h_unlensed(f)` |
| `w` | dimensionless frequency, `w = 8*pi*G*M_L*f/c^3` (i.e. `w = 4 G M_L omega /c^3`, `omega=2*pi*f`) |

`F` depends on the lens only through `w` and `y` for a point lens (no explicit
`M_L`, `D_L`, `D_S`, `D_LS` separately) — this collapse is itself one of the
first things checked in `tests/`.

## Scalar-field reduction

For wavelengths long compared to the lens's Schwarzschild radius but the lens
still treated in the geometric-optics limit of *its own* background spacetime
(i.e. thin-lens, weak-field lens — always true here, `R_E >> R_Sch(M_L)`), each
GW polarization amplitude obeys the same scalar wave (Helmholtz) equation as a
massless scalar field propagating on that background, and both polarizations
pick up the *identical* scalar `F(f,y)` — lensing does not mix or rotate
polarizations in this regime. This is derived in `theory.tex` Sec. 2, following
Peters (1974) and Takahashi & Nakamura (2003); it is the reason this whole
project can legitimately stop at "scalar field" and call the lensing problem
solved for the GW case.

## Orbit conventions

- Outer orbit: Keplerian, elements `(a_out, e_out, i_out, Omega_out, omega_out,
  t_0)` in the standard astronomical sense; `i_out` measured from the plane of
  the sky (`i_out=0` face-on ⇒ no lensing modulation from orbital motion in the
  y-direction alone — a deliberate degenerate case used as a check).
- Inner binary phase/frequency evolution: Case A's unlensed waveform is
  IMRPhenomD (Khan et al. 2016, via `pycbc`/LALSimulation where available --
  WSL2/Linux, see `README.md`), a full inspiral-merger-ringdown NR-calibrated
  model; `src/gwlens/taylorf2.py`'s leading-order (Newtonian/quadrupole,
  restricted) 2PN inspiral-only `df/dt` is the automatic fallback when
  `pycbc` is not importable (native Windows), and is what
  `src/gwlens/chirp.py` and Case B's quasi-monochromatic treatment use
  throughout (no merger is ever reached there) — logged as a choice, see
  `wiki/log.md`.

## File-naming / provenance pattern

Followed from `GW-AI-course/day5/exercise/{b,c}` (course scaffold used with
permission of its structure, not its content — no solved exercise from the
course is copied into this repo):

- Each `cases/case_*/` carries `run.py`, `RESULTS.md`, and
  `provenance/{claims.yaml,numbers.json}`.
- `claims.yaml`: one entry per stated result, each with `evidence` pointing at
  `file.py::function` and/or a literature reference, plus which `numbers` in
  `numbers.json` back it.
- Nothing in `RESULTS.md` may state a number that is not also in some
  `numbers.json`, written by the script that computed it (never hand-typed).
