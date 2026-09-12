# Conventions

The project's notation and units contract, in short form.

`theory/theory.pdf` is the source of truth for every *physical* convention:
signature, units, the Fourier sign, the reference phase, the meaning of `w`
and `y`. This page states them compactly so they can be checked without
opening the text, and adds the conventions that are the repository's own —
constants, file layout, provenance — which the text does not cover. **If this
page and `theory.pdf` ever disagree, `theory.pdf` wins and this page is
wrong.**

## The system

One hierarchical triple, followed across its own evolution, not two unrelated
examples. Pinned once in `src/gwlens/system.py`.

- **Inner binary** (the GW *source*): compact objects of masses `m1, m2`,
  total `Mb = m1+m2`, inspiralling towards merger.
- **Outer body** (the GW *lens*): a third, more massive point mass `M_L`,
  in the intermediate-mass range. The value is chosen so that the wave-optics
  parameter `w` is neither extreme at the two epochs studied, not because it
  is astrophysically typical.
- **Hierarchy**: `a_in << a_out`. Theory §1.2 splits this into four separate
  physical requirements, all checked in `tests/test_system.py`: dynamical
  stability (Mardling & Aarseth 2001, margin 93x), the inner binary inside its
  Hill radius about the lens (292x), the outer orbit's own GW decay time
  (4.0e7 yr, of which 1.6e-8 elapses before the inner merger), and Kozai-Lidov
  (t_KL = 16 yr, 24x the time to merger, and further quenched by the inner
  binary's own GR periastron precession, faster by 3.4e4).
- **Two epochs of the same system**:
  - **Case A ("chirp")**: late inspiral. The band of a ground-based detector
    is swept in seconds, against an outer period of days, so the lens is
    frozen at a single impact parameter for the whole signal. `w = 62` to
    `778`: geometric optics.
  - **Case B ("quasi-monochromatic")**: early inspiral, `f_B = 0.05 Hz`, in
    the LISA band and nearly constant over an observation spanning several
    outer periods. `y(t)` sweeps and modulates the pattern. `w_B = 0.31`:
    full diffraction.
- The triple is Galactic, so `D_L`, `D_S`, `D_LS` are ordinary Euclidean
  distances and the lens redshift is `z_L = 0`. Note that the orbital Doppler
  shift plays the same role the `(1+z_L)` factor plays in the cosmological
  literature (theory §4.4).

## Units

- **In `theory.tex`**: `G` and `c` explicit in every formula, so that each
  expression is directly evaluable in SI.
- **In code**: solar masses, seconds, Hz, AU/pc. Constants pinned once in
  `src/gwlens/units.py` from CODATA 2018 (`G`, `c`) and the IAU 2015 nominal
  solar mass parameter (`GMsun`). No `astropy` dependency.
- `GMsun/c^2 = 1.4766... km` and `GMsun/c^3 = 4.9255... us`, both derived in
  code from the pinned constants, never hand-typed twice.

## Fourier convention

```
h(f) = INT h(t) e^{-2*pi*i*f*t} dt        h(t) = INT h(f) e^{+2*pi*i*f*t} df
```

Same sign as `numpy.fft.rfft`/`irfft`, which is why the case scripts can
multiply by `F(f)` and inverse-FFT with no extra bookkeeping.

In this convention a later-arriving image carries a `-2*pi*i*f*Delta_t` phase.
Most of the lensing literature (Takahashi & Nakamura 2003 included) uses the
opposite sign, so **a formula copied from a paper must be conjugated before it
lands in this repo**. This cannot be caught by comparing two evaluators of `F`
against each other, because `|F|` is conjugation-invariant; it needs a
time-domain causality check
(`tests/test_waveoptics.py::check_causality_of_lensed_pulse`).

## Reference phase

`F(w,y)` is referenced to the arrival of the first image, so it carries zero
phase there and only the physical delay `Delta_T(y)` between images survives.
This is needed because `psi(x) = ln|x|` solves the deflection Poisson equation
only up to an additive constant, which would otherwise leak into the absolute
phase of `F` — never into `|F|` or `Delta_T`, but very much into any absolute
lensed-vs-unlensed timing read off a time-domain reconstruction. Derived in
theory §2.5.

## Lensing / wave-optics notation

| symbol | meaning |
|---|---|
| `D_L, D_S, D_LS` | observer-lens, observer-source, lens-source distances. `D_LS` is the *projection of the orbital separation on the line of sight*, not the separation itself (theory §4.2) |
| `xi_0 = R_E` | length scale on the lens plane, fixed to the Einstein radius |
| `eta_0 = xi_0 D_S/D_L` | corresponding scale on the source plane |
| `x = xi/xi_0` | dimensionless lens-plane position (2-vector) |
| `y = eta/eta_0 = beta/theta_E` | dimensionless impact parameter. The two forms are the same quantity: `eta` is the physical source-plane offset and `beta` the same offset as an angle on the sky, with `theta_E = xi_0/D_L` |
| `psi(x)` | dimensionless deflection potential, `ln\|x\|` for a point lens, defined up to an additive constant |
| `t_d(x,y)` | time delay (geometric + Shapiro) of the ray through `x`, relative to no lens |
| `F(w,y)` | dimensionless complex amplification factor, `h_L(f) = F * h(f)` |
| `w` | dimensionless frequency, `w = 8*pi*G*M_L*f/c^3` |

For a point lens `F` depends on the lens only through `w` and `y` — not on
`M_L`, `D_L`, `D_S`, `D_LS` separately. That collapse is itself checked in
`tests/`.

## Scalar-field reduction

In the short-wavelength, weak-field regime each GW polarization amplitude
obeys the same scalar wave equation as a massless scalar field on the lens
background, and both polarizations pick up the *identical* `F`: lensing
neither mixes nor rotates them. This is what lets the whole project stop at a
scalar diffraction problem. Derived in theory §2.2; the error it incurs is
quantified in theory §4.5.

## Source kinematics

The outer orbit does two independent things to the observed waveform, and the
code keeps them in separate modules because they are separate physics:
**optics** (`waveoptics.py`, driven by `geometry.py`'s `y(t)`) and
**kinematics** (`doppler.py`). The second dominates the *phase*: at Case B's
`f_B` the Roemer delay writes 90.5 GW cycles against the lens's 0.017. It does
not touch amplitudes.

For a circular orbit the frequency relation splits exactly:

```
f_obs = f_em / [ (1 + beta_los(t)) * (1 + z_orb) ]
```

- The **first** bracket is applied by evaluating the emitted waveform at the
  retarded emission time solved from `t_obs = t_em + z_src(t_em)/c`. Nothing
  multiplies a Doppler factor anywhere: `dt_em/dt_obs = 1/(1+beta_los(t_em))`
  identically, so the substitution already *is* the shift. Evaluate `beta` at
  `t_em`, never at `t_obs`: they differ by up to `a_out/c ~ 900 s`, over which
  `beta` changes by ~2.7e-4, numerically indistinguishable from `beta^2` and
  easy to mistake for a missing `gamma`.
- The **second** bracket is the gravitational redshift in the lens potential
  plus the transverse Doppler shift, which for a circular geodesic combine
  into `1+z_orb = (1 - 3GM_L/(a_out c^2))^{-1/2} = 1 + 4.08e-4`. Being
  constant it is exactly degenerate with the chirp mass and with `t_c`, so it
  is **reported as a bias, never applied**. It stops being constant for an
  eccentric outer orbit (`todo.md`).
- Do not apply the special-relativistic `D = [gamma(1+beta)]^-1` on top: the
  `1/gamma` is already inside `z_orb`, so it would be double-counted and the
  gravitational piece would still be missing.

Case A uses the same module and applies nothing: over 15 s of a 4-day orbit a
constant delay is reabsorbed by `t_c` and a constant slope by the chirp mass,
leaving only the line-of-sight acceleration, which is 0.0068 rad at `f_isco`.

## Orbit conventions

- Outer orbit: Keplerian, elements `(a_out, e_out, i_out, Omega_out,
  omega_out, t_0)` in the standard astronomical sense, with `i_out` measured
  from the plane of the sky. `i_out = 0` (face-on) gives no lensing at all and
  is used as a degenerate check.
- Inner binary: Case A's unlensed waveform is IMRPhenomD (Khan et al. 2016,
  via `pycbc`/LALSimulation where available — WSL2/Linux, see `README.md`).
  `src/gwlens/taylorf2.py`'s 2PN inspiral-only waveform is the automatic
  fallback when `pycbc` is not importable, and is what `src/gwlens/chirp.py`
  and Case B use throughout, since Case B never reaches merger.

## Provenance pattern

This is what the assignment asks for under "document the provenance of each
result" and "record how each result was verified".

- Each `cases/case_*/` carries `run.py`, `RESULTS.md` and
  `provenance/{claims.yaml,numbers.json}`.
- `claims.yaml`: one entry per stated result, each with `evidence` pointing at
  `file.py::function` and/or a literature reference, plus which entries of
  `numbers.json` back it.
- **No number may appear in a `RESULTS.md` that is not also in some
  `numbers.json`, written by the script that computed it.** Never hand-typed.
