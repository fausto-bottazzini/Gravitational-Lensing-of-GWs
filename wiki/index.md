# Index

Read this first. Everything below is reachable from here.

## The project

Gravitational lensing of gravitational waves by a hierarchical triple, in wave
optics, treated as a scalar field. One system, two epochs: a chirp (static
lens) and a quasi-monochromatic wide binary (moving lens). See
[[conventions]] for the exact setup and notation, [[log]] for how we got here,
[[todo]] for what is left.

## Canonical owners

| topic | canonical page / file |
|---|---|
| notation, units, physical scenario | `wiki/conventions.md` |
| decision history (every bug caught, and how) | `wiki/log.md` |
| open items | `wiki/todo.md` |
| full derivation (the "book chapter") | `theory/theory.tex` / `theory/theory.pdf` |
| paraxial/thin-lens validity, quantified | `theory/paraxial_validity.py` (theory.tex Sec. 4.3) |
| wave-optics core code (point lens `F(w,y)`) | `src/gwlens/waveoptics.py` |
| the one pinned system (m1, m2, M_L, orbit) | `src/gwlens/system.py` |
| inner-binary waveforms | `src/gwlens/chirp.py` (leading order), `src/gwlens/taylorf2.py` (2PN inspiral-only, fallback), `src/gwlens/imr_waveform.py` (picks pycbc's IMRPhenomD when available, else the fallback) |
| outer-orbit / D_LS(t) geometry | `src/gwlens/geometry.py` |
| correctness checks | `tests/` |
| Case A (chirp, static lens) | `cases/case_A_chirp/` |
| Case B (monochromatic, moving lens) | `cases/case_B_monochromatic/` |
| final presentation | `report/report.html`, `report/report.pdf` |
| how to reproduce all of it | `README.md`, `reproduce.sh` |
