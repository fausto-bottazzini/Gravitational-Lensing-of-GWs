# theory.tex outline

A textbook-chapter-style derivation, matching the numbering used in
`report/report.tex`'s references back into it. Sections map 1:1 onto
`wiki/conventions.md`.

1. **Setup and scope.** The physical picture (hierarchical triple, source =
   inner binary, lens = outer body), the hierarchy `a_in<<a_out`, why
   cosmological redshift is dropped, the two epochs.
2. **From the wave equation to a scalar field.** Linearized GW propagation on
   a weak-field lens background; each polarization amplitude obeys the same
   scalar wave (Helmholtz) equation; statement (with reference, not full
   proof) that lensing does not mix polarizations here. This is the
   "reaching scalar field" milestone the user asked for explicitly.
3. **The diffraction integral.** Kirchhoff/Fresnel diffraction from the
   Helmholtz equation, thin-lens approximation, the dimensionless (w,y,x)
   reduction, definition of F(w,y).
4. **Point-mass lens.** Fermat potential, image positions, magnifications
   (Jacobian derivation), time delay; geometric-optics limit as a
   stationary-phase evaluation of the same integral (Morse-index argument);
   the closed-form confluent-hypergeometric F(w,y); the Bessel-reduced 1D
   integral; three-way numerical cross-check (mirrors tests/test_waveoptics.py,
   with the actual numbers and the sign-error story from wiki/log.md as a
   worked "how do you actually verify this" example).
5. **The hierarchical triple as a lens system.** Outer Keplerian orbit,
   rotation into the observer frame, the D_LS(t)=z_los(t) subtlety and the
   front/back "repeated lensing" split (derivation, not just code comments).
6. **The inner binary as a source.** Leading-order quadrupole chirp,
   restricted amplitude, why full PN is out of scope.
7. **Putting it together: two regimes.** Static-lens vs quasi-monochromatic
   conditions derived from the timescale hierarchy; how Case A and Case B in
   `cases/` instantiate them.
8. **What a detector would see.** The idealized simple-detector model used
   in `src/gwlens/detector.py` and its limitations, stated explicitly.

References (BibTeX in theory/refs.bib, shared with report/report.tex):
Takahashi & Nakamura 2003; Deguchi & Watson 1986; Peters 1964; Peters &
Mathews 1963; Paczynski 1986; Schneider, Ehlers & Falco 1992; D'Orazio &
Loeb 2020; Ulmer & Goodman 1995; Maggiore "Gravitational Waves" Vol. 1.
