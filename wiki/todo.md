# Open items

What is still open. Everything that is done and why it was done that way is in
`log.md`; this page only lists work that is not.

## Pending

- [ ] **The validity numbers belong in `report/`, not in `theory/`.**
  `theory.pdf` states each condition and says the expansion parameter is
  small, which is what a textbook should do; whether it is small for THIS
  system is a statement about the repo's numbers. Three such numbers are
  currently asserted only inside `tests/test_system.py` and need to reach
  the report: the wave-optics point-source ratio `a_in/eta_0` = 6.4e-3
  (30x tighter than the geometric `a_in/a_out` it used to stand in for),
  the adiabatic margin of 5844x where `|F|>1.01`, and the weak-deflection
  bound -- no ray comes closer than the Einstein ring at 121 r_g, where the
  second-order Schwarzschild deflection is 2.4% (1.2% at the y actually
  used).
- [ ] Cross-references to `theory.tex` that moved with the rewrite:
  `report/report.tex:251` cites "§4.2--4.3" for the paraxial numbers, now
  §4.5, and `:305` cites a section by a name that changed;
  `theory/paraxial_validity.py`'s docstring says Sec. 4.3, now 4.5;
  `README.md` describes `tests/` as "every correctness check referenced from
  `theory.tex`", and the text no longer references any test.

## Deliberately not done

Each of these is a scope decision, not an oversight; the reasoning is in
`log.md` under "Stated limits".

- [ ] Higher than 2PN in `taylorf2.py`, absent a reference implementation to
  check the extra coefficients against.
- [ ] Eccentric outer orbit (`e_out > 0`). The lensing pulses would acquire
  the asymmetric shape D'Orazio & Loeb discuss, and — the more interesting
  half — `z_orb` would stop being constant, so it would have to be applied to
  the waveform instead of reported as an `M_chirp` bias.
  `doppler.los_velocity` already carries the eccentric branch and is exercised
  at `e=0.4` in `tests/test_doppler.py`; the redshift side is not.
- [ ] Feeding the observed frequency back into `F` for Case B rather than
  evaluating at the nominal `w_B`. Quantified as 1.0e-4 in `|F|` at the pulse,
  3.5e-3 at the worst point of the lensed half, where nothing is read off.
- [ ] A genuinely non-paraxial treatment of the ~60 s window around each
  `D_LS = 0` crossing (`theory.pdf` §4.5). Quantified as not worth the
  complexity, since `F -> 1` there anyway, but an open problem for anyone who
  wants to push further.
