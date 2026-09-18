# Open items

What is still open. Everything that is done and why it was done that way is in
`log.md`; this page only lists work that is not.

## Pending

- [ ] `README.md` describes `tests/` as "every correctness check referenced
  from `theory.tex` and the case `RESULTS.md` files", and the theory text no
  longer references any test. To be fixed when `README.md` gets its own pass,
  which is the last stage of the tidy-up. That same pass should say that the
  cheap way to read the theory is `theory/theory.tex` and not `theory.pdf`
  (`theory/OUTLINE.md` explains why; a `theory.md` was created and deleted on
  the strength of exactly that reasoning).

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
  3.5e-3 at the worst point of the lensed half, and 2.5e-5 / 5.6e-4 cycles
  for `arg F`, at points where nothing is read off.
- [ ] A genuinely non-paraxial treatment of the ~60 s window around each
  `D_LS = 0` crossing (`theory.pdf` §4.5). Quantified as not worth the
  complexity, since `F -> 1` there anyway, but an open problem for anyone who
  wants to push further.
