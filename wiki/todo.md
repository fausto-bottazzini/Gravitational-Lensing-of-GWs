# Open items

What is still open. Everything that is done and why it was done that way is in
`log.md`; this page only lists work that is not.

## Pending

Nothing. The tidy-up pass is finished: `theory/`, `src/`+`tests/`, the case
figures, `report/` (both the deck and the written report), the provenance and
`RESULTS.md` files, `wiki/`, `README.md` and `reproduce.sh` have each had a
pass, and the two items that used to live here — the validity numbers reaching
the report, and the README's stale description of `tests/` — are done.

What is left is the list below, which is deliberate, plus whatever the next
fresh-agent pass turns up.

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
  evaluating at the nominal `w_B`. Quantified as 1.1e-4 in `|F|` at the pulse,
  3.5e-3 at the worst point of the lensed half, and 2.5e-5 / 5.6e-4 cycles
  for `arg F`, at points where nothing is read off.
- [ ] A genuinely non-paraxial treatment of the ~60 s window around each
  `D_LS = 0` crossing (`theory.pdf` §4.5). Quantified as not worth the
  complexity, since `F -> 1` there anyway, but an open problem for anyone who
  wants to push further.
