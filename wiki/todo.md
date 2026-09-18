# Open items

What is still open. Everything that is done and why it was done that way is in
`log.md`; this page only lists work that is not.

## Pending

The folder-by-folder tidy-up pass is finished: `theory/`, `src/`+`tests/`, the
case figures, `report/` (both the deck and the written report), the provenance
and `RESULTS.md` files, `wiki/`, `README.md` and `reproduce.sh` have each had
one, and the two items that used to live here — the validity numbers reaching
the report, and the README's stale description of `tests/` — are done.

What follows is what the six-reviewer round of 2026-09-18 (`log.md`) turned up
and this pass did **not** act on. Each was checked and is real; each was left
because it is presentation work rather than correctness, or because the fix is
larger than the defect. The measurements in brackets are that round's, not
re-verified here.

### `report/report.html`

- [ ] The ring rendering on slide 9 keeps full per-pixel resolution down to
  `per > 2*pix`, which is exactly Nyquist rather than the comfortable margin
  its comment claims. In the middle of the sweep that draws beat structure
  instead of either the rings or their mean; `per > 4*pix` would close the
  band. [Around w=139: 132 luminance extrema measured across 160 px against
  175 predicted — more than the pixels can represent.]
- [ ] The schematic's labels ("binaria interna", "i = 87°", "lente 5x10^4
  M_sun") superimpose at conjunction, which is the one orbital phase the
  panel exists to be read at.
- [ ] Slide 8's `arg F(f)` axis carries no tick values, beside text that says
  "397 vueltas". The panel is a 10-13 Hz zoom showing about eleven of them,
  and nothing on it says so.
- [ ] The slide-4 schematic runs on its own idealized clock — conjunction at
  day 2.0 of 4 — while its readout reads "día X de 4" and slide 12, which is
  on the real ephemeris, puts conjunction at day 1.0006.
- [ ] `history.replaceState` writes `#N` on every slide change but nothing
  listens for `hashchange`, so a deep link works on first open and Back,
  Forward or a hand-edited hash silently do nothing.

### `cases/case_B_monochromatic/caseB_one_period.html`

- [ ] Wheel-zoom and drag on the detail canvas divide by the full canvas width
  instead of the plotted region, so the point under the cursor does not stay
  put. `genCv`'s click handler already corrects for the margins, so this is an
  oversight rather than a convention. [~2300 s of drift over 30 notches
  anchored at the right edge of the plot area.]
- [ ] The readout calls a full turn "17416 ciclos", from the instantaneous
  carrier at conjunction, where provenance's `one_period_carrier_cycles` is
  17280, from `f = 0.05` Hz exactly. Both are defensible; they are the same
  quantity 0.8% apart, and only one of them is in `numbers.json`.

### `src/` and `cases/`

- [ ] A two-image reconstruction check. `RESULTS.md` used to assert that
  rebuilding `h_lensed` from nothing but two delayed, rescaled copies of the
  unlensed waveform agrees to `1.1e-13`; no code in the repo produced that
  number, so the sentence was removed rather than left asserted. The check is
  worth having for real — it is the one end-to-end statement that Case A is
  pure geometric optics, as opposed to two checks of the switch that reaches
  it.
- [ ] `F_hybrid(w, np.inf)` returns NaN with only a numpy warning.
  `geometry.impact_parameter_of_time` returns `inf` on the unlensed half by
  design and the `F = 1` substitution is the caller's job; no caller forgets
  it today, but if one did the failure would be silent.
- [ ] Dead code in Case A's `run.py`: `h_diff` and `zmask_idx` are computed and
  never used, and the twenty-line comment above `h_diff` still describes the
  difference panel that the two ratio panels replaced.
- [ ] `w_B` — Case B's headline number, quoted in nearly every claim statement
  — appears in no claim's `numbers:` list, and Case A leaves thirteen keys
  uncited the same way, several of them quoted in `RESULTS.md`.

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
- [ ] Making slides 8, 9, 12 and 14 of `report.html` fit a short viewport.
  They grow past it and the page scrolls rather than the slide, despite the
  `overflow-y` the CSS sets — `.slide` has `min-height:100vh` and no maximum.
  [At 1024x768: slide 8 +80 px, 9 +333, 12 +301, 14 +370.] Left as it is on
  purpose: the deck is not projected, the browser scrolls, and making the
  four fit would mean cutting content to win back a scrollbar.
