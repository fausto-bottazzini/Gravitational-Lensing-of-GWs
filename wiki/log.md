# Log

What was decided and why, what broke and what caught it. Internal working
notes, kept so that any claim in the deliverables can be defended by someone
who did not watch them being made. Condensed 2026-09-11 from a much longer
append-only version; nothing here is a result, and every number quoted is
traceable to a `numbers.json` or a test.

## Decisions

**The assignment is in the repo.** `wiki/final-project.html` is a byte-identical
copy of the course brief (`matiaszaldarriaga/GW-AI-course`), verbatim, because
it is the specification this work is measured against and a paraphrase of a
specification is not one. It used to be reachable only as
`../GW-AI-course/final-project.html`, a sibling checkout outside this
repository — so the one document stating the test the repo is held to, *hand it
to a fresh agent that knows nothing*, was the one document that agent could not
read.

**One system, two epochs.** Case A (chirp, static lens) and Case B
(quasi-monochromatic, moving lens) are the *same* hierarchical triple at two
points of one inspiral, not two unrelated examples. Pinned in
`src/gwlens/system.py`.

**Lens mass.** `M_L = 5e4 Msun`, chosen so the wave-optics parameter `w` is
not extreme at either epoch, not because it is an astrophysically typical
third body. It puts Case B at `w_B = 0.31` (full diffraction) and Case A at
`w = 62` to `778` (geometric optics), i.e. on either side of the transition.
Case A was originally hoped to be `w = O(1)` as well; when it turned out not
to be, the case was reframed from "diffraction pattern" to "dense
interference fringes approaching the geometric-optics limit", which is the
better story for a chirp anyway.

**Waveform.** Case A uses IMRPhenomD via `pycbc`/LALSimulation. `pip install
pycbc` hangs indefinitely on native Windows resolving `lalsuite`, but installs
cleanly in WSL2 (~2 min, prebuilt wheels), which is what generated the
committed Case A figures. `src/gwlens/taylorf2.py` (2PN, inspiral-only) is the
automatic fallback when `pycbc` is not importable, so `reproduce.sh` runs
end to end anywhere. Truncated at 2PN deliberately: the 2.5-3.5PN coefficients
carry enough rational and log terms that reproducing them without a reference
implementation to check against was judged a real risk of a silent error.

**No `astropy`.** Constants pinned once in `src/gwlens/units.py` from CODATA
2018 and IAU 2015, to keep the environment small.

**`theory/` is the anchor.** Since the 2026-09-08..11 rewrite, `theory.pdf` is
the source of truth for every physical convention; `wiki/conventions.md`
exposes them in short form and the rest of the repo follows it, not the other
way round.

**Language.** The presented artefacts — `theory.pdf`, `report.pdf`,
`report.html`, and every figure inside them — are in formal Spanish.
`README.md` stays in English, together with everything else that exists for
whoever works on the repo: this wiki, code comments, `RESULTS.md`,
`claims.yaml`.

## Bugs caught, and what caught them

The pattern worth noticing: every bug below that survived longest was one
where the checks compared the code against *itself*.

**The reference phase of the strong image.** `F_geometric_optics` matched
`F_point_lens` in `|F|` at large `w` but not in phase, off by `w*phi_+(y)`.
The closed form carries the minimum image's own Fermat phase, which a naive
stationary-phase writeup drops. Diagnosed from the *ratio* of the two
evaluators. Guarded by `test_waveoptics.py::test_geometric_optics_limit`.

**The global Fourier convention.** The most serious one, and it was found by
an independent review, not by the suite. All four `F(w,y)` evaluators had been
derived self-consistently in Takahashi & Nakamura's convention and then
multiplied into a waveform inverted with `numpy.fft.irfft`, whose convention
is the opposite: Case A's saddle-point image arrived 3.43 s *before* the
merger. Every existing check compared evaluators of `F` against each other, so
all of them shared the error. `|F|` is conjugation-invariant, so nothing in
Case B or in any diffraction figure was affected; the time-domain
reconstruction was. Fixed by conjugating each evaluator. Guarded by
`test_waveoptics.py::check_causality_of_lensed_pulse`, which builds a pulse,
lenses it exactly as the case script does, and checks that the weak image
lands *after* the strong one. That is the only class of check that can catch a
global sign error.

**The same class of error, twice more.** `taylorf2.py`'s SPA phase used the
literature's sign, and the reconstructed chirp peak drifted with the amount of
zero-padding (17 s at 2x, 245 s at 16x) instead of sitting at `t_c`. And
`pycbc`'s frequency-domain waveform is not delivered aligned, so an unshifted
`irfft` put the merger at the very edge of the window. Both fixed with an
explicit frequency-domain time shift; the first is guarded by
`check_ifft_reconstructs_chirp_at_tc`, which checks the peak both lands at
`t_c` and stops moving with padding.

**The outer orbit's kinematics were missing entirely.** A grep for "doppler"
or "roemer" over the repo returned nothing. The source orbits the lens at
`beta_los = 0.0165` with `a_out/c = 907 s`, so the Roemer modulation is 1810 s
peak to peak — **90.5 GW cycles** of phase at Case B's `f_B`, against the
**0.017 cycles** the lens itself writes. Case B's headline, that repeated
lensing pulses are the outer orbit's signature on the waveform, was wrong by
orders of magnitude: the pulses are real and correctly computed, but they are
not the largest thing the orbit does. D'Orazio & Loeb (2020), cited here *for*
repeated lensing, model the Doppler boost alongside it. Three independent
reviews had missed it, each having checked the lensing calculation against
itself rather than asking what else the same orbit does. Added
`src/gwlens/doppler.py`; applied in Case B (pure timing, so no amplitude
result changed and three figures came back byte-identical) and shown
negligible in Case A (0.0068 rad).

**`reproduce.sh` never used the venv it created.** On Windows a venv only ever
creates `python.exe`, never `python3`, so a bare `python3` after activation
silently ran the system interpreter: every `pip install` and every script ran
outside the isolation the script claimed to provide. Fixed by locating the
venv's own interpreter explicitly.

**Sample rate.** Case A used `fs = 4*f_isco = 502.5 Hz`, chosen when only an
inspiral-only waveform existed; real IMRPhenomD ringdown carries power out to
~640 Hz. Raised to 2048 Hz.

**The echo's arrival time.** Adding `Delta_T` to the *unlensed* merger time
misses the second image's echo by ~0.6 s, because the fixed geometric delay
runs from the strong image's own arrival while the combined near-merger peak
is itself displaced by interference between the two images — by almost exactly
the same 0.6 s. `case_A_chirp/run.py` now searches for the actual local
maximum and logs both the observed time and its offset from the naive
prediction, rather than assuming a reference point.

**The animation marker.** `y` diverges at each `D_LS = 0` crossing, sending
the marker's CSS position to values like `191%`, silently clipped by
`overflow:hidden`. Clamped to the visible stage edge.

## Errors found in `theory.tex` (2026-09-08..11)

The text was rewritten from a report into a textbook chapter over four days of
annotate-and-correct rounds on the PDF. Errors of substance found and fixed in
the process, all verified numerically before and after:

- **The Mardling & Aarseth stability margin was quoted with `q_out`
  inverted.** Their Eq. 90 defines `q_out = m3/(m1+m2)`; the reciprocal turns
  a threshold of 51 into 2.8 and a margin of 93x into 1694x. Resolved
  2026-09-12: the CODE was right all along (`M_LENS/MTOT`, and it reports
  93x); what carried the wrong number was `check_hierarchy`'s own docstring,
  which described `q_out` backwards and quoted 1700x. Docstring corrected,
  and it now names 1700x explicitly as the reading to watch for if this ever
  regresses.
- **The Kozai-Lidov comparison was stated backwards.** `t_KL = 16 yr` against
  240 d to merger is a factor 24 of *safety*; the text said it was "4% of the
  time to merger" and the only condition not negligible on its own. Also
  resolved 2026-09-12, the same way: the test's message was already reading
  it correctly, and only the prose around it was stale.
- `|F(w,0)|^2` was printed as `pi*w/2 / (1-e^{-pi*w/2})`; the identity
  `|Gamma(1+iz)|^2 = pi*z/sinh(pi*z)` gives `pi*w/(1-e^{-pi*w})`. The numbers
  quoted beside it were right; the formula was not.
- The coefficient relating `w` to frequency was printed as `1.24 f` instead of
  `6.19 f`.
- Several factors of `c` were missing: in the Shapiro term, in the geometric
  term of the time delay and its whole dimensionless chain, and in the
  curvature radius. The chapter declared `c=1` and then wrote `c` explicitly
  elsewhere. All restored; the text now carries `G` and `c` everywhere.
- The thin-lens approximation was justified against `D_L` and `D_S`, which are
  kiloparsec, when the comparison that binds is against `D_{LS}`, which is
  1.81 AU. Same cosmological reflex the chapter warns about.
- Kirchhoff's integral theorem was stated with the outward normal but with
  Born & Wolf's ordering, giving `-phi(P)`. Verified numerically on a plane
  wave: the printed form returned exactly `-1` where it must return `+1`.
- The radial form of `F` was attributed to Takahashi & Nakamura's Eq. 19,
  which is their singular-isothermal-sphere case.

Also added in the rewrite: an appendix deriving every result step by step,
including the closed-form point-lens `F` (Weber-Sonine, Kummer), which had
been cited and not derived; `y(theta)` in closed form; and the error budget of
the whole treatment, which is 1-2%.

## The behavioural audit of `src/` and `tests/` (2026-09-12)

The question was whether the code does what it says, not whether it says
something sensible, so the method was to probe behaviour and then to plant
bugs. 21 deliberate faults in `src/` -- the conjugation dropped from
`F_point_lens`, the literature sign on the SPA phase, `Delta_T` 1% long,
`beta` read at `t_obs` instead of `t_em`, the chirp-mass exponent, `D_LS` set
to `|z_los|`, and so on -- and then: which checks notice?

20 of 21 were caught. The survivor was dropping the 3 from
`orbital_redshift_factor`'s `(1-3GM/rc^2)^{-1/2}`, i.e. keeping the
gravitational redshift and losing the transverse-Doppler piece. Nothing
noticed, and the function's docstring claimed the split was "checked in
tests/test_doppler.py" when it was only printed there. Closed with a sixth
check that rebuilds the factor from `2GM/rc^2` and `(r Omega/c)^2` with
`Omega^2 = GM/r^3`, so the 3 is a result rather than a literal. No other check
could ever have caught it: a constant redshift produces no waveform feature.

Two validity conditions `theory.pdf` states as governing had nothing behind
them either, and both now do:

- **Point source, in the wave-optics sense.** Sec. 1.3 derives `a_in << b`,
  then says in as many words that this is "the weaker of the two" and that
  what matters in wave optics is the source size against the Einstein radius.
  It asserts the stronger condition holds "with room to spare" and gives no
  number. It is `a_in/eta_0 = 6.4e-3` at Case B, against `a_in/a_out = 2.1e-4`
  -- 30x tighter than what `check_hierarchy` was testing in its place.
- **The adiabatic condition**, which is Case B's whole licence. An external
  review reported it failing by a factor of 3. It does not: `|dln y/dt|^-1`
  goes to zero at both ends of the lensed window for any system whatsoever
  (`D_LS -> 0` there, so `y -> infinity`), and minimised over the whole lensed
  half the number is the sampling grid's, not the system's -- 43 s at 2e4
  samples per period, 0.26 s at 2e5, 0.03 s at 2e6. Evaluated where the
  prescription does something, `|F| > 1.01`, the margin is 5844x, and the
  check verifies it is the same at two grid densities.

A third was added when the same question was asked of the deflection itself:
no ray comes closer than the Einstein ring, `xi_0 = 121 r_g`, whatever `y`
does -- as `y -> 0` the images move onto the ring, they do not fall inward --
so the second-order Schwarzschild term is a fixed 2.4% there and 1.2% at the
`y` actually used. It scales as `M_L^(-1/3)`, so a heavier lens is worse.

`F_point_lens` was validated from first principles at exactly one point of the
`(w,y)` plane, and `w=1` is the one place where "error `= O(eta)`" and the true
"error `= O(eta w x_+^2)`" cannot be told apart. `F_radial_1d`'s docstring
claimed the default `eta=0.01` gives ~1.3%; at `w=25` the same call is 49%
off. The convergence check now also runs at `w=25`, just under `F_hybrid`'s
`w<30` cut, which is the demanding end of the range `F_point_lens` is actually
used over.

## The observing baseline, and what it is limited by (2026-09-12)

`T_OBS_B_S` went from six outer periods to three (24 d -> 12 d). The limit is
not the lens, it is the waveform: Case B uses the leading-order quadrupole
model, and the phase a 0PN model omits accumulates with the window while the
phase the lens writes does not. At `f_B` the source is at `v/c = 0.030`, so
the 1PN term is 5.8e-3 of the leading one -- but across tens of thousands of
carrier cycles that is still large next to `arg F = 0.018 rad`, even after
absorbing `phi_c`, `t_c` and `M_chirp`. Hence the framing: **Case B is an
amplitude result.** `|F|` is untouched by phase truncation.

The system itself was interrogated at the same time and left alone. Masses,
lens mass, outer period and inclination are each pinned by constraints that
pull against one another: lowering `f_B` kills the lensing (17.7% modulation
to 2.4%) unless `M_L` rises in proportion, which worsens the weak-deflection
term; growing the orbit improves nearly everything structural but lengthens
the observation linearly and narrows the usable inclination window from 3.0 to
1.0 degrees. The whole system lives inside 3 degrees of inclination.

## The figures (2026-09-12..13)

Redrawn so that each shows what its caption says. The ones that were not:

- `caseB_detector_view.png`'s zoom panel was captioned as showing the lensed
  and unlensed waveforms "visibly dephase over just a few carrier cycles".
  They do not. Both are evaluated at the same retarded time, so the 90-cycle
  Roemer delay is common to them and cancels, leaving only `F` -- whose phase
  at the pulse peak is 1.01 degrees, a 0.056 s shift on a 20 s carrier, under
  a tenth of a pixel. What the panel shows is amplitude: crests 17.7% taller.
- `caseA_F_of_f.png`'s top panel shaded the envelope across the whole band
  with no curve inside it, which rendered as a featureless rectangle: it
  asserted that the envelope is constant instead of showing it. Replaced by
  two resolved windows at `w=62` and `w=778`.
- `caseA_strain_time.png` drew the lensed waveform over the unlensed one,
  which simply hid the lower curve.

Three decisions worth keeping, because each was reached by trying the
alternative and looking at it:

- **Envelope and carrier cannot share a panel in Case B.** A pulse is 6653 s
  wide (FWHM) against a 20 s carrier: 333 cycles. Across ~1360 usable pixels,
  +-300 s resolves the cycles (45 px each) but `|F|` moves 0.4% of its peak
  height, so any envelope drawn on it is a horizontal rule; +-2500 s makes the
  envelope move but drops the cycles to 5 px and the carrier fills in as a
  block. Both were built. Shipped is the narrow one.
- **The orbit track came off the diffraction pattern.** It is a near-straight
  line crossing rings whose spacing is the entire content of the figure. A
  single marker at closest approach carries what it carried.
- **Titles are labels, not captions.** Several had grown into two-line
  explanations inside the image. The explanations moved to `cases/FIGURES.md`,
  which now has one entry per figure.

`F(f)` is drawn on its own dense grid rather than the waveform's FFT grid: the
fringe period is 0.29 Hz and the bin spacing 0.04 Hz, so the curves were being
drawn from ~7 points per fringe and came out polygonal. Same fix, same reason,
for the one-pulse zoom in `caseB_repeated_pulses.png`.

Colour bars inside a single decade were rendering "1.1 x 10^0", "6 x 10^-1" --
scientific notation whose exponent is identical on every tick. Both cases now
choose plain decimals with explicit ticks when the span is narrow.

## The deliverables pass (2026-09-13..17)

The figures, the deck and the report, gone over one at a time with the same
rule each time: render it, look at it, and check every number against
`numbers.json` instead of carrying it over.

**Case A's unlensed envelope was ringing**, and it was the band edge, not the
Hilbert transform and not the padding (both ruled out by measurement). A 0.5 Hz
half-cosine taper at `f_lower` fixed it. The width was chosen by measuring the
scatter of the envelope against its PN prediction at six widths: 4.0% with no
taper, 2.4% at 0.3 Hz, **2.0% at 0.5**, 5.0% at 1, 22% at 2, 34% at 3.5. It
also moved two end-to-end numbers, which is why they are re-measured and not
carried: the Paczynski agreement went from 3.3e-4 to **1.6e-4**, and the
TaylorF2-fallback one from 9.6e-9 to **7.8e-8**.

**`caseA_F_of_f.png` was one curve drawn three times.** Three windows of
`|F(f)|` and `arg F(f)` at different points of the band, which measurement
showed to be identical to four decimals -- the fringe period is 0.29120 Hz at
10 Hz, at 60 and at 123. Replaced by the thing those three windows were
evidence *for*: the circle `F` traces in the complex plane, one turn per
fringe, from which the modulus bounds, the phase bound and the period are all
read at once.

**The ratio panel's shrinking fringes were not a Hilbert artefact.** The panel
used to be truncated at -3.22 s, where `Delta_T (df/dt)/f` reaches 0.4, on the
argument that an analytic envelope cannot follow modulation faster than its own
carrier. That argument is wrong here. The lensed signal is two copies and the
analytic signal is linear, so the envelope ratio is exactly
`|sqrt(mu_+) + sqrt(|mu_-|) q(t) exp(i dphi)|` with `q = A(t-dT)/A(t)`: the
modulation DEPTH is `sqrt(|mu_-|) q(t)`, and `q -> 0` at both ends of the
window -- at the start the delayed copy has not arrived, towards the merger
`A(t)` runs away while `A(t-dT)` is still inspiral. Predicted against measured,
the depth agrees to **1.3%** from -8 s to 60 ms before the merger, which is
deep inside the region the old text called untrustworthy. The panel now runs to
the merger with those bounds drawn over it. Caught by a reader asking why a
ratio tending to 1 should be surprising.

**"Roemer/Doppler" as one label put the effect at the wrong orbital phase.**
Delay and shift are the same quantity a derivative apart, so where one is
extremal the other vanishes: at conjunction the delay is at **+904.9 s of
+-904.9** and `|v_los|/c` is **1.3e-6** against an orbital maximum of
**1.645e-2**. The lensing pulse coincides with the delay. The figure now draws
both on twin axes and the crossing pattern says it; `theory.tex` gained the
half of the quadrature argument it was missing (it had "separation goes as
sin, velocity as cos" but never said that the DELAY is therefore extremal at
conjunction, which is the half that is easy to lose in front of a plot).

**A phase bound for Case B that did not exist.** `|F|` had one for the
fixed-`w` approximation; `arg F` did not, and `arg F` is what
`caseB_doppler_vs_lensing.png`'s lower panel draws.
`arg_F_cycles_err_from_fixed_w_max` = **5.6e-4 cycles** turns that panel into
two things: the two large excursions at `y ~ 3` are 21x it and are a result,
the far ringing sits inside it and is not.

**Case B's detector view was showing 30 carrier cycles**, at which density the
two curves cross sixty times and read as one sinusoid drawn twice -- the figure
looked like it was *understating* the 17.7%. Narrowed to five. Nothing was lost:
`|F|` moves 0.064% of its peak across +-300 s and 0.002% across +-50 s. Also
recorded there, because nothing explained it: the top panel's baseline drifts up
**1.29%**, which is the source still inspiralling and matches
`(t_c/(t_c-T))^(1/4)` to five digits.

**One outer turn needed to be a page, not a picture.** 17280 carrier cycles per
turn, so no fixed figure shows the turn and the wave together. The first attempt
was a three-panel PNG that was correct and unusable -- the pulse it exists to
show occupies 1.9% of the axis. `caseB_one_period.html` generates the signal in
the browser instead of carrying it: `F(t)` and the retarded-time map tabulated
every 69 s, the carrier from closed forms. The check that those agree with the
repository's own phase, `one_period_closed_form_phase_err_cycles`, is what
caught the first version having the prefactor wrong by a factor of two --
`8pi/5` instead of `16pi/5`. It read **4407 cycles**; corrected, **5.8e-10**.
Locally the wrong version looked fine, which is the whole argument for checking
the generator against the repository rather than against the picture.

A reader then hit two things in that page immediately. Panning off the pulse
appeared to zoom, because each window was auto-scaled to its own maximum -- the
axis meant something different at every position, which destroys the one
comparison the page exists for. Fixed globally. And fixing your eye on a crest
while switching stages shows it move, which invites reading that as the Roemer
delay: it is not, since two windows a quarter turn apart are separated by
thousands of accumulated cycles of which Roemer contributes about 45. The
underlying intuition -- the wave is stretched, crests drag -- is right but
happens *within* a window, so the page gained a fixed-period comb and a drag
readout: **-0.001 cycles at the pulse, -0.089 at quadrature**. A comb measures
`d(tau)/dt`, not `tau`, which is why the drag vanishes exactly where the delay
is largest.

**The report was rewritten as a report.** It read as a lab notebook -- its own
changelog, citations to `wiki/log.md`, assertions that regenerated figures came
out byte-identical. All of that belongs here, not there. It is now organised
around the fact that the two epochs sit on opposite sides of the
diffraction/geometric-optics transition, and the validity numbers that used to
live only in `tests/test_system.py` are in it.

**Scope, twice.** The report's caption for `caseA_strain_time.png` was edited
without being asked -- `report/` is its own stage -- and reverted. And a
`theory.md` that had appeared as "the token-cheap version of the PDF" was
deleted: measured, it was 126915 characters against the `.tex`'s 133091, a 5%
saving in exchange for hand-maintained synchronisation that had already failed
once, leaving a file that said "generated from theory.tex" and was not. A
line-by-line dump is not a summary. `theory/OUTLINE.md` now says to read the
`.tex`.

**The repo's own entry points, last (2026-09-18).** `README.md` and
`reproduce.sh`, which are what a fresh agent hits first and which had drifted
behind everything else.

`reproduce.sh` never ran `theory/pulsos_esquema.py`. Once `theory.tex` started
including its output as figure 4.1, the document had a build input that no
reproduce path rebuilt -- the same failure `paraxial_validity.py` had until
2026-09-07, and the second time this shape of bug has appeared, which is worth
noting: a file a document depends on, generated by a script nothing calls. It
runs in `run_theory()` now, and `./reproduce.sh theory` was run end to end to
confirm it regenerates the figure and rebuilds the 41 pages.

`./reproduce.sh checks` on native Windows was then run deliberately, to see
what the README's own caution actually does. It does what it says:
`CHECKS_imr_waveform.json`'s pycbc check turned into `"passed": null,
"skipped": true`, because that venv cannot have pycbc. Documented, but silent
-- the script wrote over what the repository ships and said nothing. It now
prints a warning naming the files and the `git checkout` that undoes them.
Nothing refuses to run; the fallback is legitimate.

`README.md` had: `tests/` described as holding "every correctness check
referenced from `theory.tex`", which the theory text has not done since the
rewrite; `bibliography/`, `checks/` and `reproduce.sh` missing from the map; a
clone command naming a directory this repository is not called; and Case B
summarised with "Roemer/Doppler" as one thing, which is exactly the conflation
corrected everywhere else this week. Added: the echo as Case A's headline
result, `caseB_one_period.html` as a figure that is not an image, what
`provenance/` holds, and a line telling agents to read `theory/theory.tex`
rather than the PDF.

## Independent reviews

Three fresh-agent reviews were run during the build, plus several more over
`theory.tex` during the rewrite. They were working documents, never part of the
deliverable: `checks/independent_review/` is `.gitignore`d, and the three
review files were deleted on 2026-09-17 once everything actionable in them had
been acted on. The folder is kept for the next one.

What they were good for: the Fourier-convention bug, the `reproduce.sh` venv
bug, the aliased Case A animation panel, and most of the `theory.tex` errors
above. What they were not good for: the missing source kinematics, which all
three missed; and, in the later rounds, a steadily rising share of findings
that were matters of taste, or were reported against a stale copy of the file.
A review returns a list whether or not there is anything left to fix, so the
list needs to be filtered before it is acted on.

### A fourth pass, after the fact (2026-09-18)

One more review was run on the finished repository rather than during the
build. It arrived as a working tree of uncommitted edits, so it was triaged
against `origin/master` instead of accepted. It found four real errors, every
one of them in text that nothing executes:

- `theory.tex` Eq. (A.36) carried the `5/3` exponent twice. The line
  immediately below it already reads `m1 m2 Mb^(-1/3) = Mchirp^(5/3)`, so the
  equation contradicted its own next sentence. `chirp.py` was always right and
  no number moved.
- `report.tex` loaded neither `babel` nor `fontenc`. The Spanish report said
  "Figure 1..8" and hyphenated accented words as if they were English: eight
  captions wrong on a delivered PDF, for the sake of two preamble lines.
- Eq. (7.5) wrote the analytic signal as `h(t) e^(i phi)` rather than
  `A(t) e^(i phi)`.
- `update_report_data.py` carried a comment asserting that `json.loads`
  rejects `NaN` and `Infinity`. It does not; Python accepts both. The guard
  that comment described did not exist. Rejection now happens on the way out,
  with `allow_nan=False`, and `</script` is matched case-insensitively.

Plus four smaller ones: a stale "6 full outer periods" in `doppler.py` (it is
3), a `\section` title whose bare math broke its PDF bookmark, a citation
sending the paraxial-validity claim to chapter 5 when it is in chapter 4, and
`Maggiore2008` sitting in `refs.bib` uncited by anything.

The rest was rejected. A rewrite of `theory/OUTLINE.md` cut it in half and
dropped the measured argument for reading the `.tex` over the PDF; the same
paragraph was deleted from `README.md`; a scope caveat was inserted in eight
places at once; and a pass of hedging over `report.tex` and `report.html`
replaced physical statements with a reviewer's register. That caveat was
itself real and was kept, once in `theory.tex` 7.2 where `y(t)` is introduced,
and once in the report's limitations list.

Two things are worth remembering from it. Its own delivery check asserted that
every `refs.bib` entry is cited in `theory.tex`, ignoring `report.tex`, which
shares the same `.bib`; it then "found" two orphans that were not orphans and
added citations to satisfy itself. And the pattern from the earlier rounds held
again: everything it caught was in prose, equations and preambles, and nothing
it caught was in the code. Filtering mattered more than finding, four fixes out
of seventeen touched files.

### Six passes at once, one per part (2026-09-18)

Six reviewers were run in parallel over the finished repo, one per part
(theory, src+tests, cases, report.tex, the two HTML pages, and the entry
points), each told to report rather than edit, to say how it verified every
finding, and to mark its own confidence. Everything below was re-checked here
before being acted on; roughly a fifth of what came back was rejected.

The four that mattered:

- **report.tex Eq. (2) was the complex conjugate of what the code computes.**
  `waveoptics.F_geometric_optics` builds exactly that expression as `F_raw`
  and returns `np.conjugate(F_raw)`; theory's `eq:Fgeo` has `+i` with
  `e^{-i w DeltaT}`. `arg F` had the wrong sign at every frequency in the
  delivered PDF. This is the Fourier-convention bug's family, three years of
  project-time later and in a place no test looks.
- **theory 4.4 claimed the source's kinematics bias the lens mass** the way a
  cosmological `(1+z_L)` would. They do not: lens and observer share one
  static background, the frequency is conserved along the ray, and the
  kinematic shift happens upstream at emission, so it is already inside
  `f_obs` — which is what `w` is built from. A time-dependent `beta_los` could
  not have been a bias on a constant anyway.
- **theory 7.2 contradicted itself** about what separates the lens from the
  kinematics, opening and closing on "orbital phase" with a paragraph in
  between correctly establishing that the Roemer delay is extremal exactly
  where the pulse peaks.
- **Three checks could not fail.** `check_leading_order_group_delay` says it
  validates `spa_phase` and never calls it, leaving every TaylorF2 phase
  coefficient untested; the envelope check asserted only monotonicity, which
  any positive power satisfies; and the causality check's weak-image window
  was one-sided and 5x wide, so it passed with the exact normalization bug its
  own docstring says it rules out. Found by mutation, not by reading.

Below those, about twenty stale numbers and pointers, nearly all in the
provenance layer: a `claims.yaml` entry describing a figure that had been
replaced, a peak amplification of 1.065 that exists in no file (and whose
claim, "peak and RMS differ", is false — they agree to 3e-5), an echo bound
violated by its own number, a ring spacing quoted as `2*pi/(w*y)` when the
number beside it comes from `2*pi/(w*sqrt(y^2+4))`, and `report.html`'s
"con lensing" badge keyed to the alignment window, marking 9.5% of the orbit
where `fraction_of_time_lensed` is 0.4994 and the next slide says half.

**What was rejected, and why it is the interesting part.** Two reviewers
independently reported a units error on the same pair of numbers, and both
were wrong: `arg F` at the pulse is 0.0176 rad and the peak-to-peak lens phase
is 0.0174 cycles. Two different quantities that look alike, and each reviewer
grabbed the wrong twin. Both texts were already right, and so is the 5206
ratio, which compares peak-to-peak with peak-to-peak. A third reported three
orphan `refs.bib` entries; two were cited in `report.tex`, which shares the
same `.bib` — its check only looked at `theory.tex`, so it found orphans that
were not, and its own fix was to add citations satisfying itself.

The pattern from every earlier round held again, harder: everything found was
in prose, equations, docstrings, provenance and preambles. Nothing was in the
physics the tests cover — except where the tests turned out not to cover it.

What the round found and this pass did not act on is in `todo.md` rather than
here: four presentation defects in the deck, two in the one-turn page, and
four in `src/`/`cases/`. They are listed there with what was measured, so the
round's output survives having been filtered.

## Stated limits

Deliberately not done, and why:

- Higher than 2PN in `taylorf2.py`, absent a reference implementation to check
  the coefficients against.
- Eccentric outer orbit. The kinematics side is ready (`doppler.los_velocity`
  carries the eccentric branch and is exercised at `e=0.4`), but `z_orb` stops
  being constant and would have to be applied rather than reported.
- Feeding the observed frequency back into `F` for Case B: quantified as
  1.1e-4 in `|F|` at the pulse and 3.5e-3 at the worst point of the lensed
  half, and 2.5e-5 / 5.6e-4 cycles for `arg F`, at points where nothing is
  read off.
- A non-paraxial treatment of the ~60 s window around each `D_LS = 0`
  crossing, where the approximation is marginal and `F -> 1` anyway.
- The idealized detector view carries no noise curve, antenna pattern or
  matched filtering: it is an illustration, not a detectability claim.
