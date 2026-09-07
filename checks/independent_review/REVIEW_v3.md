# Independent review v3

Third independent review. Fresh session, no prior context, no communication with whoever
built this repo. Run 2026-09-07 on a Windows 10 checkout with WSL2 Ubuntu available, in a
git worktree at `HEAD = f2a0248`. `REVIEW.md` and `REVIEW_v2.md` were read as background
only; nothing either of them passed was taken on trust, and I re-derived the pieces that
mattered from scratch.

**I was able to run the `pycbc`/IMRPhenomD path**, which `REVIEW_v2.md` could not: the
repo's own `.venv-wsl/` (pycbc 2.11.0) is reachable from WSL, so every Case A number and
figure in this review was measured on the *committed* waveform, not the TaylorF2 fallback.

## Verdict

The reproducibility machinery is in excellent shape and the `REVIEW_v2.md` sign-convention
fix is correctly implemented — I re-derived the Fourier delay theorem independently, and
the conjugation in `waveoptics.py` is right. Running Case A from WSL reproduced
`cases/case_A_chirp/provenance/numbers.json` **digit-for-digit** (every field except the
`F_evaluation_seconds` timer) and every Case A PNG byte-for-byte; Case B reproduced
bit-identically too; `build_report.py` rebuilt `report.html` byte-identical; 26/26 checks
pass. The closed form `F_point_lens`, the Paczynski cross-check, the point-lens Jacobian,
`ΔT(y)`, the TaylorF2 2PN coefficients and its sign flip, `f_isco`, the Kepler geometry,
the paraxial-validity table and the `units.py` constants all check out — I verified each
numerically or by re-derivation and found no error in any of them.

What I found instead is a cluster of errors in the **exposition of the Case A time-domain
result**, all traceable to one thing the project never noticed: **`F` as coded omits
Takahashi & Nakamura's `exp(-i w φ_m(y))` reference-phase normalisation, so the whole
lensed waveform carries an arbitrary, convention-dependent time offset of
`4GM_L φ₊(y_A)/c³ = −0.6027 s`.** The repo observed that offset
(`peak_time_lensed_minus_unlensed_s = -0.6025`), asked "does that make sense?"
(`wiki/log.md`), and answered with a *physical* story — "an interference effect between two
images, by an amount unrelated to ΔT" — that is simply not what is happening. It is not
interference; it is a single number, `4GM_L φ₊/c³`, reproduced to better than half a
sample. I confirmed this by re-running Case A with the standard normalisation restored:
**the 0.60 s shift becomes exactly 0.000000 s, and the "naive prediction" the project says
misses by 0.6 s becomes exactly right.** That wrong explanation is currently repeated,
verbatim, in `RESULTS.md`, `claims.yaml`, `report.tex`, `report.html`, `wiki/log.md` and
`theory.tex` — six documents. The same block of prose also states the echo is "~2 orders of
magnitude" / "~100×" fainter than the merger; it is **4.1×** fainter, and its measured
amplitude ratio matches `√|μ₋|` to 0.7% — the single best quantitative confirmation in the
whole project, and it is misreported by a factor of 25 instead of being claimed.

Separately, `caseA_diffraction_pattern.png`'s left panel — the figure whose caption says the
central spike is "finite in wave optics where the geometric-optics magnification formally
diverges" — is computed **entirely from `F_geometric_optics`**, which diverges there. Its
central value (1123) is set by the innermost radial grid sample and is 5.8× the true finite
value (194.45). The figure demonstrates the opposite of its own caption, and cites
`check_einstein_ring_regularization` as evidence for it.

And the regression test added to catch the v2 sign bug,
`check_causality_of_lensed_pulse`, does not exercise `F_point_lens` at all: 1.6×10⁻³⁴ of
its test pulse's power sits below the `w_geo_threshold=30` switch, so only
`F_geometric_optics` is tested. The conjugation in the evaluator that *all* of Case B and
the sub-5 Hz part of Case A actually use is still untested.

On reproducibility I would pass this without hesitation — it is better than most published
code. On scientific rigour I would dock it: the headline Case A narrative explains a
normalisation artifact as physics, one figure illustrates the negation of its caption, and
the number that would have been the project's cleanest validation (`√|μ₋|`) is quoted
wrong by 25×. Several `REVIEW_v2.md` findings also remain unfixed (7, 8, 9 below).

## Reproduction log

Environment: Windows 10 Pro; repo `.venv/` (numpy 1.26.4, scipy 1.12.0, matplotlib 3.8.3,
mpmath 1.3.0) for native runs; repo `.venv-wsl/` under WSL2 Ubuntu (pycbc 2.11.0, numpy 2.x)
for the IMRPhenomD path. Worktree clean at `f2a0248` when I started.

1. Read `wiki/index.md` → `conventions.md` → `README.md`, then the source, in that order.
2. Ran all six `tests/test_*.py` natively: 26 `[PASS]`, 0 `[FAIL]`, 1 of the 26 a skip
   (`pycbc_merger_time_alignment`, still reported as a pass — v2 Finding 16, unfixed).
   Every `tests/CHECKS_*.json` came back byte-identical to the committed one except
   `CHECKS_imr_waveform.json`, which flipped to `SKIPPED` exactly as `README.md` warns.
3. Ran `cases/case_A_chirp/run.py` under WSL with pycbc. `numbers.json` reproduced
   **identically in every physical field**; all four PNGs byte-identical.
4. Ran `cases/case_B_monochromatic/run.py` both natively and under WSL. `numbers.json`
   identical in every physical field both times. The PNGs are byte-identical only from
   WSL — the native-Windows matplotlib produces visually identical but byte-different PNGs,
   so Case B's *committed figures* also came from WSL, which `README.md` does not say
   (it implies only Case A needs it).
5. `report/build_report.py` rebuilt `report.html` **byte-identical**. Good signal.
6. Re-derived independently, rather than trusting the tests: the Fourier delay theorem and
   the conjugation fix; `μ(x)=x⁴/(x⁴−1)` from the axisymmetric Jacobian; `ΔT(y)` against
   the closed form `y√(y²+4)/2 + ln[(√(y²+4)+y)/(√(y²+4)−y)]` (3.48594 at `y_A`, matching
   `time_delay_difference` to 12 digits); the TaylorF2 stationary-phase group delay against
   `chirp.time_to_merger` analytically; `|F(w,0)| = √(πw/2 / sinh(πw/2))·e^{πw/4} → √(πw)`
   against `check_einstein_ring_regularization`'s printed 1.812/3.963/7.927; the Case A
   band's `|F|` extrema against `√μ₊ ∓ √|μ₋| = 0.788704 / 1.267903` (exactly
   `min_abs_F_in_band` / `max_abs_F_in_band`); and `w = 8πGM_Lf/c³` against
   `4GM_Lω/c³`.
7. Rebuilt the whole Case A lensed waveform from a **pure two-image geometric-optics
   superposition** (two delayed, scaled copies of `H` at `4GM_Lφ₊/c³` and
   `4GM_L(φ₊+ΔT)/c³`) and compared against `run.py`'s `h_lensed`:
   `max|h_l − h_geo| / max|h_l| = 1.1×10⁻¹³`. That is what made Findings 1–3 findable.
8. Restored the worktree, then made the code simplifications in the section at the end and
   re-ran everything (native + WSL) to confirm identical output.

## Findings

### MAJOR

**1. [major] The "0.60 s early merger peak" is not interference — it is the un-subtracted
Fermat-potential reference phase, i.e. a convention artifact with no physical content.**

`src/gwlens/waveoptics.py:110-134` (`F_geometric_optics`) deliberately keeps the minimum
image's own, un-subtracted Fermat phase `exp(i w φ(x₊,y))`, and `F_point_lens`
(`waveoptics.py:216-237`) likewise omits Takahashi & Nakamura's `−i(w/2)·2φ_m(y)` term
(their Eq. 17). Both docstrings say so, on purpose, so that the evaluators agree with each
other. The consequence, which nobody traced: after conjugation, that factor is
`exp(−i·2πf·[4GM_L φ₊(y_A)/c³])`, i.e. **a rigid time shift of the entire lensed
waveform**.

```
φ₊(y_A) = ½(x₊−y_A)² − ln x₊ = −0.6118454438
4GM_L/c³                     =  0.9850981897 s
φ₊ · 4GM_L/c³                = -0.6027278391 s
numbers.json:
  "peak_time_lensed_minus_unlensed_s": -0.6025390625     (= -1234 samples at fs=2048)
  -0.6027278391 / (1/2048)  = -1234.39 samples  ->  nearest sample  -0.6025390625  ✔
```

It agrees to 1.9×10⁻⁴ s — *less than half a sample*. There is no residual for
"interference" to explain. `cases/case_A_chirp/RESULTS.md:96-99` says:

> "Near merger the GW frequency sweeps fast enough that F(f)'s phase varies rapidly across
> the band, so the *interference* between the two images can shift where their *sum*'s
> envelope peaks, by an amount unrelated to ΔT"

That is wrong on both counts: the shift is not from interference (the lensed waveform is
two cleanly separated copies, item 7 of the reproduction log), and the amount is not
"unrelated to ΔT" — it is `4GM_L φ₊/c³`, from the same Fermat potential ΔT comes from.
The same is true of `echo_peak_time_minus_prediction_s = -0.6031`: since the echo lands at
`t₊ + ΔT` and the "naive prediction" uses `t_unlensed + ΔT`, the difference is
*identically* `φ₊·4GM_L/c³`. `RESULTS.md:106-110` calls the agreement between the two
numbers "the same 0.6 s ... to three digits" and treats it as evidence for the
interference story; it is a tautology — both numbers *are* `4GM_L φ₊/c³`.

**Worse, the sign of the offset is arbitrary.** `ψ(x) = ln|x|` is a deflection potential,
defined only up to an additive constant (fixed here by the choice `ξ₀ = R_E`), so
`φ₊(y)` — and therefore the absolute arrival time of the lensed signal relative to
"no lens" — is *not a physical observable* in this formalism at all. Only `ΔT` between
images is. The correct answer to `wiki/log.md`'s "tiene sentido eso?" is not "yes, it's
interference"; it is "the question is not well posed as coded — the zero of the lensing
time delay is a convention, and this code happens to pick one in which the strong image
arrives 0.60 s 'early'." Takahashi & Nakamura subtract `φ_m(y)` precisely to remove this
ambiguity (it makes `F → √μ₊` with no phase, putting the strong image at `t = 0`).

I verified this directly by re-running the committed Case A pipeline with the normalisation
restored (`F → F·exp(+i w φ₊)`, everything else untouched, pycbc/IMRPhenomD):

| | as coded | TN2003-normalised |
|---|---|---|
| `peak_time_lensed_minus_unlensed_s` | −0.602539 s | **+0.000000 s** |
| echo time − unlensed peak | 2.8311 s | **3.4341 s = ΔT exactly** |
| echo time − lensed peak | 3.4336 s | 3.4341 s |
| `peak_strain_amplification` | 1.0907770 | **1.0646908** |
| `rms_strain_amplification` | 1.0554959 | 1.0554959 (Parseval, unchanged) |

So the entire "note on causality" — the thing the last commit was written to add — is
explaining away an artifact, and the headline number `peak_strain_amplification = 1.091` is
itself normalisation-dependent at the 2.5% level (a sub-sample shift landing differently on
the FFT grid). The `RESULTS.md:74` claim that it is "the single largest |h(t)| sample" is
still true; the claim that it is a *physical* 1.091 is not stable under a change of an
arbitrary constant.

Affected: `cases/case_A_chirp/RESULTS.md:86-117`;
`cases/case_A_chirp/provenance/claims.yaml:44-60` (`causality_and_echo`);
`report/report.tex:133-169`; `report/report_template.html:186-193` (and `report.html`);
`wiki/log.md:241-269`; `theory/theory.tex` is not wrong here but inherits the omission via
Eq. `eq:Fgeo`, which prints `e^{-iwφ(x₊,y)}` with no remark that this is a *choice* whose
value is arbitrary.

Recommended fix: either (a) adopt the TN2003 normalisation in `waveoptics.py` (multiply
every evaluator by `exp(+i w φ_m(y))` after conjugating, which keeps all four evaluators
mutually consistent and makes `check_geometric_optics_limit` still close), regenerate, and
delete the causality paragraph — the story becomes far cleaner: *the lensed merger peaks at
exactly the unlensed time, and a second copy arrives exactly ΔT later*; or (b) keep the
current convention but say plainly in one sentence that the absolute lensed-vs-unlensed
timing is a gauge choice and only `ΔT` is physical. Either way the current explanation must
go.

**2. [major] The second-image echo is ~4× fainter than the merger, not "~2 orders of
magnitude" / "~100×" — and its true amplitude is a `√|μ₋|` check the project never makes.**

Measured on the committed pycbc run:

```
unlensed peak envelope        = 3.99874e-20
echo_peak_env_lensed          = 9.65235e-21   (numbers.json, unchanged)
ratio                         = 0.24138       vs  sqrt|mu_-|      = 0.23960   (0.7% agreement)
echo / lensed peak envelope   = 0.23005       vs  sqrt|mu_-/mu_+| = 0.23300   (1.3% agreement)
```

That is **a factor of 4.1**, i.e. 0.6 dex. The claim appears as:

- `cases/case_A_chirp/RESULTS.md:101-102` — "a second, genuinely separate,
  ~2-orders-of-magnitude-fainter copy"
- `cases/case_A_chirp/provenance/claims.yaml:50-51` — "a second, ~2-orders-fainter ... copy"
- `cases/case_A_chirp/provenance/claims.yaml:143` — "the echo is ~2 orders of magnitude
  below the merger peak and invisible on a linear scale"
- `cases/case_A_chirp/run.py:333` (comment) — "the echo is ~1-2 orders of magnitude below
  the merger peak"
- `report/report.tex:150` — "~2 órdenes de magnitud más débil"
- `report/report_template.html:190` and `report.html` — "una segunda copia, **~100× más
  débil**"
- `wiki/log.md:256-257, 272-273`

The "~60x above it" statement immediately following in `RESULTS.md:116` *is* correct
(`9.65e-21 / 1.549e-22 = 62.3`) — but that compares the echo to the unlensed **FFT leakage
floor**, not to the merger peak. My best guess is that the two comparisons got conflated.
`caseA_detector_envelope.png` shows the truth plainly: on its log axis the echo sits about
a quarter of the way down from the merger peak, not two decades.

This matters beyond the wrong adjective: `0.2414` vs `√|μ₋| = 0.2396` is the sharpest
independent confirmation in the project that `F` is right in *amplitude and phase together*
— a prediction from pure geometric optics landing within 0.7% of a full FFT/IFFT
reconstruction of an NR-calibrated waveform. Neither the ratio nor `√|μ₋|` is in any
`numbers.json`. Log it and quote it; it is worth more than the sentence it would replace.

**3. [major] `caseA_diffraction_pattern.png`'s left panel demonstrates the opposite of its
own caption: its "Einstein-ring diffraction spike" is a divergent geometric-optics caustic,
cut off only by the grid.**

`cases/case_A_chirp/run.py:125` — `F_r = np.array([wo.F_hybrid(w, r) for r in r_grid])`
with `w = w_start = 61.9`. Because `61.9 > w_geo_threshold = 30`, **every** point of that
map comes from `F_geometric_optics`, and `magnification(x₊) → ∞` as `y → 0`. With
`r_grid = np.linspace(1e-3, ...)` (`run.py:124`) the innermost sample is `y = 10⁻³`:

```
plotted centre  |F_geo(61.9, 1e-3)|^2 = 1123.48        <- diverges as 1/y; halve r_grid[0], it doubles
true            |F(61.9, 0)|^2 = pi*w = 194.45         (= |F_point_lens(61.9,0)|^2, exactly)
overstatement                            5.78x
```

and the interpolator's `fill_value=(mag2[0], ...)` (`run.py:127`) paints that same 1123 over
the whole central disc. The colourbar in the committed PNG runs past 10³, confirming it.
Yet `cases/case_A_chirp/RESULTS.md:122-126` says:

> "a central bright spike at y=0 (the Einstein-ring diffraction spike — **finite in wave
> optics where the geometric-optics magnification formally diverges**, see
> `tests/test_waveoptics.py::check_einstein_ring_regularization`; ... the GW analogue of
> the Arago/Poisson spot)"

and `report/report.tex:97-101` repeats it in Spanish. The cited check is correct and passes
— but it tests `F_point_lens`, which this figure does not use. The figure shows the
*unregularised* caustic and is offered as proof of regularisation.

The rest of that panel is fine (the geometric-optics error is <2% for `y ≳ 0.2` and <0.1%
for `y ≳ 0.5` at this `w`; I checked), so the fix is small: evaluate the inner region with
`F_point_lens` regardless of `w`, e.g. `wo.F_point_lens(w, r) if r < 0.5 else
wo.F_hybrid(w, r)`, or simply state in the caption that the map is the geometric-optics
approximation and is not valid inside `y ≲ 0.2`.

**4. [major] Three documents claim `caseA_F_of_f.png` shows `arg F`. It does not — both
panels plot `|F(f)|`.**

`cases/case_A_chirp/run.py:238` and `:244` both plot `np.abs(F[...])`; there is no phase
panel anywhere in the figure (I opened the committed PNG to be sure — top: `|F(f)|` over
10–125.6 Hz; bottom: `|F(f)|` over 10–13 Hz). The claims:

- `cases/case_A_chirp/RESULTS.md:80-82` — "`caseA_F_of_f.png`, bottom sub-panel: arg F
  sweeps through a full 2*pi cycle roughly once per fringe"
- `cases/case_A_chirp/provenance/claims.yaml:42` — `evidence: caseA_F_of_f.png (bottom
  sub-panel: arg F)`
- `cases/case_A_chirp/provenance/claims.yaml:92` — `shows: |F(f)| and arg F(f) across the
  whole chirp band`
- `report/report.tex:126-127` — "(figura anterior, subpanel inferior: `arg F` barre un ciclo
  completo de 2π aproximadamente una vez por franja)" — and in `report.tex` the "previous
  figure" at that point is `caseA_strain_time.png`, not even `caseA_F_of_f.png`.

The underlying physical statement is true (`arg F` does advance ~2π per fringe), but it is
sourced to a panel that does not exist, in a repo whose stated rule is that nothing is
asserted without evidence. Either add the panel (`np.angle(F)` on a second axis — three
lines) or drop the citation.

**5. [major] `check_causality_of_lensed_pulse` — the regression test written to guard the
v2 sign bug — never exercises `F_point_lens`, and its search window is 0.02 s from
clipping the image it is looking for.**

`tests/test_waveoptics.py:127-179`. The test pulse is a σ=0.05 s Gaussian at 40 Hz;
with `M_L = 5×10⁴ M☉` the `w = 30` hybrid switch sits at `f = 4.847 Hz`. Measured fraction
of the pulse's `|H(f)|²` below that: **1.6×10⁻³⁴**. So `F_hybrid` returns
`F_geometric_optics` for every frequency that carries any power, and the conjugation in
`F_point_lens` — the evaluator **all** of Case B uses, and Case A uses below 4.85 Hz — has
no time-domain test at all. Delete the `conjugate()` from `F_point_lens` today and the full
suite still passes 26/26. Fix: add a low-`w` variant (e.g. `M_L = 5×10⁴`, pulse at
`f₀ = 1 Hz`, which gives `w ≈ 6.2`), or force `F_point_lens` explicitly.

Two smaller defects in the same test, both symptoms of Finding 1:

- The main pulse does **not** land at `t_ref = 6.0 s`. It lands at
  `t_ref + φ₊·4GM/c³ = 5.7141 s` (measured: 5.7141, amplitude 1.0820 = `√μ₊`). The test
  never checks the main image's position, and builds its search windows from `t_ref ± ΔT`.
  The weak image is therefore at 7.7698 s while the window is `[7.7495, 8.3495]` — inside
  by **0.0203 s**, i.e. 0.4σ of the pulse it is trying to find. Any change to `y`, `M_L` or
  the window width could silently move it out, and the test would report `amp_after ≈ 0`
  and fail for a reason that has nothing to do with causality.
- `expected = np.sqrt(mu_minus / mu_plus)` (`line 173`) is the wrong normalisation for what
  is measured. `local_max_near` returns the absolute amplitude of `lensed`, whose unlensed
  pulse peak is 1.0, so the prediction is `√|μ₋| = 0.4133`, not `√(|μ₋|/μ₊) = 0.3820`. The
  committed `CHECKS.json` prints "= 0.4101 (expect ~0.3820)" — 7.4% off — when the correct
  expectation would have matched to 0.8%. The comment on line 173 ("relative to the main
  image's own amplitude scale") describes a normalisation the code does not apply.

**6. [major] Case A contains no wave optics beyond the two-image geometric-optics limit,
but is repeatedly described as if it did.**

`cases/case_A_chirp/RESULTS.md:59-66`:

> "This is still genuinely wave optics (interference, not the smooth `w->infinity`
> envelope)"

and `report/report.tex:88-91`: *"así que esto es genuinamente interferencia de ondas"*. But
`F_hybrid` returns `F_geometric_optics` — literally the `w → ∞` asymptotic formula — for
every in-band frequency (`w ≥ 61.9`), and the IMRPhenomD waveform has no support below
10 Hz. I reconstructed `h_lensed` from nothing but two delayed, rescaled copies of the
unlensed waveform (`√μ₊` at `4GMφ₊/c³`, `−i√|μ₋|` at `4GM(φ₊+ΔT)/c³`, conjugated) and got

```
max |h_lensed - h_two_image| / max |h_lensed| = 1.1e-13
```

i.e. agreement to machine precision. The interference fringes and the echo are real and
worth showing — but they are *geometric-optics* two-image interference, and the exact
diffraction integral contributes nothing measurable. `RESULTS.md:52-56` already says every
in-band `F` value is the geometric-optics formula; the "genuinely wave optics" sentence
four lines later contradicts it. (Case B, at `w_B = 0.309`, *is* genuinely wave-optical and
uses `F_point_lens` throughout — that half of the story is sound.)

### MEDIUM

**7. [medium] `theory.tex` §7.2's adiabaticity argument is off by five orders of magnitude
and compares the wrong timescales.** `theory/theory.tex:916-919`:

> "`y(t)` varía en la escala temporal orbital externa (`P_out ~ días`), mientras que la GW
> misma oscila en una escala temporal `1/f_B ~ 20 s` — **nueve órdenes de magnitud más
> rápido**."

`P_out / (1/f_B) = 345600 / 20 = 1.73×10⁴` — **four** orders of magnitude, not nine. The
conclusion (adiabatic multiplication is fine) is correct and in fact more robustly so than
the argument given: the relevant timescale is not the GW period at all but the lens's own
response time, the image delay `ΔT = 3.43 s`, over which `Δy/y ~ 10⁻⁴`
(`ΔT/P_out = 1.0×10⁻⁵`). Fixing the number and the comparison would strengthen the section.

**8. [medium] `theory.tex` Ch. 5's retro-lensing argument still misuses `D_LS/(2R_Sch)`, the
half-fix of `REVIEW_v2.md` Finding 5.** `theory/theory.tex:735-741` now cites 919 instead of
the unsourced 1800, which is an improvement — but it labels it:

> "`D_{LS}/(2R_Sch) ~ 919` en el punto de lensing más profundo (**un parámetro de impacto
> típico de la órbita, en unidades del radio de Schwarzschild**)"

`D_LS` is the *line-of-sight* lens–source separation, not an impact parameter, and the
denominator is `2R_Sch`, not `R_Sch`. The actual transverse impact parameter at deepest
lensing is `ρ = 0.0951 AU = 96 R_Sch` (or `ξ₀ = R_E = 0.0599 AU = 61 R_Sch`) — an order of
magnitude below the number quoted, so "casi tres órdenes de magnitud por debajo [de unos
pocos R_Sch]" should read closer to one and a half. And, as `REVIEW_v2.md` already noted, a
*length* ratio is still not the *solid-angle* fraction the sentence claims to bound; the
solid-angle version is `~(b_c/D_LS)²/4 ≈ 5×10⁻⁷`, which supports the conclusion far more
strongly. The conclusion (retro-lensing negligible) is certainly right; the argument still
is not.

**9. [medium] Five `REVIEW_v2.md` findings are still unfixed in the current tree.** I
re-checked each against `HEAD`:

- **v2 #6** — `theory/theory.tex:751` ("Capítulo~5", inside Chapter 5, a self-reference) and
  `theory/theory.tex:899` ("con `h_no lensado(f)` del Capítulo~5"). Both should be
  Chapter 6 ("La binaria interna como fuente"). Still hard-coded, still wrong.
- **v2 #21** — `theory/theory.tex:248`: `t_d = (1+z_L)·(D_S/(D_L·D_LS))·[½|x−y|² − ψ(x)]`.
  Still missing `ξ₀²`; still dimensionally an inverse length; still does not produce the `w`
  defined two equations later.
- **v2 #22** — `theory/theory.tex:204`: `ḡ_μν = η_μν + 2U δ_μν`. Still the wrong sign
  (`η_μν − 2U δ_μν` gives the standard `g₀₀ = −(1+2U)`, `g_ij = (1−2U)δ_ij`, and is what
  the very next line's `n(r) = 1 − 2U` follows from).
- **v2 #23** — `report/report.tex:61-62` still cites "`theory.pdf` Cap.~1, §6.1" for the
  lens mass and the hierarchy; §6.1 is "El generador preferido: IMRPhenomD".
- **v2 #16** — `tests/test_imr_waveform.py:75` still `return True, "SKIPPED: ..."`, while
  the module docstring (line 5) still says "Skipped checks are reported as such, not
  silently passed". Running the suite natively today writes `"passed": true` for a check
  that did not run.
- **v2 #11 / #10** (cosmetic, unfixed) — `report/report.tex:188` still says the unlensed
  track is "verde **punteado**" where `run.py:182` draws `"--"` (dashed);
  `report_template.html:220` still tells the reader to look for a "verde/gris" marker whose
  CSS is cyan/grey.
- **v2 #15** (cosmetic, unfixed) — `caseA_diffraction_pattern.png`'s `suptitle` is still
  clipped at both ends in the committed PNG ("...int lens: ... y_A = 1.589, marke").

**10. [medium] `caseB_animation_data.json` is not valid JSON — it contains bare `NaN` and
`Infinity` — because `run.py` divides by `θ_E = 0` at each crossing.**
`cases/case_B_monochromatic/run.py:168-174` deliberately uses `|z_los|` (not the masked
version) to build a continuous track, so at the two `z_los = 0` samples `theta_E_traj = 0`
and `y1_traj = x_sky/D_L/0 = ±inf`, `y2_traj = 0/0 = nan`. Every run prints two
`RuntimeWarning: divide by zero` / `invalid value` lines (visible in the WSL run above).
The values then go straight into `json.dump`, which emits the non-standard `NaN`/`Infinity`
literals. It happens to work because `build_report.py` inlines the file as *JavaScript*
(where those are valid identifiers) rather than parsing it as JSON — but `JSON.parse` or any
strict parser rejects the committed file, and in the page `clampPct(NaN)` yields
`left: "NaN%"`, so the marker silently freezes for those frames. One-line fix: clip
`D_LS_traj_pc` away from zero, or mask those samples out of the export.

**11. [medium] `y_max` in Case B's provenance is stale, and half the exported orbit is still
off-frame.** `cases/case_B_monochromatic/provenance/claims.yaml:62` still states the choice
as *"y_max = 8 * min_y"*; `run.py:141` hard-codes `y_max = 25.0` (the v2 #3 fix — 8·min_y
would be 12.71). With the current value I measure, on the committed
`caseB_animation_data.json`: **458 of 900 orbit frames (51%) have `|y₁|` or `|y₂| > y_max`
and are clamped to the frame edge**, plus 1 non-finite frame (Finding 10); the finite `y₁`
range is `[−1028, +727]`. That is a genuine improvement over v2's 78%, and as v2 said no
finite `y_max` can fix it — but the provenance file should say what the code does, and the
panel header should say the marker parks at the edge for about half the orbit.

**12. [medium] `|z_los|/D_L ~ 10⁻⁸` is quoted in two places; the actual value is
1.8×10⁻⁹.** `src/gwlens/geometry.py:97` ("|z_los|/D_L ~ 1e-8 here") and
`theory/theory.tex:547` ("`|z_los|/D_L ~ 10^{-8}`"). With `a_out = 8.810×10⁻⁶ pc` and
`D_L = 5000 pc`, `a_out·sin i / D_L = 1.76×10⁻⁹`. Nothing depends on it (the approximation
is *more* justified than claimed), but it is the wrong order of magnitude in a repo that
otherwise quotes these carefully.

**13. [medium] `theory.tex` Eq. `eq:Fradial` (3.5) was not updated by the sign fix and, as
printed, equals `conj(F)`.** Eq. `eq:Fdimless` (2.5), `eq:Fgeo` (3.3) and `eq:Fclosed` (3.4)
were all rewritten in the corrected convention (`eq:Fdimless` even grows an explicit
`\overline{...}` bar). `theory/theory.tex:447-451` still prints

```
F(w,y) = -i w e^{iwy^2/2} \int_0^\infty dx  x J_0(wxy) exp[i w (x^2/2 - ln x)]
```

with no bar — which is exactly what `F_radial_1d` computes *before* its
`np.conjugate(F_raw)` on `waveoptics.py:213`. So the one displayed equation for the
Bessel-reduced integral is in the old (wrong) convention while the three around it are in
the new one.

**14. [medium] `F_radial_1d`'s docstring contradicts its own code comment (and
`theory.tex`) about the iε sign.** `src/gwlens/waveoptics.py:178-180`:

> "regularized with a small positive `eta` added to the exponent's imaginary part
> (x^2 -> x^2*(1-i eta)), i.e. **w -> w(1-i*eta)**"

but `waveoptics.py:195-197` says, emphatically, "the `(1+i*eta)` factor — **NOT**
`(1-i*eta)` — is what damps at large x; the sign was wrong in an earlier version and blew
the integral up", and `theory/theory.tex:454` also says `w → w(1+iη)`. The code is right
(`i·w(1+iη)x²/2 = i w x²/2 − ηwx²/2`, decaying); the docstring reinstates the exact bug the
comment three lines below warns about.

### MINOR

**15. [minor] `orbital_plane_to_sky` and `impact_parameter_of_time` disagree on which way
`z` points.** `src/gwlens/geometry.py:66-67` — *"z = line of sight (**toward the
observer**)"*; `geometry.py:108` — *"`z_los(t) < 0` (source instantaneously nearer Earth
than the lens, i.e. in FRONT of it) -> D_LS<0"*, and `geometry.py:126` sets
`lensed_mask = z_los > 0`. Those are opposite conventions. For the circular orbit actually
used it makes no numerical difference (the lensed half is symmetric, the choice just shifts
which half by `P/2`), which is why no test caught it — but there is no test of the `z` sign
at all, and with `e_out ≠ 0` (listed as a `todo`) it would put the pulses at the wrong
orbital phase.

**16. [minor] Case B's `y_min` "cross-check" is not independent.**
`cases/case_B_monochromatic/RESULTS.md:30-35` — *"found by two independently-written search
routines"* — and `claims.yaml:23-24` — *"found by independent search code in each case
script"*. Both are an `argmin` over the output of the *same* function,
`geometry.impact_parameter_of_time`, with the same parameters. The agreement to 4 s.f.
(1.588949 vs 1.589005) is a statement about grid resolution (2000 vs 3600 samples per the
respective scans), not an independent confirmation of anything.

**17. [minor] `report.html`'s inference from peak-vs-RMS is not valid.**
`report/report_template.html:184-185` — *"Amplificación de la muestra pico: 1.091.
Amplificación RMS: 1.056. La diferencia entre ambas es la señal de que F(f) desfasa, no sólo
amplifica."* A purely real, frequency-dependent `|F(f)|` would also make peak and RMS
amplification differ. The dephasing claim is true; this is not evidence for it. (The
`claims.yaml:31-39` version of the same claim is more carefully worded.)

**18. [minor] The "diffraction ringing" in Case B is real and correctly computed, but not
legible in the figure it is cited from.** I re-evaluated `|F(w_B, y(t))|²` on a grid 40×
finer than `run.py`'s 600-samples-per-period one across a full pulse: the two agree
point-for-point (peak 1.38453 both; the ringing minima at 0.9111 both), so the 600/period
sampling is adequate and the ringing is genuine diffraction, not aliasing —
`claims.yaml:35` and `report.tex:196-197` are correct. But at
`caseB_repeated_pulses.png`'s 24-day full width the entire ±0.6-day ringing structure is
compressed into ~25 px and renders as a solid smear; nothing about "anillado visible a
ambos lados" is actually visible. A single-period inset would show it. (The first minimum,
`|F|² = 0.911` at `y ≈ 4.2` — a 9% *dimming* flanking each pulse — is a nice result and is
in no `numbers.json`.)

**19. [minor] `wiki/conventions.md` is stale in two places, and it declares itself the tie-
breaker.** Line 15-16 puts the lens in the *"IMBH range, ~10^2-10^4 Msun"* while
`system.py:23` uses `M_LENS_MSUN = 5.0e4`, outside it (`wiki/log.md:56-59` acknowledges the
drift; `conventions.md` was never updated). Lines 122-124 still say the inner-binary
waveform is *"leading-order (Newtonian/quadrupole, restricted) df/dt ... explicitly not full
PN"*, which has not been true for Case A since the IMRPhenomD integration. Given the page's
own opening rule ("If code, `theory.tex` or `report` ever disagree with this page, this page
wins until corrected here first"), these matter more than ordinary doc drift.

**20. [minor] `wiki/todo.md` still says `theory.pdf` is 18 pages (it is 23) and points at a
file that does not exist.** `todo.md:4` — "(18 pages)"; I counted 23 `/Type /Page` objects
in the committed `theory.pdf`. `todo.md:35` cites *"`wiki/log.md` and
`wiki/running_on_wsl.md` for the documented workflow"* — `wiki/running_on_wsl.md` does not
exist anywhere in the repo (the WSL workflow lives in `README.md`).

**21. [minor] `claims.yaml` points at a function that does not exist.**
`cases/case_A_chirp/provenance/claims.yaml:8` — `evidence: run.py::build_waveform`. The
function is `build_waveform_fd` (`run.py:41`).

**22. [minor] `reproduce.sh all` leaves `report.html` out of sync, and the HTML build is
still gated behind LaTeX.** The header comment was fixed since v2, but the behaviour was
not: `run_cases` regenerates `cases/*/caseX_animation_data.json` while `report.html` (which
inlines them) is only rebuilt by `run_report`, which under `set -e` runs the four-step
`pdflatex`/`bibtex` block first — so `./reproduce.sh report` still cannot produce
`report.html` on a machine without LaTeX, and `./reproduce.sh` (`all`) still leaves the
committed HTML describing data it no longer contains. Moving the `build_report.py` call into
`run_cases`, or into its own `html` step, closes both.

**23. [minor] `tests/CHECKS_*.json` message strings are numpy-version dependent.** Re-running
`tests/test_taylorf2.py` under WSL (numpy 2.x) reproduces every *number* identically but
writes `"[np.float64(15.17110500274658), ...]"` where the committed file (numpy 1.26) has
`"[15.17110500274658, ...]"`. Same for `pn_terms_are_small_corrections`. Harmless, but it
means these two committed artifacts cannot be reproduced byte-identically across the repo's
own two supported environments. `float(...)`/`f"{x:.6f}"` in the message would fix it.

### NITPICK

- `README.md:137` says `checks/independent_review/` is "run once already"; there are two
  reviews there (three with this one).
- `src/gwlens/waveoptics.py:264` (`F_point_lens_array`) and `src/gwlens/geometry.py:24`
  (`y_of_angular_offset`) are defined and never called anywhere in the repo.
- `cases/case_A_chirp/run.py:117` sets `half_width = y_max` inside the `half_width is None`
  branch; the variable is never read after that point.
- `caseA_animation_data.json` still carries `h_unlensed` and `h_lensed` (69 KB of 190 KB)
  that `report_template.html` no longer reads — the v2 #2 fix switched the panel to the
  envelopes but left the raw arrays in the payload. The animation window also stops at
  `t_peak + 0.3 s` (`run.py:352`), so the echo — the most striking thing in the case — never
  appears in the interactive panel at all.
- `cases/case_A_chirp/RESULTS.md:36` gives the chirp duration as 15.18 s (`t_end`) but the
  4.4e-5 duty-cycle figure beside it comes from `check_static_lens_regime`, which uses
  `t_c = 15.199`. A 0.1% mismatch, stated as if one followed from the other.
- `RESULTS.md:56` says every plotted `F(f)` "is accurate to <4% by that test". The 4% comes
  from `check_hybrid_matches_at_threshold`'s worst case, which is `w = 30, y = 0.3` —
  neither of which Case A uses. At the actual `(w ≥ 61.9, y_A = 1.589)` I measure
  `|F_geo − F_exact|/|F| ≤ 3.9×10⁻⁴` across the whole band. The bound is honest but
  ~100× pessimistic; quoting the relevant number would read better.

## Areas where I found nothing wrong

Stated plainly, because most of the repo is in this category:

- **`src/gwlens/units.py`** — `G`, `c`, `GM_sun`, `pc`, `AU` all match CODATA 2018 /
  IAU 2015 / IAU 2012 to the digits given; `GM_sun/c² = 1.47662 km` and
  `GM_sun/c³ = 4.92549 µs` derived, never re-typed.
- **`src/gwlens/waveoptics.py`'s physics** — `image_positions`, `magnification`,
  `fermat_potential`, `time_delay_difference`, `total_magnification_paczynski`,
  `F_point_lens` (verified against `|F(w,0)| = e^{πw/4}√(π(w/2)/sinh(πw/2))` analytically),
  `F_radial_1d`'s Bessel reduction (re-derived), `F_bruteforce_2d`, and the conjugation fix
  itself. `w = 8πGM_Lf/c³` is consistent with `wφ = 2πf·t_d` given `t_d = (4GM_L/c³)φ`.
- **The Fourier-convention fix.** I re-derived the delay theorem from scratch, checked the
  code's `exp(-i2πfΔt)` delay in `imr_waveform.py:57`, checked the TaylorF2 sign flip
  against the stationary-phase condition analytically, and confirmed all four `F` evaluators
  are now mutually consistent *and* consistent with `numpy.fft`. `REVIEW_v2.md` Finding 1 is
  correctly and completely fixed.
- **`taylorf2.py`** — 2PN coefficients as cited; the overall sign; `spa_amplitude` identical
  to `chirp.restricted_pn_amplitude_fd`; group delay `(1/2π)dΨ/df = −(t_c − τ(f))`
  re-derived by hand and matching `chirp.time_to_merger` exactly.
- **`chirp.py`** — `f_isco`, `time_to_merger`, `freq_of_time` (algebraic inversion verified),
  chirp mass, the td-vs-fd amplitude distinction.
- **`geometry.py`'s orbit mechanics** — Kepler solver, true anomaly, the
  `R_z(Ω)R_x(i)R_z(ω)` projection, and the `D_LS(t) = z_los(t)` insight, which is the
  physically interesting part of this project and is right.
- **`theory/paraxial_validity.py` and Table 4.1** — all six numbers recomputed and matching:
  `R_Sch = 9.8706×10⁻⁴ AU`, `D_LS = 1.8147 AU`, ratio 919.2, `ξ₀/D_LS = 0.032983` (which is
  also exactly the deflection angle `2R_Sch/ξ₀`, a nice identity the doc could mention),
  59.8 s, `1.73×10⁻⁴`.
- **Parseval.** `rms_strain_amplification = 1.0554959020` computed from the time series
  matches `√(Σ|H·F|²/Σ|H|²) = 1.0554959020` to 10 digits — the `RESULTS.md:76-78` claim that
  it is phase-convention-independent is exactly right.
- **Every number in both `RESULTS.md` files traces to its `numbers.json`** at the stated
  precision. The only prose numbers that are *not* in a `numbers.json` are 0.79/1.27
  (`REVIEW_v2.md` Finding 8, still open but analytically exact), the "~2 orders of
  magnitude" of Finding 2 above (which is not merely unlogged but wrong), and the "~60x"
  ratio, which is a correct division of two logged numbers.

## Code simplifications applied

All of these are behaviour-preserving. After making them I re-ran the full suite plus both
case scripts, natively and under WSL with pycbc, and confirmed that every
`tests/CHECKS_*.json`, both `provenance/numbers.json` (every field except the
`F_evaluation_seconds` wall-clock timer), all eight case PNGs, both
`caseX_animation_data.json` and `report/report.html` come back byte-identical to `f2a0248`.

1. **`cases/case_B_monochromatic/run.py:82-84` → 82-83** — `chirp.phase_of_time(...)` was
   called twice with identical arguments, each call throwing away the half of the return
   value the other kept. Replaced with one call unpacking both. (`REVIEW_v2.md` Finding 24.)
   The function is deterministic, so the result is identical; it also halves the cost of the
   200 000-point integration.

2. **`src/gwlens/chirp.py:76`** — `phase_of_time(t, f_func, t0, t1, n=200_000)`'s first
   parameter `t` was never read in the body (the grid is built from `t0`, `t1`, `n`).
   Dropped it. The only caller in the repo is item 1 above, updated in the same commit.

3. **`cases/case_A_chirp/run.py:70`** — `n_pad = 1 << (2 * n - 1).bit_length()` →
   `n_pad = 2 * n`. `n` is set to a power of two on the previous line, so for any `n = 2^k`
   the old expression evaluates to `1 << (k+1) = 2n` identically (65536 both ways here).

4. **`tests/test_taylorf2.py:63`** — `n_pad = 1 << (n * pad_mult - 1).bit_length()` →
   `n_pad = n * pad_mult`, for the same reason (`n` a power of two, `pad_mult ∈ {2,4,8}`).
   `check_ifft_reconstructs_chirp_at_tc` reports the same three peak times to all digits.

5. **`cases/case_A_chirp/run.py:132`** — `np.hypot(Y1 - 0.0, Y2 - 0.0)` → `np.hypot(Y1, Y2)`.

6. **`cases/case_A_chirp/run.py:186`** — `t[i_peak_lensed] - t[np.argmax(np.abs(h_unlensed))]`
   → `t[i_peak_lensed] - NUMBERS["peak_time_unlensed_s"]`. The logged value *is*
   `t[np.argmax(np.abs(h_unlensed))]` (line 164), so this is the same float, without
   recomputing the argmax over 65 536 samples.

7. **`cases/case_A_chirp/run.py:336`** — the Figure-4 mask `view4` was rebuilt from the same
   expression as `view` (line 264), with the same `t_peak`, `t_end_seconds` and `t_hi`.
   Reuses `view`.

8. **`src/gwlens/imr_waveform.py:45-46`** — `H = np.zeros_like(freqs, dtype=complex)`
   followed by `H[:] = _taylorf2_banded(...)` → `H = _taylorf2_banded(...)`.
   `_taylorf2_banded` already returns a complex array of exactly that shape and dtype
   (`imr_waveform.py:96`), so the pre-allocation and copy were a no-op.

Deliberately **not** changed (noted as suggestions only, since each alters an API or a
plotted value): removing the unused `F_point_lens_array` / `y_of_angular_offset`; removing
the dead `half_width = y_max` in `make_ring_pattern`; dropping the unused `h_unlensed`/
`h_lensed` arrays from `caseA_animation_data.json`; and everything under Findings above.

## What I did NOT check

- I did not rebuild `theory.pdf` or `report.pdf` (no LaTeX run in this session). I verified
  by `git log` that each PDF was last committed in the same commit as its `.tex`, and read
  the `.tex` sources directly; every `.tex` line I quote is at `HEAD`.
- I did not open `report.html` in a browser. `REVIEW_v2.md` did that thoroughly, and I
  confirmed its two rendering findings were fixed *in the code* (the Case A panel now draws
  `env_unlensed`/`env_lensed`, `report_template.html:373-386`; the Case B marker now uses
  `θ_E(t)`, `case_B_monochromatic/run.py:168-174`) rather than on screen. My Finding 11
  numbers come from the exported JSON, not from pixels.
- I did not verify the IMRPhenomD phenomenological coefficients, or `pycbc`'s output, against
  any external reference — the project correctly declines to re-derive them, and I have no
  independent implementation here.
- I did not check the TaylorF2 2.5PN–3.5PN terms (not present) or re-verify the 2PN
  coefficients against the Buonanno et al. (2009) paper itself; I checked them against my own
  recollection of the standard form, as the two prior reviews did. No network access.
- I did not re-derive `theory.tex` Ch. 2's linearised-wave-equation → Helmholtz reduction
  from first principles. I checked the metric sign (Finding 9), the `ξ₀²` omission
  (Finding 9), the dimensionless reduction, and that `n(r) = 1 − 2U` and
  `(∇²+ω²)ψ = 4ω²Uψ` are the standard forms.
- I did not evaluate whether `M_L = 5×10⁴ M☉` at `a_out = 1.82 AU` with a 20+15 M☉ inner
  binary is a dynamically plausible or long-lived configuration; the project explicitly
  disclaims astrophysical typicality, and I checked only that the stated regime consequences
  follow, which they do.
- I did not measure figure claims off the PNGs pixel by pixel. Findings 3 and 4 were
  established by recomputing what the plotting code feeds `pcolormesh`/`plot`, then opening
  the committed PNGs to confirm the result is what is actually rendered.
- I have no way to verify authorship or timeline claims in `wiki/log.md`; I checked only
  that the current code matches each account it gives, which — with the exception of the
  causality narrative in Finding 1 — it does.

## Recommended order of fixes

1. **Finding 1** — decide the normalisation (subtract `φ_m(y)`, or state that absolute
   timing is a gauge choice), regenerate Case A from WSL, and rewrite the causality
   paragraph in all six places. If you take the subtraction, the story gets *better*, not
   worse: the lensed merger peaks at exactly the unlensed time and the echo lands at exactly
   `ΔT`.
2. **Finding 2** — replace "~2 orders of magnitude" / "~100×" with the measured 0.24 in all
   seven places, log `echo_env / unlensed_peak_env` in `numbers.json`, and claim the
   `√|μ₋|` agreement (0.7%) as the check it is.
3. **Findings 3 and 4** — use the exact evaluator inside `y ≲ 0.5` for the left panel of
   `caseA_diffraction_pattern.png` (or fix the caption), and either add the `arg F` panel or
   stop citing it.
4. **Finding 5** — extend `check_causality_of_lensed_pulse` to a low-`w` pulse so
   `F_point_lens`'s conjugation is actually guarded, assert the *main* image's position, and
   fix the `expected` normalisation.
5. **Findings 7–9, 12–14** — the `theory.tex`/`report.tex` corrections, each a line or two.
6. **Findings 10, 11, 19–23** — the export/provenance/staleness fixes.
