# Independent review v2

Second independent review. Fresh session, no prior context, no communication with whoever
built this repo. Run 2026-09-07 on a Windows 10 checkout (Git Bash, MiKTeX, Node.js, Chrome
and Edge installed). `checks/independent_review/REVIEW.md` (the first review) was read as
background only; nothing in it was taken on trust.

Working tree was clean at `abc8055` when I started, and is clean again now — I restored
every file `reproduce.sh` overwrote (`git checkout -- cases/ report/report.pdf
theory/theory.pdf`) and deleted the four temporary HTML harness files I created under
`report/`.

## Verdict

The reproducibility machinery still works — `./reproduce.sh`, `./reproduce.sh theory` and
`./reproduce.sh report` all ran end to end with zero intervention, 26/26 checks passed,
`report.html` rebuilt byte-identical, and every number quoted in `RESULTS.md`,
`report.tex`, `report.html` and `theory.tex`'s Table 4.1 traces correctly to its
`numbers.json` (the two rounding/prose nitpicks the first review found have been fixed).
But this review found something the first one did not: **the lensing amplification factor
`F` is written in Takahashi & Nakamura's Fourier convention (a delayed image carries
`+i·w·ΔT`) and is then multiplied into a waveform that is inverse-transformed with
`numpy.fft.irfft`, whose convention is the opposite** — so in Case A the faint saddle-point
image arrives **3.43 s *before* the merger instead of 3.43 s after it**, and the overall
lensing time shift has the wrong sign. I verified this with a decisive numerical experiment
(a Gaussian wave packet through the repo's own `F_geometric_optics` and `numpy.fft`
round-trip: the two images come out at `t₀+0.604 s` with amplitude `√μ₊=1.0283` and
`t₀−2.832 s` with amplitude `√|μ₋|=0.2396` — magnitudes exactly right, ordering acausal).
`wiki/conventions.md` §"Fourier convention" states the pairing that produces this: it
declares `h(f)=∫h(t)e^{−2πift}dt` and then asserts that *that* convention is the one in
which a later image gets a `+` phase. It is not; it is the one in which a later image gets
a `−` phase. The repo applied its own stated rule ("any formula copied from a paper using
the opposite convention must have its `i` sign flipped") to `taylorf2.py` and to
`imr_waveform.py`'s time shift, but not to the lensing `F` — which is the formula the whole
project is about. `|F|` is unaffected, so every Case B number, every diffraction-pattern
figure and every `|F(f)|` fringe plot is untouched, and `rms_strain_amplification` is
protected by Parseval; but `peak_strain_amplification=0.962` and the headline Case A story
built on it ("the lensed merger comes out broader and lower") rest on a phase that is
conjugated relative to the transform used. Beyond that, the new material is a mixed bag:
the pycbc integration is well engineered and honestly labelled, but the committed evidence
set is internally inconsistent (the check that validates the pycbc path is committed in its
*skipped* state, recorded as `"passed": true`); the interactive HTML works, is attractive,
and its lensed/unlensed colour coding is genuinely correct — but the Case A waveform panel
is decimated below Nyquist for 3.5% of its samples and renders the merger as an aliased
smear, and the Case B marker is clamped to the frame edge for 78% of the orbit; the new
Regge–Wheeler–Zerilli chapter's *formulae* are correct but its retro-lensing argument
forward-references an estimate that does not exist anywhere in the document. The Spanish is
fluent and technically competent but terminologically inconsistent (four different
hispanizations of "lensed" coexist). I would still pass this on reproducibility. On
scientific rigour I would now dock it materially: the sign-convention error is exactly the
class of bug this project's own log is proudest of catching twice, and it survived in the
one place it mattered most.

## What changed since the first review, and how it held up

### The pycbc / IMRPhenomD integration

**Design: good.** `src/gwlens/imr_waveform.py` is a clean try-pycbc-then-fall-back wrapper,
the fallback is a thin band+taper around the already-validated `taylorf2.htilde` (and
`check_fallback_matches_taylorf2_directly` proves it to 0.00e+00 relative error), and
`waveform_source` in `numbers.json` genuinely records which path ran. I confirmed that:
running `./reproduce.sh` here fell back to TaylorF2 and `numbers.json` correctly flipped to
`"TaylorF2 2PN inspiral-only (pycbc not available in this environment -- fallback...)"`.
So the labelling contract holds.

**But the committed evidence is inconsistent.** `cases/case_A_chirp/provenance/numbers.json`
is the WSL/pycbc artefact; `tests/CHECKS_imr_waveform.json` is the Windows artefact and
records `"pycbc_merger_time_alignment": {"passed": true, "message": "SKIPPED: pycbc not
importable in this environment"}`. So the one check that validates the alignment fix the
committed Case A results depend on has never been committed in a state showing it passed,
and `claims.yaml`'s `imr_waveform_and_bugs_fixed` cites it as evidence anyway. Worse,
`tests/test_imr_waveform.py`'s module docstring says *"Skipped checks are reported as such,
not silently passed"* while `check_pycbc_merger_time_alignment` literally
`return True, "SKIPPED: ..."` — printing `[PASS]`, writing `"passed": true`, and letting
the runner print `ALL CHECKS PASSED`. That is silently passing a skip.

**The Windows path silently overwrites the committed results.** `HEAD` is literally
*"Restore pycbc-generated Case A results after a reproduce.sh test run on native Windows
overwrote them"*, and my run reproduced exactly that: `peak_strain_amplification` went
0.9617 → 1.0178 and four PNGs plus `caseA_animation_data.json` were rewritten. `README.md`
warns about the waveform difference, but nothing warns that the default reproduce path
*destroys* the repository's own committed figures. A grader running the documented command
would end up with a repo whose `report.html` claims IMRPhenomD over TaylorF2 data (the
"IMRPhenomD real (inspiral+merger+ringdown)" sentence and the `MERGER / RINGDOWN` badge are
hard-coded in `report_template.html`, not driven by `waveform_source`).

### The Spanish translation

Fluent, idiomatic Rioplatense, technically competent, and it preserves the substance of the
English `wiki/log.md` accounts I compared it against (the TaylorF2 sign-bug checkbox in
`theory.tex` §6.2 is a faithful and complete rendering of the log entry, padding numbers and
all). I found no mistranslation that changes a physical meaning. What I did find is
terminological drift and a few weak word choices — see Findings 9–11. Note also that every
figure is still in English while the surrounding prose is Spanish; defensible (code stays
English) but visually jarring in `report.pdf`, where an English figure title sits directly
under a Spanish caption.

### The new Chapter 5 (Regge–Wheeler–Zerilli / Teukolsky)

The physics *formulae* are correct. I checked each claim:

- Eq. (5.1), the Regge–Wheeler potential `V_l = (1−2GM/c²r)[l(l+1)/r² − 6GM/(c²r³)]`, is
  right: the general spin-`s` RW potential is `f(r)[l(l+1)/r² + (1−s²)·2M/r³]`, and `s=2`
  gives `(1−4)·2M/r³ = −6M/r³`. ✔ Citation (Regge & Wheeler 1957) correct. Tortoise
  coordinate correct.
- Zerilli (1970) for polar parity, with the isospectrality claim ("distinto potencial pero
  el mismo espectro de transmisión/reflexión … dualidad axial-polar de Chandrasekhar") —
  correct; the standard attribution is Chandrasekhar & Detweiler (1975), and citing the
  1983 book is acceptable.
- Potential barrier peaking near `r=3GM/c²` (photon sphere) — correct.
- Teukolsky (1973), Weyl scalars `ψ₀,ψ₄` in the Newman–Penrose basis, master equation
  separable with spin-weighted spheroidal harmonics — correct.

The **retro-lensing argument that ties the chapter to the rest of the document does not
hold up** (Finding 5): it is a circular forward/backward reference to an estimate that is
never made, and the number it quotes (`~1800`) is a length ratio being offered as evidence
about a solid-angle fraction. Also, inserting this chapter broke two hard-coded chapter
cross-references (Finding 6).

## Did you manage to actually view report.html?

**Yes.** No browser extension was connected (`list_connected_browsers` returned `[]`), but
Chrome is installed at `C:\Program Files\Google\Chrome\Application\chrome.exe` and headless
screenshotting works:

```
chrome.exe --headless=new --disable-gpu --no-sandbox --hide-scrollbars \
  --virtual-time-budget=6000 --window-size=1400,14000 \
  --screenshot=full.png "file:///E:/CURSO_GW/Gravitational%20Lensing%20of%20GWs/report/report.html"
```

I rendered the full page (content is ~3960 px tall at 1400 px wide) and read it in three
tiles. To exercise the animations I wrote four temporary copies of `report.html` into
`report/` (so the `../cases/*.png` relative paths still resolved), each with an injected
`window.load` handler that set `#scrubA`/`#scrubB` and dispatched an `input` event, and
screenshotted each. I drove Case A to `frac = 0.5, 0.93, 0.985, 1.0` and Case B to
`frac = 0.08, 0.25, 0.42, 0.60`. All four temporary files were deleted afterwards
(`git status` clean).

**What actually works:** the page renders correctly in dark mode; both canvases draw; the
scrub handler fires and repaints both panels; the Case A fringe plot's red frequency cursor
moves; the `MERGER / RINGDOWN` badge appears at `t = 27.74 s`; the Case B readout counts
days and orbital periods correctly; and — the thing I most wanted to verify — **the Case B
marker's colour genuinely tracks `z_los > 0`**: cyan with the badge "detrás del lente" at
day 1.92, grey with "delante del lente — sin lensing" at day 14.4. That part of the physics
is correctly represented.

**What does not work well:** see Findings 2, 3 and 4. The merger/ringdown distinction is
*not* legible in the Case A panel, the Case B marker is parked on the frame edge most of
the time, and the panel text points the reader at a green marker that is actually cyan.

## Reproduction log

Environment: Windows 10 Pro, Git Bash, MiKTeX (`pdflatex`, `bibtex` on PATH), Node 22,
Chrome. Repo at `abc8055`, working tree clean, `.venv/` and `.venv-wsl/` both present
(gitignored).

1. Read `README.md` as a stranger would. Clear, honest, and the "read this before Case A
   surprises you" section is exactly the right thing to put there.
2. `bash reproduce.sh` — completed in ~90 s, no intervention.
   - All 5 `tests/test_*.py` printed `ALL CHECKS PASSED`: 26 `[PASS]`, 0 `[FAIL]`, 1 of the
     26 actually a skip (see above).
   - **Case A fell back to TaylorF2** (`pycbc` not importable in the Python that actually
     ran). `numbers.json`'s `waveform_source` correctly recorded the fallback. ✔
   - Diff vs. committed: Case A `peak_strain_amplification` 0.961714 → 1.017834,
     `peak_time_unlensed_s` 27.70955 → 27.71751, `abs_F_at_f_isco` differing in the 16th
     digit, `waveform_source` flipped; 4 PNGs and `caseA_animation_data.json` rewritten.
     Case B was **bit-identical** except `F_evaluation_seconds`. That is a strong
     determinism result for Case B.
3. `bash reproduce.sh theory` — clean, `theory.pdf`, **23 pages** (was 18 at the first
   review; the new chapter accounts for the growth). Only routine over/underfull-hbox
   warnings.
4. `bash reproduce.sh report` — clean, `report.pdf`, 7 pages (was 6), then
   `build_report.py` rewrote `report/report.html` **byte-identical** to the committed one.
   Excellent reproducibility signal for the HTML build.
5. Cross-checked every number in both `RESULTS.md` against `provenance/numbers.json`
   programmatically: **all match** to the stated precision, including the two the first
   review flagged (Case B now says `0.0520 Hz`, correct; `F_hybrid`'s docstring now says
   `w_geo_threshold=30` with `3.1% at w=30`, matching the code and the test output).
   One exception — see Finding 8.
6. Cross-checked `report.html`, `report.tex` and `theory.tex` Table 4.1 against
   `numbers.json` / `paraxial_validity_numbers.json`: all six paraxial numbers match
   (`9.87e-4`, `1.81`, `919`, `0.033`, `60 s`, `1.7e-4`), and every figure quoted in the
   Spanish report matches its source.
7. Restored the working tree.

**One reproducibility defect I hit along the way (Finding 12): `reproduce.sh` never uses
the venv it creates on Windows.** `source .venv/Scripts/activate` succeeds, but a Windows
venv creates only `python.exe`/`pythonw.exe` — no `python3.exe` — so every `python3` in
`reproduce.sh` falls through to whatever `python3` is on PATH. Here that was the Microsoft
Store Python 3.11, and `python3 -m pip install --quiet -r requirements.txt` installed the
pinned dependencies into **that interpreter's user site-packages**, not the venv.
Confirmed: `.venv/Lib/site-packages/` contains only `pip`, `setuptools` and
`pkg_resources` — nothing was ever installed into it — while
`python3 -c "import sys; print(sys.prefix)"` inside the activated venv prints
`C:\Program Files\WindowsApps\PythonSoftwareFoundation.Python.3.11...`. So `README.md`'s
claim that `reproduce.sh` "creates its own venv and installs the four pinned dependencies
(`requirements.txt`) into it" is false on Windows, the environment isolation the README
promises does not happen, and the script quietly mutates the user's global Python. It still
*works* (pip installs the pins wherever `python3` points), which is why nobody noticed.

## Findings

### CRITICAL

**1. [critical] Lensing `F(f)` is applied in the opposite Fourier convention from the
inverse transform, so the second lensing image arrives *before* the first.**
`src/gwlens/waveoptics.py:102` — `envelope = np.sqrt(mu_plus) - 1j*np.sqrt(np.abs(mu_minus))
* np.exp(1j * w * dT)` — and `waveoptics.py:103`'s `np.exp(1j * w * phi_plus)`. This is
Takahashi & Nakamura's form, which is written in the convention
`h̃(f) = ∫ h(t) e^{+2πift} dt` (the one in which a delay `t_d` contributes `e^{+2πif t_d}`;
it is also the convention forced by the diffraction integral `F = (w/2πi)∫d²x e^{+iwφ}` in
`F_bruteforce_2d`, `waveoptics.py:130`, and `theory.tex` Eq. (2.6), where `wφ = 2πf·t_d`
with `t_d ≥ 0`). But `cases/case_A_chirp/run.py:154-155` does
`np.fft.irfft(H_unlensed * F, n=n_pad)`, and numpy's forward transform is
`e^{−2πift}` — the opposite. Decisive test (Gaussian packet at `t₀=20 s`, `M_L=5e4`,
`y=1.589`, through the repo's own `F_geometric_optics` and `numpy.fft`):

```
primary image   t = 20.604 s  rel. amp 1.0283  (= sqrt(mu_+),  correct)
saddle image    t = 17.168 s  rel. amp 0.2396  (= sqrt|mu_-|,  correct)
separation      3.436 s       (= 4GM/c^3 * DeltaT = 3.434 s, correct magnitude)
```

Amplitudes and separation are exactly right; the **ordering is inverted**. The physically
required behaviour is the weaker saddle image arriving `+3.43 s` *after* the minimum image
(`ΔT = φ(x₋) − φ(x₊) > 0` — `theory.tex` §3.1 says so explicitly: *"la imagen de punto
silla `x₋` llega más tarde"*). The overall `e^{iwφ₊}` offset is likewise flipped, shifting
the whole lensed waveform `+0.603 s` later instead of `−0.603 s` earlier.

`wiki/conventions.md:51-59` is the root cause and states the error outright: it declares
`h(f) = ∫ h(t) e^{-2πift} dt` and then says *"This sign convention is the one that makes
the … geometric-optics limit of `F(f)` come out as a time delay `+i·2πf·Δt` in the exponent
(a later-arriving image gets a `+` phase), matching Takahashi & Nakamura (2003)."* Both
halves are true individually and false together: T&N use the **other** convention. The repo
already knows this in two other places — `taylorf2.py:66` deliberately returns
`-psi_literature_sign` with a long docstring explaining the flip, and
`imr_waveform.py:57` delays with `np.exp(-1j*2*np.pi*freqs*delta)` — so the codebase
contains a direct internal contradiction: `+i·w·ΔT` for a delay in one file, `−i·2πf·Δt`
for a delay in another, both fed to the same `irfft`.

**Impact.** `|F|` is invariant under conjugation, so *nothing* in Case B changes (Case B
uses only `|F|` for its envelope, `|F|²` and every reported number), and no
diffraction-pattern figure, `|F(f)|` fringe plot, or `rms_strain_amplification` changes
(Parseval — I verified `rms = 1.055849` identically with `F` and `conj(F)`). What does
change is every phase-sensitive Case A time-domain result. Measured on the TaylorF2 path
(the only one I can run here):

| | as coded | convention-consistent (`conj(F)`) |
|---|---|---|
| `peak_strain_amplification` | 1.0178 | **1.0599** |
| `t` of lensed peak | 28.320 s | 27.115 s |
| peak ratio inside the plotted/animated window | 0.493 | 1.060 |
| fraction of lensed power inside that window | 84.8% | 97.4% |

So the headline Case A claim — `RESULTS.md:74` *"Peak-sample amplification: 0.962 (the
single largest |h(t)| sample is smaller lensed than unlensed)"*, repeated in
`report.tex:116` and `report.html` — is computed from a conjugated filter, and in the
TaylorF2 analogue the fix moves it from below to further above 1. The qualitative
"smeared merger" story (`RESULTS.md:83-100`, `report.tex:119-132`, `wiki/log.md:156-162`)
is also partly a windowing artefact of the same bug: as coded, 12.7% of the lensed power
lands in `t ∈ [28.1, 30) s`, entirely outside the `[t_peak − 1.15·t_end, t_peak + 0.3]`
window that `caseA_strain_time.png`, `caseA_detector_envelope.png` and the HTML animation
all use — so the lensed waveform *looks* lower partly because a chunk of it is off-screen.
Fixing the sign puts 97.4% of the lensed power back inside the window.

**Secondary consequence worth noting:** with the correct sign, the saddle image lands
`3.43 s` *after* the merger — a genuine, physically real, delayed "echo" burst at 24% of
the primary amplitude, well inside the 32.6 s FFT buffer. That is arguably the single most
striking observable feature of this system, and the project currently plots it 3.43 s
before the merger, buried inside the inspiral, and never mentions it.

**Why no test caught it:** every check in `tests/test_waveoptics.py` compares evaluators of
`F` *against each other* (`F_point_lens` ↔ `F_radial_1d` ↔ `F_bruteforce_2d` ↔
`F_geometric_optics`). They all share the same convention, so they agree perfectly and the
error is invisible. There is no test anywhere that checks the *time-domain* lensing
response — e.g. "a wave packet through `F` must produce a *later*, weaker second image."
That is exactly the test `wiki/log.md:108-116` says it learned to write for TaylorF2
(`check_ifft_reconstructs_chirp_at_tc`, and its explicit note that the algebraic
self-consistency check "could not have caught this on its own"). The same lesson was not
transferred to the lensing filter.

### MAJOR

**2. [major] The Case A animation panel is decimated below Nyquist and renders the merger
as aliased noise — while the caption asks the reader to read the physics off it.**
`cases/case_A_chirp/run.py:279` — `stride = max(1, len(idx) // 1400)`. The window holds
~8920 samples at `fs = 502.5 Hz`, so `stride = 6` and the exported `caseA_animation_data
.json` has an effective sample rate of **83.76 Hz → Nyquist 41.9 Hz**, against a signal that
sweeps to `f_isco = 125.6 Hz`. I measured: **3.5% of the exported samples are above
Nyquist**, starting at `t = 27.40 s` — i.e. the entire merger and ringdown, the whole point
of the panel. Consequences:

- The decimated arrays give `max|h_lensed| / max|h_unlensed| = 0.602`, while the same
  quantity in `numbers.json` (computed at full resolution) is `0.962`. The picture and the
  number printed 40 px below it disagree by 36%.
- `report_template.html:186` tells the reader *"visible directamente arriba: el merger
  lensado sale más ancho y más bajo que el no lensado"* — pointing at an aliased trace.
  In the rendered screenshots the merger is an indistinguishable wall of cyan; no
  merger/ringdown structure is legible at all.
- `run.py:302-303` already exports `env_unlensed` / `env_lensed` (full-resolution Hilbert
  envelopes, alias-robust, 68 KB of the payload) — and `report_template.html`'s JS
  **never reads them**. Drawing the envelopes, or a rolling zoom window, would have fixed
  this for free.

**3. [major] The Case B marker is clamped to the edge of the pattern for 78% of the orbit,
so the animation's core visual claim rarely happens on screen.**
`report_template.html:433-435` clamps the marker to `[1.5%, 98.5%]`, with a comment
attributing this to the `D_LS → 0` divergence (`theory.tex` §4.2). That is not the dominant
cause. The exported orbit spans `y1 ∈ [−36.10, +36.10]` while the pattern PNG covers only
`±y_max = ±12.71` (`caseB_animation_data.json`), so **77.7% of the orbit samples fall
outside the plotted pattern entirely** and the marker sits pinned to the frame edge. I
confirmed this visually at four scrub positions: at day 1.92 and day 14.4 the marker is on
the left edge in both. The panel is simply zoomed ~2.8× too far in for the orbit it is
supposed to display; `y_max = 8.0 * min_y` (`cases/case_B_monochromatic/run.py:133`) was
chosen for the *static* figure, where the track is drawn as a line and clipping is
invisible, and reused unchanged for the animation, where it is not.

**4. [major] The animated marker's radius is not `y(t)`.**
`cases/case_B_monochromatic/run.py:147-149` builds the plotted track with a single
`theta_E_typ = np.median(theta_E[lensed_mask])`, but the physics uses the time-dependent
`θ_E(t) ∝ √z_los(t)` (`geometry.py:129-136`). So the marker's distance from the lens is
`ρ(t)/(D_L·θ_E_typ)`, not `y(t) = ρ(t)/(D_L·θ_E(t))`. I measured the discrepancy: at the
moment of deepest lensing the true `y = 1.589` but the marker is plotted at radius
**1.889 — 19% too far out**, where `|F|² = 1.317` rather than the 1.385 the readout
simultaneously displays. The axes are labelled `y₁, y₂`, which by `wiki/conventions.md`
means the dimensionless impact parameter, so the label is wrong for the track. Some fixed-
`θ_E` projection is unavoidable (a fixed `|F|²(y)` map cannot host a track whose `θ_E`
varies), but it is nowhere stated — neither in `caseB_pattern_with_orbit.png`'s title nor
in `RESULTS.md:50-55` nor in the HTML panel header — and it silently desynchronises the
animation's two panels at exactly the moment of interest.

**5. [major] The retro-lensing argument in the new Chapter 5 forward-references an estimate
that does not exist.** `theory.tex:660-665` (Ch. 5) says: *"Esta es la versión rigurosa y
cuantitativa del argumento que descarta la retrodispersión ('retro-lensing') … (§4.3): ahí
se estimó que la fracción de ángulo sólido de la radiación de la fuente que pasa lo bastante
cerca del lente como para retrodispersarse es despreciable, comparando `a_out` con
`R_Sch(M_L)` (razón ~1800)."* §4.3 (`theory.tex:515-592`) contains no such estimate — it
only *defers to Chapter 5*, at line 588-592: *"El Capítulo 5 discute … y por qué ese mismo
formalismo también resolvería la pregunta relacionada de si parte de la radiación … puede
retrodispersarse."* The two sections point at each other and neither makes the argument.
I grepped: `1800` appears exactly once in `theory.tex`, in the Chapter 5 sentence, and in
no `numbers.json` — so this is the one number in the whole document with no provenance,
in a repo whose stated rule is that no stated number lacks a source. (It is arithmetically
right: `a_out / R_Sch = 1841`. But it is a *length* ratio being offered as evidence about a
*solid-angle* fraction; the actual solid-angle fraction subtended by the strong-deflection
cross-section `b_c ≈ 2.6 R_Sch` at `a_out` is `~(b_c/a_out)²/4 ≈ 5×10⁻⁷`, which supports the
conclusion far more strongly than `1/1800` does.) The conclusion — retro-lensing is
negligible — is certainly correct; the argument as written is not one.

**6. [major] `theory.tex`'s new chapter broke two hard-coded chapter cross-references.**
Inserting "Más allá de este proyecto" as Chapter 5 pushed "La binaria interna como fuente"
to Chapter 6, but two references to it are hard-coded, not `\ref`:
- `theory.tex:677` (inside Chapter 5 itself): *"distinto del ringdown del remanente de la
  propia binaria fuente (Capítulo~5)"* → renders as a **self-reference to Chapter 5**;
  should be Chapter 6. Verified in `theory.pdf` p.14.
- `theory.tex:825` (§7.1): *"con `h_no lensado(f)` del Capítulo~5"* → should be Chapter 6.
  Verified in `theory.pdf` p.15.

Every other chapter reference in the document uses `\ref` and is correct.

**7. [major] `theory.tex` §2.1's stated validity criterion is violated by Case B's own
design.** `theory.tex:153-158`: *"el tratamiento de rayos/lente delgado de abajo es válido
cuando `λ_GW ≪ R`. Para un lente puntual `R ~ R_E` (el radio de Einstein)."* I computed the
numbers for this system: at the deepest-lensing phase `R_E = 8.95×10⁹ m`, and at
`f_B = 0.05 Hz`, `λ_GW = 6.00×10⁹ m`, so **`λ_GW / R_E = 0.67` — not ≪ 1**. This is not an
accident; it is the *point*. `w_B = 0.31 = O(1)` is precisely the statement that the GW
wavelength is comparable to the Einstein scale, and `wiki/log.md` records the lens mass
being chosen to make it so. The physical conclusion survives — the curvature-coupling term
is controlled by `λ/R_curv`, not `λ/R_E`, and `R_curv ≈ R_E√(R_E/R_g) = 7.5×10¹⁰ m` gives
`(λ/R_curv)² ≈ 0.4%` — but the criterion as written identifies the curvature scale with the
diffraction scale, which is the one identification this project cannot make. A GR-literate
grader will notice that the document's stated justification for its central approximation
fails by its own numbers.

### MINOR

**8. [minor] `RESULTS.md` states two numbers that are not in any `numbers.json`, in
violation of its own opening sentence.** `cases/case_A_chirp/RESULTS.md:68` — *"|F(f)|
itself stays modest: it oscillates between about 0.79 and 1.27 throughout"*. Neither 0.79
nor 1.27 appears in `provenance/numbers.json` (I grepped). The file's own line 3-5 says
*"every number quoted here is in `provenance/numbers.json` … not hand-typed"*, and
`wiki/conventions.md:115-116` states the rule globally. The values are analytically correct
(`√μ₊ ∓ √|μ₋| = 1.0283 ∓ 0.2396 = 0.7887 / 1.2679`, and they appear in
`caseA_animation_data.json`'s fringe array as 0.78928 / 1.26748) — but they are not logged,
so the provenance chain is broken exactly where the project claims it never is.

**9. [minor] Four different Spanish renderings of "lensed" coexist, sometimes in the same
paragraph.** `lenteada` (`report.tex:27`, `theory.tex:39`), `lensado/lensada`
(`report.tex:110,116,122`; `report_template.html:158` "lensada vs. no lensada"),
`lenteado/no lenteada` (`report.tex:166`; `report_template.html:219`), plus the untranslated
English `lensing` throughout (`pulsos de lensing`, `geometría de lensing`, `sin lensing`).
`report.tex` alone uses three of the four. For a document whose whole subject is the noun
being conjugated, one choice should have been made and enforced.

**10. [minor] `report.html`'s prose describes the wrong marker colour.**
`report_template.html:220` — *"se ve arriba como el marcador saliendo del patrón
(verde/gris) durante la mitad 'de adelante' de cada vuelta"*. The marker's CSS is
`--lensed: #17b4d4` (cyan) / `--unlensed: #9a9a9a` (grey) (`report_template.html:11`,
`73-78`); I confirmed cyan-and-grey in the rendered screenshots. "Verde" is left over from
the *static* figure `caseB_pattern_with_orbit.png`, where the unlensed track is drawn in
`lime` (`cases/case_B_monochromatic/run.py:157`). The reader is told to look for a green
marker that does not exist.

**11. [minor] Two Spanish word choices weaken the key physics statement.** `report.tex:123`
renders "smeared out" as *"esparcido"* (scattered / spread about) — `difuminado`,
`emborronado` or `ensanchado` carry the temporal-broadening meaning; `esparcido` does not.
`report.tex:129` renders "scrambles the coherent buildup" as *"revuelve la acumulación
coherente"* — `revuelve` is "stirs"; `destruye`/`aleatoriza la fase de` would be precise.
Also `report.tex:151` says the unlensed track is *"verde punteado"* (dotted) where the code
draws `"--"` (dashed) and the English `RESULTS.md:53` correctly says "dashed green".

**12. [minor] `reproduce.sh` creates a venv it never uses on Windows.** See the
Reproduction log above. `reproduce.sh:20-28` — `python3` does not exist inside a Windows
venv's `Scripts/`, so every `python3` invocation escapes to the system interpreter and
`pip install -r requirements.txt` populates the user's *global* site-packages.
`README.md:31-34` promises isolation that does not occur. A one-line fix (`PY=python` after
activation, or `python -m venv` + `"$VIRTUAL_ENV/Scripts/python"`) would close it.

**13. [minor] `reproduce.sh`'s own header comment about `report.html` is stale, and `all`
leaves the HTML inconsistent with `cases/`.** `reproduce.sh:8-10` says *"report/report.html
needs no rebuild step, it just references the same PNGs cases/ writes"* — false since the
animation JSON started being inlined; `run_report` now calls `build_report.py` precisely
because a rebuild *is* needed. Two consequences: (a) `./reproduce.sh` (`all`) regenerates
`caseA_animation_data.json` but not `report.html`, so the two silently diverge; (b) because
`build_report.py` runs *after* the `pdflatex` block inside the same `run_report()` under
`set -e`, **you cannot rebuild `report.html` at all without a LaTeX installation**, even
though the HTML needs none.

**14. [minor] `wiki/todo.md` and `wiki/index.md` were not updated across the last three
commits.** `todo.md:4` says `theory.pdf` is "(18 pages)" — it is now 23. `todo.md:31-33`
still lists as *open*: *"`pycbc`/LALSimulation as a cross-check … if a working Windows
install path is ever found (the attempt this session hung indefinitely; **abandoned**)"* —
directly contradicting `README.md`, `RESULTS.md`, `theory.tex` §6.1 and `log.md`, all of
which say pycbc is integrated and generated the committed results. `todo.md:34-36` says an
animation was *"dropped in favor of the static pattern+orbit-track figure"* — it was
subsequently built. `index.md:24` lists the inner-binary waveform modules as
`chirp.py` and `taylorf2.py` and **omits `imr_waveform.py`**, now the primary one. The
first review flagged stale-comment drift as its main finding; the category got worse, not
better.

**15. [minor] `caseA_diffraction_pattern.png`'s title is clipped in the committed PNG.**
The second line of the `suptitle` renders as `int lens: axisymmetric rings; … y_A = 1.589,
marke` — the leading `(po` and trailing `d)` are cut off by the figure edge
(`cases/case_A_chirp/run.py:252-256`: `fig.suptitle(...)` followed by `fig.tight_layout()`,
which does not reserve width for an over-long suptitle). This is baked into the committed
file, so it appears in `report.pdf` (p.3), `report.html`, and anywhere else the figure is
shown. Purely cosmetic, but it is the most prominent figure in Case A.

**16. [minor] A skipped check is recorded as a pass.** `tests/test_imr_waveform.py:75` —
`return True, "SKIPPED: pycbc not importable in this environment"` — prints `[PASS]`,
writes `"passed": true` to `CHECKS_imr_waveform.json`, and lets `main()` print
`ALL CHECKS PASSED`, while the module docstring (line 5) claims *"Skipped checks are
reported as such, not silently passed."* A separate `"skipped"` state would cost three
lines and would make the committed `CHECKS_imr_waveform.json` honest.

### NITPICK

**17. [nitpick] The title reads "de Ondas de Ondas".** *"Lente Gravitacional en Óptica de
Ondas de Ondas Gravitacionales"* (`report_template.html:111`, `report.tex:15`,
`theory.tex:27`). Grammatically fine, but the stutter is the first thing a reader sees on
all three documents. *"Lentes gravitacionales de ondas gravitacionales en óptica de ondas"*
would avoid it.

**18. [nitpick] ~37% of `report.html`'s 348 KB is data the page never reads.**
`caseA_animation_data.json`'s `env_unlensed` + `env_lensed` (68.4 KB) and `t_end`, and
`caseB_animation_data.json`'s `envelope.unlensed` + `envelope.lensed` (41.2 KB) plus the
unused second half of the orbit arrays (~19 KB — `orbit` covers 1.5 periods, `frameB` uses
`tCur % P_days` so only the first is reachable). The Case A envelopes are the ones that
would have fixed Finding 2.

**19. [nitpick] The Case A readout exposes raw FFT-buffer time.** `t = 10.25 s … 28.01 s`
with the page simultaneously telling the reader the chirp lasts 15.2 s, and the readout's
first frame showing `f ≈ 9.5 Hz` when the stated band starts at 10 Hz
(`caseA_animation_data.json`'s `f_of_t[0] = 9.490`). Subtracting the window start would
make the axis mean something to a reader.

**20. [nitpick] The Case A fringe cursor is frozen for roughly half the animation.** The
fringe panel covers 10–13 Hz, which the chirp crosses at `t ≈ 7.3 s` of its 15.2 s; the
cursor is clamped to `fringeXs[last]` (`report_template.html:390`) for the rest, under a
header that says *"el marcador sigue la frecuencia instantánea"*. Confirmed at scrub 0.93,
0.985 and 1.0 — cursor pinned to the right edge in all three.

**21. [nitpick] Dimensionally inconsistent equation in `theory.tex` §2.2.**
`theory.tex:233` writes `t_d = (1+z_L)·(D_S/(D_L·D_LS))·[½|x−y|² − ψ(x)]`. The bracket is
dimensionless and the prefactor is `1/length`, so `t_d` comes out as an inverse length. The
standard form (Takahashi & Nakamura Eq. 8) carries `ξ₀²`: `t_d = (1+z_L)·(D_S ξ₀²/(D_L
D_LS))·[…]`. With `ξ₀ = R_E`, that prefactor is `4GM_L/c²`, which is exactly what makes
`w = 2πf·4GM_L/c³ = 8πGM_L f/c³` (Eq. 2.8) come out right — so the omission is a typo, not
a propagated error, but the displayed equation as printed does not produce the `w` two lines
below it.

**22. [nitpick] Weak-field metric sign in `theory.tex` Eq. (2.3).** `theory.tex:189` writes
`ḡ_μν = η_μν + 2U δ_μν` with `U = −GM_L/r`. That gives `g₀₀ = −1 + 2U = −(1 + 2GM/r)` and
`g_ij = (1 + 2U)δ_ij`, i.e. both signs opposite to the standard isotropic weak-field metric
(`g₀₀ = −(1+2U)`, `g_ij = (1−2U)δ_ij`, i.e. `η_μν − 2U δ_μν`). The very next lines are
correct — Eq. (2.5) `(∇²+ω²)ψ = 4ω²Uψ` and `n(r) = 1 − 2U` both match the literature and
follow from the *correct* metric — so this is an isolated sign slip in one displayed
equation with nothing downstream depending on it.

**23. [nitpick] `report.tex:61-62` cites `theory.pdf` "Cap.~1, §6.1"** to justify the lens
mass and the hierarchy. §6.1 is "El generador preferido: IMRPhenomD" and discusses neither.
Chapter 1's `choicebox` does justify the lens mass; the `§6.1` half of the citation appears
to be collateral damage from the same chapter renumbering as Finding 6.

**24. [nitpick] Duplicated work in `cases/case_B_monochromatic/run.py:82-83`.**
`chirp.phase_of_time(...)` is called twice with identical arguments, each call discarding
the other's return value. Harmless, just 2× the cost and confusing to read.

**25. [nitpick] `wiki/log.md`'s bug accounts all check out.** I spot-checked five against
the current code and every one is accurate: the restored `exp(i·w·φ₊)` factor
(`waveoptics.py:103`), the TaylorF2 sign flip (`taylorf2.py:66`), the pycbc alignment shift
`H·exp(−i2πfΔ)` (`imr_waveform.py:57`), `make_ring_pattern`'s `center`/`half_width` display-
grid fix (`case_A_chirp/run.py:91-113`), and Case B's dedicated fine grid for the raw-
waveform panel (`case_B_monochromatic/run.py:186-198`). Noted here as a positive: the log
is a genuinely trustworthy record, which is rarer than it should be. The irony of
Finding 1 is that the log describes catching precisely this class of error, twice, in
adjacent files.

## What I did NOT have time/ability to check

- **I could not run the pycbc/IMRPhenomD path.** `pycbc` is not importable from any Python
  on this machine's PATH and I did not enter the `.venv-wsl/` WSL environment. Everything I
  measured about Case A's time-domain behaviour (Finding 1's table, Finding 2's aliasing
  numbers) was measured on the TaylorF2 fallback. The convention error in Finding 1 is
  waveform-independent — it is a property of `F` and `numpy.fft`, and I proved it with a
  synthetic wave packet, not with either waveform — but **the exact corrected value of
  `peak_strain_amplification` for the committed IMRPhenomD run has to be recomputed in WSL.**
  I can say the number will change and that the TaylorF2 analogue moves 1.018 → 1.060; I
  cannot say what 0.962 becomes.
- I did not re-derive `theory.tex` Chapter 2 (linearised wave equation → Helmholtz) from
  first principles; I checked Eq. (2.5), Eq. (2.6), Eq. (2.8) and the `ξ₀²`/metric-sign
  issues against the standard forms and the internal consistency of the chain, nothing more.
  The `−2R_μανβh^αβ` sign in Eq. (2.2) is convention-dependent and I did not pin down which
  Riemann convention is intended; since the term is dropped, it does not matter here.
- I did not verify the IMRPhenomD phenomenological coefficients or `pycbc`'s output against
  anything — the project explicitly and correctly declines to re-derive them, and I have no
  independent reference here.
- I did not check the 2PN TaylorF2 coefficients against a second source beyond my own
  recollection of Buonanno et al. (2009); no network access.
- I did not test `report.html` in Firefox or Safari, on a phone, at the `max-width: 720px`
  breakpoint, or in light mode — only headless Chrome at 1400 px in the default (dark) OS
  theme. `prefers-color-scheme: light` is the `:root` default and looks structurally sound
  in the CSS, but I did not render it. I also did not test the `▶ reproducir` button or the
  `velocidad` button by real clicking (no browser automation); I exercised the scrub input's
  `input` handler, which is the same `setFrac`/`onFrame` path both buttons drive.
- I did not audit `src/gwlens/chirp.py`, `system.py` or `units.py` line by line — the first
  review did, their outputs are consistent with everything downstream, and my time went to
  the new material and the HTML.
- I did not measure figure claims off the PNGs pixel by pixel (e.g. that
  `caseA_F_of_f.png`'s band really spans 0.79–1.27); I verified that analytically instead.
- I have no way to verify authorship or timeline claims in `wiki/log.md`; I verified only
  that the current code matches every account it gives.

## Recommended order of fixes

1. Finding 1 — conjugate `F` where it is applied to numpy-convention data (or flip the two
   `exp(+i…)` in `F_geometric_optics` and the equivalent in `F_point_lens`/`F_radial_1d`/
   `F_bruteforce_2d`, keeping all four consistent), correct `wiki/conventions.md`'s Fourier
   paragraph, **add a time-domain regression test** ("second image arrives later, at
   `4GM/c³·ΔT`, with amplitude `√|μ₋|`"), regenerate Case A in WSL, and update every
   document that quotes `peak_strain_amplification` or the "smeared merger" story.
2. Findings 2–4 — the HTML panels: export at full resolution or draw the already-exported
   envelopes; widen `y_max` for the animation; either use `θ_E(t)` for the marker or state
   the fixed-`θ_E` projection in the panel header and the figure title.
3. Findings 5–7 — write the retro-lensing estimate that Chapter 5 claims §4.3 contains (or
   drop the cross-reference), replace the two hard-coded `Capítulo 5`s with `\ref`, and fix
   §2.1's curvature-scale criterion.
4. Findings 12–14, 16 — the reproducibility and staleness fixes; each is a few lines.
