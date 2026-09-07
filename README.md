# Wave-optics gravitational lensing of gravitational waves in a hierarchical triple

Final project for *Gravitational Waves and AI-Assisted Research* (course
repository: [matiaszaldarriaga/GW-AI-course](https://github.com/matiaszaldarriaga/GW-AI-course)).

**The question**: an inner compact binary (the GW source) orbits a third,
more massive body (the lens) far enough away that tidal effects are
negligible. How is the GW lensed, in wave optics, and how does the answer
differ between the late inspiral (a fast chirp, lens effectively frozen) and
an earlier, wider phase of the same binary (quasi-monochromatic, lens
genuinely moving)?

**Where to start reading**: [`report/report.html`](report/report.html) or
[`report/report.pdf`](report/report.pdf) — the results, with figures, meant
to be read start to finish. For the full derivation (wave equation to scalar
field, the point-lens closed form and its checks, the triple's geometry),
read [`theory/theory.pdf`](theory/theory.pdf), written as a self-contained
textbook chapter. Everything below this line is either a map of the repo or
about *reproducing* those two documents, not about the physics itself.

## Repository map

| | |
|---|---|
| `report/` | the results: `report.html` (interactive, presentation-style) and `report.pdf` (full detail), the two documents this project is presented from. `report_template.html` + `build_report.py` generate `report.html` by inlining `cases/*/caseX_animation_data.json` |
| `theory/` | the derivation, as a textbook chapter: `theory.tex` -> `theory.pdf` |
| `src/gwlens/` | the physics library: wave optics, orbits, the inner-binary waveform (`imr_waveform.py`, `taylorf2.py`, `chirp.py`), the one pinned system |
| `cases/case_A_chirp/` | the static-lens case: script, figures, `RESULTS.md`, `provenance/` |
| `cases/case_B_monochromatic/` | the moving-lens case: script, figures, `RESULTS.md`, `provenance/` |
| `tests/` | every correctness check referenced from `theory.tex` and the case `RESULTS.md` files |
| `wiki/` | this project's own working notes: `conventions.md` (notation/units, read this before the code), `log.md` (append-only decision log, including every bug caught and how), `todo.md` |

## What each case actually shows

Both cases are the *same* hierarchical triple (`src/gwlens/system.py`: a
20+15 M<sub>&#8857;</sub> inner binary, a 5&times;10<sup>4</sup> M<sub>&#8857;</sub>
lens, a 4-day, near-edge-on outer orbit), at two epochs of its inspiral:

- **Case A** (`cases/case_A_chirp/`): late inspiral, the chirp sweeps 10 Hz
  to merger in 15 s — far too fast for the outer orbit to move, so the lens
  is frozen at one impact parameter for the whole signal. Shows the
  frequency-domain interference fringes across the band, the extended
  amplification pattern on the source plane (and how finely it's resolved
  changes with frequency), and the lensed vs. unlensed chirp.
- **Case B** (`cases/case_B_monochromatic/`): earlier, wider inspiral,
  f=0.05 Hz and nearly constant over a 24-day (6 outer period) observation.
  The outer orbital motion *is* resolved: the impact parameter sweeps
  through the same diffraction pattern once per period, producing periodic
  "repeated lensing" amplification pulses — and exactly half of each orbit
  turns out to be unlensed entirely (the source is in front of the lens,
  not behind it).

Each `RESULTS.md` states every number with a pointer to the function that
produced it and the check behind it; nothing there is asserted without both.

## Reproduce everything

```
git clone <this repo>
cd gw-lensing-hierarchical-triple
./reproduce.sh          # checks + both cases (~1-2 min, creates .venv)
./reproduce.sh theory   # rebuilds theory/theory.pdf (needs a LaTeX install)
./reproduce.sh report   # rebuilds report/report.pdf (needs a LaTeX install)
```

No step needs anything pre-installed beyond Python 3.11+ and (for the
`theory`/`report` steps only) a working LaTeX distribution (MiKTeX/TeX Live,
with `pdflatex` and `bibtex` on `PATH`) — `reproduce.sh` creates its own venv
and installs the four pinned dependencies (`requirements.txt`) into it.
Nothing in `tests/` needs LaTeX.

Each case script (`cases/case_A_chirp/run.py`,
`cases/case_B_monochromatic/run.py`) is also runnable on its own and takes
well under a minute; each `tests/test_*.py` is runnable on its own too and
prints a PASS/FAIL line per check.

### One waveform, two possible sources — read this before Case A surprises you

Case A's unlensed gravitational waveform (`src/gwlens/imr_waveform.py`)
tries, in order:

1. **`pycbc`'s `IMRPhenomD`** (Khan et al. 2016) — a published,
   NR-calibrated approximant with a real merger and ringdown, exactly what
   LIGO/Virgo pipelines use. This is what the figures and numbers
   *committed in this repository* were generated with.
2. If `pycbc` cannot be imported: **`src/gwlens/taylorf2.py`**, a
   hand-built, from-scratch, restricted 2PN inspiral — validated on its own
   terms (`tests/test_taylorf2.py`) but with no merger or ringdown (cut off
   at the ISCO frequency, where the post-Newtonian approximation itself
   stops being valid).

**`pycbc` is not in `requirements.txt`** and is not expected to install on
native Windows in any reasonable time (`lalsuite`'s dependency resolution
hangs there — see `wiki/log.md`). It installs cleanly with a plain
`pip install pycbc` on Linux and macOS, **including WSL2** (this is exactly
how this repository's own Case A results were produced, from a Windows
machine, without needing native Windows support at all):

```
# inside WSL2 Ubuntu (or any Linux/macOS shell):
cd gw-lensing-hierarchical-triple      # the SAME repo, e.g. via /mnt/c/... on WSL
python3 -m venv .venv-wsl
source .venv-wsl/bin/activate
pip install pycbc mpmath                # do NOT use -r requirements.txt here:
                                         # those numpy/scipy versions are pinned
                                         # for the Windows fallback venv and can
                                         # be too old for pycbc's own build/wheel
                                         # requirements on a newer Python (as on
                                         # a fresh WSL2 Ubuntu install) -- pip
                                         # resolves compatible modern versions on
                                         # its own; see wiki/log.md
python3 cases/case_A_chirp/run.py      # now uses IMRPhenomD automatically
```

Nothing else changes: the same script, the same command, on Linux/macOS/WSL2
picks up `pycbc` automatically and produces the full inspiral-merger
-ringdown result; on native Windows without `pycbc` it falls back to
TaylorF2 automatically, with a clear label
(`numbers.json`'s `waveform_source` field) recording which one actually
ran. Case B, `tests/`, and everything else never touch `pycbc` at all.

**Caution:** afterwards, running `reproduce.sh checks` (or `cases`/`all`) on
*native* Windows will silently regenerate `tests/CHECKS_imr_waveform.json`
and Case A's own `numbers.json`/figures with the TaylorF2 fallback, since
that venv has no `pycbc` — overwriting the WSL-produced, pycbc-based results
this repository actually ships. If you want to reproduce the exact committed
Case A numbers and figures (including the ~3.4 s second-image echo), do it
from the WSL venv above, and do it last, right before comparing against what
is committed.

## Reproducibility test this repo was held to

Per the assignment (`GW-AI-course/final-project.html`): *hand it to a fresh
agent that knows nothing, ask it to reproduce a result, and see how far it
gets.* This repo was checked that way three times over (independent
fresh-agent passes with no prior context, each fixing what the last one
found); the same test can be repeated by anyone with a copy of this repo and
an agent.
