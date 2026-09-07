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
textbook chapter. Everything below this line is about *reproducing* those
two documents, not about the physics itself.

## Reproduce everything

```
git clone <this repo>
cd gw-lensing-hierarchical-triple
./reproduce.sh          # checks + both cases (~1-2 min, creates .venv)
./reproduce.sh theory   # rebuilds theory/theory.pdf (needs a LaTeX install)
```

No step needs anything pre-installed beyond Python 3.11+ and (for the
`theory` step only) a working LaTeX distribution (MiKTeX/TeX Live, with
`pdflatex` and `bibtex` on `PATH`) — `reproduce.sh` creates its own venv and
installs the four pinned dependencies (`requirements.txt`) into it. Nothing
in `cases/` or `tests/` needs LaTeX.

Each case script (`cases/case_A_chirp/run.py`,
`cases/case_B_monochromatic/run.py`) is also runnable on its own and takes
well under a minute; each `tests/test_*.py` is runnable on its own too and
prints a PASS/FAIL line per check.

## Repository map

| | |
|---|---|
| `report/` | the results: `report.html` and `report.pdf`, the two documents this project is presented from |
| `theory/` | the derivation, as a textbook chapter: `theory.tex` -> `theory.pdf` |
| `src/gwlens/` | the physics library: wave optics, orbits, the inner-binary waveform, the one pinned system |
| `cases/case_A_chirp/` | the static-lens case: script, figures, `RESULTS.md`, `provenance/` |
| `cases/case_B_monochromatic/` | the moving-lens case: script, figures, `RESULTS.md`, `provenance/` |
| `tests/` | every correctness check referenced from `theory.tex` and the case `RESULTS.md` files |
| `checks/independent_review/` | a fresh agent's independent check of this repository, run with no prior context |
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

## Reproducibility test this repo was held to

Per the assignment (`GW-AI-course/final-project.html`): *hand it to a fresh
agent that knows nothing, ask it to reproduce a result, and see how far it
gets.* `checks/independent_review/` is exactly that, run once already; the
same test can be repeated by anyone with a copy of this repo and an agent.
