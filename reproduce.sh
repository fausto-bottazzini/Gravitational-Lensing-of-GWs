#!/usr/bin/env bash
# Reproduce every result in this repository from a clean checkout.
#
#   ./reproduce.sh          # everything: checks, both cases
#   ./reproduce.sh checks   # just the correctness checks (fast, ~30s)
#   ./reproduce.sh cases    # just the two case studies (figures + numbers)
#   ./reproduce.sh theory   # rebuild theory/theory.pdf (needs a LaTeX install)
#   ./reproduce.sh report   # rebuild report/report.pdf and report.html (needs
#                           #   a LaTeX install for the PDF; the HTML step
#                           #   needs no LaTeX, just the venv's python)
#
# Nothing here is required to already exist: a venv is created, pinned
# dependencies (requirements.txt) are installed into it, then everything
# runs from there via the venv's OWN python executable (found explicitly
# below, not assumed to be on PATH after "activation" -- on Windows a venv
# only ever gets a `python.exe`, never a `python3.exe`/`python3`, so a bare
# `python3` call after activating still silently runs the SYSTEM Python
# instead of the venv's, defeating the isolation this script claims to
# give you; caught by an independent review, see wiki/log.md).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

STEP="${1:-all}"

if [ ! -d .venv ]; then
  echo "== creating .venv =="
  python3 -m venv .venv
fi

# Find the venv's own interpreter directly -- do not rely on `python3`/
# `python` resolving correctly on PATH after sourcing an activate script;
# that is exactly what silently broke isolation on Windows (see above).
if [ -x .venv/bin/python3 ]; then
  PY=.venv/bin/python3
elif [ -x .venv/Scripts/python.exe ]; then
  PY=.venv/Scripts/python.exe
elif [ -x .venv/bin/python ]; then
  PY=.venv/bin/python
else
  echo "ERROR: could not find a python executable inside .venv/" >&2
  exit 1
fi
echo "== using venv interpreter: $PY =="

echo "== installing pinned dependencies =="
"$PY" -m pip install --quiet -r requirements.txt

run_checks() {
  echo "== running correctness checks (tests/) =="
  for f in tests/test_*.py; do
    echo "--- $f ---"
    "$PY" "$f"
  done
}

run_cases() {
  echo "== Case A: chirp, static lens =="
  "$PY" cases/case_A_chirp/run.py
  echo "== Case B: quasi-monochromatic, moving lens =="
  "$PY" cases/case_B_monochromatic/run.py
}

run_theory() {
  echo "== rebuilding theory/theory.pdf (requires pdflatex + bibtex) =="
  ( cd theory && \
    pdflatex -interaction=nonstopmode -halt-on-error theory.tex && \
    bibtex theory && \
    pdflatex -interaction=nonstopmode -halt-on-error theory.tex && \
    pdflatex -interaction=nonstopmode -halt-on-error theory.tex && \
    rm -f *.aux *.bbl *.blg *.log *.out *.toc )
}

run_report() {
  echo "== rebuilding report/report.pdf (requires pdflatex + bibtex; needs cases/ figures to already exist) =="
  ( cd report && \
    pdflatex -interaction=nonstopmode -halt-on-error report.tex && \
    bibtex report && \
    pdflatex -interaction=nonstopmode -halt-on-error report.tex && \
    pdflatex -interaction=nonstopmode -halt-on-error report.tex && \
    rm -f *.aux *.bbl *.blg *.log *.out *.toc )
  echo "== rebuilding report/report.html (inlines cases/*/caseX_animation_data.json) =="
  "$PY" report/build_report.py
}

case "$STEP" in
  checks) run_checks ;;
  cases)  run_cases ;;
  theory) run_theory ;;
  report) run_report ;;
  all)    run_checks; run_cases ;;
  *) echo "usage: $0 [checks|cases|theory|report|all]"; exit 1 ;;
esac

echo "== done =="
