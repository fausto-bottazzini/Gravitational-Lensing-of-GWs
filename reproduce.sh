#!/usr/bin/env bash
# Reproduce every result in this repository from a clean checkout.
#
#   ./reproduce.sh          # everything: checks, both cases
#   ./reproduce.sh checks   # just the correctness checks (fast, ~30s)
#   ./reproduce.sh cases    # just the two case studies (figures + numbers)
#   ./reproduce.sh theory   # rebuild theory/theory.pdf (needs a LaTeX install)
#   ./reproduce.sh report   # rebuild report/report.pdf (needs a LaTeX install;
#                           #   report/report.html needs no rebuild step, it
#                           #   just references the same PNGs cases/ writes)
#
# Nothing here is required to already exist: a venv is created, pinned
# dependencies (requirements.txt) are installed into it, then everything
# runs from there. See README.md for what each step produces and where.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

STEP="${1:-all}"

if [ ! -d .venv ]; then
  echo "== creating .venv =="
  python3 -m venv .venv
fi
# shellcheck disable=SC1091
source .venv/bin/activate 2>/dev/null || source .venv/Scripts/activate

echo "== installing pinned dependencies =="
python3 -m pip install --quiet -r requirements.txt

run_checks() {
  echo "== running correctness checks (tests/) =="
  for f in tests/test_*.py; do
    echo "--- $f ---"
    python3 "$f"
  done
}

run_cases() {
  echo "== Case A: chirp, static lens =="
  python3 cases/case_A_chirp/run.py
  echo "== Case B: quasi-monochromatic, moving lens =="
  python3 cases/case_B_monochromatic/run.py
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
  python3 report/build_report.py
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
