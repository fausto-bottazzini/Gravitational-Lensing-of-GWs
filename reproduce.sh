#!/usr/bin/env bash
# Reproduce every result in this repository from a clean checkout.
#
#   ./reproduce.sh          # everything: checks, both cases (incl. report.html)
#   ./reproduce.sh checks   # just the correctness checks (fast, ~30s)
#   ./reproduce.sh cases    # the two case studies (figures + numbers) AND
#                           #   report.html, which inlines their output --
#                           #   needs no LaTeX, just the venv's python
#   ./reproduce.sh html     # just report.html, from whatever case data
#                           #   already exists on disk
#   ./reproduce.sh theory   # rebuild theory/theory.pdf (needs a LaTeX install)
#   ./reproduce.sh report   # rebuild report/report.pdf AND report.html
#                           #   (needs a LaTeX install, for the PDF)
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

# requirements.txt is pinned, and the pins have wheels through Python 3.12
# only: numpy 1.26.4 and scipy 1.12.0 both predate cp313. On a newer
# interpreter pip falls back to building them from source, which spends ten
# minutes on numpy and then dies in scipy for want of a system OpenBLAS -- the
# first wall a fresh agent hits on a current distro (Ubuntu 26.04 ships 3.14
# and nothing older). So ask the interpreter first, and let pip resolve modern
# versions when the pins cannot apply. That path is not a degraded fallback:
# an independent reproduction on Python 3.14 with numpy 2.3.5 / scipy 1.16.3 /
# matplotlib 3.11.2 returned 44/44 checks and Case A bit-identical, figures
# included (wiki/log.md).
PYVER=$("$PY" -c 'import sys; print("%d%02d" % sys.version_info[:2])')
PYTAG=$("$PY" -c 'import sys; print("%d.%d" % sys.version_info[:2])')
if [ "$PYVER" -ge 313 ]; then
  echo "== Python $PYTAG is newer than the pins, resolving dependencies unpinned =="
  "$PY" -m pip install --quiet numpy scipy matplotlib mpmath
else
  echo "== installing pinned dependencies =="
  "$PY" -m pip install --quiet -r requirements.txt
fi

# Case A's committed numbers and figures come from pycbc's IMRPhenomD, which
# does not install on native Windows (README.md, "One waveform, two possible
# sources"). Without it the fallback runs, silently, and overwrites what the
# repository ships with an inspiral-only result: no merger, no ringdown, and
# no second-image echo. The fallback is legitimate and the script should not
# refuse to run -- but it should not be quiet about it either, and it should
# say how to undo it.
warn_if_no_pycbc() {
  if ! "$PY" -c "import pycbc" >/dev/null 2>&1; then
    echo
    echo "!! ---------------------------------------------------------------"
    echo "!! pycbc is NOT importable in this venv, so Case A used the"
    echo "!! TaylorF2 fallback: inspiral only, no merger, no ringdown, and"
    echo "!! no second-image echo. Anything just written under"
    echo "!! cases/case_A_chirp/ or tests/CHECKS_imr_waveform.json is a"
    echo "!! fallback result and is NOT what this repository ships."
    echo "!!"
    echo "!!   to undo:        git checkout tests/ cases/case_A_chirp/"
    echo "!!   to reproduce:   see README.md, 'One waveform, two possible"
    echo "!!                   sources' -- run from a Linux/WSL2 venv."
    echo "!! ---------------------------------------------------------------"
    echo
  fi
}

run_checks() {
  echo "== running correctness checks (tests/) =="
  for f in tests/test_*.py; do
    echo "--- $f ---"
    "$PY" "$f"
  done
  # theory/paraxial_validity.py is a validity check on the system parameters
  # like everything above, not part of the physics pipeline -- but it writes
  # the numbers theory.tex's paraxial table quotes, and that table says "run
  # this script to regenerate them, not by hand". Until 2026-09-07 nothing in
  # this script ran it, so those numbers could have gone stale against
  # system.py with no reproduce path noticing. It lives here rather than in
  # run_theory() because it needs no LaTeX, and a reader without a LaTeX
  # install should still be regenerating it.
  echo "--- theory/paraxial_validity.py ---"
  "$PY" theory/paraxial_validity.py
  warn_if_no_pycbc
}

run_cases() {
  echo "== Case A: chirp, static lens =="
  "$PY" cases/case_A_chirp/run.py
  echo "== Case B: quasi-monochromatic, moving lens =="
  "$PY" cases/case_B_monochromatic/run.py
  run_html
  warn_if_no_pycbc
}

run_html() {
  # No LaTeX needed -- just the venv's python, inlining
  # cases/*/caseX_animation_data.json into report/report.html. Kept as its
  # own step (called from both `cases` and `report` below) so report.html
  # never goes stale relative to freshly-regenerated case data even on a
  # machine with no LaTeX install, which `run_report` alone requires (an
  # independent review caught `all`/`cases` leaving report.html out of sync
  # with the JSON it inlines; see wiki/log.md).
  echo "== rebuilding report/report.html (inlines cases/*/caseX_animation_data.json) =="
  "$PY" report/update_report_data.py
}

run_theory() {
  # theory.tex includes theory/pulsos_esquema.pdf, which this script generates
  # -- a build input, so it is regenerated here rather than trusted to be on
  # disk. Same failure mode paraxial_validity.py had until 2026-09-07: a file
  # theory.tex depends on, with no reproduce path that rebuilds it. It needs
  # no LaTeX, only the venv.
  echo "== regenerating theory/pulsos_esquema.pdf (figure 4.1) =="
  "$PY" theory/pulsos_esquema.py
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
  run_html
}

case "$STEP" in
  checks) run_checks ;;
  cases)  run_cases ;;
  html)   run_html ;;
  theory) run_theory ;;
  report) run_report ;;
  all)    run_checks; run_cases ;;
  *) echo "usage: $0 [checks|cases|html|theory|report|all]"; exit 1 ;;
esac

echo "== done =="
