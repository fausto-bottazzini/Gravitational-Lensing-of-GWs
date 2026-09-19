# What depends on what, and what is closed

This page exists to stop finished work from being reopened. It is the first
thing to read before editing anything in this repository.

## The rule

- A node marked **CLOSED** is not edited again. Not tidied, not improved, not
  "while I am here". Only the author reopens it, and says so explicitly.
- Closing is per **folder**, and the date is recorded below.
- If a change to an open node would force a change inside a closed one, that is
  not a licence to make it: stop and ask first.
- Reviews and audits may *report* on closed folders. Acting on those reports is
  a separate decision.
- An arrow `A --> B` means **B is built from, or asserts something about, A**.
  So a change to A can invalidate B, and the arrows are the blast radius.

## Status

| Folder | Status | Closed on | What it contains |
|---|---|---|---|
| `theory/` | **CLOSED** | 2026-09-18 | `theory.tex`/`.pdf`, `OUTLINE.md`, `refs.bib`, `paraxial_validity.py` + its JSON, `pulsos_esquema.py` + its PDF |
| `src/gwlens/` | **CLOSED** | 2026-09-19 | the physics library and the pinned system |
| `tests/` | **CLOSED** | 2026-09-19 | the checks and their `CHECKS*.json` records |
| `cases/` | **CLOSED** | 2026-09-19 | both `run.py`, the figures, `RESULTS.md`, `provenance/`, `FIGURES.md` |
| `report/` | **CLOSED** | 2026-09-19 | `report.tex`/`.pdf`/`.bib`, `report.html`, `update_report_data.py` |
| `wiki/` | open | — | these notes, plus the course brief |
| `bibliography/` | **CLOSED** | 2026-09-18 | `README.md` (mirrors `refs.bib`) |
| root | open | — | `README.md`, `reproduce.sh`, `requirements.txt` |

## The graph

```mermaid
flowchart LR
  subgraph SRC["src/gwlens · CLOSED"]
    SYS["system.py<br/><i>the one pinned system</i>"]
    LIB["waveoptics · geometry<br/>doppler · chirp<br/>taylorf2 · imr_waveform"]
  end

  subgraph TST["tests/ · CLOSED"]
    TESTS["test_*.py"]
    CHK["CHECKS*.json"]
  end

  subgraph CAS["cases/ · CLOSED"]
    RUNA["case_A_chirp/run.py"]
    RUNB["case_B_monochromatic/run.py"]
    NUM["provenance/numbers.json"]
    FIGS["the PNGs<br/>+ caseB_one_period.html"]
    ANIM["case*_animation_data.json"]
    PROV["RESULTS.md · claims.yaml<br/>FIGURES.md"]
  end

  subgraph THE["theory/ · CLOSED"]
    PVPY["paraxial_validity.py"]
    PVJS["paraxial_validity_numbers.json"]
    PUPY["pulsos_esquema.py"]
    PUPDF["pulsos_esquema.pdf"]
    BIB["refs.bib"]
    TEX["theory.tex"]
    TPDF["theory.pdf"]
    OUT["OUTLINE.md"]
  end

  subgraph REP["report/ · CLOSED"]
    UPD["update_report_data.py"]
    RTEX["report.tex"]
    RPDF["report.pdf"]
    RHTML["report.html"]
  end

  SYS --> LIB
  LIB --> TESTS --> CHK
  LIB --> RUNA & RUNB
  RUNA & RUNB --> NUM & FIGS & ANIM
  NUM --> PROV
  CHK --> PROV
  FIGS --> PROV

  SYS --> PVPY --> PVJS --> TEX
  PUPY --> PUPDF --> TEX
  BIB --> TEX --> TPDF
  TEX --> OUT

  NUM --> RTEX
  FIGS --> RTEX
  BIB --> RTEX --> RPDF
  TEX -. cited by chapter TITLE .-> RTEX

  ANIM --> UPD --> RHTML
  NUM --> RHTML
  FIGS --> RHTML

  classDef closed fill:#e8f5e9,stroke:#2e7d32,stroke-width:2px;
  class PVPY,PVJS,PUPY,PUPDF,BIB,TEX,TPDF,OUT,SYS,LIB,TESTS,CHK,
        RUNA,RUNB,NUM,FIGS,ANIM,PROV,UPD,RTEX,RPDF,RHTML closed;
```

## Edges that have actually caused damage

These are the non-obvious ones. Each is here because it broke something.

- **`theory.tex` → `report.tex`, by chapter TITLE.** The `\teoria{}` macro cites
  theory by the chapter's words, not by a label. Renaming a chapter in a closed
  `theory/` silently falsifies citations in an open `report/`, and nothing
  errors. One of them already pointed at the wrong chapter.
- **`refs.bib` is shared** by `theory.tex` *and* `report.tex`. An entry can look
  orphaned in one and be cited in the other — a review got this wrong and
  "fixed" two non-orphans.
- **`case*_animation_data.json` → `report.html`.** The page embeds the data, it
  does not read it. Re-running a case and not running
  `report/update_report_data.py` leaves the deck showing stale numbers with no
  sign of it.
- **`caseB_pattern_only.png` → `report.html`.** The deck positions a moving
  marker over it *by percentage of the image*, so that PNG is saved without
  `bbox_inches="tight"`. Change its margins and the marker drifts off the
  pattern silently.
- **`numbers.json` → everything that quotes a number.** `RESULTS.md`,
  `claims.yaml`, `FIGURES.md`, `report.tex` and `report.html` all transcribe
  values by hand. Every re-run needs those swept; this has gone stale
  repeatedly.
- **Native Windows → Case A.** Running `reproduce.sh checks`/`cases` without
  `pycbc` overwrites Case A's committed IMRPhenomD results with the TaylorF2
  fallback. Case A is regenerated from the WSL venv, never natively.
- **`system.py` is upstream of both worlds.** It feeds the cases *and*
  `theory/paraxial_validity.py`, whose JSON is transcribed into theory's Table
  4.1. Changing the pinned system reaches into a closed folder.

## Closing a folder

1. Its own files agree with each other and with what produced them.
2. Everything downstream of it in the graph is still true — check the arrows
   leaving it.
3. Record it in the status table with the date, and in `log.md`.
