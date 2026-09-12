"""Assembles report/report.html from report_template.html by inlining the
two case scripts' animation-data JSON as JS `<script>` blocks (fetch() of a
local JSON file is blocked by CORS when a page is opened directly as a
file:// URL, which this report is meant to be -- inlining sidesteps that
entirely, no local web server needed to view the page).

Run whenever cases/*/caseX_animation_data.json changes:
    python3 report/build_report.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

SOURCES = {
    "/*__CASE_A_DATA__*/": ROOT / "cases/case_A_chirp/caseA_animation_data.json",
    "/*__CASE_B_DATA__*/": ROOT / "cases/case_B_monochromatic/caseB_animation_data.json",
}


def main():
    template = (HERE / "report_template.html").read_text(encoding="utf-8")
    out = template

    for placeholder, path in SOURCES.items():
        raw = path.read_text(encoding="utf-8")

        # Validate as strict JSON before inlining, not after. Inlined into a
        # <script> block this text is parsed as JAVASCRIPT, where NaN and
        # Infinity are perfectly good identifiers -- so a case script that
        # emits them produces a page that renders fine and a .json file that
        # no JSON parser will read. That exact bug shipped once here (a
        # division by an exactly-zero D_LS at the front/back crossings; see
        # wiki/log.md), and it was found by eye rather than by anything that
        # would have caught it a second time.
        #
        # `parse_constant` is not optional here: bare `json.loads` ACCEPTS
        # NaN/Infinity/-Infinity as a documented CPython extension to the
        # spec, so the obvious version of this guard silently passes the one
        # input it exists to reject. Verified by injecting a NaN and watching
        # the build succeed before this hook was added.
        def _reject(name):
            raise ValueError(f"non-finite value {name!r} (valid JavaScript, not valid JSON)")

        try:
            json.loads(raw, parse_constant=_reject)
        except ValueError as exc:
            raise SystemExit(
                f"ERROR: {path.relative_to(ROOT)} is not valid JSON ({exc}).\n"
                "  Inlining it anyway would still produce a working-looking page "
                "(JS accepts NaN/Infinity, JSON does not).\n"
                "  Re-run the case script that writes it and fix the source of "
                "the non-finite value."
            ) from exc

        # And check the placeholder is actually there. A silent no-op replace
        # would emit `window.__CASE_A__ = ;` -- a syntax error that kills the
        # whole <script>, and so every visualization on the page, while the
        # build still reports success.
        if placeholder not in out:
            raise SystemExit(
                f"ERROR: placeholder {placeholder} not found in report_template.html.\n"
                "  Nothing would have been inlined and the page's JS would be "
                "broken, silently."
            )
        out = out.replace(placeholder, raw)

    (HERE / "report.html").write_text(out, encoding="utf-8")
    print(f"wrote {HERE / 'report.html'} ({len(out)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
