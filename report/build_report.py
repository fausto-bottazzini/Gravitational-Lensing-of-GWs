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

template = (HERE / "report_template.html").read_text(encoding="utf-8")

case_a = (ROOT / "cases/case_A_chirp/caseA_animation_data.json").read_text(encoding="utf-8")
case_b = (ROOT / "cases/case_B_monochromatic/caseB_animation_data.json").read_text(encoding="utf-8")

out = template.replace("/*__CASE_A_DATA__*/", case_a).replace("/*__CASE_B_DATA__*/", case_b)

(HERE / "report.html").write_text(out, encoding="utf-8")
print(f"wrote {HERE / 'report.html'} ({len(out)/1024:.0f} KB)")
