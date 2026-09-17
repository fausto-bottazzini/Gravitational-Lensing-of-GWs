"""Refresh the animation data embedded in report/report.html, in place.

report.html is the only HTML file: you edit it, you open it, there is no
generated copy. What this script does is narrow -- it replaces the contents of
the two data blocks

    <script id="datos-caso-A"> ... </script>
    <script id="datos-caso-B"> ... </script>

with whatever the case scripts last wrote to their `case*_animation_data.json`,
and leaves every other byte alone. Run it after re-running a case:

    python3 report/update_report_data.py

The data is embedded rather than fetched because the page is meant to open as
a `file://` URL, where `fetch()` of a local JSON file is blocked by CORS.
Embedding sidesteps that completely: no local web server needed to view it.

(This replaces an earlier arrangement where report_template.html was the
source and report.html a build artefact. Two nearly identical HTML files were
confusing and nothing was gained: the edits all happen in one of them, and the
data blocks can be swapped in place just as easily.)
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PAGE = HERE / "report.html"

SOURCES = {
    "datos-caso-A": ROOT / "cases/case_A_chirp/caseA_animation_data.json",
    "datos-caso-B": ROOT / "cases/case_B_monochromatic/caseB_animation_data.json",
}
VARIABLES = {"datos-caso-A": "__CASE_A__", "datos-caso-B": "__CASE_B__"}


def main():
    html = PAGE.read_text(encoding="utf-8")
    original = html

    for block_id, path in SOURCES.items():
        raw = path.read_text(encoding="utf-8")
        # Parsed before it is embedded, not after. Inside a <script> block this
        # text is read as JAVASCRIPT, where NaN and Infinity are ordinary
        # identifiers and would slip through silently; json.loads rejects them.
        data = json.loads(raw)
        payload = json.dumps(data, separators=(",", ":"))
        if "</script" in payload:
            sys.exit("ERROR: %s contains '</script', which would close the "
                     "block early." % path)

        pattern = re.compile(
            r'(<script id="%s">).*?(</script>)' % re.escape(block_id), re.S)
        if not pattern.search(html):
            sys.exit("ERROR: no <script id=\"%s\"> block in %s.\n"
                     "The page must keep those two blocks for this script to "
                     "have somewhere to put the data." % (block_id, PAGE.name))
        html = pattern.sub(
            lambda m: "%swindow.%s = %s;%s" % (
                m.group(1), VARIABLES[block_id], payload, m.group(2)),
            html, count=1)
        print("  %-14s <- %s (%d KB)" % (block_id, path.name, len(payload) // 1024))

    if html == original:
        print("sin cambios")
    else:
        PAGE.write_text(html, encoding="utf-8")
    print("%s (%d KB)" % (PAGE, len(html) // 1024))


if __name__ == "__main__":
    main()
