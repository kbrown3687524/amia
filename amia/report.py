"""Generate a self-contained HTML index for an AMIA pipeline run."""

import csv
import html
from pathlib import Path


REPORT_NAME = "amia_report.html"


def _escape(value):
    return html.escape(str(value), quote=True)


def _file_link(path, output_dir):
    relative_path = path.relative_to(output_dir).as_posix()
    return f'<a href="{_escape(relative_path)}">{_escape(relative_path)}</a>'


def _csv_preview(path, output_dir, limit=12):
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.reader(stream))
    except (OSError, UnicodeError, csv.Error):
        return ""
    if not rows:
        return "<p class=muted>Empty CSV file.</p>"
    header, data = rows[0], rows[1:limit + 1]
    header_html = "".join(f"<th>{_escape(value)}</th>" for value in header)
    body_html = "".join(
        "<tr>" + "".join(f"<td>{_escape(value)}</td>" for value in row) + "</tr>"
        for row in data
    )
    note = ""
    if len(rows) - 1 > limit:
        note = f"<p class=muted>Showing {limit} of {len(rows) - 1} rows.</p>"
    return f"<table><thead><tr>{header_html}</tr></thead><tbody>{body_html}</tbody></table>{note}"


def generate_report(output_dir, config_path, status, completed_steps, error=None):
    """Write and return the path to the run-level HTML report."""
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / REPORT_NAME
    files = sorted(
        path for path in output_dir.rglob("*")
        if path.is_file() and path.resolve() != report_path.resolve()
    )
    html_files = [path for path in files if path.suffix.lower() == ".html"]
    csv_files = [path for path in files if path.suffix.lower() == ".csv"]
    other_files = [path for path in files if path not in html_files and path not in csv_files]
    status_class = "success" if status == "complete" else "failed"
    error_html = f"<pre class=error>{_escape(error)}</pre>" if error else ""
    completed_html = "".join(f"<li>{_escape(step)}</li>" for step in completed_steps)

    def artifact_list(paths):
        if not paths:
            return "<p class=muted>No files in this category.</p>"
        return "<ul class=artifacts>" + "".join(
            f"<li>{_file_link(path, output_dir)} <span class=muted>({path.stat().st_size:,} bytes)</span></li>"
            for path in paths
        ) + "</ul>"

    csv_sections = "".join(
        f"<details><summary>{_file_link(path, output_dir)}</summary>{_csv_preview(path, output_dir)}</details>"
        for path in csv_files
    ) or "<p class=muted>No CSV tables were produced.</p>"
    html_sections = "".join(
        f"<li>{_file_link(path, output_dir)}</li>" for path in html_files
    ) or "<li class=muted>No standalone HTML reports were produced.</li>"

    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>AMIA Pipeline Report</title>
<style>
:root {{ color-scheme: light; --ink:#1f2933; --muted:#66727f; --line:#d9e1e8; --blue:#1769aa; --green:#18794e; --red:#b42318; --paper:#ffffff; --wash:#f3f6f8; }}
* {{ box-sizing:border-box; }} body {{ margin:0; background:var(--wash); color:var(--ink); font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif; }}
main {{ max-width:1120px; margin:0 auto; padding:32px 20px 56px; }} header {{ background:var(--ink); color:white; padding:28px; border-radius:10px; }}
h1 {{ margin:0 0 8px; font-size:clamp(1.7rem,4vw,2.5rem); }} header p {{ margin:0; color:#d7e0e8; }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(220px,1fr)); gap:16px; margin:20px 0; }} section, .metric {{ background:var(--paper); border:1px solid var(--line); border-radius:8px; padding:20px; }}
h2 {{ margin:0 0 14px; font-size:1.15rem; }} h3 {{ margin:20px 0 8px; font-size:1rem; }}
.metric strong {{ display:block; font-size:1.45rem; }} .label,.muted {{ color:var(--muted); }}
.badge {{ display:inline-block; padding:3px 10px; border-radius:999px; font-weight:700; }} .success {{ color:var(--green); }} .failed {{ color:var(--red); }}
ul {{ padding-left:20px; }} .artifacts {{ columns:2 320px; }} a {{ color:var(--blue); }}
table {{ width:100%; border-collapse:collapse; display:block; overflow:auto; }} th,td {{ border:1px solid var(--line); padding:7px 9px; text-align:left; white-space:nowrap; }} th {{ background:#eaf0f5; }}
details {{ margin:10px 0; border:1px solid var(--line); padding:10px; border-radius:6px; background:#fbfcfd; }} summary {{ cursor:pointer; }} pre {{ white-space:pre-wrap; overflow:auto; }} .error {{ color:var(--red); background:#fff1f0; padding:12px; border-radius:6px; }}
@media (max-width:600px) {{ main {{ padding:18px 12px 40px; }} header, section, .metric {{ padding:16px; }} }}
</style>
</head>
<body><main>
<header><h1>AMIA Pipeline Report</h1><p>Generated from <strong>{_escape(config_path)}</strong></p></header>
<div class=grid>
<div class=metric><span class=label>Status</span><strong class={status_class}>{_escape(status.title())}</strong></div>
<div class=metric><span class=label>Completed steps</span><strong>{len(completed_steps)}</strong></div>
<div class=metric><span class=label>Output files</span><strong>{len(files)}</strong></div>
<div class=metric><span class=label>Report generated</span><strong>{_escape(report_path.name)}</strong></div>
</div>
<section><h2>Pipeline Progress</h2><ol>{completed_html or '<li class=muted>No steps completed.</li>'}</ol>{error_html}</section>
<section><h2>Existing Reports</h2><ul>{html_sections}</ul></section>
<section><h2>CSV Results Preview</h2>{csv_sections}</section>
<section><h2>All Artifacts</h2><h3>Data tables</h3>{artifact_list(csv_files)}<h3>HTML reports</h3>{artifact_list(html_files)}<h3>Other output</h3>{artifact_list(other_files)}</section>
</main></body></html>"""
    report_path.write_text(document, encoding="utf-8")
    return report_path
