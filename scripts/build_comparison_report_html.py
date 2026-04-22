#!/usr/bin/env python3
"""构建双体系 comparison HTML 报告。"""

from __future__ import annotations

from html import escape
from pathlib import Path


def build_comparison_html_report(output_dir: Path, comparison_result: dict) -> Path:
    output_dir = Path(output_dir)
    html_path = output_dir / "comparison_report.html"

    summary = comparison_result["summary"]
    tables = comparison_result["tables"]
    artifacts = comparison_result["artifacts"]
    case_a = summary["case_a"]["label"]
    case_b = summary["case_b"]["label"]
    comparison_mode = summary.get("comparison_mode", "generic")
    comparison_context = summary.get("comparison_context", "")
    kicker = _mode_kicker(comparison_mode)
    subtitle = _mode_subtitle(comparison_mode)
    context_html = f'<p style="margin-top:10px;"><strong>Context:</strong> {escape(comparison_context)}</p>' if comparison_context else ""

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(case_a)} vs {escape(case_b)} Comparison Report</title>
  <style>
    :root {{
      --bg: #f3efe6;
      --panel: #fffdf8;
      --card: #ffffff;
      --ink: #23342f;
      --muted: #667772;
      --accent: #355f55;
      --accent-soft: #dce8e2;
      --accent-warm: #c78a52;
      --border: #d8ddd7;
      --shadow: 0 10px 24px rgba(31, 51, 45, 0.08);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: linear-gradient(180deg, #efe8d9 0%, var(--bg) 18%, #f7f4ec 100%);
      color: var(--ink);
      line-height: 1.5;
    }}
    .page {{
      width: min(1320px, calc(100vw - 48px));
      margin: 20px auto 56px;
    }}
    .nav {{
      position: sticky;
      top: 12px;
      z-index: 20;
      display: flex;
      gap: 10px;
      flex-wrap: wrap;
      background: rgba(255, 253, 248, 0.92);
      backdrop-filter: blur(12px);
      border: 1px solid rgba(216, 221, 215, 0.9);
      border-radius: 22px;
      padding: 14px;
      box-shadow: var(--shadow);
    }}
    .nav a {{
      text-decoration: none;
      color: var(--accent);
      border: 1px solid #c9d8d1;
      border-radius: 999px;
      padding: 10px 16px;
      font-weight: 700;
      background: #f8fbf9;
    }}
    .hero {{
      margin-top: 18px;
      background: linear-gradient(145deg, #244f46 0%, #446e65 100%);
      color: #f8fbf9;
      border-radius: 30px;
      padding: 30px;
      box-shadow: var(--shadow);
    }}
    .kicker {{
      display: inline-block;
      padding: 8px 14px;
      border-radius: 999px;
      background: rgba(255,255,255,0.12);
      border: 1px solid rgba(255,255,255,0.18);
      font-size: 13px;
      font-weight: 700;
      letter-spacing: 0.06em;
      text-transform: uppercase;
    }}
    .hero h1 {{
      margin: 18px 0 10px;
      font-size: clamp(40px, 5vw, 64px);
      line-height: 0.95;
    }}
    .hero p {{
      margin: 0;
      max-width: 760px;
      font-size: 18px;
      color: rgba(248, 251, 249, 0.88);
    }}
    .cards {{
      margin-top: 24px;
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 14px;
    }}
    .card {{
      background: rgba(255,255,255,0.12);
      border: 1px solid rgba(255,255,255,0.16);
      border-radius: 20px;
      padding: 16px 18px;
    }}
    .card .label {{
      font-size: 12px;
      text-transform: uppercase;
      letter-spacing: 0.05em;
      opacity: 0.78;
    }}
    .card .value {{
      margin-top: 8px;
      font-size: 26px;
      font-weight: 800;
    }}
    section {{
      margin-top: 22px;
      background: var(--panel);
      border-radius: 26px;
      padding: 24px;
      box-shadow: var(--shadow);
      border: 1px solid rgba(216, 221, 215, 0.72);
    }}
    h2 {{
      margin: 0 0 16px;
      font-size: 28px;
    }}
    .two-col {{
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 18px;
    }}
    .subcard {{
      background: var(--card);
      border: 1px solid var(--border);
      border-radius: 20px;
      padding: 18px;
    }}
    .subcard h3 {{
      margin: 0 0 10px;
      font-size: 18px;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 14px;
    }}
    th, td {{
      border-top: 1px solid var(--border);
      padding: 10px 8px;
      text-align: left;
      vertical-align: top;
    }}
    th {{
      color: var(--muted);
      font-weight: 700;
    }}
    .plot-grid {{
      display: grid;
      grid-template-columns: 1fr;
      gap: 16px;
    }}
    img {{
      width: 100%;
      border-radius: 18px;
      border: 1px solid var(--border);
      background: #ffffff;
    }}
    .downloads {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 12px;
    }}
    .downloads a {{
      display: block;
      padding: 14px 16px;
      border-radius: 16px;
      border: 1px solid var(--border);
      color: var(--accent);
      text-decoration: none;
      font-weight: 700;
      background: #fcfdfa;
    }}
    ul {{
      margin: 0;
      padding-left: 20px;
    }}
    @media (max-width: 980px) {{
      .cards, .two-col, .downloads {{ grid-template-columns: 1fr; }}
    }}
  </style>
</head>
<body>
  <div class="page">
    <nav class="nav">
      <a href="#overview">Overview</a>
      <a href="#comparability">Comparability</a>
      <a href="#identity">Identity</a>
      <a href="#quality-interface">Quality & Interface</a>
      <a href="#flexibility">Flexibility</a>
      <a href="#rrcs">RRCS</a>
      <a href="#interaction">Interaction</a>
      <a href="#downloads">Downloads</a>
    </nav>

    <section class="hero" id="overview">
      <span class="kicker">{escape(kicker)}</span>
      <h1>{escape(case_a)} vs {escape(case_b)}</h1>
      <p>{escape(subtitle)}</p>
      {context_html}
      <div class="cards">
        <div class="card"><div class="label">Comparability</div><div class="value">{escape(summary['comparability']['status'])}</div></div>
        <div class="card"><div class="label">Metrics</div><div class="value">{summary['metric_count']}</div></div>
        <div class="card"><div class="label">Condition A</div><div class="value">{escape(case_a)}</div></div>
        <div class="card"><div class="label">Condition B</div><div class="value">{escape(case_b)}</div></div>
      </div>
    </section>

    <section id="comparability">
      <h2>Comparability Check</h2>
      <ul>{''.join(f'<li>{escape(reason)}</li>' for reason in summary['comparability']['reasons'])}</ul>
      <div class="subcard" style="margin-top:16px;">{_render_metric_table(tables['comparison_rows'], {'quality', 'interface'})}</div>
    </section>

    <section id="identity">
      <h2>Identity Comparison</h2>
      <div class="two-col">
        <div class="subcard">
          <h3>{escape(case_a)}</h3>
          {_render_identity_side(tables['identity_rows'], case_a)}
        </div>
        <div class="subcard">
          <h3>{escape(case_b)}</h3>
          {_render_identity_side(tables['identity_rows'], case_b)}
        </div>
      </div>
    </section>

    <section id="quality-interface">
      <h2>Quality &amp; Interface</h2>
      <div class="plot-grid">
        <img src="{_relpath(output_dir, artifacts['quality_interface_plot'])}" alt="Quality and interface comparison">
      </div>
      <div class="subcard" style="margin-top:16px;">
        <h3>Mechanistic Summary</h3>
        <ul>{''.join(f'<li>{escape(item)}</li>' for item in summary['takeaways'])}</ul>
      </div>
    </section>

    <section id="flexibility">
      <h2>Flexibility Shift</h2>
      <img src="{_relpath(output_dir, artifacts['flexibility_plot'])}" alt="Flexibility comparison">
      <div class="subcard" style="margin-top:16px;">{_render_simple_table(tables['rmsf_region_rows'])}</div>
    </section>

    <section id="rrcs">
      <h2>RRCS Shift</h2>
      <img src="{_relpath(output_dir, artifacts['rrcs_plot'])}" alt="RRCS region comparison">
      <div class="subcard" style="margin-top:16px;">{_render_simple_table(tables['rrcs_region_rows'])}</div>
    </section>

    <section id="interaction">
      <h2>Interaction Landscape Shift</h2>
      <img src="{_relpath(output_dir, artifacts['interaction_plot'])}" alt="Interaction family comparison">
      <div class="subcard" style="margin-top:16px;">{_render_simple_table(tables['interaction_family_rows'])}</div>
    </section>

    <section id="downloads">
      <h2>Downloads</h2>
      <div class="downloads">
        {_render_download_link(output_dir, artifacts['summary_json'], 'Comparison summary JSON')}
        {_render_download_link(output_dir, artifacts['comparison_table_csv'], 'Comparison table CSV')}
        {_render_download_link(output_dir, artifacts['identity_comparison_csv'], 'Identity comparison CSV')}
        {_render_download_link(output_dir, artifacts['rmsf_region_comparison_csv'], 'RMSF region comparison CSV')}
        {_render_download_link(output_dir, artifacts['rrcs_region_comparison_csv'], 'RRCS region comparison CSV')}
        {_render_download_link(output_dir, artifacts['interaction_family_comparison_csv'], 'Interaction family comparison CSV')}
        {_render_download_link(output_dir, artifacts['quality_interface_plot'], 'Quality & interface plot')}
        {_render_download_link(output_dir, artifacts['flexibility_plot'], 'Flexibility plot')}
        {_render_download_link(output_dir, artifacts['rrcs_plot'], 'RRCS comparison plot')}
        {_render_download_link(output_dir, artifacts['interaction_plot'], 'Interaction comparison plot')}
      </div>
    </section>
  </div>
</body>
</html>
"""

    html_path.write_text(html, encoding="utf-8")
    return html_path


def _render_metric_table(rows: list[dict], categories: set[str]) -> str:
    filtered = [row for row in rows if row["category"] in categories]
    return _render_simple_table(filtered)


def _render_identity_side(rows: list[dict], label: str) -> str:
    items = [f"<p><strong>{escape(row['field'])}</strong><br>{escape(str(row.get(label, '')))}</p>" for row in rows]
    return "".join(items)


def _render_simple_table(rows: list[dict]) -> str:
    if not rows:
        return "<p>No comparable data.</p>"
    headers = list(rows[0].keys())
    thead = "".join(f"<th>{escape(str(header))}</th>" for header in headers)
    body_rows = []
    for row in rows:
        body_rows.append("<tr>" + "".join(f"<td>{escape(_format_value(row.get(header)))}</td>" for header in headers) + "</tr>")
    return f"<table><thead><tr>{thead}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def _render_download_link(output_dir: Path, target: str, label: str) -> str:
    return f'<a href="{_relpath(output_dir, target)}">{escape(label)}</a>'


def _relpath(base: Path, target: str) -> str:
    return Path(target).resolve().relative_to(base.resolve()).as_posix()


def _format_value(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def _mode_kicker(comparison_mode: str) -> str:
    return {
        "generic": "Condition comparison",
        "mutation": "Mutation comparison",
        "sampling": "Sampling comparison",
        "replicate": "Replicate comparison",
    }.get(comparison_mode, "Condition comparison")


def _mode_subtitle(comparison_mode: str) -> str:
    return {
        "generic": "Compare two analyzed Immunex conditions across quality, interface burial, flexibility, and interaction landscape.",
        "mutation": "Compare two analyzed mutation conditions across quality, interface burial, flexibility, and interaction landscape.",
        "sampling": "Compare two analyzed sampling conditions across quality, interface burial, flexibility, and interaction landscape.",
        "replicate": "Compare two analyzed replicate conditions across quality, interface burial, flexibility, and interaction landscape.",
    }.get(comparison_mode, "Compare two analyzed Immunex conditions across quality, interface burial, flexibility, and interaction landscape.")
