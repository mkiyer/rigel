"""Assemble the self-contained report HTML.

Inlines the design-system CSS, the front-end render layer, the view-model
payload, and — when ``vl-convert-python`` is installed — the Vega/Vega-Lite/
Vega-Embed runtime, so the resulting single ``.html`` renders offline with no
CDN, server, or Node dependency.
"""

from __future__ import annotations

import json
from functools import lru_cache
from importlib.resources import files

_ASSETS = files("rigel.report") / "assets"

# ES6 snippet handed to vl_convert.javascript_bundle: pulls the three libraries
# into the bundle and exposes them as globals for report.js to call.
_VEGA_EXPOSE = "window.vegaEmbed = vegaEmbed; window.vega = vega; window.vegaLite = vegaLite;"


def _asset(name: str) -> str:
    return (_ASSETS / name).read_text(encoding="utf-8")


@lru_cache(maxsize=1)
def _vega_bundle() -> str | None:
    """Return the inlinable Vega runtime, or ``None`` if vl-convert is absent.

    The bundle is identical for every report, so it is built once per process and
    cached — this dominates a single build (~0.8 s) and would otherwise repeat for
    every sample in a batch.
    """
    try:
        import vl_convert as vlc
    except ImportError:
        return None
    return vlc.javascript_bundle(_VEGA_EXPOSE)


def _json_for_script(obj) -> str:
    """JSON-encode for safe embedding inside an inline <script> element."""
    return json.dumps(obj, allow_nan=False).replace("</", "<\\/")


_SECTIONS = """
<div class="topbar"><div class="topbar-in">
  <div class="wm"><span class="dot"></span><span class="star">RIGEL</span><span class="sub">Quality&nbsp;Control</span></div>
  <div class="spacer"></div>
  <span class="chip">v__VERSION__</span>
  <span class="chip">schema&nbsp;__SCHEMA__</span>
  <button class="tbtn" id="theme-toggle" type="button" title="Toggle light / dark">◐ theme</button>
</div></div>

<div class="layout">
  <nav class="rail">
    <a href="#hero" class="active"><span class="rn">◆</span>Overview</a>
    <a href="#s-align"><span class="rn">01</span>Alignment</a>
    <a href="#s-frag"><span class="rn">02</span>Fragments</a>
    <a href="#s-strand"><span class="rn">03</span>Strand model</a>
    <a href="#s-fl"><span class="rn">04</span>Fragment length</a>
    <a href="#s-quant"><span class="rn">05</span>Quantification</a>
    <a href="#s-enrich"><span class="rn">06</span>Enrichment</a>
    <a href="#s-density"><span class="rn">07</span>gDNA density</a>
    <a href="#s-genes"><span class="rn">08</span>Gene expression</a>
    <a href="#s-config"><span class="rn">✦</span>Run configuration</a>
  </nav>
  <main class="main">

    <div class="banner" id="warn-banner" style="display:none"></div>

    <section class="hero" id="hero">
      <h1 id="hero-title">Rigel QC report</h1>
      <div class="meta" id="hero-meta"></div>
      <div class="verdicts" id="verdicts"></div>
    </section>

    <section class="panel" id="s-align">
      <div class="ph"><span class="eyebrow">01</span><h2>Alignment statistics</h2><span class="desc">Every read from the BAM, by fate</span></div>
      <div class="pb">
        <div class="kpis" id="align-kpis"></div>
        <div class="grid2 wide-left">
          <div>
            <p class="cap">Read fate — read groups by NH tag</p>
            <div class="stackbar" id="align-bar"></div>
            <div class="legend" id="align-legend"></div>
          </div>
          <div class="tw"><p class="cap">Counts (records vs read groups)</p><table class="data" id="align-table"></table></div>
        </div>
      </div>
    </section>

    <section class="panel" id="s-frag">
      <div class="ph"><span class="eyebrow">02</span><h2>Fragment composition</h2><span class="desc">Overlapping partitions — separate bars, not one pie</span></div>
      <div class="pb">
        <div class="kpis" id="frag-kpis"></div>
        <div class="grid2">
          <div>
            <p class="cap">Genomic context</p>
            <div class="stackbar" id="ctx-bar"></div><div class="legend" id="ctx-legend"></div>
            <p class="cap" style="margin-top:22px">Chimeric breakdown</p>
            <div class="stackbar" id="chim-bar"></div><div class="legend" id="chim-legend"></div>
          </div>
          <div>
            <p class="cap">Splice junctions · by class</p>
            <div class="stackbar" id="sj-bar"></div><div class="legend" id="sj-legend"></div>
          </div>
        </div>
      </div>
    </section>

    <section class="panel" id="s-strand">
      <div class="ph"><span class="eyebrow">03</span><h2>Strand model</h2><span class="desc">Library protocol &amp; strand specificity</span></div>
      <div class="pb"><div class="grid2 wide-left">
        <div>
          <p class="cap">Strand specificity</p>
          <div id="gauge"></div>
          <div class="kpis" id="strand-kpis" style="margin-top:16px"></div>
        </div>
        <div>
          <p class="cap">Read orientation</p>
          <div id="orient"></div>
          <div class="note" id="strand-note" style="margin-top:16px"></div>
        </div>
      </div></div>
    </section>

    <section class="panel" id="s-fl">
      <div class="ph"><span class="eyebrow">04</span><h2>Fragment length distributions</h2><span class="desc">Vega-Lite · hover, zoom, export</span></div>
      <div class="pb">
        <div class="grid2 wide-left">
          <div>
            <p class="cap">gDNA vs RNA — density-normalized overlay</p>
            <div class="vega-chart" id="vega-overlay"></div>
          </div>
          <div class="tw"><p class="cap">Summary statistics</p><table class="data" id="fl-table"></table></div>
        </div>
        <p class="cap" style="margin-top:24px">Per-category histograms</p>
        <div class="tw"><div class="vega-chart" id="vega-small_multiples"></div></div>
      </div>
    </section>

    <section class="panel" id="s-quant">
      <div class="ph"><span class="eyebrow">05</span><h2>Quantification</h2><span class="desc">Mature mRNA · nascent RNA · gDNA</span></div>
      <div class="pb">
        <div class="kpis" id="pool-kpis"></div>
        <div class="grid2 wide-left">
          <div>
            <p class="cap">Assigned fragment pools</p>
            <div class="stackbar" id="pool-bar" style="height:34px"></div>
            <div class="legend" id="pool-legend"></div>
            <div class="tw" style="margin-top:20px"><table class="data" id="pool-table"></table></div>
          </div>
          <div style="display:flex;align-items:center;justify-content:center"><div id="pool-donut"></div></div>
        </div>
      </div>
    </section>

    <section class="panel" id="s-enrich">
      <div class="ph"><span class="eyebrow">06</span><h2>Calibration · Enrichment</h2><span class="desc">Capture on-target enrichment</span></div>
      <div class="pb">
        <div class="kpis" id="enrich-kpis"></div>
        <div class="grid2 wide-left">
          <div>
            <p class="cap">gDNA density — by region count vs by gDNA mass</p>
            <div class="vega-chart" id="vega-capture_kde"></div>
          </div>
          <div><div class="note" id="capture-note"></div></div>
        </div>
      </div>
    </section>

    <section class="panel" id="s-density">
      <div class="ph"><span class="eyebrow">07</span><h2>Calibration · gDNA density across the genome</h2><span class="desc">Per-reference gDNA levels</span></div>
      <div class="pb">
        <div class="kpis" id="density-kpis"></div>
        <p class="cap" id="genome-cap">gDNA density across the genome</p>
        <div class="tw"><div class="vega-chart" id="vega-genome"></div></div>
        <div class="note">Log density, independent y per reference (so a high-density reference such as
          chrM can't flatten the others). Full per-base signal is in
          <span class="num">calibration_track.bedgraph</span> (IGV / UCSC).</div>

        <p class="cap" style="margin-top:22px">All references · by gDNA mass</p>
        <div class="search">
          <input id="rsearch" type="text" placeholder="reference…" aria-label="Search references"/>
          <span class="hint" id="rcount"></span>
        </div>
        <div class="tw"><table class="data" id="ref-table"></table></div>
        <div class="pager" id="ref-pager"></div>
      </div>
    </section>

    <section class="panel" id="s-genes">
      <div class="ph"><span class="eyebrow">08</span><h2>Gene expression</h2><span class="desc">Look up any gene</span></div>
      <div class="pb">
        <div class="search">
          <input id="gsearch" type="text" placeholder="gene name, ID, biotype…" aria-label="Search genes"/>
          <span class="hint" id="gcount"></span>
        </div>
        <div class="tw"><table class="data" id="gene-table"></table></div>
        <div class="pager" id="gene-pager"></div>
      </div>
    </section>

    <section class="panel" id="s-config">
      <div class="ph"><span class="eyebrow">✦</span><h2>Run configuration</h2><span class="desc">Exactly how this run was invoked</span></div>
      <div class="pb"><details class="cfg" open><summary>Resolved parameters <span class="tw-c">rigel quant</span></summary><div class="cfg-body" id="cfg-body"></div></details></div>
    </section>

    <div class="foot">Generated by <b>rigel report</b> · self-contained &amp; offline · built from the run’s <span class="num">summary.json</span> + companion tables.</div>
  </main>
</div>
"""


def render_html(model: dict, charts: dict, title: str) -> str:
    """Render the full self-contained HTML document string."""
    css = _asset("report.css")
    js = _asset("report.js")
    bundle = _vega_bundle()
    meta = model.get("meta", {})

    payload = _json_for_script({"model": model, "charts": charts})
    # Token replace (not str.format) so literal { } in the section HTML/CSS/JS are safe.
    schema = f"report/{meta.get('schema_version')}" if meta.get("schema_version") else "report"
    sections = _SECTIONS.replace("__VERSION__", str(meta.get("version") or "?")).replace(
        "__SCHEMA__", schema
    )

    vega_tag = f"<script>{bundle}</script>" if bundle else ""

    return (
        "<!doctype html>\n"
        '<html lang="en"><head>\n'
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{title}</title>\n"
        f"<style>{css}</style>\n"
        "</head><body>\n"
        f"{sections}\n"
        f'<script id="rigel-data" type="application/json">{payload}</script>\n'
        "<script>window.__RIGEL__ = JSON.parse(document.getElementById('rigel-data').textContent);</script>\n"
        f"{vega_tag}\n"
        f"<script>{js}</script>\n"
        "</body></html>\n"
    )
