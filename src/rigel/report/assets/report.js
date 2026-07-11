"use strict";
/* Rigel QC report front-end. Consumes window.__RIGEL__ = {model, fl_specs}.
   Native HTML/SVG for tiles, bars, gauge, tables; Vega-Lite (window.vegaEmbed,
   inlined at build time) for the fragment-length distribution charts. */
(function () {
  const R = window.__RIGEL__ || {};
  const M = R.model || {};
  const FL_SPECS = R.fl_specs || {};

  /* ---------- helpers ---------- */
  const grp = (n) => Math.round(n).toLocaleString("en-US");
  const si = (n) => n >= 1e6 ? (n / 1e6).toFixed(1) + "M" : n >= 1e3 ? (n / 1e3).toFixed(1) + "k" : String(Math.round(n));
  const pc = (x) => (x * 100).toFixed(1) + "%";
  const SVGNS = "http://www.w3.org/2000/svg";
  const $ = (id) => document.getElementById(id);
  function H(tag, cls, html) { const e = document.createElement(tag); if (cls) e.className = cls; if (html != null) e.innerHTML = html; return e; }
  function el(tag, attrs, kids) { const e = document.createElementNS(SVGNS, tag); for (const k in (attrs || {})) e.setAttribute(k, attrs[k]); (kids || []).forEach((c) => e.appendChild(c)); return e; }
  function txt(t) { const e = document.createElementNS(SVGNS, "tspan"); e.textContent = t; return e; }
  function svg(w, h) { return el("svg", { viewBox: `0 0 ${w} ${h}`, role: "img" }); }
  function slotColor(cls) { return cls === "fmuted" ? "var(--muted)" : `var(--s${String(cls).replace(/\D/g, "") || "1"})`; }

  /* ---------- components ---------- */
  function icon(name) {
    const m = {
      target: '<circle cx="8" cy="8" r="6.5" fill="none" stroke="currentColor" stroke-width="1.4"/><circle cx="8" cy="8" r="3" fill="none" stroke="currentColor" stroke-width="1.4"/><circle cx="8" cy="8" r="1" fill="currentColor"/>',
      arrow: '<path d="M2 8h9M8 4l4 4-4 4" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>',
      warn: '<path d="M8 1.5l6.5 12H1.5z" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linejoin="round"/><path d="M8 6v3.5" stroke="currentColor" stroke-width="1.4" stroke-linecap="round"/><circle cx="8" cy="11.6" r=".9" fill="currentColor"/>',
      check: '<path d="M2.5 8.5l3.5 3.5 7.5-8" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"/>',
    };
    return `<svg class="ic" viewBox="0 0 16 16" style="color:var(--vc,var(--accent))">${m[name] || ""}</svg>`;
  }
  function kpis(node, items) {
    if (!node) return; node.innerHTML = "";
    (items || []).forEach((d) => {
      const k = H("div", "kpi");
      k.appendChild(H("div", "l", d.l));
      k.appendChild(H("div", "v", `${d.v}${d.u ? `<small> ${d.u}</small>` : ""}`));
      node.appendChild(k);
    });
  }
  function stackBar(node, items) {
    if (!node) return;
    node.innerHTML = "";
    const tot = (items || []).reduce((s, d) => s + d.value, 0);
    if (tot <= 0) {
      // Empty state — a full-width neutral track, not a row of 2px slivers.
      const seg = H("div", "seg"); seg.style.flex = "1"; seg.style.background = "var(--surface-3)";
      const l = H("div", "pct", "none detected"); l.style.color = "var(--muted)"; l.style.textShadow = "none";
      seg.appendChild(l); node.appendChild(seg); return;
    }
    items.forEach((d) => {
      const seg = H("div", "seg"); seg.style.flex = `${d.value} 0 0`;
      seg.style.background = slotColor(d.cls);
      const frac = d.value / tot;
      if (frac > 0.06) seg.appendChild(H("div", "pct", pc(frac)));
      seg.title = `${d.label}: ${grp(d.value)} (${pc(frac)})`;
      node.appendChild(seg);
    });
  }
  function legend(node, items, withCounts) {
    if (!node) return;
    node.innerHTML = "";
    items.forEach((d) => {
      const li = H("span", "li");
      const sw = H("span", "swatch"); sw.style.background = slotColor(d.cls);
      li.appendChild(sw); li.appendChild(H("span", null, d.label));
      if (withCounts) li.appendChild(H("span", "num", ` ${grp(d.value)}`));
      node.appendChild(li);
    });
  }
  function gauge(node, val) {
    if (!node) return;
    const W = 300, Hh = 168, cx = 150, cy = 150, r = 118, s = svg(W, Hh);
    const a0 = Math.PI, pt = (a, rr) => [cx + rr * Math.cos(a), cy - rr * Math.sin(a)];
    const arc = (b, c, col, w) => {
      const large = Math.abs(b - c) > Math.PI ? 1 : 0, sweep = c < b ? 1 : 0;
      const [x0, y0] = pt(b, r), [x1, y1] = pt(c, r);
      return el("path", { d: `M ${x0} ${y0} A ${r} ${r} 0 ${large} ${sweep} ${x1} ${y1}`, fill: "none", stroke: col, "stroke-width": w, "stroke-linecap": "round" });
    };
    s.appendChild(arc(a0, 0, "var(--grid)", 16));
    for (let v = 0.5; v <= 1.001; v += 0.1) {
      const a = a0 - (v - 0.5) / 0.5 * Math.PI;
      const [x1, y1] = pt(a, r + 9), [x2, y2] = pt(a, r + 2);
      s.appendChild(el("line", { x1, y1, x2, y2, class: "ax" }));
      const [tx, ty] = pt(a, r + 20);
      s.appendChild(el("text", { class: "tk", x: tx, y: ty + 3, "text-anchor": "middle" }, [txt(v.toFixed(1))]));
    }
    const av = a0 - (Math.min(1, Math.max(0.5, val)) - 0.5) / 0.5 * Math.PI;
    s.appendChild(arc(a0, av, "var(--s1)", 16));
    const [dx, dy] = pt(av, r); s.appendChild(el("circle", { cx: dx, cy: dy, r: 6, fill: "var(--s1)", stroke: "var(--surface)", "stroke-width": 2 }));
    s.appendChild(el("text", { x: cx, y: cy - 26, "text-anchor": "middle", class: "dl", style: "font-size:34px" }, [txt(val.toFixed(3))]));
    s.appendChild(el("text", { x: cx, y: cy - 8, "text-anchor": "middle", class: "tklab" }, [txt("strand specificity")]));
    node.innerHTML = ""; node.appendChild(s);
  }
  function orient(node, read1sense) {
    if (!node) return;
    const W = 440, Hh = 116, y = 62, s = svg(W, Hh);
    s.appendChild(el("line", { x1: 30, y1: y, x2: 410, y2: y, class: "ax", "stroke-width": 2 }));
    s.appendChild(el("path", { d: `M 410 ${y} l -10 -5 l 0 10 z`, fill: "var(--baseline)" }));
    s.appendChild(el("text", { x: 30, y: y - 14, class: "tklab" }, [txt("transcript  5′ → 3′  (sense strand)")]));
    // Reads share a fixed [x1,x2] span; the arrowhead end encodes sense/antisense
    // so a leftward (antisense) arrow never runs off the canvas.
    const x1 = 120, x2 = 300, r1s = read1sense, r2s = !read1sense;
    readArrow(s, x1, x2, y - 34, "s1", `R1  ${r1s ? "sense →" : "← antisense"}`, r1s);
    readArrow(s, x1, x2, y + 24, "s2", `R2  ${r2s ? "sense →" : "← antisense"}`, r2s);
    node.innerHTML = ""; node.appendChild(s);
  }
  function readArrow(s, x1, x2, y, color, label, sense) {
    s.appendChild(el("line", { x1, y1: y, x2, y2: y, stroke: `var(--${color})`, "stroke-width": 3, "stroke-linecap": "round" }));
    const head = sense ? `M ${x2} ${y} l -9 -5 l 0 10 z` : `M ${x1} ${y} l 9 -5 l 0 10 z`;
    s.appendChild(el("path", { d: head, fill: `var(--${color})` }));
    s.appendChild(el("text", { x: x1, y: y - 9, class: "tklab", style: "font-family:var(--mono);font-size:11px" }, [txt(label)]));
  }
  function donut(node, items, centerFrac) {
    if (!node) return;
    const W = 210, Hh = 210, cx = 105, cy = 105, r = 88, rin = 52, s = svg(W, Hh);
    const tot = items.reduce((a, d) => a + d.value, 0) || 1; let a = -Math.PI / 2;
    items.forEach((d) => {
      const frac = d.value / tot; if (frac <= 0) return;
      const a2 = a + frac * 2 * Math.PI, large = frac > 0.5 ? 1 : 0;
      const x0 = cx + r * Math.cos(a), y0 = cy + r * Math.sin(a), x1 = cx + r * Math.cos(a2), y1 = cy + r * Math.sin(a2);
      const xi0 = cx + rin * Math.cos(a2), yi0 = cy + rin * Math.sin(a2), xi1 = cx + rin * Math.cos(a), yi1 = cy + rin * Math.sin(a);
      s.appendChild(el("path", { d: `M ${x0} ${y0} A ${r} ${r} 0 ${large} 1 ${x1} ${y1} L ${xi0} ${yi0} A ${rin} ${rin} 0 ${large} 0 ${xi1} ${yi1} Z`, fill: slotColor(d.cls), stroke: "var(--surface)", "stroke-width": 2 }));
      a = a2;
    });
    s.appendChild(el("text", { x: cx, y: cy - 4, "text-anchor": "middle", class: "dl", style: "font-size:24px" }, [txt(pc(centerFrac))]));
    s.appendChild(el("text", { x: cx, y: cy + 15, "text-anchor": "middle", class: "tklab" }, [txt("RNA")]));
    node.innerHTML = ""; node.appendChild(s);
  }

  /* ---------- tables ---------- */
  function fillTable(node, head, rows) {
    if (!node) return;
    node.innerHTML = `<thead><tr>${head}</tr></thead>`;
    const tb = H("tbody"); rows.forEach((r) => tb.appendChild(r)); node.appendChild(tb);
  }
  function alignTable() {
    const rows = (M.alignment.table || []).map(([k, v, f]) => {
      const tr = H("tr");
      tr.appendChild(H("td", null, k));
      tr.appendChild(H("td", "n num", grp(v)));
      const c = H("td", "n barcell");
      c.innerHTML = `<div class="fill" style="width:${(f * 100).toFixed(0)}%"></div><span class="num">${pc(f)}</span>`;
      tr.appendChild(c); return tr;
    });
    fillTable($("align-table"), `<th>Metric</th><th class="n">Reads</th><th class="n">%</th>`, rows);
  }
  function flTable() {
    const rows = (M.fl.table || []).map((r) => {
      const tr = H("tr");
      const color = r[0] === "rna" ? "var(--s1)" : r[0] === "gdna" ? "var(--s8)" : "var(--ink)";
      tr.appendChild(H("td", null, `<span class="num" style="color:${color}">${r[0]}</span>`));
      tr.appendChild(H("td", "n num", si(r[1])));
      [r[2], r[3], r[4], r[5]].forEach((v) => tr.appendChild(H("td", "n num", v)));
      tr.appendChild(H("td", "n num", r[6] + "%"));
      return tr;
    });
    fillTable($("fl-table"), `<th>Category</th><th class="n">N</th><th class="n">mean</th><th class="n">sd</th><th class="n">med</th><th class="n">mode</th><th class="n">ovf</th>`, rows);
  }
  function poolTable() {
    const rows = (M.quant.table || []).map(([k, v, f, d]) => {
      const tr = H("tr");
      tr.appendChild(H("td", null, `<b>${k}</b>`));
      tr.appendChild(H("td", "n num", grp(v)));
      tr.appendChild(H("td", "n num", pc(f)));
      tr.appendChild(H("td", null, `<span style="color:var(--ink-2);font-size:12px">${d}</span>`));
      return tr;
    });
    fillTable($("pool-table"), `<th>Pool</th><th class="n">Fragments</th><th class="n">Share</th><th>Origin</th>`, rows);
  }
  function geneTable(filter) {
    const q = (filter || "").toLowerCase();
    const all = M.genes.rows || [];
    const rows = [];
    let shown = 0;
    for (const g of all) {
      if (q && !(String(g[0]).toLowerCase().includes(q) || String(g[1]).toLowerCase().includes(q))) continue;
      shown++;
      const tr = H("tr");
      tr.appendChild(H("td", null, `<b>${g[0]}</b>`));
      tr.appendChild(H("td", null, `<span class="num" style="color:var(--muted);font-size:11.5px">${g[1]}</span>`));
      tr.appendChild(H("td", "n num", grp(g[2])));
      tr.appendChild(H("td", "n num", grp(g[3])));
      tr.appendChild(H("td", "n num", g[4].toLocaleString("en-US")));
      tr.appendChild(H("td", "n num", g[5]));
      rows.push(tr);
    }
    fillTable($("gene-table"), `<th>Gene</th><th>ID</th><th class="n">count</th><th class="n">spliced</th><th class="n">TPM</th><th class="n">n_tx</th>`, rows);
    const c = $("gcount");
    if (c) {
      const g = M.genes;
      const tail = g.truncated ? ` · top ${grp(g.shown)} by TPM (full set in gene_quant.feather)` : "";
      c.textContent = `${grp(shown)} shown · ${grp(g.total)} expressed genes${tail}`;
    }
  }
  function config() {
    const b = $("cfg-body"); if (!b) return; b.innerHTML = "";
    const cfg = M.config || {};
    for (const g in cfg) {
      b.appendChild(H("div", "cfg-group", g));
      for (const k in cfg[g]) {
        const r = H("div", "cfg-row");
        r.appendChild(H("span", "kk", k));
        r.appendChild(H("span", "vv", cfg[g][k]));
        b.appendChild(r);
      }
    }
  }

  /* ---------- hero / verdicts ---------- */
  function hero() {
    const m = M.meta || {};
    if ($("hero-title")) $("hero-title").textContent = m.sample || "Rigel QC report";
    const meta = $("hero-meta");
    if (meta) {
      const bits = [];
      if (m.sample) bits.push(`Sample <b>${m.sample}</b>`);
      if (m.created) bits.push(`Generated <b>${m.created}</b>`);
      if (m.version) bits.push(`Rigel <b>${m.version}</b>`);
      if (m.index) bits.push(`Index <b>${m.index.split("/").pop()}</b>`);
      meta.innerHTML = bits.join("");
    }
    const v = $("verdicts"); if (v) {
      v.innerHTML = "";
      (M.verdicts || []).forEach((d) => {
        const t = H("div", `vt s-${d.s}`);
        t.innerHTML = `<div class="k">${icon(d.icon)}${d.k}</div><div class="v">${d.v}</div><div class="n">${d.n}</div>`;
        v.appendChild(t);
      });
    }
    const wb = $("warn-banner");
    if (wb && m.warnings && m.warnings.length) {
      wb.innerHTML = `<b>Note.</b> ${m.warnings.join(" ")}`;
      wb.style.display = "";
    }
  }

  /* ---------- Vega-Lite (fragment length) ---------- */
  const PAL = {
    light: { cat: ["#2a78d6", "#1baf7a", "#eda100", "#008300", "#4a3aa7", "#e34948", "#e87ba4", "#eb6834"], ink: "#0d1117", ink2: "#48505e", muted: "#8b93a1", grid: "#e6eaf0", baseline: "#c7cdd8" },
    dark: { cat: ["#3987e5", "#199e70", "#c98500", "#008300", "#9085e9", "#e66767", "#d55181", "#d95926"], ink: "#eef1f6", ink2: "#aab2c0", muted: "#6c7482", grid: "#222834", baseline: "#333b48" },
  };
  function curTheme() {
    const a = document.documentElement.getAttribute("data-theme");
    if (a === "dark" || a === "light") return a;
    return window.matchMedia && matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  }
  function themeConfig(t) {
    const p = PAL[t] || PAL.light;
    const font = "system-ui,-apple-system,'Segoe UI',Roboto,sans-serif";
    return {
      background: "transparent", font,
      view: { stroke: null },
      range: { category: p.cat },
      mark: { color: p.cat[0] },
      axis: { labelColor: p.muted, titleColor: p.ink2, gridColor: p.grid, domainColor: p.baseline, tickColor: p.baseline, labelFont: font, titleFont: font, gridOpacity: 0.7, labelFontSize: 10, titleFontSize: 11 },
      legend: { labelColor: p.ink2, titleColor: p.ink2, labelFont: font, titleFont: font },
      header: { labelColor: p.ink2, titleColor: p.ink2, labelFont: font, titleFont: font },
      title: { color: p.ink, font },
    };
  }
  function embedAll() {
    if (typeof window.vegaEmbed !== "function") return;
    const t = curTheme();
    const opts = { config: themeConfig(t), renderer: "svg", actions: { export: true, source: true, compiled: false, editor: false } };
    for (const [id, spec] of Object.entries(FL_SPECS)) {
      const node = $("vega-" + id);
      if (!node) continue;
      node.innerHTML = "";
      window.vegaEmbed(node, spec, opts).catch((e) => { node.innerHTML = `<div class="chart-note">chart error: ${e}</div>`; });
    }
  }
  function initFL() {
    const hasVega = typeof window.vegaEmbed === "function";
    const hasSpecs = Object.keys(FL_SPECS).length > 0;
    if (!hasSpecs || !hasVega) {
      ["vega-overlay", "vega-small_multiples"].forEach((id) => {
        const n = $(id); if (n) n.innerHTML = `<div class="chart-note">Fragment-length histograms unavailable${hasVega ? " (no fragment_lengths.feather)" : " (Vega runtime not embedded)"}.</div>`;
      });
      return;
    }
    embedAll();
    // Re-embed on theme change (host toggle stamps data-theme; OS pref may change)
    new MutationObserver(embedAll).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
    if (window.matchMedia) {
      const mq = matchMedia("(prefers-color-scheme: dark)");
      (mq.addEventListener ? mq.addEventListener.bind(mq, "change") : mq.addListener.bind(mq))(embedAll);
    }
  }

  /* ---------- theme toggle (standalone file: no host toggle) ---------- */
  function initThemeToggle() {
    const btn = $("theme-toggle"); if (!btn) return;
    btn.addEventListener("click", () => {
      const next = curTheme() === "dark" ? "light" : "dark";
      document.documentElement.setAttribute("data-theme", next);
    });
  }

  /* ---------- scroll spy ---------- */
  function spy() {
    const links = [...document.querySelectorAll(".rail a")];
    const map = links.map((a) => ({ a, el: document.querySelector(a.getAttribute("href")) })).filter((x) => x.el);
    const obs = new IntersectionObserver((ents) => {
      ents.forEach((e) => { if (e.isIntersecting) { links.forEach((l) => l.classList.remove("active")); const m = map.find((x) => x.el === e.target); if (m) m.a.classList.add("active"); } });
    }, { rootMargin: "-30% 0px -60% 0px" });
    map.forEach((x) => obs.observe(x.el));
  }

  /* ---------- boot ---------- */
  function boot() {
    hero();
    kpis($("align-kpis"), M.alignment.kpis); stackBar($("align-bar"), M.alignment.fate); legend($("align-legend"), M.alignment.fate, true); alignTable();
    kpis($("frag-kpis"), M.fragments.kpis);
    stackBar($("ctx-bar"), M.fragments.ctx); legend($("ctx-legend"), M.fragments.ctx, true);
    stackBar($("chim-bar"), M.fragments.chim); legend($("chim-legend"), M.fragments.chim, true);
    stackBar($("sj-bar"), M.fragments.splice); legend($("sj-legend"), M.fragments.splice, true);
    gauge($("gauge"), M.strand.spec); kpis($("strand-kpis"), M.strand.kpis); orient($("orient"), M.strand.read1_sense);
    if ($("strand-note")) {
      const g = M.strand.contamination_gap, ex = M.strand.exonic_all_spec;
      $("strand-note").innerHTML = `<b>Protocol: ${M.strand.protocol}.</b> ` +
        (ex != null ? `The all-exonic model reads <span class="num">${ex.toFixed(3)}</span> vs the spliced model’s <span class="num">${M.strand.spec.toFixed(3)}</span> — a contamination gap of <span class="num">${g != null ? g.toFixed(4) : "—"}</span> (larger ⇒ more unstranded gDNA drag).` : "");
    }
    flTable();
    stackBar($("pool-bar"), M.quant.pools); legend($("pool-legend"), M.quant.pools, true); poolTable();
    kpis($("pool-kpis"), M.quant.kpis); donut($("pool-donut"), M.quant.pools, M.quant.rna_share);
    geneTable("");
    const gs = $("gsearch"); if (gs) gs.addEventListener("input", (e) => geneTable(e.target.value));
    config();
    initFL(); initThemeToggle(); spy();
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot); else boot();
})();
