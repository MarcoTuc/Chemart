// The Chemart simulation pit: pick a chemistry, run it, watch it, measure it.
// Talks to the local pit server (hub/src/chemart_hub/pit.py); runs stream in
// as newline-delimited JSON and are drawn with uPlot as they come.
(function () {
  "use strict";

  const $ = (id) => document.getElementById(id);
  const HEADERS = { "Content-Type": "application/json", "X-Chemart-Pit": "1" };
  const COLORS = ["#0d9488", "#4f46e5", "#d97706", "#db2777", "#2563eb", "#65a30d", "#9333ea",
                  "#dc2626", "#0891b2", "#ca8a04", "#7c3aed", "#059669", "#888888"];
  const TYPE_LABEL = { given: "given networks", generator: "generators", gas: "Turing gases" };
  const DIST = { constant: ["value"], uniform: ["low", "high"], loguniform: ["low", "high"],
                 lognormal: ["mean", "sigma"], exponential: ["scale"], gamma: ["shape", "scale"] };

  let chemistries = [];
  let info = null;           // describe_chemistry of the current chemistry
  let getGenerate = () => ({});
  let getEvolve = () => ({});
  let getRates = () => null;
  let getX0 = () => null;
  let controller = null;

  // ------------------------------------------------------------------ helpers
  function el(tag, attrs, ...children) {
    const node = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (k === "text") node.textContent = v;
      else if (k === "hidden") node.hidden = v;
      else if (k === "style") node.style.cssText = v;     // CSSOM, allowed by the pit's CSP
      else node.setAttribute(k, v);
    }
    for (const child of children) if (child != null) node.append(child);
    return node;
  }

  function short(text, n) {
    text = String(text);
    return text.length > n ? text.slice(0, n - 1) + "…" : text;
  }

  function fmt(v) {
    if (v === null || v === undefined) return "—";
    if (typeof v === "number") return Number.isInteger(v) ? String(v) : v.toPrecision(4);
    if (typeof v === "boolean") return v ? "yes" : "no";
    if (typeof v === "object") return Object.entries(v).map(([k, x]) => `${k}: ${fmt(x)}`).join(", ");
    return String(v);
  }

  function status(text, error) {
    const s = $("status");
    s.textContent = text;
    s.className = error ? "pit-error" : "muted";
  }

  async function getJSON(url) {
    const r = await fetch(url);
    if (!r.ok) throw new Error((await r.json()).detail || r.statusText);
    return r.json();
  }

  // POST and read newline-delimited JSON, calling onMessage for each line.
  async function stream(url, body, onMessage) {
    controller = new AbortController();
    $("stop").disabled = false;
    try {
      const r = await fetch(url, { method: "POST", headers: HEADERS, body: JSON.stringify(body),
                                   signal: controller.signal });
      const reader = r.body.getReader();
      const decoder = new TextDecoder();
      let buffer = "";
      for (;;) {
        const { value, done } = await reader.read();
        if (done) break;
        buffer += decoder.decode(value, { stream: true });
        let nl;
        while ((nl = buffer.indexOf("\n")) >= 0) {
          const line = buffer.slice(0, nl);
          buffer = buffer.slice(nl + 1);
          if (line.trim()) onMessage(JSON.parse(line));
        }
      }
    } catch (err) {
      if (err.name === "AbortError") status("Stopped.");
      else status(String(err), true);
    } finally {
      $("stop").disabled = true;
      controller = null;
    }
  }

  // ------------------------------------------------------------ parameter forms
  // One input per parameter, from the JSON Schema of describe_chemistry.
  function paramForm(container, schema) {
    container.replaceChildren();
    const inputs = {};
    for (const [name, s] of Object.entries(schema.properties || {})) {
      let input;
      if (s.enum) {
        input = el("select", {});
        for (const v of s.enum) input.append(el("option", { value: JSON.stringify(v), text: String(v) }));
        input.value = JSON.stringify(s.default);
      } else if (s.type === "boolean") {
        input = el("input", { type: "checkbox" });
        input.checked = !!s.default;
      } else if (s.type === "integer" || s.type === "number") {
        input = el("input", { type: "number", step: s.type === "integer" ? "1" : "any" });
        if (s.minimum !== undefined) input.min = s.minimum;
        if (s.maximum !== undefined) input.max = s.maximum;
        input.value = s.default ?? "";
      } else if (s.type === "array" || s.type === "object") {
        input = el("textarea", { rows: "2", spellcheck: "false" });
        input.value = JSON.stringify(s.default ?? (s.type === "array" ? [] : {}));
      } else {
        input = el("input", { type: "text" });
        input.value = s.default ?? "";
      }
      inputs[name] = { input, schema: s };
      container.append(el("label", {}, name, input, el("small", { text: short(s.description || "", 220) })));
    }
    return () => {
      const out = {};
      for (const [name, { input, schema: s }] of Object.entries(inputs)) {
        let v;
        if (s.enum) v = JSON.parse(input.value);
        else if (s.type === "boolean") v = input.checked;
        else if (s.type === "integer") v = parseInt(input.value, 10);
        else if (s.type === "number") v = parseFloat(input.value);
        else if (s.type === "array" || s.type === "object") v = JSON.parse(input.value || "null");
        else v = input.value;
        if (JSON.stringify(v) !== JSON.stringify(s.default)) out[name] = v;
      }
      return out;
    };
  }

  // A rate or initial-state editor: the chemistry's own, a constant, a
  // distribution, a table (JSON) or a file (.json or two-column .csv).
  function specEditor(container, title, keyHint) {
    container.replaceChildren();
    const mode = el("select", {});
    for (const [v, t] of [["own", "the chemistry's own"], ["constant", "one value for all"],
                          ["dist", "a distribution"], ["table", "a table (JSON)"], ["file", "a file"]]) {
      mode.append(el("option", { value: v, text: t }));
    }
    const body = el("div", {});
    const box = el("fieldset", { class: "pit-spec" }, el("legend", { text: title }), mode, body);
    container.append(box);
    let fileSpec = null;
    let get = () => null;
    mode.addEventListener("change", () => {
      body.replaceChildren();
      fileSpec = null;
      if (mode.value === "constant") {
        const v = el("input", { type: "number", step: "any", value: "1" });
        body.append(v);
        get = () => parseFloat(v.value);
      } else if (mode.value === "dist") {
        const d = el("select", {});
        for (const name of Object.keys(DIST)) d.append(el("option", { value: name, text: name }));
        const fields = el("div", { class: "pit-row" });
        const draw = () => {
          fields.replaceChildren();
          for (const p of DIST[d.value]) {
            fields.append(el("label", {}, p, el("input", { type: "number", step: "any", value: p === "high" ? "10" : p === "sigma" || p === "scale" || p === "shape" || p === "value" || p === "low" ? "1" : "0", "data-p": p })));
          }
        };
        d.addEventListener("change", draw);
        draw();
        body.append(d, fields);
        get = () => {
          const spec = { dist: d.value };
          for (const input of fields.querySelectorAll("input")) spec[input.dataset.p] = parseFloat(input.value);
          return spec;
        };
      } else if (mode.value === "table") {
        const t = el("textarea", { rows: "3", spellcheck: "false", placeholder: keyHint });
        body.append(t);
        get = () => (t.value.trim() ? JSON.parse(t.value) : null);
      } else if (mode.value === "file") {
        const f = el("input", { type: "file", accept: ".json,.csv" });
        const note = el("small", { class: "muted", text: "read in the browser and sent as a table" });
        f.addEventListener("change", async () => {
          const file = f.files[0];
          if (!file) return;
          const text = await file.text();
          if (file.name.endsWith(".json")) fileSpec = JSON.parse(text);
          else {
            fileSpec = {};
            for (const line of text.split(/\r?\n/)) {
              const i = line.lastIndexOf(",");
              if (i < 0) continue;
              const key = line.slice(0, i).trim(), value = parseFloat(line.slice(i + 1));
              if (key && !Number.isNaN(value)) fileSpec[key] = value;
            }
          }
          note.textContent = `${Object.keys(fileSpec).length} entries from ${file.name}`;
        });
        body.append(f, note);
        get = () => fileSpec;
      } else {
        get = () => null;
      }
    });
    return () => get();
  }

  // ------------------------------------------------------------------ charts
  function chart(title, labels, wide) {
    const box = el("div", { class: "pit-chart" + (wide ? " wide" : "") }, el("h4", { text: title }));
    $("charts").append(box);
    const width = Math.max(280, box.clientWidth - 20);
    const series = [{ label: "t" }].concat(labels.map((l, i) => ({
      label: short(l, 24), stroke: COLORS[i % COLORS.length], width: 1.5,
    })));
    const data = [[]].concat(labels.map(() => []));
    const plot = new uPlot({ width, height: wide ? 300 : 180, series, scales: { x: { time: false } },
                             legend: { show: labels.length > 1 } }, data, box);
    return { plot, data, labels };
  }

  // Frames can arrive faster than the screen refreshes: redraw each chart at
  // most once per animation frame.
  const dirty = new Set();
  let scheduled = false;
  const PLOT_POINTS = 3000;      // a run without end would otherwise fill the page's memory

  function push(c, t, values) {
    c.data[0].push(t);
    c.labels.forEach((l, i) => c.data[i + 1].push(values[l] ?? null));
    if (c.data[0].length > PLOT_POINTS) for (const row of c.data) row.shift();
    dirty.add(c);
    if (!scheduled) {
      scheduled = true;
      requestAnimationFrame(() => {
        scheduled = false;
        for (const d of dirty) d.plot.setData(d.data);
        dirty.clear();
      });
    }
  }

  function measuresTable(result, title) {
    const rows = Object.entries(result.values).map(([k, v]) =>
      el("tr", {}, el("td", { class: "mono", text: k }), el("td", { text: fmt(v) })));
    const why = Object.entries(result.reasons || {}).map(([k, v]) =>
      el("tr", {}, el("td", { class: "mono faint", text: k }), el("td", { class: "muted", text: v })));
    $("measures").replaceChildren(el("h4", { text: title }),
      el("table", {}, el("tbody", {}, ...rows, ...why)));
  }

  function abundance(state) {
    const total = Object.values(state).reduce((a, b) => a + b, 0) || 1;
    const top = Object.entries(state).sort((a, b) => b[1] - a[1]);
    const max = top.length ? top[0][1] : 1;
    $("abundance").replaceChildren(el("h4", { text: "Most abundant now" }), el("table", {}, el("tbody", {},
      ...top.map(([s, v]) => el("tr", {},
        el("td", { class: "pit-id", title: s, text: s }),
        el("td", { class: "num", text: fmt(v) }),
        el("td", { style: "width:40%" }, el("div", { class: "pit-bar", style: `width:${(100 * v / max).toFixed(1)}%` })),
        el("td", { class: "num muted", text: (100 * v / total).toFixed(1) + "%" }))))));
  }

  function reset() {
    $("charts").replaceChildren();
    $("abundance").replaceChildren();
    $("measures").replaceChildren();
    $("download").replaceChildren();
  }

  function downloadLink(run) {
    $("download").replaceChildren(el("a", { class: "btn btn-sm", href: `/api/runs/${run}.json`,
                                            text: "Download the trajectory (JSON)" }));
  }

  // The length parameter of the process, and what ticking "run until I stop" does to it.
  function endlessNote() {
    const box = $("endless");
    box.disabled = !info.duration;
    $("endless-note").textContent = info.duration
      ? `sets ${info.duration} = 0; Stop ends it and keeps what ran`
      : "this chemistry has no length parameter to zero";
  }

  // ------------------------------------------------------------------- runs
  function seed() {
    const v = $("seed").value;
    return v === "" ? null : parseInt(v, 10);
  }

  async function simulate() {
    reset();
    status("Generating the network…");
    let species = null;
    let tracked = null;
    const body = {
      chemistry: info.id, seed: seed(), params: getGenerate(), method: $("sim-method").value,
      t_end: parseFloat($("t-end").value), points: parseInt($("points").value, 10),
      volume: parseFloat($("volume").value), rates: getRates(), x0: getX0(),
      fill_only: $("fill-only").checked,
    };
    await stream("/api/run", body, (m) => {
      if (m.type === "network") {
        measuresTable(m.measures, `Measures of the network (${m.summary.species} species, ${m.summary.reactions} reactions)`);
        status(body.method === "ode" ? "Integrating…" : "Sampling a stochastic path…");
      } else if (m.type === "frame") {
        if (!species) {
          species = chart(body.method === "ode" ? "Amounts (rate equations)" : "Amounts (one stochastic path)",
                          Object.keys(m.state), true);
          tracked = chart("Population", Object.keys(m.measures), false);
        }
        push(species, m.t, m.state);
        push(tracked, m.t, m.measures);
      } else if (m.type === "done") {
        status(`Done: ${m.frames} frames.`);
        downloadLink(m.run);
      } else if (m.type === "error") {
        status(m.message, true);
      }
    });
  }

  async function evolve() {
    reset();
    status("Evolving…");
    const charts = {};
    const tracked = Array.from(document.querySelectorAll("#tracked .chip.on")).map((c) => c.dataset.name);
    const params = getEvolve();
    const endless = $("endless").checked && info.duration;
    if (endless) params[info.duration] = 0;
    const body = { chemistry: info.id, seed: seed(), params, method: "evolve",
                   every: parseInt($("every").value, 10), window: parseInt($("window").value, 10),
                   measures: tracked };
    let frames = 0;
    await stream("/api/run", body, (m) => {
      if (m.type === "frame") {
        frames += 1;
        const values = Object.assign({}, m.measures);
        for (const [k, v] of Object.entries(m.observables || {})) if (typeof v === "number") values[k] = v;
        for (const [k, v] of Object.entries(values)) {
          if (!charts[k]) charts[k] = chart(k, [k], false);
          push(charts[k], m.t, { [k]: v });
        }
        if (frames % 2 === 1) abundance(m.state);
        status(`Evolving… ${info.clock || "t"} = ${fmt(m.t)}${endless ? " (press Stop when you have seen enough)" : ""}`);
      } else if (m.type === "start") {
        downloadLink(m.run);
      } else if (m.type === "network") {
        measuresTable(m.measures, `Measures of the evolved network (${m.summary.species} species, ${m.summary.reactions} reactions)`);
      } else if (m.type === "done") {
        status(`Done: ${m.frames} frames. Turnover ${fmt(m.trajectory_measures.turnover)}, novelty rate ${fmt(m.trajectory_measures.novelty_rate)}.`);
        downloadLink(m.run);
      } else if (m.type === "error") {
        status(m.message, true);
      }
    });
  }

  async function measureNetwork() {
    reset();
    status("Measuring…");
    const r = await fetch("/api/measure", { method: "POST", headers: HEADERS, body: JSON.stringify({
      chemistry: info.id, seed: seed(), params: getGenerate(), cost: "moderate" }) });
    const result = await r.json();
    if (result.error) return status(result.error, true);
    measuresTable(result, "Measures of the network (cheap and moderate)");
    status("Done.");
  }

  async function sweep() {
    reset();
    const param = $("sweep-param").value;
    const values = $("sweep-values").value.split(",").map((v) => v.trim()).filter(Boolean)
      .map((v) => (Number.isNaN(Number(v)) ? v : Number(v)));
    const seeds = Array.from({ length: parseInt($("sweep-seeds").value, 10) }, (_, i) => i);
    const y = $("sweep-measure").value;
    const rows = [];
    const means = {};
    const c = chart(`${y} against ${param} (mean over ${seeds.length} seeds)`, [y], true);
    status("Sweeping…");
    await stream("/api/sweep", { chemistry: info.id, param, values, seeds, params: getGenerate() }, (m) => {
      if (m.type === "row") {
        rows.push(m.row);
        const key = m.row[param];
        (means[key] = means[key] || []).push(m.row[y]);
        if (means[key].length === seeds.length && typeof key === "number") {
          const avg = means[key].reduce((a, b) => a + b, 0) / seeds.length;
          push(c, key, { [y]: avg });
        }
        status(`Swept ${rows.length} of ${values.length * seeds.length} networks…`);
      } else if (m.type === "done") {
        status(`Done: ${rows.length} networks.`);
        const cols = Object.keys(rows[0] || {});
        $("measures").replaceChildren(el("table", {},
          el("thead", {}, el("tr", {}, ...cols.map((k) => el("th", { class: "mono", text: k })))),
          el("tbody", {}, ...rows.map((r) => el("tr", {}, ...cols.map((k) => el("td", { text: fmt(r[k]) })))))));
      } else if (m.type === "error") {
        status(m.message, true);
      }
    });
  }

  // --------------------------------------------------------------- choosing
  async function choose(id) {
    reset();
    status("Loading…");
    info = await getJSON(`/api/chemistry/${encodeURIComponent(id)}`);
    const row = chemistries.find((c) => c.id === id);
    const faces = info.faces || [];
    $("about").textContent = `${TYPE_LABEL[info.type] || info.type} · ${info.family}` +
      (info.clock ? ` · clock: ${info.clock}` : "") + ` · ${info.intuition ? short(info.intuition, 260) : ""}`;

    $("faces").replaceChildren(...faces.map((f) => el("span", { class: "chip on",
      text: f === "generate" ? "network: generate_network" : "process: chemart.evolve" })));
    $("network-panel").hidden = !faces.includes("generate");
    $("process-panel").hidden = !faces.includes("evolve");

    if (faces.includes("generate")) {
      getGenerate = paramForm($("generate-params"), info.params);
      getRates = specEditor($("rates-editor"), "Rates", '{"A + B -> C": 2.0, "*": 1.0}');
      getX0 = specEditor($("x0-editor"), "Initial amounts", '{"A": 1.0, "*": 0.0}');
      const numeric = Object.entries(info.params.properties || {})
        .filter(([, s]) => s.type === "integer" || s.type === "number");
      $("sweep").hidden = info.type !== "generator" || !numeric.length;
      $("sweep-param").replaceChildren(...numeric.map(([n]) => el("option", { value: n, text: n })));
      const scalar = info.measures.filter((m) => m.input === "network" && m.cost === "cheap");
      $("sweep-measure").replaceChildren(...scalar.map((m) => el("option", { value: m.name, text: m.name })));
    } else {
      getGenerate = () => ({});
    }
    if (faces.includes("evolve")) {
      getEvolve = paramForm($("evolve-params"), info.evolve_params || info.params);
      endlessNote();
      const choices = info.measures.filter((m) => m.input === "state" ||
        (m.input === "network" && m.cost === "cheap"));
      const defaults = new Set(["richness", "shannon", "dominance", "n_reactions"]);
      $("tracked").replaceChildren(...choices.map((m) => {
        const chip = el("a", { class: "chip" + (defaults.has(m.name) ? " on" : ""), href: "#",
                                title: m.meaning, text: m.name, "data-name": m.name });
        chip.addEventListener("click", (e) => { e.preventDefault(); chip.classList.toggle("on"); });
        return chip;
      }));
    } else {
      getEvolve = () => ({});
    }
    status(row ? `${row.name}: ready.` : "Ready.");
  }

  async function start() {
    chemistries = await getJSON("/api/chemistries");
    const select = $("chemistry");
    for (const type of ["given", "generator", "gas"]) {
      const group = el("optgroup", { label: TYPE_LABEL[type] });
      for (const c of chemistries.filter((c) => c.type === type).sort((a, b) => a.name.localeCompare(b.name))) {
        group.append(el("option", { value: c.id, text: c.name }));
      }
      select.append(group);
    }
    select.addEventListener("change", () => choose(select.value));
    $("sim-method").addEventListener("change", () => { $("volume-label").hidden = $("sim-method").value !== "ssa"; });
    $("run-sim").addEventListener("click", simulate);
    $("run-measure").addEventListener("click", measureNetwork);
    $("run-sweep").addEventListener("click", sweep);
    $("run-evolve").addEventListener("click", evolve);
    $("stop").addEventListener("click", () => controller && controller.abort());
    const first = new URLSearchParams(location.search).get("chemistry") || "brusselator";
    select.value = first;
    await choose(select.value);
  }

  start().catch((err) => status(String(err), true));
})();
