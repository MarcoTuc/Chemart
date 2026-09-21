// Chemart Hub, static edition: browse and search in the browser.
// The page is built with every repo on it; this applies the filters, the
// search and the sort from the URL, the way the live hub does on the server.
(function () {
  const grid = document.querySelector(".browse .tiles");
  const params = new URLSearchParams(location.search);
  const f = {
    q: (params.get("q") || "").trim(),
    type: params.get("type") || "",
    family: params.get("family") || "",
    tag: params.get("tag") || "",
    fidelity: params.get("fidelity") || "",
    code: params.get("code") || "",
    author: params.get("author") || "",
    sort: params.get("sort") || "updated",
    provides: params.getAll("provides").filter(Boolean),
  };
  if (!["updated", "created", "name"].includes(f.sort)) f.sort = "updated";

  // The site may live under a path (a GitHub project page); find it from the stylesheet.
  const css = document.querySelector('link[href$="/static/chemart.css"]');
  const base = css ? css.getAttribute("href").replace(/\/static\/chemart\.css$/, "") : "";
  const browse = base + "/browse/";

  function link(change) {
    const m = Object.assign({}, f, change);
    const qs = new URLSearchParams();
    for (const key of ["q", "type", "family", "tag", "fidelity", "code", "author"]) {
      if (m[key]) qs.append(key, m[key]);
    }
    if (m.sort && m.sort !== "updated") qs.append("sort", m.sort);
    for (const v of m.provides) qs.append("provides", v);
    const s = qs.toString();
    return browse + (s ? "?" + s : "");
  }

  document.querySelectorAll("a.chip[data-key]").forEach((a) => {
    const key = a.dataset.key, value = a.dataset.value;
    let on, next;
    if (key === "provides") {
      on = f.provides.includes(value);
      next = { provides: on ? f.provides.filter((v) => v !== value) : f.provides.concat([value]) };
    } else if (key === "type" || key === "code") {
      on = f[key] === value;
      next = { [key]: value };
    } else {
      on = f[key] === value;
      next = { [key]: on ? "" : value };
    }
    a.classList.toggle("on", on);
    a.href = link(next);
  });

  const select = document.getElementById("sort");
  if (select) {
    select.value = f.sort;
    const go = (event) => { if (event) event.preventDefault(); location.href = link({ sort: select.value }); };
    select.addEventListener("change", () => go());
    if (select.form) select.form.addEventListener("submit", go);
  }
  document.querySelectorAll('input[name="q"]').forEach((input) => { if (!input.value) input.value = f.q; });

  const title = document.getElementById("browse-title");
  if (title) title.textContent = f.q ? "“" + f.q + "”" : (f.author || "Browse the shelves");
  const filtered = f.q || f.type || f.family || f.tag || f.fidelity || f.code || f.author || f.provides.length;
  const clear = document.getElementById("browse-clear");
  if (clear) clear.hidden = !filtered;

  if (!grid) return;
  const tiles = Array.from(grid.querySelectorAll("a.tile"));

  function apply(texts) {
    const terms = f.q.toLowerCase().split(/\s+/).filter(Boolean);
    const shown = tiles.filter((t) => {
      const d = t.dataset;
      if (f.type && d.type !== f.type) return false;
      if (f.family && d.family !== f.family) return false;
      if (f.fidelity && d.fidelity !== f.fidelity) return false;
      if (f.code && d.code !== f.code) return false;
      if (f.author && d.author !== f.author.toLowerCase()) return false;
      if (f.tag && !d.tags.split(" ").includes(f.tag)) return false;
      const provides = d.provides.split(" ");
      if (!f.provides.every((p) => provides.includes(p))) return false;
      if (terms.length) {
        const hay = (d.id + " " + t.textContent + " " + d.tags + " " + (texts[d.id] || "")).toLowerCase();
        if (!terms.every((term) => hay.includes(term))) return false;
      }
      return true;
    });
    const key = { updated: "updated", created: "created" }[f.sort];
    shown.sort((a, b) => key
      ? (b.dataset[key] || "").localeCompare(a.dataset[key] || "")
      : a.dataset.id.split("/")[1].localeCompare(b.dataset.id.split("/")[1]));
    tiles.forEach((t) => { t.hidden = true; });
    shown.forEach((t) => { t.hidden = false; grid.appendChild(t); });

    const count = document.getElementById("browse-count");
    if (count) count.textContent = shown.length + " result" + (shown.length === 1 ? "" : "s");
    let empty = document.getElementById("browse-empty");
    if (!shown.length && !empty) {
      empty = document.createElement("div");
      empty.id = "browse-empty";
      empty.className = "empty";
      empty.textContent = "No repos match these filters.";
      grid.after(empty);
    }
    if (empty) empty.hidden = shown.length > 0;
  }

  if (!f.q) { apply({}); return; }
  fetch(base + "/api/v1/index.json")
    .then((r) => r.json())
    .then((data) => {
      const texts = {};
      for (const repo of data.repos || []) texts[repo.id] = repo.text || "";
      apply(texts);
    })
    .catch(() => apply({}));
})();
