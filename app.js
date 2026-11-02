(function () {
  "use strict";
  var $ = function (s, r) { return (r || document).querySelector(s); };
  var $$ = function (s, r) { return Array.prototype.slice.call((r || document).querySelectorAll(s)); };
  var reduce = window.matchMedia && matchMedia("(prefers-reduced-motion: reduce)").matches;
  var TOOLS = window.TT_TOOLS, CATS = window.TT_CATS, RES = window.TT_RESULTS;
  var catOf = function (id) { return CATS.filter(function (c) { return c.id === id; })[0]; };
  var esc = function (s) { return String(s).replace(/[&<>"]/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]; }); };

  /* ---------- nav ---------- */
  var burger = $("#burger"), menu = $("#menu");
  burger.addEventListener("click", function () {
    var o = menu.classList.toggle("open");
    burger.setAttribute("aria-expanded", o);
  });
  $$("a", menu).forEach(function (a) { a.addEventListener("click", function () { menu.classList.remove("open"); burger.setAttribute("aria-expanded", false); }); });

  /* ---------- hero: a matrix of genotype calls (0 / 1 / 2) ---------- */
  (function () {
    var cv = $("#matrix"), ctx = cv.getContext("2d"), cell = 26, cols, rows, grid, w, h, mouse = { x: -999, y: -999 }, raf;
    var colors = ["rgba(92,197,154,.20)", "rgba(92,197,154,.50)", "rgba(240,180,41,.75)"];
    function size() {
      var r = cv.getBoundingClientRect(), d = Math.min(window.devicePixelRatio || 1, 2);
      w = r.width; h = r.height; cv.width = w * d; cv.height = h * d; ctx.setTransform(d, 0, 0, d, 0, 0);
      cols = Math.ceil(w / cell); rows = Math.ceil(h / cell);
      grid = []; for (var i = 0; i < cols * rows; i++) grid.push(Math.random() < .55 ? 0 : (Math.random() < .6 ? 1 : 2));
    }
    function draw(t) {
      ctx.clearRect(0, 0, w, h); ctx.font = "600 12px ui-monospace,Consolas,monospace"; ctx.textAlign = "center"; ctx.textBaseline = "middle";
      var scan = ((t || 0) / 16) % (w + 300) - 150;
      for (var y = 0; y < rows; y++) for (var x = 0; x < cols; x++) {
        var px = x * cell + cell / 2, py = y * cell + cell / 2, v = grid[y * cols + x];
        var d = Math.hypot(px - mouse.x, py - mouse.y), near = Math.max(0, 1 - d / 150), s = Math.max(0, 1 - Math.abs(px - scan) / 90);
        ctx.globalAlpha = Math.min(1, .5 + near + s);
        ctx.fillStyle = (near > .35 || s > .5) ? "#f0b429" : colors[v];
        ctx.fillText(String(v), px, py);
        if (!reduce && Math.random() < .00035) grid[y * cols + x] = (v + 1 + (Math.random() < .5 ? 1 : 0)) % 3;
      }
      ctx.globalAlpha = 1;
      if (!reduce) raf = requestAnimationFrame(draw);
    }
    size(); draw(0);
    window.addEventListener("resize", function () { size(); if (reduce) draw(0); });
    cv.parentNode.addEventListener("mousemove", function (e) { var r = cv.getBoundingClientRect(); mouse.x = e.clientX - r.left; mouse.y = e.clientY - r.top; });
    cv.parentNode.addEventListener("mouseleave", function () { mouse.x = mouse.y = -999; });
    document.addEventListener("visibilitychange", function () { if (document.hidden) cancelAnimationFrame(raf); else if (!reduce) raf = requestAnimationFrame(draw); });
  })();

  /* tilt of the app window */
  (function () {
    var box = $("#tilt"), win = $(".win", box);
    if (reduce) { win.style.transform = "none"; return; }
    box.addEventListener("mousemove", function (e) {
      var r = box.getBoundingClientRect(), x = (e.clientX - r.left) / r.width - .5, y = (e.clientY - r.top) / r.height - .5;
      win.style.transform = "rotateY(" + (x * 12 - 3) + "deg) rotateX(" + (-y * 9 + 2) + "deg)";
    });
    box.addEventListener("mouseleave", function () { win.style.transform = ""; });
  })();

  /* scroll reveal */
  (function () {
    var els = $$(".sec .wrap > *, .hero-text > *");
    els.forEach(function (e) { e.classList.add("reveal"); });
    if (!("IntersectionObserver" in window)) { els.forEach(function (e) { e.classList.add("in"); }); return; }
    var io = new IntersectionObserver(function (en) { en.forEach(function (x) { if (x.isIntersecting) { x.target.classList.add("in"); io.unobserve(x.target); } }); }, { threshold: .08 });
    els.forEach(function (e) { io.observe(e); });
    setTimeout(function () { els.forEach(function (e) { e.classList.add("in"); }); }, 2500);
  })();

  /* ---------- tool explorer ---------- */
  var grid = $("#toolGrid"), chips = $("#catChips"), q = $("#q"), cat = "all", visible = TOOLS.slice();
  function chipHtml(id, label, color, n) {
    return '<button role="tab" data-c="' + id + '" aria-selected="' + (id === cat) + '"' + (color ? ' style="--c:' + color + '"' : '') + '>' + (color ? "<i></i>" : "") + esc(label) + " <b>" + n + "</b></button>";
  }
  chips.innerHTML = chipHtml("all", "All tools", "", TOOLS.length) + CATS.map(function (c) { return chipHtml(c.id, c.label, c.color, TOOLS.filter(function (t) { return t.cat === c.id; }).length); }).join("");
  grid.innerHTML = TOOLS.map(function (t) {
    return '<button class="card" data-id="' + t.id + '" style="--c:' + catOf(t.cat).color + '"><span class="num">' + t.n + '</span><h3>' + esc(t.name) + '</h3><p>' + esc(t.desc) +
      '</p><span class="row"><span>Demo: ' + esc(t.time) + '</span><span class="go">Open &rarr;</span></span></button>';
  }).join("");
  function filter() {
    var s = q.value.trim().toLowerCase();
    visible = TOOLS.filter(function (t) { return (cat === "all" || t.cat === cat) && (!s || (t.name + " " + t.desc + " " + t.see.join(" ")).toLowerCase().indexOf(s) > -1); });
    var ids = visible.map(function (t) { return t.id; });
    $$(".card", grid).forEach(function (c) { c.classList.toggle("hide", ids.indexOf(c.dataset.id) < 0); });
    $("#empty").hidden = visible.length > 0;
  }
  chips.addEventListener("click", function (e) {
    var b = e.target.closest("button"); if (!b) return;
    cat = b.dataset.c; $$("button", chips).forEach(function (x) { x.setAttribute("aria-selected", x === b); }); filter();
  });
  q.addEventListener("input", filter);

  /* tool dialog */
  var modal = $("#modal"), card = $(".m-card", modal), cur = null, lastFocus = null;
  function open(id, push) {
    var t = TOOLS.filter(function (x) { return x.id === id; })[0]; if (!t) return;
    cur = t; var c = catOf(t.cat);
    modal.style.setProperty("--c", c.color);
    $("#mImg").src = "assets/shots/" + t.id + ".jpg"; $("#mImg").alt = t.name + " results";
    $("#mCap").textContent = ["phenotypic", "mendelian", "qtl_seq", "et_gwas", "pybsa", "clustering", "adom", "breedsight", "ocs", "heterosis"].indexOf(t.id) > -1 ? "Figure written by the tool on the demo files" : "The tool in TerraTools, on the demo files";
    $("#mTag").textContent = "#" + t.n + " · " + c.label; $("#mTag").style.setProperty("--c", c.color);
    $("#mTitle").textContent = t.name; $("#mDesc").textContent = t.desc; $("#mDemo").textContent = t.demo;
    $("#mSee").innerHTML = t.see.map(function (s) { return "<li>" + esc(s) + "</li>"; }).join("");
    $("#mTime").textContent = "Run time on the demo: " + t.time; $("#mSteps").textContent = t.steps;
    $("#mZip").href = "downloads/demo/" + t.zip + ".zip";
    lastFocus = document.activeElement; modal.hidden = false; document.body.style.overflow = "hidden"; card.focus();
    if (push !== false) history.replaceState(null, "", "#tool=" + id);
  }
  function close() { modal.hidden = true; document.body.style.overflow = ""; history.replaceState(null, "", location.pathname + location.search + "#tools"); if (lastFocus) lastFocus.focus(); }
  grid.addEventListener("click", function (e) { var c = e.target.closest(".card"); if (c) open(c.dataset.id); });
  modal.addEventListener("click", function (e) { if (e.target.hasAttribute("data-close")) close(); });
  $("#mNext").addEventListener("click", function () {
    var list = visible.length ? visible : TOOLS, i = list.indexOf(cur); open(list[(i + 1) % list.length].id);
  });
  document.addEventListener("keydown", function (e) {
    if (modal.hidden) return;
    if (e.key === "Escape") close();
    if (e.key === "ArrowRight") $("#mNext").click();
    if (e.key === "Tab") { var f = $$("button,a[href]", card); if (!f.length) return; var a = f[0], z = f[f.length - 1]; if (e.shiftKey && document.activeElement === a) { e.preventDefault(); z.focus(); } else if (!e.shiftKey && document.activeElement === z) { e.preventDefault(); a.focus(); } }
  });
  var m = location.hash.match(/tool=([\w]+)/); if (m) open(m[1], false);

  /* ---------- demo list ---------- */
  $("#demoList").innerHTML = TOOLS.map(function (t) {
    return '<a class="demo-row" style="--c:' + catOf(t.cat).color + '" href="downloads/demo/' + t.zip + '.zip" download><span class="num">' + t.n + '</span><div><b>' + esc(t.name) + '</b><small>' + esc(t.time) + " &middot; " + esc(t.see[0]) + '</small></div><span class="z">ZIP</span></a>';
  }).join("");

  /* ---------- workflow ---------- */
  (function () {
    var F = window.TT_FLOW, box = $("#flow"), det = $("#flowDetail");
    var name = function (id) { return TOOLS.filter(function (t) { return t.id === id; })[0]; };
    box.innerHTML = F.map(function (f, i) { var t = name(f.id); return '<button class="fnode" data-i="' + i + '"><span class="n">Step ' + (i + 1) + '</span><h3>' + esc(t.name) + '</h3></button>'; }).join("");
    function show(i) {
      $$(".fnode", box).forEach(function (n, j) { n.classList.toggle("on", j === i); });
      var f = F[i], next = F[i + 1] && name(F[i + 1].id);
      det.innerHTML = '<div class="file">' + esc(f.out) + '</div><p>' + esc(f.note) + (next ? " It is the input of <b>" + esc(next.name) + "</b>." : " This is the end of the chain: a ranked list you can take to the lab.") + '</p>';
    }
    box.addEventListener("click", function (e) { var b = e.target.closest(".fnode"); if (b) show(+b.dataset.i); });
    show(0);
    var i = 0, t = null;
    if (!reduce) { t = setInterval(function () { i = (i + 1) % F.length; show(i); }, 3200); box.addEventListener("click", function () { clearInterval(t); }); det.addEventListener("click", function () { clearInterval(t); }); }
  })();

  /* ---------- download ---------- */
  $("#copySha").addEventListener("click", function () {
    var b = this, txt = $("#sha").textContent;
    (navigator.clipboard ? navigator.clipboard.writeText(txt) : Promise.reject()).then(function () { b.textContent = "Copied"; }, function () {
      var r = document.createRange(); r.selectNode($("#sha")); getSelection().removeAllRanges(); getSelection().addRange(r); b.textContent = "Selected";
    });
    setTimeout(function () { b.textContent = "Copy"; }, 1600);
  });
  $("#heroDl").addEventListener("click", function () { var d = $("#dlBtn"); setTimeout(function () { d.animate([{ transform: "scale(1)" }, { transform: "scale(1.04)" }, { transform: "scale(1)" }], { duration: 700, iterations: 2 }); }, 500); });

  /* counters */
  (function () {
    if (reduce) return;
    $$("[data-count]").forEach(function (el) {
      var to = +el.dataset.count, pre = el.dataset.prefix ? el.dataset.prefix.replace("&lt;", "<") : "", t0 = performance.now();
      (function tick(t) { var p = Math.min(1, (t - t0) / 1200); el.textContent = pre + Math.round(to * (1 - Math.pow(1 - p, 3))); if (p < 1) requestAnimationFrame(tick); })(t0);
    });
  })();

  /* ---------- interactive plots ---------- */
  function setup(cv) {
    var d = Math.min(window.devicePixelRatio || 1, 2), r = cv.getBoundingClientRect();
    cv.width = Math.round(r.width * d); cv.height = Math.round(r.height * d);
    var c = cv.getContext("2d"); c.setTransform(d, 0, 0, d, 0, 0); return { c: c, w: r.width, h: r.height };
  }
  function axisY(c, L, T, W, H, lo, hi, step, label) {
    c.font = "12px DM Sans, sans-serif"; c.fillStyle = "#5a6f64"; c.strokeStyle = "#e4ece7"; c.lineWidth = 1; c.textAlign = "right"; c.textBaseline = "middle";
    for (var v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) {
      var y = T + H - (v - lo) / (hi - lo) * H; c.beginPath(); c.moveTo(L, y); c.lineTo(L + W, y); c.stroke(); c.fillText(String(Math.round(v * 100) / 100), L - 8, y);
    }
    c.save(); c.translate(14, T + H / 2); c.rotate(-Math.PI / 2); c.textAlign = "center"; c.fillText(label, 0, 0); c.restore();
  }
  function placeTip(tip, wrap, x, y, html) {
    tip.innerHTML = html; tip.hidden = false;
    var w = tip.offsetWidth / 2; tip.style.left = Math.max(w, Math.min(wrap.clientWidth - w, x)) + "px"; tip.style.top = Math.max(46, y) + "px";
  }

  // ET-GWAS Manhattan
  var et = (function () {
    var cv = $("#cvEt"), wrap = cv.parentNode, tip = $("#tipEt"), pts = RES.et.snps, S, hover = -1, chrLen = {}, off = {}, tot = 0, ymax = 0;
    pts.forEach(function (p) { chrLen[p[0]] = Math.max(chrLen[p[0]] || 0, p[1]); ymax = Math.max(ymax, p[2]); });
    var chrs = Object.keys(chrLen).map(Number).sort(function (a, b) { return a - b; });
    chrs.forEach(function (c) { off[c] = tot; tot += chrLen[c] + 1.5e6; });
    var bonf = -Math.log10(0.05 / pts.length), hi = Math.ceil(ymax + 1);
    var L = 52, T = 14, B = 34, pos = [];
    function draw() {
      S = setup(cv); var c = S.c, W = S.w - L - 10, H = S.h - T - B; pos = [];
      c.clearRect(0, 0, S.w, S.h); axisY(c, L, T, W, H, 0, hi, 2, "−log₁₀(p)");
      var show = $("#etPlanted").checked;
      var yb = T + H - bonf / hi * H; c.setLineDash([6, 5]); c.strokeStyle = "#d9534f"; c.lineWidth = 1.4; c.beginPath(); c.moveTo(L, yb); c.lineTo(L + W, yb); c.stroke(); c.setLineDash([]);
      c.fillStyle = "#d9534f"; c.font = "12px DM Sans"; c.textAlign = "right"; c.fillText("Bonferroni", L + W, yb - 8);
      c.textAlign = "center"; c.textBaseline = "top"; c.fillStyle = "#5a6f64";
      chrs.forEach(function (ch) { c.fillText(ch, L + (off[ch] + chrLen[ch] / 2) / tot * W, T + H + 8); });
      pts.forEach(function (p, i) {
        var x = L + (off[p[0]] + p[1]) / tot * W, y = T + H - p[2] / hi * H, sig = p[2] > bonf;
        pos.push([x, y]);
        c.beginPath(); c.arc(x, y, sig ? 6 : 3.2, 0, 7);
        c.fillStyle = sig ? "#d9534f" : (p[0] % 2 ? "#1f7a57" : "#6aa4d9"); c.globalAlpha = sig ? 1 : .72; c.fill(); c.globalAlpha = 1;
        if (i === hover) { c.lineWidth = 2.5; c.strokeStyle = "#0a2a20"; c.stroke(); }
        if (show && p[6]) { c.beginPath(); c.arc(x, y, 12, 0, 7); c.lineWidth = 2.5; c.strokeStyle = "#f0b429"; c.stroke(); c.font = "700 12px DM Sans"; c.fillStyle = "#8a5d00"; c.textAlign = "center"; c.textBaseline = "bottom"; c.fillText(p[3] < .2 ? "rare allele" : "common allele", x, y - 16); }
      });
    }
    function move(e) {
      var r = cv.getBoundingClientRect(), mx = e.clientX - r.left, my = e.clientY - r.top, best = -1, bd = 14;
      pos.forEach(function (p, i) { var d = Math.hypot(p[0] - mx, p[1] - my); if (d < bd) { bd = d; best = i; } });
      if (best !== hover) { hover = best; draw(); }
      if (best < 0) { tip.hidden = true; return; }
      var p = pts[best];
      placeTip(tip, wrap, pos[best][0], pos[best][1], "<b>" + p[4] + "</b><br>MAF " + p[3] + " &middot; −log₁₀p " + p[2].toFixed(2) + "<br>q-value " + p[5]);
    }
    cv.addEventListener("mousemove", move); cv.addEventListener("click", move); cv.addEventListener("mouseleave", function () { hover = -1; tip.hidden = true; draw(); });
    $("#etPlanted").addEventListener("change", function () { draw(); $("#etHint").textContent = this.checked ? "Gold rings: the two SNPs the trait was built on. The tool found both." : "Hover a point to read it."; });
    return draw;
  })();

  // QTL-Seq delta SNP-index
  var qs = (function () {
    var cv = $("#cvQs"), wrap = cv.parentNode, tip = $("#tipQs"), lvl = 3, names = Object.keys(RES.qtlseq), hov = null, S, geo = [];
    var lab = { 2: "90%", 3: "95%", 4: "99%" }, maxX = 22;
    function draw() {
      S = setup(cv); var c = S.c, L = 52, T = 14, B = 34, gap = 12, pw = (S.w - L - 8 - gap * (names.length - 1)) / names.length, H = S.h - T - B; geo = [];
      c.clearRect(0, 0, S.w, S.h);
      axisY(c, L, T, S.w - L - 8, H, -1, 1, .5, "Δ SNP-index");
      var out = {};
      names.forEach(function (n, k) {
        var x0 = L + k * (pw + gap), d = RES.qtlseq[n], X = function (v) { return x0 + v / maxX * pw; }, Y = function (v) { return T + H - (v + 1) / 2 * H; };
        geo.push({ n: n, x0: x0, pw: pw });
        c.fillStyle = k % 2 ? "rgba(31,122,87,.04)" : "rgba(61,127,196,.05)"; c.fillRect(x0, T, pw, H);
        // confidence band
        c.beginPath(); d.forEach(function (p, i) { c[i ? "lineTo" : "moveTo"](X(p[0]), Y(p[lvl])); }); for (var i = d.length - 1; i >= 0; i--) c.lineTo(X(d[i][0]), Y(-d[i][lvl])); c.closePath();
        c.fillStyle = "#f7ebc8"; c.fill();
        c.strokeStyle = "rgba(217,83,79,.7)"; c.lineWidth = 1.2; c.beginPath(); d.forEach(function (p, i) { c[i ? "lineTo" : "moveTo"](X(p[0]), Y(p[lvl])); }); c.stroke();
        c.beginPath(); d.forEach(function (p, i) { c[i ? "lineTo" : "moveTo"](X(p[0]), Y(-p[lvl])); }); c.stroke();
        // data
        c.beginPath(); d.forEach(function (p, i) { c[i ? "lineTo" : "moveTo"](X(p[0]), Y(p[1])); }); c.strokeStyle = "#0a2a20"; c.lineWidth = 1.6; c.stroke();
        var n_out = 0;
        d.forEach(function (p) { if (Math.abs(p[1]) > p[lvl]) { n_out++; c.beginPath(); c.arc(X(p[0]), Y(p[1]), 2.4, 0, 7); c.fillStyle = "#d9534f"; c.fill(); } });
        out[n] = n_out;
        c.fillStyle = "#16261f"; c.font = "700 13px DM Sans"; c.textAlign = "center"; c.textBaseline = "top"; c.fillText(n.replace("chr", "Chromosome "), x0 + pw / 2, T + H + 8);
        c.fillStyle = "#5a6f64"; c.font = "11px DM Sans"; c.textAlign = "left"; c.fillText("0", x0, T + H + 22 - 2); c.textAlign = "right"; c.fillText("22 Mb", x0 + pw, T + H + 22 - 2);
        if (hov && hov.n === n) { var hx = X(hov.p[0]); c.strokeStyle = "rgba(10,42,32,.5)"; c.setLineDash([4, 4]); c.beginPath(); c.moveTo(hx, T); c.lineTo(hx, T + H); c.stroke(); c.setLineDash([]); c.beginPath(); c.arc(hx, Y(hov.p[1]), 5, 0, 7); c.fillStyle = "#f0b429"; c.fill(); c.strokeStyle = "#0a2a20"; c.lineWidth = 2; c.stroke(); }
      });
      $("#qsHint").textContent = "Windows beyond the " + lab[lvl] + " limit: " + names.map(function (n) { return n.replace("chr", "chr") + " " + out[n]; }).join(" · ");
    }
    cv.addEventListener("mousemove", function (e) {
      var r = cv.getBoundingClientRect(), mx = e.clientX - r.left, g = geo.filter(function (g) { return mx >= g.x0 && mx <= g.x0 + g.pw; })[0];
      if (!g) { hov = null; tip.hidden = true; draw(); return; }
      var mb = (mx - g.x0) / g.pw * maxX, d = RES.qtlseq[g.n], best = d[0];
      d.forEach(function (p) { if (Math.abs(p[0] - mb) < Math.abs(best[0] - mb)) best = p; });
      hov = { n: g.n, p: best }; draw();
      placeTip(tip, wrap, g.x0 + best[0] / maxX * g.pw, e.clientY - r.top, "<b>" + g.n + " &middot; " + best[0].toFixed(1) + " Mb</b><br>Δ SNP-index " + best[1].toFixed(2) + "<br>" + lab[lvl] + " limit ±" + best[lvl].toFixed(2) + (Math.abs(best[1]) > best[lvl] ? " <b>beyond</b>" : ""));
    });
    cv.addEventListener("mouseleave", function () { hov = null; tip.hidden = true; draw(); });
    $("#qsLevel").addEventListener("click", function (e) { var b = e.target.closest("button"); if (!b) return; lvl = +b.dataset.l; $$("button", this).forEach(function (x) { x.classList.toggle("on", x === b); }); draw(); });
    return draw;
  })();

  // Diversity strip plot
  var dv = (function () {
    var cv = $("#cvDv"), wrap = cv.parentNode, tip = $("#tipDv"), names = Object.keys(RES.div), hov = -1, S, rows = [], seed = 7;
    function rnd() { seed = (seed * 16807) % 2147483647; return seed / 2147483647; }
    var jit = names.map(function (n) { return RES.div[n].shannon.map(function () { return rnd() - .5; }); });
    var pal = ["#6aa4d9", "#2d9d78", "#e0a526", "#c9743a", "#a24b6b"];
    function draw() {
      S = setup(cv); var c = S.c, L = 112, R = 22, T = 14, B = 38, W = S.w - L - R, H = S.h - T - B, rh = H / names.length, lo = 0, hi = 1.3; rows = [];
      c.clearRect(0, 0, S.w, S.h); c.font = "12px DM Sans"; c.textBaseline = "top"; c.textAlign = "center"; c.strokeStyle = "#e4ece7"; c.fillStyle = "#5a6f64";
      for (var v = 0; v <= 1.2; v += .2) { var x = L + v / hi * W; c.beginPath(); c.moveTo(x, T); c.lineTo(x, T + H); c.stroke(); c.fillText(v.toFixed(1), x, T + H + 8); }
      c.fillText("Shannon diversity of each sample", L + W / 2, T + H + 22);
      names.forEach(function (n, i) {
        var d = RES.div[n], cy = T + rh * (i + .5), mean = d.shannon.reduce(function (a, b) { return a + b; }, 0) / d.shannon.length; rows.push([T + rh * i, T + rh * (i + 1), n, mean]);
        if (hov === i) { c.fillStyle = "rgba(240,180,41,.13)"; c.fillRect(0, T + rh * i, S.w, rh); }
        c.fillStyle = "#16261f"; c.font = "600 14px DM Sans"; c.textAlign = "right"; c.textBaseline = "middle"; c.fillText(n.replace("_", " "), L - 14, cy);
        d.shannon.forEach(function (s, k) { c.beginPath(); c.arc(L + s / hi * W, cy + jit[i][k] * rh * .62, 4.4, 0, 7); c.fillStyle = pal[i]; c.globalAlpha = .75; c.fill(); });
        c.globalAlpha = 1; var mx = L + mean / hi * W; c.fillStyle = "#0a2a20"; c.fillRect(mx - 1.5, cy - rh * .38, 3, rh * .76);
      });
    }
    cv.addEventListener("mousemove", function (e) {
      var r = cv.getBoundingClientRect(), my = e.clientY - r.top, i = -1; rows.forEach(function (x, k) { if (my >= x[0] && my < x[1]) i = k; });
      if (i !== hov) { hov = i; draw(); }
      if (i < 0) { tip.hidden = true; return; }
      var d = RES.div[names[i]]; placeTip(tip, wrap, e.clientX - r.left, (rows[i][0] + rows[i][1]) / 2 - 16, "<b>" + names[i].replace("_", " ") + "</b><br>mean Shannon " + rows[i][3].toFixed(2) + "<br>effective alleles " + d.eff + " &middot; heterozygosity " + d.het);
    });
    cv.addEventListener("mouseleave", function () { hov = -1; tip.hidden = true; draw(); });
    return draw;
  })();

  var draws = { et: et, qs: qs, dv: dv };
  function showPlot(p) {
    $$("#plotTabs button").forEach(function (b) { b.setAttribute("aria-selected", b.dataset.p === p); });
    ["et", "qs", "dv"].forEach(function (k) { $("#pane-" + k).hidden = k !== p; });
    draws[p]();
  }
  $("#plotTabs").addEventListener("click", function (e) { var b = e.target.closest("button"); if (b) showPlot(b.dataset.p); });
  var cp = "et";
  $$("#plotTabs button").forEach(function (b) { b.addEventListener("click", function () { cp = b.dataset.p; }); });
  window.addEventListener("resize", function () { draws[cp](); });
  showPlot("et");
  if ("fonts" in document) document.fonts.ready.then(function () { draws[cp](); });
})();
