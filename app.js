(function () {
  "use strict";
  var $ = function (s, r) { return (r || document).querySelector(s); };
  var $$ = function (s, r) { return Array.prototype.slice.call((r || document).querySelectorAll(s)); };
  var reduce = window.matchMedia && matchMedia("(prefers-reduced-motion: reduce)").matches;
  var TOOLS = window.TT_TOOLS, CATS = window.TT_CATS;
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
      '</p><span class="row"><span class="go">What it does &rarr;</span></span></button>';
  }).join("");
  function filter() {
    var s = q.value.trim().toLowerCase();
    visible = TOOLS.filter(function (t) { return (cat === "all" || t.cat === cat) && (!s || (t.name + " " + t.desc + " " + t.does.join(" ")).toLowerCase().indexOf(s) > -1); });
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
    $("#mImg").src = "assets/shots/" + t.id + ".jpg"; $("#mImg").alt = "The " + t.name + " tool in TerraTools";
    $("#mTag").textContent = "#" + t.n + " · " + c.label; $("#mTag").style.setProperty("--c", c.color);
    $("#mTitle").textContent = t.name; $("#mDesc").textContent = t.desc;
    $("#mSee").innerHTML = t.does.map(function (s) { return "<li>" + esc(s) + "</li>"; }).join("");
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
    return '<a class="demo-row" style="--c:' + catOf(t.cat).color + '" href="downloads/demo/' + t.zip + '.zip" download><span class="num">' + t.n + '</span><div><b>' + esc(t.name) + '</b><small>' + esc(t.desc) + '</small></div><span class="z">ZIP</span></a>';
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

})();
