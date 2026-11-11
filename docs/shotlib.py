"""
shotlib.py — drives the real TerraTools window (offscreen) and captures annotated screenshots
for the user manual.  Scenarios are generator functions that `yield <milliseconds>` to wait.

    from shotlib import App
    app = App(out_dir)
    app.run([scenario1, scenario2, ...])
"""
from __future__ import annotations

import json
import os
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QTWEBENGINE_CHROMIUM_FLAGS", "--no-sandbox --disable-gpu")
os.environ.setdefault("QT_SCALE_FACTOR", "1.5")

from PIL import Image, ImageDraw, ImageFont   # noqa: E402
from PySide6.QtCore import QTimer             # noqa: E402
from PySide6.QtWidgets import QApplication, QFileDialog   # noqa: E402

ACCENT = (230, 90, 20)
FONT = None
for f in (r"C:\Windows\Fonts\segoeuib.ttf", r"C:\Windows\Fonts\arialbd.ttf"):
    if Path(f).exists():
        FONT = f
        break

LIB = r"""
window.__q = function(spec){
  const root = (spec.global || spec.global_) ? document : (document.querySelector('.view.active') || document);
  let els = [...root.querySelectorAll(spec.sel || '*')].filter(e => e.getClientRects().length && e.offsetWidth > 0);
  if (spec.text) {
    const t = spec.text.toLowerCase();
    els = els.filter(e => (e.innerText || '').toLowerCase().includes(t));
    els = els.filter(e => !els.some(o => o !== e && e.contains(o)));
  }
  let el = els[spec.nth || 0];
  if (!el) return null;
  if (spec.closest) el = el.closest(spec.closest) || el;
  for (let i = 0; i < (spec.parent || 0); i++) el = el.parentElement || el;
  return el;
};
window.__rect = function(spec){
  const el = window.__q(spec); if (!el) return null;
  const r = el.getBoundingClientRect();
  return {x: r.left, y: r.top, w: r.width, h: r.height};
};
window.__click = function(spec){ const el = window.__q(spec); if (!el) return false; el.scrollIntoView({block:'center'}); el.click(); return true; };
window.__setSelect = function(label, value){
  const root = document.querySelector('.view.active') || document;
  const f = [...root.querySelectorAll('.field')].find(f => (f.querySelector('.field-title')||{innerText:''}).innerText.toLowerCase().includes(label.toLowerCase()) && f.querySelector('select'));
  if (!f) return false;
  const s = f.querySelector('select');
  Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value').set.call(s, value);
  s.dispatchEvent(new Event('change', {bubbles: true})); return true;
};
window.__setInput = function(label, value){
  const root = document.querySelector('.view.active') || document;
  const f = [...root.querySelectorAll('.field')].find(f => (f.querySelector('.field-title')||{innerText:''}).innerText.toLowerCase().includes(label.toLowerCase()) && f.querySelector('input'));
  if (!f) return false;
  const s = f.querySelector('input');
  Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set.call(s, String(value));
  s.dispatchEvent(new Event('input', {bubbles: true})); return true;
};
window.__setSlider = function(label, value){
  const root = document.querySelector('.view.active') || document;
  const s = [...root.querySelectorAll('input[type=range]')].find(i => (i.getAttribute('aria-label')||'').toLowerCase().includes(label.toLowerCase()));
  if (!s) return false;
  Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value').set.call(s, String(value));
  s.dispatchEvent(new Event('input', {bubbles: true})); s.dispatchEvent(new Event('change', {bubbles: true})); return true;
};
window.__scrollTop = function(){ const v = document.querySelector('.view.active'); if (v) v.scrollTo(0, 0); return 1; };
true;
"""


class App:
    def __init__(self, out_dir):
        self.out = Path(out_dir)
        self.out.mkdir(parents=True, exist_ok=True)
        self.qt = QApplication(sys.argv)
        self.next_open = ""
        QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: (self.next_open, ""))
        QFileDialog.getOpenFileNames = staticmethod(lambda *a, **k: ([self.next_open], ""))
        from ui.main_window import MainWindow
        self.win = MainWindow()
        self.win.resize(1440, 1000)
        self.win.show()
        self.page = self.win.browser.page()
        self.bridges = self.win.tool_bridges
        mf = self.out / "manifest.json"
        self.manifest = json.loads(mf.read_text(encoding="utf-8")) if mf.exists() else {}   # partial re-runs merge in
        self.failures = []
        self.width = 1700

    # ---------- primitives (generators: use `yield from`) ----------
    def js(self, code):
        box = {"done": False, "v": None}

        def cb(v):
            box["v"], box["done"] = v, True
        self.page.runJavaScript(code, 0, cb)
        waited = 0
        while not box["done"] and waited < 20000:
            yield 30
            waited += 30
        return box["v"]

    def lib(self):
        yield from self.js(LIB)

    def call(self, expr):
        """Evaluate a __q helper call and return its JSON-decoded value."""
        yield from self.lib()
        v = yield from self.js(f"JSON.stringify({expr})")
        return json.loads(v) if isinstance(v, str) else v

    def route(self, r, height=1000, wait=1800):
        self.win.resize(1440, height)
        yield from self.js(f"location.replace('#/{r}'); 1")
        yield wait

    def size(self, height):
        self.win.resize(1440, height)
        yield 600

    def click(self, **spec):
        ok = yield from self.call(f"window.__click({json.dumps(spec)})")
        if not ok:
            raise RuntimeError(f"click target not found: {spec}")
        yield 500

    def tab(self, text, wait=900):
        yield from self.click(sel=".tab-btn", text=text)
        yield wait

    def browse(self, nth, path, wait=1500, sel=".file-drop button"):
        self.next_open = str(path)
        yield from self.click(sel=sel, nth=nth)
        yield wait

    def select(self, label, value, wait=500):
        ok = yield from self.call(f"window.__setSelect({json.dumps(label)}, {json.dumps(value)})")
        if not ok:
            raise RuntimeError(f"select not found: {label}")
        yield wait

    def text_input(self, label, value, wait=300):
        ok = yield from self.call(f"window.__setInput({json.dumps(label)}, {json.dumps(value)})")
        if not ok:
            raise RuntimeError(f"input not found: {label}")
        yield wait

    def slider(self, label, value, wait=300):
        ok = yield from self.call(f"window.__setSlider({json.dumps(label)}, {json.dumps(value)})")
        if not ok:
            raise RuntimeError(f"slider not found: {label}")
        yield wait

    def run_and_wait(self, timeout_s=300, settle=2800):
        """Click the run button and wait until the run finishes; returns the status text."""
        yield from self.click(sel=".run-bar .btn-run")
        waited = 0
        started = False
        while waited < timeout_s * 1000:
            running = yield from self.call("(function(){const b=document.querySelector('.view.active .run-bar .btn-run'); return b ? b.classList.contains('running') : null})()")
            if running:
                started = True
            elif started or waited > 6000:
                break
            yield 700
            waited += 700
        yield settle
        status = yield from self.call("(function(){const r=document.querySelector('.view.active .run-bar'); return r ? r.innerText : ''})()")
        return status

    def stat(self, name):
        return self.bridges[name]

    # ---------- screenshots ----------
    def shot(self, name, callouts=(), title="", height=None, focus=None):
        """Capture the window; `callouts` = [{spec..., note: '...'}] drawn as numbered boxes.
        height: window height in px; focus: spec of an element to scroll to the top of the view."""
        if height:
            self.win.resize(1440, height)
            yield 1800                      # the web view needs a moment to lay out at the new size
        yield from self.lib()
        if focus:
            yield from self.js("(function(){const el = window.__q(%s); const v = document.querySelector('.view.active');"
                               "if (el && v){ v.scrollTo(0, 0); const r = el.getBoundingClientRect(); v.scrollBy(0, r.top - 120);} return 1})()" % json.dumps(focus))
        else:
            yield from self.js("window.__scrollTop()")
        yield 700
        rects = []
        for c in callouts:
            spec = {k: v for k, v in c.items() if k not in ("note", "pad")}
            r = yield from self.call(f"window.__rect({json.dumps(spec)})")
            rects.append((c, r))
        yield 250
        pm = self.win.grab()
        path = self.out / f"{name}.png"
        pm.save(str(path))
        dpr = pm.devicePixelRatio()
        im = Image.open(path).convert("RGB")
        d = ImageDraw.Draw(im)
        font = ImageFont.truetype(FONT, int(15 * dpr)) if FONT else None
        notes, n = [], 0
        for c, r in rects:
            if not r:
                self.failures.append(f"{name}: callout not found {c.get('text') or c.get('sel')}")
                continue
            n += 1
            pad = c.get("pad", 4)
            x0, y0 = (r["x"] - pad) * dpr, (r["y"] - pad) * dpr
            x1, y1 = (r["x"] + r["w"] + pad) * dpr, (r["y"] + r["h"] + pad) * dpr
            d.rounded_rectangle([x0, y0, x1, y1], radius=8 * dpr / 1.5, outline=ACCENT, width=max(3, int(2.4 * dpr)))
            rad = 12 * dpr
            cx, cy = max(x0 - 2, rad + 2), max(y0 - 2, rad + 2)
            d.ellipse([cx - rad, cy - rad, cx + rad, cy + rad], fill=ACCENT, outline=(255, 255, 255), width=max(2, int(1.5 * dpr)))
            d.text((cx, cy), str(n), fill=(255, 255, 255), font=font, anchor="mm")
            notes.append(c.get("note", ""))
        w = self.width
        im = im.resize((w, int(im.height * w / im.width)), Image.LANCZOS)
        im.save(self.out / f"{name}.jpg", quality=86, optimize=True)
        path.unlink(missing_ok=True)
        self.manifest[name] = {"title": title, "notes": notes}
        print(f"[shot] {name}  ({n} callouts)", flush=True)

    # ---------- runner ----------
    def run(self, scenarios, boot_ms=7000):
        def driver():
            try:
                yield boot_ms
                for sc in scenarios:
                    label = getattr(sc, "__name__", str(sc))
                    print(f"=== {label}", flush=True)
                    try:
                        yield from sc(self)
                    except Exception as e:             # keep going; report at the end
                        self.failures.append(f"{label}: {e}")
                        traceback.print_exc()
            finally:
                (self.out / "manifest.json").write_text(json.dumps(self.manifest, indent=1), encoding="utf-8")
                print("\nFAILURES:" if self.failures else "\nno failures", *self.failures, sep="\n  ")
                self.qt.quit()
        gen = driver()

        def tick():
            try:
                delay = next(gen)
            except StopIteration:
                return
            QTimer.singleShot(max(int(delay), 1), tick)
        QTimer.singleShot(50, tick)
        self.qt.exec()
