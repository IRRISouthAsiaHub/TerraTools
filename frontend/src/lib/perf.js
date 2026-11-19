/* Frontend diagnostics and the motion guard.
 *
 *  - diag(ev, fields): queues a structured record; at most one batch per 4 s goes to Python (TerraToolsBridge.diag_log),
 *    which writes it to %LOCALAPPDATA%\TerraTools\logs\diagnostics.log next to the Python-side records.
 *  - long tasks (> 200 ms on the page's main thread) are reported by a PerformanceObserver (no polling).
 *  - motionGuard(): after start-up, samples ~2 s of animation frames. If the median frame takes longer than 30 ms the
 *    machine cannot keep the decorative animations smooth, so <html data-motion="reduced"> stops them (styles.css).
 *    Software graphics mode (reported by Python) starts reduced.
 */
const inQt = () => Boolean(window.qt && window.qt.webChannelTransport)   // (not imported from bridge.js, which imports this file)

let queue = []
let timer = null
let sender = null

export function setDiagSender(fn) { sender = fn }

export function diag(ev, fields = {}) {
  if (!inQt()) return
  queue.push({ ev, ...fields })
  if (!timer) {
    timer = setTimeout(() => {
      timer = null
      const batch = queue.splice(0, 50)
      if (batch.length && sender) { try { sender(JSON.stringify(batch)) } catch { /* diagnostics must never break the page */ } }
    }, 4000)
  }
}

export function watchLongTasks() {
  try {
    const po = new PerformanceObserver((list) => {
      for (const e of list.getEntries()) if (e.duration > 200) diag('long_task', { ms: Math.round(e.duration), at: Math.round(e.startTime) })
    })
    po.observe({ entryTypes: ['longtask'] })
  } catch { /* longtask not supported */ }
}

export function motionGuard() {
  const html = document.documentElement
  if (html.dataset.gfx === 'software') { html.dataset.motion = 'reduced'; diag('motion', { mode: 'reduced', why: 'software graphics' }); return }
  setTimeout(() => {
    if (document.visibilityState !== 'visible') return
    const ts = []
    let n = 0
    const loop = (t) => { ts.push(t); if (++n < 90 && document.visibilityState === 'visible') requestAnimationFrame(loop); else done() }
    const done = () => {
      if (ts.length < 30) return
      const d = ts.slice(1).map((t, i) => t - ts[i]).sort((a, b) => a - b)
      const p50 = d[Math.floor(d.length / 2)]
      const reduced = p50 > 30
      if (reduced) html.dataset.motion = 'reduced'
      diag('motion', { mode: reduced ? 'reduced' : 'full', p50_ms: Math.round(p50 * 10) / 10, p95_ms: Math.round(d[Math.floor(d.length * 0.95)]), frames: ts.length })
    }
    requestAnimationFrame(loop)
  }, 3500)
}
