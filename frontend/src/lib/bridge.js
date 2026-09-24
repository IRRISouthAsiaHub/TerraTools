import { useEffect, useRef, useState, useCallback } from 'react'

/* ------------------------------------------------------------------ *
 * QWebChannel access.
 *
 * Inside TerraTools every window registers exactly one bridge object on
 * its own channel (terraTools, gwasBridge, clusteringBridge, ...). The
 * same React build is loaded by all of them, so a page asks for the
 * object it needs by name.
 *
 * In a plain browser (npm run dev) there is no Qt, so a small mock is
 * returned instead so the UI can be developed and previewed standalone.
 * ------------------------------------------------------------------ */

let channelPromise = null

function loadChannel() {
  if (channelPromise) return channelPromise
  channelPromise = new Promise((resolve) => {
    if (!window.qt || !window.qt.webChannelTransport) return resolve(null)
    const connect = () =>
      new window.QWebChannel(window.qt.webChannelTransport, (ch) => resolve(ch))
    if (window.QWebChannel) return connect()
    const s = document.createElement('script')
    s.src = 'qrc:///qtwebchannel/qwebchannel.js'
    s.onload = connect
    s.onerror = () => resolve(null)
    document.head.appendChild(s)
  })
  return channelPromise
}

export const inQt = () => Boolean(window.qt && window.qt.webChannelTransport)

/** Turn a local filesystem path returned by Python into a loadable URL. */
export function toFileUrl(p) {
  if (!p) return ''
  if (/^(data:|https?:|file:|blob:)/i.test(p)) return p
  let s = p.replace(/\\/g, '/').replace(/#/g, '%23').replace(/\?/g, '%3F')
  if (/^[A-Za-z]:/.test(s)) return 'file:///' + s
  return 'file://' + s
}

/**
 * useBridge(objectName, mockFactory?)
 * Returns { bridge, ready }. `bridge` is null until connected.
 */
export function useBridge(objectName, mockFactory) {
  const [bridge, setBridge] = useState(null)
  useEffect(() => {
    let alive = true
    loadChannel().then((ch) => {
      if (!alive) return
      if (ch && ch.objects[objectName]) setBridge(ch.objects[objectName])
      else if (!ch && mockFactory) setBridge(mockFactory())
      else console.error(`[TerraTools] bridge "${objectName}" not found`,
        ch ? Object.keys(ch.objects) : '(no channel)')
    })
    return () => { alive = false }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [objectName])
  return { bridge, ready: Boolean(bridge) }
}

/**
 * useTool(bridge, runSlot)
 * Wraps the resultReady / statusUpdate / run_* / browse_file / save_file
 * contract that every tool bridge follows.
 *
 * state: 'idle' | 'running' | 'success' | 'error'
 */
export function useTool(bridge, runSlot, opts = {}) {
  const [state, setState] = useState('idle')
  const [status, setStatus] = useState('Ready')
  const [result, setResult] = useState(null)
  const resultCb = useRef(null)
  const accept = useRef(opts.accept)
  accept.current = opts.accept
  const onResult = useCallback((cb) => { resultCb.current = cb }, [])

  useEffect(() => {
    if (!bridge) return
    const onStatus = (t) => setStatus(t)
    const onDone = (json) => {
      try {
        const r = typeof json === 'string' ? JSON.parse(json) : json
        const ok = r.status === 'success'
        if (ok || (accept.current && accept.current(r))) {
          setResult(r)
          setState(ok ? 'success' : 'error')
          if (!ok) setStatus('Error: ' + (r.message || 'Unknown error'))
          resultCb.current && resultCb.current(r)
        } else {
          setState('error')
          setStatus('Error: ' + (r.message || 'Unknown error'))
        }
      } catch (err) {
        setState('error')
        setStatus('Frontend error: ' + err.message)
      }
    }
    bridge.statusUpdate.connect(onStatus)
    bridge.resultReady.connect(onDone)
    return () => {
      try { bridge.statusUpdate.disconnect(onStatus) } catch { /* noop */ }
      try { bridge.resultReady.disconnect(onDone) } catch { /* noop */ }
    }
  }, [bridge])

  /** browse(slotName = 'browse_file', ...args) -> Promise<string path> */
  const browse = useCallback(
    (slot = 'browse_file', ...args) =>
      new Promise((res) => (bridge ? bridge[slot](...args, (p) => res(p || '')) : res(''))),
    [bridge],
  )
  const save = useCallback(
    (src, name) => bridge && bridge.save_file(src, name),
    [bridge],
  )
  const run = useCallback(
    (config) => {
      if (!bridge) return
      setState('running')
      setStatus('Starting…')
      bridge[runSlot](JSON.stringify(config))
    },
    [bridge, runSlot],
  )
  const fail = useCallback((msg) => { setState('error'); setStatus(msg) }, [])

  return { bridge, state, status, setStatus, setState, result, setResult, run, browse, save, fail, onResult }
}

/**
 * Some Python bridges push results with page.runJavaScript("showXxx({...})").
 * useGlobals({ showXxx: fn }) exposes those names on window and always calls
 * the latest handler.
 */
export function useGlobals(handlers) {
  const ref = useRef(handlers)
  ref.current = handlers
  useEffect(() => {
    const names = Object.keys(ref.current)
    names.forEach((n) => { window[n] = (...a) => ref.current[n](...a) })
    return () => names.forEach((n) => { delete window[n] })
  }, [])
}

/** Subscribe to an extra bridge signal (dataLoaded, progress, ...). */
export function useSignal(bridge, name, handler) {
  const ref = useRef(handler)
  ref.current = handler
  useEffect(() => {
    if (!bridge || !bridge[name]) return
    const f = (...a) => ref.current(...a)
    bridge[name].connect(f)
    return () => { try { bridge[name].disconnect(f) } catch { /* noop */ } }
  }, [bridge, name])
}

/** Call a bridge slot that returns a value (QWebChannel passes a callback last). */
export const callSlot = (bridge, name, ...args) =>
  new Promise((res) => (bridge ? bridge[name](...args, res) : res(null)))

/** Same, but parses a JSON string return value. */
export const callJson = async (bridge, name, ...args) => {
  const r = await callSlot(bridge, name, ...args)
  try { return typeof r === 'string' ? JSON.parse(r) : r } catch { return null }
}

/** Load a local JSON file written by Python (fetch works for file:// in QtWebEngine). */
export async function loadJson(path) {
  const r = await fetch(toFileUrl(path))
  return r.json()
}

/* ------------------------------------------------------------------ *
 * Browser-only mocks (never used inside Qt).
 * ------------------------------------------------------------------ */

function signal() {
  const fns = new Set()
  return {
    connect: (f) => fns.add(f),
    disconnect: (f) => fns.delete(f),
    emit: (...a) => fns.forEach((f) => f(...a)),
  }
}

const svg = (title, hue) =>
  'data:image/svg+xml;utf8,' +
  encodeURIComponent(
    `<svg xmlns="http://www.w3.org/2000/svg" width="900" height="420" viewBox="0 0 900 420">
      <rect width="900" height="420" fill="#f5faf7"/>
      ${Array.from({ length: 140 }, (_, i) => {
        const x = 30 + (i * 6.2) % 840
        const y = 380 - Math.abs(Math.sin(i * 1.7 + hue)) * (i % 23 === 0 ? 300 : 120)
        return `<circle cx="${x}" cy="${y}" r="3.2" fill="hsl(${150 + hue * 10},45%,${35 + (i % 5) * 6}%)"/>`
      }).join('')}
      <text x="450" y="30" text-anchor="middle" font-family="sans-serif" font-size="16" fill="#1a3a2a">${title}</text>
    </svg>`,
  )

export function mockToolBridge(kind) {
  const resultReady = signal()
  const statusUpdate = signal()
  const bridge = {
    resultReady,
    statusUpdate,
    browse_file: (cb) => setTimeout(() => cb('C:\\data\\demo_genotypes.csv'), 300),
    save_file: (src, name) => statusUpdate.emit(`Saved ${name} (mock)`),
  }
  const run = () => {
    let i = 0
    const steps = ['Reading input…', 'Filtering markers…', 'Fitting models…', 'Rendering plots…']
    const t = setInterval(() => {
      if (i < steps.length) return statusUpdate.emit(steps[i++])
      clearInterval(t)
      if (kind === 'gwas') {
        resultReady.emit(JSON.stringify({
          status: 'success',
          plots: {
            '01_Manhattan_OLS': svg('Manhattan — OLS', 1),
            '02_Manhattan_MLM': svg('Manhattan — MLM', 2),
            '03_QQ_OLS': svg('QQ — OLS', 3),
            '04_QQ_MLM': svg('QQ — MLM', 4),
          },
          excel: 'C:\\tmp\\gwas.xlsx', zip: 'C:\\tmp\\gwas.zip',
          summary: 'Models run: OLS, MLM\nSignificant markers (Bonferroni): 7\nTop marker: Chr3_1289221 (p = 2.1e-9)',
        }))
      } else {
        resultReady.emit(JSON.stringify({
          status: 'success', n_samples: 240, n_features: 3120, n_clusters_used: 5,
          n_train: 192, n_test: 48, pc1_variance: 18.4, pc2_variance: 9.7, split_method: 'stratified',
          plot_png: svg('PCA — K-Means clusters', 5),
          train_csv: 'C:\\tmp\\train.csv', test_csv: 'C:\\tmp\\test.csv', cluster_table_csv: 'C:\\tmp\\members.csv',
          cluster_table: Array.from({ length: 240 }, (_, i) => ({
            Sample: `IR${1000 + i}`, Cluster: (i % 5) + 1, PC1: +(Math.sin(i) * 4).toFixed(3), PC2: +(Math.cos(i * 1.3) * 3).toFixed(3),
            Set: i % 5 === 0 ? 'Test' : 'Train',
          })),
        }))
      }
    }, 650)
  }
  bridge.run_gwas = run
  bridge.run_clustering = run
  return bridge
}

export function mockLauncherBridge() {
  return { open_tool: (n) => console.log('[mock] open_tool', n) }
}

/**
 * Generic browser-only mock: every slot call that ends in a callback answers
 * with a demo path; every other call is logged. Used for `npm run dev`.
 */
export function mockGeneric(signals = []) {
  const sig = {}
  ;['resultReady', 'statusUpdate', ...signals].forEach((n) => { sig[n] = signal() })
  return new Proxy(sig, {
    get(t, name) {
      if (name in t) return t[name]
      return (...args) => {
        const cb = args[args.length - 1]
        console.log('[mock bridge]', String(name), args.slice(0, -1))
        if (typeof cb === 'function') setTimeout(() => cb(/browse/.test(String(name)) ? 'C:\\data\\demo_input.csv' : ''), 250)
        else if (/^run/.test(String(name))) setTimeout(() => sig.resultReady.emit(JSON.stringify({ status: 'error', message: 'Preview mode: no Python backend.' })), 700)
      }
    },
  })
}

