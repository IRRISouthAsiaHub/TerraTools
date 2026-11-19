import { useContext, useEffect, useId, useMemo, useRef, useState } from 'react'
import { createPortal } from 'react-dom'
import { AnimatePresence, motion, animate } from 'framer-motion'
import {
  FolderOpen, FileText, Play, LoaderCircle, CircleCheck, TriangleAlert,
  Search, X, Maximize2, Sprout, Download, Info,
} from 'lucide-react'
import { toFileUrl } from '../lib/bridge.js'
import { RouteContext } from '../lib/runState.js'
import { useFiles } from '../lib/files.js'
import { useRegisterImages, useRouteImages } from '../lib/imageRegistry.js'

/* ───────────── Ambient background ───────────── */

// deterministic pseudo-random so fireflies don't jump on re-render
const rnd = (i, k) => {
  const x = Math.sin(i * 12.9898 + k * 78.233) * 43758.5453
  return x - Math.floor(x)
}

// Fireflies are plain CSS animations (transform + opacity only), so Chromium runs them on the compositor thread and
// they cost the page's main thread nothing. They used to be framer-motion animations: JavaScript updating the style
// of every firefly on every frame, forever, which kept the renderer busy (a style recalculation per frame) even when
// the window was idle.
export function Background({ variant = 'photo', count = 16 }) {
  const flies = useMemo(
    () => Array.from({ length: count }, (_, i) => ({
      x: rnd(i, 1) * 100, y: 40 + rnd(i, 2) * 60,
      dx: (rnd(i, 3) - 0.5) * 160, dy: -(120 + rnd(i, 4) * 320),
      d: 9 + rnd(i, 5) * 12, delay: -rnd(i, 6) * 14, s: 0.6 + rnd(i, 7) * 0.9,
    })),
    [count],
  )
  return (
    <div className={`bg ${variant === 'photo' ? '' : 'bg-mint'}`} aria-hidden>
      {flies.map((f, i) => (
        <span
          key={i}
          className="firefly drift"
          style={{ left: `${f.x}%`, top: `${f.y}%`, '--dx': `${f.dx}px`, '--dy': `${f.dy}px`, '--s': f.s, animationDuration: `${f.d}s`, animationDelay: `${f.delay}s` }}
        />
      ))}
    </div>
  )
}

/* ───────────── Layout ───────────── */

export function LogoIcon({ icon: Icon = Sprout, size = 26 }) {
  return <div className="logo-icon"><Icon size={size} strokeWidth={2} /></div>
}

export function ToolHeader({ icon, title, subtitle }) {
  return (
    <motion.header
      className="tool-header"
      initial={{ opacity: 0, y: -14 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.5, ease: [0.22, 1, 0.36, 1] }}
    >
      <LogoIcon icon={icon} size={22} />
      <div>
        <h1>{title}</h1>
        <p>{subtitle}</p>
      </div>
    </motion.header>
  )
}

/** Tabs with a sliding highlight, animated pane switching. */
export function Tabs({ tabs, value, onChange, idPrefix = 'tabs' }) {
  return (
    <div className="tabs" role="tablist">
      {tabs.map((t) => {
        const active = t.id === value
        const Icon = t.icon
        return (
          <button
            key={t.id}
            role="tab"
            aria-selected={active}
            className={`tab-btn ${active ? 'active' : ''}`}
            onClick={() => onChange(t.id)}
          >
            {active && <motion.span layoutId={`${idPrefix}-bg`} className="tab-bg" transition={{ type: 'spring', stiffness: 500, damping: 38 }} />}
            {active && <motion.span layoutId={`${idPrefix}-ul`} className="tab-underline" transition={{ type: 'spring', stiffness: 500, damping: 38 }} />}
            {Icon && <Icon size={16} />}
            <span>{t.label}</span>
            {t.badge != null && <span className="tab-badge">{t.badge}</span>}
          </button>
        )
      })}
    </div>
  )
}

export function TabPane({ id, children }) {
  return (
    <motion.div
      key={id}
      initial={{ opacity: 0, y: 10 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.3, ease: [0.22, 1, 0.36, 1] }}
    >
      {children}
    </motion.div>
  )
}

export function SectionHeader({ title, children }) {
  return (
    <div className="section-header">
      <h2>{title}</h2>
      {children && <p>{children}</p>}
    </div>
  )
}

export function Card({ title, icon: Icon, children, className = '', delay = 0 }) {
  return (
    <motion.div
      className={`card ${className}`}
      initial={{ opacity: 0, y: 16 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.45, delay, ease: [0.22, 1, 0.36, 1] }}
    >
      {title && <h3>{Icon && <Icon size={17} />}{title}</h3>}
      {children}
    </motion.div>
  )
}

/* ───────────── Inputs ───────────── */

export function Slider({ label, value, onChange, min, max, step = 1, format }) {
  const pct = ((value - min) / (max - min)) * 100
  return (
    <>
      <div className="field-label">
        <span>{label}</span>
        <b>{format ? format(value) : value}</b>
      </div>
      <input
        type="range" className="slider" min={min} max={max} step={step} value={value}
        style={{ '--pct': `${pct}%` }}
        onChange={(e) => onChange(Number(e.target.value))}
        aria-label={label}
      />
    </>
  )
}

export function Toggle({ label, checked, onChange }) {
  return (
    <div className="toggle-row">
      <span>{label}</span>
      <label className="switch">
        <input type="checkbox" checked={checked} onChange={(e) => onChange(e.target.checked)} />
        <span className="track" /><span className="knob" />
      </label>
    </div>
  )
}

export function ChipGroup({ options, value, onChange }) {
  const toggle = (o) => onChange(value.includes(o) ? value.filter((v) => v !== o) : [...value, o])
  return (
    <div className="chip-grid">
      {options.map((o) => {
        const checked = value.includes(o)
        return (
          <motion.label key={o} className={`chip ${checked ? 'checked' : ''}`} whileTap={{ scale: 0.94 }}>
            <input type="checkbox" checked={checked} onChange={() => toggle(o)} />
            <AnimatePresence initial={false}>
              {checked && (
                <motion.span
                  initial={{ width: 0, opacity: 0 }} animate={{ width: 14, opacity: 1 }} exit={{ width: 0, opacity: 0 }}
                  style={{ display: 'inline-flex', overflow: 'hidden' }}
                >
                  <CircleCheck size={14} />
                </motion.span>
              )}
            </AnimatePresence>
            {o}
          </motion.label>
        )
      })}
    </div>
  )
}

export function FilePicker({ value, onBrowse, placeholder }) {
  return (
    <div className="file-drop">
      <div className={`file-path ${value ? 'filled' : ''}`} title={value}>
        <FileText size={17} />
        <span>{value ? `‎${value}` : placeholder}</span>
      </div>
      <button className="btn btn-secondary" onClick={onBrowse}><FolderOpen size={16} /> Browse</button>
    </div>
  )
}

export function TextField({ label, value, onChange, placeholder, hint, type = 'text', min, max, step, width }) {
  return (
    <div className="field" style={width ? { width } : undefined}>
      {label && <label className="field-title">{label}</label>}
      <input
        className="input" type={type} value={value} placeholder={placeholder} min={min} max={max} step={step}
        onChange={(e) => onChange(type === 'number' ? (e.target.value === '' ? '' : Number(e.target.value)) : e.target.value)}
      />
      {hint && <div className="field-hint">{hint}</div>}
    </div>
  )
}

export function SelectField({ label, value, onChange, options, hint }) {
  return (
    <div className="field">
      {label && <label className="field-title">{label}</label>}
      <select className="input" value={value} onChange={(e) => onChange(e.target.value)}>
        {options.map((o) => {
          const v = typeof o === 'object' ? o.value : o
          const l = typeof o === 'object' ? o.label : o
          return <option key={v} value={v}>{l}</option>
        })}
      </select>
      {hint && <div className="field-hint">{hint}</div>}
    </div>
  )
}

/** Single-choice pills (radio behaviour). */
export function RadioChips({ options, value, onChange }) {
  return (
    <div className="chip-grid">
      {options.map((o) => {
        const v = typeof o === 'object' ? o.value : o
        const l = typeof o === 'object' ? o.label : o
        const checked = v === value
        return (
          <motion.label key={v} className={`chip ${checked ? 'checked' : ''}`} whileTap={{ scale: 0.94 }}>
            <input type="radio" checked={checked} onChange={() => onChange(v)} />
            {l}
          </motion.label>
        )
      })}
    </div>
  )
}

export function ColorField({ label, value, onChange }) {
  return (
    <div className="field">
      {label && <label className="field-title">{label}</label>}
      <div className="color-row">
        <input type="color" value={/^#[0-9a-f]{6}$/i.test(value) ? value : '#000000'} onChange={(e) => onChange(e.target.value)} />
        <input className="input" value={value} onChange={(e) => onChange(e.target.value)} />
      </div>
    </div>
  )
}

export function Note({ children }) {
  return <div className="note">{children}</div>
}

/** Secondary tab row used inside a results pane. */
export function SubTabs({ tabs, value, onChange, id = 'sub' }) {
  return (
    <div className="subtabs">
      {tabs.map((t) => {
        const Icon = t.icon
        const active = t.id === value
        return (
          <button key={t.id} className={`subtab-btn ${active ? 'active' : ''}`} onClick={() => onChange(t.id)}>
            {active && <motion.span layoutId={`${id}-pill`} className="subtab-pill" transition={{ type: 'spring', stiffness: 520, damping: 40 }} />}
            {Icon && <Icon size={15} />}<span>{t.label}</span>
          </button>
        )
      })}
    </div>
  )
}

/** Embedded local HTML/plotly file. */
export function Frame({ src, height = 560, title, downloadName }) {
  const [loaded, setLoaded] = useState(false)
  const ref = useRef(null)
  // Interactive plotly pages listen for this message and download their own PNG
  // (the system Save dialog then opens).
  const downloadPng = () => ref.current && ref.current.contentWindow &&
    ref.current.contentWindow.postMessage({ type: 'tt-download-png', name: downloadName }, '*')
  return (
    <>
      <div className="frame-wrap" style={{ height }}>
        {!loaded && <div className="frame-loading"><span className="spinner" /></div>}
        <iframe ref={ref} title={title || 'plot'} src={toFileUrl(src)} onLoad={() => setLoaded(true)} style={{ opacity: loaded ? 1 : 0 }} />
      </div>
      {downloadName && (
        <div className="frame-actions">
          <button className="btn btn-secondary sm" onClick={downloadPng} disabled={!loaded}><Download size={14} /> Download PNG</button>
        </div>
      )}
    </>
  )
}

/** Monospace log; sticks to the bottom as text is appended. */
export function LogBox({ children }) {
  const ref = useRef(null)
  useEffect(() => { if (ref.current) ref.current.scrollTop = ref.current.scrollHeight }, [children])
  return <pre ref={ref} className="summary-box log">{children}</pre>
}

/** Determinate progress bar, value 0..1. */
export function ProgressBar({ value }) {
  return (
    <div className="pbar">
      <motion.div className="pbar-fill" animate={{ width: `${Math.round(Math.max(0, Math.min(1, value)) * 100)}%` }} transition={{ duration: 0.4, ease: [0.22, 1, 0.36, 1] }} />
      <span>{Math.round(Math.max(0, Math.min(1, value)) * 100)}%</span>
    </div>
  )
}

/** File picker + optional one-line preview beneath it. */
export function FileWithPreview({ label, value, placeholder, onBrowse, preview }) {
  return (
    <div className="field">
      {label && <label className="field-title">{label}</label>}
      <FilePicker value={value} placeholder={placeholder} onBrowse={onBrowse} />
      {preview && <div className="field-hint">{preview}</div>}
    </div>
  )
}

/** Server-rendered (trusted, local) HTML fragment. */
export function HtmlBox({ html, className = 'summary-box html' }) {
  return <div className={className} dangerouslySetInnerHTML={{ __html: html }} />
}

const EXT_LABEL = { csv: 'CSV', pdf: 'PDF', txt: 'Text', svg: 'SVG', eps: 'EPS', html: 'HTML', png: 'PNG', xlsx: 'Excel', zip: 'ZIP', json: 'JSON' }

/** Listing of output files with per-file save buttons. */
export function FileList({ files, onSave }) {
  if (!files || !files.length) return <div className="empty-state" style={{ padding: 30 }}>No output files found.</div>
  return (
    <div className="file-list">
      {files.map((f, i) => (
        <motion.div
          key={f.path + i} className="file-row-item"
          initial={{ opacity: 0, x: -12 }} animate={{ opacity: 1, x: 0 }} transition={{ delay: Math.min(i, 15) * 0.03 }}
        >
          <span className="fname"><FileText size={16} />{f.name}<em>{EXT_LABEL[f.ext] || f.ext}</em></span>
          <button className="btn btn-secondary sm" onClick={() => onSave(f.path, f.name)}><Download size={14} /> Save</button>
        </motion.div>
      ))}
    </div>
  )
}

/** One plot from an {html, png} pair: the PNG is shown whole in a card that opens a pop-up when clicked; the interactive page is the fallback. */
export function PlotSlot({ plot, name, title, height = 520 }) {
  const html = plot && plot.html
  const png = plot && plot.png
  const label = title || name
  if (png) return <div className="plot-block"><ImageGallery single images={[{ src: png, label, pdf: plot.pdf }]} /></div>
  return (
    <div className="plot-block">
      {title && <div className="plot-title">{title}</div>}
      {html ? <Frame src={html} height={height} title={title || name} />
        : <div className="frame-wrap"><div className="empty-state" style={{ padding: 24 }}><Info size={28} /><div>Not available for this dataset.</div></div></div>}
    </div>
  )
}

export function TableSlot({ rows, title, maxHeight = 340 }) {
  return (
    <div className="plot-block">
      {title && <div className="plot-title">{title}</div>}
      {rows && rows.length ? <DataTable rows={rows} maxHeight={maxHeight} /> : <div className="frame-wrap"><div className="empty-state" style={{ padding: 24 }}>No data.</div></div>}
    </div>
  )
}

/** Standard tool window: header + tabs + panes + run bar. */
export function ToolPage({ icon, title, subtitle, tabs, tab, setTab, id, children, run }) {
  const route = useContext(RouteContext)
  const { galleries, images } = useRouteImages(route)
  const { saveImages } = useFiles()
  return (
    <div className="view-inner">
      <ToolHeader icon={icon} title={title} subtitle={subtitle} />
      <div className="panel on-mint">
        <Tabs tabs={tabs} value={tab} onChange={setTab} idPrefix={id} />
        {galleries > 1 && images.length > 1 && (
          <div className="imgbar">
            <span><Download size={14} /> {images.length} images on this tab</span>
            <button className="btn btn-secondary sm" onClick={() => saveImages(images)}>Download all as PNG</button>
          </div>
        )}
        <div className="tab-content">
          {children}
        </div>
        {run && <RunBar {...run} />}
      </div>
    </div>
  )
}

/* ───────────── Run bar ───────────── */

export function RunBar({ label, runningLabel = 'Running…', state, status, onRun, disabled }) {
  const running = state === 'running'
  return (
    <div className="run-bar">
      <motion.button
        className={`btn-run ${running ? 'running' : ''}`}
        onClick={onRun}
        disabled={running || disabled}
        whileTap={running ? undefined : { scale: 0.985 }}
        style={{ opacity: running ? 0.9 : 1 }}
      >
        {running ? <LoaderCircle size={18} style={{ animation: 'spin 1s linear infinite' }} /> : <Play size={18} fill="currentColor" />}
        <span>{running ? runningLabel : label}</span>
      </motion.button>
      <AnimatePresence>
        {running && (
          <motion.div className="progress" initial={{ opacity: 0, scaleY: 0 }} animate={{ opacity: 1, scaleY: 1 }} exit={{ opacity: 0 }}>
            <i />
          </motion.div>
        )}
      </AnimatePresence>
      <div className={`status-bar ${state}`} aria-live="polite">
        {state === 'success' && <CircleCheck size={16} />}
        {state === 'error' && <TriangleAlert size={16} />}
        <AnimatePresence mode="wait" initial={false}>
          <motion.span key={status} initial={{ opacity: 0, x: -6 }} animate={{ opacity: 1, x: 0 }} exit={{ opacity: 0 }} transition={{ duration: 0.18 }}>
            {status}
          </motion.span>
        </AnimatePresence>
      </div>
    </div>
  )
}

/* ───────────── Results ───────────── */

export function EmptyState({ icon: Icon, children }) {
  return (
    <div className="empty-state">
      <Icon size={44} strokeWidth={1.5} />
      <div>{children}</div>
    </div>
  )
}

function CountUp({ value }) {
  const ref = useRef(null)
  useEffect(() => {
    const num = Number(value)
    if (!Number.isFinite(num)) { if (ref.current) ref.current.textContent = value; return }
    const dec = String(value).includes('.') ? String(value).split('.')[1].length : 0
    const c = animate(0, num, {
      duration: 1.1, ease: [0.22, 1, 0.36, 1],
      onUpdate: (v) => { if (ref.current) ref.current.textContent = dec ? v.toFixed(dec) : Math.round(v).toLocaleString() },
    })
    return () => c.stop()
  }, [value])
  return <span ref={ref}>0</span>
}

export function StatCard({ icon: Icon, value, label, delay = 0, className = 'stat-card' }) {
  return (
    <motion.div
      className={className}
      initial={{ opacity: 0, y: 18 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.5, delay, ease: [0.22, 1, 0.36, 1] }}
      whileHover={{ y: -4 }}
    >
      {Icon && <Icon size={26} />}
      <h3 className={className === 'stat-card' ? '' : 'val'}><CountUp value={value} /></h3>
      <p className={className === 'stat-card' ? '' : 'lbl'}>{label}</p>
    </motion.div>
  )
}

export function SummaryGrid({ items }) {
  return (
    <div className="summary-grid">
      {items.map((it, i) => (
        <motion.div
          key={it.label} className="summary-card"
          initial={{ opacity: 0, y: 14, scale: 0.96 }} animate={{ opacity: 1, y: 0, scale: 1 }}
          transition={{ delay: i * 0.06, duration: 0.4, ease: [0.22, 1, 0.36, 1] }}
        >
          <div className="val">{typeof it.value === 'number' ? <CountUp value={it.value} /> : it.value}</div>
          <div className="lbl">{it.label}</div>
        </motion.div>
      ))}
    </div>
  )
}

export function DownloadButtons({ items, onSave }) {
  return (
    <div className="downloads">
      {items.map((d, i) => (
        <motion.button
          key={d.name} className="btn btn-secondary" onClick={() => onSave(d.path, d.name)}
          initial={{ opacity: 0, y: 8 }} animate={{ opacity: 1, y: 0 }} transition={{ delay: 0.1 + i * 0.06 }}
          whileHover={{ y: -2 }} whileTap={{ scale: 0.96 }}
        >
          <d.icon size={16} /> {d.label}
        </motion.button>
      ))}
    </div>
  )
}

/** Grid (or single) of result images. Every image has a Download PNG button; click opens it full size. */
export function ImageGallery({ images, single = false }) {
  const [open, setOpen] = useState(null)
  const route = useContext(RouteContext)
  const id = useId()
  const { saveImage, saveImages } = useFiles()
  useRegisterImages(route, id, images)
  // Result files are often rewritten under the same path on every run.
  const stamp = useMemo(() => Date.now(), [images])
  const url = (p) => { const u = toFileUrl(p); return /^file:/.test(u) ? `${u}?t=${stamp}` : u }
  useEffect(() => {
    if (!open) return
    const onKey = (e) => e.key === 'Escape' && setOpen(null)
    window.addEventListener('keydown', onKey)
    return () => window.removeEventListener('keydown', onKey)
  }, [open])

  return (
    <>
      {images.length > 1 && (
        <div className="gallery-bar">
          <span>{images.length} images</span>
          <button className="btn btn-secondary sm" onClick={() => saveImages(images)}><Download size={14} /> Download all as PNG</button>
        </div>
      )}
      <div className={single ? 'plot-single' : 'results-grid'}>
        {images.map((im, i) => (
          <motion.div
            key={im.label + i} className="result-card"
            initial={{ opacity: 0, y: 20, scale: 0.97 }} animate={{ opacity: 1, y: 0, scale: 1 }}
            transition={{ delay: i * 0.08, duration: 0.5, ease: [0.22, 1, 0.36, 1] }}
          >
            <button className="img-wrap" onClick={() => setOpen(im)} aria-label={`View ${im.label} full size`}>
              <img src={url(im.src)} alt={im.label} decoding="async" loading="lazy" />
            </button>
            <div className="label">
              <span className="lbl-text">{im.label}</span>
              <span className="card-actions">
                <button className="btn btn-secondary sm" onClick={() => saveImage(im.src, im.label)}><Download size={14} /> Download PNG</button>
                {im.pdf && <button className="btn btn-secondary sm" onClick={() => saveImage(im.pdf, `${im.label}.pdf`)} title="Vector PDF"><Download size={14} /> PDF</button>}
                <button className="icon-btn" onClick={() => setOpen(im)} title="View full size" aria-label="View full size"><Maximize2 size={15} /></button>
              </span>
            </div>
          </motion.div>
        ))}
      </div>
      {createPortal(
      <AnimatePresence>
        {open && (
          <motion.div
            className="lightbox" onClick={() => setOpen(null)}
            initial={{ opacity: 0 }} animate={{ opacity: 1 }} exit={{ opacity: 0 }}
          >
            <button className="close" aria-label="Close"><X size={20} /></button>
            <motion.figure
              onClick={(e) => e.stopPropagation()}
              initial={{ scale: 0.92, y: 20 }} animate={{ scale: 1, y: 0 }} exit={{ scale: 0.95 }}
              transition={{ type: 'spring', stiffness: 300, damping: 28 }}
            >
              <img src={url(open.src)} alt={open.label} decoding="async" />
              <figcaption>
                {open.label}
                <button className="btn btn-primary sm" onClick={() => saveImage(open.src, open.label)}><Download size={14} /> Download PNG</button>
                {open.pdf && <button className="btn btn-primary sm" onClick={() => saveImage(open.pdf, `${open.label}.pdf`)}><Download size={14} /> PDF</button>}
              </figcaption>
            </motion.figure>
          </motion.div>
        )}
      </AnimatePresence>,
        document.body,
      )}
    </>
  )
}

const PAGE = 100
const CLUSTER_COLORS = ['#52b788', '#2d6a4f', '#e9c46a', '#e76f51', '#8ab17d', '#457b9d', '#b5838d', '#6d597a', '#f4a261', '#264653']

const fmt = (v) => (v === null || v === undefined ? '' : typeof v === 'number' && !Number.isInteger(v) ? v.toFixed(4) : String(v))

export function DataTable({ rows, colorColumn, highlight, renderCell, columns, maxHeight = 400 }) {
  const [q, setQ] = useState('')
  const [shown, setShown] = useState(PAGE)
  const cols = useMemo(() => columns || (rows && rows.length ? Object.keys(rows[0]) : []), [rows, columns])
  const filtered = useMemo(() => {
    if (!q.trim()) return rows || []
    const n = q.toLowerCase()
    return (rows || []).filter((r) => cols.some((c) => String(r[c] ?? '').toLowerCase().includes(n)))
  }, [rows, cols, q])
  useEffect(() => setShown(PAGE), [q, rows])

  if (!rows || !rows.length) return <div className="empty-state" style={{ padding: 24 }}>No data.</div>
  return (
    <>
      <div className="table-tools">
        <label className="search">
          <Search size={15} />
          <input value={q} onChange={(e) => setQ(e.target.value)} placeholder="Filter rows…" />
          {q && <button onClick={() => setQ('')} aria-label="Clear"><X size={14} /></button>}
        </label>
        <small>{filtered.length.toLocaleString()} of {rows.length.toLocaleString()} rows</small>
      </div>
      <div className="table-wrap" style={{ maxHeight }}>
        <table className="results-table">
          <thead><tr>{cols.map((c) => <th key={c}>{c}</th>)}</tr></thead>
          <tbody>
            {filtered.slice(0, shown).map((r, i) => (
              <tr key={i}>
                {cols.map((c) => (
                  <td key={c} className={highlight && highlight(c, r) ? 'distorted' : undefined}>
                    {c === colorColumn && (
                      <span className="cluster-dot" style={{ background: CLUSTER_COLORS[(Number(r[c]) - 1) % CLUSTER_COLORS.length] || CLUSTER_COLORS[0] }} />
                    )}
                    {renderCell ? (renderCell(c, r[c], r) ?? fmt(r[c])) : fmt(r[c])}
                  </td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      {shown < filtered.length && (
        <div className="more-row">
          <button className="btn btn-secondary" onClick={() => setShown((s) => s + PAGE)}>
            Show {Math.min(PAGE, filtered.length - shown)} more
          </button>
        </div>
      )}
    </>
  )
}

