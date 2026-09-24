import { useEffect, useMemo, useRef, useState } from 'react'
import { AnimatePresence, motion } from 'framer-motion'
import { Sprout, House, ArrowLeft, ChevronRight, PanelLeftClose, PanelLeftOpen, Search, CircleCheck, TriangleAlert, X } from 'lucide-react'
import { TOOLS, CATEGORIES } from '../lib/tools.js'
import { useRunStates } from '../lib/runState.js'

const ease = [0.22, 1, 0.36, 1]
export const toolById = (id) => TOOLS.find((t) => t.id === id)

/* ───────────── Sidebar ───────────── */

function NavItem({ active, icon: Icon, label, onClick, collapsed, running, badge }) {
  return (
    <button className={`nav-item ${active ? 'active' : ''}`} onClick={onClick} title={collapsed ? label : undefined}>
      {active && <motion.span layoutId="nav-pill" className="nav-pill" transition={{ type: 'spring', stiffness: 520, damping: 42 }} />}
      <Icon size={18} className="nav-icon" />
      {!collapsed && <span className="nav-label">{label}</span>}
      {running && <span className="run-dot" title="Running" />}
      {!running && badge && <span className={`done-dot ${badge}`} />}
    </button>
  )
}

export function Sidebar({ route, navigate, collapsed, setCollapsed }) {
  const runs = useRunStates()
  return (
    <motion.aside className={`sidebar ${collapsed ? 'collapsed' : ''}`} animate={{ width: collapsed ? 72 : 264 }} transition={{ duration: 0.3, ease }}>
      <div className="brand" onClick={() => navigate('')}>
        <div className="brand-mark"><Sprout size={20} /></div>
        {!collapsed && <div><div className="brand-name">TerraTools</div><div className="brand-sub">IRRI · SAH</div></div>}
      </div>
      <nav className="nav-scroll">
        <NavItem active={route === ''} icon={House} label="Home" onClick={() => navigate('')} collapsed={collapsed} />
        {CATEGORIES.map((c) => (
          <div key={c.id} className="nav-group">
            {collapsed ? <div className="nav-sep" /> : <div className="nav-group-label">{c.label}</div>}
            {TOOLS.filter((t) => t.cat === c.id).map((t) => (
              <NavItem key={t.id} active={route === t.id} icon={t.icon} label={t.name} collapsed={collapsed}
                onClick={() => navigate(t.id)} running={runs[t.id] === 'running'}
                badge={route !== t.id && (runs[t.id] === 'success' ? 'ok' : runs[t.id] === 'error' ? 'err' : null)} />
            ))}
          </div>
        ))}
      </nav>
      <button className="collapse-btn" onClick={() => setCollapsed(!collapsed)} title={collapsed ? 'Expand sidebar' : 'Collapse sidebar'}>
        {collapsed ? <PanelLeftOpen size={18} /> : <><PanelLeftClose size={18} /><span>Collapse</span></>}
      </button>
    </motion.aside>
  )
}

/* ───────────── Top bar ───────────── */

function QuickOpen({ navigate }) {
  const [q, setQ] = useState('')
  const [open, setOpen] = useState(false)
  const [hi, setHi] = useState(0)
  const input = useRef(null)
  const hits = useMemo(() => {
    const n = q.trim().toLowerCase()
    return n ? TOOLS.filter((t) => t.name.toLowerCase().includes(n) || t.desc.toLowerCase().includes(n)).slice(0, 7) : TOOLS.slice(0, 7)
  }, [q])
  useEffect(() => {
    const onKey = (e) => { if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === 'k') { e.preventDefault(); input.current && input.current.focus() } }
    window.addEventListener('keydown', onKey)
    return () => window.removeEventListener('keydown', onKey)
  }, [])
  const go = (t) => { navigate(t.id); setQ(''); setOpen(false); input.current && input.current.blur() }
  return (
    <div className="quick">
      <label className="quick-input">
        <Search size={15} />
        <input ref={input} value={q} placeholder="Open a tool…" onFocus={() => { setOpen(true); setHi(0) }}
          onBlur={() => setTimeout(() => setOpen(false), 120)}
          onChange={(e) => { setQ(e.target.value); setHi(0) }}
          onKeyDown={(e) => {
            if (e.key === 'ArrowDown') { e.preventDefault(); setHi((h) => Math.min(h + 1, hits.length - 1)) }
            if (e.key === 'ArrowUp') { e.preventDefault(); setHi((h) => Math.max(h - 1, 0)) }
            if (e.key === 'Enter' && hits[hi]) go(hits[hi])
            if (e.key === 'Escape') input.current.blur()
          }} />
        <kbd>Ctrl K</kbd>
      </label>
      <AnimatePresence>
        {open && hits.length > 0 && (
          <motion.div className="quick-menu" initial={{ opacity: 0, y: -6 }} animate={{ opacity: 1, y: 0 }} exit={{ opacity: 0, y: -6 }} transition={{ duration: 0.15 }}>
            {hits.map((t, i) => (
              <button key={t.id} className={i === hi ? 'hi' : ''} onMouseDown={() => go(t)} onMouseEnter={() => setHi(i)}>
                <span className="app-icon sm"><t.icon size={15} /></span>
                <span><b>{t.name}</b><small>{CATEGORIES.find((c) => c.id === t.cat).label}</small></span>
              </button>
            ))}
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  )
}

export function TopBar({ route, back, navigate }) {
  const tool = toolById(route)
  const cat = tool && CATEGORIES.find((c) => c.id === tool.cat)
  return (
    <header className="topbar">
      <motion.button className="back-btn" onClick={back} disabled={!tool} whileTap={tool ? { scale: 0.94 } : undefined} title="Back (Alt+←)">
        <ArrowLeft size={17} /><span>Back</span>
      </motion.button>
      <nav className="crumbs">
        <button onClick={() => navigate('')}>Home</button>
        <AnimatePresence initial={false}>
          {tool && (
            <motion.span key={route} className="crumb-tail" initial={{ opacity: 0, x: -8 }} animate={{ opacity: 1, x: 0 }} exit={{ opacity: 0 }} transition={{ duration: 0.2 }}>
              <ChevronRight size={14} /><span className="muted">{cat.label}</span>
              <ChevronRight size={14} /><b>{tool.name}</b>
            </motion.span>
          )}
        </AnimatePresence>
      </nav>
      <QuickOpen navigate={navigate} />
    </header>
  )
}

/* ───────────── "Finished in the background" notifications ───────────── */

export function RunToasts({ route, navigate }) {
  const runs = useRunStates()
  const prev = useRef(runs)
  const [toasts, setToasts] = useState([])
  useEffect(() => {
    const fresh = Object.entries(runs)
      .filter(([r, s]) => prev.current[r] === 'running' && (s === 'success' || s === 'error') && r !== route)
      .map(([r, s]) => ({ id: `${r}-${Date.now()}`, route: r, ok: s === 'success' }))
    prev.current = runs
    if (!fresh.length) return
    setToasts((t) => [...t, ...fresh])
    fresh.forEach((f) => setTimeout(() => setToasts((t) => t.filter((x) => x.id !== f.id)), 8000))
  }, [runs, route])
  return (
    <div className="toasts">
      <AnimatePresence>
        {toasts.map((t) => {
          const tool = toolById(t.route)
          return (
            <motion.div key={t.id} className={`toast ${t.ok ? 'ok' : 'err'}`} layout
              initial={{ opacity: 0, x: 40 }} animate={{ opacity: 1, x: 0 }} exit={{ opacity: 0, x: 40 }} transition={{ duration: 0.3, ease }}>
              {t.ok ? <CircleCheck size={18} /> : <TriangleAlert size={18} />}
              <div onClick={() => { navigate(t.route); setToasts((x) => x.filter((y) => y.id !== t.id)) }}>
                <b>{tool ? tool.name : t.route}</b>
                <span>{t.ok ? 'Analysis finished — click to view results' : 'Analysis failed — click to see why'}</span>
              </div>
              <button onClick={() => setToasts((x) => x.filter((y) => y.id !== t.id))} aria-label="Dismiss"><X size={14} /></button>
            </motion.div>
          )
        })}
      </AnimatePresence>
    </div>
  )
}
