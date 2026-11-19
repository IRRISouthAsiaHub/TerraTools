import { useMemo, useState } from 'react'
import { AnimatePresence, motion } from 'framer-motion'
import { Layers, CircleCheck, Network, Microscope, Search, X, ArrowUpRight, Leaf, LoaderCircle } from 'lucide-react'
import { Tabs, StatCard } from '../components/kit.jsx'
import { TOOLS, CATEGORIES, ALL_TAB } from '../lib/tools.js'
import { useRunStates } from '../lib/runState.js'

const ease = [0.22, 1, 0.36, 1]

// deterministic sparkle positions for the banner
const rnd = (i, k) => { const x = Math.sin(i * 12.9898 + k * 78.233) * 43758.5453; return x - Math.floor(x) }
const SPARKS = Array.from({ length: 14 }, (_, i) => ({ x: rnd(i, 1) * 100, y: 30 + rnd(i, 2) * 70, d: 7 + rnd(i, 3) * 8, delay: -rnd(i, 4) * 10, dy: -(40 + rnd(i, 5) * 90) }))

function AppCard({ tool, index, onOpen, run }) {
  const Icon = tool.icon
  // Entrance, hover lift and press are CSS (see .app-card): no per-mouse-move JavaScript and no inline transforms.
  return (
    <button
      className="app-card"
      onClick={() => onOpen(tool.id)}
      style={{ '--i': Math.min(index, 12) }}
    >
      <div className="top">
        <span className="category-label">{CATEGORIES.find((c) => c.id === tool.cat).label}</span>
        {run === 'running'
          ? <span className="status-pill running"><LoaderCircle size={12} className="spin" /> Running</span>
          : run === 'success' ? <span className="status-pill ok"><CircleCheck size={12} /> Results ready</span>
          : <span className="app-number">{tool.n}</span>}
      </div>
      <h3><span className="app-icon"><Icon size={19} /></span>{tool.name}</h3>
      <p>{tool.desc}</p>
      <div className="open-row">
        <span>{run === 'running' || run === 'success' ? 'Continue' : 'Open'}</span>
        <ArrowUpRight size={17} />
      </div>
    </button>
  )
}

export default function Dashboard({ navigate }) {
  const [tab, setTab] = useState('all')
  const [q, setQ] = useState('')
  const runs = useRunStates()

  const tabs = useMemo(
    () => [{ ...ALL_TAB, badge: TOOLS.length }, ...CATEGORIES.map((c) => ({ ...c, badge: TOOLS.filter((t) => t.cat === c.id).length }))],
    [],
  )
  const visible = useMemo(() => {
    const n = q.trim().toLowerCase()
    return TOOLS.filter((t) => (tab === 'all' || t.cat === tab) && (!n || t.name.toLowerCase().includes(n) || t.desc.toLowerCase().includes(n)))
  }, [tab, q])
  const groups = tab === 'all' ? CATEGORIES : CATEGORIES.filter((c) => c.id === tab)
  let idx = 0

  return (
    <div className="view-inner home">
      <motion.section className="banner" initial={{ opacity: 0, y: -12 }} animate={{ opacity: 1, y: 0 }} transition={{ duration: 0.6, ease }}>
        <div className="banner-photo" />
        <div className="banner-shade" />
        {SPARKS.map((s, i) => (
          <span key={i} className="firefly rise" style={{ left: `${s.x}%`, top: `${s.y}%`, '--dy': `${s.dy}px`, animationDuration: `${s.d}s`, animationDelay: `${s.delay}s` }} />
        ))}
        <div className="banner-content">
          <div className="logo-row">
            <h1>TerraTools</h1>
          </div>
          <p className="subtitle">
            A centralized platform for all IRRI-SAH data analysis applications — from trait mapping to genomic selection.
          </p>
        </div>
      </motion.section>

      <div className="stats">
        <StatCard icon={Layers} value={TOOLS.length} label="Total Applications" delay={0.1} />
        <StatCard icon={CircleCheck} value={TOOLS.length} label="Active Applications" delay={0.16} />
        <StatCard icon={Network} value={CATEGORIES.length} label="Tool Categories" delay={0.22} />
        <StatCard icon={Microscope} value="IRRI" label="International Rice Research Institute" delay={0.28} />
      </div>

      <motion.div className="panel on-mint" initial={{ opacity: 0, y: 20 }} animate={{ opacity: 1, y: 0 }} transition={{ duration: 0.6, delay: 0.15, ease }}>
        <Tabs tabs={tabs} value={tab} onChange={setTab} idPrefix="dash" />
        <div className="tab-content">
          <div className="toolbar">
            <div className="section-header" style={{ marginBottom: 0 }}>
              <h2>{tab === 'all' ? 'All Applications' : `${CATEGORIES.find((c) => c.id === tab).label}`}</h2>
              <p>Choose an application — it opens right here. Use Back or the sidebar to move between tools.</p>
            </div>
            <label className="search">
              <Search size={16} />
              <input value={q} onChange={(e) => setQ(e.target.value)} placeholder="Search applications…" />
              {q && <button onClick={() => setQ('')} aria-label="Clear search"><X size={15} /></button>}
            </label>
          </div>

          <div key={tab + '|' + q}>
            {groups.map((g) => {
              const items = visible.filter((t) => t.cat === g.id)
              if (!items.length) return null
              return (
                <section key={g.id}>
                  {tab === 'all' && <div className="group-title"><h3><g.icon size={18} />{g.label}</h3></div>}
                  <div className="app-grid" style={tab !== 'all' ? { marginTop: 24 } : undefined}>
                    {items.map((t) => <AppCard key={t.id} tool={t} index={idx++} onOpen={navigate} run={runs[t.id]} />)}
                  </div>
                </section>
              )
            })}
          </div>

          <AnimatePresence>
            {!visible.length && (
              <motion.div className="no-results" initial={{ opacity: 0 }} animate={{ opacity: 1 }} exit={{ opacity: 0 }}>
                <Leaf size={40} strokeWidth={1.5} />
                <div>No applications match “{q}”.</div>
              </motion.div>
            )}
          </AnimatePresence>
        </div>
      </motion.div>

      <footer className="app-footer">
        contact: p.sinha@cgiar.org · International Rice Research Institute
      </footer>
    </div>
  )
}
