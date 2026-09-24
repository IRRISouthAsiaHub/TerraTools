import { useMemo, useState } from 'react'
import { AnimatePresence, LayoutGroup, motion } from 'framer-motion'
import { Sprout, Layers, CircleCheck, Network, Microscope, Search, X, ArrowUpRight, LoaderCircle, Leaf } from 'lucide-react'
import { Background, LogoIcon, Tabs, StatCard } from '../components/kit.jsx'
import { useBridge, mockLauncherBridge } from '../lib/bridge.js'
import { TOOLS, CATEGORIES, ALL_TAB } from '../lib/tools.js'

const ease = [0.22, 1, 0.36, 1]

function AppCard({ tool, index, launching, onOpen, showCategory }) {
  const Icon = tool.icon
  const track = (e) => {
    const r = e.currentTarget.getBoundingClientRect()
    e.currentTarget.style.setProperty('--mx', `${e.clientX - r.left}px`)
    e.currentTarget.style.setProperty('--my', `${e.clientY - r.top}px`)
  }
  return (
    <motion.button
      layout
      className={`app-card ${launching ? 'launching' : ''}`}
      onClick={() => onOpen(tool)}
      onMouseMove={track}
      initial={{ opacity: 0, y: 26, scale: 0.97 }}
      animate={{ opacity: 1, y: 0, scale: 1 }}
      exit={{ opacity: 0, scale: 0.94, transition: { duration: 0.18 } }}
      transition={{ duration: 0.5, delay: Math.min(index, 12) * 0.05, ease, layout: { duration: 0.4, ease } }}
      whileHover={{ y: -6 }}
      whileTap={{ scale: 0.985 }}
    >
      <div className="top">
        <span className="category-label">{CATEGORIES.find((c) => c.id === tool.cat).label}</span>
        <span className="app-number">{tool.n}</span>
      </div>
      <h3><span className="app-icon"><Icon size={19} /></span>{tool.name}</h3>
      <p>{tool.desc}</p>
      <div className="open-row">
        <span>{launching ? 'Launching…' : 'Open Application'}</span>
        {launching ? <LoaderCircle size={17} /> : <ArrowUpRight size={17} />}
      </div>
    </motion.button>
  )
}

export default function Dashboard() {
  const { bridge } = useBridge('terraTools', mockLauncherBridge)
  const [tab, setTab] = useState('all')
  const [q, setQ] = useState('')
  const [launching, setLaunching] = useState(null)

  const open = (tool) => {
    if (!bridge) return console.error('TerraTools bridge is not ready')
    setLaunching(tool.id)
    bridge.open_tool(tool.id)
    setTimeout(() => setLaunching(null), 1600)
  }

  const tabs = useMemo(
    () => [{ ...ALL_TAB, badge: TOOLS.length }, ...CATEGORIES.map((c) => ({ ...c, badge: TOOLS.filter((t) => t.cat === c.id).length }))],
    [],
  )

  const visible = useMemo(() => {
    const n = q.trim().toLowerCase()
    return TOOLS.filter(
      (t) => (tab === 'all' || t.cat === tab) && (!n || t.name.toLowerCase().includes(n) || t.desc.toLowerCase().includes(n)),
    )
  }, [tab, q])

  const groups = tab === 'all' ? CATEGORIES : CATEGORIES.filter((c) => c.id === tab)
  let idx = 0

  return (
    <>
      <Background variant="photo" count={20} />
      <div className="page">
        <div className="container">
          <motion.header className="hero" initial={{ opacity: 0, y: -20 }} animate={{ opacity: 1, y: 0 }} transition={{ duration: 0.7, ease }}>
            <div className="logo-row">
              <LogoIcon icon={Sprout} />
              <h1>TerraTools</h1>
            </div>
            <p className="subtitle">
              A centralized platform for accessing all IRRI-SAH data analysis applications — from trait mapping to genomic selection.
            </p>
          </motion.header>

          <motion.div className="panel" initial={{ opacity: 0, y: 30 }} animate={{ opacity: 1, y: 0 }} transition={{ duration: 0.7, delay: 0.1, ease }}>
            <LayoutGroup>
              <Tabs tabs={tabs} value={tab} onChange={setTab} idPrefix="dash" />
            </LayoutGroup>

            <div className="tab-content">
              <div className="stats">
                <StatCard icon={Layers} value={TOOLS.length} label="Total Applications" delay={0.2} />
                <StatCard icon={CircleCheck} value={TOOLS.length} label="Active Applications" delay={0.28} />
                <StatCard icon={Network} value={CATEGORIES.length} label="Tool Categories" delay={0.36} />
                <StatCard icon={Microscope} value="IRRI" label="International Rice Research Institute" delay={0.44} />
              </div>

              <div className="toolbar">
                <div className="section-header" style={{ marginBottom: 0 }}>
                  <h2>{tab === 'all' ? 'All Applications at a Glance' : `${CATEGORIES.find((c) => c.id === tab).label} Applications`}</h2>
                  <p>Pick an application to launch it in its own window.</p>
                </div>
                <label className="search">
                  <Search size={16} />
                  <input value={q} onChange={(e) => setQ(e.target.value)} placeholder="Search applications…" />
                  {q && <button onClick={() => setQ('')} aria-label="Clear search"><X size={15} /></button>}
                </label>
              </div>

              <LayoutGroup>
                {groups.map((g) => {
                  const items = visible.filter((t) => t.cat === g.id)
                  if (!items.length) return null
                  return (
                    <section key={g.id}>
                      {tab === 'all' && (
                        <motion.div layout className="group-title">
                          <h3><g.icon size={18} />{g.label}</h3>
                        </motion.div>
                      )}
                      <motion.div layout className="app-grid" style={tab !== 'all' ? { marginTop: 24 } : undefined}>
                        <AnimatePresence mode="popLayout">
                          {items.map((t) => (
                            <AppCard key={`${tab}-${t.id}`} tool={t} index={idx++} launching={launching === t.id} onOpen={open} />
                          ))}
                        </AnimatePresence>
                      </motion.div>
                    </section>
                  )
                })}
              </LayoutGroup>

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

          <footer className="site-footer">
            <p>contact: p.sinha@cgiar.org &nbsp;|&nbsp; International Rice Research Institute</p>
            <p>This platform provides centralized access to data analysis tools for agricultural and genomic research.</p>
          </footer>
        </div>
      </div>
    </>
  )
}
