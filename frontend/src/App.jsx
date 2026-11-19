import { useCallback, useEffect, useRef, useState } from 'react'
import { Sidebar, TopBar, RunToasts } from './components/Shell.jsx'
import { Background } from './components/kit.jsx'
import { RouteContext } from './lib/runState.js'
import { useBridge } from './lib/bridge.js'
import { setDiagSender, watchLongTasks, motionGuard } from './lib/perf.js'
import Dashboard from './pages/Dashboard.jsx'
import GenofileConverterTool from './pages/GenofileConverterTool.jsx'
import RecombinationTool from './pages/RecombinationTool.jsx'
import MendelianTool from './pages/MendelianTool.jsx'
import LinkageTool from './pages/LinkageTool.jsx'
import GenotypeViewerTool from './pages/GenotypeViewerTool.jsx'
import QtlTool from './pages/QtlTool.jsx'
import GeneExtractorTool from './pages/GeneExtractorTool.jsx'
import GenePrioritizationTool from './pages/GenePrioritizationTool.jsx'
import PybsaseqTool from './pages/PybsaseqTool.jsx'
import QtlSeqTool from './pages/QtlSeqTool.jsx'
import InSilicoBsaTool from './pages/InSilicoBsaTool.jsx'
import DiversityTool from './pages/DiversityTool.jsx'
import GwasTool from './pages/GwasTool.jsx'
import EtGwasTool from './pages/EtGwasTool.jsx'
import ClusteringTool from './pages/ClusteringTool.jsx'
import AdditiveDominanceTool from './pages/AdditiveDominanceTool.jsx'
import GenomicPredictionTool from './pages/GenomicPredictionTool.jsx'
import CrossSelectionTool from './pages/CrossSelectionTool.jsx'
import HeterosisTool from './pages/HeterosisTool.jsx'
import PhenotypicTool from './pages/PhenotypicTool.jsx'
import FavourableAlleleTool from './pages/FavourableAlleleTool.jsx'
import HelpPage from './pages/HelpPage.jsx'

// Hash route -> page. Keys match the tool ids in lib/tools.js.
const ROUTES = {
  '': Dashboard,
  genofile: GenofileConverterTool,
  mendelian: MendelianTool,
  recombination: RecombinationTool,
  linkage: LinkageTool,
  genotype_viewer: GenotypeViewerTool,
  qtl: QtlTool,
  gene_extractor: GeneExtractorTool,
  gene_prioritization: GenePrioritizationTool,
  pybsaseq: PybsaseqTool,
  qtl_seq: QtlSeqTool,
  insilico_bsa: InSilicoBsaTool,
  diversity: DiversityTool,
  gwas: GwasTool,
  et_gwas: EtGwasTool,
  clustering: ClusteringTool,
  additive_dominance: AdditiveDominanceTool,
  genomic_prediction: GenomicPredictionTool,
  cross_selection: CrossSelectionTool,
  heterosis: HeterosisTool,
  phenotypic: PhenotypicTool,
  favourable_allele: FavourableAlleleTool,
  help: HelpPage,
}

const parse = () => {
  const r = window.location.hash.replace(/^#\/?/, '')
  return r in ROUTES ? r : ''
}

let collapsedPref = false
try { collapsedPref = localStorage.getItem('tt.sidebar') === '1' } catch { /* storage unavailable */ }

export default function App() {
  const [route, setRoute] = useState(parse)
  // Diagnostics: page reports (long tasks, slow runs) go to the Python log; the active tool is named in stall reports.
  const { bridge: tt } = useBridge('terraTools')
  useEffect(() => { watchLongTasks(); motionGuard() }, [])
  useEffect(() => { if (tt && tt.diag_log) setDiagSender((s) => tt.diag_log(s)) }, [tt])
  useEffect(() => { if (tt && tt.set_active_tool) tt.set_active_tool(route || 'home') }, [tt, route])
  // Tool pages stay mounted once opened, so a running analysis, its inputs
  // and its results survive going Back and coming in again.
  const [opened, setOpened] = useState(() => [...new Set(['', parse()])])
  const [collapsed, setCollapsedState] = useState(collapsedPref)
  // Our own back stack. Hash changes use location.replace, so the web view's
  // history never grows and nothing else can "go back" behind our back.
  const stack = useRef([])

  useEffect(() => {
    const on = () => {
      const r = parse()
      setRoute(r)
      setOpened((o) => (o.includes(r) ? o : [...o, r]))
    }
    window.addEventListener('hashchange', on)
    return () => window.removeEventListener('hashchange', on)
  }, [])

  const navigate = useCallback((r) => {
    const cur = parse()
    if (r === cur) return
    stack.current.push(cur)
    window.location.replace('#/' + r)
  }, [])

  const back = useCallback(() => {
    const cur = parse()
    if (!cur) return
    let prev = stack.current.pop()
    while (prev === cur) prev = stack.current.pop()
    window.location.replace('#/' + (prev ?? ''))
  }, [])

  // Alt+Left, and the mouse "back" side button.
  useEffect(() => {
    const key = (e) => { if (e.altKey && e.key === 'ArrowLeft') { e.preventDefault(); back() } }
    const mouse = (e) => { if (e.button === 3) { e.preventDefault(); back() } }
    window.addEventListener('keydown', key)
    window.addEventListener('mouseup', mouse)
    return () => { window.removeEventListener('keydown', key); window.removeEventListener('mouseup', mouse) }
  }, [back])

  const setCollapsed = (v) => {
    setCollapsedState(v)
    try { localStorage.setItem('tt.sidebar', v ? '1' : '0') } catch { /* storage unavailable */ }
  }

  return (
    <div className="app-shell">
      <Sidebar route={route} navigate={navigate} collapsed={collapsed} setCollapsed={setCollapsed} />
      <div className="app-main">
        <Background variant="mint" count={8} />
        <TopBar route={route} back={back} navigate={navigate} />
        <div className="views">
          {opened.map((r) => {
            const Page = ROUTES[r]
            const active = r === route
            return (
              <RouteContext.Provider key={r || 'home'} value={r}>
                <section className={`view ${active ? 'active' : ''}`} inert={!active} aria-hidden={!active}>
                  <Page navigate={navigate} />
                </section>
              </RouteContext.Provider>
            )
          })}
        </div>
      </div>
      <RunToasts route={route} navigate={navigate} />
    </div>
  )
}
