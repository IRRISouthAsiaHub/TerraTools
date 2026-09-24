import { useEffect, useState } from 'react'
import { Network, SlidersHorizontal, LayoutGrid, FileText, Group, ChartColumn, Grid3x3, Leaf, GitBranch, Table2, Info } from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, FilePicker, Note, SubTabs, Frame,
  EmptyState, DownloadButtons, DataTable, HtmlBox,
} from '../components/kit.jsx'
import { useBridge, useTool, mockGeneric, loadJson } from '../lib/bridge.js'

const TABS = [
  { id: 'input', label: 'Input & Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]
const SUB = [
  { id: 'block', label: 'Block Analysis', icon: ChartColumn },
  { id: 'heatmap', label: 'Metrics Heatmap', icon: Grid3x3 },
  { id: 'variety', label: 'Variety Comparison', icon: Leaf },
  { id: 'pca', label: 'PCA & Clustering', icon: Network },
  { id: 'tree', label: 'NJ Tree', icon: GitBranch },
  { id: 'table', label: 'Results Table', icon: Table2 },
]

const NA = () => <div className="empty-state" style={{ padding: 30 }}><Info size={30} /><div>Not applicable for this dataset.</div></div>
const PlotOr = ({ src, height }) => (src ? <Frame src={src} height={height} /> : <div className="frame-wrap"><NA /></div>)

export default function DiversityTool() {
  const { bridge } = useBridge('diversityBridge', () => mockGeneric())
  const tool = useTool(bridge, 'run_diversity')
  const [tab, setTab] = useState('input')
  const [sub, setSub] = useState('block')
  const [file, setFile] = useState('')
  const [k, setK] = useState(5)
  const [top, setTop] = useState(25)
  const [rows, setRows] = useState(null)
  tool.onResult(() => { setTab('results'); setSub('block') })

  const r = tool.result
  useEffect(() => {
    setRows(null)
    if (r && r.results_json) loadJson(r.results_json).then(setRows).catch((e) => console.error('[DIVERSITY] table load failed', e))
  }, [r])

  const run = () => {
    if (!file) return tool.fail('Please select an input file.')
    tool.run({ input_file: file, n_clusters: k, top_blocks: top })
  }
  const p = (r && r.plots) || {}

  return (
    <ToolPage
      id="div" icon={Network} tabs={TABS} tab={tab} setTab={setTab}
      title="Diversity Analysis"
      subtitle="Haplotypic diversity metrics, clustering, PCA and neighbour-joining tree across varieties"
      run={{ label: 'Run Full Analysis', runningLabel: 'Analysing…', state: tool.state, status: tool.status, onRun: run }}
    >
      {tab === 'input' && (
        <TabPane id="input" key="input">
          <SectionHeader title="Input Data & Parameters">Upload your haplotype matrix and configure the analysis.</SectionHeader>
          <Card title="Haplotype Matrix" icon={FileText}>
            <FilePicker value={file} onBrowse={async () => { const x = await tool.browse(); if (x) setFile(x) }}
              placeholder="Select CSV (Column 1: Variety, Column 2: Sample ID, Columns 3+: haplotype blocks)" />
            <Note>Column 1: Variety/Cultivar name • Column 2: Sample ID • Columns 3+: haplotype block values (numeric)</Note>
          </Card>
          <div className="row">
            <Card title="Number of Clusters" icon={Group} delay={0.05}>
              <Slider label="Groups for hierarchical clustering" value={k} onChange={setK} min={2} max={20} />
            </Card>
            <Card title="Top Variable Blocks" icon={ChartColumn} delay={0.1}>
              <Slider label="Most informative blocks to display" value={top} onChange={setTop} min={5} max={50} step={5} />
            </Card>
          </div>
        </TabPane>
      )}
      {tab === 'results' && (
        <TabPane id="results" key="results">
          <SectionHeader title="Results">Diversity visualizations, summary statistics, and downloadable results table.</SectionHeader>
          {!r ? <EmptyState icon={Network}>Run the analysis to see results here.</EmptyState> : (
            <>
              <HtmlBox html={r.summary_html} />
              <div style={{ height: 18 }} />
              <DownloadButtons onSave={tool.save} items={[{ icon: FileText, label: 'Results Table (CSV)', path: r.results_csv, name: 'diversity_results.csv' }]} />
              <SubTabs tabs={SUB} value={sub} onChange={setSub} id="divsub" />
              {sub === 'block' && <PlotOr src={p.block_variance} height={560} />}
              {sub === 'heatmap' && <PlotOr src={p.heatmap} height={560} />}
              {sub === 'variety' && (
                <div className="plot-grid-2"><PlotOr src={p.violin} height={540} /><PlotOr src={p.radar} height={620} /></div>
              )}
              {sub === 'pca' && (
                <div className="plot-grid-2"><PlotOr src={p.pca} height={540} /><PlotOr src={p.clusters_pca} height={540} /></div>
              )}
              {sub === 'tree' && <PlotOr src={p.nj_tree} height={750} />}
              {sub === 'table' && (rows ? <DataTable rows={rows} /> : <div className="empty-state"><span className="spinner" /></div>)}
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
