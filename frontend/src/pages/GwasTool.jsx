import { useState } from 'react'
import { AnimatePresence } from 'framer-motion'
import {
  Globe, Database, SlidersHorizontal, ChartArea, LayoutGrid, FileText, Dna, Ruler,
  Scale, Waypoints, Link2, FlaskConical, FileSpreadsheet, FileArchive,
} from 'lucide-react'
import {
  Background, ToolHeader, Tabs, TabPane, SectionHeader, Card, Slider, Toggle, ChipGroup,
  FilePicker, RunBar, EmptyState, ImageGallery, DownloadButtons,
} from '../components/kit.jsx'
import { useBridge, useTool, mockToolBridge } from '../lib/bridge.js'

const MODELS = ['OLS', 'MLM', 'EMMAX', 'FarmCPU', 'GEMMA', 'BLINK', 'mrMLM', 'FASTmrMLM']

const TABS = [
  { id: 'data', label: 'Data & Models', icon: Database },
  { id: 'advanced', label: 'Advanced', icon: SlidersHorizontal },
  { id: 'plots', label: 'Plot Settings', icon: ChartArea },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]

const DEFAULTS = {
  models: ['OLS', 'MLM'],
  maf: 0.05, permutations: 100, top_n: 15, minimum_pve: 0,
  mlm_pcs: 3, gemma_pcs: 3, farmcpu_iterations: 10, farmcpu_bin_size: 10, blink_ld_r2: 0.7,
  show_bonferroni: true, show_suggestive: true, show_lod: true, show_pve: false,
}

export default function GwasTool() {
  const { bridge } = useBridge('gwasBridge', () => mockToolBridge('gwas'))
  const tool = useTool(bridge, 'run_gwas')
  const [tab, setTab] = useState('data')
  const [file, setFile] = useState('')
  const [c, setC] = useState(DEFAULTS)
  const set = (k) => (v) => setC((p) => ({ ...p, [k]: v }))

  tool.onResult(() => setTab('results'))

  const browse = async () => { const p = await tool.browse(); if (p) setFile(p) }
  const run = () => {
    if (!file) return tool.fail('Please select an input file.')
    if (!c.models.length) return tool.fail('Please select at least one model.')
    tool.run({ input_file: file, ...c })
  }

  const r = tool.result
  const plots = r
    ? Object.entries(r.plots).map(([name, src]) => ({ src, label: name.replace(/^\d+_/, '').replace(/_/g, ' ') }))
    : []

  return (
    <>
      <Background variant="mint" count={10} />
      <div className="page">
        <div className="container narrow">
          <ToolHeader
            icon={Globe}
            title="GWAS Analysis"
            subtitle="Genome-wide association study — genotype–phenotype analysis across multiple models"
          />
          <div className="panel on-mint">
            <Tabs tabs={TABS} value={tab} onChange={setTab} idPrefix="gwas" />
            <div className="tab-content">
              <AnimatePresence mode="wait">
                {tab === 'data' && (
                  <TabPane id="data" key="data">
                    <SectionHeader title="Input Data & Models">
                      Select your genotype/phenotype file and choose which GWAS models to run.
                    </SectionHeader>
                    <Card title="Input Data" icon={FileText}>
                      <FilePicker value={file} onBrowse={browse} placeholder="Select genotype/phenotype CSV, TSV or TXT file" />
                    </Card>
                    <Card title="GWAS Models" icon={Dna} delay={0.05}>
                      <ChipGroup options={MODELS} value={c.models} onChange={set('models')} />
                    </Card>
                    <Card title="Basic Parameters" icon={Ruler} delay={0.1}>
                      <div className="row" style={{ gap: 36 }}>
                        <div style={{ flex: 1, minWidth: 260 }}>
                          <Slider label="Minimum MAF" value={c.maf} onChange={set('maf')} min={0} max={0.5} step={0.01} format={(v) => v.toFixed(2)} />
                          <Slider label="Permutations" value={c.permutations} onChange={set('permutations')} min={20} max={200} step={10} />
                        </div>
                        <div style={{ flex: 1, minWidth: 260 }}>
                          <Slider label="Top N" value={c.top_n} onChange={set('top_n')} min={1} max={100} />
                          <Slider label="Minimum PVE (%)" value={c.minimum_pve} onChange={set('minimum_pve')} min={0} max={100} step={0.5} />
                        </div>
                      </div>
                    </Card>
                  </TabPane>
                )}
                {tab === 'advanced' && (
                  <TabPane id="advanced" key="advanced">
                    <SectionHeader title="Advanced Model Parameters">Fine-tune the mixed-model and iterative-model settings.</SectionHeader>
                    <div className="row">
                      <Card title="Mixed Model Parameters" icon={Scale}>
                        <Slider label="MLM PCs" value={c.mlm_pcs} onChange={set('mlm_pcs')} min={0} max={50} />
                        <Slider label="GEMMA PCs" value={c.gemma_pcs} onChange={set('gemma_pcs')} min={0} max={50} />
                      </Card>
                      <Card title="FarmCPU" icon={Waypoints} delay={0.06}>
                        <Slider label="Iterations" value={c.farmcpu_iterations} onChange={set('farmcpu_iterations')} min={1} max={100} />
                        <Slider label="Bin size" value={c.farmcpu_bin_size} onChange={set('farmcpu_bin_size')} min={1} max={1000} />
                      </Card>
                      <Card title="BLINK" icon={Link2} delay={0.12}>
                        <Slider label="LD r² threshold" value={c.blink_ld_r2} onChange={set('blink_ld_r2')} min={0} max={1} step={0.01} format={(v) => v.toFixed(2)} />
                      </Card>
                    </div>
                  </TabPane>
                )}
                {tab === 'plots' && (
                  <TabPane id="plots" key="plots">
                    <SectionHeader title="Plot Options">Control which reference lines and filters appear on the Manhattan plot.</SectionHeader>
                    <Card>
                      <Toggle label="Show Bonferroni threshold" checked={c.show_bonferroni} onChange={set('show_bonferroni')} />
                      <Toggle label="Show suggestive threshold" checked={c.show_suggestive} onChange={set('show_suggestive')} />
                      <Toggle label="Show permutation/LOD threshold" checked={c.show_lod} onChange={set('show_lod')} />
                      <Toggle label="Apply PVE filtering" checked={c.show_pve} onChange={set('show_pve')} />
                    </Card>
                  </TabPane>
                )}
                {tab === 'results' && (
                  <TabPane id="results" key="results">
                    <SectionHeader title="Results">Generated plots, downloadable Excel report, and plot archive.</SectionHeader>
                    {!r ? (
                      <EmptyState icon={FlaskConical}>Run the pipeline to see results here.</EmptyState>
                    ) : (
                      <>
                        <ImageGallery images={plots} />
                        <DownloadButtons onSave={tool.save} items={[
                          { icon: FileSpreadsheet, label: 'Excel Report', path: r.excel, name: 'gwas_results.xlsx' },
                          { icon: FileArchive, label: 'All Plots (ZIP)', path: r.zip, name: 'gwas_plots.zip' },
                        ]} />
                        {r.summary && <pre className="summary-box">{r.summary}</pre>}
                      </>
                    )}
                  </TabPane>
                )}
              </AnimatePresence>
            </div>
            <RunBar label="Run GWAS" runningLabel="Running GWAS…" state={tool.state} status={tool.status} onRun={run} />
          </div>
        </div>
      </div>
    </>
  )
}
