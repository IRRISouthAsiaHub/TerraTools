import { useState } from 'react'
import { Dna, SlidersHorizontal, LayoutGrid, FileText, Sprout, Ruler, Eraser } from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, Toggle, FilePicker, RadioChips, TextField,
  EmptyState, ImageGallery, DownloadButtons, DataTable, HtmlBox,
} from '../components/kit.jsx'
import { useBridge, useTool, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'input', label: 'Input & Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]

export default function MendelianTool() {
  const { bridge } = useBridge('mendelianBridge', () => mockGeneric())
  const tool = useTool(bridge, 'run_mendelian')
  const [tab, setTab] = useState('input')
  const [file, setFile] = useState('')
  const [pop, setPop] = useState('RIL')
  const [pthr, setPthr] = useState(0.05)
  const [missing, setMissing] = useState('NA,.,-9')
  const [dropUnplaced, setDropUnplaced] = useState(true)
  tool.onResult(() => setTab('results'))

  const run = () => {
    if (!file) return tool.fail('Please select an input file.')
    tool.run({ input_file: file, pop, pthr, missing_codes: missing, drop_unplaced: dropUnplaced })
  }
  const r = tool.result

  return (
    <ToolPage
      id="mend" icon={Dna} tabs={TABS} tab={tab} setTab={setTab}
      title="Mendelian Segregation"
      subtitle="Chi-square test for expected segregation ratios in RIL (1:1) and F2 (1:2:1) populations"
      run={{ label: 'Run Segregation Test', runningLabel: 'Testing…', state: tool.state, status: tool.status, onRun: run }}
    >
      {tab === 'input' && (
        <TabPane id="input" key="input">
          <SectionHeader title="Input Data & Settings">Select your genotype CSV and configure the segregation test.</SectionHeader>
          <Card title="Genotype Data" icon={FileText}>
            <FilePicker value={file} onBrowse={async () => { const p = await tool.browse(); if (p) setFile(p) }}
              placeholder="Select genotype CSV (Marker, Chromosome, then genotype columns)" />
          </Card>
          <div className="row">
            <Card title="Population Type" icon={Sprout} delay={0.05}>
              <RadioChips value={pop} onChange={setPop} options={[{ value: 'RIL', label: 'RIL (1:1)' }, { value: 'F2', label: 'F2 (1:2:1)' }]} />
            </Card>
            <Card title="Test Threshold" icon={Ruler} delay={0.1}>
              <Slider label="p-value threshold" value={pthr} onChange={setPthr} min={0.001} max={0.2} step={0.001} format={(v) => v.toFixed(3)} />
            </Card>
          </div>
          <Card title="Missing Data" icon={Eraser} delay={0.15}>
            <TextField label="Treat these values as missing (comma-separated)" value={missing} onChange={setMissing} />
            <div style={{ marginTop: 10 }}>
              <Toggle label="Keep only markers with a known chromosome" checked={dropUnplaced} onChange={setDropUnplaced} />
            </div>
          </Card>
        </TabPane>
      )}
      {tab === 'results' && (
        <TabPane id="results" key="results">
          <SectionHeader title="Results">Marker-level chi-square test results, distortion plots, and downloads.</SectionHeader>
          {!r ? <EmptyState icon={Dna}>Run the segregation test to see results here.</EmptyState> : (
            <>
              <ImageGallery images={[
                { src: r.manhattan_png, label: '-log10(p) by marker' },
                { src: r.rateplot_png, label: 'Distortion rate by chromosome' },
              ]} />
              <DownloadButtons onSave={tool.save} items={[
                { icon: FileText, label: 'Filtered Markers (CSV)', path: r.results_filtered_csv, name: 'markers_filtered.csv' },
                { icon: FileText, label: 'Full Results (CSV)', path: r.results_full_csv, name: 'markers_results.csv' },
              ]} />
              <pre className="summary-box">{
                `Population type: ${r.pop}\np-value threshold: ${r.pthr}\n` +
                `Markers total: ${r.n_markers_total}  |  Markers tested: ${r.n_markers_kept}\n` +
                `Distorted markers (p < threshold): ${r.n_distorted}\n` +
                `Non-distorted markers (filtered set): ${r.n_filtered_nondistorted}`}</pre>
              <div style={{ height: 20 }} />
              <DataTable rows={r.table} highlight={(c, row) => c === 'p' && row.p != null && row.p < r.pthr} />
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
