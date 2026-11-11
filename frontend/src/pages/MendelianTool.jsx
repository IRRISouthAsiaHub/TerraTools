import { useEffect, useRef, useState } from 'react'
import { Dna, SlidersHorizontal, LayoutGrid, FileText, Sprout, Ruler, Eraser, FileSpreadsheet } from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, Toggle, FilePicker, RadioChips, SelectField, TextField, Note,
  EmptyState, ImageGallery, DownloadButtons, DataTable, HtmlBox,
} from '../components/kit.jsx'
import { useBridge, useTool, mockGeneric } from '../lib/bridge.js'
import { usePipe } from '../lib/pipeline.js'

const TABS = [
  { id: 'input', label: 'Input & Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]

const TABLE_COLS = ['Marker', 'Chromosome', 'n', 'c0', 'c1', 'c2', 'het_rate', 'chisq', 'df', 'p', 'p_adj', 'distorted', 'low_power', '-log10p']
// p-values keep their precision (0.0000 would hide the strongest signals); booleans read as yes / no
const cell = (c, v) => {
  if (v === null || v === undefined) return ''
  if (c === 'p' || c === 'p_adj') return v < 0.0001 ? Number(v).toExponential(2) : Number(v).toFixed(4)
  if (typeof v === 'boolean') return v ? 'yes' : 'no'
  return undefined
}

export default function MendelianTool({ navigate }) {
  const { bridge } = useBridge('mendelianBridge', () => mockGeneric())
  const tool = useTool(bridge, 'run_mendelian')
  const [tab, setTab] = useState('input')
  const [file, setFile] = useState('')
  const [pop, setPop] = useState('RIL')
  const [pthr, setPthr] = useState(0.05)
  const [mtc, setMtc] = useState('fdr')
  const [missing, setMissing] = useState('NA,.,-9,-1')
  const [dropUnplaced, setDropUnplaced] = useState(true)
  tool.onResult(() => setTab('results'))

  // "Open in Mendelian Segregation" from the Genofile Converter
  const pipe = usePipe()
  const lastStamp = useRef(0)
  useEffect(() => {
    const g = pipe.genofile
    if (g && g.target === 'mendelian' && g.stamp !== lastStamp.current) { lastStamp.current = g.stamp; setFile(g.mendelian) }
  }, [pipe.genofile])
  const run = () => {
    if (!file) return tool.fail('Please select an input file.')
    tool.run({ input_file: file, pop, pthr, mtc_method: mtc, missing_codes: missing, drop_unplaced: dropUnplaced })
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
            <Note>
              Genotypes from a genotyping company (HapMap, VCF, A/G calls)? <a href="#/genofile" onClick={(e) => { e.preventDefault(); navigate && navigate('genofile') }}>Convert them with the Genofile Converter</a> first:
              it writes this CSV.
            </Note>
          </Card>
          <div className="row">
            <Card title="Population Type" icon={Sprout} delay={0.05}>
              <RadioChips value={pop} onChange={setPop} options={[{ value: 'RIL', label: 'RIL (1:1)' }, { value: 'F2', label: 'F2 (1:2:1)' }]} />
            </Card>
            <Card title="Significance" icon={Ruler} delay={0.1}>
              <Slider label="Significance threshold (α)" value={pthr} onChange={setPthr} min={0.001} max={0.2} step={0.001} format={(v) => v.toFixed(3)} />
              <SelectField label="Multiple-testing correction" value={mtc} onChange={setMtc} options={[
                { value: 'fdr', label: 'FDR (Benjamini-Hochberg)' },
                { value: 'bonferroni', label: 'Bonferroni' },
                { value: 'none', label: 'None (raw p, not recommended)' },
              ]} hint="Hundreds of markers tested at a raw α of 0.05 flag many as distorted by chance alone. The correction is applied on top of α." />
            </Card>
          </div>
          <Note>Coding assumed: 0 = AA, 1 = Aa, 2 = aa. In a RIL population the heterozygotes (1) are ignored in the test.</Note>
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
              ].filter((x) => x.src)} />
              <DownloadButtons onSave={tool.save} items={[
                { icon: FileText, label: 'Filtered Markers (CSV)', path: r.results_filtered_csv, name: `markers_filtered_${r.pop.toLowerCase()}_${r.mtc_method}${r.pthr}.csv` },
                { icon: FileText, label: 'Full Results (CSV)', path: r.results_full_csv, name: `markers_results_${r.pop.toLowerCase()}.csv` },
                { icon: FileSpreadsheet, label: 'LM-ready File (XLSX)', path: r.lm_ready_xlsx, name: `markers_LMready_${r.pop.toLowerCase()}.xlsx` },
              ]} />
              <Note>The LM-ready file holds the markers that pass the quality check (Marker, Chromosome, genotypes): load it directly in Linkage Mapping.</Note>
              <pre className="summary-box">{
                `Population type: ${r.pop}
α = ${r.pthr}, correction: ${r.mtc_label}
` +
                `Markers in file: ${r.n_markers_total}  |  Markers tested: ${r.n_tested}
` +
                `Flagged distorted: ${r.n_distorted}
` +
                `Passing the quality check (kept for mapping): ${r.n_pass}
` +
                `Expected cell count < 5 (read their p-values with caution): ${r.n_low_power}`}</pre>
              <div style={{ height: 20 }} />
              <DataTable rows={r.table} columns={TABLE_COLS} renderCell={cell} highlight={(c, row) => row.distorted === true && (c === 'p_adj' || c === 'distorted')} />
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
