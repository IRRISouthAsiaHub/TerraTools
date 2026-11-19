import { useEffect, useRef, useState } from 'react'
import { Dna, SlidersHorizontal, LayoutGrid, FileText, Sprout, Ruler, Eraser, FileSpreadsheet } from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, Toggle, FilePicker, SelectField, TextField, Note,
  EmptyState, ImageGallery, DownloadButtons, DataTable, HtmlBox,
} from '../components/kit.jsx'
import { useBridge, useTool, mockGeneric } from '../lib/bridge.js'
import { usePipe } from '../lib/pipeline.js'

const TABS = [
  { id: 'input', label: 'Input & Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]

const TABLE_COLS = ['Marker', 'Chromosome', 'n', 'c0', 'c1', 'c2', 'het_rate', 'test', 'chisq', 'df', 'p', 'p_adj', 'call', '-log10p']

// Every population has an explicit model: the classes that are tested and their expected ratio (the engine rejects any other type)
const POPS = [
  { value: 'RIL', label: 'RIL (recombinant inbred line)', model: 'Classes 0 and 2 in 1:1 (exact binomial test). Heterozygous calls (1) are not part of the test.' },
  { value: 'DH', label: 'DH (doubled haploid)', model: 'Classes 0 and 2 in 1:1 (exact binomial test). Heterozygous calls (1) are not part of the test.' },
  { value: 'BC1P1', label: 'BC1 to parent 1', model: 'Classes 0 and 1 in 1:1 (exact binomial test). Class 2 cannot occur in this backcross.' },
  { value: 'BC1P2', label: 'BC1 to parent 2', model: 'Classes 1 and 2 in 1:1 (exact binomial test). Class 0 cannot occur in this backcross.' },
  { value: 'F2', label: 'F2', model: 'Classes 0, 1 and 2 in 1:2:1 (exact multinomial test).' },
  { value: 'F3', label: 'F3', model: 'Classes 0, 1 and 2 in 3:2:3 (exact multinomial test).' },
  { value: 'F4', label: 'F4', model: 'Classes 0, 1 and 2 in 7:2:7 (exact multinomial test).' },
]
// p-values keep their precision (0.0000 would hide the strongest signals); booleans read as yes / no
const cell = (c, v) => {
  if (v === null || v === undefined) return ''
  if (c === 'het_rate') return Number(v).toFixed(3)
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
      subtitle="Exact tests of the expected segregation ratio of RIL, DH, backcross (BC1) and F2 to F4 populations"
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
              <SelectField label="Population" value={pop} onChange={setPop} options={POPS.map(({ value, label }) => ({ value, label }))}
                hint={POPS.find((x) => x.value === pop)?.model} />
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
          <Note>Coding assumed: 0 = parent 1 (AA), 1 = heterozygous (Aa), 2 = parent 2 (aa). Every population is tested against its own model, and a type without one (F1, for example) is refused. A marker with too few calls to ever reach α is reported as <b>inconclusive</b>, not as passed.</Note>
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
          <SectionHeader title="Results">Marker-level exact test results, distortion plots, and downloads.</SectionHeader>
          {!r ? <EmptyState icon={Dna}>Run the segregation test to see results here.</EmptyState> : (
            <>
              <ImageGallery images={[
                { src: r.manhattan_png, label: 'Segregation plot: -log10(p) by marker' },
                { src: r.hetplot_png, label: 'Heterozygosity plot: % heterozygous calls by marker' },
                { src: r.rateplot_png, label: 'Distortion summary: % distorted markers by chromosome' },
              ].filter((x) => x.src)} />
              <DownloadButtons onSave={tool.save} items={[
                { icon: FileText, label: 'Filtered Markers (CSV)', path: r.results_filtered_csv, name: `markers_filtered_${r.pop.toLowerCase()}_${r.mtc_method}${r.pthr}.csv` },
                { icon: FileText, label: 'Full Results (CSV)', path: r.results_full_csv, name: `markers_results_${r.pop.toLowerCase()}.csv` },
                { icon: FileSpreadsheet, label: 'LM-ready File (XLSX)', path: r.lm_ready_xlsx, name: `markers_LMready_${r.pop.toLowerCase()}.xlsx` },
              ]} />
              <Note>The LM-ready file holds the markers that pass the quality check (Marker, Chromosome, genotypes): load it directly in Linkage Mapping.</Note>
              <pre className="summary-box">{
                `Population type: ${r.pop_label || r.pop}
α = ${r.pthr}, correction: ${r.mtc_label}
` +
                `Markers in file: ${r.n_markers_total}  |  Markers tested: ${r.n_tested}  |  Inconclusive: ${r.n_inconclusive}
` +
                `Flagged distorted: ${r.n_distorted}
` +
                `Passing the quality check (kept for mapping): ${r.n_pass}` +
                (r.n_inconclusive ? `
Inconclusive (too few calls to reach α, neither distorted nor passed): ${r.n_inconclusive}` : '') +
                (r.n_unexpected_calls ? `
Calls outside the model, not used in the tests: ${r.n_unexpected_calls}` : '')}</pre>
              <div style={{ height: 20 }} />
              {r.table_truncated && <Note>The table lists the {r.table.length.toLocaleString()} most significant of {r.table_total.toLocaleString()} markers. Every marker is in the Full Results (CSV).</Note>}
              <DataTable rows={r.table} columns={TABLE_COLS} renderCell={cell} highlight={(c, row) => row.call === 'distorted' && (c === 'p_adj' || c === 'call')} />
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
