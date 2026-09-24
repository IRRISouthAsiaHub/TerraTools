import { useState } from 'react'
import { Dna, SlidersHorizontal, LayoutGrid, FolderOpen, Info, PaintBucket, Grid3x3, ChartColumn, Table2, Download } from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, SelectField, ChipGroup, FilePicker, Note, SubTabs,
  EmptyState, StatCard, ImageGallery, DataTable,
} from '../components/kit.jsx'
import { useBridge, useTool, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'input', label: 'Input & Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]
const SUB = [
  { id: 'heat', label: 'Heatmaps', icon: Grid3x3 },
  { id: 'hist', label: 'Histogram', icon: ChartColumn },
  { id: 'table', label: 'Data Table', icon: Table2 },
  { id: 'dl', label: 'Downloads', icon: Download },
]
const MATRICES = { Additive: 'additive', Dominance: 'dominance' }

const NotComputed = ({ label }) => (
  <div className="result-card" style={{ cursor: 'default' }}>
    <div className="empty-state" style={{ padding: '40px 10px' }}><Info size={30} /><div>Not computed.</div></div>
    <div className="label">{label}</div>
  </div>
)

function Pair({ items }) {
  const have = items.filter((i) => i.src)
  return (
    <div className="plot-grid-2">
      {items.map((i) => (i.src
        ? <ImageGallery key={i.label} single images={[i]} />
        : <NotComputed key={i.label} label={i.label} />))}
      {have.length === 0 && null}
    </div>
  )
}

export default function AdditiveDominanceTool() {
  const { bridge } = useBridge('additiveDominanceBridge', () => mockGeneric())
  const tool = useTool(bridge, 'run_matrices')
  const [tab, setTab] = useState('input')
  const [sub, setSub] = useState('heat')
  const [file, setFile] = useState('')
  const [impute, setImpute] = useState('mean')
  const [types, setTypes] = useState(['Additive', 'Dominance'])
  const [which, setWhich] = useState('A')
  tool.onResult((r) => {
    setTab('results'); setSub('heat')
    setWhich(r.table_additive && r.table_additive.length ? 'A' : 'D')
  })

  const run = () => {
    if (!file) return tool.fail('Please select a SNP file first.')
    if (!types.length) return tool.fail('Select at least one matrix type.')
    tool.run({ snp_file: file, impute_method: impute, matrix_types: types.map((t) => MATRICES[t]) })
  }
  const r = tool.result
  const s = (r && r.summary) || {}
  const msg = r
    ? `Done. Built ${(s.matrix_types || []).join(' + ')} matrix for ${s.n_samples} samples × ${s.n_snps} markers.` +
      (s.mono_removed ? ` (${s.mono_removed} monomorphic marker(s) removed)` : '')
    : 'Run the analysis to see results here.'
  const rows = r && (which === 'A' ? r.table_additive : r.table_dominance)

  return (
    <ToolPage
      id="ad" icon={Dna} tabs={TABS} tab={tab} setTab={setTab}
      title="Genomic Relationship Matrix Builder"
      subtitle="Additive (rrBLUP A.mat) and dominance (sommer D.mat) genomic relationship matrices from SNP dosage data"
      run={{ label: 'Compute Matrices', runningLabel: 'Computing…', state: tool.state, status: tool.status, onRun: run }}
    >
      {tab === 'input' && (
        <TabPane id="input" key="input">
          <SectionHeader title="SNP Data" />
          <Card title="SNP File" icon={FolderOpen}>
            <FilePicker value={file} placeholder="Select a CSV or TSV file…"
              onBrowse={async () => { const p = await tool.browse('browse_snp_file'); if (p) setFile(p) }} />
            <Note>
              Rows = samples, columns = SNP markers • First column = sample IDs<br />
              Marker encoding: 0 (AA) / 1 (Aa) / 2 (aa) • Missing values are imputed automatically
            </Note>
          </Card>
          <SectionHeader title="Parameters" />
          <div className="row">
            <Card title="Imputation" icon={PaintBucket} delay={0.05}>
              <SelectField label="Missing value method" value={impute} onChange={setImpute}
                options={[{ value: 'mean', label: 'Mean imputation' }, { value: 'median', label: 'Median imputation' }, { value: 'zero', label: 'Zero fill' }]} />
            </Card>
            <Card title="Matrix Type" icon={Grid3x3} delay={0.1}>
              <ChipGroup options={Object.keys(MATRICES)} value={types} onChange={setTypes} />
            </Card>
          </div>
        </TabPane>
      )}
      {tab === 'results' && (
        <TabPane id="results" key="results">
          <SectionHeader title="Results">{msg}</SectionHeader>
          {!r ? <EmptyState icon={Dna}>Run the analysis to see results here.</EmptyState> : (
            <>
              <div className="summary-grid">
                <StatCard className="summary-card" value={s.n_samples ?? 0} label="Samples" />
                <StatCard className="summary-card" value={s.n_snps ?? 0} label="SNP Markers" delay={0.06} />
                <StatCard className="summary-card" value={Number(s.pct_missing ?? 0)} label="Missing (%)" delay={0.12} />
                <StatCard className="summary-card" value={s.maf_pass ?? 0} label="MAF > 0.05" delay={0.18} />
              </div>
              <SubTabs tabs={SUB} value={sub} onChange={setSub} id="adsub" />
              {sub === 'heat' && <Pair items={[{ src: r.heatmap_additive, label: 'Additive Matrix' }, { src: r.heatmap_dominance, label: 'Dominance Matrix' }]} />}
              {sub === 'hist' && <Pair items={[{ src: r.hist_additive, label: 'Additive — off-diagonal distribution' }, { src: r.hist_dominance, label: 'Dominance — off-diagonal distribution' }]} />}
              {sub === 'table' && (
                <>
                  <div style={{ maxWidth: 240 }}>
                    <SelectField label="View matrix" value={which} onChange={setWhich} options={[{ value: 'A', label: 'Additive' }, { value: 'D', label: 'Dominance' }]} />
                  </div>
                  <div style={{ height: 14 }} />
                  {rows && rows.length ? <DataTable rows={rows} maxHeight={480} /> : <div className="empty-state">Not computed.</div>}
                </>
              )}
              {sub === 'dl' && (
                <div className="downloads">
                  {r.additive_csv && <button className="btn btn-secondary" onClick={() => bridge.download_result('additive')}><Download size={16} /> Additive Matrix (CSV)</button>}
                  {r.dominance_csv && <button className="btn btn-secondary" onClick={() => bridge.download_result('dominance')}><Download size={16} /> Dominance Matrix (CSV)</button>}
                  {!r.additive_csv && !r.dominance_csv && <EmptyState icon={Download}>Nothing to download.</EmptyState>}
                </div>
              )}
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
