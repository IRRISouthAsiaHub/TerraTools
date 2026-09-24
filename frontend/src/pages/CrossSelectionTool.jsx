import { useState } from 'react'
import { GitBranch, FolderOpen, ChartColumn, Table2, Trophy, Dna, Settings2, Database, Search, Filter, Download, FileText, Star } from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, FilePicker, ColorField, SelectField, TextField, Note,
  EmptyState, StatCard, ImageGallery, DataTable,
} from '../components/kit.jsx'
import { useBridge, useTool, callJson, useSignal, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'configure', label: 'Configure', icon: FolderOpen },
  { id: 'dashboard', label: 'Dashboard', icon: ChartColumn },
  { id: 'allcrosses', label: 'All Crosses', icon: Table2 },
  { id: 'topcrosses', label: 'Top Crosses', icon: Trophy },
  { id: 'parents', label: 'Parent Summary', icon: Dna },
]
const REL = ['All', 'Unrelated', 'Distantly_Related', 'Moderately_Related', 'Closely_Related']
const PERF = ['All', 'Low', 'Medium', 'High', 'Very_High']

const badge = (c, v) => {
  if (c === 'Relatedness_Category' || c === 'Relatedness') return v ? <span className={`badge rel-${v}`}>{String(v).replace(/_/g, ' ')}</span> : '—'
  if (v === null || v === undefined) return '—'
  return undefined
}

export default function CrossSelectionTool() {
  const { bridge } = useBridge('crossBridge', () => mockGeneric(['errorOccurred']))
  const tool = useTool(bridge, 'run_cross_selection')
  const [tab, setTab] = useState('configure')
  const [criterion, setCriterion] = useState('')
  const [kmat, setKmat] = useState('')
  const [validation, setValidation] = useState(null)
  const [threshold, setThreshold] = useState(0.15)
  const [maxPer, setMaxPer] = useState(5)
  const [povColor, setPovColor] = useState('#3A7D44')
  const [kinColor, setKinColor] = useState('#C8760A')
  const [dataset, setDataset] = useState({})
  const [rows, setRows] = useState([])
  const [low, setLow] = useState(null)
  const [med, setMed] = useState(null)
  const [flt, setFlt] = useState({ relatedness: 'All', performance: 'All', kinship_max: 1 })

  useSignal(bridge, 'errorOccurred', (m) => tool.fail('Error: ' + m))

  const band = async (kinship_min, kinship_max) => {
    const r = await callJson(bridge, 'filter_crosses', JSON.stringify({ relatedness: 'All', performance: 'All', kinship_min, kinship_max }))
    return r && r.status === 'success' ? r.rows.slice(0, 50) : []
  }
  tool.onResult(async (r) => {
    const s = r.summary || {}
    setDataset({ rows: s.n_criterion_rows, genos: s.n_kmat_genotypes, overlap: s.n_overlap })
    setRows(r.top_crosses || [])
    setTab('dashboard')
    setLow(await band(0, 0.125))
    setMed(await band(0.125, 0.25))
  })

  const preview = async () => {
    if (!criterion && !kmat) return tool.setStatus('Select at least one file to preview.')
    const r = await callJson(bridge, 'preview_data', criterion, kmat)
    setValidation(r)
    if (r && r.status !== 'error') setDataset({ rows: r.criterion && r.criterion.rows, genos: r.kmat && r.kmat.rows, overlap: r.overlap })
  }
  const run = () => {
    if (!criterion || !kmat) return tool.fail('Please select both a criterion file and a kinship matrix.')
    tool.run({ criterion_file: criterion, kmat_file: kmat, threshold, max_per_cluster: maxPer, plot: { pov_color: povColor, kinship_color: kinColor } })
  }
  const applyFilter = async () => {
    if (!tool.result) return tool.setStatus('Run an analysis first.')
    const r = await callJson(bridge, 'filter_crosses', JSON.stringify(flt))
    if (!r) return
    if (r.status === 'error') return tool.setStatus(r.message)
    setRows(r.rows)
  }

  const r = tool.result
  const s = (r && r.summary) || {}
  const pl = (r && r.plots) || {}
  const dl = (type) => bridge && bridge.download_result(type)
  const fmt = (v) => v ?? 0

  const topTable = (data, title, icon, exportKey) => (
    <Card title={title} icon={icon}>
      {data ? <DataTable rows={data} columns={['Parent1', 'Parent2', 'Predicted_Offspring_Value', 'Kinship_Coefficient', 'Relatedness_Category']} renderCell={badge} maxHeight={420} />
        : <div className="empty-state" style={{ padding: 24 }}>No results yet.</div>}
      <div style={{ marginTop: 12 }}><button className="btn btn-secondary" onClick={() => bridge.export_top_crosses(exportKey)}><FileText size={16} /> Export (CSV)</button></div>
    </Card>
  )

  return (
    <ToolPage
      id="cross" icon={GitBranch} tabs={TABS} tab={tab} setTab={setTab}
      title="Cross Selection"
      subtitle="Half-diallel cross optimization with kinship-aware parent selection"
      run={tab === 'configure' ? { label: 'Run Cross Selection', runningLabel: 'Optimizing crosses…', state: tool.state, status: tool.status, onRun: run } : undefined}
    >
      {tab === 'configure' && (
        <TabPane id="configure" key="configure">
          <div className="grid-2" style={{ gridTemplateColumns: 'repeat(auto-fit, minmax(320px, 1fr))', alignItems: 'start' }}>
            <Card title="Input Files" icon={FolderOpen}>
              <div className="field">
                <label className="field-title">Criterion CSV</label>
                <FilePicker value={criterion} placeholder="Select criterion CSV…" onBrowse={async () => { const p = await tool.browse('browse_criterion'); if (p) setCriterion(p) }} />
                <div className="field-hint">Expected columns: Genotype, Criterion (the selection index or breeding value used to rank parents).</div>
              </div>
              <div className="field">
                <label className="field-title">Kinship Matrix CSV</label>
                <FilePicker value={kmat} placeholder="Select kinship matrix CSV…" onBrowse={async () => { const p = await tool.browse('browse_kmat'); if (p) setKmat(p) }} />
                <div className="field-hint">Square matrix, genotype names in the first column and as the header row.</div>
              </div>
              <button className="btn btn-secondary" style={{ width: '100%', marginTop: 16 }} onClick={preview}><Search size={16} /> Preview & Validate</button>
              {validation && (validation.status === 'error'
                ? <div className="info-line warn">{validation.message}</div>
                : <>
                    {validation.criterion && <div className="info-line">Criterion loaded: {validation.criterion.rows} rows, {validation.criterion.unique_genotypes} unique genotypes</div>}
                    {validation.kmat && <div className="info-line">Kinship matrix loaded: {validation.kmat.rows} × {validation.kmat.cols}</div>}
                    {validation.overlap !== undefined && <div className={`info-line ${validation.overlap === 0 ? 'warn' : ''}`}>Overlap: {validation.overlap} genotypes match between files</div>}
                  </>)}
            </Card>
            <Card title="Thinning & Cross Settings" icon={Settings2} delay={0.06}>
              <Slider label="Kinship Threshold" value={threshold} onChange={setThreshold} min={0.05} max={0.5} step={0.01} format={(v) => v.toFixed(2)} />
              <Slider label="Max per Cluster" value={maxPer} onChange={setMaxPer} min={1} max={20} />
              <div className="grid-2" style={{ marginTop: 6 }}>
                <ColorField label="Offspring Value Color" value={povColor} onChange={setPovColor} />
                <ColorField label="Kinship Color" value={kinColor} onChange={setKinColor} />
              </div>
            </Card>
            <Card title="Dataset Information" icon={Database} delay={0.12}>
              <div className="summary-grid" style={{ gridTemplateColumns: 'repeat(3, 1fr)' }}>
                {[['Criterion Rows', dataset.rows], ['Kinship Genotypes', dataset.genos], ['Overlap', dataset.overlap]].map(([l, v]) => (
                  <div className="summary-card" key={l}><div className="val">{v ?? '—'}</div><div className="lbl">{l}</div></div>
                ))}
              </div>
              <Note>
                Genotypes are first thinned so no cluster of closely related lines contributes more than the chosen maximum, then every remaining
                pair is crossed (half-diallel). Predicted Offspring Value is the parents' mid-point criterion value; Kinship Coefficient comes from
                the matrix you supplied.
              </Note>
            </Card>
          </div>
        </TabPane>
      )}
      {tab === 'dashboard' && (
        <TabPane id="dashboard" key="dashboard">
          <SectionHeader title="Analysis Summary" />
          {!r ? <EmptyState icon={ChartColumn}>Run an analysis to see this dashboard.</EmptyState> : (
            <>
              <div className="summary-grid" style={{ gridTemplateColumns: 'repeat(auto-fit, minmax(130px, 1fr))' }}>
                <StatCard className="summary-card" value={fmt(s.n_crosses)} label="Total Crosses" />
                <StatCard className="summary-card" value={fmt(s.n_parents)} label="Parents" delay={0.05} />
                <StatCard className="summary-card" value={Number(fmt(s.mean_pov))} label="Mean POV" delay={0.1} />
                <StatCard className="summary-card" value={Number(fmt(s.max_pov))} label="Max POV" delay={0.15} />
                <StatCard className="summary-card" value={fmt(s.n_low_kinship)} label="K < 0.125" delay={0.2} />
                <StatCard className="summary-card" value={Number(fmt(s.mean_kinship))} label="Mean Kinship" delay={0.25} />
              </div>
              <div className="downloads">
                {[['hist_pov', 'POV Histogram'], ['hist_kinship', 'Kinship Histogram'], ['bar_relatedness', 'Relatedness Bar'], ['scatter_vk', 'Value/Kinship Scatter']].map(([k, l]) => (
                  <button key={k} className="btn btn-secondary" onClick={() => dl(k)}><Download size={15} /> {l}</button>
                ))}
                <button className="btn btn-secondary" onClick={() => dl('crosses_csv')}><FileText size={15} /> Full Crosses CSV</button>
              </div>
              <ImageGallery images={[
                pl.hist_pov && { src: pl.hist_pov, label: 'Predicted Offspring Value Distribution' },
                pl.hist_kinship && { src: pl.hist_kinship, label: 'Kinship Coefficient Distribution' },
                pl.bar_relatedness && { src: pl.bar_relatedness, label: 'Crosses by Relatedness Category' },
                pl.scatter_vk && { src: pl.scatter_vk, label: 'Value vs Kinship Scatter' },
              ].filter(Boolean)} />
            </>
          )}
        </TabPane>
      )}
      {tab === 'allcrosses' && (
        <TabPane id="allcrosses" key="allcrosses">
          <SectionHeader title="Half-Diallel Cross Table" />
          {!r ? <EmptyState icon={Table2}>No crosses yet. Run an analysis first.</EmptyState> : (
            <>
              <div className="grid-4" style={{ marginBottom: 14 }}>
                <SelectField label="Filter Relatedness" value={flt.relatedness} onChange={(v) => setFlt({ ...flt, relatedness: v })} options={REL} />
                <SelectField label="Filter Performance" value={flt.performance} onChange={(v) => setFlt({ ...flt, performance: v })} options={PERF} />
                <TextField label="Max Kinship" type="number" min={0} max={1} step={0.01} value={flt.kinship_max} onChange={(v) => setFlt({ ...flt, kinship_max: v })} />
                <button className="btn btn-secondary" style={{ height: 41 }} onClick={applyFilter}><Filter size={15} /> Apply Filter</button>
              </div>
              <div style={{ marginBottom: 12 }}>
                <button className="btn btn-secondary" onClick={() => bridge.export_filtered_crosses(JSON.stringify({ ...flt, kinship_min: 0 }))}><FileText size={15} /> Export Filtered Table (CSV)</button>
              </div>
              <DataTable rows={rows} renderCell={badge} maxHeight={520}
                columns={['Parent1', 'Parent2', 'P1_Criterion_Value', 'P2_Criterion_Value', 'Predicted_Offspring_Value', 'Kinship_Coefficient', 'Relatedness_Category', 'Performance_Category']} />
            </>
          )}
        </TabPane>
      )}
      {tab === 'topcrosses' && (
        <TabPane id="topcrosses" key="topcrosses">
          {!r ? <EmptyState icon={Trophy}>No results yet.</EmptyState> : (
            <div className="grid-2" style={{ gridTemplateColumns: 'repeat(auto-fit, minmax(420px, 1fr))', alignItems: 'start' }}>
              {topTable(low, 'Top 50 — Low Kinship (K < 0.125)', Star, 'low')}
              {topTable(med, 'Top 50 — Medium Kinship (0.125 ≤ K < 0.25)', Star, 'medium')}
            </div>
          )}
        </TabPane>
      )}
      {tab === 'parents' && (
        <TabPane id="parents" key="parents">
          {!r ? <EmptyState icon={Dna}>No results yet.</EmptyState> : (
            <div className="grid-2" style={{ gridTemplateColumns: 'repeat(auto-fit, minmax(380px, 1fr))', alignItems: 'start' }}>
              <Card title="Parent Usage & Criterion Values" icon={Dna}>
                <DataTable rows={r.parent_table} columns={['Genotype', 'Times_Used', 'Criterion_Value']} maxHeight={480} />
                <div style={{ marginTop: 12 }}><button className="btn btn-secondary" onClick={() => dl('parent_csv')}><FileText size={15} /> Download Parent Table</button></div>
              </Card>
              <Card title="Parent Criterion Values (Ranked)" icon={ChartColumn}>
                {pl.parent_bar ? <ImageGallery single images={[{ src: pl.parent_bar, label: 'Parent criterion values' }]} /> : <div className="empty-state">No chart.</div>}
              </Card>
            </div>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
