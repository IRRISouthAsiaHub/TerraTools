import { useState } from 'react'
import {
  Waypoints, FolderOpen, ChartColumn, Palette, Scale, Target, Settings2, Database, Download, ChartArea, Map as MapIcon,
  List, Rocket, Search, Table2, FileText, Image,
} from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, SelectField, TextField, ColorField, ChipGroup, FilePicker, Note,
  EmptyState, StatCard, ImageGallery, DataTable, RunBar,
} from '../components/kit.jsx'
import { useBridge, useTool, useGlobals, useSignal, callJson, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'analysis', label: 'Analysis', icon: FolderOpen },
  { id: 'results', label: 'Results', icon: ChartColumn },
  { id: 'plot', label: 'Plot Settings', icon: Palette },
  { id: 'compare', label: 'Compare Models', icon: Scale },
  { id: 'common', label: 'Common QTL', icon: Target },
]
const METHODS = [
  { value: 'SIM', label: 'Simple Interval Mapping (SIM)' }, { value: 'CIM', label: 'Composite Interval Mapping (CIM)' },
  { value: 'MQM', label: 'Multiple QTL Model (MQM)' }, { value: 'RF', label: 'Random Forest QTL' },
  { value: 'SVR', label: 'Support Vector Regression' }, { value: 'XGBOOST', label: 'XGBoost QTL' },
]
const MODEL_DESC = [
  ['SIM', 'Single-marker regression'], ['CIM', 'Composite interval mapping'], ['MQM', 'Multiple QTL model'],
  ['RF', 'Random Forest importance'], ['SVR', 'Kernel-based regression'], ['XGBOOST', 'Gradient boosting'],
]
const QTL_COLS = ['Marker', 'Chromosome', 'Position', 'LOD', 'P_Value', 'PVE_Marker_Percent', 'Additive_Effect', 'Interval_Size_cM', 'Method']
const CMP_COLS = ['Marker', 'Chromosome', 'Position', 'Models_List', 'Mean_LOD', 'Max_LOD', 'Mean_PVE']
const COMMON_COLS = ['Marker', 'Chromosome', 'Position', 'N_Models_Detected', 'Models_List', 'Mean_LOD', 'Max_LOD', 'Mean_PVE']

const PLOT_DEFAULTS = {
  manhattan_height: 800, manhattan_width: 1400, nonsig_point_size: 6, sig_point_size: 10,
  nonsig_color: '#3A7D44', sig_color: '#CC0000', threshold_width: 1.5, threshold_style: 'dash',
  ideogram_height: 900, ideogram_width: 1400, chr_dark: '#2E6B35', chr_light: '#6AAB6E', chr_border: '#1A4A1F',
  qtl_marker_size_min: 0.7, qtl_marker_size_max: 2,
}

const Pic = ({ src, label }) => (src ? <ImageGallery single images={[{ src, label }]} /> : null)

export default function QtlTool() {
  const { bridge } = useBridge('qtlBridge', () => mockGeneric(['errorOccurred']))
  const tool = useTool(bridge, 'run_qtl')
  const [tab, setTab] = useState('analysis')
  const [geno, setGeno] = useState('')
  const [pheno, setPheno] = useState('')
  const [traits, setTraits] = useState(null)
  const [trait, setTrait] = useState('')
  const [method, setMethod] = useState('CIM')
  const [lod, setLod] = useState(2.5)
  const [valFrac, setValFrac] = useState(0.2)
  const [cof, setCof] = useState(5)
  const [win, setWin] = useState(10)
  const [lodDrop, setLodDrop] = useState(1)
  const [plot, setPlot] = useState(PLOT_DEFAULTS)
  const [models, setModels] = useState(['SIM', 'CIM', 'MQM'])
  const [cmp, setCmp] = useState(null)
  const [cmpStatus, setCmpStatus] = useState('Select at least two models to compare.')
  const [cmpBusy, setCmpBusy] = useState(false)
  const [minModels, setMinModels] = useState(2)
  const [common, setCommon] = useState(null)
  const [commonStatus, setCommonStatus] = useState('')
  const setP = (k) => (v) => setPlot((p) => ({ ...p, [k]: v }))

  tool.onResult(() => setTab('results'))
  useSignal(bridge, 'errorOccurred', (m) => tool.fail('Error: ' + m))
  useGlobals({
    showCompareResult: (r) => {
      setCmpBusy(false)
      if (!r) return
      setCmpStatus(r.status === 'success' ? 'Model comparison completed.' : r.message || 'Comparison failed.')
      if (r.status === 'success') setCmp(r)
    },
    showCommonResult: (r) => {
      if (!r) return
      setCommonStatus(r.status === 'error' ? r.message || 'Failed.' : 'Common QTL found.')
      if (r.status !== 'error') setCommon(r)
    },
  })

  const pickPheno = async () => {
    const p = await tool.browse('browse_phenotype')
    if (!p) return
    setPheno(p)
    const cols = await callJson(bridge, 'get_phenotype_columns', p)
    setTraits(cols || [])
    setTrait((cols && cols[0]) || '')
    if (cols && cols.length) tool.setStatus('Phenotype loaded. Select a trait.')
  }

  const run = () => {
    if (!geno || !pheno) return tool.fail('Please select both genotype and phenotype files.')
    if (!trait) return tool.fail('Please select a phenotype/trait.')
    tool.run({
      genotype_file: geno, phenotype_file: pheno, phenotype_column: trait, mapping_method: method,
      lod_threshold: lod, validation_fraction: valFrac, n_cofactors: cof, window_size: win, lod_drop: lodDrop, plot,
    })
  }
  const compare = () => {
    if (models.length < 2) return setCmpStatus('Select at least two models.')
    setCmpStatus('Running model comparison...'); setCmpBusy(true)
    bridge.run_compare(JSON.stringify({
      genotype_file: geno, phenotype_file: pheno, phenotype_column: trait, models, lod_threshold: lod, lod_drop: lodDrop, plot,
    }))
  }
  const findCommon = () => { setCommonStatus('Searching…'); bridge.find_common(JSON.stringify({ min_models: minModels, lod_threshold: lod })) }
  const savePlot = () => { bridge.save_plot_settings(JSON.stringify(plot)); tool.setStatus('Plot settings saved.') }
  const dl = (t) => bridge && bridge.download_result(t)

  const r = tool.result
  const s = (r && r.summary) || {}

  return (
    <ToolPage
      id="qtl" icon={Waypoints} tabs={TABS} tab={tab} setTab={setTab}
      title="QTL Mapping"
      subtitle="Quantitative Trait Locus analysis and multi-model comparison"
      run={tab === 'analysis' ? { label: 'Run QTL Analysis', runningLabel: 'Mapping QTL…', state: tool.state, status: tool.status, onRun: run } : undefined}
    >
      {tab === 'analysis' && (
        <TabPane id="analysis" key="analysis">
          <div className="grid-2" style={{ gridTemplateColumns: 'repeat(auto-fit, minmax(330px, 1fr))', alignItems: 'start' }}>
            <Card title="Input Files" icon={FolderOpen}>
              <div className="field">
                <label className="field-title">Genotype CSV</label>
                <FilePicker value={geno} placeholder="Select genotype CSV…" onBrowse={async () => { const p = await tool.browse('browse_genotype'); if (p) setGeno(p) }} />
                <div className="field-hint">Expected format: row 1 = marker names, row 2 = chromosome, row 3 = position, followed by sample genotype rows.</div>
              </div>
              <div className="field">
                <label className="field-title">Phenotype CSV</label>
                <FilePicker value={pheno} placeholder="Select phenotype CSV…" onBrowse={pickPheno} />
              </div>
              <SelectField label="Trait / Phenotype" value={trait} onChange={setTrait}
                options={traits && traits.length ? traits : [{ value: '', label: traits ? 'No traits found' : 'Load phenotype file first' }]} />
            </Card>
            <Card title="Analysis Settings" icon={Settings2} delay={0.06}>
              <SelectField label="Mapping Method" value={method} onChange={setMethod} options={METHODS} />
              <Slider label="LOD Threshold" value={lod} onChange={setLod} min={0.5} max={10} step={0.1} format={(v) => v.toFixed(1)} />
              <Slider label="Validation Fraction" value={valFrac} onChange={setValFrac} min={0.1} max={0.4} step={0.05} format={(v) => v.toFixed(2)} />
              {method === 'CIM' && (
                <div className="grid-2" style={{ marginTop: 6 }}>
                  <TextField label="CIM cofactors (K)" type="number" min={1} max={20} value={cof} onChange={setCof} />
                  <TextField label="Window (± markers)" type="number" min={3} max={30} value={win} onChange={setWin} />
                </div>
              )}
              <TextField label="LOD Drop — Support Interval" type="number" min={0.3} max={2} step={0.1} value={lodDrop} onChange={setLodDrop} />
            </Card>
            <Card title="Dataset Information" icon={Database} delay={0.12}>
              <div className="summary-grid" style={{ gridTemplateColumns: 'repeat(3, 1fr)' }}>
                {[['Samples', s.samples], ['Markers', s.markers], ['Chromosomes', s.chromosomes]].map(([l, v]) => (
                  <div className="summary-card" key={l}><div className="val">{v ?? '—'}</div><div className="lbl">{l}</div></div>
                ))}
              </div>
              <div className="plot-title">Available Models</div>
              <div className="model-grid">
                {MODEL_DESC.map(([k, d]) => <div className="model-card" key={k}><strong>{k}</strong><p>{d}</p></div>)}
              </div>
              <Note>QTL Mapper supports six mapping approaches and chromosome-wise visualization of detected QTL.</Note>
            </Card>
          </div>
        </TabPane>
      )}
      {tab === 'results' && (
        <TabPane id="results" key="results">
          <SectionHeader title="Results & Downloads" />
          {!r ? <EmptyState icon={ChartArea}>Run a QTL analysis to display the results.</EmptyState> : (
            <>
              <div className="downloads">
                <button className="btn btn-secondary" onClick={() => dl('manhattan')}><Image size={15} /> Manhattan PNG</button>
                <button className="btn btn-secondary" onClick={() => dl('ideogram')}><Image size={15} /> Ideogram PNG</button>
                <button className="btn btn-secondary" onClick={() => dl('qtl_csv')}><FileText size={15} /> QTL Table CSV</button>
                <button className="btn btn-secondary" onClick={() => dl('full_csv')}><FileText size={15} /> Full Scan CSV</button>
              </div>
              <div className="plot-title">Chromosome-wise Manhattan Plot</div>
              <Pic src={r.manhattan} label="QTL Manhattan Plot" />
              <div className="plot-title">Chromosome Ideogram</div>
              <Pic src={r.ideogram} label="QTL Chromosome Ideogram" />
              <div className="plot-title">Significant QTL</div>
              <DataTable rows={r.qtl_table} columns={QTL_COLS} maxHeight={420} />
            </>
          )}
        </TabPane>
      )}
      {tab === 'plot' && (
        <TabPane id="plot" key="plot">
          <SectionHeader title="Plot Settings">Applied the next time you run an analysis or comparison.</SectionHeader>
          <div className="row">
            <Card title="Manhattan Plot" icon={ChartArea}>
              <div className="plot-title">Plot Dimensions</div>
              <div className="grid-2">
                <TextField label="Height (px)" type="number" min={400} max={1200} step={50} value={plot.manhattan_height} onChange={setP('manhattan_height')} />
                <TextField label="Width (px)" type="number" min={800} max={2000} step={100} value={plot.manhattan_width} onChange={setP('manhattan_width')} />
              </div>
              <div className="plot-title" style={{ marginTop: 16 }}>Point Styling</div>
              <div className="grid-2">
                <TextField label="Non-significant size" type="number" min={2} max={12} value={plot.nonsig_point_size} onChange={setP('nonsig_point_size')} />
                <TextField label="Significant size" type="number" min={4} max={20} value={plot.sig_point_size} onChange={setP('sig_point_size')} />
                <ColorField label="Non-significant color" value={plot.nonsig_color} onChange={setP('nonsig_color')} />
                <ColorField label="Significant color" value={plot.sig_color} onChange={setP('sig_color')} />
              </div>
              <div className="plot-title" style={{ marginTop: 16 }}>Threshold</div>
              <div className="grid-2">
                <TextField label="Line width" type="number" min={0.5} max={4} step={0.5} value={plot.threshold_width} onChange={setP('threshold_width')} />
                <SelectField label="Line style" value={plot.threshold_style} onChange={setP('threshold_style')}
                  options={[{ value: 'solid', label: 'Solid' }, { value: 'dash', label: 'Dashed' }, { value: 'dot', label: 'Dotted' }, { value: 'dashdot', label: 'Dash-Dot' }]} />
              </div>
            </Card>
            <Card title="Ideogram" icon={MapIcon} delay={0.06}>
              <div className="plot-title">Plot Dimensions</div>
              <div className="grid-2">
                <TextField label="Height (px)" type="number" min={400} max={1200} step={50} value={plot.ideogram_height} onChange={setP('ideogram_height')} />
                <TextField label="Width (px)" type="number" min={800} max={2000} step={100} value={plot.ideogram_width} onChange={setP('ideogram_width')} />
              </div>
              <div className="plot-title" style={{ marginTop: 16 }}>Chromosome Colors</div>
              <div className="grid-2">
                <ColorField label="Dark band" value={plot.chr_dark} onChange={setP('chr_dark')} />
                <ColorField label="Light band" value={plot.chr_light} onChange={setP('chr_light')} />
                <ColorField label="Border" value={plot.chr_border} onChange={setP('chr_border')} />
              </div>
              <div className="plot-title" style={{ marginTop: 16 }}>QTL Marker Sizing</div>
              <div className="grid-2">
                <TextField label="Minimum size" type="number" min={0.3} max={1.5} step={0.1} value={plot.qtl_marker_size_min} onChange={setP('qtl_marker_size_min')} />
                <TextField label="Maximum size" type="number" min={1} max={3} step={0.1} value={plot.qtl_marker_size_max} onChange={setP('qtl_marker_size_max')} />
              </div>
              <button className="btn btn-secondary" style={{ marginTop: 18 }} onClick={savePlot}>Save Plot Settings</button>
            </Card>
          </div>
        </TabPane>
      )}
      {tab === 'compare' && (
        <TabPane id="compare" key="compare">
          <Card title="Compare QTL Mapping Models" icon={Scale}>
            <div className="plot-title">Select Models</div>
            <ChipGroup options={METHODS.map((m) => m.value)} value={models} onChange={setModels} />
            <button className="btn btn-secondary" style={{ marginTop: 16 }} onClick={compare} disabled={cmpBusy}><Rocket size={16} /> {cmpBusy ? 'Running…' : 'Run Comparison'}</button>
            <div className="info-line">{cmpStatus}</div>
          </Card>
          <div className="plot-grid-2">
            <Card title="Multi-Model Ideogram" icon={MapIcon}>{cmp && cmp.ideogram ? <Pic src={cmp.ideogram} label="Comparison Ideogram" /> : <div className="empty-state" style={{ padding: 20 }}>Comparison ideogram will appear here.</div>}</Card>
            <Card title="QTL Counts per Model" icon={ChartColumn}>{cmp && cmp.bar ? <Pic src={cmp.bar} label="QTL Comparison" /> : <div className="empty-state" style={{ padding: 20 }}>Comparison chart will appear here.</div>}</Card>
          </div>
          <Card title="Comparison Table" icon={Table2}><DataTable rows={cmp && cmp.table} columns={CMP_COLS} maxHeight={380} /></Card>
        </TabPane>
      )}
      {tab === 'common' && (
        <TabPane id="common" key="common">
          <Card title="Common QTL Across Models" icon={Target}>
            <div style={{ maxWidth: 450 }}>
              <Slider label="Minimum Models Detecting QTL" value={minModels} onChange={setMinModels} min={2} max={6} />
            </div>
            <button className="btn btn-secondary" style={{ marginTop: 16 }} onClick={findCommon}><Search size={16} /> Find Common QTL</button>
            {commonStatus && <div className="info-line">{commonStatus}</div>}
          </Card>
          <div className="plot-grid-2">
            <Card title="Common QTL Ideogram" icon={MapIcon}>{common && common.ideogram ? <Pic src={common.ideogram} label="Common QTL Ideogram" /> : <div className="empty-state" style={{ padding: 20 }}>Common-QTL ideogram will appear here.</div>}</Card>
            <Card title="Common QTL Counts" icon={ChartColumn}>{common && common.bar ? <Pic src={common.bar} label="Common QTL Counts" /> : <div className="empty-state" style={{ padding: 20 }}>Common-QTL counts will appear here.</div>}</Card>
          </div>
          <Card title="Common QTL Table" icon={Table2}><DataTable rows={common && common.table} columns={COMMON_COLS} maxHeight={380} /></Card>
        </TabPane>
      )}
    </ToolPage>
  )
}
