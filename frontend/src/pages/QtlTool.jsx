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
const QTL_COLS = ['QTL_ID', 'Chromosome', 'Left_Marker', 'Peak_Marker', 'Right_Marker', 'Left_Position', 'Peak_Position', 'Right_Position',
  'Interval_Size_cM', 'CI_Start_cM', 'CI_End_cM', 'LOD', 'P_Value', 'Genomewide_P', 'PVE_Marker_Percent', 'PVE_CIM_Percent',
  'Additive_Effect', 'SE_Additive', 'Dominance_Effect', 'Dominance_Class', 'Method']
const CMP_COLS = ['Region_ID', 'Chromosome', 'Left_Marker', 'Peak_Marker', 'Right_Marker', 'Left_Position', 'Right_Position', 'Region_Size_cM',
  'Left_Mb', 'Peak_Mb', 'Right_Mb', 'Order_Check', 'N_Models_Detected', 'Models_List', 'Peak_Range_cM', 'Mean_LOD', 'Max_LOD', 'Mean_PVE']
const COMMON_COLS = CMP_COLS

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
  const [valFrac, setValFrac] = useState(0)
  const [cof, setCof] = useState(5)
  const [win, setWin] = useState(10)
  const [lodDrop, setLodDrop] = useState(1.5)
  const [peakDrop, setPeakDrop] = useState(2)
  const [pThr, setPThr] = useState(0.05)
  const [thrMode, setThrMode] = useState('lod')
  const [nPerm, setNPerm] = useState(1000)
  const [permAlpha, setPermAlpha] = useState(0.05)
  const [plot, setPlot] = useState(PLOT_DEFAULTS)
  const [models, setModels] = useState(['SIM', 'CIM', 'MQM'])
  const [cmp, setCmp] = useState(null)
  const [cmpStatus, setCmpStatus] = useState('Select at least two models to compare.')
  const [cmpBusy, setCmpBusy] = useState(false)
  const [minModels, setMinModels] = useState(2)
  const [common, setCommon] = useState(null)
  const [commonStatus, setCommonStatus] = useState('')
  const [commonBusy, setCommonBusy] = useState(false)
  const setP = (k) => (v) => setPlot((p) => ({ ...p, [k]: v }))

  tool.onResult(() => setTab('results'))
  useSignal(bridge, 'errorOccurred', (m) => tool.fail('Error: ' + m))
  useGlobals({
    showCompareResult: (r) => {
      setCmpBusy(false)
      if (!r) return
      setCmpStatus(r.message || (r.status === 'success' ? 'Model comparison completed.' : 'Comparison failed.'))
      if (r.status === 'success') setCmp(r)
    },
    showCommonResult: (r) => {
      setCommonBusy(false)
      if (!r) return
      setCommonStatus(r.status === 'error' ? r.message || 'Failed.' : `Found ${(r.table || []).length} common QTL region(s).`)
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

  const settings = () => ({
    threshold_mode: thrMode, lod_threshold: lod, n_perm: nPerm, perm_alpha: permAlpha, p_threshold: pThr,
    lod_drop: lodDrop, peak_drop: peakDrop, validation_fraction: valFrac, n_cofactors: cof, window_size: win,
  })
  const run = () => {
    if (!geno || !pheno) return tool.fail('Please select both genotype and phenotype files.')
    if (!trait) return tool.fail('Please select a phenotype/trait.')
    tool.run({ genotype_file: geno, phenotype_file: pheno, phenotype_column: trait, mapping_method: method, ...settings(), plot })
  }
  const compare = () => {
    if (models.length < 2) return setCmpStatus('Select at least two models.')
    setCmpStatus('Running model comparison...'); setCmpBusy(true)
    bridge.run_compare(JSON.stringify({
      genotype_file: geno, phenotype_file: pheno, phenotype_column: trait, models, ...settings(), plot,
    }))
  }
  const findCommon = () => { setCommonStatus('Searching…'); setCommonBusy(true); bridge.find_common(JSON.stringify({ min_models: minModels, lod_threshold: lod })) }
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
              <SelectField label="Significance threshold" value={thrMode} onChange={setThrMode}
                options={[{ value: 'lod', label: 'Fixed LOD threshold' }, { value: 'permutation', label: 'Genome-wide permutation test (SIM / CIM / MQM)' }]} />
              {thrMode === 'lod' ? (
                <Slider label="LOD Threshold" value={lod} onChange={setLod} min={0.5} max={10} step={0.1} format={(v) => v.toFixed(1)} />
              ) : (
                <div className="grid-2">
                  <TextField label="Permutations" type="number" min={100} max={10000} step={100} value={nPerm} onChange={setNPerm} />
                  <TextField label="Genome-wide α" type="number" min={0.01} max={0.2} step={0.01} value={permAlpha} onChange={setPermAlpha} />
                </div>
              )}
              {thrMode === 'permutation' && (method === 'CIM' || method === 'MQM') && (
                <Note>
                  Every permutation repeats the whole {method} analysis on the shuffled phenotype (the cofactors are selected again, then the scan), as Doerge &amp; Churchill (1996) advise;
                  a single-marker null would give a threshold that is too low. It takes longer than SIM: about 0.15 s per permutation for 250 lines and 4,000 markers.
                  RF, SVR and XGBoost use the typed LOD threshold.
                </Note>
              )}
              <div className="grid-2" style={{ marginTop: 6 }}>
                <TextField label="Peak p-value cut-off" type="number" min={0.000001} max={1} step={0.01} value={pThr} onChange={setPThr} />
                <TextField label="Peak separation (LOD drop)" type="number" min={0.5} max={5} step={0.1} value={peakDrop} onChange={setPeakDrop} />
              </div>
              <TextField label="Support interval — LOD drop (1.5 ≈ 95%)" type="number" min={0.3} max={3} step={0.1} value={lodDrop} onChange={setLodDrop} />
              {method === 'CIM' && (
                <div className="grid-2" style={{ marginTop: 6 }}>
                  <TextField label="CIM cofactors (max)" type="number" min={1} max={20} value={cof} onChange={setCof} />
                  <TextField label="CIM window (cM)" type="number" min={1} max={50} value={win} onChange={setWin} />
                </div>
              )}
              <Slider label="Hold-out fraction (0 = use all individuals)" value={valFrac} onChange={setValFrac} min={0} max={0.4} step={0.05} format={(v) => v.toFixed(2)} />
              <Note>
                A QTL is a region: one peak (local LOD maximum above the threshold and p-value cut-off) with the flanking markers
                where the LOD falls more than the “LOD drop” below the peak. Two peaks on a chromosome are separate QTL only if the
                curve dips by at least “Peak separation” between them. RF / SVR / XGBoost scores are importance values on a LOD-like
                scale, not true LOD.
              </Note>
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
              <div className="plot-title">QTL regions</div>
              <div className="info-line">
                {s.significant_qtl ?? 0} QTL · LOD threshold {s.threshold_lod != null ? Number(s.threshold_lod).toFixed(2) : '—'}
                {s.threshold_note ? ` (${s.threshold_note})` : ''} · support interval {s.lod_drop}-LOD drop
                {s.imputed_genotypes ? ` · ${s.imputed_genotypes} missing genotypes filled from flanking markers` : ''}
                {s.cofactors && s.cofactors.length ? ` · cofactors: ${s.cofactors.join(', ')}` : ''}
              </div>
              <DataTable rows={r.qtl_table} columns={QTL_COLS} maxHeight={420} />
              <div className="field-hint" style={{ marginTop: 8 }}>
                Additive effect (a) = half the difference between the two homozygote classes, in trait units; positive means the allele
                of the higher-coded parent (genotype 2 / “b”) raises the trait. For CIM / MQM it is estimated with the background
                cofactors in the model. PVE_Marker = single-marker R²; the full CSV has every column (dominance, SE, 95% CI, Bayes interval).
              </div>
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
          <Card title="Comparison Table" icon={Table2}>
            {cmp && cmp.csv && <button className="btn btn-secondary" style={{ marginBottom: 12 }} onClick={() => dl('compare_csv')}><FileText size={15} /> Download table (CSV)</button>}
            <DataTable rows={cmp && cmp.table} columns={CMP_COLS} maxHeight={380} />
          </Card>
        </TabPane>
      )}
      {tab === 'common' && (
        <TabPane id="common" key="common">
          <Card title="Common QTL Across Models" icon={Target}>
            <div style={{ maxWidth: 450 }}>
              <Slider label="Minimum Models Detecting QTL" value={minModels} onChange={setMinModels} min={2} max={6} />
            </div>
            <button className="btn btn-secondary" style={{ marginTop: 16 }} onClick={findCommon} disabled={commonBusy}><Search size={16} /> {commonBusy ? 'Searching…' : 'Find Common QTL'}</button>
            {commonStatus && <div className="info-line">{commonStatus}</div>}
          </Card>
          <div className="plot-grid-2">
            <Card title="Common QTL Ideogram" icon={MapIcon}>{common && common.ideogram ? <Pic src={common.ideogram} label="Common QTL Ideogram" /> : <div className="empty-state" style={{ padding: 20 }}>Common-QTL ideogram will appear here.</div>}</Card>
            <Card title="Common QTL Counts" icon={ChartColumn}>{common && common.bar ? <Pic src={common.bar} label="Common QTL Counts" /> : <div className="empty-state" style={{ padding: 20 }}>Common-QTL counts will appear here.</div>}</Card>
          </div>
          <Card title="Common QTL Table" icon={Table2}>
            {common && common.csv && <button className="btn btn-secondary" style={{ marginBottom: 12 }} onClick={() => dl('common_csv')}><FileText size={15} /> Download table (CSV)</button>}
            <DataTable rows={common && common.table} columns={COMMON_COLS} maxHeight={380} />
          </Card>
        </TabPane>
      )}
    </ToolPage>
  )
}
