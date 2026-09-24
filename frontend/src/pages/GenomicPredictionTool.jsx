import { useState } from 'react'
import {
  Eye, SlidersHorizontal, LayoutGrid, FileText, Dna, Layers, Brain, TreePine, Settings2, Terminal,
  Table2, ListOrdered, ChartLine, Download, TriangleAlert, Info, Sprout,
} from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, Toggle, TextField, FilePicker, Note, SubTabs,
  EmptyState, ImageGallery, DataTable, DownloadButtons, LogBox,
} from '../components/kit.jsx'
import { useBridge, useTool, useSignal, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'input', label: 'Input & Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]
const SUB = [
  { id: 'metrics', label: 'Metrics & Interpretation', icon: Table2 },
  { id: 'preds', label: 'Test Predictions', icon: ListOrdered },
  { id: 'plots', label: 'Plots', icon: ChartLine },
  { id: 'dl', label: 'Downloads', icon: Download },
]
const PLOT_LABELS = {
  regression: 'Best-Fold Regression + Residuals + Distribution (BreedSight)',
  comparison: 'Train vs Validation R² and Pearson r — All Models',
  folds: 'R² and Pearson r per Fold (95% CI)',
  heatmap: 'Performance Heatmap (RMSE & R² × Folds)',
  test_grid: 'Test Predictions — All Models',
  distribution: 'Phenotype Distribution — True vs Predicted',
}
const FILES = [
  { key: 'train_pheno_file', slot: 'browse_train_pheno', title: 'Training Phenotypic Data', icon: FileText, ph: 'sample_id, phenotypes' },
  { key: 'train_geno_file', slot: 'browse_train_geno', title: 'Training Genotypic Data', icon: Dna, ph: 'sample_id, SNP1, SNP2, …' },
  { key: 'test_pheno_file', slot: 'browse_test_pheno', title: 'Testing Phenotypic Data', icon: FileText, ph: 'sample_id, [phenotypes optional]' },
  { key: 'test_geno_file', slot: 'browse_test_geno', title: 'Testing Genotypic Data', icon: Dna, ph: 'sample_id, SNP1, SNP2, …' },
]
const DEFAULTS = {
  outer_n_splits: 5, heritability: 0.82, epochs: 300, batch_size: 32, learning_rate: 0.001, alpha: 0.5,
  l2_reg: 0.01, dropout_rate: 0.3, rf_n_estimators: 200, rf_max_depth: 15,
  use_pca: true, feature_selection: true, use_raw_genotypes: false, n_components: 50, rfe_n_features: 100,
}

/** Tiny **bold** + line-break renderer (the interpretation text is simple markdown). */
function MdLite({ text }) {
  if (!text) return null
  return (
    <div className="interp-box">
      {text.split('\n').map((line, i) => (
        <p key={i} style={{ minHeight: line ? undefined : '0.7em' }}>
          {line.split(/\*\*(.+?)\*\*/g).map((seg, j) => (j % 2 ? <strong key={j}>{seg}</strong> : seg))}
        </p>
      ))}
    </div>
  )
}

export default function GenomicPredictionTool() {
  const { bridge } = useBridge('genomicPredictionBridge', () => mockGeneric(['logUpdate']))
  const tool = useTool(bridge, 'run_prediction')
  const [tab, setTab] = useState('input')
  const [sub, setSub] = useState('metrics')
  const [files, setFiles] = useState({})
  const [c, setC] = useState(DEFAULTS)
  const [log, setLog] = useState('')
  const set = (k) => (v) => setC((p) => ({ ...p, [k]: v }))
  useSignal(bridge, 'logUpdate', (line) => setLog((l) => l + line + '\n'))
  tool.onResult(() => { setTab('results'); setSub('metrics') })

  const run = () => {
    if (!FILES.every((f) => files[f.key])) return tool.fail('Please select all four input files.')
    setLog('')
    tool.run({ ...Object.fromEntries(FILES.map((f) => [f.key, files[f.key]])), ...c })
  }
  const r = tool.result
  const dl = (r && r.downloads) || {}
  const plots = r ? Object.entries(r.plots || {}).filter(([, p]) => p) : []

  return (
    <ToolPage
      id="gp" icon={Eye} tabs={TABS} tab={tab} setTab={setTab}
      title="BreedSight Genomic Prediction"
      subtitle="Hybrid Deep Neural Network + Random Forest ensemble for genomic selection, benchmarked against Lasso/RBLUP/GBLUP"
      run={{ label: 'Run BreedSight Analysis', runningLabel: 'Training models…', state: tool.state, status: tool.status, onRun: run }}
    >
      {tab === 'input' && (
        <TabPane id="input" key="input">
          <div className="warning-box"><TriangleAlert size={18} />
            <span>This tool requires TensorFlow and scikit-learn to be installed. First runs and large datasets can take a long time — training happens once per cross-validation fold, plus a final full retrain.</span>
          </div>
          <SectionHeader title="Data Files">Upload training and testing phenotypic/genotypic data.</SectionHeader>
          <div className="row">
            {FILES.map((f, i) => (
              <Card key={f.key} title={f.title} icon={f.icon} delay={i * 0.05}>
                <FilePicker value={files[f.key]} placeholder={f.ph}
                  onBrowse={async () => { const p = await tool.browse(f.slot); if (p) setFiles((x) => ({ ...x, [f.key]: p })) }} />
              </Card>
            ))}
          </div>
          <Card title="Expected CSV Format" icon={Info}>
            <div className="field-hint" style={{ fontSize: '.84rem', lineHeight: 1.8 }}>
              <strong>Phenotypic:</strong> sample_id, phenotypes<br />
              <strong>Genotypic:</strong> sample_id, SNP1, SNP2, …<br />
              Marker encoding: 0 (AA) / 1 (Aa) / 2 (aa) • No missing values — impute before upload
            </div>
          </Card>
          <SectionHeader title="Model Parameters" />
          <div className="row">
            <Card title="Cross-Validation & Heritability" icon={Layers}>
              <Slider label="Number of CV folds" value={c.outer_n_splits} onChange={set('outer_n_splits')} min={2} max={10} />
              <Slider label="Trait heritability h²" value={c.heritability} onChange={set('heritability')} min={0.1} max={1} step={0.01} format={(v) => v.toFixed(2)} />
            </Card>
            <Card title="BreedSight — DNN + RF Ensemble" icon={Brain} delay={0.05}>
              <Slider label="Training epochs" value={c.epochs} onChange={set('epochs')} min={50} max={500} step={50} />
              <Slider label="Batch size" value={c.batch_size} onChange={set('batch_size')} min={16} max={128} step={16} />
              <TextField label="Learning rate" type="number" step={0.0001} value={c.learning_rate} onChange={set('learning_rate')} />
              <Slider label="DNN weight α (RF weight = 1−α)" value={c.alpha} onChange={set('alpha')} min={0} max={1} step={0.05} format={(v) => v.toFixed(2)} />
            </Card>
            <Card title="Regularization & RF" icon={TreePine} delay={0.1}>
              <Slider label="L2 regularization" value={c.l2_reg} onChange={set('l2_reg')} min={0} max={0.1} step={0.005} format={(v) => v.toFixed(3)} />
              <Slider label="Dropout rate" value={c.dropout_rate} onChange={set('dropout_rate')} min={0} max={0.7} step={0.05} format={(v) => v.toFixed(2)} />
              <Slider label="RF trees" value={c.rf_n_estimators} onChange={set('rf_n_estimators')} min={50} max={500} step={50} />
              <Slider label="RF max depth (0 = None)" value={c.rf_max_depth} onChange={set('rf_max_depth')} min={0} max={30} />
            </Card>
          </div>
          <Card title="Feature Engineering" icon={Settings2}>
            <div className="row" style={{ gap: 36 }}>
              <div style={{ flex: 1, minWidth: 260 }}>
                <Toggle label="Use PCA on GRM features" checked={c.use_pca} onChange={set('use_pca')} />
                <Toggle label="RFE feature selection" checked={c.feature_selection} onChange={set('feature_selection')} />
                <Toggle label="Use raw genotypes (skip GRM)" checked={c.use_raw_genotypes} onChange={set('use_raw_genotypes')} />
              </div>
              <div style={{ flex: 1, minWidth: 260 }}>
                <Slider label="PCA components" value={c.n_components} onChange={set('n_components')} min={10} max={200} step={10} />
                <Slider label="RFE features to select" value={c.rfe_n_features} onChange={set('rfe_n_features')} min={20} max={500} step={10} />
              </div>
            </div>
          </Card>
          <Card title="Live Log" icon={Terminal}>
            <LogBox>{log || 'Waiting to start…'}</LogBox>
          </Card>
        </TabPane>
      )}
      {tab === 'results' && (
        <TabPane id="results" key="results">
          <SectionHeader title="Results">
            {r ? `Best model: ${r.best_model} — Val R² = ${Number(r.best_val_r2).toFixed(4)}` : 'Run the analysis to see results here.'}
          </SectionHeader>
          {!r ? <EmptyState icon={Sprout}>Run the analysis to see results here.</EmptyState> : (
            <>
              <SubTabs tabs={SUB} value={sub} onChange={setSub} id="gpsub" />
              {sub === 'metrics' && (
                <>
                  <DataTable rows={r.metrics_table} maxHeight={360} />
                  <div style={{ height: 16 }} />
                  <MdLite text={r.interpretation_model} />
                  <MdLite text={r.interpretation_heritability} />
                </>
              )}
              {sub === 'preds' && <DataTable rows={r.test_predictions_table} maxHeight={520} />}
              {sub === 'plots' && (plots.length
                ? <ImageGallery images={plots.map(([k, p]) => ({ src: p, label: PLOT_LABELS[k] || k }))} />
                : <EmptyState icon={Info}>No plots were generated.</EmptyState>)}
              {sub === 'dl' && (
                <DownloadButtons onSave={tool.save} items={[
                  dl.metrics_csv && { icon: FileText, label: 'Model Metrics (CSV)', path: dl.metrics_csv, name: 'model_metrics.csv' },
                  dl.test_csv && { icon: FileText, label: 'Test Predictions (CSV)', path: dl.test_csv, name: 'test_predictions.csv' },
                  dl.val_csv && { icon: FileText, label: 'Validation Predictions (CSV)', path: dl.val_csv, name: 'validation_predictions.csv' },
                ].filter(Boolean)} />
              )}
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
