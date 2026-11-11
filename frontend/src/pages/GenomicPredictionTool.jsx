import { useEffect, useState } from 'react'
import {
  Eye, SlidersHorizontal, LayoutGrid, FileText, Dna, Layers, Brain, TreePine, Settings2, Terminal,
  Table2, ListOrdered, Download, Info, Sprout,
} from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, Toggle, TextField, FilePicker, Note, SubTabs,
  EmptyState, ImageGallery, DataTable, DownloadButtons, LogBox, ChipGroup, SummaryGrid,
} from '../components/kit.jsx'
import { useBridge, useTool, useSignal, mockGeneric } from '../lib/bridge.js'
import { usePipe } from '../lib/pipeline.js'

const TABS = [
  { id: 'input', label: 'Input & Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]
const SUB = [
  { id: 'summary', label: 'Summary', icon: Table2 },
  { id: 'preds', label: 'Test Predictions', icon: ListOrdered },
  { id: 'dl', label: 'Downloads', icon: Download },
]
// the validation plots shown under the summary
const VALIDATION_PLOTS = [
  ['regression', 'Validation, best fold: true vs predicted, residuals and distribution (BreedSight)'],
  ['folds', 'Validation R² and Pearson r per fold'],
  ['comparison', 'Training vs validation R² and Pearson r'],
]
const FILES = [
  { key: 'train_pheno_file', slot: 'browse_train_pheno', title: 'Training Phenotypic Data', icon: FileText, ph: 'sample_id, phenotypes' },
  { key: 'train_geno_file', slot: 'browse_train_geno', title: 'Training Additive Matrix', icon: Dna, ph: 'training_additive.csv: sample_id + one column per training line' },
  { key: 'test_pheno_file', slot: 'browse_test_pheno', title: 'Testing Phenotypic Data', icon: FileText, ph: 'sample_id, [phenotypes optional]' },
  { key: 'test_geno_file', slot: 'browse_test_geno', title: 'Testing Additive Matrix', icon: Dna, ph: 'test_additive.csv: sample_id + the same training-line columns' },
]
// the built-in comparison models (BreedSight alone is run unless some are ticked)
const COMPARE = ['GBLUP', 'RBLUP', 'Lasso']
const DEFAULTS = {
  outer_n_splits: 5, heritability: 0.82, epochs: 300, batch_size: 32, learning_rate: 0.001, alpha: 0.5,
  l2_reg: 0.01, dropout_rate: 0.3, rf_n_estimators: 200, rf_max_depth: 15,
  use_pca: true, feature_selection: true, use_raw_genotypes: false, n_components: 50, rfe_n_features: 100,
}

const num = (v, d = 3) => (v === null || v === undefined || Number.isNaN(Number(v)) ? '—' : Number(v).toFixed(d))

/** Small read-only table (no search box): the summary has only a few rows. */
function SimpleTable({ head, rows }) {
  return (
    <div className="table-wrap">
      <table className="results-table">
        <thead><tr>{head.map((h) => <th key={h}>{h}</th>)}</tr></thead>
        <tbody>{rows.map((row, i) => <tr key={i}>{row.map((c, j) => <td key={j}>{c}</td>)}</tr>)}</tbody>
      </table>
    </div>
  )
}

export default function GenomicPredictionTool() {
  const { bridge } = useBridge('genomicPredictionBridge', () => mockGeneric(['logUpdate']))
  const tool = useTool(bridge, 'run_prediction')
  const [tab, setTab] = useState('input')
  const [sub, setSub] = useState('summary')
  const [files, setFiles] = useState({})
  const [c, setC] = useState(DEFAULTS)
  const [log, setLog] = useState('')
  const [compare, setCompare] = useState([])
  const pipe = usePipe()
  const [fromTools, setFromTools] = useState('')
  // the Additive Dominance tool has just built the matrices for this analysis
  useEffect(() => {
    const m = pipe.matrices
    if (!m) return
    setFiles({ train_pheno_file: m.train_pheno_file, train_geno_file: m.train_geno_file, test_pheno_file: m.test_pheno_file, test_geno_file: m.test_geno_file })
    setFromTools(m.kind)
  }, [pipe.matrices])
  const set = (k) => (v) => setC((p) => ({ ...p, [k]: v }))
  useSignal(bridge, 'logUpdate', (line) => setLog((l) => l + line + '\n'))
  tool.onResult(() => { setTab('results'); setSub('summary') })

  const run = () => {
    if (!FILES.every((f) => files[f.key])) return tool.fail('Please select all four input files.')
    setLog('')
    tool.run({ ...Object.fromEntries(FILES.map((f) => [f.key, files[f.key]])), ...c, baselines: compare })
  }
  const r = tool.result
  const dl = (r && r.downloads) || {}
  const plots = (r && r.plots) || {}
  const metrics = (r && r.metrics_table) || []
  const mine = metrics.find((m) => m.Model === 'BreedSight') || metrics[0]
  const valPlots = VALIDATION_PLOTS.filter(([k]) => plots[k]).map(([k, label]) => ({ src: plots[k], label }))

  return (
    <ToolPage
      id="gp" icon={Eye} tabs={TABS} tab={tab} setTab={setTab}
      title="BreedSight Genomic Prediction"
      subtitle="Hybrid Deep Neural Network + Random Forest ensemble for genomic selection, optionally compared with GBLUP, RBLUP and Lasso"
      run={{ label: 'Run BreedSight Analysis', runningLabel: 'Training models…', state: tool.state, status: tool.status, onRun: run }}
    >
      {tab === 'input' && (
        <TabPane id="input" key="input">
          <SectionHeader title="Data Files">The phenotypic files and the additive (relationship) matrices of the training and testing lines.</SectionHeader>
          {fromTools && (
            <Note>
              The four files are filled in from <b>Additive Dominance</b> ({fromTools} matrix), which used the training and testing files of Clustering.
              Choose another file in any box to change it.
            </Note>
          )}
          <div className="row">
            {FILES.map((f, i) => (
              <Card key={f.key} title={f.title} icon={f.icon} delay={i * 0.05}>
                <FilePicker value={files[f.key]} placeholder={f.ph}
                  onBrowse={async () => { const p = await tool.browse(f.slot); if (p) { setFiles((x) => ({ ...x, [f.key]: p })); setFromTools('') } }} />
              </Card>
            ))}
          </div>
          <Card title="Expected CSV Format" icon={Info}>
            <div className="field-hint" style={{ fontSize: '.84rem', lineHeight: 1.8 }}>
              <strong>Phenotypic:</strong> sample_id, phenotypes (any further columns, such as markers or Cluster, are ignored; the testing file may leave phenotypes out)<br />
              <strong>Training additive matrix:</strong> sample_id, then one column per training line — rows = training lines, columns = the same training lines<br />
              <strong>Testing additive matrix:</strong> sample_id, then the <em>same</em> training-line columns, in the same order — rows = testing lines<br />
              Rows of each phenotypic file and its matrix must be the same lines in the same order • No missing values<br />
              Build all four with <strong>Clustering &amp; Genomic Selection</strong> (training / testing files, with your phenotype file) and
              <strong> Additive Dominance</strong> (Training + Testing mode), which compute the matrix once over all lines.
              With “Use raw genotypes” on, the matrix files hold marker genotypes instead, with the same markers in both.
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
          <Card title="Compare with other models (optional)" icon={Layers}>
            <ChipGroup options={COMPARE} value={compare} onChange={setCompare} />
            <Note>
              BreedSight alone is run unless you tick a model here. Each ticked model is trained on the same cross-validation folds,
              so its accuracy can be put side by side with BreedSight; it also adds a prediction column for the testing lines.
              GBLUP (RBF kernel ridge), RBLUP (ridge) and Lasso use fixed penalties on the BreedSight features.
            </Note>
          </Card>
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
            {r ? `${mine ? mine.N_Folds : ''}-fold cross-validation on the training lines${r.test_metrics_table ? ', and the hold-out testing lines' : ''}.` : 'Run the analysis to see results here.'}
          </SectionHeader>
          {!r ? <EmptyState icon={Sprout}>Run the analysis to see results here.</EmptyState> : (
            <>
              <SubTabs tabs={SUB} value={sub} onChange={setSub} id="gpsub" />
              {sub === 'summary' && mine && (
                <>
                  <SummaryGrid items={[
                    { label: 'Validation R²', value: num(mine.Avg_Val_R2) },
                    { label: 'Validation Pearson r', value: num(mine.Avg_Val_Pearson_r) },
                    { label: 'Validation RMSE', value: num(mine.Avg_Val_RMSE, 2) },
                    { label: 'Training R²', value: num(mine.Avg_Train_R2) },
                  ]} />
                  <h4 style={{ margin: '16px 0 8px' }}>Cross-validation (validation folds)</h4>
                  <SimpleTable
                    head={['Model', 'Val R²', 'SD', 'Val Pearson r', 'Val RMSE', 'Train R²', 'Train Pearson r']}
                    rows={metrics.map((m) => [m.Model, num(m.Avg_Val_R2), num(m.Std_Val_R2), num(m.Avg_Val_Pearson_r), num(m.Avg_Val_RMSE, 2), num(m.Avg_Train_R2), num(m.Avg_Train_Pearson_r)])} />
                  {r.test_metrics_table && r.test_metrics_table.length > 0 && (
                    <>
                      <h4 style={{ margin: '16px 0 8px' }}>Hold-out testing lines (final models trained on all training lines)</h4>
                      <SimpleTable
                        head={['Model', 'Test R²', 'Test Pearson r', 'Test RMSE', 'Lines']}
                        rows={r.test_metrics_table.map((m) => [m.Model, num(m.Test_R2), num(m.Test_Pearson_r), num(m.Test_RMSE, 2), m.N_test])} />
                    </>
                  )}
                  {valPlots.length > 0 && (
                    <>
                      <h4 style={{ margin: '20px 0 8px' }}>Validation plots</h4>
                      <ImageGallery images={valPlots} />
                    </>
                  )}
                </>
              )}
              {sub === 'preds' && (
                <>
                  <DataTable rows={r.test_predictions_table} maxHeight={520} />
                  {plots.test_grid && <ImageGallery single images={[{ src: plots.test_grid, label: 'Testing lines: true vs predicted' }]} />}
                </>
              )}
              {sub === 'dl' && (
                <DownloadButtons onSave={tool.save} items={[
                  dl.metrics_csv && { icon: FileText, label: 'Model Metrics (CSV)', path: dl.metrics_csv, name: 'model_metrics.csv' },
                  dl.test_csv && { icon: FileText, label: 'Test Predictions (CSV)', path: dl.test_csv, name: 'test_predictions.csv' },
                  dl.val_csv && { icon: FileText, label: 'Validation Predictions (CSV)', path: dl.val_csv, name: 'validation_predictions.csv' },
                  dl.test_metrics_csv && { icon: FileText, label: 'Test-Set Metrics (CSV)', path: dl.test_metrics_csv, name: 'test_metrics.csv' },
                  dl.summary_txt && { icon: FileText, label: 'Analysis Summary (TXT)', path: dl.summary_txt, name: 'analysis_summary.txt' },
                ].filter(Boolean)} />
              )}
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
