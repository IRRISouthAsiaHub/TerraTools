import { useState } from 'react'
import {
  Sprout, SlidersHorizontal, LayoutGrid, ChartColumn, CircleHelp, FolderOpen, Settings2, Terminal, Download,
  Trophy, Star, FlaskConical, ChartLine, FileText, BookOpen,
} from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, TextField, FileWithPreview, Note, EmptyState,
  StatCard, ImageGallery, DataTable, LogBox, ProgressBar,
} from '../components/kit.jsx'
import { useBridge, useTool, useSignal, callJson, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'setup', label: 'Setup & Run', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
  { id: 'plots', label: 'Plots', icon: ChartColumn },
  { id: 'help', label: 'Help', icon: CircleHelp },
]

const FILES = [
  { key: 'train_pheno_file', label: 'Training Phenotype CSV (col 0 = ID, col 1 = trait value)', req: true },
  { key: 'train_geno_file', label: 'Training Genotype CSV (row index = sample ID)', req: true },
  { key: 'parent_geno_file', label: 'Parent Genotype CSV (row index = parent ID)', req: true },
  { key: 'test_geno_file', label: 'Test Genotype CSV (optional)', req: false },
]

export default function HeterosisTool() {
  const { bridge } = useBridge('heterosisBridge', () => mockGeneric(['progressUpdate', 'errorOccurred']))
  const tool = useTool(bridge, 'run_heterosis')
  const [tab, setTab] = useState('setup')
  const [paths, setPaths] = useState({})
  const [previews, setPreviews] = useState({})
  const [c, setC] = useState({ num_F1: 5, num_P: 5, k_folds: 5, epochs: 100, batch_size: 16, dropout_rate: 0.2, l1_reg: 0.001, l2_reg: 0.001 })
  const [progress, setProgress] = useState(0)
  const [log, setLog] = useState('')
  const set = (k) => (v) => setC((p) => ({ ...p, [k]: v }))

  useSignal(bridge, 'progressUpdate', (fraction, message) => {
    if (fraction !== null && fraction >= 0) setProgress(fraction)
    if (message) setLog((l) => (l ? l + '\n' : '') + message)
  })
  useSignal(bridge, 'errorOccurred', (m) => tool.fail('Error: ' + m))
  tool.onResult((r) => {
    setProgress(1)
    if (r.logs) setLog(r.logs)
    setTab('results')
  })

  const pick = (key) => async () => {
    const path = await tool.browse('browse_file', key)
    if (!path) return
    setPaths((x) => ({ ...x, [key]: path }))
    setPreviews((x) => ({ ...x, [key]: 'Loading preview…' }))
    const pv = await callJson(bridge, 'get_file_preview', path)
    setPreviews((x) => ({ ...x, [key]: pv && pv.status === 'success' ? `${pv.n_rows} rows × ${pv.n_cols} columns` : `Could not preview file: ${(pv && pv.message) || ''}` }))
  }
  const ready = paths.train_pheno_file && paths.train_geno_file && paths.parent_geno_file

  const run = () => {
    setLog(''); setProgress(0)
    tool.run(c)
  }
  const r = tool.result
  const cv = (r && r.cv_summary) || {}

  return (
    <ToolPage
      id="het" icon={Sprout} tabs={TABS} tab={tab} setTab={setTab}
      title="Heterosis Prediction"
      subtitle="Deep learning GEBV prediction & hybrid screening"
      run={tab === 'setup' ? {
        label: 'Run Prediction', runningLabel: 'Training…', state: tool.state, disabled: !ready,
        status: !ready && tool.state === 'idle' ? 'Upload the three required files first.' : tool.state === 'idle' ? 'Ready. Adjust parameters and click Run Prediction.' : tool.status,
        onRun: run,
      } : undefined}
    >
      {tab === 'setup' && (
        <TabPane id="setup" key="setup">
          <div className="grid-2" style={{ gridTemplateColumns: 'repeat(auto-fit, minmax(380px, 1fr))', alignItems: 'start' }}>
            <Card title="Input Files" icon={FolderOpen}>
              {FILES.map((f) => (
                <FileWithPreview key={f.key} label={f.label} value={paths[f.key]} placeholder="Select file…" onBrowse={pick(f.key)} preview={previews[f.key]} />
              ))}
              <Note>Genotype files are numeric marker matrices (samples × markers). Parent and training genotype files must share the same marker columns.</Note>
            </Card>
            <Card title="Parameters" icon={Settings2} delay={0.06}>
              <Slider label="Top F1 hybrids to report" value={c.num_F1} onChange={set('num_F1')} min={1} max={50} />
              <Slider label="Top parents to report" value={c.num_P} onChange={set('num_P')} min={1} max={30} />
              <Slider label="K-Fold CV splits" value={c.k_folds} onChange={set('k_folds')} min={2} max={10} />
              <Slider label="Epochs" value={c.epochs} onChange={set('epochs')} min={10} max={500} step={10} />
              <Slider label="Batch size" value={c.batch_size} onChange={set('batch_size')} min={8} max={128} step={8} />
              <Slider label="Dropout rate" value={c.dropout_rate} onChange={set('dropout_rate')} min={0} max={0.5} step={0.05} format={(v) => v.toFixed(2)} />
              <div className="grid-2" style={{ marginTop: 6 }}>
                <TextField label="L1 regularization" type="number" step={0.0001} value={c.l1_reg} onChange={set('l1_reg')} />
                <TextField label="L2 regularization" type="number" step={0.0001} value={c.l2_reg} onChange={set('l2_reg')} />
              </div>
            </Card>
          </div>
          <Card title="Training Log" icon={Terminal} delay={0.12}>
            <ProgressBar value={progress} />
            <LogBox>{log || 'Logs will appear here after training…'}</LogBox>
          </Card>
        </TabPane>
      )}
      {tab === 'results' && (
        <TabPane id="results" key="results">
          <SectionHeader title="Results" />
          {!r ? <EmptyState icon={FlaskConical}>Run a prediction to see results.</EmptyState> : (
            <>
              <Card title="Downloads" icon={Download}>
                <div className="downloads" style={{ margin: 0 }}>
                  {[['f1_hybrids', 'Top F1 Hybrids'], ['top_parents', 'Top Parents'], ['all_hybrids', 'All Hybrid Results']].map(([k, l]) => (
                    <button key={k} className="btn btn-secondary" onClick={() => bridge.download_result(k)}><FileText size={15} /> {l}</button>
                  ))}
                  <button className="btn btn-secondary" disabled={!paths.test_geno_file} onClick={() => bridge.download_result('test_predictions')}><FileText size={15} /> Test Predictions</button>
                </div>
              </Card>
              <Card title="Cross-Validation Summary" icon={ChartLine}>
                <div className="summary-grid">
                  <StatCard className="summary-card" value={cv.mse ?? '—'} label="MSE" />
                  <StatCard className="summary-card" value={cv.rmse ?? '—'} label="RMSE" delay={0.05} />
                  <StatCard className="summary-card" value={cv.r2 ?? '—'} label="R²" delay={0.1} />
                  <StatCard className="summary-card" value={cv.pearson_r ?? '—'} label="Pearson r" delay={0.15} />
                </div>
                <div className="plot-title">Per-Fold Metrics</div>
                <DataTable rows={r.fold_metrics} columns={['Fold', 'MSE', 'RMSE', 'R2', 'Pearson_r']} maxHeight={260} />
              </Card>
              <Card title="Top F1 Hybrid Combinations" icon={Trophy}>
                <DataTable rows={r.f1_hybrids} columns={r.f1_hybrids_columns} maxHeight={340} />
              </Card>
              <Card title="Top Parent Lines by Mean GEBV" icon={Star}>
                <DataTable rows={r.top_parents} columns={r.top_parents_columns} maxHeight={340} />
              </Card>
              <Card title="Test Genotype Predictions" icon={FlaskConical}>
                {r.test_predictions && r.test_predictions.length
                  ? <DataTable rows={r.test_predictions} columns={r.test_predictions_columns} maxHeight={340} />
                  : <div className="empty-state" style={{ padding: 24 }}>No test genotype file was provided.</div>}
              </Card>
            </>
          )}
        </TabPane>
      )}
      {tab === 'plots' && (
        <TabPane id="plots" key="plots">
          <SectionHeader title="Hybrid Performance Visualization" />
          {r && r.plot_path ? (
            <>
              <ImageGallery single images={[{ src: r.plot_path, label: 'Hybrid Performance' }]} />
              <button className="btn btn-secondary" onClick={() => bridge.download_result('plot')}><Download size={16} /> Download Plot PNG</button>
            </>
          ) : <EmptyState icon={ChartColumn}>Run a prediction first.</EmptyState>}
        </TabPane>
      )}
      {tab === 'help' && (
        <TabPane id="help" key="help">
          <SectionHeader title="Application Guide" />
          <Card title="Required Input Format" icon={BookOpen}>
            <ul className="plain-list">
              <li>Training Phenotype CSV — column 0: sample ID, column 1: numeric trait value</li>
              <li>Training Genotype CSV — row index: sample ID, columns: numeric marker calls</li>
              <li>Parent Genotype CSV — row index: parent ID, columns: same markers as training genotype</li>
              <li>Test Genotype CSV (optional) — same marker columns, for out-of-sample GEBV prediction</li>
            </ul>
          </Card>
          <Card title="Workflow" icon={BookOpen} delay={0.06}>
            <ol className="plain-list">
              <li>Upload the three required CSV files (and an optional test genotype file) in Setup & Run.</li>
              <li>Adjust cross-validation and network hyperparameters as needed.</li>
              <li>Click Run Prediction — this trains K models for cross-validation, then a final model on the full dataset.</li>
              <li>All pairwise parent combinations are scored to identify the top-performing F1 hybrids.</li>
              <li>Review CV metrics, hybrid/parent rankings, and predictions in Results; export CSVs as needed.</li>
              <li>Inspect the hybrid ranking and GEBV distribution charts in Plots.</li>
            </ol>
          </Card>
          <Card title="Model & Metrics" icon={BookOpen} delay={0.12}>
            <Note><b>Model architecture:</b> 8-layer DNN (512→512→256→128→64→32→16→8) with batch normalization, LeakyReLU, dropout, and L1/L2 regularization.</Note>
            <Note><b>Metrics:</b> MSE · RMSE · R² · Pearson r, evaluated via K-fold cross-validation.</Note>
            <Note><b>Heterosis indices:</b> Mid-Parent Heterosis (MPH) and Better-Parent Heterosis (BPH), both computed relative to the training population.</Note>
          </Card>
        </TabPane>
      )}
    </ToolPage>
  )
}
