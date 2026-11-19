import { useEffect, useState } from 'react'
import { Dna, SlidersHorizontal, LayoutGrid, FolderOpen, Info, PaintBucket, Grid3x3, ChartColumn, Table2, Download, FileText, ArrowRight, TriangleAlert } from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, SelectField, ChipGroup, RadioChips, FilePicker, Note, SubTabs,
  EmptyState, StatCard, ImageGallery, DataTable,
} from '../components/kit.jsx'
import { useBridge, useTool, mockGeneric } from '../lib/bridge.js'
import { setPipe, usePipe } from '../lib/pipeline.js'

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
const DOM_METHODS = [
  { value: 'nishio', label: 'Nishio & Satoh 2014 (sommer D.mat) — for prediction' },
  { value: 'vitezica', label: 'Vitezica et al. 2013 (orthogonal) — for variance partition' },
]
const MODES = [
  { value: 'single', label: 'One file — full matrix' },
  { value: 'split', label: 'Training + Testing files' },
]

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

export default function AdditiveDominanceTool({ navigate }) {
  const { bridge } = useBridge('additiveDominanceBridge', () => mockGeneric())
  const tool = useTool(bridge, 'run_matrices')
  const pipe = usePipe()
  const [tab, setTab] = useState('input')
  const [sub, setSub] = useState('heat')
  const [mode, setMode] = useState('single')
  const [file, setFile] = useState('')
  const [trainFile, setTrainFile] = useState('')
  const [testFile, setTestFile] = useState('')
  const [fromClu, setFromClu] = useState(false)
  const [impute, setImpute] = useState('mean')
  const [domMethod, setDomMethod] = useState('nishio')
  const [types, setTypes] = useState(['Additive', 'Dominance'])
  const [which, setWhich] = useState('A')

  // the Clustering tool has just produced a training and a testing file: use them
  useEffect(() => {
    const c = pipe.clustering
    if (c && c.train_csv && c.test_csv) { setMode('split'); setTrainFile(c.train_csv); setTestFile(c.test_csv); setFromClu(true) }
  }, [pipe.clustering])

  tool.onResult((r) => {
    setTab('results'); setSub('heat')
    setWhich(r.table_additive && r.table_additive.length ? 'A' : 'D')
  })

  const run = () => {
    if (!types.length) return tool.fail('Select at least one matrix type.')
    const m = types.map((t) => MATRICES[t])
    if (mode === 'split') {
      if (!trainFile || !testFile) return tool.fail('Please select both the training and the testing file.')
      return tool.run({ train_file: trainFile, test_file: testFile, impute_method: impute, matrix_types: m, dominance_method: domMethod })
    }
    if (!file) return tool.fail('Please select a SNP file first.')
    tool.run({ snp_file: file, impute_method: impute, matrix_types: m, dominance_method: domMethod })
  }
  const r = tool.result
  const s = (r && r.summary) || {}
  const isSplit = s.mode === 'split'
  const sp = s.split || {}
  const msg = r
    ? `Done. Built ${(s.matrix_types || []).join(' + ')} matrix for ${s.n_samples} samples × ${s.n_snps} markers.` +
      (s.mono_removed ? ` (${s.mono_removed} monomorphic marker(s) removed)` : '')
    : 'Run the analysis to see results here.'
  const rows = r && (which === 'A' ? r.table_additive : r.table_dominance)

  // hand the phenotyped training / testing files and one matrix pair on to genomic prediction
  const toPrediction = (kind) => {
    setPipe({
      matrices: {
        kind,
        train_pheno_file: trainFile, test_pheno_file: testFile,
        train_geno_file: kind === 'additive' ? r.train_additive_csv : r.train_dominance_csv,
        test_geno_file: kind === 'additive' ? r.test_additive_csv : r.test_dominance_csv,
      },
    })
    navigate && navigate('genomic_prediction')
  }
  const blocks = r ? [
    { key: 'train_additive', path: r.train_additive_csv, label: `Training additive (${sp.n_train} × ${sp.n_train})` },
    { key: 'test_additive', path: r.test_additive_csv, label: `Testing additive (${sp.n_test} × ${sp.n_train})` },
    { key: 'train_dominance', path: r.train_dominance_csv, label: `Training dominance (${sp.n_train} × ${sp.n_train})` },
    { key: 'test_dominance', path: r.test_dominance_csv, label: `Testing dominance (${sp.n_test} × ${sp.n_train})` },
  ].filter((b) => b.path) : []

  return (
    <ToolPage
      id="ad" icon={Dna} tabs={TABS} tab={tab} setTab={setTab}
      title="Genomic Relationship Matrix Builder"
      subtitle="Additive (VanRaden / rrBLUP A.mat) and dominance (Nishio & Satoh or Vitezica) genomic relationship matrices from SNP dosage data"
      run={{ label: 'Compute Matrices', runningLabel: 'Computing…', state: tool.state, status: tool.status, onRun: run }}
    >
      {tab === 'input' && (
        <TabPane id="input" key="input">
          <SectionHeader title="What to compute">
            One matrix of a single file, or the training and testing matrices of the two files made by the Clustering tool.
          </SectionHeader>
          <Card title="Mode" icon={Grid3x3}>
            <RadioChips options={MODES} value={mode} onChange={setMode} />
          </Card>
          {mode === 'single' ? (
            <>
              <SectionHeader title="SNP Data" />
              <Card title="SNP File" icon={FolderOpen}>
                <FilePicker value={file} placeholder="Select a CSV or TSV file…"
                  onBrowse={async () => { const p = await tool.browse('browse_snp_file'); if (p) setFile(p) }} />
                <Note>
                  Rows = samples, columns = SNP markers • First column = sample IDs<br />
                  Marker encoding: 0 (AA) / 1 (Aa) / 2 (aa) • Missing values (NA, blank or -1) are imputed automatically
                </Note>
              </Card>
            </>
          ) : (
            <>
              <SectionHeader title="Training and Testing files">
                The two files written by Clustering &amp; Genomic Selection (the phenotypes column and the Cluster column are not treated as markers).
              </SectionHeader>
              {fromClu && pipe.clustering && trainFile === pipe.clustering.train_csv && testFile === pipe.clustering.test_csv && (
                <Note>Filled in from the Clustering run you have just made. Choose other files below if you want different ones.</Note>
              )}
              <div className="row">
                <Card title="Training file" icon={FolderOpen} delay={0.04}>
                  <FilePicker value={trainFile} placeholder="Training_data_Cluster.csv"
                    onBrowse={async () => { const p = await tool.browse('browse_train_file'); if (p) { setTrainFile(p); setFromClu(false) } }} />
                </Card>
                <Card title="Testing file" icon={FolderOpen} delay={0.08}>
                  <FilePicker value={testFile} placeholder="Testing_data_Cluster.csv"
                    onBrowse={async () => { const p = await tool.browse('browse_test_file'); if (p) { setTestFile(p); setFromClu(false) } }} />
                </Card>
              </div>
              <Note>
                The allele frequencies, the centring, the imputation and the monomorphic-marker filter all come from the <b>training lines</b> only, for the
                additive and the dominance matrix, so the testing lines never influence the training matrix. The additive matrix is scaled as in VanRaden (2008),
                <code>X<sub>c</sub>X<sub>c</sub>′ / 2Σp(1−p)</code> (diagonal about 1, as <code>rrBLUP::A.mat</code>), and written as <b>training × training</b> and
                <b>testing × training</b> — the training lines are the columns of both files. Missing calls (-1 or blank) are imputed from the training lines.
              </Note>
            </>
          )}
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
          {types.includes('Dominance') && (
            <Card title="Dominance matrix" icon={Grid3x3} delay={0.14}>
              <SelectField label="Dominance relationship matrix" value={domMethod} onChange={setDomMethod} options={DOM_METHODS} />
              <Note>
                <b>Nishio &amp; Satoh</b> centres the heterozygosity of every marker; it is the matrix of <code>sommer::D.mat</code> and works for genomic
                prediction. It is not orthogonal to the additive effects, so variance components estimated with it put too much into dominance and too little into
                additive (Vitezica et al. 2013). <b>Vitezica</b> uses the dominance deviation of each genotype, orthogonal to the additive effect, and is the one to
                use when you want to partition the variance. Both assume genotype frequencies close to Hardy–Weinberg, and neither is informative for RIL or DH lines.
              </Note>
            </Card>
          )}
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
                {isSplit && <StatCard className="summary-card" value={sp.n_train ?? 0} label="Training lines" delay={0.24} />}
                {isSplit && <StatCard className="summary-card" value={sp.n_test ?? 0} label="Testing lines" delay={0.3} />}
              </div>
              {(Array.isArray(s.warnings) ? s.warnings : []).filter((w) => typeof w === 'string').map((w) => (
                <div className="warning-box" key={w} style={{ marginBottom: 12 }}><TriangleAlert size={18} /><span>{w}</span></div>
              ))}
              {typeof s.dominance_method === 'string' && s.dominance_method && <Note>Additive: {s.additive_method}. Dominance: {s.dominance_method}. Heterozygous calls: {s.het_rate}%.</Note>}
              {isSplit && (
                <div className="note" style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', gap: 16, flexWrap: 'wrap' }}>
                  <span><b>Next step:</b> genomic prediction, with the training and testing files and this matrix. They are filled in for you.</span>
                  <span style={{ display: 'flex', gap: 10, flexWrap: 'wrap' }}>
                    {r.train_additive_csv && <button className="btn btn-primary" onClick={() => toPrediction('additive')}>Predict with additive <ArrowRight size={16} /></button>}
                    {r.train_dominance_csv && <button className="btn btn-primary" onClick={() => toPrediction('dominance')}>Predict with dominance <ArrowRight size={16} /></button>}
                  </span>
                </div>
              )}
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
                  {isSplit
                    ? blocks.map((b) => (
                      <button key={b.key} className="btn btn-secondary" onClick={() => bridge.download_result(b.key)}><FileText size={16} /> {b.label}</button>
                    ))
                    : (
                      <>
                        {r.additive_csv && <button className="btn btn-secondary" onClick={() => bridge.download_result('additive')}><Download size={16} /> Additive Matrix (CSV)</button>}
                        {r.dominance_csv && <button className="btn btn-secondary" onClick={() => bridge.download_result('dominance')}><Download size={16} /> Dominance Matrix (CSV)</button>}
                      </>
                    )}
                  {!r.additive_csv && !r.dominance_csv && !blocks.length && <EmptyState icon={Download}>Nothing to download.</EmptyState>}
                </div>
              )}
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
