import { useState } from 'react'
import { Group, SlidersHorizontal, LayoutGrid, Dna, Percent, Download, FileText, Sprout, ArrowRight, TriangleAlert } from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, Toggle, FilePicker, Note,
  EmptyState, SummaryGrid, ImageGallery, DownloadButtons, DataTable,
} from '../components/kit.jsx'
import { useBridge, useTool, mockToolBridge } from '../lib/bridge.js'
import { setPipe } from '../lib/pipeline.js'

const TABS = [
  { id: 'input', label: 'Input & Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]

export default function ClusteringTool({ navigate }) {
  const { bridge } = useBridge('clusteringBridge', () => mockToolBridge('clustering'))
  const tool = useTool(bridge, 'run_clustering')
  const [tab, setTab] = useState('input')
  const [file, setFile] = useState('')
  const [pheno, setPheno] = useState('')
  const [keepTest, setKeepTest] = useState(false)
  const [k, setK] = useState(5)
  const [testSize, setTestSize] = useState(0.2)

  // the files go on to the next tool (Additive Dominance)
  tool.onResult((res) => {
    setTab('results')
    setPipe({ clustering: { train_csv: res.train_csv, test_csv: res.test_csv, train_has_pheno: res.train_has_pheno, test_has_pheno: res.test_has_pheno } })
  })

  const browse = async () => { const p = await tool.browse(); if (p) setFile(p) }
  const browsePheno = async () => { const p = await tool.browse('browse_pheno_file'); if (p) setPheno(p) }
  const run = () => {
    if (!file) return tool.fail('Please select an input file.')
    tool.run({ input_file: file, n_clusters: k, test_size: testSize, pheno_file: pheno, keep_test_pheno: keepTest })
  }

  const r = tool.result

  return (
    <ToolPage
      id="clu" icon={Group} tabs={TABS} tab={tab} setTab={setTab}
      title="Genotypic Clustering & Genomic Selection"
      subtitle="PCA + K-Means clustering with stratified train/test split for genomic selection"
      run={{ label: 'Run Clustering & Split', runningLabel: 'Clustering…', state: tool.state, status: tool.status, onRun: run }}
    >
      {tab === 'input' && (
        <TabPane id="input" key="input">
          <SectionHeader title="Genotype Data">
            Upload SNP genotype data (CSV or TSV). Missing values coded as -1 are handled automatically.
          </SectionHeader>
          <Card title="Genotype File" icon={Dna}>
            <FilePicker value={file} onBrowse={browse} placeholder="Select CSV/TSV (sample IDs in first column, SNP columns after)" />
          </Card>
          <Card title="Phenotype File (optional)" icon={Sprout} delay={0.04}>
            <FilePicker value={pheno} onBrowse={browsePheno} placeholder="Select CSV / TSV / Excel (sample IDs in first column, trait in the second)" />
            <Note>
              With a phenotype file the <b>training</b> file also gets a <code>phenotypes</code> column, so it can be used as a phenotyped training set.
              Lines are matched by sample ID. A genotyped line without a phenotype cannot be trained on, so it is put in the test set.
            </Note>
            <Toggle label="Also put the phenotypes in the test file (lets a model be checked against the test lines)" checked={keepTest} onChange={setKeepTest} />
          </Card>
          <div className="row">
            <Card title="Number of Clusters" icon={Group} delay={0.06}>
              <Slider label="K" value={k} onChange={setK} min={2} max={20} />
            </Card>
            <Card title="Test Size" icon={Percent} delay={0.12}>
              <Slider label="Fraction held out for testing" value={testSize} onChange={setTestSize}
                min={0.1} max={0.8} step={0.05} format={(v) => v.toFixed(2)} />
            </Card>
          </div>
        </TabPane>
      )}
      {tab === 'results' && (
        <TabPane id="results" key="results">
          <SectionHeader title="Results">
            {r ? `${r.n_samples} samples, ${r.n_features} markers, K=${r.n_clusters_used} (${r.split_method} split)` : 'Run the analysis to see results here.'}
          </SectionHeader>
          {!r ? (
            <EmptyState icon={Group}>Run the analysis to see results here.</EmptyState>
          ) : (
            <>
              <SummaryGrid items={[
                { label: 'Samples', value: r.n_samples },
                { label: 'Markers', value: r.n_features },
                { label: 'Clusters', value: r.n_clusters_used },
                { label: 'Train samples', value: r.n_train },
                { label: 'Test samples', value: r.n_test },
                { label: 'PC1 / PC2 variance', value: `${r.pc1_variance}% / ${r.pc2_variance}%` },
              ]} />
              {r.pheno_used && (
                <Note>
                  Phenotypes ({r.pheno_column}) were matched for <b>{r.n_pheno_matched}</b> lines and added to the training file
                  {r.test_has_pheno ? ' and the test file' : ''}.
                  {r.n_geno_without_pheno > 0 && <> {r.n_geno_without_pheno} genotyped line(s) have no phenotype (in the test set).</>}
                  {r.n_pheno_without_geno > 0 && <> {r.n_pheno_without_geno} phenotyped line(s) are not in the genotype file and were left out.</>}
                </Note>
              )}
              {(r.notes || []).filter((n) => /test file has no phenotypes/.test(n)).map((n) => (
                <div className="warning-box" key={n} style={{ marginTop: 12, marginBottom: 0 }}><TriangleAlert size={18} /><span>{n}</span></div>
              ))}
              <ImageGallery single images={[{ src: r.plot_png, label: 'PCA — K-Means clusters' }]} />
              <DownloadButtons onSave={tool.save} items={[
                { icon: Download, label: 'Plot (PNG)', path: r.plot_png, name: 'cluster_plot.png' },
                { icon: FileText, label: 'Training Data (CSV)', path: r.train_csv, name: 'Training_data_Cluster.csv' },
                { icon: FileText, label: 'Testing Data (CSV)', path: r.test_csv, name: 'Testing_data_Cluster.csv' },
                { icon: FileText, label: 'Cluster Membership (CSV)', path: r.cluster_table_csv, name: 'cluster_membership.csv' },
              ]} />
              <div className="note" style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', gap: 16, flexWrap: 'wrap' }}>
                <span><b>Next step:</b> build the genomic relationship matrices of these training and testing files. They are filled in for you.</span>
                <button className="btn btn-primary" onClick={() => navigate && navigate('additive_dominance')}>
                  Relationship matrices <ArrowRight size={16} />
                </button>
              </div>
              <DataTable rows={r.cluster_table} colorColumn="Cluster" />
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
