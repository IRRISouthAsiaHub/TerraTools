import { useState } from 'react'
import { AnimatePresence } from 'framer-motion'
import { Group, SlidersHorizontal, LayoutGrid, Dna, Percent, Download, FileText } from 'lucide-react'
import {
  Background, ToolHeader, Tabs, TabPane, SectionHeader, Card, Slider, FilePicker,
  RunBar, EmptyState, SummaryGrid, ImageGallery, DownloadButtons, DataTable,
} from '../components/kit.jsx'
import { useBridge, useTool, mockToolBridge } from '../lib/bridge.js'

const TABS = [
  { id: 'input', label: 'Input & Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]

export default function ClusteringTool() {
  const { bridge } = useBridge('clusteringBridge', () => mockToolBridge('clustering'))
  const tool = useTool(bridge, 'run_clustering')
  const [tab, setTab] = useState('input')
  const [file, setFile] = useState('')
  const [k, setK] = useState(5)
  const [testSize, setTestSize] = useState(0.2)

  tool.onResult(() => setTab('results'))

  const browse = async () => { const p = await tool.browse(); if (p) setFile(p) }
  const run = () => {
    if (!file) return tool.fail('Please select an input file.')
    tool.run({ input_file: file, n_clusters: k, test_size: testSize })
  }

  const r = tool.result

  return (
    <>
      <Background variant="mint" count={10} />
      <div className="page">
        <div className="container narrow">
          <ToolHeader
            icon={Group}
            title="Genotypic Clustering & Genomic Selection"
            subtitle="PCA + K-Means clustering with stratified train/test split for genomic selection"
          />
          <div className="panel on-mint">
            <Tabs tabs={TABS} value={tab} onChange={setTab} idPrefix="clu" />
            <div className="tab-content">
              <AnimatePresence mode="wait">
                {tab === 'input' && (
                  <TabPane id="input" key="input">
                    <SectionHeader title="Genotype Data">
                      Upload SNP genotype data (CSV or TSV). Missing values coded as -1 are handled automatically.
                    </SectionHeader>
                    <Card title="Genotype File" icon={Dna}>
                      <FilePicker value={file} onBrowse={browse} placeholder="Select CSV/TSV (sample IDs in first column, SNP columns after)" />
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
                        <ImageGallery single images={[{ src: r.plot_png, label: 'PCA — K-Means clusters' }]} />
                        <DownloadButtons onSave={tool.save} items={[
                          { icon: Download, label: 'Plot (PNG)', path: r.plot_png, name: 'cluster_plot.png' },
                          { icon: FileText, label: 'Training Data (CSV)', path: r.train_csv, name: 'Training_data_Cluster.csv' },
                          { icon: FileText, label: 'Testing Data (CSV)', path: r.test_csv, name: 'Testing_data_Cluster.csv' },
                          { icon: FileText, label: 'Cluster Membership (CSV)', path: r.cluster_table_csv, name: 'cluster_membership.csv' },
                        ]} />
                        <DataTable rows={r.cluster_table} colorColumn="Cluster" />
                      </>
                    )}
                  </TabPane>
                )}
              </AnimatePresence>
            </div>
            <RunBar label="Run Clustering & Split" runningLabel="Clustering…" state={tool.state} status={tool.status} onRun={run} />
          </div>
        </div>
      </div>
    </>
  )
}
