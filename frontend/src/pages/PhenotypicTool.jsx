import { useState } from 'react'
import {
  ChartLine, SlidersHorizontal, LayoutGrid, FileSpreadsheet, Dna, Globe, Group, Palette, List, Table2, ArrowDownUp,
  Activity, ChartColumn, Calculator, Waypoints, Braces, Grid3x3, Scale, Download, FileArchive,
} from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, SelectField, Toggle, ChipGroup, FilePicker, SubTabs,
  EmptyState, PlotSlot, TableSlot, DownloadButtons,
} from '../components/kit.jsx'
import { useBridge, useTool, useSignal, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'input', label: 'Input & Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]

const p = (name, title, height) => ({ k: 'plot', name, title, height })
const t = (name, title) => ({ k: 't', name, title })

// Each sub-tab: rows of blocks; a row with 2+ blocks becomes a two-column grid.
const SECTIONS = [
  { id: 'env1', label: 'Env: Summary', icon: List, rows: [[t('env_desc_stats', 'Descriptive statistics by environment')], [p('env_boxplot', 'Boxplot')], [p('env_violin_plot', 'Violin plot')]] },
  { id: 'env2', label: 'Env: ANOVA & Means', icon: Table2, rows: [[t('env_anova', 'ANOVA')], [p('env_means_plot', 'Means plot')], [t('env_means', 'Environment means')]] },
  { id: 'env3', label: 'Env: Ranks & Correlation', icon: ArrowDownUp, rows: [[p('env_rank_plot', 'Rank plot')], [t('env_rank', 'Ranks')], [p('env_cor_heatmap', 'Correlation heatmap')], [t('env_cor', 'Correlation table')]] },
  { id: 'env4', label: 'Env: BLUEs', icon: Activity, rows: [[p('env_blues_plot', 'BLUEs plot')], [t('env_blues', 'BLUEs')]] },
  { id: 'desc', label: 'Descriptives', icon: ChartColumn, rows: [[t('descriptive_stats', 'Descriptive statistics')], [p('distribution_plot', 'Distribution')]] },
  { id: 'anova', label: 'ANOVA', icon: Calculator, rows: [[t('anova_table', 'ANOVA table')], [p('residual_plot', 'Residuals'), p('qq_plot', 'Normal Q-Q')]] },
  { id: 'mm', label: 'Mixed Models', icon: ChartLine, rows: [[p('blues_plot', 'BLUEs')], [t('mixed_model_blues', 'Mixed-model BLUEs')]] },
  { id: 'herit', label: 'Heritability', icon: Dna, rows: [[t('variance_components', 'Variance components'), t('heritability', 'Heritability')], [p('heritability_plot', 'Heritability')]] },
  { id: 'corr', label: 'Correlations', icon: Waypoints, rows: [[p('correlation_heatmap', 'Correlation heatmap'), t('correlation_table', 'Correlation table')], [p('correlation_network', 'Correlation network')]] },
  { id: 'clust', label: 'Clustering', icon: Group, rows: [[p('cluster_plot', 'Clusters')], [p('dendrogram', 'Dendrogram'), t('cluster_summary', 'Cluster summary')]] },
  { id: 'pca', label: 'PCA', icon: Braces, rows: [[p('pca_plot', 'PCA'), p('pca_variance_plot', 'Variance explained')], [t('pca_loadings', 'Loadings'), t('pca_variance', 'Variance')]] },
  { id: 'stab', label: 'Stability', icon: Scale, rows: [[p('stability_plot', 'Stability'), p('ammi_plot', 'AMMI')], [t('stability_parameters', 'Stability parameters')]] },
  { id: 'dl', label: 'Downloads', icon: Download, rows: null },
]

const COLOR_SCHEMES = [
  { value: 'default', label: 'Default (Green/Red)' }, { value: 'blue', label: 'Blue Gradient' },
  { value: 'viridis', label: 'Viridis' }, { value: 'plasma', label: 'Plasma' }, { value: 'magma', label: 'Magma' },
]

export default function PhenotypicTool() {
  const { bridge } = useBridge('phenotypicBridge', () => mockGeneric(['dataLoaded']))
  const tool = useTool(bridge, 'run_analysis')
  const [tab, setTab] = useState('input')
  const [sub, setSub] = useState('env1')
  const [file, setFile] = useState('')
  const [preview, setPreview] = useState(null)
  const [trait, setTrait] = useState('')
  const [interaction, setInteraction] = useState(false)
  const [mode, setMode] = useState('all')
  const [envs, setEnvs] = useState([])
  const [k, setK] = useState(3)
  const [method, setMethod] = useState('kmeans')
  const [scheme, setScheme] = useState('default')
  tool.onResult(() => { setTab('results'); setSub('env1') })

  useSignal(bridge, 'dataLoaded', (json) => {
    try {
      const pv = JSON.parse(json)
      if (pv.status !== 'ok') return tool.setStatus('Error reading file: ' + pv.message)
      setPreview(pv)
      setTrait(pv.traits[0] || '')
      setEnvs(pv.environments)
      tool.setStatus(`Loaded ${pv.n_obs} observations, ${pv.traits.length} traits, ${pv.n_env} environments.`)
    } catch (e) { console.error('[PHENO] dataLoaded', e) }
  })

  const run = () => {
    if (!file) return tool.fail('Please select an input file.')
    if (!trait) return tool.fail('Please select a trait.')
    tool.run({
      input_file: file, trait, include_interaction: interaction, analysis_mode: mode,
      selected_environments: mode === 'specific' ? envs : [], n_clusters: k, clustering_method: method, color_scheme: scheme,
    })
  }
  const r = tool.result
  const sec = SECTIONS.find((s) => s.id === sub)

  return (
    <ToolPage
      id="phen" icon={ChartLine} tabs={TABS} tab={tab} setTab={setTab}
      title="Phenotypic Data Analysis"
      subtitle="Descriptive stats, ANOVA, mixed models & BLUEs, heritability, correlations, PCA, clustering and G×E stability"
      run={{ label: 'Run Full Analysis', runningLabel: 'Analysing…', state: tool.state, status: tool.status, onRun: run }}
    >
      {tab === 'input' && (
        <TabPane id="input" key="input">
          <SectionHeader title="Data Upload & Parameters">
            Requires columns: Genotype, Environment, Rep, Block, then trait columns from column 10 onward.
          </SectionHeader>
          <Card title="Phenotypic Data File" icon={FileSpreadsheet}>
            <FilePicker value={file} placeholder="Select .xlsx, .xls or .csv file"
              onBrowse={async () => { const x = await tool.browse(); if (x) setFile(x) }} />
            {preview && (
              <div className="summary-grid" style={{ marginTop: 16, marginBottom: 0 }}>
                {[['Variables', preview.n_vars], ['Observations', preview.n_obs], ['Genotypes', preview.n_geno ?? 'N/A'], ['Environments', preview.n_env]].map(([l, v]) => (
                  <div key={l} className="summary-card"><div className="val">{v}</div><div className="lbl">{l}</div></div>
                ))}
              </div>
            )}
          </Card>
          <div className="row">
            <Card title="Trait & Model" icon={Dna} delay={0.05}>
              <SelectField label="Trait to analyze" value={trait} onChange={setTrait}
                options={preview && preview.traits.length ? preview.traits : [{ value: '', label: 'Upload a file first' }]} />
              <div style={{ marginTop: 10 }}><Toggle label="Include G×E interaction" checked={interaction} onChange={setInteraction} /></div>
            </Card>
            <Card title="Environment Selection" icon={Globe} delay={0.1}>
              <SelectField label="Analysis mode" value={mode} onChange={setMode}
                options={[{ value: 'all', label: 'All Environments' }, { value: 'specific', label: 'Specific Environment(s)' }]} />
              {mode === 'specific' && preview && (
                <div className="field">
                  <label className="field-title">Select environment(s)</label>
                  <ChipGroup options={preview.environments} value={envs} onChange={setEnvs} />
                </div>
              )}
            </Card>
          </div>
          <div className="row">
            <Card title="Clustering" icon={Group} delay={0.15}>
              <Slider label="Number of clusters" value={k} onChange={setK} min={2} max={10} />
              <SelectField label="Method" value={method} onChange={setMethod} options={[{ value: 'kmeans', label: 'K-means' }, { value: 'hierarchical', label: 'Hierarchical' }]} />
            </Card>
            <Card title="Color Scheme" icon={Palette} delay={0.2}>
              <SelectField value={scheme} onChange={setScheme} options={COLOR_SCHEMES} />
            </Card>
          </div>
        </TabPane>
      )}
      {tab === 'results' && (
        <TabPane id="results" key="results">
          <SectionHeader title="Results">
            {r ? `Trait: ${r.trait}` + (Array.isArray(r.active_envs) && r.active_envs.length ? ` — Environments: ${r.active_envs.join(', ')}` : ' — All Environments') : 'Run the analysis to see results here.'}
          </SectionHeader>
          {!r ? <EmptyState icon={ChartLine}>Run the analysis to see results here.</EmptyState> : (
            <>
              <SubTabs tabs={SECTIONS} value={sub} onChange={setSub} id="phsub" />
              {sec.rows ? sec.rows.map((row, i) => (
                <div key={sub + i} className={row.length > 1 ? 'plot-grid-2' : undefined}>
                  {row.map((b) => b.k === 'plot'
                    ? <PlotSlot key={b.name} name={b.name} title={b.title} plot={(r.plots || {})[b.name]} onSave={tool.save} />
                    : <TableSlot key={b.name} title={b.title} rows={(r.tables_data || {})[b.name]} />)}
                </div>
              )) : (
                <Card title="Excel Workbooks" icon={FileArchive}>
                  <DownloadButtons onSave={tool.save} items={[
                    r.excel_all && { icon: FileSpreadsheet, label: 'Complete Analysis (xlsx)', path: r.excel_all, name: 'Phenotypic_Analysis.xlsx' },
                    r.excel_env && { icon: FileSpreadsheet, label: 'Environment Analysis (xlsx)', path: r.excel_env, name: 'Environment_Analysis.xlsx' },
                  ].filter(Boolean)} />
                </Card>
              )}
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
