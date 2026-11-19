import { useEffect, useState } from 'react'
import {
  Globe, Database, SlidersHorizontal, ChartArea, LayoutGrid, FileText, Dna,
  Scale, Waypoints, Link2, FlaskConical, FileSpreadsheet, FileArchive, ShieldCheck, Network,
} from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, Toggle, ChipGroup,
  FilePicker, EmptyState, ImageGallery, DownloadButtons, SelectField, TextField, Note,
} from '../components/kit.jsx'
import { useBridge, useTool, mockToolBridge, callJson } from '../lib/bridge.js'

const SOURCES = [
  { value: '3k', label: 'Rice 3K panel — phenotype file only' },
  { value: 'vcf', label: 'My own genotypes — VCF + phenotype file' },
  { value: 'file', label: 'My own genotypes — ready-made GWAS file' },
]

const ENCODINGS = [
  { value: 'auto', label: 'Detect (stops and asks if the values do not decide it)' },
  { value: 'dosage', label: '0 / 1 / 2 — copies of the counted allele' },
  { value: 'centred', label: '−1 / 0 / 1 — centred genotypes' },
]
const KINSHIPS = [
  { value: 'vanraden', label: 'VanRaden, all QC-passed markers' },
  { value: 'ld_pruned', label: 'VanRaden, LD-pruned markers' },
]


const MODELS = ['OLS', 'MLM', 'EMMAX', 'FarmCPU', 'GEMMA', 'BLINK', 'mrMLM', 'FASTmrMLM']

const TABS = [
  { id: 'data', label: 'Data & Models', icon: Database },
  { id: 'advanced', label: 'Advanced', icon: SlidersHorizontal },
  { id: 'plots', label: 'Plot Settings', icon: ChartArea },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]

const DEFAULTS = {
  models: ['OLS', 'MLM'],
  genotype_encoding: 'auto', missing_codes: '',
  maf: 0, max_snp_missing: 1, max_sample_missing: 1, drop_duplicate_positions: false, imputation: 'mean',
  kinship_method: 'vanraden', ld_prune_r2: 0.2, loco: false, n_pcs: 3, pca_ld_prune: false, compare_pcs: '',
  significance_rule: 'bonferroni', alpha: 0.05, permutations: 1000, permutation_null: 'residual_permutation', permutation_seed: 1,
  top_n: 15, minimum_pve: 0,
  ols_pcs: 0, mlm_pcs: 0, emmax_pcs: 0, gemma_pcs: 0, gemma_test: 'wald', compare_gemma: false, gemma_binary: '',
  farmcpu_iterations: 10,
  mrmlm_radius_kb: 20, mrmlm_lod: 3, mrmlm_n_pca: 3,
  show_bonferroni: true, show_suggestive: true, show_pve: false,
}

export default function GwasTool() {
  const { bridge } = useBridge('gwasBridge', () => mockToolBridge('gwas'))
  const tool = useTool(bridge, 'run_gwas')
  const [tab, setTab] = useState('data')
  const [mode, setMode] = useState('3k')
  const [file, setFile] = useState('')
  const [vcf, setVcf] = useState('')
  const [pheno, setPheno] = useState('')
  const [traits, setTraits] = useState(null)
  const [trait, setTrait] = useState('')
  const [info3k, setInfo3k] = useState(null)
  const [c, setC] = useState(DEFAULTS)
  const set = (k) => (v) => setC((p) => ({ ...p, [k]: v }))

  tool.onResult(() => setTab('results'))
  useEffect(() => {
    let live = true
    if (bridge && bridge.rice3k_info) callJson(bridge, 'rice3k_info').then((i) => { if (live && i && i.n_samples) setInfo3k(i) })
    return () => { live = false }
  }, [bridge])

  const browse = async () => { const p = await tool.browse(); if (p) setFile(p) }
  const browseVcf = async () => { const p = await tool.browse('browse_vcf'); if (p) setVcf(p) }
  const browsePheno = async () => {
    const p = await tool.browse('browse_phenotype')
    if (!p) return
    setPheno(p)
    const cols = await callJson(bridge, 'get_phenotype_columns', p)
    setTraits(cols || [])
    setTrait((cols && cols[0]) || '')
    tool.setStatus(cols && cols.length ? 'Phenotype file loaded. Choose the trait.' : 'No numeric trait columns found in the phenotype file.')
  }
  const run = () => {
    if (!c.models.length) return tool.fail('Please select at least one model.')
    if (mode === 'file') {
      if (!file) return tool.fail('Please select an input file.')
      return tool.run({ input_mode: 'file', input_file: file, ...c })
    }
    if (!pheno || !trait) return tool.fail('Please select the phenotype file and the trait.')
    if (mode === 'vcf' && !vcf) return tool.fail('Please select the VCF file.')
    tool.run({ input_mode: mode, phenotype_file: pheno, trait, vcf_file: mode === 'vcf' ? vcf : '', ...c })
  }

  const r = tool.result
  const plots = r
    ? Object.entries(r.plots).map(([name, src]) => ({ src, label: name.replace(/^\d+_/, '').replace(/_/g, ' ') }))
    : []
  const pct = (v) => `${Math.round(v * 100)}%`

  return (
    <ToolPage
      id="gwas" icon={Globe} tabs={TABS} tab={tab} setTab={setTab}
      title="GWAS Analysis"
      subtitle="Genome-wide association study — genotype–phenotype analysis across multiple models"
      run={{ label: 'Run GWAS', runningLabel: 'Running GWAS…', state: tool.state, status: tool.status, onRun: run }}
    >
      {tab === 'data' && (
        <TabPane id="data" key="data">
          <SectionHeader title="Input Data & Models">
            Choose the genotypes, then the GWAS models to run.
          </SectionHeader>
          <Card title="Genotype source" icon={Database}>
            <SelectField label="Use" value={mode} onChange={setMode} options={SOURCES} />
          </Card>
          <Card title="Input Data" icon={FileText} delay={0.03}>
            {mode === 'file' ? (
              <>
                <FilePicker value={file} onBrowse={browse} placeholder="Select genotype/phenotype CSV, TSV or TXT file" />
                <div className="field-hint">Columns: ID, phenotype, then one column per marker. Marker positions come from the marker names (Chr3_123456) or from a position row and a chromosome row under the header.</div>
              </>
            ) : (
              <>
                {mode === 'vcf' && (
                  <div className="field">
                    <label className="field-title">VCF file</label>
                    <FilePicker value={vcf} onBrowse={browseVcf} placeholder="Select variants.vcf or variants.vcf.gz" />
                  </div>
                )}
                <div className="field">
                  <label className="field-title">Phenotype file (CSV, TSV or Excel)</label>
                  <FilePicker value={pheno} onBrowse={browsePheno} placeholder="Select phenotype file — ID column first, then traits" />
                </div>
                <SelectField label="Trait to analyse" value={trait} onChange={setTrait}
                  options={traits && traits.length ? traits : [{ value: '', label: traits ? 'No trait columns found' : 'Load the phenotype file first' }]} />
                {mode === '3k' && info3k && info3k.template && (
                  <DownloadButtons onSave={tool.save} items={[{ icon: FileText, label: 'Download ID template (CSV)', path: info3k.template, name: 'rice3k_phenotype_template.csv' }]} />
                )}
              </>
            )}
          </Card>
          {mode === 'file' && (
            <Card title="Genotype coding" icon={ShieldCheck} delay={0.04}>
              <SelectField label="Coding of the genotype matrix" value={c.genotype_encoding} onChange={set('genotype_encoding')} options={ENCODINGS} />
              <TextField label="Extra values to read as missing (optional)" value={c.missing_codes} onChange={set('missing_codes')} placeholder="e.g. 3, -5" />
            </Card>
          )}
          <Card title="GWAS Models" icon={Dna} delay={0.05}>
            <ChipGroup options={MODELS} value={c.models} onChange={set('models')} />
          </Card>
        </TabPane>
      )}
      {tab === 'advanced' && (
        <TabPane id="advanced" key="advanced">
          <SectionHeader title="Advanced Model Parameters">Relationship matrix, population structure and per-model settings.</SectionHeader>
          <div className="row">
            <Card title="Relationship matrix" icon={Network}>
              <SelectField label="Kinship" value={c.kinship_method} onChange={set('kinship_method')} options={KINSHIPS} />
              {c.kinship_method === 'ld_pruned' && (
                <Slider label="LD pruning r² threshold" value={c.ld_prune_r2} onChange={set('ld_prune_r2')} min={0.05} max={1} step={0.05} format={(v) => v.toFixed(2)} />
              )}
              <Toggle label="Leave-one-chromosome-out (the tested marker's chromosome is left out of K)" checked={c.loco} onChange={set('loco')} />
            </Card>
          </div>
          <div className="row">
            <Card title="Covariates per model" icon={Scale} delay={0.06}>
              <Slider label="OLS PCs (0 = uncorrected)" value={c.ols_pcs} onChange={set('ols_pcs')} min={0} max={30} />
              <Slider label="MLM PCs (0 = kinship only)" value={c.mlm_pcs} onChange={set('mlm_pcs')} min={0} max={30} />
              <Slider label="EMMAX PCs" value={c.emmax_pcs} onChange={set('emmax_pcs')} min={0} max={30} />
              <Slider label="GEMMA-style PCs" value={c.gemma_pcs} onChange={set('gemma_pcs')} min={0} max={30} />
              <SelectField label="GEMMA-style test" value={c.gemma_test} onChange={set('gemma_test')} options={[{ value: 'wald', label: 'Wald (REML fits)' }, { value: 'lrt', label: 'Likelihood ratio (ML fits)' }]} />
              <Toggle label="Compare the GEMMA-style scan with the GEMMA program (if installed)" checked={c.compare_gemma} onChange={set('compare_gemma')} />
              {c.compare_gemma && <TextField label="GEMMA executable (leave empty to search the PATH)" value={c.gemma_binary} onChange={set('gemma_binary')} placeholder="C:\tools\gemma.exe" />}
            </Card>
            <Card title="FarmCPU" icon={Waypoints} delay={0.08}>
              <Slider label="Iterations" value={c.farmcpu_iterations} onChange={set('farmcpu_iterations')} min={1} max={100} />
            </Card>
          </div>
          <div className="row">
            <Card title="mrMLM / FASTmrMLM" icon={Dna} delay={0.14}>
              <Slider label="Search radius (kb)" value={c.mrmlm_radius_kb} onChange={set('mrmlm_radius_kb')} min={1} max={500} />
              <Slider label="Critical LOD (QTN rule)" value={c.mrmlm_lod} onChange={set('mrmlm_lod')} min={1} max={10} step={0.5} format={(v) => v.toFixed(1)} />
              <Slider label="Structure PCs" value={c.mrmlm_n_pca} onChange={set('mrmlm_n_pca')} min={0} max={20} />
            </Card>
          </div>
        </TabPane>
      )}
      {tab === 'plots' && (
        <TabPane id="plots" key="plots">
          <SectionHeader title="Plot and table options">Reference lines on the Manhattan plot and an optional PVE filter.</SectionHeader>
          <Card>
            <Toggle label="Show Bonferroni threshold" checked={c.show_bonferroni} onChange={set('show_bonferroni')} />
            <Toggle label="Show suggestive threshold (1 / number of SNPs)" checked={c.show_suggestive} onChange={set('show_suggestive')} />
            <Toggle label="List only genome-wide / QTN markers with at least this PVE" checked={c.show_pve} onChange={set('show_pve')} />
            {c.show_pve && <Slider label="Minimum PVE (%)" value={c.minimum_pve} onChange={set('minimum_pve')} min={0} max={100} step={0.5} />}
          </Card>
        </TabPane>
      )}
      {tab === 'results' && (
        <TabPane id="results" key="results">
          <SectionHeader title="Results">Generated plots, downloadable Excel report, and plot archive.</SectionHeader>
          {!r ? (
            <EmptyState icon={FlaskConical}>Run the pipeline to see results here.</EmptyState>
          ) : (
            <>
              <ImageGallery images={plots} />
              <DownloadButtons onSave={tool.save} items={[
                { icon: FileSpreadsheet, label: 'Excel Report', path: r.excel, name: 'gwas_results.xlsx' },
                { icon: FileArchive, label: 'All Plots (ZIP)', path: r.zip, name: 'gwas_plots.zip' },
                r.gwas_input && { icon: FileText, label: 'GWAS-ready input (CSV)', path: r.gwas_input, name: 'gwas_ready_input.csv' },
              ].filter(Boolean)} />
              {r.message && <Note><strong>{r.message}.</strong></Note>}
              {Object.keys(r.errors || {}).length > 0 && (
                <Note>{Object.entries(r.errors).map(([k, v]) => <div key={k}><strong>{k}</strong> did not run: {v}</div>)}</Note>
              )}
              {r.summary && <pre className="summary-box">{r.summary}</pre>}
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
