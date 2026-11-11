import { useEffect, useState } from 'react'
import {
  Layers, Database, SlidersHorizontal, LayoutGrid, FileText, Ruler, Scale, FlaskConical,
  FileSpreadsheet, FileArchive, ChartArea,
} from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, Toggle, FilePicker, EmptyState, ImageGallery,
  DownloadButtons, DataTable, SelectField, TextField, Note,
} from '../components/kit.jsx'
import { useBridge, useTool, mockGeneric, callJson } from '../lib/bridge.js'

const SOURCES = [
  { value: '3k', label: 'Rice 3K panel — phenotype file only' },
  { value: 'vcf', label: 'My own genotypes — VCF + phenotype file' },
  { value: 'file', label: 'My own genotypes — ready-made file' },
]

const CALL_BY = [
  { value: 'fdr', label: 'FDR (q-value) below the level' },
  { value: 'suggestive', label: 'Suggestive: p-value below the threshold' },
]
const SUGGESTIVE = [
  { value: 1e-3, label: '1 × 10⁻³' }, { value: 1e-4, label: '1 × 10⁻⁴' }, { value: 5e-5, label: '5 × 10⁻⁵ (paper)' },
  { value: 1e-5, label: '1 × 10⁻⁵' }, { value: 1e-6, label: '1 × 10⁻⁶' },
]

const TABS = [
  { id: 'data', label: 'Data & Bulks', icon: Database },
  { id: 'advanced', label: 'Statistics', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]

const MTA_COLS = ['SNP', 'Chr', 'Position', 'MAF', 'AF_low', 'AF_random', 'AF_high', 'Delta_AF_high_minus_low', 'P_value', 'Q_value', 'PVE_%', 'Counted_allele_enriched_in']
const REGION_COLS = ['Region', 'Chr', 'Start', 'End', 'Size_bp', 'N_SNPs', 'Peak_SNP', 'Peak_P', 'Enriched_in', 'Max_PVE_%']
const BULK_COLS = ['Bulk', 'Lines', 'Mean', 'Median', 'Min', 'Max']

const DEFAULTS = {
  bulk_pct: 20, maf: 0.05, min_alleles: 0, seed: 1,
  use_gc: true, call_by: 'fdr', fdr: 0.05, suggestive: 5e-5, merge_kb: 100,
}

export default function EtGwasTool() {
  const { bridge } = useBridge('etGwasBridge', () => mockGeneric())
  const tool = useTool(bridge, 'run_et_gwas')
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
    const cfg = { ...c, min_alleles: Number(c.min_alleles) || 0, bulk_pct: Number(c.bulk_pct), seed: Number(c.seed) || 1 }
    if (mode === 'file') {
      if (!file) return tool.fail('Please select the input file.')
      return tool.run({ input_mode: 'file', input_file: file, ...cfg })
    }
    if (!pheno || !trait) return tool.fail('Please select the phenotype file and the trait.')
    if (mode === 'vcf' && !vcf) return tool.fail('Please select the VCF file.')
    tool.run({ input_mode: mode, phenotype_file: pheno, trait, vcf_file: mode === 'vcf' ? vcf : '', ...cfg })
  }

  const r = tool.result
  const plots = r && r.plots
    ? Object.entries(r.plots).map(([name, src]) => ({ src, label: name.replace(/^\d+_/, '').replace(/_/g, ' ') }))
    : []

  return (
    <ToolPage
      id="et_gwas" icon={Layers} tabs={TABS} tab={tab} setTab={setTab}
      title="ET-GWAS — Extreme Trait GWAS"
      subtitle="Compares the allele frequencies of low, random and high phenotype bulks to find associations, including rare variants"
      run={{ label: 'Run ET-GWAS', runningLabel: 'Running ET-GWAS…', state: tool.state, status: tool.status, onRun: run }}
    >
      {tab === 'data' && (
        <TabPane id="data" key="data">
          <SectionHeader title="Input Data & Bulks">
            Extreme Trait GWAS (Gnanapragasam et al. 2024, Life Science Alliance, doi 10.26508/lsa.202302352): the lines with the lowest and the
            highest phenotypes form two bulks, an equal number of lines from the middle form a random bulk, and every SNP is tested for a
            steady change of its allele frequency from the low to the random to the high bulk.
          </SectionHeader>
          <Card title="Genotype source" icon={Database}>
            <SelectField label="Use" value={mode} onChange={setMode} options={SOURCES} />
            {mode === '3k' && (
              <Note>
                Phenotype file only — the rice 3K genotypes ({info3k ? `${info3k.n_samples.toLocaleString()} accessions, ${info3k.n_markers.toLocaleString()} SNPs` : '3,024 accessions, 52,600 SNPs'}) are built in.
                The first column must hold the accession IDs{info3k ? <> (e.g. <code>{info3k.examples.slice(0, 3).join(', ')}</code>)</> : null} or the accession names; the other columns are traits.
                The method needs a panel with a wide phenotype range: at least about 100 phenotyped accessions, ideally several hundred.
              </Note>
            )}
            {mode === 'vcf' && (
              <Note>A VCF (.vcf or .vcf.gz, GT field) and a phenotype file whose first column holds the sample IDs exactly as in the VCF header. Only biallelic SNPs on numbered chromosomes are used.</Note>
            )}
            {mode === 'file' && (
              <Note>A ready-made file: first column = ID, second = phenotype, then one column per SNP named <code>Chr1_123456</code> (chromosome and position) with 0/1/2 dosages.</Note>
            )}
          </Card>
          <Card title="Input Data" icon={FileText} delay={0.03}>
            {mode === 'file' ? (
              <FilePicker value={file} onBrowse={browse} placeholder="Select genotype + phenotype CSV, TSV or TXT file" />
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
          <Card title="Bulks and filters" icon={Ruler} delay={0.06}>
            <div className="row" style={{ gap: 36 }}>
              <div style={{ flex: 1, minWidth: 260 }}>
                <Slider label="Bulk size (% of the lines in each extreme bulk)" value={c.bulk_pct} onChange={set('bulk_pct')} min={5} max={40} step={1} format={(v) => `${v}%`} />
                <Note>The paper compared 10, 15 and 20% and found 20% best: a diversity panel is highly segregating, so small bulks miss many SNPs.</Note>
              </div>
              <div style={{ flex: 1, minWidth: 260 }}>
                <Slider label="Minimum MAF" value={c.maf} onChange={set('maf')} min={0} max={0.5} step={0.01} format={(v) => v.toFixed(2)} />
                <TextField label="Minimum called alleles per bulk (0 = the bulk size)" type="number" min={0} value={c.min_alleles} onChange={set('min_alleles')}
                  hint="SNPs seen in too few lines of a bulk are dropped. Missingness is otherwise not filtered, so that rare variants stay." />
              </div>
            </div>
          </Card>
        </TabPane>
      )}
      {tab === 'advanced' && (
        <TabPane id="advanced" key="advanced">
          <SectionHeader title="Statistics">How associations are called and grouped into regions.</SectionHeader>
          <div className="row">
            <Card title="Significance" icon={Scale}>
              <SelectField label="An SNP is an association (MTA) when" value={c.call_by} onChange={set('call_by')} options={CALL_BY} />
              <Slider label="FDR level" value={c.fdr} onChange={set('fdr')} min={0.01} max={0.2} step={0.01} format={(v) => v.toFixed(2)} />
              <SelectField label="Suggestive threshold (p-value)" value={c.suggestive} onChange={(v) => set('suggestive')(Number(v))} options={SUGGESTIVE} />
              <Toggle label="Divide the statistics by the inflation factor (genomic control)" checked={c.use_gc} onChange={set('use_gc')} />
            </Card>
            <Card title="Regions and random bulk" icon={ChartArea} delay={0.06}>
              <Slider label="Merge MTAs closer than (kb)" value={c.merge_kb} onChange={set('merge_kb')} min={10} max={1000} step={10} />
              <Note>The paper treats MTAs in the same 100 kb bin as the same locus.</Note>
              <TextField label="Random seed (random bulk)" type="number" min={1} value={c.seed} onChange={set('seed')}
                hint="Change it to see how much the random bulk matters; the extreme bulks do not change." />
            </Card>
          </div>
        </TabPane>
      )}
      {tab === 'results' && (
        <TabPane id="results" key="results">
          <SectionHeader title="Results">Associations, regions, the bulks that were used, plots and downloads.</SectionHeader>
          {!r ? (
            <EmptyState icon={FlaskConical}>Run the analysis to see results here.</EmptyState>
          ) : (
            <>
              <DownloadButtons onSave={tool.save} items={[
                { icon: FileSpreadsheet, label: 'Excel Report', path: r.excel, name: 'et_gwas_results.xlsx' },
                { icon: FileText, label: 'MTAs (CSV)', path: r.mta_csv, name: 'et_gwas_mtas.csv' },
                { icon: FileText, label: 'All SNPs (CSV)', path: r.all_csv, name: 'et_gwas_all_snps.csv' },
                { icon: FileText, label: 'Bulk members (CSV)', path: r.bulks_csv, name: 'et_gwas_bulks.csv' },
                { icon: FileArchive, label: 'All Plots (ZIP)', path: r.zip, name: 'et_gwas_plots.zip' },
              ]} />
              {r.summary && <pre className="summary-box">{r.summary}</pre>}
              <div className="plot-title">Marker–trait associations (MTAs)</div>
              {(r.mta_table || []).length ? (
                <DataTable rows={r.mta_table} columns={MTA_COLS} maxHeight={380} />
              ) : (
                <EmptyState icon={FlaskConical}>No SNP passed the significance rule. Try the suggestive threshold, a bulk size of 20%, or check the QQ plot and the inflation factor.</EmptyState>
              )}
              {(r.region_table || []).length > 0 && (
                <>
                  <div className="plot-title">Regions</div>
                  <DataTable rows={r.region_table} columns={REGION_COLS} maxHeight={260} />
                </>
              )}
              <div className="plot-title">Bulks</div>
              <DataTable rows={r.bulk_table || []} columns={BULK_COLS} maxHeight={220} />
              <div className="plot-title">Plots</div>
              <ImageGallery images={plots} />
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
