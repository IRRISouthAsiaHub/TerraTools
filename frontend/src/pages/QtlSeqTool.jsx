import { useEffect, useMemo, useState } from 'react'
import {
  Dna, FileInput, SlidersHorizontal, LayoutGrid, ChartLine, FlaskConical, ShieldCheck, Download, FolderOpen, Ban, TriangleAlert,
  CircleCheck, CircleX, Info, Users, Settings2, FileText, Wrench, Activity,
} from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, FilePicker, TextField, SelectField, RadioChips, Toggle, EmptyState, SummaryGrid,
  ProgressBar, LogBox, SubTabs,
} from '../components/kit.jsx'
import { TrackPlot, LevelLegend, NullPlot, DepthLimitPlot } from '../components/snpCharts.jsx'
import { useBridge, useTool, useSignal, callJson, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'input', label: 'Input', icon: FileInput },
  { id: 'settings', label: 'Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
  { id: 'chromosome', label: 'Chromosomes', icon: ChartLine },
  { id: 'simulation', label: 'Simulation', icon: FlaskConical },
  { id: 'validation', label: 'Validation', icon: ShieldCheck },
]

const DEFAULTS = {
  mode: 'bam', high_name: 'High bulk', low_name: 'Low bulk', parent1_name: 'Parent 1', parent2_name: 'Parent 2', index_parent: 'parent2',
  population: 'F2', n_high: 20, n_low: 20, bc_recurrent: 'non_index',
  min_depth: 7, max_depth: 0, min_snp_index: 0.3, apply_snp_index_filter: true, min_parent_depth: 7, min_parent_gq: 30, min_qual: 0,
  require_pass: false, min_mapq: 20, min_baseq: 20, min_alt_reads: 2, min_alt_fraction: 0.05, base_error: 0.01,
  window_size: 2000000, window_step: 10000, min_snps_per_window: 10,
  replicates: 10000, seed: 12345, depth_mode: 'auto', fixed_depth: 20, ci_method: 'snp_average',
  qtl_min_level: 90, merge_gap: 100000, interval_basis: 'top_level', engine: 'auto', chromosomes: '',
  high_sample: '', low_sample: '', parent1_sample: '', parent2_sample: '',
}
const NUMERIC = ['n_high', 'n_low', 'min_depth', 'max_depth', 'min_snp_index', 'min_parent_depth', 'min_parent_gq', 'min_qual', 'min_mapq', 'min_baseq',
  'min_alt_reads', 'min_alt_fraction', 'base_error', 'window_size', 'window_step', 'min_snps_per_window', 'replicates', 'seed', 'fixed_depth',
  'qtl_min_level', 'merge_gap']

const EMPTY_FILES = {
  ref_fasta: '', high_bam: '', low_bam: '', parent1_bam: '', parent2_bam: '', high_vcf: '', low_vcf: '', parent_vcf: '', annotation: '', output_dir: '',
}

const fmtBp = (v) => (v === null || v === undefined ? '' : Number(v).toLocaleString())
const f3 = (v) => (typeof v === 'number' ? v.toFixed(3) : v)

function mockQtlSeq() {
  const m = mockGeneric(['progressUpdate', 'inputsChecked', 'validationReady'])
  let last = null
  const json = (u) => fetch(u).then((r) => r.json())
  return new Proxy(m, {
    get(t, name) {
      if (name === 'get_defaults') return (cb) => setTimeout(() => cb(JSON.stringify({ defaults: DEFAULTS, output_default: 'C:\\TerraTools\\results\\qtl_seq' })), 50)
      if (name === 'get_dependencies') return (cb) => setTimeout(() => cb(JSON.stringify([
        { name: 'samtools', found: true, where: 'inside WSL', message: 'samtools found (inside WSL).', needed_for: 'reading BAM files quickly (samtools mpileup)' },
        { name: 'bcftools', found: true, where: 'inside WSL', version: 'bcftools 1.22', message: 'bcftools 1.22 found (inside WSL).', needed_for: 'indexing and checking the filtered VCF' },
        { name: 'pysam', found: true, where: 'inside WSL', version: 'pysam 0.24.1', message: 'pysam 0.24.1 found (inside WSL).', needed_for: 'counting reads of BAM files' },
      ])), 50)
      if (name === 'get_validation_experiments') return (cb) => json('/qtlseq_demo_validation.json').then((v) => cb(JSON.stringify(v.refs))).catch(() => cb('[]'))
      if (name === 'check_inputs') return () => setTimeout(() => t.inputsChecked.emit(JSON.stringify({ status: 'success', ok: true, errors: [], warnings: ['No parental data: parental polymorphism cannot be validated and the analysis is a bulk-only allele-frequency comparison.'], files: {} })), 400)
      if (name === 'run_qtlseq') {
        return () => {
          let p = 0
          const id = setInterval(() => { p += 12; t.progressUpdate.emit(JSON.stringify({ pct: Math.min(p, 95), message: p < 40 ? 'Reading the high-bulk VCF' : p < 80 ? 'Simulating the null distribution of the Delta SNP-index' : 'Sliding-window averages' })) }, 250)
          json('/qtlseq_demo_result.json').then((r) => { last = r; setTimeout(() => { clearInterval(id); t.resultReady.emit(JSON.stringify(r)) }, 2200) })
            .catch(() => { clearInterval(id); t.resultReady.emit(JSON.stringify({ status: 'error', message: 'Preview mode: there is no Python backend.' })) })
        }
      }
      if (name === 'get_chromosome') return (c, cb) => setTimeout(() => { const x = last && last.plot.chromosomes.find((k) => k.name === c); cb(JSON.stringify(x ? { status: 'success', ...x, qtls: last.qtls.filter((q) => q.Chromosome === c) } : { status: 'error', message: 'Run the analysis first.' })) }, 80)
      if (name === 'run_validation') return (k) => json('/qtlseq_demo_validation.json').then((v) => setTimeout(() => t.validationReady.emit(JSON.stringify(v.result)), 1200))
      return t[name]
    },
  })
}

function Badge({ ok, children }) {
  return <span className={`badge ${ok ? '' : 'warn'}`}>{children}</span>
}

function Messages({ errors = [], warnings = [] }) {
  return (
    <>
      {errors.map((e) => <div className="info-line warn" key={e} style={{ display: 'flex', gap: 8 }}><CircleX size={16} style={{ flex: 'none', marginTop: 2 }} /><span>{e}</span></div>)}
      {warnings.map((w) => <div className="info-line" key={w} style={{ display: 'flex', gap: 8, background: 'rgba(233,196,106,.2)', color: '#7a5b0b' }}><TriangleAlert size={16} style={{ flex: 'none', marginTop: 2 }} /><span>{w}</span></div>)}
    </>
  )
}

export default function QtlSeqTool() {
  const { bridge } = useBridge('qtlSeqBridge', () => mockQtlSeq())
  const tool = useTool(bridge, 'run_qtlseq', { accept: (r) => r.status === 'cancelled' })
  const [tab, setTab] = useState('input')
  const [files, setFiles] = useState(EMPTY_FILES)
  const [s, setS] = useState(DEFAULTS)
  const [advanced, setAdvanced] = useState(false)
  const [deps, setDeps] = useState(null)
  const [check, setCheck] = useState(null)
  const [outDefault, setOutDefault] = useState('')
  const [prog, setProg] = useState({ pct: 0, message: '' })
  const [log, setLog] = useState([])
  const [view, setView] = useState('delta')
  const [chrom, setChrom] = useState('')
  const [chromData, setChromData] = useState(null)
  const [chromErr, setChromErr] = useState('')
  const [pickedQtl, setPickedQtl] = useState('')
  const [valKey, setValKey] = useState('blast_ril')
  const [valSeed, setValSeed] = useState(1)
  const [valRefs, setValRefs] = useState([])
  const [valRes, setValRes] = useState(null)
  const [valRunning, setValRunning] = useState(false)

  useEffect(() => {
    if (!bridge) return
    callJson(bridge, 'get_defaults').then((d) => { if (d && d.defaults) { setS((x) => ({ ...x, ...d.defaults })); setOutDefault(d.output_default || '') } })
    callJson(bridge, 'get_dependencies').then((d) => Array.isArray(d) && setDeps(d))
    callJson(bridge, 'get_validation_experiments').then((d) => Array.isArray(d) && setValRefs(d))
  }, [bridge])
  useSignal(bridge, 'progressUpdate', (j) => {
    try { const p = JSON.parse(j); setProg(p); setLog((l) => (l.length && l[l.length - 1] === p.message ? l : [...l.slice(-300), p.message])) } catch { /* ignore */ }
  })
  useSignal(bridge, 'inputsChecked', (j) => { try { setCheck(JSON.parse(j)) } catch { /* ignore */ } })
  useSignal(bridge, 'validationReady', (j) => { try { setValRes(JSON.parse(j)) } catch { /* ignore */ } setValRunning(false) })

  tool.onResult((r) => {
    if (r.status === 'cancelled') { tool.setState('idle'); tool.setStatus('Cancelled.'); return }
    if (r.status === 'success') {
      setView('delta'); setChrom(''); setChromData(null); setPickedQtl('')
      setTab('results')
    }
  })

  const r = tool.result && tool.result.status === 'success' ? tool.result : null
  const running = tool.state === 'running'
  const set = (k) => (v) => setS((x) => ({ ...x, [k]: v }))
  const setFile = (k, slot) => async () => { const p = await tool.browse(slot); if (p) { setFiles((f) => ({ ...f, [k]: p })); setCheck(null) } }
  const clearFile = (k) => () => { setFiles((f) => ({ ...f, [k]: '' })); setCheck(null) }
  const isBam = s.mode === 'bam'
  const parentCount = isBam ? [files.parent1_bam, files.parent2_bam].filter(Boolean).length : (files.parent_vcf ? 1 : 0)

  const config = () => {
    const num = Object.fromEntries(NUMERIC.map((k) => [k, Number(s[k])]))
    return { ...s, ...num, ...files, output_dir: files.output_dir || '' }
  }

  const problems = () => {
    const p = []
    if (isBam) {
      if (!files.ref_fasta) p.push('Select the reference genome (FASTA).')
      if (!files.high_bam || !files.low_bam) p.push('Select the BAM files of the high and the low bulk.')
    } else if (!files.high_vcf || !files.low_vcf) p.push('Select the VCF of the high bulk and the VCF of the low bulk.')
    if (!(Number(s.n_high) >= 1 && Number(s.n_low) >= 1)) p.push('Enter the number of plants in each bulk.')
    return p
  }

  const run = () => {
    const p = problems()
    if (p.length) return tool.fail(p[0])
    setProg({ pct: 0, message: 'Starting…' }); setLog([])
    tool.run(config())
  }
  const runCheck = () => {
    const p = problems()
    if (p.length) return setCheck({ ok: false, errors: p, warnings: [] })
    setCheck({ pending: true })
    bridge.check_inputs(JSON.stringify(config()))
  }
  const browseOut = async () => { const p = await tool.browse('browse_folder'); if (p) setFiles((f) => ({ ...f, output_dir: p })) }

  /* ------------------------------------------------ chromosomes */
  const chromList = r ? r.plot.chromosomes : []
  const openChrom = (name, qtlId = '') => { setChrom(name); setPickedQtl(qtlId); setTab('chromosome') }
  useEffect(() => {
    if (tab !== 'chromosome' || !r || !bridge) return
    const name = chrom || (r.qtls.length ? r.qtls[0].Chromosome : (chromList[0] && chromList[0].name))
    if (!name) return
    if (name !== chrom) setChrom(name)
    setChromErr('')
    callJson(bridge, 'get_chromosome', name).then((d) => {
      if (d && d.status === 'success') setChromData(d)
      else { setChromData(null); setChromErr((d && d.message) || 'The chromosome could not be loaded.') }
    })
  }, [tab, chrom, r, bridge]) // eslint-disable-line react-hooks/exhaustive-deps

  const download = (key, fmt = '') => bridge && bridge.download_result(key, fmt)
  const genes = useMemo(() => (r ? r.genes : []), [r])

  const modeNote = useMemo(() => {
    if (parentCount >= 2 && isBam) return 'Both parents: SNPs are kept only where the parents are confidently homozygous for different alleles.'
    if (parentCount === 1 && isBam) return 'One parent: that parent must be homozygous, and the other allele is taken from the bulks. Parental polymorphism cannot be completely validated.'
    if (!isBam && parentCount === 1) return 'A parental VCF with two samples (parent 1, parent 2) gives full parental filtering; with one sample the polymorphism cannot be completely validated.'
    return 'No parents: the analysis compares the allele frequencies of the two bulks only. Parental polymorphism cannot be checked, so sites that are not polymorphic between the parents and sequencing errors can remain.'
  }, [parentCount, isBam])

  const progressCard = (running || log.length > 0 && tool.state !== 'success') && (
    <Card title="Progress" icon={Activity}>
      <ProgressBar value={(prog.pct || 0) / 100} />
      <div style={{ margin: '4px 0 10px' }}>{prog.message}</div>
      <LogBox>{log.join('\n')}</LogBox>
      {running && <div style={{ marginTop: 12 }}><button className="btn btn-secondary" onClick={() => bridge && bridge.cancel()}><Ban size={15} /> Cancel</button></div>}
    </Card>
  )

  /* ------------------------------------------------ validation */
  const runValidation = () => { setValRes(null); setValRunning(true); bridge && bridge.run_validation(valKey, Number(valSeed) || 1) }
  const vref = valRefs.find((v) => v.key === valKey)

  return (
    <ToolPage
      id="qtlseq" icon={Dna} tabs={TABS} tab={tab} setTab={setTab}
      title="QTL-seq"
      subtitle="Bulked segregant analysis with the Δ SNP-index, simulated confidence limits and QTL intervals"
      run={tab === 'input' || tab === 'settings' ? { label: 'Run QTL-seq', runningLabel: 'Analysing…', state: tool.state, status: tool.status, onRun: run } : undefined}
    >
      {/* ========================================================= INPUT */}
      {tab === 'input' && (
        <TabPane id="input" key="input">
          <SectionHeader title="Input">Two bulks of plants with a high and a low trait value, and optionally their parents.</SectionHeader>
          <div className="note" style={{ marginTop: 0, marginBottom: 18 }}>
            <p><b>BAM files</b> (aligned reads, coordinate-sorted) or <b>VCF files</b> with allele depths (AD). In a VCF the genotype (GT) of a pooled bulk is not used: the read counts of each allele are.</p>
            <p>The SNP-index of a bulk is the share of its reads that carry the <i>index allele</i>; the Δ SNP-index is high bulk minus low bulk, so positive values mean the index allele is enriched in the high bulk.</p>
          </div>
          {progressCard}
          <RadioChips value={s.mode} onChange={(v) => { set('mode')(v); setCheck(null) }}
            options={[{ value: 'bam', label: 'BAM files' }, { value: 'vcf', label: 'VCF files' }]} />

          <div className="row" style={{ marginTop: 16 }}>
            <Card title={isBam ? 'Reference and bulks' : 'Bulks'} icon={Dna}>
              {isBam && (
                <div className="field"><label className="field-title">Reference genome (FASTA)</label>
                  <FilePicker value={files.ref_fasta} placeholder="No file selected" onBrowse={setFile('ref_fasta', 'browse_fasta')} /></div>
              )}
              <div className="field"><label className="field-title">{s.high_name || 'High bulk'} — {isBam ? 'BAM (required)' : 'VCF (required)'}</label>
                <FilePicker value={isBam ? files.high_bam : files.high_vcf} placeholder="No file selected"
                  onBrowse={setFile(isBam ? 'high_bam' : 'high_vcf', isBam ? 'browse_bam' : 'browse_vcf')} /></div>
              <div className="field"><label className="field-title">{s.low_name || 'Low bulk'} — {isBam ? 'BAM (required)' : 'VCF (required)'}</label>
                <FilePicker value={isBam ? files.low_bam : files.low_vcf} placeholder="No file selected"
                  onBrowse={setFile(isBam ? 'low_bam' : 'low_vcf', isBam ? 'browse_bam' : 'browse_vcf')} /></div>
              {!isBam && <div className="field-hint">A jointly called VCF with both bulks can be given twice; name the two samples under Advanced settings.</div>}
            </Card>
            <Card title="Parents (optional)" icon={Users} delay={0.05}>
              {isBam ? (
                <>
                  <div className="field"><label className="field-title">{s.parent1_name || 'Parent 1'} — BAM</label>
                    <FilePicker value={files.parent1_bam} placeholder="Not used" onBrowse={setFile('parent1_bam', 'browse_bam')} />
                    {files.parent1_bam && <button className="btn btn-secondary sm" style={{ marginTop: 6 }} onClick={clearFile('parent1_bam')}>Remove</button>}</div>
                  <div className="field"><label className="field-title">{s.parent2_name || 'Parent 2'} — BAM</label>
                    <FilePicker value={files.parent2_bam} placeholder="Not used" onBrowse={setFile('parent2_bam', 'browse_bam')} />
                    {files.parent2_bam && <button className="btn btn-secondary sm" style={{ marginTop: 6 }} onClick={clearFile('parent2_bam')}>Remove</button>}</div>
                </>
              ) : (
                <div className="field"><label className="field-title">Parental VCF (one or two samples)</label>
                  <FilePicker value={files.parent_vcf} placeholder="Not used" onBrowse={setFile('parent_vcf', 'browse_vcf')} />
                  {files.parent_vcf && <button className="btn btn-secondary sm" style={{ marginTop: 6 }} onClick={clearFile('parent_vcf')}>Remove</button>}</div>
              )}
              <div className={`info-line ${parentCount === 0 || parentCount === 1 && isBam ? 'warn' : ''}`}>{modeNote}</div>
            </Card>
          </div>

          <div className="row">
            <Card title="Experiment" icon={Settings2} delay={0.08}>
              <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fill,minmax(210px,1fr))', gap: '2px 14px' }}>
                <TextField label="High phenotype bulk name" value={s.high_name} onChange={set('high_name')} />
                <TextField label="Low phenotype bulk name" value={s.low_name} onChange={set('low_name')} />
                <TextField label="Plants in the high bulk" type="number" min={1} value={s.n_high} onChange={set('n_high')} />
                <TextField label="Plants in the low bulk" type="number" min={1} value={s.n_low} onChange={set('n_low')} />
                <SelectField label="Population type" value={s.population} onChange={set('population')}
                  options={[{ value: 'F2', label: 'F2' }, { value: 'RIL', label: 'RIL (recombinant inbred lines)' }, { value: 'BC', label: 'Backcross' }]} />
                {s.population === 'BC' && (
                  <SelectField label="Recurrent parent" value={s.bc_recurrent} onChange={set('bc_recurrent')}
                    options={[{ value: 'non_index', label: 'Parent without the index allele' }, { value: 'index', label: 'Parent with the index allele' }]} />
                )}
                <TextField label="Parent 1 name" value={s.parent1_name} onChange={set('parent1_name')} />
                <TextField label="Parent 2 name" value={s.parent2_name} onChange={set('parent2_name')} />
                <SelectField label="SNP-index counts the allele of" value={s.index_parent} onChange={set('index_parent')}
                  options={[{ value: 'parent2', label: s.parent2_name || 'Parent 2' }, { value: 'parent1', label: s.parent1_name || 'Parent 1' }]}
                  hint="Positive Δ then means this parent's allele is enriched in the high bulk." />
              </div>
            </Card>
            <Card title="Annotation and output" icon={FileText} delay={0.1}>
              <div className="field"><label className="field-title">Gene annotation (GFF3 / GTF, optional)</label>
                <FilePicker value={files.annotation} placeholder="Not used" onBrowse={setFile('annotation', 'browse_annotation')} />
                {files.annotation && <button className="btn btn-secondary sm" style={{ marginTop: 6 }} onClick={clearFile('annotation')}>Remove</button>}
                <div className="field-hint">With an annotation the genes inside each QTL interval are listed. They are candidates by position only.</div></div>
              <div className="field"><label className="field-title">Output folder</label>
                <FilePicker value={files.output_dir || outDefault} placeholder="Default results folder" onBrowse={browseOut} />
                <div className="field-hint">A new folder QTLseq_results is created inside it for every analysis.</div></div>
            </Card>
          </div>

          <Card title="Check before running" icon={ShieldCheck} delay={0.12}>
            {deps && (
              <div style={{ display: 'flex', flexWrap: 'wrap', gap: 10, marginBottom: 12 }}>
                {deps.map((d) => (
                  <span key={d.name} className="chip checked" title={d.message} style={{ cursor: 'default', background: d.found ? undefined : 'var(--text-muted)' }}>
                    {d.found ? <CircleCheck size={14} /> : <Wrench size={14} />} {d.name}{d.version ? ` ${d.version.replace(/^\S+\s/, '')}` : ''}: {d.found ? d.where : 'not found'}
                  </span>
                ))}
              </div>
            )}
            {deps && deps.filter((d) => !d.found).map((d) => <div className="info-line" key={d.name}>{d.message}</div>)}
            <button className="btn btn-secondary" style={{ marginTop: 10 }} onClick={runCheck} disabled={running || (check && check.pending)}><ShieldCheck size={15} /> Check inputs</button>
            {check && check.pending && <div className="info-line">Reading the files…</div>}
            {check && !check.pending && (
              <div style={{ marginTop: 6 }}>
                {check.ok && <div className="info-line" style={{ display: 'flex', gap: 8 }}><CircleCheck size={16} /> The inputs look usable.</div>}
                <Messages errors={check.errors} warnings={check.warnings} />
              </div>
            )}
          </Card>
        </TabPane>
      )}

      {/* ========================================================= SETTINGS */}
      {tab === 'settings' && (
        <TabPane id="settings" key="settings">
          <SectionHeader title="Analysis settings">The defaults are the values normally used for this analysis. Everything used is saved with the results.</SectionHeader>
          {progressCard}
          <div className="row">
            <Card title="SNP filters" icon={SlidersHorizontal}>
              <TextField label="Minimum read depth (each bulk)" type="number" min={1} value={s.min_depth} onChange={set('min_depth')} hint="SNPs covered by fewer reads in either bulk are removed." />
              <TextField label="Minimum SNP-index" type="number" min={0} max={1} step={0.05} value={s.min_snp_index} onChange={set('min_snp_index')}
                hint="A SNP is kept when at least one bulk reaches this SNP-index; the other bulk's value is kept as it is, even below the threshold." />
              <Toggle label="Apply the SNP-index filter" checked={s.apply_snp_index_filter} onChange={set('apply_snp_index_filter')} />
              <TextField label="Maximum read depth (0 = no limit)" type="number" min={0} value={s.max_depth} onChange={set('max_depth')} hint="Removes repeats and collapsed duplications with very high coverage." />
            </Card>
            <Card title="Sliding window" icon={ChartLine} delay={0.05}>
              <TextField label="Window size (bp)" type="number" min={1000} step={100000} value={s.window_size} onChange={set('window_size')} />
              <TextField label="Window step (bp)" type="number" min={100} step={1000} value={s.window_step} onChange={set('window_step')} />
              <TextField label="Minimum SNPs per window" type="number" min={1} value={s.min_snps_per_window} onChange={set('min_snps_per_window')}
                hint="Windows with fewer SNPs are left out of the curve." />
            </Card>
          </div>
          <div className="row">
            <Card title="Confidence limits (simulation)" icon={FlaskConical} delay={0.08}>
              <TextField label="Simulation replicates" type="number" min={100} step={1000} value={s.replicates} onChange={set('replicates')} hint="10,000 replicates are used for every read depth." />
              <TextField label="Random seed" type="number" value={s.seed} onChange={set('seed')} hint="The same seed and settings always give the same limits." />
              <SelectField label="Read depth in the simulation" value={s.depth_mode} onChange={set('depth_mode')}
                options={[{ value: 'auto', label: 'Automatic — the depth of every SNP' }, { value: 'fixed', label: 'One fixed depth for all SNPs' }]} />
              {s.depth_mode === 'fixed' && <TextField label="Sequencing depth" type="number" min={1} value={s.fixed_depth} onChange={set('fixed_depth')} />}
              <SelectField label="Limits of a window" value={s.ci_method} onChange={set('ci_method')}
                options={[{ value: 'snp_average', label: 'Average of the limits of its SNPs (conservative)' }, { value: 'window', label: 'Whole window simulated (more SNPs narrow the limits)' }]}
                hint="The first follows each SNP's own read depth; the second also takes the number of SNPs in the window into account." />
            </Card>
            <Card title="QTL intervals" icon={Dna} delay={0.1}>
              <SelectField label="Call a QTL where the window mean passes the" value={String(s.qtl_min_level)} onChange={(v) => set('qtl_min_level')(Number(v))}
                options={[{ value: '90', label: '90% limit (and stricter)' }, { value: '95', label: '95% limit (and stricter)' }, { value: '99', label: '99% limit only' }]} />
              <TextField label="Merge significant windows closer than (bp)" type="number" min={0} step={10000} value={s.merge_gap} onChange={set('merge_gap')}
                hint="Windows of the same sign with less than this distance of non-significant windows between them are one QTL." />
              <SelectField label="Interval start and end" value={s.interval_basis} onChange={set('interval_basis')}
                options={[{ value: 'top_level', label: 'Windows at the strictest limit the QTL reaches' }, { value: 'window_extent', label: 'All merged windows (outer ends)' }, { value: 'window_midpoints', label: 'Window midpoints (where the curve is outside the limit)' }]} />
            </Card>
          </div>
          <Toggle label="Show advanced settings" checked={advanced} onChange={setAdvanced} />
          {advanced && (
            <div className="row" style={{ marginTop: 12 }}>
              <Card title="Parents and quality" icon={Settings2}>
                <TextField label="Minimum parent read depth" type="number" min={1} value={s.min_parent_depth} onChange={set('min_parent_depth')} />
                <TextField label="Minimum parent genotype quality (Phred)" type="number" min={0} value={s.min_parent_gq} onChange={set('min_parent_gq')} hint="Confidence that a parent is homozygous." />
                <TextField label="Read error rate of the genotype caller" type="number" min={0.0001} max={0.2} step={0.005} value={s.base_error} onChange={set('base_error')} />
                <TextField label="Minimum VCF QUAL (0 = off)" type="number" min={0} value={s.min_qual} onChange={set('min_qual')} />
                <Toggle label="VCF FILTER must be PASS" checked={s.require_pass} onChange={set('require_pass')} />
              </Card>
              <Card title="BAM reading" icon={Wrench} delay={0.04}>
                <SelectField label="Read counting" value={s.engine} onChange={set('engine')}
                  options={[{ value: 'auto', label: 'Automatic (samtools, then pysam, then built-in)' }, { value: 'samtools', label: 'samtools mpileup' },
                    { value: 'pysam', label: 'pysam' }, { value: 'builtin', label: 'Built-in reader (slower)' }]} />
                <TextField label="Minimum mapping quality" type="number" min={0} value={s.min_mapq} onChange={set('min_mapq')} />
                <TextField label="Minimum base quality" type="number" min={0} value={s.min_baseq} onChange={set('min_baseq')} />
                <TextField label="Minimum reads of the other allele (no / one parent)" type="number" min={1} value={s.min_alt_reads} onChange={set('min_alt_reads')} />
                <TextField label="Minimum share of the other allele (no / one parent)" type="number" min={0} max={1} step={0.01} value={s.min_alt_fraction} onChange={set('min_alt_fraction')} />
                <TextField label="Chromosomes (comma separated, empty = all)" value={s.chromosomes} onChange={set('chromosomes')} placeholder="chr1, chr2, chr3" />
              </Card>
              <Card title="VCF sample names" icon={FileText} delay={0.08}>
                <div className="field-hint" style={{ marginTop: 0 }}>Only needed for multi-sample VCF files. Empty = the first sample (the parental VCF: first = parent 1, second = parent 2).</div>
                <TextField label="High bulk sample" value={s.high_sample} onChange={set('high_sample')} />
                <TextField label="Low bulk sample" value={s.low_sample} onChange={set('low_sample')} />
                <TextField label="Parent 1 sample" value={s.parent1_sample} onChange={set('parent1_sample')} />
                <TextField label="Parent 2 sample" value={s.parent2_sample} onChange={set('parent2_sample')} />
              </Card>
            </div>
          )}
        </TabPane>
      )}

      {/* ========================================================= RESULTS */}
      {tab === 'results' && (
        <TabPane id="results" key="results">
          <SectionHeader title="Results">{r ? `${r.summary.population} · ${r.summary.n_high} + ${r.summary.n_low} plants · ${r.summary.marker_text}` : 'Results of the last analysis.'}</SectionHeader>
          {!r ? <EmptyState icon={Info}>No results yet. Choose the input files and click Run QTL-seq.</EmptyState> : (
            <>
              <SummaryGrid items={[
                { label: 'Informative SNPs', value: r.summary.informative_snps },
                { label: 'Significant QTLs', value: r.summary.qtls },
                { label: 'QTLs at 99%', value: r.summary.qtl_99 },
                { label: 'QTLs at 95%', value: r.summary.qtl_95 },
                { label: 'QTLs at 90%', value: r.summary.qtl_90 },
                { label: 'Windows', value: r.summary.windows },
              ]} />
              <div className="info-line">SNP-index counts the allele of <b>{r.summary.index_text}</b>; Δ = {r.summary.high_name} − {r.summary.low_name}. A QTL that reaches 99% is also counted at 95% and 90%. Windows are compared with the limits one by one; the limits are not corrected for the number of windows.</div>
              <Messages warnings={r.warnings} />

              <h3 style={{ margin: '22px 0 6px', fontSize: '1rem' }}>SNP filtering</h3>
              <div className="table-wrap"><table className="results-table"><tbody>
                {[['SNPs identified', r.ledger.total_snps], ['Passing the parental polymorphism filter', r.ledger.passing_parental],
                  ['Removed for quality', r.ledger.removed_quality], ['Removed for low read depth', r.ledger.removed_low_depth],
                  ['Removed for high read depth', r.ledger.removed_high_depth], ['Removed for SNP-index below the threshold', r.ledger.removed_snp_index],
                  ['SNPs retained for the analysis', r.ledger.retained]].map(([k, v]) => <tr key={k}><td>{k}</td><td style={{ textAlign: 'right' }}>{fmtBp(v)}</td></tr>)}
              </tbody></table></div>

              <h3 style={{ margin: '24px 0 4px', fontSize: '1rem' }}>Genome-wide plots</h3>
              <SubTabs id="gw" value={view} onChange={setView} tabs={[
                { id: 'delta', label: 'Δ SNP-index' }, { id: 'high', label: r.summary.high_name }, { id: 'low', label: r.summary.low_name }]} />
              <TrackPlot chroms={r.plot.chromosomes} y={view} qtls={r.qtls} height={view === 'delta' ? 380 : 280}
                onPickChrom={(c) => openChrom(c)} title={`Genome-wide ${view}`} />
              {view === 'delta' && <LevelLegend />}
              <div style={{ display: 'flex', gap: 10, flexWrap: 'wrap', margin: '10px 0 0' }}>
                {[['genomewide', 'Δ SNP-index'], ['high_plot', 'High bulk'], ['low_plot', 'Low bulk']].map(([k, l]) => (
                  <span key={k} style={{ display: 'inline-flex', gap: 6 }}>
                    <button className="btn btn-secondary sm" onClick={() => download(k, 'png')}><Download size={13} /> {l} PNG</button>
                    <button className="btn btn-secondary sm" onClick={() => download(k, 'pdf')}>PDF</button>
                  </span>
                ))}
              </div>

              <h3 style={{ margin: '26px 0 4px', fontSize: '1rem' }}>Significant QTLs</h3>
              {r.qtls.length === 0 ? (
                <div className="info-line">No sliding window crossed the limit of the selected level, so no QTL is reported. A Δ of ±1 is not a threshold: only the simulated limits decide.</div>
              ) : (
                <>
                  <div className="field-hint" style={{ marginTop: 0 }}>Click a row to open its chromosome with the candidate genes.</div>
                  <div className="table-wrap" style={{ maxHeight: 360 }}>
                    <table className="results-table">
                      <thead><tr>{['QTL', 'Chromosome', 'Start', 'End', 'Peak', 'Peak Δ', 'Mean Δ', 'Min Δ', 'Max Δ', 'SNPs', 'Level', 'Limit', 'Direction'].map((c) => <th key={c}>{c}</th>)}</tr></thead>
                      <tbody>
                        {r.qtls.map((q) => (
                          <tr key={q.QTL_ID} onClick={() => openChrom(q.Chromosome, q.QTL_ID)} style={{ cursor: 'pointer' }}>
                            <td>{q.QTL_ID}</td><td>{q.Chromosome}</td><td>{fmtBp(q.Start)}</td><td>{fmtBp(q.End)}</td><td>{fmtBp(q.Peak_Position)}</td>
                            <td><Badge ok={q.Peak_Delta_SNP_index >= 0}>{f3(q.Peak_Delta_SNP_index)}</Badge></td>
                            <td>{f3(q.Mean_Delta_SNP_index)}</td><td>{f3(q.Minimum_Delta_SNP_index)}</td><td>{f3(q.Maximum_Delta_SNP_index)}</td>
                            <td>{fmtBp(q.SNP_Count)}</td><td>{q.Significance_Level}</td><td>{f3(q.Confidence_Threshold)}</td>
                            <td>{q.Peak_Delta_SNP_index >= 0 ? `Enriched in ${r.summary.high_name}` : `Enriched in ${r.summary.low_name}`}</td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                </>
              )}

              <h3 style={{ margin: '26px 0 8px', fontSize: '1rem' }}>Downloads</h3>
              <div style={{ display: 'flex', gap: 10, flexWrap: 'wrap' }}>
                <button className="btn btn-primary" onClick={() => download('qtls')}><Download size={16} /> QTL table (CSV)</button>
                <button className="btn btn-secondary" onClick={() => download('windows')}><Download size={16} /> All windows (CSV)</button>
                <button className="btn btn-secondary" onClick={() => download('snp_index')}><Download size={16} /> SNP-index per position (CSV)</button>
                <button className="btn btn-secondary" onClick={() => download('delta')}><Download size={16} /> Δ SNP-index (CSV)</button>
                {r.files.genes && <button className="btn btn-secondary" onClick={() => download('genes')}><Download size={16} /> Candidate genes (CSV)</button>}
                <button className="btn btn-secondary" onClick={() => download('vcf')}><Download size={16} /> Filtered SNPs (VCF.gz)</button>
                {r.files.vcf_index && <button className="btn btn-secondary" onClick={() => download('vcf_index')}><Download size={16} /> VCF index (.tbi)</button>}
                <button className="btn btn-secondary" onClick={() => download('report')}><Download size={16} /> Report (PDF)</button>
                <button className="btn btn-secondary" onClick={() => bridge && bridge.download_all('')}><Download size={16} /> Copy the whole result folder</button>
                <button className="btn btn-secondary" onClick={() => bridge && bridge.open_folder(r.out_dir)}><FolderOpen size={16} /> Open the folder</button>
              </div>
              <div className="field-hint" style={{ marginTop: 8 }}>Results folder: {r.out_dir}. It also holds the simulation table, the parameters (JSON), the log, and a PNG and a PDF of every plot.{r.summary.vcf_check ? ` ${r.summary.vcf_check}` : ''}</div>
              {log.length > 0 && <details style={{ marginTop: 14 }}><summary style={{ cursor: 'pointer' }}>Log of the run</summary><LogBox>{log.join('\n')}</LogBox></details>}
            </>
          )}
        </TabPane>
      )}

      {/* ========================================================= CHROMOSOMES */}
      {tab === 'chromosome' && (
        <TabPane id="chromosome" key="chromosome">
          <SectionHeader title="Chromosome view">SNP-index of both bulks and the Δ SNP-index of one chromosome, with the limits, the QTL intervals and the peaks.</SectionHeader>
          {!r ? <EmptyState icon={Info}>No results yet. Run the analysis first.</EmptyState> : (
            <>
              <SelectField label="Chromosome" value={chrom} onChange={(v) => { setChrom(v); setPickedQtl('') }}
                options={chromList.map((c) => ({ value: c.name, label: `${c.name} — ${c.n_snps.toLocaleString()} SNPs${r.qtls.some((q) => q.Chromosome === c.name) ? ' — QTL' : ''}` }))} />
              {chromData ? (
                <>
                  <h3 style={{ margin: '14px 0 2px', fontSize: '.95rem' }}>SNP-index of {r.summary.high_name}</h3>
                  <TrackPlot chroms={[chromData]} y="high" single qtls={chromData.qtls} height={190} title={`${chromData.name} ${r.summary.high_name}`} />
                  <h3 style={{ margin: '14px 0 2px', fontSize: '.95rem' }}>SNP-index of {r.summary.low_name}</h3>
                  <TrackPlot chroms={[chromData]} y="low" single qtls={chromData.qtls} height={190} title={`${chromData.name} ${r.summary.low_name}`} />
                  <h3 style={{ margin: '14px 0 2px', fontSize: '.95rem' }}>Δ SNP-index ({r.summary.high_name} − {r.summary.low_name})</h3>
                  <TrackPlot chroms={[chromData]} y="delta" single qtls={chromData.qtls} height={340} title={`${chromData.name} Δ SNP-index`} />
                  <LevelLegend />
                  <div style={{ display: 'flex', gap: 8, marginTop: 8 }}>
                    <button className="btn btn-secondary sm" onClick={() => bridge && bridge.download_chromosome_plot(chromData.name, 'png')}><Download size={13} /> Chromosome plot PNG</button>
                    <button className="btn btn-secondary sm" onClick={() => bridge && bridge.download_chromosome_plot(chromData.name, 'pdf')}>PDF</button>
                  </div>
                  {chromData.qtls.length > 0 && (
                    <>
                      <h3 style={{ margin: '22px 0 6px', fontSize: '1rem' }}>QTLs on {chromData.name}</h3>
                      <div className="table-wrap"><table className="results-table">
                        <thead><tr>{['QTL', 'Start', 'End', 'Peak', 'Peak Δ', 'Level', 'Limit', 'SNPs'].map((c) => <th key={c}>{c}</th>)}</tr></thead>
                        <tbody>{chromData.qtls.map((q) => (
                          <tr key={q.QTL_ID} onClick={() => setPickedQtl(q.QTL_ID)} style={{ cursor: 'pointer', fontWeight: pickedQtl === q.QTL_ID ? 700 : 400 }}>
                            <td>{q.QTL_ID}</td><td>{fmtBp(q.Start)}</td><td>{fmtBp(q.End)}</td><td>{fmtBp(q.Peak_Position)}</td>
                            <td>{f3(q.Peak_Delta_SNP_index)}</td><td>{q.Significance_Level}</td><td>{f3(q.Confidence_Threshold)}</td><td>{fmtBp(q.SNP_Count)}</td>
                          </tr>))}</tbody>
                      </table></div>
                    </>
                  )}
                  <h3 style={{ margin: '22px 0 6px', fontSize: '1rem' }}>Genes overlapping the QTL {pickedQtl ? `(${pickedQtl})` : 'intervals on this chromosome'}</h3>
                  {(() => {
                    const rows = genes.filter((g) => g.Chromosome === chromData.name && (!pickedQtl || g.QTL_ID === pickedQtl))
                    if (!r.files.genes) return <div className="info-line">No annotation file was given, so no genes are listed. Add a GFF3 or GTF on the Input tab to list the genes inside each interval.</div>
                    if (!rows.length) return <div className="info-line">No annotated gene overlaps {pickedQtl || 'a QTL interval on this chromosome'}.</div>
                    return (
                      <>
                        <div className="table-wrap" style={{ maxHeight: 340 }}><table className="results-table">
                          <thead><tr>{['QTL', 'Gene', 'Name', 'Start', 'End', 'Strand', 'Annotation', 'Distance from peak (bp)'].map((c) => <th key={c}>{c}</th>)}</tr></thead>
                          <tbody>{rows.map((g) => (
                            <tr key={g.QTL_ID + g.Gene_ID}><td>{g.QTL_ID}</td><td>{g.Gene_ID}</td><td>{g.Gene_Name || '—'}</td><td>{fmtBp(g.Gene_Start)}</td><td>{fmtBp(g.Gene_End)}</td>
                              <td>{g.Strand}</td><td style={{ maxWidth: 320, overflow: 'hidden', textOverflow: 'ellipsis' }} title={g.Annotation}>{g.Annotation || '—'}</td>
                              <td>{g.Gene_Contains_Peak ? '0 (contains the peak)' : fmtBp(g.Distance_From_Peak)}</td></tr>))}</tbody>
                        </table></div>
                        <div className="field-hint">These genes lie inside the interval; that alone does not make one of them the causal gene.</div>
                      </>
                    )
                  })()}
                </>
              ) : <div className={`info-line ${chromErr ? 'warn' : ''}`}>{chromErr || 'Loading the chromosome…'}</div>}
            </>
          )}
        </TabPane>
      )}

      {/* ========================================================= SIMULATION */}
      {tab === 'simulation' && (
        <TabPane id="simulation" key="simulation">
          <SectionHeader title="Simulation of the confidence limits">How far the Δ SNP-index can move by chance when there is no QTL.</SectionHeader>
          {!r ? <EmptyState icon={Info}>No results yet. Run the analysis first.</EmptyState> : (
            <>
              <div className="info-line">{r.plot.meta.population} · {r.plot.meta.n_high} + {r.plot.meta.n_low} plants · {r.plot.meta.replicates.toLocaleString()} replicates for each of {r.plot.meta.pairs.toLocaleString()} read-depth pair(s) · seed {r.plot.meta.seed} · limit method “{r.summary.ci_method}”.</div>
              <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit,minmax(340px,1fr))', gap: 20, marginTop: 14 }}>
                <div><h3 style={{ margin: '0 0 6px', fontSize: '.95rem' }}>Null distribution</h3><NullPlot nul={r.plot.null} /></div>
                <div><h3 style={{ margin: '0 0 6px', fontSize: '.95rem' }}>Limits by read depth</h3><DepthLimitPlot data={r.plot.limits_by_depth} /></div>
              </div>
              <LevelLegend extra={false} mean={false} />
              <h3 style={{ margin: '22px 0 6px', fontSize: '1rem' }}>Assumptions</h3>
              <div className="note"><ul style={{ margin: 0, paddingLeft: 18 }}>{r.assumptions.map((a) => <li key={a}>{a}</li>)}</ul></div>
              <div style={{ display: 'flex', gap: 8, marginTop: 12 }}>
                <button className="btn btn-secondary" onClick={() => download('simulation_plot', 'png')}><Download size={15} /> Simulation plot PNG</button>
                <button className="btn btn-secondary" onClick={() => download('simulation_plot', 'pdf')}>PDF</button>
                <button className="btn btn-secondary" onClick={() => download('simulation')}><Download size={15} /> Limits table (CSV)</button>
              </div>
            </>
          )}
        </TabPane>
      )}

      {/* ========================================================= VALIDATION */}
      {tab === 'validation' && (
        <TabPane id="validation" key="validation">
          <SectionHeader title="Validation">Checks that the pipeline finds a QTL of a reported size and position. It runs on simulated data.</SectionHeader>
          <div className="note" style={{ marginTop: 0 }}>
            <p>The reference values below are findings reported for two rice experiments. They are <b>expectations to compare with</b>, not results the program returns. The original sequencing reads are not part of TerraTools, so each experiment is <b>simulated</b> with the same design (population, plants per bulk, QTL position and size) and analysed with the same filters, simulation, windows and QTL rules as real data. A pass shows that the pipeline recovers such a QTL; it does not replicate the published numbers, which needs the original reads.</p>
          </div>
          <RadioChips value={valKey} onChange={(v) => { setValKey(v); setValRes(null) }} options={valRefs.length ? valRefs.map((v) => ({ value: v.key, label: v.title })) : [{ value: 'blast_ril', label: 'Experiment 1' }, { value: 'vigour_f2', label: 'Experiment 2' }]} />
          {vref && (
            <Card title="Reference design and findings" icon={FileText}>
              <div className="table-wrap"><table className="results-table"><tbody>
                <tr><td>Cross</td><td>{vref.cross}</td></tr>
                <tr><td>Population</td><td>{vref.population} ({vref.generation}), {vref.total_lines} plants in total</td></tr>
                <tr><td>Bulks</td><td>{vref.n_high} + {vref.n_low} plants</td></tr>
                {vref.qtls.map((q) => <tr key={q.name}><td>{q.name}</td><td>{q.chromosome}: {(q.start / 1e6).toFixed(2)}–{(q.end / 1e6).toFixed(2)} Mb · maximum Δ ≈ {q.max_delta} · {q.p}</td></tr>)}
                {vref.linkage_lod && <tr><td>Linkage mapping</td><td>LOD {vref.linkage_lod}</td></tr>}
              </tbody></table></div>
            </Card>
          )}
          <div style={{ display: 'flex', gap: 14, alignItems: 'flex-end', flexWrap: 'wrap', margin: '14px 0' }}>
            <TextField label="Random seed of the simulated data" type="number" min={1} value={valSeed} onChange={setValSeed} width={200} />
            <button className="btn btn-primary" onClick={runValidation} disabled={valRunning}><ShieldCheck size={16} /> {valRunning ? 'Running…' : 'Run validation'}</button>
          </div>
          {valRes && valRes.status === 'error' && <div className="info-line warn">{valRes.message}</div>}
          {valRes && valRes.status === 'success' && (
            <>
              <div className="info-line" style={{ display: 'flex', gap: 8 }}>{valRes.passed ? <CircleCheck size={16} /> : <TriangleAlert size={16} />}
                {valRes.passed ? 'All expected QTLs were recovered at the reported place and with the reported significance.' : 'At least one expected QTL was not recovered as expected (see below).'}</div>
              <div className="table-wrap"><table className="results-table">
                <thead><tr>{['Expected', 'Region (reported)', 'Δ (reported)', 'Called interval', 'Peak Δ', 'Level', 'Result'].map((c) => <th key={c}>{c}</th>)}</tr></thead>
                <tbody>{valRes.checks.map((c) => (
                  <tr key={c.expected}><td>{c.expected}</td><td>{c.chromosome ? `${c.chromosome}: ${c.expected_interval}` : c.expected_interval}</td>
                    <td>{c.expected_max_delta ?? ''}</td><td>{c.called_interval || ''}</td><td>{c.peak_delta ?? ''}</td>
                    <td>{c.level ? `${c.level}% (reported ${c.expected_level}%)` : ''}</td>
                    <td>{c.info ? <span title={c.note}>{c.note}</span> : <Badge ok={c.passed}>{c.passed ? 'recovered' : 'not recovered'}</Badge>}{!c.info && c.note ? ` — ${c.note}` : ''}</td></tr>))}</tbody>
              </table></div>
              <h3 style={{ margin: '20px 0 4px', fontSize: '1rem' }}>Δ SNP-index of the simulated experiment</h3>
              <TrackPlot chroms={valRes.plot.chromosomes} y="delta" qtls={valRes.qtls} height={340} />
              <LevelLegend />
              <div className="field-hint">{valRes.note}</div>
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
