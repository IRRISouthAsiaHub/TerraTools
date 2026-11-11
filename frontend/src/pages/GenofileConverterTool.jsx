import { useState } from 'react'
import { motion } from 'framer-motion'
import {
  FileCog, SlidersHorizontal, LayoutGrid, BookOpen, FileInput, FileSpreadsheet, FileText, Columns3, Users, Info, Ban, ListChecks,
  Table2, ArrowRight, Dna, Link2, Settings2,
} from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, SelectField, RadioChips, Toggle, FilePicker, Note, StatCard, DataTable, DownloadButtons, EmptyState,
} from '../components/kit.jsx'
import { useBridge, useTool, callJson, mockGeneric } from '../lib/bridge.js'
import { setPipe } from '../lib/pipeline.js'
import { MultiSelect } from './favallele/Pickers.jsx'

/*
 * Genofile Converter — the start of the input-file chain.
 * Takes a genotype file as delivered by a genotyping company (HapMap, VCF / VCF.GZ or a raw allele sheet) and writes the
 * 0 / 1 / 2 / −1 coded genotype file that the other analyses read (tools/genofile_converter/backend/genofile.py).
 */

const TABS = [
  { id: 'input', label: 'Input & Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
  { id: 'formats', label: 'Formats Guide', icon: BookOpen },
]

const FORMATS = [
  { value: 'hapmap', label: 'HapMap' },
  { value: 'vcf', label: 'VCF / VCF.GZ' },
  { value: 'raw', label: 'Raw allele calls' },
]
const SEPARATORS = [
  { value: ':', label: 'Colon  —  A:G' },
  { value: '/', label: 'Slash  —  A/G' },
  { value: '|', label: 'Pipe  —  A|G' },
  { value: 'none', label: 'None  —  AG' },
]
const NONE = '__none__'
const REASON_TEXT = {
  'Invalid parent': 'A parent call is missing or unreadable',
  'Non-homozygous parent': 'A parent is heterozygous',
  'Non-segregating': 'Both parents carry the same allele',
  'Unplaced contig': 'Unplaced contig (chromosome is not a number)',
}
const CODES = [
  { code: '0', cls: 'c0', title: 'Homozygous Parent 1', example: 'same as Parent 1 (A/A)' },
  { code: '1', cls: 'c1', title: 'Heterozygous', example: 'one allele of each parent (A/G)' },
  { code: '2', cls: 'c2', title: 'Homozygous Parent 2', example: 'same as Parent 2 (G/G)' },
  { code: '-1', cls: 'cm', title: 'Missing / unexpected', example: 'no call (N, ./.) or an allele neither parent has' },
]

const codeCell = (meta) => (c, v) => {
  if (meta.includes(c)) return undefined
  const k = String(v)
  const cls = k === '0' ? 'c0' : k === '1' ? 'c1' : k === '2' ? 'c2' : 'cm'
  return <span className={`code-chip ${cls}`}>{k}</span>
}

function Bars({ items, total }) {
  return (
    <div className="hbars">
      {items.map((it, i) => (
        <div key={it.label} className="hbar">
          <div className="hbar-head"><span>{it.label}</span><b>{it.value.toLocaleString()}</b></div>
          <div className="hbar-track">
            <motion.div className={`hbar-fill ${it.cls || ''}`} initial={{ width: 0 }}
              animate={{ width: `${total ? (100 * it.value) / total : 0}%` }}
              transition={{ duration: 0.7, delay: 0.1 + i * 0.08, ease: [0.22, 1, 0.36, 1] }} />
          </div>
        </div>
      ))}
    </div>
  )
}

function Steps({ step }) {
  const labels = ['Choose the file', 'Set the parents', 'Convert', 'Use the files']
  return (
    <div className="conv-steps">
      {labels.map((l, i) => (
        <div key={l} className={`conv-step ${i < step ? 'done' : ''} ${i === step ? 'current' : ''}`}>
          <span className="num">{i + 1}</span><span>{l}</span>
          {i < labels.length - 1 && <span className="line" />}
        </div>
      ))}
    </div>
  )
}

const lastName = (p) => (p.split(/[\\/]/).pop() || '')
const stemOf = (p) => lastName(p).replace(/(\.hmp)?(\.txt)?(\.vcf)?(\.gz)?$/i, '').replace(/\.[^.]+$/, '') || 'genotypes'

export default function GenofileConverterTool({ navigate }) {
  const { bridge } = useBridge('genofileBridge', () => withMock(mockGeneric()))
  const tool = useTool(bridge, 'convert_genotypes')
  const [tab, setTab] = useState('input')
  const [file, setFile] = useState('')
  const [info, setInfo] = useState(null)
  const [fmt, setFmt] = useState('')
  const [sheet, setSheet] = useState('')
  // HapMap / VCF
  const [p1, setP1] = useState('')
  const [p2, setP2] = useState('')
  const [exclude, setExclude] = useState([])
  const [naming, setNaming] = useState('original')
  const [chromStyle, setChromStyle] = useState('number')
  const [skipUnplaced, setSkipUnplaced] = useState(true)
  // raw
  const [cols, setCols] = useState({ marker: '', chromosome: '', p1: '', p2: '' })
  const [sep, setSep] = useState(':')

  tool.onResult(() => setTab('results'))

  const inspect = async (path, format, sh) => {
    tool.setState('running'); tool.setStatus('Reading the file…')
    const r = await callJson(bridge, 'inspect_file', path, format || '', sh || '')
    if (!r || r.status !== 'success') {
      setInfo(null)
      tool.setState('error'); tool.setStatus((r && r.message) || 'Could not read the file.')
      return
    }
    setInfo(r); setFmt(r.format); setSheet(r.sheet || ''); setExclude([])
    if (r.format === 'raw') {
      setSep(r.separator || ':')
      setCols({ marker: r.guesses.marker, chromosome: r.guesses.chromosome || NONE, p1: r.guesses.p1, p2: r.guesses.p2 })
      setChromStyle('keep')
      tool.setStatus(`${r.format_label} file: ${r.columns.length} columns found. Check the column roles, then convert.`)
    } else {
      setP1(r.guesses.p1 || ''); setP2(r.guesses.p2 || '')
      setChromStyle('number')
      tool.setStatus(`${r.format_label} file with ${r.n_samples} samples. ${r.guesses.p1 && r.guesses.p2 ? 'Check the two parents, then convert.' : 'Choose the two parents, then convert.'}`)
    }
    tool.setState('idle')
  }
  const browse = async () => {
    const path = await tool.browse('browse_file')
    if (!path) return
    setFile(path); tool.setResult(null)
    inspect(path, '', '')
  }

  const raw = fmt === 'raw'
  const meta = [cols.marker, cols.chromosome, cols.p1, cols.p2].filter((c) => c && c !== NONE)
  const nOffspring = !info ? 0 : raw ? info.columns.filter((c) => !meta.includes(c) && !exclude.includes(c)).length
    : info.samples.filter((s) => s !== p1 && s !== p2 && !exclude.includes(s)).length
  const problem = !info ? 'Choose a genotype file first.'
    : raw ? (!cols.marker ? 'Choose the Marker column.'
      : !cols.p1 || !cols.p2 ? 'Choose both parent columns (P1 and P2).'
      : cols.p1 === cols.p2 ? 'P1 and P2 must be different columns.'
      : new Set(meta).size !== meta.length ? 'Each column can only have one role.'
      : !nOffspring ? 'No offspring columns left to code.' : null)
    : (!p1 || !p2 ? 'Choose both parents (Parent 1 and Parent 2).'
      : p1 === p2 ? 'Parent 1 and Parent 2 must be different samples.'
      : !nOffspring ? 'No offspring samples left to code.' : null)

  const convert = () => {
    if (problem) return tool.fail(problem)
    const base = { input_file: file, format: fmt, exclude_samples: exclude }
    tool.run(raw
      ? { ...base, sheet, separator: sep, marker_col: cols.marker, chromosome_col: cols.chromosome === NONE ? '' : cols.chromosome, p1_col: cols.p1, p2_col: cols.p2, chrom_style: chromStyle }
      : { ...base, p1_sample: p1, p2_sample: p2, marker_naming: naming, chrom_style: chromStyle, skip_unplaced: skipUnplaced })
    tool.setStatus(`Coding ${nOffspring.toLocaleString()} individuals against the two parents…`)
  }

  const r = tool.result
  const step = r ? 3 : info ? (problem ? 1 : 2) : 0
  const outName = (suffix, ext) => `${stemOf(file)}_${suffix}.${ext}`
  const toTool = (id) => {
    setPipe({ genofile: { mendelian: r.coded_csv, linkage: r.coded_file, name: lastName(file), target: id, stamp: Date.now() } })
    navigate && navigate(id)
  }
  const set = (k) => (v) => setCols((c) => ({ ...c, [k]: v }))
  const colOptions = info && raw ? info.columns : []
  const sampleOptions = info && !raw ? info.samples : []
  const pick = (list) => [{ value: '', label: 'Select…' }, ...list.map((s) => ({ value: s, label: s }))]

  return (
    <ToolPage
      id="geno" icon={FileCog} tabs={TABS} tab={tab} setTab={setTab}
      title="Genofile Converter"
      subtitle="Turn the genotype files of a genotyping company (HapMap, VCF or raw alleles) into the input files of the other TerraTools analyses"
      run={tab === 'input' ? {
        label: 'Convert Genotype File', runningLabel: 'Converting…', state: tool.state, disabled: !info,
        status: !info && tool.state === 'idle' ? 'Choose a genotype file to begin.' : tool.status, onRun: convert,
      } : undefined}
    >
      {tab === 'input' && (
        <TabPane id="input" key="input">
          <SectionHeader title="Convert a Genotype File">
            Genotyping companies send alleles (<b>A/G</b>, <b>AG</b>, <b>0/1</b>). The analyses read codes <b>0 / 1 / 2 / −1</b> measured against the two parents.
            Convert your file here first: the coded file you get is the starting point for the other analyses.
          </SectionHeader>
          <Steps step={step} />

          <div className="grid-2" style={{ gridTemplateColumns: 'minmax(320px, 1.1fr) minmax(300px, 1fr)', alignItems: 'start' }}>
            <Card title="Genotype file" icon={FileInput}>
              <FilePicker value={file} onBrowse={browse} placeholder="HapMap (.hmp.txt), VCF (.vcf, .vcf.gz) or an Excel / CSV sheet of allele calls" />
              {info && (
                <div className="field" style={{ marginTop: 14 }}>
                  <label className="field-title">File format {info.format === fmt ? '(detected)' : ''}</label>
                  <RadioChips options={FORMATS} value={fmt} onChange={(v) => { setFmt(v); setInfo(null); inspect(file, v, '') }} />
                </div>
              )}
              {info && raw && info.sheets.length > 1 && (
                <SelectField label="Sheet" value={sheet} onChange={(s) => { setSheet(s); tool.setResult(null); inspect(file, 'raw', s) }} options={info.sheets} />
              )}
              <Note>
                The format is detected from the file; change it above if it is wrong. HapMap and VCF files list the <b>parents as two of their samples</b>;
                a raw sheet holds them in two columns, <b>P1</b> and <b>P2</b>.
              </Note>
            </Card>

            <Card title="How calls are coded" icon={Info} delay={0.05}>
              <div className="code-legend">
                {CODES.map((c) => (
                  <div key={c.code} className="code-legend-row">
                    <span className={`code-chip lg ${c.cls}`}>{c.code}</span>
                    <div><b>{c.title}</b><small>{c.example}</small></div>
                  </div>
                ))}
              </div>
              <div className="legend-excl">
                <Ban size={14} />
                <span>Markers are <b>left out</b> when a parent is missing, a parent is heterozygous,
                or both parents carry the same allele — they cannot tell which parent an allele came from.</span>
              </div>
            </Card>
          </div>

          {info && !raw && (
            <>
              <Card title="Parents and samples" icon={Users} delay={0.03}>
                <div className="grid-4">
                  <SelectField label="Parent 1 (coded 0)" value={p1} onChange={setP1} options={pick(sampleOptions)} />
                  <SelectField label="Parent 2 (coded 2)" value={p2} onChange={setP2} options={pick(sampleOptions)} />
                  <MultiSelect label="Samples to leave out (optional)" allLabel="None — use all" clearLabel="Clear" noun="samples" search
                    options={sampleOptions.filter((s) => s !== p1 && s !== p2).map((s) => ({ value: s, label: s }))} value={exclude} onChange={setExclude} />
                </div>
                <div className={`info-line ${problem ? 'warn' : ''}`}>
                  {problem || <>The file has <b>{info.n_samples}</b> samples: the two parents and <b>{nOffspring}</b> individuals that are coded against them.</>}
                </div>
              </Card>
              <Card title="Output settings" icon={Settings2} delay={0.06}>
                <div className="grid-3">
                  <SelectField label="Marker names" value={naming} onChange={setNaming} options={[
                    { value: 'original', label: 'As in the file (chr_position if empty)' },
                    { value: 'chrpos', label: 'chromosome_position  (chr01_1044946)' }]}
                    hint="Names like chr01_1044946 carry the physical position of the marker, which some analyses use." />
                  <SelectField label="Chromosome labels" value={chromStyle} onChange={setChromStyle} options={[
                    { value: 'number', label: 'Numbers only  (Chr01 → 1)' },
                    { value: 'keep', label: 'As in the file' }]}
                    hint={info.example_chromosomes && info.example_chromosomes.length ? `In your file: ${info.example_chromosomes.join(', ')}` : undefined} />
                  <div className="field">
                    <label className="field-title">Unplaced contigs</label>
                    <Toggle label="Leave out markers on scaffolds" checked={skipUnplaced && chromStyle === 'number'} onChange={setSkipUnplaced} />
                    <div className="field-hint">Only with numbered chromosomes: scaffolds and organelle contigs are not mapped.</div>
                  </div>
                </div>
              </Card>
            </>
          )}

          {info && raw && (
            <Card title="Match the columns" icon={Columns3} delay={0.05}>
              <div className="grid-4">
                <SelectField label="Marker" value={cols.marker} onChange={set('marker')} options={colOptions} />
                <SelectField label="Chromosome" value={cols.chromosome} onChange={set('chromosome')}
                  options={[{ value: NONE, label: '(none — infer from marker names)' }, ...colOptions]} />
                <SelectField label="Parent 1 (P1)" value={cols.p1} onChange={set('p1')} options={pick(colOptions)} />
                <SelectField label="Parent 2 (P2)" value={cols.p2} onChange={set('p2')} options={pick(colOptions)} />
                <SelectField label="Allele separator" value={sep} onChange={setSep} options={SEPARATORS} />
              </div>
              <div className={`info-line ${problem ? 'warn' : ''}`}>
                {problem || <>All other <b>{nOffspring}</b> columns are treated as individuals (columns with no genotype calls are skipped automatically).</>}
              </div>
            </Card>
          )}

          {info && (
            <Card title="File preview" icon={Table2} delay={0.08}>
              <DataTable rows={info.preview_rows} columns={info.preview_columns} maxHeight={260} />
            </Card>
          )}
        </TabPane>
      )}

      {tab === 'results' && (
        <TabPane id="results" key="results">
          {!r ? (
            <EmptyState icon={FileCog}>Convert a genotype file first.</EmptyState>
          ) : (
            <>
              <SectionHeader title="Conversion Result">
                {r.parents[0]} (Parent 1) × {r.parents[1]} (Parent 2) — {r.n_offspring.toLocaleString()} individuals coded.
              </SectionHeader>
              <div className="summary-grid">
                <StatCard className="summary-card" value={r.total_markers} label="Markers in file" />
                <StatCard className="summary-card" value={r.informative} label="Coded (kept)" delay={0.05} />
                <StatCard className="summary-card" value={r.excluded} label="Left out" delay={0.1} />
                <StatCard className="summary-card" value={r.n_offspring} label="Individuals" delay={0.15} />
                <StatCard className="summary-card" value={r.missing_pct} label="Missing / unexpected (%)" delay={0.2} />
              </div>

              <Card title="Your converted file" icon={FileText}>
                {r.informative === 0 && <div className="info-line warn">No marker could be coded. Check that the two parents are the right samples.</div>}
                {r.duplicate_names > 0 && (
                  <div className="info-line warn">{r.duplicate_names} marker name{r.duplicate_names > 1 ? 's occur' : ' occurs'} more than once in the file; use <b>chromosome_position</b> marker names if the analyses complain.</div>
                )}
                <p className="geno-lead">
                  One table with <b>a row per marker</b> ({r.informative.toLocaleString()}) and <b>a column per individual</b> ({r.n_offspring.toLocaleString()}). The first two columns are
                  <b> Marker</b> and <b>Chromosome</b>; every other column holds the genotype of one individual, coded against the parents
                  {' '}<b>{r.parents[0]}</b> and <b>{r.parents[1]}</b>:
                </p>
                <div className="geno-codes">
                  {CODES.map((c) => <span key={c.code}><span className={`code-chip ${c.cls}`}>{c.code}</span> {c.title}</span>)}
                </div>
                <p className="geno-lead">The parents themselves are not in the table. Markers that cannot be coded are listed below with the reason.</p>
                <ul className="geno-list">
                  <li><b>CSV</b> — just this table.</li>
                  <li><b>Excel workbook</b> — {r.big ? 'the summary and the marker reports (the table itself is too big for a sheet; use the CSV)' : 'the same table on the first sheet, then'}: <i>Marker_Details</i> (missing / unexpected calls per marker), <i>Excluded_Markers</i> (every marker left out and why) and <i>Summary</i> (the settings and counts).</li>
                </ul>
                <DownloadButtons onSave={tool.save} items={[
                  { icon: FileText, label: 'Coded genotypes (CSV)', path: r.coded_csv, name: outName('coded', 'csv') },
                  { icon: FileSpreadsheet, label: r.big ? 'Reports (Excel)' : 'Coded genotypes + reports (Excel)', path: r.big ? r.report_xlsx : r.coded_file, name: outName(r.big ? 'report' : 'coded', 'xlsx') },
                ]} />
                <div className="geno-next">
                  <span>Continue with this file in</span>
                  <button type="button" className="btn btn-secondary sm" disabled={!r.informative} onClick={() => toTool('mendelian')}><Dna size={14} /> Mendelian Segregation</button>
                  <button type="button" className="btn btn-secondary sm" disabled={!r.informative} onClick={() => toTool('linkage')}><Link2 size={14} /> Linkage Mapping</button>
                </div>
              </Card>

              <div className="grid-2">
                <Card title="Why markers were left out" icon={Ban} delay={0.05}>
                  <Bars total={Math.max(r.excluded, 1)}
                    items={Object.entries(r.reasons).filter(([k, v]) => v > 0 || k !== 'Unplaced contig').map(([k, v]) => ({ label: REASON_TEXT[k] || k, value: v, cls: 'warn' }))} />
                </Card>
                <Card title="Genotype codes across all kept markers" icon={ListChecks} delay={0.1}>
                  <Bars total={Object.values(r.code_counts).reduce((a, b) => a + b, 0)}
                    items={CODES.map((c) => ({ label: `${c.code}  ·  ${c.title}`, value: r.code_counts[c.code] || 0, cls: c.cls }))} />
                </Card>
              </div>
              {r.skipped_columns.length > 0 && <div className="info-line">Skipped (no genotype calls): {r.skipped_columns.join(', ')}</div>}

              <Card title={`Coded markers (first ${r.preview_rows.length} of ${r.informative.toLocaleString()})`} icon={Table2}>
                <DataTable rows={r.preview_rows} columns={r.preview_columns} renderCell={codeCell(['Marker', 'Chromosome'])} maxHeight={380} />
              </Card>
              {r.excluded > 0 && (
                <Card title={`Markers left out (${r.excluded.toLocaleString()}${r.excluded > r.excluded_rows.length ? `, first ${r.excluded_rows.length} shown` : ''})`} icon={Ban}>
                  <DataTable rows={r.excluded_rows} columns={['Marker', 'Chromosome', 'P1', 'P2', 'Reason']} maxHeight={320} />
                </Card>
              )}
            </>
          )}
        </TabPane>
      )}

      {tab === 'formats' && (
        <TabPane id="formats" key="formats">
          <SectionHeader title="Supported Formats">What the three kinds of genotype files look like, and what is made of them.</SectionHeader>
          <div className="grid-2">
            <Card title="HapMap" icon={FileText}>
              <p className="geno-lead">Tab-separated, one row per marker. The 11 standard columns, then one column per sample.</p>
              <pre className="geno-code">{'rs#  alleles  chrom  pos  strand … QCcode  P1  P2  F2_1  F2_2\nm1   A/G      chr01  717702  +  …  NA      AA  GG  AG    AA'}</pre>
              <p className="geno-lead">Calls are two letters (<b>AA</b>, <b>AG</b>), a single letter (<b>A</b>, or the IUPAC code <b>R</b> for A/G) or <b>NN</b> for missing.</p>
            </Card>
            <Card title="VCF / VCF.GZ" icon={FileText} delay={0.04}>
              <p className="geno-lead">Variant Call Format, also gzipped. The GT field of every sample is read.</p>
              <pre className="geno-code">{'#CHROM  POS     ID   REF ALT … FORMAT  P1   P2   F2_1\nchr01   717702  m1   A   G   … GT      0/0  1/1  0/1'}</pre>
              <p className="geno-lead"><b>0/0</b>, <b>0/1</b>, <b>1/1</b> (also with <b>|</b>) and <b>./.</b> for missing; the genotype is turned into the REF / ALT alleles first, so multi-allelic markers work.</p>
            </Card>
            <Card title="Raw allele calls" icon={FileSpreadsheet} delay={0.08}>
              <p className="geno-lead">An Excel or CSV sheet, one row per marker; the parents have a column each.</p>
              <pre className="geno-code">{'Marker  Chromosome  P1   P2   F2_1  F2_2\nm1      1           A:A  G:G  A:G   A:A'}</pre>
              <p className="geno-lead">Alleles may be written <b>A:G</b>, <b>A/G</b>, <b>A|G</b> or <b>AG</b>; the separator is detected.</p>
            </Card>
            <Card title="What you get" icon={ListChecks} delay={0.12}>
              <p className="geno-lead">One coded table, as a CSV and as an Excel workbook:</p>
              <ul className="geno-list">
                <li><b>Rows</b> are markers, <b>columns</b> are <i>Marker</i>, <i>Chromosome</i>, then one column per individual.</li>
                <li><b>Cells</b> are <b>0</b> (Parent 1 allele), <b>1</b> (heterozygous), <b>2</b> (Parent 2 allele) or <b>−1</b> (missing / unexpected).</li>
                <li>The <b>Excel workbook</b> adds the sheets <i>Marker_Details</i>, <i>Excluded_Markers</i> and <i>Summary</i>.</li>
              </ul>
              <Note>The parents are not part of the table; they define what 0 and 2 mean. Any TerraTools analysis that reads a 0 / 1 / 2 genotype table can use it.</Note>
            </Card>
          </div>
        </TabPane>
      )}
    </ToolPage>
  )
}

/* Browser-preview (npm run dev) stand-ins for the converter slots. */
function withMock(base) {
  const samples = ['P1', 'P2', ...Array.from({ length: 12 }, (_, i) => `F2_${i + 1}`)]
  const calls = (i, j) => ['0/0', '0/1', '1/1', '0/1', './.'][(i * 3 + j) % 5]
  const rows = Array.from({ length: 6 }, (_, i) => ({ Marker: `chr01_${717702 + i * 1900}`, Chromosome: 'chr01', Position: 717702 + i * 1900, REF: 'A', ALT: 'G', P1: '0/0', P2: '1/1', F2_1: calls(i, 0), F2_2: calls(i, 1), F2_3: calls(i, 2), F2_4: calls(i, 3) }))
  const code = { '0/0': 0, '0/1': 1, '1/1': 2, './.': -1 }
  return new Proxy(base, {
    get(t, k) {
      if (k === 'browse_file') return (cb) => setTimeout(() => cb('C:\\data\\genotypes.vcf.gz'), 200)
      if (k === 'inspect_file') return (_p, _f, _s, cb) => setTimeout(() => cb(JSON.stringify({
        status: 'success', format: 'vcf', format_label: 'VCF', file: 'genotypes.vcf.gz', samples, n_samples: samples.length, guesses: { p1: 'P1', p2: 'P2' },
        preview_columns: ['Marker', 'Chromosome', 'Position', 'REF', 'ALT', 'P1', 'P2', 'F2_1', 'F2_2', 'F2_3', 'F2_4'], preview_rows: rows, example_chromosomes: ['chr01', 'chr02'],
      })), 300)
      if (k === 'convert_genotypes') return () => setTimeout(() => t.resultReady.emit(JSON.stringify({
        status: 'success', format: 'vcf', coded_csv: 'C:\tmp\genotypes_coded.csv', coded_file: 'C:\tmp\genotypes_coded.xlsx', report_xlsx: 'C:\tmp\genotypes_coded.xlsx', big: false,
        total_markers: 1840, informative: 1512, excluded: 328, n_offspring: 12, missing_pct: 6.4, skipped_columns: [], duplicate_names: 0, parents: ['P1', 'P2'],
        reasons: { 'Invalid parent': 41, 'Non-homozygous parent': 97, 'Non-segregating': 190, 'Unplaced contig': 0 },
        code_counts: { 0: 4610, 1: 8920, 2: 3450, '-1': 1164 },
        preview_columns: ['Marker', 'Chromosome', 'F2_1', 'F2_2', 'F2_3', 'F2_4'],
        preview_rows: rows.map((r) => ({ Marker: r.Marker, Chromosome: 1, F2_1: code[r.F2_1], F2_2: code[r.F2_2], F2_3: code[r.F2_3], F2_4: code[r.F2_4] })),
        excluded_rows: [{ Marker: 'chr01_1385', Chromosome: 1, P1: '0/0', P2: '0/0', Reason: 'Non-segregating' }],
      })), 1200)
      return t[k]
    },
  })
}
