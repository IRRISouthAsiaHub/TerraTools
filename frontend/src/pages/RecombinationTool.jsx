import { useEffect, useRef, useState } from 'react'
import {
  Shuffle, SlidersHorizontal, LayoutGrid, BookOpen, FileInput, FileText, FileSpreadsheet, FileArchive, Users, Settings2, Table2,
  Info, Image as ImageIcon, TriangleAlert, ListChecks, Dna,
} from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, SelectField, TextField, Slider, Toggle, FilePicker, Note, StatCard, DataTable,
  DownloadButtons, EmptyState, ImageGallery,
} from '../components/kit.jsx'
import { useBridge, useTool, callJson, mockGeneric } from '../lib/bridge.js'
import { usePipe } from '../lib/pipeline.js'
import { MultiSelect } from './favallele/Pickers.jsx'

/*
 * Recombination Analysis — crossovers, donor content and recombination hot / cold regions of a mapping population
 * (backend: tools/recombination/backend/recombination.py).
 */

const TABS = [
  { id: 'input', label: 'Input & Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
  { id: 'method', label: 'Method Guide', icon: BookOpen },
]
const NONE = '__none__'
const SMOOTHING = [
  { value: 'legacy', label: 'Original script (first short run takes the value of the run before it)' },
  { value: 'robust', label: 'Shortest run first, into its longer neighbour' },
  { value: 'none', label: 'No smoothing (every isolated call counts)' },
]
const lastName = (p) => (p.split(/[\\/]/).pop() || '')

export default function RecombinationTool() {
  const { bridge } = useBridge('recombinationBridge', () => withMock(mockGeneric()))
  const tool = useTool(bridge, 'run_recombination')
  const [tab, setTab] = useState('input')
  const [file, setFile] = useState('')
  const [info, setInfo] = useState(null)
  const [parent, setParent] = useState('')
  const [donor, setDonor] = useState(NONE)
  const [exclude, setExclude] = useState([])
  const [smoothing, setSmoothing] = useState('legacy')
  const [pct, setPct] = useState(5)
  const [floor, setFloor] = useState(3)
  const [minDonor, setMinDonor] = useState(1)
  const [bin, setBin] = useState(0.5)
  const [cen, setCen] = useState(true)
  tool.onResult(() => setTab('results'))

  const inspect = async (path) => {
    tool.setState('running'); tool.setStatus('Reading the file…')
    const r = await callJson(bridge, 'inspect_file', path, '')
    if (!r || r.status !== 'success') {
      setInfo(null); tool.setState('error'); tool.setStatus((r && r.message) || 'Could not read the file.')
      return
    }
    setInfo(r); setParent(r.guess_parent || NONE); setDonor(NONE); setExclude([])
    tool.setStatus(`${r.n_markers.toLocaleString()} markers and ${r.n_samples} samples. ${r.guess_parent ? `Is ${r.guess_parent} the recurrent parent? ` : ''}Check the settings, then run.`)
    tool.setState('idle')
  }
  const browse = async () => {
    const path = await tool.browse('browse_file')
    if (!path) return
    setFile(path); tool.setResult(null); inspect(path)
  }

  // "Open in Recombination Analysis" from the Genofile Converter
  const pipe = usePipe()
  const lastStamp = useRef(0)
  useEffect(() => {
    const g = pipe.genofile
    if (bridge && g && g.target === 'recombination' && g.stamp !== lastStamp.current) {
      lastStamp.current = g.stamp; setFile(g.recombination); tool.setResult(null); inspect(g.recombination)
    }
  }, [pipe.genofile, bridge]) // eslint-disable-line react-hooks/exhaustive-deps

  const samples = info ? info.samples : []
  const par = parent === NONE ? '' : parent
  const don = donor === NONE ? '' : donor
  const nLines = info ? samples.filter((s) => s !== par && s !== don && !exclude.includes(s)).length : 0
  const problem = !info ? 'Choose a genotype file first.'
    : par && par === don ? 'The parent and the donor must be different samples.'
    : nLines < 5 ? 'At least 5 lines are needed.' : null

  const run = () => {
    if (problem) return tool.fail(problem)
    tool.run({
      input_file: file, parent: par, donor: don, exclude_samples: exclude, smoothing, min_block_pct: pct, min_block_min: Number(floor) || 3,
      min_donor_rils: Number(minDonor) || 1, bin_mb: Number(bin) || 0.5, centromeres: cen,
    })
    tool.setStatus(`Detecting crossovers in ${nLines} lines…`)
  }

  const r = tool.result
  const stem = lastName(file).replace(/\.[^.]+$/, '') || 'recombination'
  const parentOptions = [{ value: NONE, label: '(none — the file is already coded: 0 = Parent 1, 2 = Parent 2)' }, ...samples.map((s) => ({ value: s, label: s }))]
  const donorOptions = [{ value: NONE, label: '(none — the donor allele is the other homozygous allele)' }, ...samples.filter((s) => s !== par).map((s) => ({ value: s, label: s }))]

  return (
    <ToolPage
      id="recomb" icon={Shuffle} tabs={TABS} tab={tab} setTab={setTab}
      title="Recombination Analysis"
      subtitle="Locate the crossovers of a mapping population, the donor content of every line and the recombination hot and cold regions"
      run={tab === 'input' ? {
        label: 'Run Analysis', runningLabel: 'Analysing…', state: tool.state, disabled: !info,
        status: !info && tool.state === 'idle' ? 'Choose a genotype file to begin.' : tool.status, onRun: run,
      } : undefined}
    >
      {tab === 'input' && (
        <TabPane id="input" key="input">
          <SectionHeader title="Recombination in a Mapping Population">
            From the genotypes of RILs or backcross-derived lines, TerraTools finds every <b>crossover</b> (a change between the
            recurrent-parent and the donor allele), measures the <b>donor genome</b> of each line and shows where along the chromosomes recombination happens.
          </SectionHeader>
          <div className="grid-2" style={{ gridTemplateColumns: 'minmax(320px, 1.1fr) minmax(300px, 1fr)', alignItems: 'start' }}>
            <Card title="Genotype file" icon={FileInput}>
              <FilePicker value={file} onBrowse={browse} placeholder="CSV / Excel table: one row per marker, one column per line (0 / 1 / 2 / −1)" />
              {info && (
                <div className="info-line" style={{ marginTop: 12 }}>
                  <b>{info.n_markers.toLocaleString()}</b> markers on <b>{info.chromosomes.length}</b> chromosomes, <b>{info.n_samples}</b> samples;
                  {' '}{info.het_pct}% heterozygous and {info.missing_pct}% missing calls. Positions from {info.pos_source}.
                </div>
              )}
              <Note>
                Needs the <b>position of every marker in base pairs</b> (a <i>POS</i> column, or marker names like <i>chr1_134521</i>) and genotypes coded
                <b> 0 / 1 / 2</b> with <b>−1</b> for missing. The recurrent parent can be one of the columns, or the file can come already coded
                from the Genofile Converter.
              </Note>
            </Card>
            <Card title="What is counted" icon={Info} delay={0.05}>
              <div className="code-legend">
                <div className="code-legend-row"><span className="code-chip lg c0">A</span><div><b>Recurrent parent allele</b><small>the line carries the parent's allele</small></div></div>
                <div className="code-legend-row"><span className="code-chip lg c2">B</span><div><b>Donor allele</b><small>the other homozygous allele</small></div></div>
                <div className="code-legend-row"><span className="code-chip lg c1">H</span><div><b>Heterozygous</b><small>kept, but never used to place a crossover</small></div></div>
                <div className="code-legend-row"><span className="code-chip lg cm">U</span><div><b>Missing</b><small>skipped; the crossover is placed between the flanking called markers</small></div></div>
              </div>
            </Card>
          </div>

          {info && (
            <>
              <Card title="Parent and lines" icon={Users} delay={0.03}>
                <div className="grid-3">
                  <SelectField label="Recurrent parent (allele A)" value={parent} onChange={setParent} options={parentOptions}
                    hint="The parent whose allele is called A. Markers where it is not homozygous are left out." />
                  <SelectField label="Donor sample (optional)" value={donor} onChange={setDonor} options={donorOptions}
                    hint="If the donor is in the file, markers where it carries the parent's allele are left out too." />
                  <MultiSelect label="Samples to leave out (optional)" allLabel="None — use all" clearLabel="Clear" noun="samples" search
                    options={samples.filter((s) => s !== par && s !== don).map((s) => ({ value: s, label: s }))} value={exclude} onChange={setExclude} />
                </div>
                <div className={`info-line ${problem ? 'warn' : ''}`}>
                  {problem || <><b>{nLines}</b> lines are analysed{par ? <> against <b>{par}</b></> : <> (coded file: 0 = Parent 1 = A, 2 = Parent 2 = B)</>}.</>}
                </div>
              </Card>
              <Card title="Settings" icon={Settings2} delay={0.06}>
                <div className="grid-2">
                  <div>
                    <SelectField label="Smoothing of isolated calls" value={smoothing} onChange={setSmoothing} options={SMOOTHING}
                      hint="A run of A or B calls shorter than the minimum block is treated as a genotyping error, not as a double crossover." />
                    {smoothing !== 'none' && (
                      <>
                        <Slider label="Minimum block (% of the markers of a chromosome)" value={pct} onChange={setPct} min={0} max={20} step={0.5} format={(v) => `${v} %`} />
                        <div className="field-hint" style={{ marginBottom: 10 }}>5 % of 300 markers = a block of 15 markers. A larger value removes more short donor segments.</div>
                      </>
                    )}
                  </div>
                  <div className="grid-2">
                    {smoothing !== 'none' && <TextField label="Block of at least (markers)" type="number" min={1} max={50} value={floor} onChange={setFloor} />}
                    <TextField label="Lines that must carry the other allele" type="number" min={1} max={50} value={minDonor} onChange={setMinDonor}
                      hint="A marker with fewer such lines is left out as monomorphic." />
                    <TextField label="Bin size for the rate plot (Mb)" type="number" min={0.05} max={20} step={0.05} value={bin} onChange={setBin} />
                    <div className="field">
                      <label className="field-title">Centromeres</label>
                      <Toggle label="Shade the rice centromeres" checked={cen} onChange={setCen} />
                      <div className="field-hint">IRGSP-1.0 positions of the 12 rice chromosomes.</div>
                    </div>
                  </div>
                </div>
              </Card>
              <Card title="File preview" icon={Table2} delay={0.08}>
                <DataTable rows={info.preview_rows} columns={info.preview_columns} maxHeight={260} />
              </Card>
            </>
          )}
        </TabPane>
      )}

      {tab === 'results' && (
        <TabPane id="results" key="results">
          {!r ? (
            <EmptyState icon={Shuffle}>Run the analysis first.</EmptyState>
          ) : (
            <>
              <SectionHeader title="Recombination Results">{r.message}. {r.counts.polymorphic.toLocaleString()} of {r.counts.markers.toLocaleString()} markers were informative.</SectionHeader>
              <div className="summary-grid">
                <StatCard className="summary-card" value={r.counts.lines} label="Lines" />
                <StatCard className="summary-card" value={r.counts.polymorphic} label="Polymorphic markers" delay={0.05} />
                <StatCard className="summary-card" value={r.summary['Total crossovers']} label="Crossovers" delay={0.1} />
                <StatCard className="summary-card" value={r.summary['Mean crossovers per line']} label="Mean per line" delay={0.15} />
                <StatCard className="summary-card" value={r.summary['Mean donor % (mean of the chromosomes)']} label="Mean donor (%)" delay={0.2} />
                {r.has_genes && <StatCard className="summary-card" value={r.summary['Intragenic %']} label="Intragenic (%)" delay={0.25} />}
              </div>

              {r.notes.length > 0 && (
                <Card title="Please note" icon={TriangleAlert}>
                  <ul className="geno-list">{r.notes.map((n, i) => <li key={i}>{n}</li>)}</ul>
                </Card>
              )}

              <Card title="Figures" icon={ImageIcon} delay={0.04}>
                <ImageGallery images={r.plots.map((p) => ({ src: p.path, label: p.label }))} />
              </Card>

              <Card title="Downloads" icon={FileText} delay={0.06}>
                <p className="geno-lead">
                  The Excel workbook holds every table (<i>Summary, Settings, Per_chromosome, Per_RIL, Crossovers, CO_rate_by_bin{r.has_genes ? ', Genes' : ''}</i>); the ZIP has the same tables as CSV files and all figures (PNG).
                </p>
                <DownloadButtons onSave={tool.save} items={[
                  { icon: FileSpreadsheet, label: 'All tables (Excel)', path: r.files.workbook, name: `${stem}_recombination.xlsx` },
                  { icon: FileArchive, label: 'Everything (ZIP)', path: r.files.zip, name: `${stem}_recombination.zip` },
                  { icon: FileText, label: 'Crossovers (CSV)', path: r.files.crossovers, name: `${stem}_crossovers.csv` },
                  { icon: FileText, label: 'Per line (CSV)', path: r.files.per_ril, name: `${stem}_per_line.csv` },
                  { icon: FileText, label: 'Per chromosome (CSV)', path: r.files.per_chr, name: `${stem}_per_chromosome.csv` },
                  { icon: FileText, label: 'Rate by bin (CSV)', path: r.files.bins, name: `${stem}_rate_by_bin.csv` },
                ]} />
              </Card>

              <Card title="By chromosome" icon={Table2} delay={0.08}><DataTable rows={r.per_chr} maxHeight={360} /></Card>
              <div className="grid-2">
                <Card title="Summary" icon={ListChecks} delay={0.1}><DataTable rows={r.summary_rows} maxHeight={420} /></Card>
                <Card title="Settings used" icon={Settings2} delay={0.12}><DataTable rows={r.settings} maxHeight={420} /></Card>
              </div>
              <Card title={`Lines (${r.per_ril.length} of ${r.per_ril_total}, most crossovers first)`} icon={Users}><DataTable rows={r.per_ril} maxHeight={380} /></Card>
              {r.has_genes && r.genes.length > 0 && (
                <Card title="Genes with the most intragenic crossovers" icon={Dna}><DataTable rows={r.genes} maxHeight={360} /></Card>
              )}
              <Card title={`Crossovers (first ${r.crossovers.length} of ${r.crossovers_total.toLocaleString()}; all are in the files)`} icon={Table2}>
                <DataTable rows={r.crossovers} maxHeight={420} />
              </Card>
            </>
          )}
        </TabPane>
      )}

      {tab === 'method' && (
        <TabPane id="method" key="method">
          <SectionHeader title="How the Analysis Works">Five steps, each reported in the results.</SectionHeader>
          <div className="grid-2">
            <Card title="1 · Informative markers" icon={ListChecks}>
              <p className="geno-lead">A marker is kept when the recurrent parent is <b>homozygous</b> and at least one line carries the <b>other homozygous allele</b>
                (and, if you name a donor, when the donor carries it). Markers that are the same in the parent and every line cannot show a crossover.</p>
            </Card>
            <Card title="2 · Calls" icon={Dna} delay={0.04}>
              <p className="geno-lead">Each genotype becomes <b>A</b> (as the parent), <b>B</b> (other homozygote), <b>H</b> (heterozygous) or <b>U</b> (missing).
                A file already coded by the Genofile Converter has no parent column: 0 is A and 2 is B.</p>
            </Card>
            <Card title="3 · Smoothing" icon={SlidersHorizontal} delay={0.08}>
              <p className="geno-lead">A run of A or B calls shorter than the <b>minimum block</b> (a percentage of the markers of that chromosome) is taken for a genotyping error and
                overwritten. <b>This setting changes the crossover counts a lot</b>: compare the methods and block sizes on your data before you report numbers.</p>
            </Card>
            <Card title="4 · Crossovers" icon={Shuffle} delay={0.12}>
              <p className="geno-lead">A crossover is a change between A and B in two neighbouring called markers. It lies somewhere in the interval between them;
                the table gives the interval (<i>BP_start</i>–<i>BP_end</i>), its middle and the number of markers inside it.</p>
            </Card>
            <Card title="5 · Rate, donor and genes" icon={Info} delay={0.16}>
              <p className="geno-lead"><b>Rate</b> = crossovers per line in each bin. <b>Donor %</b> = B calls among the A and B calls, per chromosome and averaged over the chromosomes.
                <b> Intragenic</b> = both flanking markers belong to the same gene, read from the marker name (<i>chr1_134521_LOC_Os01g01280</i>); it needs at least two informative markers in a gene.</p>
            </Card>
            <Card title="Reading the results" icon={TriangleAlert} delay={0.2}>
              <p className="geno-lead">Backcross-derived lines (mean donor around 20 %) recombine only inside their donor segments, so the crossovers per chromosome reflect where the donor was kept as much as the recombination rate.
                Lines with many missing calls have few locatable crossovers; leave them out.</p>
            </Card>
          </div>
        </TabPane>
      )}
    </ToolPage>
  )
}

/* Browser-preview (npm run dev) stand-ins. */
function withMock(base) {
  const samples = ['PARENT', ...Array.from({ length: 12 }, (_, i) => `RIL_${i + 1}`)]
  return new Proxy(base, {
    get(t, k) {
      if (k === 'browse_file') return (cb) => setTimeout(() => cb('C:\\data\\mapping_population.csv'), 200)
      if (k === 'inspect_file') return (_p, _s, cb) => setTimeout(() => cb(JSON.stringify({
        status: 'success', file: 'mapping_population.csv', layout: 'marker rows with a parent column', samples, n_samples: samples.length, n_markers: 4651,
        chromosomes: [1, 2, 3], pos_source: "column 'POS'", guess_parent: 'PARENT', het_pct: 1.5, missing_pct: 4.2, notes: [],
        preview_columns: ['Marker', 'Chromosome', 'Position', 'PARENT', 'RIL_1'], preview_rows: [{ Marker: 'chr1_134521_LOC_Os01g01280', Chromosome: 1, Position: 134521, PARENT: 0, RIL_1: 2 }],
      })), 300)
      return t[k]
    },
  })
}
