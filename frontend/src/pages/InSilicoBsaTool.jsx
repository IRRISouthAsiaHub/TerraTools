import { useState } from 'react'
import {
  FlaskConical, SlidersHorizontal, LayoutGrid, BookOpen, FileInput, FileText, FileSpreadsheet, FileArchive, Settings2, Table2,
  Info, Image as ImageIcon, TriangleAlert, ListChecks, Dna, Users, Sigma, Search,
} from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, SelectField, TextField, Slider, FilePicker, Note, StatCard, DataTable,
  DownloadButtons, EmptyState, ImageGallery,
} from '../components/kit.jsx'
import { useBridge, useTool, callJson, mockGeneric } from '../lib/bridge.js'

/*
 * In silico BSA — bulked segregant analysis from the genotypes and phenotypes of every line
 * (backend: tools/insilico_bsa/backend/insilico_bsa.py). Statistic: Gw of dQTG-seq (Li et al. 2022, Plant Communications 3:100319).
 */

const TABS = [
  { id: 'input', label: 'Input & Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
  { id: 'method', label: 'Method Guide', icon: BookOpen },
]
const POPULATIONS = [
  { value: '', label: 'Choose the population type…' },
  { value: 'f2', label: 'F2 / F3 / IMF2 (three genotype classes, 1:2:1)' },
  { value: 'bc1', label: 'Backcross BC1 (two genotype classes, 1:1)' },
  { value: 'ril', label: 'RIL (inbred lines, one allele per line)' },
  { value: 'dh', label: 'DH (doubled haploid, one allele per line)' },
]
const POP_NAME = { f2: 'F2 / IMF2', bc1: 'backcross', ril: 'RIL / DH', dh: 'RIL / DH' }
const PERMS = [200, 500, 1000, 2000, 5000].map((v) => ({ value: String(v), label: `${v.toLocaleString()} permutations` }))
const ALPHAS = [
  { value: '0.1', label: '0.10 (suggested for dQTG-seq1-type data)' },
  { value: '0.05', label: '0.05 (genome-wide)' },
  { value: '0.01', label: '0.01 (stringent)' },
]
const PRIMARY = [{ value: '250', label: '250 kb' }, { value: '500', label: '500 kb' }, { value: '1000', label: '1 Mb' }]
const REGION_COLS = ['Region', 'Window', 'Chromosome', 'Left_marker', 'Left_bp', 'Peak_marker', 'Peak_bp', 'Right_marker', 'Right_bp', 'P_adj_Gw_at_peak', 'LOD', 'Delta_SNP_index']
const WINDOW_COLS = ['Window', 'Chromosome', 'Window_start', 'Window_end', 'Scored_markers', 'Gw', 'G1', 'abs_dAF', 'P_adj_Gw', 'Sig_Gw']
const lociCols = (comp) => ['Locus', 'Chromosome', 'Start_bp', 'End_bp', 'Span_kb', 'Peak_marker', 'Windows_that_found_it', 'Peak_spread_kb',
  ...comp.flatMap((c) => [`Length_kb_${c.Window}`, `P_adj_${c.Window}`])]
const NEAR_COLS = ['Region', 'Window', 'Chromosome', 'Left_marker', 'Left_bp', 'Peak_marker', 'Peak_bp', 'Right_marker', 'Right_bp', 'P_adj_Gw_at_peak', 'LOD', 'Fraction_of_threshold', 'Delta_SNP_index']
const lastName = (p) => (p.split(/[\\/]/).pop() || '')

export default function InSilicoBsaTool() {
  const { bridge } = useBridge('inSilicoBsaBridge', () => withMock(mockGeneric()))
  const tool = useTool(bridge, 'run_insilico_bsa')
  const [tab, setTab] = useState('input')
  const [geno, setGeno] = useState('')
  const [pheno, setPheno] = useState('')
  const [info, setInfo] = useState(null)
  const [trait, setTrait] = useState('')
  const [population, setPopulation] = useState('')
  const [mergeKb, setMergeKb] = useState(1000)
  const [primary, setPrimary] = useState('500')
  const [minWin, setMinWin] = useState(0)
  const [annWin, setAnnWin] = useState('')
  const [gff, setGff] = useState('')
  const [flank, setFlank] = useState(100)
  const [genes, setGenes] = useState(null)
  const [annBusy, setAnnBusy] = useState(false)
  const [annErr, setAnnErr] = useState('')
  const [frac, setFrac] = useState(20)
  const [minCalls, setMinCalls] = useState(5)
  const [maf, setMaf] = useState(0.05)
  const [perms, setPerms] = useState('1000')
  const [alpha, setAlpha] = useState('0.05')
  tool.onResult(() => { setGenes(null); setAnnErr(''); setAnnWin(''); setTab('results') })

  const inspect = async (g, p) => {
    tool.setState('running'); tool.setStatus('Reading the files…')
    const r = await callJson(bridge, 'inspect_inputs', g, p || '')
    if (!r || r.status !== 'success') {
      setInfo(null); tool.setState('error'); tool.setStatus((r && r.message) || 'Could not read the files.')
      return
    }
    setInfo(r)
    setTrait((t) => (r.traits.includes(t) ? t : (r.traits[0] || '')))
    tool.setStatus(`${r.n_markers.toLocaleString()} markers and ${r.n_samples} lines${r.traits.length ? '' : '; choose the phenotype file'}.`)
    tool.setState('idle')
  }
  const browseGeno = async () => {
    const path = await tool.browse('browse_geno')
    if (!path) return
    setGeno(path); tool.setResult(null); inspect(path, pheno)
  }
  const browsePheno = async () => {
    const path = await tool.browse('browse_pheno')
    if (!path) return
    setPheno(path); tool.setResult(null); if (geno) inspect(geno, path)
  }

  const ready = info && trait && population
  const k = info ? Math.max(2, Math.round((frac / 100) * (pheno ? info.n_matched : info.n_samples))) : 0
  const problem = !info ? 'Choose a genotype file first.'
    : !info.traits.length ? 'Choose a phenotype file: the genotype file has no trait.'
    : (pheno ? info.n_matched : info.n_samples) < 20 ? 'At least 20 lines with a genotype and a trait value are needed.'
    : k < Number(minCalls) ? `Each pool would hold only ${k} lines; raise the pool fraction.`
    : !population ? 'Choose the population type: it decides how the heterozygous calls are used.' : null

  const browseGff = async () => {
    const path = await tool.browse('browse_gff')
    if (path) { setGff(path); setGenes(null); setAnnErr('') }
  }
  const annotate = async () => {
    if (!gff) { setAnnErr('Choose the annotation file (GFF3 / GTF) first.'); return }
    setAnnBusy(true); setAnnErr('')
    const out = await callJson(bridge, 'annotate_regions', JSON.stringify({ regions_file: tool.result.files.regions, window: winSel, gff_file: gff, flank_kb: Number(flank) || 0 }))
    setAnnBusy(false)
    if (!out || out.status !== 'success') { setGenes(null); setAnnErr((out && out.message) || 'The annotation failed.'); return }
    setGenes(out)
  }

  const run = () => {
    if (problem) return tool.fail(problem)
    tool.run({
      input_file: geno, pheno_file: pheno, trait, population, merge_gap_kb: Number(mergeKb) >= 0 ? Number(mergeKb) : 0, primary_window_kb: Number(primary), min_window_markers: Number(minWin) || 0, pool_fraction: frac / 100, min_calls: Number(minCalls) || 5,
      min_maf: Number(maf), permutations: Number(perms), alpha: Number(alpha),
    })
    tool.setStatus(`Pooling the ${k} lowest and ${k} highest lines and running ${Number(perms).toLocaleString()} permutations for the 250 kb, 500 kb and 1 Mb windows…`)
  }

  const r = tool.result
  const winSel = annWin || ((r && r.comparison && r.comparison.find((c) => Number(c.Window_kb) === Number(primary))) || {}).Window || ''
  const stem = lastName(geno).replace(/\.[^.]+$/, '') || 'insilico_bsa'

  return (
    <ToolPage
      id="insilico_bsa" icon={FlaskConical} tabs={TABS} tab={tab} setTab={setTab}
      title="In silico BSA"
      subtitle="Bulked segregant analysis from the genotypes and the trait of every line: the extreme lines form the pools, no pooled sequencing needed"
      run={tab === 'input' ? {
        label: 'Run Analysis', runningLabel: 'Analysing…', state: tool.state, disabled: !ready,
        status: !info && tool.state === 'idle' ? 'Choose a genotype file to begin.' : tool.status, onRun: run,
      } : undefined}
    >
      {tab === 'input' && (
        <TabPane id="input" key="input">
          <SectionHeader title="Bulked Segregant Analysis Without Pooled Sequencing">
            TerraTools sorts your lines by the trait, takes the lowest and the highest fraction as the <b>low and high pool</b>, and compares the
            genotypes of the two pools at every marker. It uses the <b>Gw statistic of dQTG-seq</b>, which counts alleles <i>and</i> genotypes, so it also
            finds <b>over-dominant</b> loci that allele-frequency methods (G′, ED, Δ SNP-index) cannot see.
          </SectionHeader>
          <div className="grid-2" style={{ gridTemplateColumns: 'minmax(320px, 1fr) minmax(320px, 1fr)', alignItems: 'start' }}>
            <Card title="Genotype file" icon={FileInput}>
              <FilePicker value={geno} onBrowse={browseGeno} placeholder="CSV / Excel: one row per marker, one column per line (0 / 1 / 2 / −1)" />
              {info && (
                <div className="info-line" style={{ marginTop: 12 }}>
                  <b>{info.n_markers.toLocaleString()}</b> markers on <b>{info.chromosomes.length}</b> chromosomes, <b>{info.n_samples}</b> lines;
                  {' '}{info.het_pct}% heterozygous and {info.missing_pct}% missing calls (looks like {POP_NAME[info.suggest_population] || 'an F2'}; choose the population type below).
                  {info.pos_source ? <> Positions are read from the {info.pos_source}.</> : null}
                </div>
              )}
              <Note>
                Needs a <b>CHROM</b> and a <b>POS</b> column; positions are read from them, never from the marker name (a <b>Marker</b> column is only the name). Marker names such as
                <i> chr1_134521</i> are used for the position only when there is no POS column. Genotypes are coded <b>0 / 1 / 2</b>
                with <b>−1</b> for missing: 0 and 2 are the two parents, 1 is heterozygous. The coded files of the Genofile Converter work. A GWAS-ready file
                (one row per line with a trait column) works too, and then needs no phenotype file.
              </Note>
            </Card>
            <Card title="Phenotype file" icon={FileText} delay={0.05}>
              <FilePicker value={pheno} onBrowse={browsePheno} placeholder="CSV / Excel: sample name in the first column, then the traits" />
              {info && info.traits.length > 0 && (
                <div className="info-line" style={{ marginTop: 12 }}>
                  Trait from the {info.pheno_from}: <b>{info.traits.length}</b> column{info.traits.length > 1 ? 's' : ''}
                  {pheno ? <>; <b>{info.n_matched}</b> of {info.n_samples} line names match</> : null}.
                </div>
              )}
              {info && info.traits.length > 1 && (
                <SelectField label="Trait" value={trait} onChange={setTrait} options={info.traits.map((t) => ({ value: t, label: t }))} />
              )}
              <Note>The first column holds the line names, exactly as they are written in the genotype file. Lines without a trait value are left out.</Note>
            </Card>
          </div>

          {info && (
            <>
              <Card title="Pools and population" icon={Users} delay={0.03}>
                <div className="grid-2">
                  <div>
                    <Slider label="Lines in each pool (% of all lines)" value={frac} onChange={setFrac} min={5} max={50} step={1} format={(v) => `${v} %`} />
                    <div className={`info-line ${problem ? 'warn' : ''}`}>
                      {problem || <><b>{k}</b> lowest and <b>{k}</b> highest lines of {pheno ? info.n_matched : info.n_samples} form the pools.</>}
                    </div>
                  </div>
                  <SelectField label="Population type" value={population} onChange={setPopulation} options={POPULATIONS}
                    hint={`Choose it yourself. Heterozygous calls are used for F2 / IMF2 and BC1 (the genotype counts are what lets Gw see over-dominance); in RIL and DH every line carries one allele, heterozygous calls are treated as missing and Gw equals G1.${info.suggest_population ? ` The data look like ${POP_NAME[info.suggest_population]} (${info.het_pct}% heterozygous).` : ''}`} />
                </div>
              </Card>
              <Card title="Windows and threshold" icon={Settings2} delay={0.06}>
                <p className="geno-lead">
                  The genome is cut into non-overlapping windows of <b>250 kb, 500 kb and 1 Mb</b>, all three analysed in the same run and compared. A window holds the mean Gw of its scored markers;
                  its threshold comes from permutations (random low and high pools of the same size, the largest window value of the genome each time).
                </p>
                <div className="grid-3">
                  <SelectField label="Primary window" value={primary} onChange={setPrimary} options={PRIMARY}
                    hint="Used for the headline threshold, the genotype figure and the G1 / ΔAF comparison. Every window size is still reported." />
                  <TextField label="Minimum scored markers per window" type="number" min={0} max={100} value={minWin} onChange={setMinWin}
                    hint="0 = automatic: 3, or fewer when the map is sparse (a 300-marker map has about one marker per window, so a fixed 3 would test almost nothing). Windows below the minimum are not tested." />
                  <TextField label="Merge windows closer than (kb)" type="number" min={0} max={50000} step={100} value={mergeKb} onChange={setMergeKb}
                    hint="Significant windows of a chromosome form one region when the gap between them is at most this (0: touching windows only). Choose it from marker density and LD decay; the Settings sheet lists the median marker spacing." />
                  <TextField label="Minimum called lines per pool" type="number" min={2} max={200} value={minCalls} onChange={setMinCalls}
                    hint="Markers with fewer calls in a pool are not scored." />
                  <SelectField label="Permutations" value={perms} onChange={setPerms} options={PERMS}
                    hint="Two pools of the same size are drawn at random from all lines each time; the same draws serve the three window sizes." />
                  <SelectField label="Significance level" value={alpha} onChange={setAlpha} options={ALPHAS} />
                  <TextField label="Minimum minor allele frequency" type="number" min={0} max={0.5} step={0.01} value={maf} onChange={setMaf}
                    hint="Rare markers are left out before the pools are compared." />
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
            <EmptyState icon={FlaskConical}>Run the analysis first.</EmptyState>
          ) : (
            <>
              <SectionHeader title="In silico BSA Results">{r.message}. {r.counts.markers.toLocaleString()} markers, {r.counts.lines} lines, {r.counts.permutations.toLocaleString()} permutations.</SectionHeader>
              <div className="summary-grid">
                <StatCard className="summary-card" value={r.counts.lines} label="Lines" />
                <StatCard className="summary-card" value={r.counts.pool} label="Lines per pool" delay={0.05} />
                <StatCard className="summary-card" value={r.counts.markers} label="Markers" delay={0.1} />
                <StatCard className="summary-card" value={r.loci.length} label="Loci" delay={0.15} />
                <StatCard className="summary-card" value={r.counts.regions} label="Regions (all sizes)" delay={0.18} />
                <StatCard className="summary-card" value={String(r.summary['Gw threshold'])} label="Gw threshold (primary)" delay={0.2} />
              </div>

              {r.notes.length > 0 && (
                <Card title="Please note" icon={TriangleAlert}>
                  <ul className="geno-list">{r.notes.map((n, i) => <li key={i}>{n}</li>)}</ul>
                </Card>
              )}

              <Card title="Loci across the window sizes" icon={Dna} delay={0.025}>
                {r.loci.length === 0
                  ? <p className="geno-lead">No window reached the threshold at any size. Try a larger pool fraction or the 0.10 level.</p>
                  : (
                    <>
                      <p className="geno-lead">
                        A locus is a group of regions from the different window sizes that overlap. <i>Windows_that_found_it</i> says how many sizes reach the threshold there; a locus found by all three
                        is the most reliable. <i>Peak_spread_kb</i> is how far apart the peak markers of the sizes lie, and the length and adjusted p-value of each size show the resolution it gives.
                      </p>
                      <DataTable rows={r.loci} columns={lociCols(r.comparison)} maxHeight={320} />
                    </>
                  )}
              </Card>

              <Card title="Significant regions of every window size" icon={Dna} delay={0.03}>
                {r.regions.length === 0
                  ? <p className="geno-lead">No region reached the threshold at any window size. The regions closest to significance are listed below.</p>
                  : (
                    <>
                      <p className="geno-lead">
                        One row per region and window size (the <i>Window</i> column says which size the region comes from). The left and right marker are the first and last scored marker of the region, the peak marker is the strongest one; each is followed by its position in bp.
                        <i>P_adj_Gw_at_peak</i> is the permutation-adjusted p-value of the peak window and <i>LOD</i> is its Gw / ln(10) (G here is the log likelihood ratio, so this is a LOD score). <i>Delta_SNP_index</i> is the allele-frequency difference at the peak marker, high pool minus low pool, for the allele of parent 2 (positive: that allele is more frequent in the high pool).
                        Positions, lengths, the genotype pattern, allele frequencies and the G1 / ΔAF check are in the Excel report.
                      </p>
                      <DataTable rows={r.regions} columns={REGION_COLS} maxHeight={320} />
                    </>
                  )}
              </Card>

              {r.near_regions && r.near_regions.length > 0 && (
                <Card title="Closest to significant (not significant)" icon={Search} delay={0.032}>
                  <p className="geno-lead">
                    {r.regions.length === 0 ? 'No region reached the threshold. ' : 'Some window sizes found no significant region. '}
                    For each window size without one, these are the regions <b>closest to significance</b>, in the same form as the table above (left, peak and right marker with their positions,
                    adjusted p-value, LOD and ΔSNP-index). <i>Fraction_of_threshold</i> is the window Gw divided by the threshold of that window size; 1 would be on the threshold.
                    They are <b>not significant</b>: treat them as leads to follow up, not as findings.
                  </p>
                  <DataTable rows={r.near_regions} columns={NEAR_COLS} maxHeight={320} />
                </Card>
              )}

              {r.regions.length > 0 && (
                <Card title="Candidate genes (separate step)" icon={Search} delay={0.03}>
                  <p className="geno-lead">
                    Lists the genes of a genome annotation that overlap each region of one window size, widened by a flank, nearest to the peak marker first. This gives positional candidates, not the causal gene.
                    Use the annotation (GFF3 / GTF) of the same genome assembly as the marker positions.
                  </p>
                  <div className="grid-3" style={{ alignItems: 'end' }}>
                    <SelectField label="Regions of" value={winSel} onChange={(v) => { setAnnWin(v); setGenes(null) }}
                      options={r.comparison.map((c) => ({ value: c.Window, label: `${c.Window} windows (${c.Regions} regions)` }))} />
                    <FilePicker value={gff} onBrowse={browseGff} placeholder="GFF3 / GTF of the reference genome (.gff3, .gtf, .gz)" />
                    <Slider label="Flank on each side of a region" value={flank} onChange={setFlank} min={0} max={1000} step={10} format={(v) => `${v} kb`} />
                  </div>
                  <button className="btn btn-primary" style={{ marginTop: 12 }} onClick={annotate} disabled={annBusy}>
                    <Search size={15} /> {annBusy ? 'Reading the annotation…' : 'Annotate regions'}
                  </button>
                  {annErr && <div className="info-line warn" style={{ marginTop: 12 }}>{annErr}</div>}
                  {genes && (
                    <div style={{ marginTop: 14 }}>
                      <div className="info-line"><b>{genes.message}.</b> {genes.genes_total > genes.genes.length ? `The first ${genes.genes.length} are shown; all are in the downloads.` : ''}</div>
                      <DataTable rows={genes.per_region} columns={['Region', 'Chromosome', 'Start_bp', 'End_bp', 'Peak_marker', 'Genes', 'Inside']} maxHeight={220} />
                      <div style={{ height: 10 }} />
                      <DataTable rows={genes.genes} columns={['Region', 'Gene_ID', 'Gene_Name', 'Gene_Start', 'Gene_End', 'Strand', 'Location', 'Distance_to_peak_bp', 'Description']} maxHeight={360} />
                      <DownloadButtons onSave={tool.save} items={[
                        { icon: FileSpreadsheet, label: 'Candidate genes (Excel)', path: genes.files.xlsx, name: `${stem}_candidate_genes.xlsx` },
                        { icon: FileText, label: 'Candidate genes (CSV)', path: genes.files.csv, name: `${stem}_candidate_genes.csv` },
                      ]} />
                    </div>
                  )}
                </Card>
              )}

              <Card title="Figures" icon={ImageIcon} delay={0.04}>
                <ImageGallery images={r.plots.map((p) => ({ src: p.path, label: p.label }))} />
              </Card>

              <Card title="Downloads" icon={FileText} delay={0.06}>
                <p className="geno-lead">
                  The Excel workbook holds <i>Summary, Settings, Window comparison, Loci, Regions, Thresholds, Windows, Pools</i> and <i>Markers</i> (every statistic of every marker); the ZIP has the same tables as CSV files and all figures (PNG).
                </p>
                <DownloadButtons onSave={tool.save} items={[
                  { icon: FileSpreadsheet, label: 'All tables (Excel)', path: r.files.workbook, name: `${stem}_insilico_bsa.xlsx` },
                  { icon: FileArchive, label: 'Everything (ZIP)', path: r.files.zip, name: `${stem}_insilico_bsa.zip` },
                  { icon: FileText, label: 'Window comparison (CSV)', path: r.files.comparison, name: `${stem}_window_comparison.csv` },
                  { icon: FileText, label: 'Loci (CSV)', path: r.files.loci, name: `${stem}_loci.csv` },
                  ...(r.files.near ? [{ icon: FileText, label: 'Closest to significant (CSV)', path: r.files.near, name: `${stem}_closest_to_significant.csv` }] : []),
                  { icon: FileText, label: 'Regions (CSV)', path: r.files.regions, name: `${stem}_regions.csv` },
                  { icon: FileText, label: 'Window statistics (CSV)', path: r.files.windows, name: `${stem}_window_statistics.csv` },
                  { icon: FileText, label: 'Marker statistics (CSV)', path: r.files.markers, name: `${stem}_marker_statistics.csv` },
                  { icon: FileText, label: 'Pool members (CSV)', path: r.files.pools, name: `${stem}_pool_members.csv` },
                ]} />
              </Card>

              <div className="grid-2">
                <Card title="Summary" icon={ListChecks} delay={0.08}><DataTable rows={r.summary_rows} maxHeight={420} /></Card>
                <Card title="Thresholds (permutation, per window size)" icon={Sigma} delay={0.1}><DataTable rows={r.thresholds} columns={['Statistic', 'Alpha', 'Threshold']} maxHeight={420} /></Card>
              </div>
              <Card title="Windows with the highest Gw" icon={Table2} delay={0.12}><DataTable rows={r.windows_top} columns={WINDOW_COLS} maxHeight={360} /></Card>
              <Card title="Markers with the highest Gw (single markers, for reference)" icon={Table2} delay={0.13}><DataTable rows={r.top_markers} columns={['Marker', 'Chromosome', 'Position', 'Gw', 'G1', 'G2', 'dAF', 'Het_low', 'Het_high']} maxHeight={320} /></Card>
              <div className="grid-2">
                <Card title="Settings used" icon={Settings2}><DataTable rows={r.settings} maxHeight={360} /></Card>
                <Card title={`Lines and pools (${r.pools.length} of ${r.pools_total})`} icon={Users}><DataTable rows={r.pools} columns={['Sample', 'Trait', 'Pool']} maxHeight={360} /></Card>
              </div>
            </>
          )}
        </TabPane>
      )}

      {tab === 'method' && (
        <TabPane id="method" key="method">
          <SectionHeader title="How the Analysis Works">Six steps; comparing three window sizes shows which findings do not depend on the window you chose.</SectionHeader>
          <div className="grid-2">
            <Card title="1 · Pools" icon={Users}>
              <p className="geno-lead">The lines are sorted by the trait. The lowest and the highest fraction (20 % by default) are the <b>low and high pool</b>. Because the genotype of every
                line is known, allele and genotype counts are counts of lines, not sequencing reads (the way Li et al. analysed a real rice IMF2 population).</p>
            </Card>
            <Card title="2 · Counts" icon={Dna} delay={0.04}>
              <p className="geno-lead">At every marker the two pools give a <b>pool × allele</b> table (2 × 2) and a <b>pool × genotype</b> table (2 × 3). Missing calls are skipped;
                markers with too few calls in a pool, or a rare minor allele, are not scored. A marker that cannot be scored gets <b>no statistic</b> (it is not counted as zero) and does not enter any window.</p>
            </Card>
            <Card title="3 · Statistics" icon={Sigma} delay={0.08}>
              <p className="geno-lead"><b>G1</b> is the G statistic of the allele table, <b>G2</b> that of the genotype table. <b>Gw</b> = G1/(G1+G2)·G1 + G2/(G1+G2)·G2.
                A locus with <i>no</i> allele-frequency difference (fully over-dominant) has G1 near 0 but a large G2, so only Gw finds it. TerraTools also reports G1 and the
                allele-frequency difference (the Δ SNP-index) next to Gw.</p>
            </Card>
            <Card title="4 · Windows" icon={SlidersHorizontal} delay={0.12}>
              <p className="geno-lead">The genome is cut into <b>non-overlapping physical windows of 250 kb, 500 kb and 1 Mb</b>. The value of a window is the <b>mean of the marker statistic over the scored markers</b>
                it holds; a window with fewer than 3 scored markers is not tested. There is no sliding kernel: a marker belongs to exactly one window of each size, and the three sizes are
                reported next to each other.</p>
            </Card>
            <Card title="5 · Permutation threshold per size" icon={Info} delay={0.16}>
              <p className="geno-lead">Gw has no known distribution, so two pools of the same size are drawn <b>at random from all lines</b> many times, the statistics are recomputed and the
                <b> largest window value of the genome</b> is kept, T<sub>b</sub> = max<sub>w</sub> Gw<sub>wb</sub>, separately for each window size (the same draws serve all three). The threshold is the (1 − α)
                quantile of the T<sub>b</sub>, and P<sub>adj</sub> = (1 + number of permutations with T<sub>b</sub> ≥ Gw) / (B + 1) controls the family-wise error across the windows of that size.</p>
            </Card>
            <Card title="6 · Regions, loci and peaks" icon={TriangleAlert} delay={0.2}>
              <p className="geno-lead">Significant windows that touch (or lie within the merge distance) form a <b>region</b>; its peak window has the highest value and its peak marker is the strongest marker in it.
                Regions of the three sizes that overlap are one <b>locus</b>, and the table shows which sizes found it, how long each interval is and how far apart the peaks lie. A locus found by all
                sizes does not depend on the window; one found only by the smallest window rests on few markers, one found only by the largest may be too coarse to place.
                The <i>Pattern</i> column reads the genotype classes at the peak: an <b>allele-frequency shift</b> (additive or dominant), or <b>heterozygotes enriched</b> in one pool without an allele
                shift (over-dominant). Candidate genes are a separate step: choose the window size, the GFF3 / GTF of your reference genome and the flank.</p>
            </Card>
            <Card title="When it is not suitable" icon={TriangleAlert} delay={0.24}>
              <p className="geno-lead">The permutation test treats the lines as exchangeable, which fits an F2, IMF2, backcross or RIL population from one cross. In a diversity panel with
                family structure, related lines in the same pool give false signals everywhere; use GWAS with a kinship matrix there. A very polygenic trait gives many regions; compare
                with QTL Mapping before you pick candidates.</p>
            </Card>
          </div>
          <Card title="Reference" icon={BookOpen} delay={0.22}>
            <p className="geno-lead">Li P, Li G, Zhang Y-W, Zuo J-F, Liu J-Y, Zhang Y-M (2022). A combinatorial strategy to identify various types of QTLs for quantitative traits using
              extreme phenotype individuals in an F2 population. <i>Plant Communications</i> 3:100319. doi:10.1016/j.xplc.2022.100319. Magwene PM, Willis JH, Kelly JK (2011). The statistics of
              bulk segregant analysis using next generation sequencing. <i>PLoS Computational Biology</i> 7:e1002255.</p>
          </Card>
        </TabPane>
      )}
    </ToolPage>
  )
}

/* Browser-preview (npm run dev) stand-ins. */
function withMock(base) {
  return new Proxy(base, {
    get(t, k) {
      if (k === 'browse_geno') return (cb) => setTimeout(() => cb('C:\\data\\genotypes_F2.csv'), 200)
      if (k === 'browse_gff') return (cb) => setTimeout(() => cb('C:\\data\\genome.gff3'), 200)
      if (k === 'annotate_regions') return (_c, cb) => setTimeout(() => cb(JSON.stringify({ status: 'error', message: 'Preview mode: no annotation is read.' })), 200)
      if (k === 'browse_pheno') return (cb) => setTimeout(() => cb('C:\\data\\phenotypes.csv'), 200)
      if (k === 'run_insilico_bsa') return (_c, _cb) => setTimeout(() => t.resultReady.emit(JSON.stringify(MOCK_RESULT)), 600)
      if (k === 'inspect_inputs') return (_g, _p, cb) => setTimeout(() => cb(JSON.stringify({
        status: 'success', file: 'genotypes_F2.csv', layout: 'marker rows', n_samples: 300, n_markers: 900, chromosomes: [1, 2, 3], pos_source: "column 'POS'",
        het_pct: 49.1, missing_pct: 3.0, suggest_population: 'f2', notes: [], traits: ['PlantHeight'], pheno_from: 'phenotype file', n_matched: 300,
        preview_columns: ['Marker', 'Chromosome', 'Position', 'F2-001'], preview_rows: [{ Marker: 'chr1_100000', Chromosome: 1, Position: 100000, 'F2-001': 1 }],
      })), 300)
      return t[k]
    },
  })
}

const MOCK_RESULT = {
  status: 'success', message: '2 loci for PlantHeight (250kb: 2, 500kb: 2, 1Mb: 1 regions; 60 lines per pool)',
  counts: { lines: 300, markers: 900, pool: 60, regions: 5, permutations: 1000 }, notes: [],
  summary: { 'Gw threshold': 5.4437 },
  comparison: [
    { Window: '250kb', Window_kb: 250, Windows_tested: 240, Median_markers_per_window: 4, Gw_threshold: 5.9279, Significant_windows: 31, Regions: 2, Total_region_length_Mb: 8.5, Median_region_length_kb: 4250, Smallest_P_adj_Gw: 0.005 },
    { Window: '500kb', Window_kb: 500, Windows_tested: 120, Median_markers_per_window: 7.5, Gw_threshold: 5.4437, Significant_windows: 17, Regions: 2, Total_region_length_Mb: 8.5, Median_region_length_kb: 4250, Smallest_P_adj_Gw: 0.005 },
    { Window: '1Mb', Window_kb: 1000, Windows_tested: 60, Median_markers_per_window: 15, Gw_threshold: 5.2707, Significant_windows: 9, Regions: 1, Total_region_length_Mb: 5, Median_region_length_kb: 5000, Smallest_P_adj_Gw: 0.005 },
  ],
  loci: [
    { Locus: 'L1', Chromosome: 3, Start_bp: 7000000, End_bp: 12000000, Span_kb: 5000, Peak_marker: 'chr3_9100000', Windows_that_found_it: '3 of 3', Peak_spread_kb: 0,
      'Length_kb_250kb': 5000, 'P_adj_250kb': 0.005, 'Length_kb_500kb': 5000, 'P_adj_500kb': 0.005, 'Length_kb_1Mb': 5000, 'P_adj_1Mb': 0.005 },
    { Locus: 'L2', Chromosome: 1, Start_bp: 1000000, End_bp: 4500000, Span_kb: 3500, Peak_marker: 'chr1_2300000', Windows_that_found_it: '2 of 3', Peak_spread_kb: 120,
      'Length_kb_250kb': 3500, 'P_adj_250kb': 0.009, 'Length_kb_500kb': 3500, 'P_adj_500kb': 0.005, 'Length_kb_1Mb': null, 'P_adj_1Mb': null },
  ],
  regions: [
    { Region: '500kb-R1', Window: '500kb', Chromosome: 3, Start_bp: 7000000, End_bp: 12000000, Length_kb: 5000, Windows_above_threshold: 10, Peak_marker: 'chr3_9100000', Peak_bp: 9100000, Peak_value: 23.1, Markers: 74, P_adj_Gw_at_peak: 0.005, Pattern: 'allele frequency shift (additive or dominant)' },
  ],
  windows_top: [{ Window: '500kb', Chromosome: 3, Window_start: 9000000, Window_end: 9500000, Scored_markers: 8, Gw: 23.1, G1: 22.0, abs_dAF: 0.41, P_adj_Gw: 0.005, Sig_Gw: true }],
  top_markers: [{ Marker: 'chr3_9100000', Chromosome: 3, Position: 9100000, Gw: 26.2, G1: 24.8, G2: 1.4, dAF: 0.43, Het_low: 0.4, Het_high: 0.5 }],
  thresholds: [{ Statistic: 'Gw (alleles + genotypes), 500kb windows', Alpha: 0.05, Threshold: 5.4437 }],
  summary_rows: [{ Metric: 'Loci', Value: 2 }], settings: [{ Setting: 'Window sizes', Value: '250kb, 500kb, 1Mb' }],
  pools: [{ Sample: 'F2-001', Trait: 80.1, Pool: 'Low' }], pools_total: 120, plots: [], files: { regions: 'C:\tmp\r.csv', workbook: '', zip: '', markers: '', pools: '', comparison: '', loci: '', windows: '' },
}
