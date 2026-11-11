import { useState } from 'react'
import { GitBranch, FolderOpen, ChartColumn, Table2, Dna, Database, Search, Download, FileText, Target } from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, FilePicker, TextField, Note,
  EmptyState, StatCard, ImageGallery, DataTable, RadioChips, SelectField,
} from '../components/kit.jsx'
import { useBridge, useTool, callJson, useSignal, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'configure', label: 'Configure', icon: FolderOpen },
  { id: 'plan', label: 'Optimal Plan', icon: Target },
]

const MAP_KINDS = [
  { value: 'map', label: 'Genetic map (cM)' },
  { value: 'vcf', label: 'No map: use a VCF' },
]

const dash = (c, v) => (v === null || v === undefined ? '—' : undefined)

export default function CrossSelectionTool() {
  const { bridge } = useBridge('crossBridge', () => mockGeneric(['errorOccurred']))
  const tool = useTool(bridge, 'run_cross_selection')
  const [tab, setTab] = useState('configure')
  const [geno, setGeno] = useState('')
  const [effects, setEffects] = useState('')
  const [gmap, setGmap] = useState('')
  const [mapKind, setMapKind] = useState('map')
  const [cmPerMb, setCmPerMb] = useState('')
  const [validation, setValidation] = useState(null)
  const [nCross, setNCross] = useState(20)
  const [maxUse, setMaxUse] = useState(4)
  const [angle, setAngle] = useState(30)
  const [cross, setCross] = useState({ cross_type: 'DH', generation: 1, prop_sel: 0.05, pool_size: 60 })
  const [dataset, setDataset] = useState({})
  const setC = (k) => (v) => setCross((c) => ({ ...c, [k]: v }))

  useSignal(bridge, 'errorOccurred', (m) => tool.fail('Error: ' + m))

  tool.onResult((r) => {
    const s = (r.ocs && r.ocs.summary) || {}
    setDataset({ rows: s.n_markers, genos: s.n_genotypes, overlap: s.n_candidates })
    setTab('plan')
  })

  const preview = async () => {
    const r = await callJson(bridge, 'preview_markers', geno, effects, gmap, mapKind, String(cmPerMb ?? ''))
    setValidation(r)
    if (r && r.status !== 'error') setDataset({ rows: r.markers, genos: r.genotypes })
  }
  const run = () => {
    if (!geno || !effects || !gmap) return tool.fail(`Please select the allele dosage matrix, the marker effects and the ${mapKind === 'vcf' ? 'VCF' : 'genetic map'}.`)
    if (mapKind === 'vcf' && !(Number(cmPerMb) > 0)) return tool.fail('Enter the recombination rate (cM per Mb) used to turn the VCF positions into a genetic map.')
    tool.run({
      geno_file: geno, effects_file: effects, map_file: gmap, map_kind: mapKind, cm_per_mb: mapKind === 'vcf' ? Number(cmPerMb) : null,
      ...cross, n_crosses: Number(nCross), max_crosses_per_parent: maxUse, penalty_degree: angle,
    })
  }

  const r = tool.result
  const o = (r && r.ocs) || null
  const os = (o && o.summary) || {}
  const dl = (type) => bridge && bridge.download_result(type)
  const fmt = (v) => v ?? 0
  const pickTo = (slot, set) => async () => { const p = await tool.browse(slot); if (p) set(p) }

  return (
    <ToolPage
      id="cross" icon={GitBranch} tabs={TABS} tab={tab} setTab={setTab}
      title="Cross Selection"
      subtitle="Optimal cross selection — the crosses with the highest usefulness, balanced against diversity"
      run={tab === 'configure' ? { label: 'Run Cross Selection', runningLabel: 'Optimizing crosses…', state: tool.state, status: tool.status, onRun: run } : undefined}
    >
      {tab === 'configure' && (
        <TabPane id="configure" key="configure">
          <div className="grid-2" style={{ gridTemplateColumns: 'repeat(auto-fit, minmax(320px, 1fr))', alignItems: 'start' }}>
            <Card title="Input Files" icon={FolderOpen}>
              <div className="field">
                <label className="field-title">Allele Dosage Matrix CSV</label>
                <FilePicker value={geno} placeholder="Select the dosage matrix…" onBrowse={pickTo('browse_geno', setGeno)} />
                <div className="field-hint">Genotypes in the rows, markers in the columns, values 0 / 1 / 2 (inbred parents: 0 and 2; a heterozygous call is read as missing). The relationship matrix is calculated from it.</div>
              </div>
              <div className="field">
                <label className="field-title">Marker Effects CSV</label>
                <FilePicker value={effects} placeholder="Select the marker effects…" onBrowse={pickTo('browse_effects', setEffects)} />
                <div className="field-hint">Columns Marker and Effect (the additive effect of each marker, for example from Genomic Prediction).</div>
              </div>
              <div className="field">
                <label className="field-title">Genetic positions of the markers</label>
                <RadioChips options={MAP_KINDS} value={mapKind} onChange={(v) => { setMapKind(v); setGmap(''); setValidation(null) }} />
              </div>
              {mapKind === 'map' ? (
                <div className="field">
                  <label className="field-title">Genetic Map CSV</label>
                  <FilePicker value={gmap} placeholder="Select the genetic map…" onBrowse={pickTo('browse_map', setGmap)} />
                  <div className="field-hint">Columns Marker, Chromosome and cM. Needed to predict the variance of the progeny from the recombination between markers.</div>
                </div>
              ) : (
                <>
                  <div className="field">
                    <label className="field-title">VCF file</label>
                    <FilePicker value={gmap} placeholder="Select the VCF (.vcf or .vcf.gz)…" onBrowse={pickTo('browse_vcf', setGmap)} />
                    <div className="field-hint">Its CHROM, POS and ID columns give the position of every marker; the ID (or CHROM_POS) must be the marker name used in the other files.</div>
                  </div>
                  <TextField label="Recombination rate (cM per Mb)" type="number" min={0.01} max={50} step={0.1} value={cmPerMb} onChange={setCmPerMb}
                    placeholder="for example 4"
                    hint="Base-pair positions are turned into cM at this one rate for the whole genome (roughly 4 in rice, 1 in maize, 0.2 in wheat; use the value of your species). A real genetic map is more accurate." />
                </>
              )}
              <button className="btn btn-secondary" style={{ width: '100%', marginTop: 16 }} onClick={preview}><Search size={16} /> Preview & Validate</button>
              {validation && (validation.status === 'error'
                ? <div className="info-line warn">{validation.message}</div>
                : <>
                    <div className="info-line">Dosage matrix: {validation.genotypes} genotypes × {validation.markers} markers on {validation.chromosomes} chromosomes · relationship matrix {validation.kinship}</div>
                    {(validation.notes || []).map((n) => <div className="info-line warn" key={n}>{n}</div>)}
                  </>)}
            </Card>
            <Card title="Optimal Cross Selection" icon={Target} delay={0.03}>
              <TextField label="Number of crosses" type="number" min={1} max={500} value={nCross} onChange={setNCross} />
              <Slider label="Max crosses per parent" value={maxUse} onChange={setMaxUse} min={1} max={8} />
              <Slider label="Penalty degree (usefulness ↔ diversity)" value={angle} onChange={setAngle} min={0} max={90} step={1} format={(v) => `${v}°`} />
              <SelectField label="Population made from each cross" value={cross.cross_type} onChange={setC('cross_type')}
                options={[{ value: 'DH', label: 'Doubled haploids (DH)' }, { value: 'RIL', label: 'Recombinant inbred lines (RIL)' }]} />
              <TextField label="Generation" type="number" min={1} max={20} value={cross.generation} onChange={setC('generation')}
                hint="1 = DH made from the F1, or RILs taken from the F2. Use 11 or more for fully inbred RILs." />
              <TextField label="Selected proportion" type="number" min={0.001} max={0.9} step={0.01} value={cross.prop_sel} onChange={setC('prop_sel')}
                hint="The share of the progeny you will keep (0.05 = the best 5 %); it sets the selection intensity in the usefulness." />
              <TextField label="Candidate pool" type="number" min={4} max={300} value={cross.pool_size} onChange={setC('pool_size')}
                hint="The usefulness is calculated for every cross among this many candidates, taken in order of genomic value." />
              <Note>
                The plan of crosses maximises the <b>usefulness</b> of its crosses — the mean of the progeny plus the selection intensity times its standard deviation, from the
                marker effects and the recombination along the genetic map — while keeping the group coancestry and the inbreeding of the progeny down.
                0° looks at usefulness only, 90° at diversity only; the Optimal Plan tab shows the trade-off at every angle.
              </Note>
            </Card>
            <Card title="Dataset Information" icon={Database} delay={0.06}>
              <div className="summary-grid" style={{ gridTemplateColumns: 'repeat(3, 1fr)' }}>
                {[['Markers', dataset.rows], ['Genotypes', dataset.genos], ['Candidates', dataset.overlap]].map(([l, v]) => (
                  <div className="summary-card" key={l}><div className="val">{v ?? '—'}</div><div className="lbl">{l}</div></div>
                ))}
              </div>
              <Note>
                Candidates are the best genotypes by genomic value (dosage × marker effects). Usefulness, cross mean and progeny SD come from the genetic map; the
                Kinship Coefficient comes from the relationship matrix calculated from the dosage matrix.
              </Note>
            </Card>
          </div>
        </TabPane>
      )}
      {tab === 'plan' && (
        <TabPane id="plan" key="plan">
          <SectionHeader title="Optimal Crossing Plan" />
          {!o ? <EmptyState icon={Target}>Run an analysis to get the optimal plan.</EmptyState> : (
            <>
              <div className="summary-grid" style={{ gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))' }}>
                <StatCard className="summary-card" value={fmt(os.n_crosses)} label="Crosses" />
                <StatCard className="summary-card" value={fmt(os.n_parents)} label={`Parents (of ${fmt(os.n_candidates)})`} delay={0.05} />
                <StatCard className="summary-card" value={Number(fmt(os.expected_gain))} label="Expected usefulness (cross mean)" delay={0.1} />
                <StatCard className="summary-card" value={Number(fmt(os.group_coancestry))} label="Group coancestry" delay={0.15} />
                <StatCard className="summary-card" value={Number(fmt(os.progeny_inbreeding))} label="Progeny inbreeding" delay={0.2} />
              </div>
              <div className="info-line">
                Truncation selection (the best parents, mated at random) would give usefulness {fmt(os.truncation_gain)}, group coancestry {fmt(os.truncation_coancestry)}{' '}
                and progeny inbreeding {fmt(os.truncation_inbreeding)}. Penalty degree used: {fmt(os.penalty_degree)}°; mean usefulness of all possible crosses: {fmt(os.mean_usefulness)}.
                {' '}Population: {os.cross_type}, generation {os.generation}, best {Math.round((os.prop_sel || 0) * 100 * 10) / 10} % selected; {os.n_markers} markers;
                relationship matrix {os.kinship_source}; map: {os.map_source}.
              </div>
              {os.engine && <div className="info-line">Method: {os.engine}.</div>}
              {os.note && <div className="info-line">{os.note}</div>}
              <div className="downloads">
                <button className="btn btn-secondary" onClick={() => dl('ocs_crosses_csv')}><FileText size={15} /> Crossing Plan CSV</button>
                <button className="btn btn-secondary" onClick={() => dl('ocs_contributions_csv')}><FileText size={15} /> Parent Contributions CSV</button>
                <button className="btn btn-secondary" onClick={() => dl('ocs_frontier_csv')}><FileText size={15} /> Usefulness / Diversity Table CSV</button>
                <button className="btn btn-secondary" onClick={() => dl('ocs_frontier_plot')}><Download size={15} /> Frontier Plot</button>
                <button className="btn btn-secondary" onClick={() => dl('ocs_usefulness_csv')}><FileText size={15} /> Usefulness of All Crosses CSV</button>
                <button className="btn btn-secondary" onClick={() => dl('ocs_kinship_csv')}><FileText size={15} /> Relationship Matrix CSV</button>
              </div>
              <div className="grid-2" style={{ gridTemplateColumns: 'repeat(auto-fit, minmax(420px, 1fr))', alignItems: 'start' }}>
                <Card title="Crosses to make" icon={Table2}>
                  <DataTable rows={o.crosses} renderCell={dash} maxHeight={420}
                    columns={['Cross_ID', 'Parent1', 'Parent2', 'Usefulness', 'Cross_Mean', 'Progeny_SD', 'Kinship_Coefficient']} />
                </Card>
                <Card title="Parent contributions" icon={Dna} delay={0.06}>
                  <DataTable rows={o.contributions} maxHeight={420} columns={['Genotype', 'Genomic_Value', 'N_Crosses', 'Contribution']} />
                </Card>
              </div>
              <div className="grid-2" style={{ gridTemplateColumns: 'repeat(auto-fit, minmax(420px, 1fr))', alignItems: 'start', marginTop: 16 }}>
                <Card title="Usefulness versus diversity" icon={ChartColumn}>
                  {o.frontier_plot ? <ImageGallery single images={[{ src: o.frontier_plot, label: 'Usefulness versus diversity' }]} /> : <div className="empty-state">No chart.</div>}
                </Card>
                <Card title="Every penalty degree" icon={Table2} delay={0.06}>
                  <DataTable rows={o.frontier} maxHeight={420}
                    columns={['Penalty_Degree', 'Expected_Gain', 'Group_Coancestry', 'Progeny_Inbreeding', 'N_Parents']} />
                </Card>
              </div>
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
