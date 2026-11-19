import { useRef, useState } from 'react'
import {
  ChartLine, SlidersHorizontal, LayoutGrid, FileSpreadsheet, Dna, Globe, Group, Palette, Table2, Activity, ChartColumn,
  Calculator, Waypoints, Braces, Scale, Download, FileArchive, ClipboardCheck, LayoutPanelTop, ListChecks, Settings2, Medal,
} from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, SelectField, Toggle, ChipGroup, RadioChips, FilePicker, SubTabs,
  EmptyState, PlotSlot, DataTable, Note,
} from '../components/kit.jsx'
import { useBridge, useTool, useSignal, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'data', label: 'Data & Design', icon: SlidersHorizontal },
  { id: 'analyses', label: 'Analyses', icon: LayoutGrid },
  { id: 'downloads', label: 'Downloads', icon: Download },
]

// One tab per analysis. multi: the analysis takes several traits (min = how many it needs); otherwise one trait.
// plots / tables: what the analysis writes (id -> label), in the order shown.
const MODS = [
  { id: 'descriptive', label: 'Descriptives', title: 'Descriptive statistics', icon: ChartColumn, multi: true, min: 1,
    plots: { env_boxplot: 'Boxplot per environment', env_violin_plot: 'Violin plot per environment', distribution_plot: 'Trait distribution' },
    tables: { descriptive_stats: 'Descriptive statistics of every trait (all environments, then each environment)' } },
  { id: 'environments', label: 'Environments', title: 'Environments', icon: Globe, multi: false,
    plots: { env_means_plot: 'Accession means per environment', env_rank_plot: 'Accession ranks across environments', env_cor_heatmap: 'Environment correlation' },
    tables: { env_means: 'Means and ranks of every accession per environment', env_cor: 'Correlation between environments' } },
  { id: 'blues', label: 'BLUEs', title: 'BLUEs', icon: Activity, multi: false,
    plots: { env_blues_plot: 'Within-environment BLUEs', blues_plot: 'Mixed-model BLUEs' },
    tables: { blues: 'BLUEs per environment and over all environments' } },
  { id: 'anova', label: 'ANOVA', title: 'ANOVA', icon: Calculator, multi: true, min: 1,
    plots: { trait_variance_plot: 'Where the variation comes from, per trait and season', residual_plot: 'Residuals vs fitted', qq_plot: 'Normal Q-Q' },
    tables: {} },
  { id: 'herit', label: 'Heritability', title: 'Heritability', icon: Dna, multi: true, min: 1,
    plots: { heritability_plot: 'Entry-mean heritability per trait' },
    tables: { heritability: 'Variance components and heritability per trait and season' } },
  { id: 'blups', label: 'BLUPs', title: 'BLUPs', icon: Medal, multi: true, min: 1,
    plots: { blup_plot: 'BLUPs of the accessions, best first' },
    tables: { blups: 'BLUP of every accession per trait and environment (Pooled = all environments together)',
              blup_summary: 'Variance components and heritability of the BLUP models' } },
  { id: 'correlation', label: 'Correlation', title: 'Correlation', icon: Waypoints, multi: true, min: 2,
    plots: { correlation_heatmap: 'Trait correlation', correlation_network: 'Correlation circles' },
    tables: { correlation_table: 'Correlation between traits (with P values)' } },
  { id: 'pca', label: 'PCA', title: 'PCA', icon: Braces, multi: true, min: 2,
    plots: { pca_plot: 'PCA biplot', pca_variance_plot: 'Variance explained' },
    tables: { pca: 'PCA: trait loadings, eigenvalues and variance explained' } },
  { id: 'cluster', label: 'Clustering', title: 'Clustering', icon: Group, multi: true, min: 2,
    plots: { cluster_plot: 'Clusters (PCA space)', dendrogram: 'Dendrogram' },
    tables: { clusters: 'Clusters: cluster means, then every accession with its principal component scores and trait means' } },
  { id: 'stability', label: 'Stability & G×E', title: 'Stability & G×E', icon: Scale, multi: false,
    plots: { stability_plot: 'Stability (Wricke ecovalence vs mean)', ammi_plot: 'G×E interaction (means)', ammi_biplot: 'AMMI1 biplot' },
    tables: { stability: 'Stability, AMMI scores and G×E effects of every accession (sorted by stability rank)' } },
]

// Settings an analysis has besides its traits
const EXTRA_DEFAULTS = {
  descriptive: { plot_trait: '', by_env: true },
  blups: { plot_trait: '', by_env: true },
  blues: { interaction: true },
  anova: { diag_trait: '' },
  correlation: { method: 'pearson' },
  pca: { scale: true },
  cluster: { k: 3, method: 'kmeans' },
}

const HERIT_INTRO = 'H² entry-mean: heritability of an accession mean (use this for selection). H² plot basis: of a single plot. ' +
  'Pooled = all environments together. In an augmented design the error comes from the check plots only, so single-environment ' +
  'values can swing a lot; 0 = no genotypic variance detected.'

const BLUP_INTRO = 'BLUP = best linear unbiased prediction of an accession: its deviation from the mean, shrunk towards 0 the less it was measured. Rank accessions on it. ' +
  'Predicted value = BLUP + the average of the fixed part (environment, check) over the plots of the accessions. SE and Reliability (1 − PEV / Vg, 0 to 1) come from the ' +
  'mixed-model equations, so the uncertainty of the fixed effects is included. H² Cullis = 1 − (mean variance of a difference between two BLUPs) / (2 Vg) (Cullis et al. 2006): ' +
  'the heritability of the BLUPs, the one to report for unbalanced trials; Mean reliability = 1 − mean(PEV) / Vg is a slightly different number.'

const COL_GENOTYPE = 'Name of the variety or line. Checks and new entries both go here; a genotype that appears more than once is treated as replicated.'
const COL_ENV = 'Location, season or year label. At least 2 different environments are needed.'
const DESIGNS = [
  { value: 'augmented', label: 'Augmented RCBD', required: ['Genotype', 'Environment', 'Rep', 'Block'],
    hint: 'Checks repeated in every block, new entries once. Blocks are fixed and their effects come from the checks; new entries are adjusted to their block. A block is Block, or Rep + Block when every such group holds a check.',
    info: { Genotype: COL_GENOTYPE, Environment: COL_ENV,
      Rep: 'Replicate number. Use 1 in every row if the trial is not replicated.',
      Block: 'Block number. Every block holds all the checks plus some new entries.' },
    example: { cols: ['Genotype', 'Environment', 'Rep', 'Block', 'Plot', 'Yield', 'PlantHeight'],
      rows: [['Check1', 'Env1', 1, 1, 101, 5.4, 98], ['New07', 'Env1', 1, 1, 102, 6.1, 105], ['Check2', 'Env1', 1, 1, 103, 5.3, 102], ['New21', 'Env1', 1, 2, 201, 5.8, 110]] } },
  { value: 'rcbd', label: 'RCBD (randomized complete block)', required: ['Genotype', 'Environment', 'Rep'],
    hint: 'Every genotype once in each replicate. Rep is the block; a Block column is not needed.',
    info: { Genotype: COL_GENOTYPE, Environment: COL_ENV,
      Rep: 'Replicate (complete block) number. Every genotype appears once in each Rep.' },
    example: { cols: ['Genotype', 'Environment', 'Rep', 'Plot', 'Yield', 'PlantHeight'],
      rows: [['Geno05', 'Env1', 1, 101, 6.8, 106], ['Geno03', 'Env1', 1, 102, 7.8, 112], ['Geno05', 'Env1', 2, 201, 6.5, 104], ['Geno03', 'Env1', 2, 202, 7.9, 113]] } },
  { value: 'alpha', label: 'Alpha-lattice', required: ['Genotype', 'Environment', 'Rep', 'Block'],
    hint: 'Incomplete blocks inside each replicate. Block is numbered inside its Rep, so block 1 of Rep 1 and block 1 of Rep 2 are different blocks.',
    info: { Genotype: COL_GENOTYPE, Environment: COL_ENV,
      Rep: 'Replicate number. Every genotype appears once in each Rep.',
      Block: 'Incomplete block number, counted inside its Rep (1, 2, 3 … again in every Rep).' },
    example: { cols: ['Genotype', 'Environment', 'Rep', 'Block', 'Plot', 'Yield', 'PlantHeight'],
      rows: [['Geno05', 'Env1', 1, 1, 101, 6.8, 106], ['Geno03', 'Env1', 1, 1, 102, 7.8, 112], ['Geno12', 'Env1', 1, 2, 105, 5.9, 99], ['Geno05', 'Env1', 2, 3, 201, 6.5, 104]] } },
]

const COLOR_SCHEMES = [
  { value: 'default', label: 'Default (Green/Red)' }, { value: 'blue', label: 'Blue Gradient' },
  { value: 'viridis', label: 'Viridis' }, { value: 'plasma', label: 'Plasma' }, { value: 'magma', label: 'Magma' },
]

/* ───────────── ANOVA table (publication style: rules, no vertical clutter) ───────────── */

const ANOVA_COL = /^(.*): (DF|SS|MS|F|P|Sig\.|F vs G×E|P vs G×E)$/
const ANOVA_FOOT = ['Plots', 'Mean', 'CV (%)']
const ANOVA_HEAD = { DF: 'DF', SS: 'SS', MS: 'MS', F: 'F', P: 'P', 'Sig.': '', 'F vs G×E': 'F (vs G×E)', 'P vs G×E': 'P (vs G×E)' }
const has = (v) => v !== null && v !== undefined

function anovaFmt(x, lab) {
  if (!has(x) || x === '') return ''
  if (lab === 'DF') return String(x)
  if (lab === 'SS' || lab === 'MS') return Number(x).toLocaleString('en-US', { minimumFractionDigits: 2, maximumFractionDigits: 2 })
  if (lab === 'F' || lab === 'F vs G×E') return Number(x).toFixed(2)
  if (lab === 'P' || lab === 'P vs G×E') return x < 0.001 ? '<0.001' : Number(x).toFixed(3)
  return String(x)
}

function AnovaOne({ rows }) {
  const cols = Object.keys(rows[0]).filter((k) => k !== 'Trait' && k !== 'Source')
    .map((key) => { const m = ANOVA_COL.exec(key); return m ? { key, grp: m[1], lab: m[2] } : null })
    .filter((c) => c && rows.some((r) => has(r[c.key])))
  const groups = [...new Set(cols.map((c) => c.grp))]
  const body = rows.filter((r) => !ANOVA_FOOT.includes(r.Source))
  const foot = rows.filter((r) => ANOVA_FOOT.includes(r.Source))
  const first = (i) => cols.findIndex((c) => c.grp === cols[i].grp) === i
  const footText = (src, x) => (!has(x) ? '–' : src === 'Plots' ? String(x) : src === 'Mean' ? String(Number(Number(x).toPrecision(5))) : Number(x).toFixed(2))
  const gxe = cols.some((c) => c.lab === 'F vs G×E')
  return (
    <>
      <div className="table-wrap">
        <table className="anova-tbl">
          <thead>
            <tr>
              <th rowSpan={2} className="src">Source</th>
              {groups.map((g) => <th key={g} colSpan={cols.filter((c) => c.grp === g).length} className="grp center">{g === 'Pooled' ? 'Pooled over environments' : g}</th>)}
            </tr>
            <tr>{cols.map((c, i) => <th key={c.key} className={`end ${first(i) ? 'grp' : ''}`}>{ANOVA_HEAD[c.lab]}</th>)}</tr>
          </thead>
          <tbody>
            {body.map((r, ri) => {
              const sub = String(r.Source).startsWith(' ')
              return (
                <tr key={ri} className={r.Source === 'Total' ? 'total' : sub ? 'muted' : ''}>
                  <td className={sub ? 'sub' : ''}>{String(r.Source).trim()}</td>
                  {cols.map((c, i) => {
                    const x = r[c.key]
                    if (c.lab === 'Sig.') {
                      return <td key={c.key} className={`sig ${x === 'ns' ? 'ns' : x ? 'yes' : ''}`}>{has(x) ? x : ''}</td>
                    }
                    return <td key={c.key} className={`end ${first(i) ? 'grp' : ''}`}>{anovaFmt(x, c.lab)}</td>
                  })}
                </tr>
              )
            })}
            {foot.map((r, ri) => (
              <tr key={'f' + ri} className={`foot ${ri === 0 ? 'first' : ''}`}>
                <td>{r.Source}</td>
                {groups.map((g) => <td key={g} colSpan={cols.filter((c) => c.grp === g).length} className="grp center">{footText(r.Source, r[`${g}: DF`])}</td>)}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <div className="field-hint" style={{ marginTop: 8 }}>
        DF degrees of freedom; SS sum of squares; MS mean square; CV = √MS(error) / mean × 100. *** P&lt;0.001, ** P&lt;0.01, * P&lt;0.05, ns not significant.
        {gxe && ' “vs G×E”: genotype tested against the genotype × environment mean square (environments as a sample).'}
      </div>
    </>
  )
}

function AnovaTables({ rows, first }) {
  const traits = [...new Set(rows.map((r) => r.Trait))]
  const ordered = first && traits.includes(first) ? [first, ...traits.filter((t) => t !== first)] : traits
  return (
    <div className="anova-list">
      {ordered.map((tr, i) => (
        <details key={tr} open={i === 0} className="anova-item">
          <summary>{tr}</summary>
          <AnovaOne rows={rows.filter((r) => r.Trait === tr)} />
        </details>
      ))}
    </div>
  )
}

/* ───────────── Result pieces ───────────── */

const cell = (c, v) => {
  if (c === 'Sig.' && v) return <b style={{ color: v === 'ns' ? '#9ca3af' : '#15803d' }}>{v}</b>
  if (typeof v === 'number') return String(v)
  return undefined
}

function TableCard({ id, title, rows, path, onSave }) {
  return (
    <div className="plot-block">
      <div className="plot-title">
        <span>{title}</span>
        {path && <button className="btn btn-secondary sm" style={{ marginLeft: 'auto' }} onClick={() => onSave(path, `${id}.csv`)}><Download size={14} /> CSV</button>}
      </div>
      <DataTable rows={rows} renderCell={cell} maxHeight={340} />
    </div>
  )
}

function Messages({ msgs }) {
  const of = (type) => (msgs || []).filter((m) => m.type === type).map((m) => m.text)
  const box = (type, cls, title) => {
    const t = of(type)
    return t.length ? (
      <div className={`info-line ${cls}`}>
        <b>{title}</b>
        <ul style={{ margin: '4px 0 0 18px', padding: 0 }}>{t.map((x, i) => <li key={i}>{x}</li>)}</ul>
      </div>
    ) : null
  }
  return <>{box('error', 'warn', 'This analysis could not be completed')}{box('warning', 'warn', 'Warnings')}{box('note', '', 'Notes')}</>
}

function ModuleResults({ mod, r, anovaFirst, onSave }) {
  const msgs = (r.msgs || []).filter((m) => m.module === mod.id)
  const plots = Object.entries(mod.plots).filter(([id]) => r.plots && r.plots[id])
  const tables = Object.entries(mod.tables).filter(([id]) => r.tables_data && r.tables_data[id])
  const wide = mod.id === 'anova' && r.anova_wide_data && r.anova_wide_data.length
  return (
    <>
      <Messages msgs={msgs} />
      {mod.id === 'herit' && <div className="field-hint" style={{ marginTop: 10 }}>{HERIT_INTRO}</div>}
      {mod.id === 'blups' && <div className="field-hint" style={{ marginTop: 10 }}>{BLUP_INTRO}</div>}
      {wide ? (
        <div className="plot-block" style={{ marginTop: 14 }}>
          <div className="plot-title">
            <span>ANOVA (one table per trait; seasons and pooled analysis side by side)</span>
            {r.tables && r.tables.anova && <button className="btn btn-secondary sm" style={{ marginLeft: 'auto' }} onClick={() => onSave(r.tables.anova, 'anova.csv')}><Download size={14} /> CSV</button>}
          </div>
          <AnovaTables rows={r.anova_wide_data} first={anovaFirst} />
        </div>
      ) : null}
      {plots.length > 0 && (
        <div className="plot-grid-2" style={{ marginTop: 14 }}>
          {plots.map(([id, label]) => <PlotSlot key={id} name={id} title={label} plot={r.plots[id]} />)}
        </div>
      )}
      {tables.map(([id, label]) => (
        <div key={id} style={{ marginTop: 14 }}>
          <TableCard id={id} title={label} rows={r.tables_data[id]} path={r.tables && r.tables[id]} onSave={onSave} />
        </div>
      ))}
      {!wide && plots.length === 0 && tables.length === 0 && (
        <div className="empty-state" style={{ padding: 24 }}>Nothing was produced. See the messages above.</div>
      )}
    </>
  )
}

/* ───────────── The page ───────────── */

const initCfg = (traits, prev = {}, fresh = true) => {
  const out = {}
  MODS.forEach((m) => {
    const old = prev[m.id] || {}
    const base = { ...(EXTRA_DEFAULTS[m.id] || {}), ...old }
    if (m.multi) {
      let t = fresh ? traits : (old.traits || []).filter((x) => traits.includes(x))
      if (!t.length) t = traits
      out[m.id] = { ...base, traits: t }
    } else {
      out[m.id] = { ...base, trait: traits.includes(old.trait) ? old.trait : traits[0] || '' }
    }
  })
  return out
}

export default function PhenotypicTool() {
  const { bridge } = useBridge('phenotypicBridge', () => mockGeneric(['dataLoaded', 'exportReady']))
  const tool = useTool(bridge, 'run_module')
  const [tab, setTab] = useState('data')
  const [sub, setSub] = useState('descriptive')
  const [file, setFile] = useState('')
  const [info, setInfo] = useState(null)
  const [loading, setLoading] = useState(false)
  const [design, setDesign] = useState('augmented')
  const [mode, setMode] = useState('all')
  const [envs, setEnvs] = useState([])
  const [scheme, setScheme] = useState('default')
  const [zero, setZero] = useState(true)
  const [cfg, setCfg] = useState(initCfg([]))
  const [results, setResults] = useState({})
  const [outdir, setOutdir] = useState('')
  const [resSig, setResSig] = useState('')
  const [exporting, setExporting] = useState('')
  const req = useRef(0)
  const fresh = useRef(true)
  const timer = useRef(null)
  const pending = useRef(null)

  const sigOf = (o) => JSON.stringify([o.file, o.design, o.mode, o.mode === 'specific' ? [...o.envs].sort() : [], o.scheme, o.zero])
  const curSig = sigOf({ file, design, mode, envs, scheme, zero })

  // ── Reading the data file (again whenever a data setting that changes what is analysed changes) ──
  const inspect = (o, delay = 0) => {
    clearTimeout(timer.current)
    const go = () => {
      if (!bridge || !o.file) return
      const id = ++req.current
      setLoading(true)
      bridge.inspect_data(JSON.stringify({
        req: id, input_file: o.file, design: o.design, analysis_mode: o.mode,
        selected_environments: o.mode === 'specific' ? o.envs : [], zero_as_na: o.zero,
      }))
    }
    if (delay) timer.current = setTimeout(go, delay); else go()
  }
  const now = (over = {}) => ({ file, design, mode, envs, zero, ...over })

  useSignal(bridge, 'dataLoaded', (json) => {
    let pv
    try { pv = JSON.parse(json) } catch { return }
    if (pv.req !== undefined && pv.req !== req.current) return          // an older request
    setLoading(false)
    if (pv.status !== 'success') { setInfo(null); return tool.setStatus('Error reading file: ' + pv.message) }
    const L = (v) => (Array.isArray(v) ? v : v === undefined || v === null || v === '' ? [] : [v])   // a one-item list must still be a list
    for (const k of ['columns', 'environments', 'traits', 'empty_traits', 'errors', 'warnings', 'notes', 'design_notes']) pv[k] = L(pv[k])
    setInfo(pv)
    const isNew = fresh.current
    if (isNew) {
      fresh.current = false
      setDesign(pv.design.code)
      setEnvs(pv.environments)
      setMode('all')
    }
    setCfg((prev) => initCfg(pv.traits, prev, isNew))
  })

  const onBrowse = async () => {
    const p = await tool.browse()
    if (!p) return
    setFile(p)
    fresh.current = true
    setInfo(null)
    inspect({ file: p, design: 'auto', mode: 'all', envs: [], zero })
  }
  const onDesign = (v) => { setDesign(v); inspect(now({ design: v })) }
  const onMode = (v) => { setMode(v); inspect(now({ mode: v })) }
  const onEnvs = (v) => { setEnvs(v); if (mode === 'specific') inspect(now({ envs: v }), 700) }
  const onZero = (v) => { setZero(v); inspect(now({ zero: v })) }

  // ── Running one analysis ──
  const mod = MODS.find((m) => m.id === sub)
  const c = cfg[sub] || {}
  const traitsOf = (m) => (m.multi ? (cfg[m.id].traits || []) : [])

  tool.onResult((r) => {
    if (!r.module) return
    const p = pending.current
    setOutdir(r.outdir)
    setResSig(p ? p.sig : curSig)
    setResults((prev) => ({ ...prev, [r.module]: r }))
  })

  const run = () => {
    if (!file || !info) return tool.fail('Select a data file first.')
    if (info.errors && info.errors.length) return tool.fail('The data cannot be analysed: ' + info.errors.join('; '))
    if (mode === 'specific' && !envs.length) return tool.fail('Choose at least one environment.')
    if (mod.multi && (c.traits || []).length < mod.min) {
      return tool.fail(`Select at least ${mod.min} trait${mod.min > 1 ? 's' : ''} for ${mod.title}.`)
    }
    if (!mod.multi && !c.trait) return tool.fail('Select a trait.')

    const settings = {}
    if (mod.multi) settings.traits = c.traits; else settings.trait = c.trait
    if (mod.id === 'descriptive') { settings.plot_trait = c.traits.includes(c.plot_trait) ? c.plot_trait : c.traits[0]; settings.by_env = !!c.by_env }
    if (mod.id === 'blups') { settings.plot_trait = c.traits.includes(c.plot_trait) ? c.plot_trait : c.traits[0]; settings.by_env = !!c.by_env }
    if (mod.id === 'blues') settings.interaction = !!c.interaction
    if (mod.id === 'anova') settings.diag_trait = c.traits.includes(c.diag_trait) ? c.diag_trait : c.traits[0]
    if (mod.id === 'correlation') settings.method = c.method
    if (mod.id === 'pca') settings.scale = !!c.scale
    if (mod.id === 'cluster') { settings.k = c.k; settings.method = c.method }

    // new data settings: the earlier results no longer match and are started afresh
    let od = outdir
    if (od && resSig !== curSig) od = ''
    if (!od) { setResults({}); setOutdir('') }
    pending.current = { sig: curSig }
    tool.run({
      input_file: file, design, analysis_mode: mode, selected_environments: mode === 'specific' ? envs : [],
      color_scheme: scheme, zero_as_na: zero, module: sub, settings, outdir: od,
    })
  }

  // ── Downloads ──
  useSignal(bridge, 'exportReady', (json) => {
    let r
    try { r = JSON.parse(json) } catch { return }
    setExporting('')
    if (r.status !== 'success') return tool.fail(r.message || 'Could not prepare the file.')
    tool.save(r.path, r.name)
  })
  const exportFile = (kind) => {
    if (!bridge || !outdir) return
    setExporting(kind)
    tool.setStatus('Preparing the file…')
    bridge.export_files(JSON.stringify({ outdir, kind }))
  }

  const dz = DESIGNS.find((d) => d.value === design) || DESIGNS[0]
  const ran = MODS.filter((m) => results[m.id])
  const stale = ran.length > 0 && resSig !== curSig
  const dc = (resSig === curSig && results[ran.length ? ran[ran.length - 1].id : ''] && results[ran[ran.length - 1].id].tables_data.design_check) || (info && info.design_check)
  const setC = (patch) => setCfg((prev) => ({ ...prev, [sub]: { ...prev[sub], ...patch } }))
  const allTraits = (info && info.traits) || []
  const traitChoices = mod.multi ? (c.traits || []) : []

  return (
    <ToolPage
      id="phen" icon={ChartLine} tabs={TABS} tab={tab} setTab={setTab}
      title="Phenotypic Data Analysis"
      subtitle="Descriptive stats, environments, BLUEs, ANOVA, heritability, BLUPs, correlation, PCA, clustering and G×E stability"
      run={tab === 'analyses' ? {
        label: `Run ${mod.title}`, runningLabel: 'Analysing…', state: tool.state, status: tool.status, onRun: run, disabled: !info || loading,
      } : null}
    >
      {tab === 'data' && (
        <TabPane id="data" key="data">
          <SectionHeader title="Data & Design">
            Load the file once; every analysis tab uses it. Columns are found by their name, in any order: {dz.required.join(', ')} for this design. Every other numeric column is a trait.
          </SectionHeader>
          <Card title="Phenotypic Data File" icon={FileSpreadsheet}>
            <FilePicker value={file} placeholder="Select .xlsx, .xls or .csv file" onBrowse={onBrowse} />
            {loading && <div className="info-line">Reading the data file…</div>}
            {info && (
              <div className="summary-grid" style={{ marginTop: 16, marginBottom: 0 }}>
                {[['Design', info.design.label], ['Plots', info.n_plots], ['Accessions', info.n_genotypes], ['Environments', info.n_environments], ['Traits', info.traits.length]].map(([l, v]) => (
                  <div key={l} className="summary-card"><div className="val" style={l === 'Design' ? { fontSize: '1.05rem' } : undefined}>{v}</div><div className="lbl">{l}</div></div>
                ))}
              </div>
            )}
            {info && info.errors.length > 0 && <div className="info-line warn"><b>The data cannot be analysed: </b>{info.errors.join('; ')}. Columns found: {info.columns.join(', ')}.</div>}
            {info && info.warnings.map((w, i) => <div key={'w' + i} className="info-line warn">{w}</div>)}
            {info && info.empty_traits.length > 0 && <div className="info-line warn">Left out (no data{zero ? ', zeros count as missing' : ''}): {info.empty_traits.join(', ')}</div>}
            {info && info.notes.map((w, i) => <div key={'n' + i} className="info-line">{w}</div>)}
            {info && info.design_notes.map((w, i) => <div key={'d' + i} className="info-line warn">{w}</div>)}
          </Card>
          <div className="row">
            <Card title="Experimental Design" icon={LayoutPanelTop} delay={0.03}>
              <SelectField label="Design of the trial" value={design} onChange={onDesign} options={DESIGNS.map(({ value, label }) => ({ value, label }))} />
              <Note>{dz.hint}</Note>
            </Card>
            <Card title="Environments" icon={Globe} delay={0.05}>
              <RadioChips value={mode} onChange={onMode} options={[{ value: 'all', label: 'All environments' }, { value: 'specific', label: 'Choose environments' }]} />
              {mode === 'specific' && info && (
                <div className="field" style={{ marginTop: 12 }}>
                  <ChipGroup options={info.environments} value={envs} onChange={onEnvs} />
                </div>
              )}
              <div style={{ marginTop: 12 }}>
                <Toggle label="Treat 0 as missing (not measured)" checked={zero} onChange={onZero} />
              </div>
            </Card>
            <Card title="Color Scheme" icon={Palette} delay={0.08}>
              <SelectField value={scheme} onChange={setScheme} options={COLOR_SCHEMES} />
              <div className="field-hint">Used by the figures of the analyses run from now on.</div>
            </Card>
          </div>
          {info && info.errors.length === 0 && (
            <>
              <Card title="Traits" icon={ListChecks} delay={0.1}>
                <DataTable rows={info.trait_stats} renderCell={(col, v) => (typeof v === 'number' ? String(v) : undefined)} maxHeight={260} />
              </Card>
              <Card title="Design Check" icon={ClipboardCheck} delay={0.12}>
                <DataTable rows={dc || []} renderCell={(col, v) => (typeof v === 'number' ? String(v) : undefined)} maxHeight={260} />
              </Card>
            </>
          )}
          <Card title="Input File Format" icon={Table2} delay={0.14}>
            <div className="table-wrap">
              <table className="results-table">
                <thead><tr><th>Column name</th><th>Needed</th><th>What it holds</th></tr></thead>
                <tbody>
                  {dz.required.map((col) => <tr key={col}><td><b>{col}</b></td><td>Yes</td><td>{dz.info[col]}</td></tr>)}
                  <tr><td><b>Any other numeric column</b>, e.g. Yield, PlantHeight</td><td>Traits</td>
                    <td>Leave a cell empty when a trait was not measured. A 0 is also read as “not measured” (switch this off under Environments).</td></tr>
                  <tr><td>Plot, Row, Column, Range, Year, Season, Entry, Line.type, Remarks</td><td>No</td>
                    <td>Layout and label columns. They are not treated as traits.</td></tr>
                </tbody>
              </table>
            </div>
            <div className="field-hint" style={{ marginTop: 12 }}>One row per plot. Example for the {dz.label} design:</div>
            <div className="table-wrap">
              <table className="results-table">
                <thead><tr>{dz.example.cols.map((col) => <th key={col}>{col}</th>)}</tr></thead>
                <tbody>{dz.example.rows.map((row, i) => <tr key={i}>{row.map((v, j) => <td key={j}>{v}</td>)}</tr>)}</tbody>
              </table>
            </div>
          </Card>
        </TabPane>
      )}

      {tab === 'analyses' && (
        <TabPane id="analyses" key="analyses">
          <SectionHeader title="Analyses">
            Each analysis has its own traits and settings and its own Run button; running one does not run the others. Results stay until the data settings change.
          </SectionHeader>
          {!info ? <EmptyState icon={ChartLine}>Load a data file in Data &amp; Design first.</EmptyState> : (
            <>
              <SubTabs tabs={MODS.map(({ id, label, icon }) => ({ id, label, icon }))} value={sub} onChange={setSub} id="phsub" />
              {stale && <div className="info-line warn">The data settings changed after these results were produced. Running an analysis starts afresh and clears them.</div>}
              <div className="row" style={{ marginTop: 14 }}>
                <Card title={mod.multi ? 'Traits to analyse' : 'Trait'} icon={Settings2}>
                  {mod.multi ? (
                    <>
                      <ChipGroup options={allTraits} value={c.traits || []} onChange={(v) => setC({ traits: v })} />
                      <div style={{ marginTop: 10, display: 'flex', gap: 8 }}>
                        <button className="btn btn-secondary sm" onClick={() => setC({ traits: allTraits })}>Select all</button>
                        <button className="btn btn-secondary sm" onClick={() => setC({ traits: [] })}>Clear</button>
                      </div>
                      <div className="field-hint">{mod.min > 1 ? `At least ${mod.min} traits are needed.` : 'Every ticked trait is analysed.'}</div>
                    </>
                  ) : (
                    <SelectField value={c.trait || ''} onChange={(v) => setC({ trait: v })} options={allTraits} />
                  )}
                </Card>
                {['descriptive', 'blues', 'blups', 'anova', 'correlation', 'pca', 'cluster'].includes(mod.id) && (
                  <Card title="Options" icon={Settings2} delay={0.04}>
                    {mod.id === 'descriptive' && (
                      <>
                        <SelectField label="Trait for the distribution plots" value={traitChoices.includes(c.plot_trait) ? c.plot_trait : traitChoices[0] || ''} onChange={(v) => setC({ plot_trait: v })} options={traitChoices} />
                        <Toggle label="Also show each environment separately" checked={!!c.by_env} onChange={(v) => setC({ by_env: v })} />
                      </>
                    )}
                    {mod.id === 'blups' && (
                      <>
                        <SelectField label="Trait for the figure" value={traitChoices.includes(c.plot_trait) ? c.plot_trait : traitChoices[0] || ''} onChange={(v) => setC({ plot_trait: v })} options={traitChoices} />
                        <Toggle label="Also compute the BLUPs of each environment separately" checked={!!c.by_env} onChange={(v) => setC({ by_env: v })} />
                      </>
                    )}
                    {mod.id === 'blues' && <Toggle label="Include Genotype × Environment (BLUEs per environment)" checked={!!c.interaction} onChange={(v) => setC({ interaction: v })} />}
                    {mod.id === 'anova' && <SelectField label="Trait for the residual plots" value={traitChoices.includes(c.diag_trait) ? c.diag_trait : traitChoices[0] || ''} onChange={(v) => setC({ diag_trait: v })} options={traitChoices} />}
                    {mod.id === 'correlation' && (
                      <SelectField label="Correlation method" value={c.method} onChange={(v) => setC({ method: v })}
                        options={[{ value: 'pearson', label: 'Pearson' }, { value: 'spearman', label: 'Spearman (rank)' }, { value: 'kendall', label: 'Kendall' }]} />
                    )}
                    {mod.id === 'pca' && <Toggle label="Standardise the traits first (correlation-matrix PCA)" checked={!!c.scale} onChange={(v) => setC({ scale: v })} />}
                    {mod.id === 'cluster' && (
                      <>
                        <Slider label="Number of clusters" value={c.k} onChange={(v) => setC({ k: v })} min={2} max={8} />
                        <RadioChips value={c.method} onChange={(v) => setC({ method: v })} options={[{ value: 'kmeans', label: 'K-means' }, { value: 'hierarchical', label: 'Hierarchical (Ward)' }]} />
                      </>
                    )}
                  </Card>
                )}
              </div>
              <div style={{ marginTop: 18 }}>
                {results[sub]
                  ? <ModuleResults mod={mod} r={results[sub]} anovaFirst={c.diag_trait} onSave={tool.save} />
                  : <div className="empty-state" style={{ padding: 24 }}>Not run yet. Choose the traits and settings above, then press Run {mod.title}.</div>}
              </div>
            </>
          )}
        </TabPane>
      )}

      {tab === 'downloads' && (
        <TabPane id="downloads" key="downloads">
          <SectionHeader title="Downloads">Everything that has been run so far in this session.</SectionHeader>
          {!outdir ? <EmptyState icon={ChartLine}>Run an analysis first.</EmptyState> : (
            <Card title="Files" icon={FileArchive}>
              <div className="downloads">
                {[
                  { kind: 'all', icon: FileSpreadsheet, label: 'All tables (.xlsx)' },
                  results.anova && { kind: 'anova', icon: FileSpreadsheet, label: 'ANOVA tables (.xlsx)' },
                  results.blups && { kind: 'blups', icon: FileSpreadsheet, label: 'BLUPs (.xlsx)' },
                  { kind: 'figures', icon: FileArchive, label: 'All figures: 300-dpi PNG + vector PDF (.zip)' },
                  { kind: 'everything', icon: FileArchive, label: 'Everything (.zip)' },
                ].filter(Boolean).map((d) => (
                  <button key={d.kind} className="btn btn-secondary" disabled={!!exporting} onClick={() => exportFile(d.kind)}>
                    <d.icon size={16} /> {exporting === d.kind ? 'Preparing…' : d.label}
                  </button>
                ))}
              </div>
              <Note>Analyses run so far: {ran.map((m) => m.title).join(', ')}.</Note>
            </Card>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
