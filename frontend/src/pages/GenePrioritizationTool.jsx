import { useEffect, useMemo, useState } from 'react'
import {
  ListOrdered, Dna, FileSpreadsheet, ListChecks, Info, Download, Search, Wifi, WifiOff, X, ExternalLink,
  ChartColumn, Settings2, ArrowUp, ArrowDown, Ban, Network,
} from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, FilePicker, TextField, SelectField, EmptyState, SummaryGrid, ProgressBar, Toggle,
} from '../components/kit.jsx'
import { ScoreBars, ContributionBars, PathwayMatrix, Legend, CAT_COLORS } from '../components/geneCharts.jsx'
import { useBridge, useTool, useSignal, callJson, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'input', label: 'Input', icon: ListChecks },
  { id: 'summary', label: 'Summary', icon: Info },
  { id: 'ranking', label: 'Ranked Genes', icon: ListOrdered },
  { id: 'evidence', label: 'Gene Evidence', icon: Dna },
  { id: 'charts', label: 'Charts', icon: ChartColumn },
]

// Used when the Python side is not there (browser preview) and as a fallback if the list cannot be read.
const SPECIES_FALLBACK = [
  { value: 'oryza_sativa', label: 'Oryza sativa (Rice)', enabled: true, id_help: 'RAP-DB / Ensembl Plants IDs such as Os02g0766700 (MSU LOC_Os IDs are translated for genes listed in Oryzabase)', id_example: 'Os02g0766700' },
  { value: 'triticum_aestivum', label: 'Triticum aestivum (Wheat)', enabled: true, id_help: 'Ensembl Plants IDs such as TraesCS1A02G000100', id_example: 'TraesCS1A02G000100' },
  { value: 'zea_mays', label: 'Zea mays (Maize)', enabled: true, id_help: 'Ensembl Plants IDs such as Zm00001eb000010', id_example: 'Zm00001eb000010' },
]

const CAT_ORDER = ['association', 'function', 'pathway', 'ortholog', 'literature']
const STATE_TEXT = { ok: 'connected', partial: 'partly answering', unavailable: 'not reachable', 'not used': 'not used yet' }

const pre = { margin: '8px 0', padding: '10px 14px', borderRadius: 8, background: 'rgba(0,0,0,0.06)', fontSize: 13, overflowX: 'auto' }

const f1 = (v) => (v === null || v === undefined ? '—' : Number(v).toFixed(1))
const dist = (v) => (v === null || v === undefined ? '' : `${Number(v).toLocaleString()} bp`)

function mockPrioritization() {
  const m = mockGeneric(['progressUpdate', 'apiStatus'])
  return new Proxy(m, {
    get(t, name) {
      if (name === 'get_species') return (cb) => setTimeout(() => cb(JSON.stringify(SPECIES_FALLBACK.map((s) => ({ ...s, reason: '' })))), 100)
      if (name === 'open_link') return (u) => window.open(u, '_blank')
      if (name === 'run_prioritization') {
        return () => {
          fetch('/gp_demo_result.json').then((r) => r.json())
            .then((r) => setTimeout(() => t.resultReady.emit(JSON.stringify(r)), 600))
            .catch(() => t.resultReady.emit(JSON.stringify({ status: 'error', message: 'Preview mode: there is no Python backend, so nothing can be retrieved.' })))
        }
      }
      return t[name]
    },
  })
}

function Ext({ href, bridge, children }) {
  if (!href) return <>{children}</>
  return (
    <a href={href} className="ext-link" onClick={(e) => { e.preventDefault(); bridge && bridge.open_link ? bridge.open_link(href) : window.open(href, '_blank') }}>
      {children}<ExternalLink size={11} style={{ marginLeft: 3, verticalAlign: -1 }} />
    </a>
  )
}

function ScoreCell({ v, color }) {
  if (v === null || v === undefined) return <span style={{ color: 'var(--text-muted)' }}>—</span>
  return (
    <span style={{ display: 'inline-flex', alignItems: 'center', gap: 8 }}>
      <span style={{ width: 70, height: 8, borderRadius: 4, background: 'rgba(82,183,136,.15)', display: 'inline-block', overflow: 'hidden' }}>
        <i style={{ display: 'block', height: '100%', width: `${Math.max(0, Math.min(100, v))}%`, background: color }} />
      </span>
      {f1(v)}
    </span>
  )
}

function ApiChips({ items }) {
  if (!items || !items.length) return null
  return (
    <div style={{ display: 'flex', flexWrap: 'wrap', gap: 10 }}>
      {items.map((a) => {
        const good = a.state === 'ok'
        const bad = a.state === 'unavailable'
        return (
          <span key={a.source} className="chip checked" title={a.error || ''}
            style={{ cursor: 'default', background: good ? undefined : bad ? 'var(--danger)' : 'var(--text-muted)' }}>
            {bad ? <WifiOff size={14} /> : <Wifi size={14} />} {a.source}: {STATE_TEXT[a.state] || a.state}
            {a.calls || a.cached ? ` · ${a.calls} request${a.calls === 1 ? '' : 's'}${a.cached ? `, ${a.cached} cached` : ''}` : ''}
          </span>
        )
      })}
    </div>
  )
}

export default function GenePrioritizationTool() {
  const { bridge } = useBridge('genePrioritizationBridge', () => mockPrioritization())
  const tool = useTool(bridge, 'run_prioritization', { accept: (r) => r.status === 'cancelled' })
  const [tab, setTab] = useState('input')
  const [species, setSpecies] = useState(SPECIES_FALLBACK)
  const [f, setF] = useState({
    species: 'oryza_sativa', trait: '', gene_text: '', gene_file: '', refresh_cache: false,
  })
  const [prog, setProg] = useState({ pct: 0, message: '', api: null })
  const [api, setApi] = useState(null)
  const [selected, setSelected] = useState(null)
  const [sort, setSort] = useState({ key: 'Rank', dir: 1 })
  const [q, setQ] = useState('')
  const [shown, setShown] = useState(100)

  useEffect(() => {
    if (!bridge) return
    callJson(bridge, 'get_species').then((s) => { if (Array.isArray(s) && s.length) setSpecies(s) })
  }, [bridge])
  useSignal(bridge, 'progressUpdate', (j) => { try { const p = JSON.parse(j); setProg((o) => ({ ...o, ...p, api: p.api || o.api })) } catch { /* ignore */ } })
  useSignal(bridge, 'apiStatus', (j) => { try { setApi(JSON.parse(j)) } catch { /* ignore */ } })

  tool.onResult((r) => {
    if (r.status === 'cancelled') { tool.setState('idle'); tool.setStatus('Cancelled.'); return }
    if (r.api_status) setApi(r.api_status)
    if (r.status === 'success') {
      setSelected(r.genes.length ? r.genes[0].Gene_ID : null)
      setSort({ key: 'Rank', dir: 1 }); setQ(''); setShown(100)
      setTab('summary')
    }
  })

  const set = (key) => (v) => setF((x) => ({ ...x, [key]: v }))
  const sp = species.find((s) => s.value === f.species) || species[0]
  const r = tool.result && tool.result.status === 'success' ? tool.result : null
  const running = tool.state === 'running'

  const browseGenes = async () => { const p = await tool.browse('browse_gene_list'); if (p) setF((x) => ({ ...x, gene_file: p })) }
  const run = () => {
    if (!sp || !sp.enabled) return tool.fail((sp && sp.reason) || 'Choose a supported species.')
    if (f.trait.trim().length < 3) return tool.fail('Enter the trait of interest, for example Drought tolerance.')
    if (!f.gene_text.trim() && !f.gene_file) return tool.fail('Provide the genes to rank: paste gene IDs or choose a gene ID list file.')
    setProg({ pct: 0, message: 'Starting…', api: null })
    setApi(null)
    tool.run({
      species: f.species, trait: f.trait.trim(), gene_text: f.gene_text, gene_file: f.gene_file, refresh_cache: f.refresh_cache,
    })
  }

  const rows = useMemo(() => {
    if (!r) return []
    const n = q.trim().toLowerCase()
    let out = r.genes
    if (n) out = out.filter((g) => [g.Gene_ID, g.Symbol, g.Description, g.Evidence_Level].some((x) => String(x || '').toLowerCase().includes(n)))
    const { key, dir } = sort
    return [...out].sort((a, b) => {
      const x = a[key], y = b[key]
      if (x === null || x === undefined || x === '') return 1
      if (y === null || y === undefined || y === '') return -1
      return (typeof x === 'number' && typeof y === 'number' ? x - y : String(x).localeCompare(String(y))) * dir
    })
  }, [r, q, sort])
  const sortBy = (key) => setSort((s) => (s.key === key ? { key, dir: -s.dir } : { key, dir: key === 'Rank' || typeof r.genes[0][key] === 'string' ? 1 : -1 }))
  const openGene = (id) => { setSelected(id); setTab('evidence') }
  const download = (t) => bridge && bridge.download_result(t)

  const cols = r ? [
    { key: 'Rank', label: 'Rank' }, { key: 'Gene_ID', label: 'Gene ID' }, { key: 'Symbol', label: 'Gene Symbol' },
    { key: 'Description', label: 'Description' }, { key: 'Final_Score', label: 'Final Score' },
    { key: 'Evidence_Level', label: 'Evidence Level' },
  ] : []
  const gd = r && selected ? r.detail[selected] : null

  return (
    <ToolPage
      id="geneprio" icon={ListOrdered} tabs={TABS} tab={tab} setTab={setTab}
      title="Gene Prioritization"
      subtitle="Rank candidate genes for a trait from Ensembl Plants, Gramene, Plant Reactome, PubMed and (rice) Oryzabase / RAP-DB"
      run={tab === 'input' ? { label: 'Prioritize Genes', runningLabel: 'Retrieving evidence…', state: tool.state, status: tool.status, onRun: run } : undefined}
    >
      {/* ========================================================= INPUT */}
      {tab === 'input' && (
        <TabPane id="input" key="input">
          <SectionHeader title="Instructions">Choose the species and the trait, then give the genes to rank.</SectionHeader>
          <div className="note" style={{ marginTop: 0, marginBottom: 20 }}>
            <p>1. Choose the <b>species</b> and type the <b>trait</b> (for example “Drought tolerance”, “Salt tolerance”, “Grain yield”). Related words are added automatically (drought: water deficit, osmotic stress …).</p>
            <p>2. Give the <b>candidate genes</b>: paste gene IDs (one per line, or separated by commas, spaces or tabs) or choose a TXT, CSV or Excel file whose first column holds the IDs — for example the Gene_ID_list.txt written by Genes Extraction for the QTL. Gene symbols and UniProt IDs are recognised too.</p>
            <p>3. Click “Prioritize Genes”. Every number comes from the databases queried: when a service cannot be reached the evidence is shown as unavailable, never estimated. The score ranks genes by the evidence found, it is not a probability that a gene is causal.</p>
          </div>

          {running && (
            <Card title="Progress" icon={Network}>
              <ProgressBar value={(prog.pct || 0) / 100} />
              <div style={{ margin: '4px 0 12px' }}>{prog.message}</div>
              <ApiChips items={prog.api} />
              <div style={{ marginTop: 14 }}><button className="btn btn-secondary" onClick={() => bridge && bridge.cancel()}><Ban size={15} /> Cancel</button></div>
            </Card>
          )}

          <div className="row">
            <Card title="Species and trait" icon={Dna}>
              <SelectField label="Species" value={f.species} onChange={set('species')}
                options={species.map((s) => ({ value: s.value, label: s.enabled ? s.label : `${s.label} — not available` }))}
                hint={sp ? (sp.enabled ? `${sp.id_help}` : sp.reason) : ''} />
              <TextField label="Trait" value={f.trait} onChange={set('trait')} placeholder="Drought tolerance" />
            </Card>
            <Card title="Gene ID list" icon={ListChecks} delay={0.06}>
              <div className="field">
                <label className="field-title">Paste gene IDs</label>
                <textarea className="input" rows={5} value={f.gene_text} onChange={(e) => set('gene_text')(e.target.value)}
                  placeholder={sp ? `${sp.id_example}\n…` : ''} style={{ resize: 'vertical', fontFamily: 'monospace' }} />
              </div>
              <div className="field">
                <label className="field-title">or a gene list file (TXT, CSV, TSV, XLSX)</label>
                <FilePicker value={f.gene_file} placeholder="No file selected" onBrowse={browseGenes} />
                {f.gene_file && <button className="btn btn-secondary sm" style={{ marginTop: 8 }} onClick={() => set('gene_file')('')}><X size={13} /> Remove</button>}
              </div>
            </Card>
          </div>

          <div className="row">
            <Card title="Options" icon={Settings2} delay={0.1}>
              <Toggle label="Refresh cached data (ask the databases again)" checked={f.refresh_cache} onChange={set('refresh_cache')} />
              <div className="field-hint">Responses are kept on this computer for 30 days, so repeating a run is fast. Switch this on to ignore them.</div>
            </Card>
          </div>

          <Card title="Database connections" icon={Wifi} delay={0.18}>
            <ApiChips items={api} />
            <div style={{ marginTop: api ? 12 : 0 }}>
              <button className="btn btn-secondary" onClick={() => bridge && bridge.check_apis()} disabled={running}><Wifi size={15} /> Check connections</button>
            </div>
          </Card>
        </TabPane>
      )}

      {/* ========================================================= SUMMARY */}
      {tab === 'summary' && (
        <TabPane id="summary" key="summary">
          <SectionHeader title="Summary">{r ? `${r.species} · ${r.trait} · ${r.run_date}` : 'Results of the last run.'}</SectionHeader>
          {!r ? <EmptyState icon={Info}>No results yet. Run the prioritization first.</EmptyState> : (
            <>
              <SummaryGrid items={[
                { label: 'Genes in the input list', value: r.summary.input_genes },
                { label: 'Candidate genes', value: r.summary.candidate_genes },
                { label: 'Successfully annotated', value: r.summary.annotated },
                { label: 'With a trait association', value: r.summary.with_association },
                { label: 'With literature support', value: r.summary.with_literature },
                { label: 'In a pathway', value: r.summary.with_pathway },
                { label: 'IDs not matched', value: r.summary.unresolved },
              ]} />
              <div className="info-line">Genes are ranked by the <b>{r.rank_by.toLowerCase()}</b>.
Trait words searched: {[...r.trait_terms.direct, ...r.trait_terms.related].join(', ')}.</div>
              {r.warnings.map((w) => <div className="info-line warn" key={w}>{w}</div>)}
              {r.unresolved.length > 0 && (
                <div className="note">
                  <p><b>Identifiers that were not matched to a gene ({r.unresolved.length}):</b></p>
                  {r.unresolved.slice(0, 30).map((u) => <p key={u.input}>{u.input} — {u.reason}</p>)}
                  {r.unresolved.length > 30 && <p>… and {r.unresolved.length - 30} more (see the Excel file).</p>}
                </div>
              )}
              <h3 style={{ margin: '22px 0 8px', fontSize: '1rem' }}>Top genes</h3>
              {r.charts.bars.some((b) => b.final)
                ? <ScoreBars bars={r.charts.bars.slice(0, 10)} onSelect={openGene} />
                : <div className="info-line">No gene has a score above 0, so there is nothing to chart. See the note above about the words searched.</div>}
              <h3 style={{ margin: '22px 0 8px', fontSize: '1rem' }}>Services</h3>
              <ApiChips items={r.api_status} />
              <div style={{ marginTop: 22, display: 'flex', gap: 12, flexWrap: 'wrap' }}>
                <button className="btn btn-primary" onClick={() => download('csv')}><Download size={16} /> Download Ranked List (CSV)</button>
                <button className="btn btn-secondary" onClick={() => download('xlsx')}><Download size={16} /> Download Workbook (Excel)</button>
                <button className="btn btn-secondary" onClick={() => download('html')}><Download size={16} /> Download Report (HTML)</button>
              </div>
              <div className="field-hint" style={{ marginTop: 8 }}>The workbook has the ranked genes, every piece of evidence with its source, link and retrieval date, publications, pathways, GO terms and the run settings.</div>
            </>
          )}
        </TabPane>
      )}

      {/* ========================================================= RANKING */}
      {tab === 'ranking' && (
        <TabPane id="ranking" key="ranking">
          <SectionHeader title="Ranked Genes">Click a column to sort, a row to see the evidence behind its score.</SectionHeader>
          {!r ? <EmptyState icon={Info}>No results yet. Run the prioritization first.</EmptyState> : (
            <>
              <div className="table-tools">
                <label className="search">
                  <Search size={15} />
                  <input value={q} onChange={(e) => { setQ(e.target.value); setShown(100) }} placeholder="Search gene, symbol, description…" />
                  {q && <button onClick={() => setQ('')} aria-label="Clear"><X size={14} /></button>}
                </label>
                <small>{rows.length.toLocaleString()} of {r.genes.length.toLocaleString()} genes</small>
              </div>
              <div className="table-wrap" style={{ maxHeight: 520 }}>
                <table className="results-table">
                  <thead>
                    <tr>{cols.map((c) => (
                      <th key={c.key} onClick={() => sortBy(c.key)} style={{ cursor: 'pointer' }}>
                        {c.label}{sort.key === c.key && (sort.dir > 0 ? <ArrowUp size={12} /> : <ArrowDown size={12} />)}
                      </th>
                    ))}</tr>
                  </thead>
                  <tbody>
                    {rows.slice(0, shown).map((g) => (
                      <tr key={g.Gene_ID} onClick={() => openGene(g.Gene_ID)} style={{ cursor: 'pointer' }}>
                        <td>{g.Rank}</td>
                        <td>{g.Gene_ID}</td>
                        <td>{g.Symbol || '—'}</td>
                        <td style={{ maxWidth: 380, overflow: 'hidden', textOverflow: 'ellipsis' }} title={g.Description}>{g.Description || '—'}</td>
                        <td><ScoreCell v={g.Final_Score} color={CAT_COLORS.association} /></td>
                        <td>{g.Evidence_Level}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              {shown < rows.length && <div className="more-row"><button className="btn btn-secondary" onClick={() => setShown((s) => s + 100)}>Show {Math.min(100, rows.length - shown)} more</button></div>}
              <div className="field-hint">“—” means the score could not be calculated because no evidence of that kind was available — it is not a score of zero. Evidence levels: Direct association · Multiple lines · Two lines · Single line · No trait evidence.</div>
            </>
          )}
        </TabPane>
      )}

      {/* ========================================================= EVIDENCE */}
      {tab === 'evidence' && (
        <TabPane id="evidence" key="evidence">
          <SectionHeader title="Gene Evidence">Everything behind the score of one gene, with its sources.</SectionHeader>
          {!r ? <EmptyState icon={Info}>No results yet. Run the prioritization first.</EmptyState> : (
            <>
              <SelectField label="Gene" value={selected || ''} onChange={setSelected}
                options={r.genes.map((g) => ({ value: g.Gene_ID, label: `${g.Rank}. ${g.Symbol ? g.Symbol + ' — ' : ''}${g.Gene_ID}` }))} />
              {gd && <EvidencePanel d={gd} r={r} bridge={bridge} />}
            </>
          )}
        </TabPane>
      )}

      {/* ========================================================= CHARTS */}
      {tab === 'charts' && (
        <TabPane id="charts" key="charts">
          <SectionHeader title="Charts">Top-ranked genes. Click a bar to open the evidence of the gene.</SectionHeader>
          {!r ? <EmptyState icon={Info}>No results yet. Run the prioritization first.</EmptyState> : (
            <>
              <h3 style={{ margin: '4px 0 8px', fontSize: '1rem' }}>Score of the top genes</h3>
              <ScoreBars bars={r.charts.bars} onSelect={openGene} />
              <h3 style={{ margin: '26px 0 8px', fontSize: '1rem' }}>What the final score is made of <small style={{ fontWeight: 400, color: 'var(--text-muted)' }}>(genes with database evidence)</small></h3>
              {r.charts.contributions.length
                ? <ContributionBars rows={r.charts.contributions} weights={r.weights} onSelect={openGene} />
                : <div className="info-line">No gene has database evidence for this trait (final score above 0), so there is nothing to break down.</div>}
              <Legend items={CAT_ORDER.map((k) => ({ label: `${r.weights[k].label} (${Math.round(r.weights[k].weight * 100)}%)`, color: CAT_COLORS[k] }))} />
              <div className="field-hint">Each colour is the share of the final score supplied by one evidence category. A category with no data does not add points and does not count against the gene: the weights of the other categories are scaled up to 100%.</div>
              <h3 style={{ margin: '26px 0 8px', fontSize: '1rem' }}>Genes and pathways <small style={{ fontWeight: 400, color: 'var(--text-muted)' }}>(best-ranked genes that belong to a pathway)</small></h3>
              {r.charts.pathways.genes.length
                ? <PathwayMatrix data={r.charts.pathways} onSelect={openGene} />
                : <div className="info-line">No candidate gene is a member of a Plant Reactome / Gramene pathway.</div>}
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}

/* ──────────────────────────── evidence panel ──────────────────────────── */

function Block({ title, children }) {
  return (
    <div style={{ margin: '18px 0 0' }}>
      <h4 style={{ margin: '0 0 6px', color: 'var(--primary)', fontSize: '.92rem' }}>{title}</h4>
      {children}
    </div>
  )
}

function EvidencePanel({ d, r, bridge }) {
  const L = (href, text) => <Ext href={href} bridge={bridge}>{text}</Ext>
  return (
    <div className="card" style={{ marginTop: 14 }}>
      <h3 style={{ fontSize: '1.15rem' }}>{d.symbol || d.id} <span style={{ fontWeight: 400, color: 'var(--text-muted)', marginLeft: 8 }}>{d.id}</span></h3>
      <div style={{ display: 'flex', flexWrap: 'wrap', gap: '4px 18px', color: 'var(--text-muted)', fontSize: '.85rem' }}>
        <span>Rank <b>{d.rank}</b></span>
        <span>{d.evidence_level}</span>
        {d.location && <span>{d.location}</span>}
        {d.biotype && <span>{d.biotype}</span>}
        <span>retrieved {d.retrieved}</span>
        {Object.entries(d.links).map(([k, u]) => <span key={k}>{L(u, k)}</span>)}
      </div>
      <p style={{ margin: '10px 0 0' }}>{d.description || <i>No description available from the databases.</i>}</p>
      {Object.keys(d.unavailable || {}).length > 0 && (
        <div className="info-line warn">Not available: {Object.entries(d.unavailable).map(([k, v]) => `${k} (${v})`).join(', ')}.</div>
      )}

      <Block title="Scores">
        <div className="table-wrap" style={{ maxHeight: 'none' }}>
          <table className="results-table">
            <thead><tr><th>Evidence category</th><th>Weight</th><th>Score (0–1)</th><th>Points</th></tr></thead>
            <tbody>
              {CAT_ORDER.filter((k) => d.scores[k].items.length > 0).map((k) => {
                const s = d.scores[k]
                return (
                  <tr key={k}>
                    <td><i style={{ display: 'inline-block', width: 10, height: 10, borderRadius: 2, background: CAT_COLORS[k], marginRight: 8 }} />{s.label}</td>
                    <td>{Math.round(s.weight * 100)}%</td>
                    <td>{s.score === null ? <span title={s.note}>unavailable</span> : s.score.toFixed(2)}</td>
                    <td>{s.score === null ? '—' : s.contribution.toFixed(1)}</td>
                  </tr>
                )
              })}
              <tr><td><b>Final score</b></td><td /><td /><td><b>{f1(d.final)}</b></td></tr>
            </tbody>
          </table>
        </div>
        <div className="field-hint">Only the evidence that was found is listed. Categories where nothing was found for this gene add no points and are not shown.</div>
        {!CAT_ORDER.some((k) => d.scores[k].items.length > 0) && <div className="info-line">No evidence for this trait was found for this gene in the databases queried.</div>}
      </Block>

      {CAT_ORDER.map((k) => {
        const s = d.scores[k]
        if (!s.items.length) return null
        return (
          <Block key={k} title={s.label}>
              <ul style={{ margin: '0 0 0 18px' }}>
                {s.items.map((it, i) => (
                  <li key={i} style={{ marginBottom: 3 }}>
                    {it.text} <span style={{ color: 'var(--text-muted)', fontSize: '.8rem' }}>— {it.source}{it.kind ? ` · ${it.kind}` : ''}</span> {it.url && L(it.url, 'source')}
                  </li>
                ))}
              </ul>
          </Block>
        )
      })}

      {d.go.length > 0 && (
        <Block title={`Gene Ontology (${d.go.length})`}>
          <div style={{ fontSize: '.86rem', lineHeight: 1.7 }}>
            {d.go.slice(0, 40).map((t) => <span key={t.id} title={`${t.namespace || ''} ${t.evidence_code || ''}`} style={{ display: 'inline-block', marginRight: 14 }}>{t.name || t.id} <span style={{ color: 'var(--text-muted)' }}>{t.id}</span></span>)}
            {d.go.length > 40 && <span style={{ color: 'var(--text-muted)' }}>… {d.go.length - 40} more</span>}
          </div>
        </Block>
      )}
      {d.domains.length > 0 && (
        <Block title="Protein domains">
          <div style={{ fontSize: '.86rem' }}>{d.domains.slice(0, 20).map((x) => x.name || x.id).join(' · ')}</div>
        </Block>
      )}
      {(d.to.length > 0 || d.qtl_to.length > 0) && (
        <Block title="Trait ontology terms (Gramene)">
          <div style={{ fontSize: '.86rem' }}>
            {d.to.map((t) => `${t.name} (${t.id})`).join(' · ')}
            {d.qtl_to.length > 0 && <div style={{ color: 'var(--text-muted)' }}>Positional, from QTL covering the gene: {d.qtl_to.slice(0, 12).map((t) => t.name).join(' · ')}</div>}
          </div>
        </Block>
      )}
      {d.oryzabase_entry && (
        <Block title="Oryzabase (curated rice gene)">
          <div style={{ fontSize: '.86rem' }}>
            <div><b>{d.oryzabase_entry.symbol}</b> {d.oryzabase_entry.name} {d.oryzabase_entry.url && L(d.oryzabase_entry.url, 'Oryzabase')}</div>
            {d.oryzabase_entry.classes.length > 0 && <div>Trait classes: {d.oryzabase_entry.classes.join(' · ')}</div>}
            {d.oryzabase_entry.to.length > 0 && <div>Trait Ontology: {d.oryzabase_entry.to.map((t) => `${t.name} (${t.id})`).join(' · ')}</div>}
            {d.oryzabase_entry.synonyms.length > 0 && <div style={{ color: 'var(--text-muted)' }}>Also known as: {d.oryzabase_entry.synonyms.join(', ')}</div>}
            {d.oryzabase_entry.explanation && <div style={{ color: 'var(--text-muted)' }}>{d.oryzabase_entry.explanation}</div>}
          </div>
        </Block>
      )}
      {d.pathways.length > 0 && (
      <Block title={`Pathways (${d.pathways.length})`}>
          <ul style={{ margin: '0 0 0 18px' }}>{d.pathways.map((p) => <li key={p.id}>{p.name} {p.url && L(p.url, p.id)}</li>)}</ul>
      </Block>
      )}
      {d.orthologs.length > 0 && (
      <Block title={`Orthologs (${d.orthologs.length})`}>
          <ul style={{ margin: '0 0 0 18px' }}>
            {d.orthologs.map((o) => <li key={o.id}>{o.name ? `${o.name} · ` : ''}{o.url ? L(o.url, o.id) : o.id} <span style={{ color: 'var(--text-muted)', fontSize: '.8rem' }}>— {o.species}{o.description ? ` · ${o.description}` : ''}</span></li>)}
          </ul>
      </Block>
      )}
      {d.publications.length > 0 && (
      <Block title={`Publications (${d.publications.length})`}>
          <>
            <ul style={{ margin: '0 0 0 18px' }}>
              {d.publications.map((p) => (
                <li key={p.pmid} style={{ marginBottom: 3 }}>
                  {p.title} <span style={{ color: 'var(--text-muted)', fontSize: '.8rem' }}>— {p.authors}, {p.journal} {p.year} · {p.relation}</span> {L(p.url, `PMID ${p.pmid}`)}
                </li>
              ))}
            </ul>
            <div className="field-hint">“Co-mention” means the gene and the trait are named in the same paper; it does not show that the paper demonstrates a relationship.</div>
          </>
      </Block>
      )}
    </div>
  )
}
