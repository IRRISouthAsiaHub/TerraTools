import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import {
  Fingerprint, SlidersHorizontal, LayoutDashboard, Dna, Table2, CircleHelp, Filter, Wheat, Download, ArrowRight, Settings2,
  ListChecks, Sprout, TriangleAlert, Grip, FileUp, Info,
} from 'lucide-react'
import { ToolPage, TabPane, SectionHeader, Card, RadioChips, Note, DataTable, EmptyState, FilePicker, SummaryGrid } from '../components/kit.jsx'
import { useFiles } from '../lib/files.js'
import { useBridge, useTool, callJson, mockGeneric } from '../lib/bridge.js'
import {
  SNPS, VARIETIES, LOST_GENES, loadDataset, DEFAULT_MODE, TRAITS, STATUSES, STATUS_COLOURS, STATUS_SHORT, REGION_LABEL,
  traitColour, geneList, snpsFor, tally, genesOf, placeGenes, statusOf, rawCall,
} from '../lib/favAllele.js'
import './favallele/favallele.css'
import RicePlant from './favallele/RicePlant.jsx'
import { TraitBars, GeneStack, Heatmap, StatusLegend, seriesColour } from './favallele/Charts.jsx'
import { MultiSelect, VarietyPicker, GeneDialog } from './favallele/Pickers.jsx'
import HeatmapTab from './favallele/HeatmapTab.jsx'

const TABS = [
  { id: 'setup', label: 'Your 4K File', icon: SlidersHorizontal },
  { id: 'overview', label: 'Overview', icon: LayoutDashboard },
  { id: 'genes', label: 'Genes', icon: Dna },
  { id: 'table', label: 'SNP Table', icon: Table2 },
  { id: 'heatmap', label: 'Heatmap', icon: Grip },
  { id: 'help', label: 'About', icon: CircleHelp },
]
const VIEWS = [
  { value: 'all', label: 'All statuses' }, { value: 'favourable', label: 'Favourable' }, { value: 'unfavourable', label: 'Unfavourable' },
  { value: 'heterozygous', label: 'Heterozygous' }, { value: 'missing', label: 'Missing' },
]
const CARDS = [
  ['Fav SNPs', 'Homozygous favourable', '#2f9e5f'], ['Unfav SNPs', 'Homozygous unfavourable', '#e5484d'],
  ['Fav genes', 'Genes with ≥1 fav allele', '#12a394'], ['Unfav genes', 'Genes with ≥1 unfav allele', '#d1603d'],
  ['Het', 'Heterozygous', '#d9980b'], ['Missing', 'No call', '#7f8a84'],
]
const GENE_STATUS_COLOUR = { Favourable: '#2f9e5f', Unfavourable: '#e5484d', Mixed: '#5b6ee1', Heterozygous: '#f2b01e', Missing: '#b4beb8' }
const geneStatus = (c) => (c.Favourable && c.Unfavourable ? 'Mixed' : c.Favourable ? 'Favourable' : c.Unfavourable ? 'Unfavourable' : c.Heterozygous ? 'Heterozygous' : 'Missing')
const firstThree = () => VARIETIES.slice(0, 3).map((v) => v.id)

export default function FavourableAlleleTool() {
  const [tab, setTab] = useState('setup')
  const [ds, setDs] = useState(0)            // grows each time a file is loaded; the data itself lives in lib/favAllele.js
  const [file, setFile] = useState('')
  const [match, setMatch] = useState(null)
  const [report, setReport] = useState([])
  const [fileMsg, setFileMsg] = useState('')
  const [info, setInfo] = useState(null)
  const [traits, setTraits] = useState([])
  const [genes, setGenes] = useState([])
  const [view, setView] = useState('all')
  const mode = DEFAULT_MODE
  const [selected, setSelected] = useState([])
  const [active, setActive] = useState(null)
  const [metric, setMetric] = useState('Favourable')
  const [dotBy, setDotBy] = useState('trait')
  const [focusTrait, setFocusTrait] = useState(null)
  const [hover, setHover] = useState(null)
  const [dialog, setDialog] = useState(null)
  const [topGenes, setTopGenes] = useState('25')
  const [heatGenes, setHeatGenes] = useState('40')
  const [showLost, setShowLost] = useState(false)
  const plantRef = useRef(null)
  const [hmGenes, setHmGenes] = useState([])
  const [hmLines, setHmLines] = useState([])
  const [hmDpi, setHmDpi] = useState(300)
  const { bridge: hmBridge } = useBridge('favourableAlleleBridge', () => mockGeneric())
  const hm = useTool(hmBridge, 'run_heatmap')
  const { saveImage } = useFiles()

  /* ───── derived data ───── */
  /* every memo that reads the loaded data lists `ds` so that it is recomputed when a new file is loaded */
  const byId = useMemo(() => new Map(VARIETIES.map((v) => [v.id, v])), [ds]) // eslint-disable-line react-hooks/exhaustive-deps
  const selVarieties = useMemo(() => selected.map((id) => byId.get(id)).filter(Boolean), [selected, byId])
  const activeV = byId.get(active) || selVarieties[0] || null
  const filtered = useMemo(() => snpsFor(traits, genes), [traits, genes, ds]) // eslint-disable-line react-hooks/exhaustive-deps
  const geneOptions = useMemo(() => geneList(traits).map(([n, t]) => ({ value: n, label: n, hint: t })), [traits, ds]) // eslint-disable-line react-hooks/exhaustive-deps
  const hmGeneTrait = useMemo(() => new Map(geneList([]).map(([n, t]) => [n, t])), [ds]) // eslint-disable-line react-hooks/exhaustive-deps
  const hmGeneOptions = useMemo(() => geneList([]).map(([n, t]) => ({ value: n, label: n, hint: t })), [ds]) // eslint-disable-line react-hooks/exhaustive-deps
  const hmLineOptions = useMemo(() => VARIETIES.map((v) => ({ value: v.id, label: v.name })), [ds]) // eslint-disable-line react-hooks/exhaustive-deps
  const traitOptions = useMemo(() => TRAITS.map((t) => ({ value: t, label: t })), [ds]) // eslint-disable-line react-hooks/exhaustive-deps
  const loaded = ds > 0
  const lostGenes = info ? info.lost_genes : LOST_GENES

  const t = useMemo(() => (activeV ? tally(filtered, activeV, mode, view) : null), [filtered, activeV, mode, view])
  const geneObjs = useMemo(() => {
    if (!activeV) return []
    const want = view === 'all' ? null : STATUSES.find((s) => s.toLowerCase() === view)
    return genesOf(filtered, activeV, mode).filter((g) => !want || g.counts[want] > 0)
  }, [filtered, activeV, mode, view])
  const placed = useMemo(() => placeGenes(geneObjs), [geneObjs])
  const dots = useMemo(() => placed.map(({ gene, x, y }) => {
    const colour = dotBy === 'trait' ? traitColour(gene.trait) : GENE_STATUS_COLOUR[geneStatus(gene.counts)]
    return { gene, x, y, colour, dim: dotBy === 'trait' && focusTrait && gene.trait !== focusTrait }
  }), [placed, dotBy, focusTrait])
  const plantTraits = useMemo(() => [...new Set(geneObjs.map((g) => g.trait))].sort(), [geneObjs])
  const hoverGene = hover ? geneObjs.find((g) => g.name === hover) : null

  const traitChart = useMemo(() => {
    const names = [...new Set(filtered.map((s) => s.traitLabel))].sort()
    const series = selVarieties.map((v, i) => {
      const counts = Object.fromEntries(names.map((n) => [n, 0]))
      filtered.forEach((s) => { if (statusOf(s, v, mode).status === metric) counts[s.traitLabel]++ })
      return { id: v.id, name: v.name, colour: seriesColour(i), values: names.map((n) => counts[n]) }
    })
    return { names, series }
  }, [filtered, selVarieties, mode, metric])

  const geneBars = useMemo(() => {
    if (!t) return []
    const n = topGenes === 'all' ? Infinity : Number(topGenes)
    return [...t.byGene].filter((g) => g.total > 0).sort((a, b) => (b.Favourable - a.Favourable) || (a.Unfavourable - b.Unfavourable)).slice(0, n)
  }, [t, topGenes])

  const heatRows = useMemo(() => {
    const by = new Map()
    filtered.forEach((s) => { if (!by.has(s.geneLabel)) by.set(s.geneLabel, []); by.get(s.geneLabel).push(s) })
    const ranked = [...by.entries()].map(([gene, snps]) => {
      let typed = 0
      snps.forEach((s) => selVarieties.forEach((v) => { if (statusOf(s, v, mode).status !== 'Missing') typed++ }))
      return { gene, snps, typed }
    }).sort((a, b) => b.typed - a.typed)
    const n = heatGenes === 'all' ? Infinity : Number(heatGenes)
    return ranked.slice(0, n).flatMap(({ gene, snps }) => snps.map((snp, k) => ({ gene, snp, n: snps.length, k })))
  }, [filtered, selVarieties, mode, heatGenes])

  const tableRows = useMemo(() => (t ? t.rows.map(({ snp, status, call, gt }) => ({
    Gene: snp.geneLabel, 'Gene ID': snp.gene, Trait: snp.traitLabel, Chr: snp.chrom, Position: snp.pos,
    REF: snp.ref || '-', ALT: snp.alt || '-', Favourable: snp.fav || '-', GT: gt < 0 ? '-' : gt, Call: call, Status: status,
  })) : []), [t])

  /* ───── balloon heatmap (drawn by the R script in tools/favourable_allele/backend) ───── */
  const runHeatmap = () => {
    if (!hmGenes.length) return hm.fail('Choose at least one gene for the heatmap.')
    if (!hmLines.length) return hm.fail('Choose at least one line for the heatmap.')
    const rank = new Map(hmGenes.map((g, i) => [g, i]))
    const snps = SNPS.filter((s) => rank.has(s.geneLabel)).sort((a, b) => (rank.get(a.geneLabel) - rank.get(b.geneLabel)) || (a.i - b.i))
    const lines = hmLines.map((id) => byId.get(id)).filter(Boolean)
    hm.run({
      dpi: hmDpi,
      snps: snps.map((s) => ({ gene: s.gene, name: s.geneLabel, trait: s.traitLabel, sub: s.sub, marker: s.marker, chrom: s.chrom, pos: s.pos, ref: s.ref, alt: s.alt, fav: s.fav, unfav: s.unfav })),
      varieties: lines.map((v) => ({ name: v.name, group: 'Variety' })), group_label: 'Lines',
      calls: lines.map((v) => snps.map((s) => { const g = rawCall(v.i, s.i); return g < 0 ? 'x' : String(g) }).join('')),
    })
  }

  /* ───── actions ───── */
  useEffect(() => {
    if (hmBridge) callJson(hmBridge, 'get_info').then(setInfo).catch(() => {})
  }, [hmBridge])
  /** Choose the user's file; the backend matches its marker IDs against the built-in trait-gene markers. */
  const browseFile = async () => {
    const p = await hm.browse()
    if (!p) return
    setFile(p); setFileMsg('Matching the marker IDs…')
    let r
    try { r = await callJson(hmBridge, 'load_genotypes', p) } catch { r = { status: 'error', message: 'The backend did not answer (preview mode has no Python backend).' } }
    if (!r || r.status !== 'success') { setMatch(null); setFileMsg((r && r.message) || 'The file could not be read.'); return }
    loadDataset(r)
    const ids = firstThree()
    setMatch(r.match); setReport(r.report); setFileMsg('')
    setTraits([]); setGenes([]); setView('all'); setSelected(ids); setActive(ids[0] || null)
    setHmGenes([]); setHmLines([]); hm.setResult(null); hm.setState('idle'); hm.setStatus('Ready')
    setDs((n) => n + 1)
  }
  const toggleVariety = (id) => {
    const next = selected.includes(id) ? selected.filter((x) => x !== id) : [...selected, id]
    setSelected(next)
    if (!next.includes(active)) setActive(next[0] || null)
  }
  const setMany = (ids) => { setSelected(ids); setActive(ids[0] || null) }
  const changeTraits = (v) => {
    setTraits(v)
    const valid = new Set(geneList(v).map(([n]) => n))
    setGenes((g) => g.filter((x) => valid.has(x)))
  }
  const openGene = useCallback((name) => {
    const g = geneObjs.find((x) => x.name === name)
    if (g) setDialog(g)
  }, [geneObjs])

  const savePlantPng = () => {
    const svg = plantRef.current
    if (!svg) return
    const xml = new XMLSerializer().serializeToString(svg)
    const url = URL.createObjectURL(new Blob([xml], { type: 'image/svg+xml;charset=utf-8' }))
    const img = new Image()
    img.onload = () => {
      const c = document.createElement('canvas')
      c.width = 1500; c.height = 2040
      const ctx = c.getContext('2d')
      ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, c.width, c.height)
      ctx.drawImage(img, 0, 0, c.width, c.height)
      URL.revokeObjectURL(url)
      saveImage(c.toDataURL('image/png'), `${activeV ? activeV.name : 'line'}_gene_map.png`)
    }
    img.src = url
  }
  const saveMatches = () => {
    const cols = ['file_id', 'marker', 'gene', 'trait', 'method', 'ref', 'alt', 'favourable']
    const head = ['ID in your file', 'Marker matched', 'Gene', 'Trait', 'Matched by', 'REF', 'ALT', 'Favourable allele']
    const esc = (v) => { const t = String(v ?? ''); return /[",\n]/.test(t) ? `"${t.replace(/"/g, '""')}"` : t }
    const csv = [head.join(','), ...report.map((r) => cols.map((c) => esc(r[c])).join(','))].join('\n')
    saveImage(`data:text/csv;charset=utf-8,${encodeURIComponent(`\ufeff${csv}`)}`, 'matched_marker_IDs.csv')
  }
  const saveCsv = () => {
    if (!tableRows.length) return
    const cols = Object.keys(tableRows[0])
    const esc = (v) => { const s = String(v ?? ''); return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s }
    const csv = [cols.join(','), ...tableRows.map((r) => cols.map((c) => esc(r[c])).join(','))].join('\n')
    saveImage(`data:text/csv;charset=utf-8,${encodeURIComponent(`﻿${csv}`)}`, `${activeV ? activeV.name : 'line'}_SNP_status.csv`)
  }

  /* ───── pieces ───── */
  const filterSummary = [
    traits.length ? `${traits.length} trait${traits.length > 1 ? 's' : ''}` : 'all traits',
    genes.length ? `${genes.length} gene${genes.length > 1 ? 's' : ''}` : 'all genes',
    view !== 'all' ? `${VIEWS.find((v) => v.value === view).label} only` : null,
  ].filter(Boolean).join(' · ')

  const bar = !loaded ? null : (
    <div className="fa-bar">
      <div className="fa-bar-pills">
        {selVarieties.map((v, i) => (
          <button key={v.id} type="button" className={`fa-bpill ${v.id === active ? 'active' : ''}`} style={{ '--c': seriesColour(i) }} onClick={() => setActive(v.id)}
            title="Show this line on the plant, in the gene chart and in the table">
            <i />{v.name}
          </button>
        ))}
        {!selVarieties.length && <span className="fa-hint">No line selected.</span>}
      </div>
      <button type="button" className="btn btn-secondary sm" onClick={() => setTab('setup')}><Settings2 size={14} /> Change</button>
      <div className="fa-bar-sub"><Filter size={13} /> {filterSummary} · {filtered.length} SNPs in {new Set(filtered.map((s) => s.geneLabel)).size} genes</div>
    </div>
  )

  const noVariety = (
    <EmptyState icon={Fingerprint}>
      {loaded ? 'Select at least one line first.' : 'Upload your 4K haplopanel file first.'}<br />
      <button type="button" className="btn btn-secondary" style={{ marginTop: 14 }} onClick={() => setTab('setup')}>{loaded ? 'Go to Your 4K File' : 'Upload the file'}</button>
    </EmptyState>
  )

  return (
    <ToolPage
      id="favallele" icon={Fingerprint} tabs={TABS.map((x) => (x.id === 'overview' && selVarieties.length ? { ...x, badge: selVarieties.length } : x))} tab={tab} setTab={setTab}
      title="Favourable Allele Analysis"
      subtitle="Favourable alleles of 166 known trait genes in your rice 4K haplopanel genotypes"
      run={tab === 'heatmap' ? { label: 'Draw Balloon Heatmap', runningLabel: 'Drawing…', state: hm.state, status: hm.status, onRun: runHeatmap } : undefined}
    >
      {tab === 'setup' && (
        <TabPane id="setup" key="setup">
          <SectionHeader title="Your 4K haplopanel file">
            Upload the genotypes of your lines. TerraTools matches the marker IDs against its own list of trait-gene markers and shows the favourable alleles at the markers it finds.
          </SectionHeader>
          <Note>
            <b>Scope.</b> This tool is for the <b>rice 4K haplopanel</b> only, for now. It is built on {info ? info.n_snps : 327} SNPs and indels in {info ? info.n_genes : 166} known trait genes
            that are kept inside TerraTools; it has <b>not been scaled yet</b> to other marker panels, other crops or other genes. Markers of your file that are not among them are ignored.
          </Note>
          <Card title="Genotype file" icon={FileUp} delay={0.03}>
            <FilePicker value={file} onBrowse={browseFile} placeholder="Select your 4K haplopanel genotype file (CSV, TXT, VCF or VCF.GZ)" />
            {fileMsg && <div className={`fa-filemsg ${match ? '' : 'err'}`}>{fileMsg}</div>}
            <ul className="fa-list" style={{ marginTop: 14 }}>
              <li><b>Layout:</b> one row per marker, one column per line. The first row holds the line names.</li>
              <li><b>Marker IDs</b> are matched as <span className="fa-mono">chr1_7268879</span> (chromosome_position), also when followed by the gene (<span className="fa-mono">chr1_7268879_LOC_Os01g13060</span>); chromosome and position columns are used as a second chance.</li>
              <li><b>Genotypes:</b> <span className="fa-mono">0</span> = REF/REF, <span className="fa-mono">1</span> = heterozygous, <span className="fa-mono">2</span> = ALT/ALT, <span className="fa-mono">-1</span> = missing. VCF genotypes (<span className="fa-mono">0/1</span>, <span className="fa-mono">./.</span>) and nucleotides (<span className="fa-mono">A/G</span>) are converted to these codes.</li>
            </ul>
          </Card>

          {!loaded ? (
            <EmptyState icon={Fingerprint}>Select your file to see which of its markers match the trait genes.</EmptyState>
          ) : (
            <>
              <SummaryGrid items={[
                { label: 'Markers in your file', value: match.file_markers.toLocaleString() },
                { label: 'Matched IDs', value: `${match.matched_markers} / ${match.panel_markers}` },
                { label: 'Trait genes found', value: `${match.genes_found} / ${match.genes_total}` },
                { label: 'Lines', value: match.lines },
              ]} />
              <Card title="Matched marker IDs" icon={ListChecks} delay={0.05}>
                <Note>
                  <b>{match.matched_markers}</b> of the {match.panel_markers} trait-gene markers were found in your file (matched by {Object.entries(match.by_method).map(([m, n]) => `${m}: ${n}`).join(', ')});
                  {' '}<b>{match.unmatched_file_markers.toLocaleString()}</b> markers of your file are not trait-gene markers and are ignored. IDs were read from the column “{match.id_column}”.
                  {' '}{match.calls_made.toLocaleString()} of {match.calls_total.toLocaleString()} genotype calls ({Math.round((100 * match.calls_made) / Math.max(1, match.calls_total))} %) are present.
                  {match.genes_missing.length > 0 && <> No marker was found for {match.genes_missing.length} gene{match.genes_missing.length > 1 ? 's' : ''}: {match.genes_missing.slice(0, 12).join(', ')}{match.genes_missing.length > 12 ? ', …' : ''}.</>}
                  {match.n_allele_mismatch > 0 && <> {match.n_allele_mismatch} marker{match.n_allele_mismatch > 1 ? 's were' : ' was'} left out because the REF/ALT alleles in your file disagree with ours (e.g. {match.allele_mismatch[0].file_id}).</>}
                </Note>
                <DataTable rows={report.map((r) => ({ 'ID in your file': r.file_id, 'Marker matched': r.marker, Gene: r.gene, Trait: r.trait, 'Matched by': r.method }))} maxHeight={320} />
                <div className="fa-quick-row">
                  <button type="button" className="btn btn-secondary sm" onClick={saveMatches}><Download size={14} /> Download matched IDs (CSV)</button>
                </div>
              </Card>

              <div className="grid-2">
                <Card title="Filters" icon={Filter} className="fa-filters-card">
                  <MultiSelect label="Traits" allLabel="All traits" noun="traits" options={traitOptions} value={traits} onChange={changeTraits} dot={traitColour} />
                  <MultiSelect label="Genes" allLabel="All genes" noun="genes" options={geneOptions} value={genes} onChange={setGenes} search />
                  <div className="field"><label className="field-title">View status</label><RadioChips options={VIEWS} value={view} onChange={setView} /></div>
                </Card>
                <Card title="Genotype codes" icon={Settings2}>
                  <ul className="fa-list">
                    {[['0', 'homozygous reference (REF/REF)'], ['1', 'heterozygous'], ['2', 'homozygous alternate (ALT/ALT)'], ['−1', 'missing (no call)']].map(([c, d]) => (
                      <li key={c}><b className="fa-mono">{c}</b> = {d}</li>
                    ))}
                  </ul>
                  <Note>
                    A SNP is <b>favourable</b> when the line is homozygous for the favourable allele of the literature, <b>unfavourable</b> when homozygous for the other allele,
                    <b> heterozygous</b> when it carries both, and <b>missing</b> without a call.
                  </Note>
                </Card>
              </div>
              <Card title="Select lines to compare" icon={ListChecks} delay={0.05}>
                <VarietyPicker varieties={VARIETIES} selected={selected} active={active} onToggle={toggleVariety} onActive={setActive}
                  onClear={() => setMany(selected.slice(0, 1))} onMany={setMany} />
              </Card>
              <div className="fa-go">
                <button type="button" className="btn-run" disabled={!selVarieties.length} onClick={() => setTab('overview')}>
                  <span>Show results for {selVarieties.length} {selVarieties.length === 1 ? 'line' : 'lines'}</span><ArrowRight size={18} />
                </button>
              </div>
            </>
          )}
        </TabPane>
      )}

      {tab === 'overview' && (
        <TabPane id="overview" key="overview">
          {bar}
          {!activeV ? noVariety : (
            <>
              <div className="fa-stats">
                {CARDS.map(([label, sub, c], i) => {
                  const v = [t.counts.Favourable, t.counts.Unfavourable, t.favGenes, t.unfavGenes, t.counts.Heterozygous, t.counts.Missing][i]
                  return (
                    <div key={label} className="fa-stat" style={{ '--c': c }}>
                      <span className="l">{label}</span><b>{v}</b><small>{sub}</small>
                    </div>
                  )
                })}
              </div>
              <div className="fa-main">
                <Card title={`Line profile — ${activeV.name}`} icon={Sprout}>
                  <div className="fa-plant-opts">
                    <span>Colour dots by</span>
                    <RadioChips options={[{ value: 'trait', label: 'Trait' }, { value: 'status', label: 'Allele status' }]} value={dotBy} onChange={setDotBy} />
                  </div>
                  <div className="fa-plant-scene">
                    <RicePlant dots={dots} hover={hover} onHover={setHover} onSelect={(g) => setDialog(g)} svgRef={plantRef} />
                  </div>
                  <div className="fa-tip">
                    {hoverGene ? (
                      <>
                        <b>{hoverGene.name}</b>
                        <span style={{ color: traitColour(hoverGene.trait) }}>{hoverGene.trait} · {REGION_LABEL[hoverGene.region]} · {hoverGene.snps.length} SNP{hoverGene.snps.length > 1 ? 's' : ''}</span>
                        <span>
                          {STATUSES.map((s) => <em key={s}><i style={{ background: STATUS_COLOURS[s] }} />{hoverGene.counts[s]} {STATUS_SHORT[s]}</em>)}
                        </span>
                      </>
                    ) : <span>Hover a dot for the gene · click it for every SNP and every selected line.</span>}
                  </div>
                  {dotBy === 'trait' ? (
                    <div className="fa-trait-legend">
                      {plantTraits.map((tr) => (
                        <button key={tr} type="button" className={focusTrait === tr ? 'on' : ''} onClick={() => setFocusTrait(focusTrait === tr ? null : tr)}
                          title="Click to highlight this trait">
                          <i style={{ background: traitColour(tr) }} />{tr}
                        </button>
                      ))}
                    </div>
                  ) : (
                    <div className="fa-legend">
                      {Object.entries(GENE_STATUS_COLOUR).map(([s, c]) => <span key={s}><i style={{ background: c }} />{s === 'Mixed' ? 'Mixed (fav + unfav)' : s}</span>)}
                    </div>
                  )}
                  <div className="fa-actions"><button type="button" className="btn btn-secondary sm" onClick={savePlantPng}><Download size={14} /> Download PNG</button></div>
                </Card>
                <Card title="SNPs by trait" icon={Wheat} delay={0.05}>
                  <div className="fa-plant-opts">
                    <span>Count</span>
                    <RadioChips options={STATUSES} value={metric} onChange={setMetric} />
                  </div>
                  <p className="fa-hint" style={{ marginBottom: 8 }}>One bar per selected line.</p>
                  <TraitBars traits={traitChart.names} series={traitChart.series} metric={metric} />
                  <div className="fa-legend">{selVarieties.map((v, i) => <span key={v.id}><i style={{ background: seriesColour(i) }} />{v.name}</span>)}</div>
                </Card>
              </div>
            </>
          )}
        </TabPane>
      )}

      {tab === 'genes' && (
        <TabPane id="genes" key="genes">
          {bar}
          {!activeV ? noVariety : (
            <>
              <Card title={`Allele status by gene — ${activeV.name}`} icon={Dna}>
                <div className="fa-plant-opts">
                  <span>Show</span>
                  <RadioChips options={[{ value: '15', label: 'Top 15' }, { value: '25', label: 'Top 25' }, { value: '50', label: 'Top 50' }, { value: 'all', label: 'All genes' }]} value={topGenes} onChange={setTopGenes} />
                </div>
                <p className="fa-hint" style={{ marginBottom: 8 }}>Sorted by favourable SNP count; bar length = SNPs typed in the gene. Click a bar for its SNPs.</p>
                <GeneStack genes={geneBars} onGene={openGene} />
                <StatusLegend />
              </Card>
              <Card title="Allele map of the selected lines" icon={Table2} delay={0.05}>
                <div className="fa-plant-opts">
                  <span>Genes (best typed first)</span>
                  <RadioChips options={[{ value: '20', label: '20' }, { value: '40', label: '40' }, { value: '80', label: '80' }, { value: 'all', label: 'All' }]} value={heatGenes} onChange={setHeatGenes} />
                </div>
                <p className="fa-hint" style={{ marginBottom: 8 }}>Rows are the selected varieties, columns the SNPs of each gene. Click a cell or a gene name for the details.</p>
                <Heatmap rows={heatRows} varieties={selVarieties} mode={mode} onGene={openGene} />
                <StatusLegend />
              </Card>
            </>
          )}
        </TabPane>
      )}

      {tab === 'table' && (
        <TabPane id="table" key="table">
          {bar}
          {!activeV ? noVariety : (
            <Card title={`SNP-level detail — ${activeV.name}`} icon={Table2}>
              <div className="fa-actions top">
                <button type="button" className="btn btn-secondary sm" onClick={saveCsv}><Download size={14} /> Download CSV</button>
              </div>
              <DataTable rows={tableRows} maxHeight={520}
                renderCell={(c, v) => {
                  if (c === 'Status') return <span className="fa-chip" style={{ '--c': STATUS_COLOURS[v] }}>{v}</span>
                  if (c === 'Favourable') return <b style={{ color: STATUS_COLOURS.Favourable }}>{v}</b>
                  if (['REF', 'ALT', 'Call', 'GT', 'Position'].includes(c)) return <span className="fa-mono">{c === 'Position' ? Number(v).toLocaleString() : v}</span>
                  return undefined
                }} />
            </Card>
          )}
        </TabPane>
      )}

      {tab === 'heatmap' && (
        <TabPane id="heatmap" key="heatmap">
          <SectionHeader title="Balloon heatmap">
            Pick the genes and the lines you want to compare. For every line and gene a balloon shows the fraction of the gene’s markers at which the line is homozygous for the favourable allele.
          </SectionHeader>
          {!loaded ? noVariety : <HeatmapTab
            geneOptions={hmGeneOptions} geneTrait={hmGeneTrait} genes={hmGenes} setGenes={setHmGenes}
            lineOptions={hmLineOptions} lines={hmLines} setLines={setHmLines} dpi={hmDpi} setDpi={setHmDpi}
            filterGenes={[...new Set(filtered.map((s) => s.geneLabel))]} filterLines={selected} tool={hm}
          />}
        </TabPane>
      )}
      {tab === 'help' && (
        <TabPane id="help" key="help">
          <SectionHeader title="About this analysis">Which favourable alleles of 166 known trait genes do the lines of your rice 4K haplopanel file carry?</SectionHeader>
          <div className="grid-2">
            <Card title="Reading the results" icon={CircleHelp}>
              <ul className="fa-list">
                <li><b>Your 4K File</b> — upload your genotype file; the marker IDs are matched against the trait-gene markers kept inside TerraTools and the matched IDs are listed. Then restrict the traits and genes and tick the lines to compare. The first one is shown in detail; click another name at the top of any page to switch.</li>
                <li><b>Overview</b> — the counts for the shown line, the rice plant with one dot per gene (placed by the trait’s plant part), and the SNP counts by trait for every selected line.</li>
                <li><b>Genes</b> — a bar for each gene with its SNP statuses, and an allele map (lines × SNPs) to see at a glance where the lines differ.</li>
                <li><b>SNP Table</b> — every SNP with its REF, ALT, favourable allele, stored code, call and status; filter it and download it as CSV.</li>
              </ul>
            </Card>
            <Card title="Status of a SNP" icon={ListChecks}>
              <ul className="fa-list">
                {[['Favourable', 'The line is homozygous for the favourable allele.'], ['Unfavourable', 'Homozygous for the other allele.'],
                  ['Heterozygous', 'Carries both alleles.'], ['Missing', 'No call, or the favourable allele matches neither REF nor ALT.']].map(([s, d]) => (
                  <li key={s}><i className="fa-sw" style={{ background: STATUS_COLOURS[s] }} /><b>{s}</b> — {d}</li>
                ))}
              </ul>
              <Note>
                The favourable allele is compared with REF and ALT; when it only matches the reverse complement of one of them, the annotation is taken to be on the opposite
                strand. For A/T and C/G SNPs the strand cannot be told, so the allele is used as written. “Fav genes” are the genes with at least one favourable SNP, so a gene
                with both kinds of SNP is counted in both.
              </Note>
            </Card>
          </div>
          <Card title={`${lostGenes.length} genes without a scorable marker`} icon={TriangleAlert} delay={0.05}>
            <Note>
              Every marker of these genes is an indel whose favourable allele cannot be matched to REF or ALT, so they are not part of the analysis.
              {' '}<button type="button" className="fa-link" onClick={() => setShowLost((s) => !s)}>{showLost ? 'Hide the list' : 'Show the list'}</button>
            </Note>
            {showLost && <div className="fa-lost">{lostGenes.map((g) => <span key={g}>{g}</span>)}</div>}
          </Card>
        </TabPane>
      )}

      <GeneDialog gene={dialog} varieties={selVarieties} mode={mode} onClose={() => setDialog(null)} />
    </ToolPage>
  )
}
