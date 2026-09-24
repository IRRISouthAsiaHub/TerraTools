import { useEffect, useRef, useState } from 'react'
import {
  Dna, Upload, SlidersHorizontal, ClipboardCheck, Table2, ChartColumn, ChartArea, CircleHelp, FolderOpen,
  Settings2, BookOpen, Trash2, Download, FileSpreadsheet, RefreshCw, Camera, Map as MapIcon, Thermometer, ChartLine,
} from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Slider, SelectField, TextField, Toggle, ChipGroup, FilePicker, Note,
  EmptyState, ImageGallery, DataTable, LogBox,
} from '../components/kit.jsx'
import { useBridge, useTool, useGlobals, useSignal, callJson, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'upload', label: 'Upload', icon: Upload },
  { id: 'parameters', label: 'Parameters', icon: SlidersHorizontal },
  { id: 'qc', label: 'QC Report', icon: ClipboardCheck },
  { id: 'results', label: 'Results', icon: Table2 },
  { id: 'mapplot', label: 'Map Plot', icon: ChartColumn },
  { id: 'density', label: 'Density & RF', icon: ChartArea },
  { id: 'cyto', label: 'Cytogenetic', icon: Dna },
  { id: 'help', label: 'Help', icon: CircleHelp },
]
const PALETTES = ['Set1', 'Set2', 'Set3', 'Paired', 'Dark2', 'Accent', 'Spectral', 'RdYlBu', 'PRGn']
const CYTO_SCHEMES = {
  classic: { a: '#3B82F6', b: '#EF4444', ab: '#F59E0B', miss: 'rgba(148,163,184,0.4)', bg: '#0F172A' },
  blue_orange: { a: '#2563EB', b: '#EA580C', ab: '#FACC15', miss: 'rgba(148,163,184,0.4)', bg: '#0F172A' },
  green_purple: { a: '#16A34A', b: '#9333EA', ab: '#EAB308', miss: 'rgba(148,163,184,0.4)', bg: '#0F172A' },
  red_blue: { a: '#DC2626', b: '#2563EB', ab: '#D97706', miss: 'rgba(148,163,184,0.4)', bg: '#0F172A' },
  pastel: { a: '#93C5FD', b: '#FCA5A5', ab: '#FDE68A', miss: 'rgba(203,213,225,0.5)', bg: '#1E293B' },
}

/**
 * R/jsonlite tables can arrive as row objects, row arrays, or a column-oriented
 * object; normalise to row objects keyed by `columns`.
 */
function normRows(rows, columns) {
  if (!rows) return []
  let out = []
  if (Array.isArray(rows)) out = rows
  else if (typeof rows === 'object') {
    const keys = Object.keys(rows)
    if (keys.some((k) => Array.isArray(rows[k]))) {
      const n = Math.max(...keys.map((k) => (Array.isArray(rows[k]) ? rows[k].length : 0)))
      for (let i = 0; i < n; i++) out.push(Object.fromEntries(keys.map((k) => [k, Array.isArray(rows[k]) ? rows[k][i] : rows[k]])))
    } else out = [rows]
  }
  if (!columns) return out
  const cols = columns.map(String)
  return out.map((row) => Object.fromEntries(cols.map((c, i) => {
    if (Array.isArray(row)) return [c, row[i]]
    let v = row[c]
    if (v === undefined) { const k = Object.keys(row).find((x) => x.toLowerCase() === c.toLowerCase()); if (k !== undefined) v = row[k] }
    return [c, v]
  })))
}
const dash = (_c, v) => (v === null || v === undefined || v === '' ? '—' : undefined)

function Tbl({ rows, columns, maxHeight = 380 }) {
  const data = normRows(rows, columns)
  return data.length ? <DataTable rows={data} columns={columns && columns.map(String)} renderCell={dash} maxHeight={maxHeight} />
    : <div className="empty-state" style={{ padding: 20 }}>No results.</div>
}

/** Cytogenetic ideogram drawn on a canvas, with hover tooltip. */
function Ideogram({ data, s }) {
  const canvasRef = useRef(null)
  const hits = useRef([])
  const [tip, setTip] = useState(null)

  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas || !data) return
    const ctx = canvas.getContext('2d')
    const chroms = Object.values(data)
    hits.current = []
    const scheme = CYTO_SCHEMES[s.scheme] || CYTO_SCHEMES.classic
    const z = s.zoom
    const w = s.width * z, gap = s.spacing * z, ppc = s.pxPerCm * z, fs = s.fontSize
    const left = s.showScale ? 70 : 30, top = 50, bottom = 60
    const maxPos = Math.max(0, ...chroms.map((c) => c.max_pos || 0))
    const mapH = Math.max(200, Math.round(maxPos * ppc))
    canvas.width = Math.max(left + chroms.length * (w + gap) + 40, 400)
    canvas.height = Math.max(top + mapH + bottom, 300)
    ctx.fillStyle = scheme.bg
    ctx.fillRect(0, 0, canvas.width, canvas.height)

    chroms.forEach((chr, ci) => {
      const x0 = left + ci * (w + gap), y0 = top
      const chrH = Math.round((chr.max_pos || 0) * ppc)
      const mk = chr.markers || []
      ctx.fillStyle = scheme.miss
      ctx.fillRect(x0, y0, w, chrH)
      mk.forEach((m, i) => {
        const yT = y0 + Math.round(m.position * ppc)
        const yB = i < mk.length - 1 ? y0 + Math.round(mk[i + 1].position * ppc) : y0 + chrH
        const h = Math.max(yB - yT, 1)
        let dom = 'miss', best = 0
        if (m.allele_a > best) { best = m.allele_a; dom = 'a' }
        if (m.allele_b > best) { best = m.allele_b; dom = 'b' }
        if (m.allele_ab > best) { best = m.allele_ab; dom = 'ab' }
        ctx.fillStyle = scheme[dom] || scheme.miss
        ctx.fillRect(x0, yT, w, h)
        if (s.showNames && h >= fs + 3) {
          ctx.fillStyle = dom === 'a' ? '#FFFFFF' : '#0F172A'
          ctx.font = `bold ${fs - 1}px monospace`
          ctx.textAlign = 'left'
          ctx.fillText(String(m.marker).substring(0, 14), x0 + 4, yT + h / 2 + (fs - 1) / 3)
        }
        hits.current.push({ x: x0, y: yT, w, h, m, chr: chr.chromosome })
      })
      ctx.strokeStyle = 'rgba(149,213,178,0.85)'
      ctx.lineWidth = Math.max(1.2, z * 1.2)
      ctx.strokeRect(x0, y0, w, chrH)
      ctx.textAlign = 'center'
      ctx.fillStyle = '#b7e4c7'; ctx.font = `bold ${fs}px sans-serif`
      ctx.fillText('Chr ' + chr.chromosome, x0 + w / 2, y0 + chrH + fs + 10)
      ctx.fillStyle = '#95d5b2'; ctx.font = `${fs - 1}px sans-serif`
      ctx.fillText('n=' + chr.n_markers, x0 + w / 2, y0 + chrH + fs * 2 + 14)
      ctx.fillStyle = '#CBD5E1'
      ctx.fillText((chr.max_pos || 0).toFixed(1) + ' cM', x0 + w / 2, y0 - 10)
    })

    if (s.showScale && maxPos > 0) {
      let step = 10
      for (const st of [5, 10, 20, 25, 50, 100, 150, 200]) if (mapH / (maxPos / st) >= 22) { step = st; break }
      const sx = left - 14
      ctx.strokeStyle = '#52b788'; ctx.lineWidth = 1.5; ctx.fillStyle = '#CBD5E1'
      ctx.font = `${fs - 1}px monospace`; ctx.textAlign = 'right'
      ctx.beginPath(); ctx.moveTo(sx, top); ctx.lineTo(sx, top + mapH); ctx.stroke()
      for (let p = 0; p <= maxPos; p += step) {
        const ty = top + Math.round(p * ppc)
        ctx.beginPath(); ctx.moveTo(sx - 5, ty); ctx.lineTo(sx, ty); ctx.stroke()
        ctx.fillText(String(p), sx - 7, ty + 4)
      }
    }
  }, [data, s])

  const onMove = (e) => {
    const c = canvasRef.current
    const r = c.getBoundingClientRect()
    const x = (e.clientX - r.left) * (c.width / r.width), y = (e.clientY - r.top) * (c.height / r.height)
    const f = hits.current.find((h) => x >= h.x && x <= h.x + h.w && y >= h.y && y <= h.y + h.h)
    setTip(f ? { ...f, cx: e.clientX, cy: e.clientY } : null)
  }
  const pct = (v) => ((v || 0) * 100).toFixed(1) + '%'

  return (
    <>
      <div className="cyto-scroll"><canvas ref={canvasRef} onMouseMove={onMove} onMouseLeave={() => setTip(null)} /></div>
      {tip && (
        <div className="cyto-tip" style={{ left: tip.cx + 16, top: tip.cy - 12 }}>
          <b>{tip.m.marker}</b><br />Chr {tip.chr} | {Number(tip.m.position).toFixed(3)} cM<br />
          Parent A: {pct(tip.m.allele_a)}<br />Parent B: {pct(tip.m.allele_b)}<br />Het (ab): {pct(tip.m.allele_ab)}<br />Missing: {pct(tip.m.allele_miss)}
        </div>
      )}
    </>
  )
}

export default function LinkageTool() {
  const { bridge } = useBridge('linkageBridge', () => mockGeneric(['errorOccurred', 'cytoDataReady']))
  const tool = useTool(bridge, 'run_linkage')
  const [tab, setTab] = useState('upload')
  const [geno, setGeno] = useState('')
  const [pop, setPop] = useState('F2')
  const [upload, setUpload] = useState(null)
  const [uploadLog, setUploadLog] = useState('No file uploaded yet.')
  const [p, setPr] = useState({
    missing_threshold: 0.3, recombination_threshold: 0.45, enable_chi_square: false, chi_square_alpha: 0.05,
    use_bonferroni: false, max_total_cm: 1500, scaling_method: 'proportional',
  })
  const setP = (k) => (v) => setPr((x) => ({ ...x, [k]: v }))
  const [chrs, setChrs] = useState([])
  const [mapSel, setMapSel] = useState([])
  const [denSel, setDenSel] = useState([])
  const [heatSel, setHeatSel] = useState([])
  const [cytoSel, setCytoSel] = useState([])
  const [mp, setMp] = useState({ theme: 'dark', color_palette: 'Set1', show_labels: false, label_size: 2.5, chrom_width: 0.4, point_size: 2.5, chrom_spacing: 2.5, plot_height: 700 })
  const setM = (k) => (v) => setMp((x) => ({ ...x, [k]: v }))
  const [bin, setBin] = useState(10)
  const [topN, setTopN] = useState(60)
  const [plots, setPlots] = useState({})
  const [busy, setBusy] = useState({})
  const [cyto, setCyto] = useState(null)
  const [cs, setCs] = useState({ zoom: 1, width: 28, spacing: 60, pxPerCm: 5, fontSize: 11, scheme: 'classic', showScale: true, showNames: false })
  const setC = (k) => (v) => setCs((x) => ({ ...x, [k]: v }))

  tool.onResult((r) => {
    const c = (r.chromosomes || []).map(String)
    setChrs(c); setMapSel(c); setDenSel(c); setHeatSel(c); setCytoSel(c)
    setPlots({ map: r.map_plot, density: r.density_plot, heatmap: r.heatmap_plot })
    setTab('results')
  })
  useSignal(bridge, 'errorOccurred', (m) => tool.fail('Error: ' + m))
  useSignal(bridge, 'cytoDataReady', (json) => {
    try {
      const d = JSON.parse(json)
      if (!d || d.status === 'error') return tool.setStatus('Cytogenetic data error: ' + ((d && d.message) || ''))
      setCyto(d.chromosomes); setCs((x) => ({ ...x, zoom: 1 }))
    } catch (e) { console.error('[LINKAGE] cyto', e) }
  })
  const plotHandler = (key, field) => (r) => {
    setBusy((b) => ({ ...b, [key]: false }))
    if (!r) return
    if (r.status === 'error') return setPlots((x) => ({ ...x, [key + 'Err']: r.message || 'Could not render plot.' }))
    if (r[field]) setPlots((x) => ({ ...x, [key]: r[field], [key + 'Err']: null }))
  }
  useGlobals({
    showMapPlot: plotHandler('map', 'map_plot'),
    showDensityPlot: plotHandler('density', 'density_plot'),
    showHeatmapPlot: plotHandler('heatmap', 'heatmap_plot'),
  })

  const browse = async () => {
    const path = await tool.browse('browse_genotype')
    if (!path) return
    setGeno(path); setUploadLog('Loading file…')
    const pv = await callJson(bridge, 'get_upload_preview', path)
    if (!pv || pv.status !== 'success') { setUpload(null); return setUploadLog('Error: ' + ((pv && pv.message) || 'Could not read file.')) }
    setUpload(pv); setUploadLog(pv.log || 'File loaded.')
    tool.setStatus('Ready. Configure parameters and click Run Analysis.')
  }
  const run = () => {
    if (!upload) return tool.fail('Upload a genotype file first.')
    tool.run({ population_type: pop, ...p })
  }
  const render = (key, slot, cfg) => { setBusy((b) => ({ ...b, [key]: true })); bridge[slot](JSON.stringify(cfg)) }
  const dl = (t) => bridge && bridge.download_result(t)
  const saveCyto = () => {
    const c = document.querySelector('.cyto-scroll canvas')
    if (!c) return
    const a = document.createElement('a')
    a.download = `cytogenetic_map_${new Date().toISOString().slice(0, 10)}.png`
    a.href = c.toDataURL('image/png', 1.0)
    a.click()
  }

  const r = tool.result
  const noRun = <EmptyState icon={ChartLine}>Run analysis first.</EmptyState>
  const plotOr = (key, label) => plots[key + 'Err'] ? <div className="info-line warn">Error: {plots[key + 'Err']}</div>
    : plots[key] ? <ImageGallery single images={[{ src: plots[key], label }]} /> : <div className="empty-state" style={{ padding: 30 }}>No plot yet.</div>
  const scheme = CYTO_SCHEMES[cs.scheme]

  return (
    <ToolPage
      id="link" icon={Dna} tabs={TABS} tab={tab} setTab={setTab}
      title="Linkage Map Creator"
      subtitle="Genetic linkage mapping and cytogenetic visualization · Kosambi cM"
      run={tab === 'parameters' ? { label: 'Run Analysis', runningLabel: 'Building linkage map…', state: tool.state, disabled: !upload,
        status: !upload && tool.state === 'idle' ? 'Upload a genotype file first.' : tool.status, onRun: run } : undefined}
    >
      {tab === 'upload' && (
        <TabPane id="upload" key="upload">
          <div className="grid-2" style={{ gridTemplateColumns: 'minmax(300px, 1fr) minmax(300px, 2fr)', alignItems: 'start' }}>
            <Card title="Upload Genotype Data" icon={FolderOpen}>
              <div className="field">
                <label className="field-title">Genotype Excel File (.xlsx)</label>
                <FilePicker value={geno} placeholder="Select genotype file…" onBrowse={browse} />
              </div>
              <SelectField label="Population Type" value={pop} onChange={setPop} options={['F2', 'RIL']} />
              <Note>
                <b>Features:</b> high-resolution linkage map plot · recombination fraction heatmap · marker density distribution ·
                cytogenetic ideogram with hover tooltips · allele-frequency banding per chromosome · full Excel export
              </Note>
            </Card>
            <Card title="Upload Status" icon={ClipboardCheck} delay={0.06}><LogBox>{uploadLog}</LogBox></Card>
          </div>
          {upload && (
            <div className="grid-2" style={{ gridTemplateColumns: 'minmax(300px, 2fr) minmax(240px, 1fr)', alignItems: 'start', marginTop: 18 }}>
              <Card title="Data Preview (first 10 rows)" icon={Table2}><Tbl rows={upload.preview_rows} columns={upload.preview_columns} maxHeight={340} /></Card>
              <Card title="Chromosome Assignments" icon={Dna}><Tbl rows={upload.chr_summary} columns={['Chromosome', 'N_Markers']} maxHeight={340} /></Card>
            </div>
          )}
        </TabPane>
      )}
      {tab === 'parameters' && (
        <TabPane id="parameters" key="parameters">
          <div className="grid-2" style={{ gridTemplateColumns: 'repeat(auto-fit, minmax(340px, 1fr))', alignItems: 'start' }}>
            <Card title="Analysis Parameters" icon={Settings2}>
              <Slider label="Missing Data Threshold" value={p.missing_threshold} onChange={setP('missing_threshold')} min={0} max={1} step={0.05} format={(v) => v.toFixed(2)} />
              <Slider label="Max Recombination Fraction (r)" value={p.recombination_threshold} onChange={setP('recombination_threshold')} min={0} max={0.5} step={0.01} format={(v) => v.toFixed(2)} />
              <div className="plot-title" style={{ marginTop: 18 }}>Chi-Square Filter</div>
              <Toggle label="Enable segregation distortion filter" checked={p.enable_chi_square} onChange={setP('enable_chi_square')} />
              {p.enable_chi_square ? (
                <>
                  <TextField label="Significance level (α)" type="number" min={0.001} max={0.2} step={0.005} value={p.chi_square_alpha} onChange={setP('chi_square_alpha')} />
                  <Toggle label="Apply Bonferroni correction" checked={p.use_bonferroni} onChange={setP('use_bonferroni')} />
                  <div className="info-line warn">Chi-square filtering removes markers with skewed allele ratios.</div>
                </>
              ) : <div className="field-hint">Chi-square filter is OFF.</div>}
              <div className="plot-title" style={{ marginTop: 18 }}>Map Scaling</div>
              <TextField label="Max Total Map Length (cM)" type="number" min={100} max={10000} step={100} value={p.max_total_cm} onChange={setP('max_total_cm')} />
              <SelectField label="Scaling Method" value={p.scaling_method} onChange={setP('scaling_method')}
                options={[{ value: 'proportional', label: 'Proportional' }, { value: 'uniform', label: 'Uniform' }]} />
            </Card>
            <Card title="Parameter Guide" icon={BookOpen} delay={0.06}>
              <ul className="plain-list">
                <li><b>Missing Data Threshold:</b> markers with more than this fraction missing are excluded.</li>
                <li><b>Max Recombination Fraction:</b> marker pairs with r ≥ threshold get a 50 cM gap.</li>
                <li><b>Monomorphic Removal:</b> automatic — fixed-allele markers are excluded.</li>
                <li><b>Chi-Square Filter:</b> optional — removes segregation-distorted markers.</li>
                <li><b>Map Scaling:</b> proportional scales all chromosomes uniformly; uniform equalizes chromosome lengths.</li>
              </ul>
              <div className="plot-title" style={{ marginTop: 14 }}>Scaling Summary</div>
              <LogBox>{(r && r.scaling_summary) || 'Run analysis to see scaling summary.'}</LogBox>
              <div className="grid-2" style={{ marginTop: 12 }}>
                <div><div className="plot-title">Before Scaling</div><LogBox>{(r && r.before_scaling_stats) || '—'}</LogBox></div>
                <div><div className="plot-title">After Scaling</div><LogBox>{(r && r.after_scaling_stats) || '—'}</LogBox></div>
              </div>
            </Card>
          </div>
        </TabPane>
      )}
      {tab === 'qc' && (
        <TabPane id="qc" key="qc">
          <Card title="Quality Control Report" icon={ClipboardCheck}><LogBox>{(r && r.qc_report) || 'Run analysis to see the QC report.'}</LogBox></Card>
          <Card title="Removed Markers — Full Detail" icon={Trash2} delay={0.06}>
            {r ? <Tbl rows={r.removal_log} columns={['Marker', 'Chromosome', 'Reason', 'Detail']} /> : <div className="empty-state" style={{ padding: 20 }}>No removal log yet.</div>}
            <button className="btn btn-secondary" style={{ marginTop: 14 }} onClick={() => dl('removal_log')} disabled={!r}><FileSpreadsheet size={15} /> Download Removal Log (.xlsx)</button>
          </Card>
        </TabPane>
      )}
      {tab === 'results' && (
        <TabPane id="results" key="results">
          {!r ? <EmptyState icon={Table2}>Upload data and run analysis to see results.</EmptyState> : (
            <>
              <Card title="Downloads" icon={Download}>
                <div className="downloads" style={{ margin: 0 }}>
                  <button className="btn btn-secondary" onClick={() => dl('full_results')}><FileSpreadsheet size={15} /> Full Excel Results</button>
                  <button className="btn btn-secondary" onClick={() => dl('all_markers')}><FileSpreadsheet size={15} /> All Markers + Positions</button>
                  <button className="btn btn-secondary" onClick={() => dl('retained_markers')}><FileSpreadsheet size={15} /> Retained Markers List</button>
                </div>
              </Card>
              <Card title="Retained Markers" icon={Table2}><Tbl rows={r.retained_markers} columns={r.retained_markers_columns} /></Card>
              <Card title="All Markers — Status + Positions" icon={Table2}><Tbl rows={r.all_markers} columns={r.all_markers_columns} /></Card>
              <Card title="Linkage Map Table" icon={MapIcon}><Tbl rows={r.linkage_map} columns={r.linkage_map_columns} /></Card>
              <Card title="Summary Statistics" icon={ChartLine}><LogBox>{r.summary_stats || '—'}</LogBox></Card>
            </>
          )}
        </TabPane>
      )}
      {tab === 'mapplot' && (
        <TabPane id="mapplot" key="mapplot">
          {!r ? noRun : (
            <div className="grid-2" style={{ gridTemplateColumns: 'minmax(280px, 1fr) minmax(320px, 2.2fr)', alignItems: 'start' }}>
              <Card title="Plot Options" icon={Settings2}>
                <div className="plot-title">Chromosomes</div>
                <ChipGroup options={chrs} value={mapSel} onChange={setMapSel} />
                <SelectField label="Plot Theme" value={mp.theme} onChange={setM('theme')}
                  options={[{ value: 'dark', label: 'Dark' }, { value: 'blue', label: 'Blue/White' }, { value: 'minimal_white', label: 'Minimal White' }]} />
                <SelectField label="Color Palette" value={mp.color_palette} onChange={setM('color_palette')} options={PALETTES} />
                <div style={{ marginTop: 8 }}><Toggle label="Show Marker Labels" checked={mp.show_labels} onChange={setM('show_labels')} /></div>
                <Slider label="Label Size" value={mp.label_size} onChange={setM('label_size')} min={1.5} max={6} step={0.5} format={(v) => v.toFixed(1)} />
                <Slider label="Chromosome Width" value={mp.chrom_width} onChange={setM('chrom_width')} min={0.1} max={1.2} step={0.05} format={(v) => v.toFixed(2)} />
                <Slider label="Marker Point Size" value={mp.point_size} onChange={setM('point_size')} min={0.5} max={6} step={0.25} format={(v) => v.toFixed(2)} />
                <Slider label="LG Spacing" value={mp.chrom_spacing} onChange={setM('chrom_spacing')} min={1} max={6} step={0.25} format={(v) => v.toFixed(2)} />
                <Slider label="Plot Height (px)" value={mp.plot_height} onChange={setM('plot_height')} min={400} max={2000} step={50} />
                <button className="btn btn-secondary" style={{ width: '100%', marginTop: 16 }} disabled={busy.map}
                  onClick={() => render('map', 'render_map_plot', { chromosomes: mapSel, ...mp })}>
                  <RefreshCw size={15} className={busy.map ? 'spin' : ''} /> {busy.map ? 'Rendering…' : 'Update Plot'}
                </button>
                <button className="btn btn-secondary" style={{ width: '100%', marginTop: 8 }} onClick={() => dl('map_plot')}><Download size={15} /> Download PNG (300 dpi)</button>
              </Card>
              <Card title="Genetic Linkage Map" icon={MapIcon} delay={0.06}>{plotOr('map', 'Genetic Linkage Map')}</Card>
            </div>
          )}
        </TabPane>
      )}
      {tab === 'density' && (
        <TabPane id="density" key="density">
          {!r ? noRun : (
            <>
              <Card title="Marker Density per Chromosome" icon={ChartArea}>
                <div className="grid-2" style={{ alignItems: 'end' }}>
                  <div><div className="plot-title">Chromosomes</div><ChipGroup options={chrs} value={denSel} onChange={setDenSel} /></div>
                  <div>
                    <Slider label="Bin Width (cM)" value={bin} onChange={setBin} min={2} max={30} step={2} />
                    <div className="downloads" style={{ marginTop: 14, marginBottom: 0 }}>
                      <button className="btn btn-secondary" disabled={busy.density} onClick={() => render('density', 'render_density_plot', { chromosomes: denSel, bin_width: bin })}>
                        <RefreshCw size={15} className={busy.density ? 'spin' : ''} /> Update</button>
                      <button className="btn btn-secondary" onClick={() => dl('density_plot')}><Download size={15} /> Download</button>
                    </div>
                  </div>
                </div>
                <div style={{ marginTop: 16 }}>{plotOr('density', 'Marker Density')}</div>
              </Card>
              <Card title="Adjacent Marker Recombination Heatmap" icon={Thermometer} delay={0.06}>
                <div className="grid-2" style={{ alignItems: 'end' }}>
                  <div><div className="plot-title">Chromosomes</div><ChipGroup options={chrs} value={heatSel} onChange={setHeatSel} /></div>
                  <div>
                    <Slider label="Max Pairs to Show" value={topN} onChange={setTopN} min={10} max={200} step={10} />
                    <div className="downloads" style={{ marginTop: 14, marginBottom: 0 }}>
                      <button className="btn btn-secondary" disabled={busy.heatmap} onClick={() => render('heatmap', 'render_heatmap_plot', { chromosomes: heatSel, top_n: topN })}>
                        <RefreshCw size={15} className={busy.heatmap ? 'spin' : ''} /> Update</button>
                      <button className="btn btn-secondary" onClick={() => dl('heatmap_plot')}><Download size={15} /> Download</button>
                    </div>
                  </div>
                </div>
                <div style={{ marginTop: 16 }}>{plotOr('heatmap', 'Recombination Heatmap')}</div>
              </Card>
            </>
          )}
        </TabPane>
      )}
      {tab === 'cyto' && (
        <TabPane id="cyto" key="cyto">
          {!r ? noRun : (
            <Card title="High-Resolution Cytogenetic Ideogram" icon={Dna}>
              <div className="grid-2" style={{ alignItems: 'end' }}>
                <div>
                  <div className="plot-title">Chromosomes</div>
                  <ChipGroup options={chrs} value={cytoSel} onChange={setCytoSel} />
                  <button className="btn btn-secondary" style={{ marginTop: 12 }} onClick={() => bridge.get_cyto_data(JSON.stringify({ chromosomes: cytoSel }))}>
                    <RefreshCw size={15} /> Load Chromosomes</button>
                </div>
                <div>
                  <Slider label="Zoom" value={cs.zoom} onChange={setC('zoom')} min={0.5} max={3} step={0.1} format={(v) => v.toFixed(1) + '×'} />
                  <button className="btn btn-secondary" style={{ marginTop: 12 }} onClick={saveCyto} disabled={!cyto}><Camera size={15} /> Save as PNG</button>
                </div>
              </div>
              <div style={{ marginTop: 16 }}>
                {cyto ? <Ideogram data={cyto} s={cs} /> : <div className="empty-state" style={{ padding: 30 }}>Click “Load Chromosomes” to draw the ideogram.</div>}
              </div>
              <div className="cyto-legend">
                <span><i style={{ background: scheme.a }} />Parent A (aa)</span>
                <span><i style={{ background: scheme.b }} />Parent B (bb)</span>
                <span><i style={{ background: scheme.ab }} />Heterozygous (ab)</span>
                <span><i style={{ background: 'rgba(148,163,184,.4)', border: '1px dashed #94A3B8' }} />Missing / unknown</span>
              </div>
              <div className="grid-2" style={{ gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', marginTop: 10 }}>
                <div><Slider label="Chromosome Width" value={cs.width} onChange={setC('width')} min={10} max={80} step={2} /></div>
                <div><Slider label="Spacing" value={cs.spacing} onChange={setC('spacing')} min={20} max={150} step={5} /></div>
                <div><Slider label="px per cM" value={cs.pxPerCm} onChange={setC('pxPerCm')} min={1} max={25} step={0.5} /></div>
                <div><Slider label="Font Size" value={cs.fontSize} onChange={setC('fontSize')} min={8} max={20} /></div>
              </div>
              <div className="grid-2" style={{ gridTemplateColumns: 'repeat(auto-fit, minmax(220px, 1fr))', marginTop: 8, alignItems: 'end' }}>
                <SelectField label="Color Scheme" value={cs.scheme} onChange={setC('scheme')} options={[
                  { value: 'classic', label: 'Dark Ideogram' }, { value: 'blue_orange', label: 'Genomics Blue' },
                  { value: 'green_purple', label: 'Forest' }, { value: 'red_blue', label: 'Red/Blue' }, { value: 'pastel', label: 'Pastel' }]} />
                <Toggle label="Show cM scale bar" checked={cs.showScale} onChange={setC('showScale')} />
                <Toggle label="Show marker names (≤200 markers)" checked={cs.showNames} onChange={setC('showNames')} />
              </div>
            </Card>
          )}
        </TabPane>
      )}
      {tab === 'help' && (
        <TabPane id="help" key="help">
          <SectionHeader title="Application Guide" />
          <Card title="Required Input Format" icon={BookOpen}>
            <ul className="plain-list">
              <li>Excel (.xlsx): Column 1 = "Marker", Column 2 = "Chromosome"</li>
              <li>Genotype columns 3+: individual sample genotype calls</li>
              <li>Genotype codes are auto-converted to a / b / ab / − (missing)</li>
              <li>Marker order in the input file is preserved during map construction</li>
            </ul>
          </Card>
          <Card title="Workflow" icon={BookOpen} delay={0.06}>
            <ol className="plain-list">
              <li>Upload a genotype Excel file and confirm the population type (F2 or RIL) in the Upload tab.</li>
              <li>Set thresholds and optional filters in the Parameters tab, then click Run Analysis.</li>
              <li>Review excluded markers in the QC Report tab.</li>
              <li>Browse tables and download Excel exports in the Results tab.</li>
              <li>Customize and export the linkage map figure in the Map Plot tab.</li>
              <li>Inspect marker density and adjacent-marker recombination in the Density & RF tab.</li>
              <li>Explore the interactive cytogenetic ideogram in the Cytogenetic tab.</li>
            </ol>
          </Card>
        </TabPane>
      )}
    </ToolPage>
  )
}
