import { useMemo, useState } from 'react'

/* Hand-drawn SVG plots of the SNP-index analysis (the app has no chart library). */

export const LEVEL_COLOR = { 90: '#f4a261', 95: '#2a6f97', 99: '#c1121f' }
const LEVELS = [90, 95, 99]
const INK = '#1a2e22'
const MUTED = '#5a7566'
const GRID = 'rgba(82,183,136,0.22)'
const POINT = ['#7a9e7e', '#3d5a40']

const niceStep = (span) => {
  const raw = span / 7
  const p = Math.pow(10, Math.floor(Math.log10(raw)))
  const f = raw / p
  return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p
}

const shortName = (n) => String(n).replace(/^(chromosome|chrom|chr)[_ ]?/i, '') || String(n)

export function LevelLegend({ extra = true, mean = true }) {
  const item = (color, label, round) => (
    <span key={label} style={{ display: 'inline-flex', alignItems: 'center', gap: 6 }}>
      <i style={{ width: round ? 12 : 18, height: round ? 12 : 3, borderRadius: round ? 3 : 0, background: color, display: 'inline-block' }} />
      {label}
    </span>
  )
  return (
    <div style={{ display: 'flex', flexWrap: 'wrap', gap: '6px 18px', fontSize: 12.5, color: MUTED, margin: '8px 0 2px' }}>
      {mean && item(INK, 'Sliding-window mean')}
      {LEVELS.map((l) => item(LEVEL_COLOR[l], `${l}% confidence limit`))}
      {extra && item('rgba(82,183,136,0.35)', 'QTL interval', true)}
    </div>
  )
}

/**
 * One track of the genome (or of one chromosome).
 *   chroms  [{ name, length, snp:{pos,high,low,delta}, win:{mid,delta,high,low,ci90u,...} }]
 *   y       'delta' | 'high' | 'low'
 *   single  true when `chroms` has one chromosome: the x axis is in Mb
 */
export function TrackPlot({ chroms, y = 'delta', qtls = [], single = false, height = 300, onPickChrom, title, showPoints = true }) {
  const W = 1000, L = 54, R = 12, T = 12, B = single ? 36 : 30
  const H = height
  const [hover, setHover] = useState(null)
  const isDelta = y === 'delta'
  const y0 = isDelta ? -1 : 0, y1 = 1
  const py = (v) => T + (H - T - B) * (1 - (v - y0) / (y1 - y0))

  const layout = useMemo(() => {
    const gap = 0.012 * chroms.reduce((a, c) => a + c.length, 0)
    const off = {}
    let x = 0
    chroms.forEach((c) => { off[c.name] = x; x += c.length + gap })
    return { off, total: Math.max(x - gap, 1) }
  }, [chroms])
  const px = (c, v) => L + ((W - L - R) * (layout.off[c] + v)) / layout.total

  const paths = useMemo(() => chroms.map((c, ci) => {
    const pts = c.snp.pos.map((p, i) => `M${px(c.name, p).toFixed(1)} ${py(c.snp[y][i]).toFixed(1)}h.01`).join('')
    const w = c.win
    const line = (arr) => w.mid.map((m, i) => `${i ? 'L' : 'M'}${px(c.name, m).toFixed(1)} ${py(arr[i]).toFixed(1)}`).join('')
    const sig = (lv) => {                                       // the stretches of the curve outside the limits
      if (!isDelta) return ''
      let d = '', open = false
      w.mid.forEach((m, i) => {
        const out = w.delta[i] > w[`ci${lv}u`][i] || w.delta[i] < w[`ci${lv}l`][i]
        if (out) { d += `${open ? 'L' : 'M'}${px(c.name, m).toFixed(1)} ${py(w.delta[i]).toFixed(1)}`; open = true } else open = false
      })
      return d
    }
    return {
      name: c.name, color: POINT[ci % 2], pts, mean: line(w[y]),
      ci: isDelta ? LEVELS.map((lv) => ({ lv, up: line(w[`ci${lv}u`]), lo: line(w[`ci${lv}l`]) })) : [],
      sig: isDelta ? LEVELS.map((lv) => ({ lv, d: sig(lv) })) : [],
    }
  }), [chroms, y, layout]) // eslint-disable-line react-hooks/exhaustive-deps

  const yTicks = isDelta ? [-1, -0.5, 0, 0.5, 1] : [0, 0.25, 0.5, 0.75, 1]
  let xTicks = []
  if (single) {
    const c = chroms[0]
    const step = niceStep(c.length / 1e6)
    for (let v = 0; v * 1e6 <= c.length; v += step) xTicks.push({ x: px(c.name, v * 1e6), label: Number.isInteger(step) ? v : v.toFixed(1) })
  } else {
    xTicks = chroms.map((c) => ({ x: px(c.name, c.length / 2), label: shortName(c.name), name: c.name }))
  }

  const onMove = (e) => {
    const r = e.currentTarget.getBoundingClientRect()
    const sx = ((e.clientX - r.left) / r.width) * W
    if (sx < L || sx > W - R) return setHover(null)
    const gx = ((sx - L) / (W - L - R)) * layout.total
    const c = chroms.find((k) => gx >= layout.off[k.name] && gx <= layout.off[k.name] + k.length)
    if (!c || !c.win.mid.length) return setHover(null)
    const pos = gx - layout.off[c.name]
    let best = 0, bd = Infinity
    c.win.mid.forEach((m, i) => { const d = Math.abs(m - pos); if (d < bd) { bd = d; best = i } })
    setHover({ x: px(c.name, c.win.mid[best]), c, i: best })
  }

  const tip = hover && (() => {
    const w = hover.c.win, i = hover.i
    const v = w[y][i]
    return `${hover.c.name} · ${(w.mid[i] / 1e6).toFixed(2)} Mb · ${y === 'delta' ? 'mean Δ' : 'mean SNP-index'} ${v.toFixed(3)}`
      + (isDelta ? ` · 95% limit ${w.ci95l[i].toFixed(2)} … ${w.ci95u[i].toFixed(2)} · 99% ${w.ci99l[i].toFixed(2)} … ${w.ci99u[i].toFixed(2)}` : '')
      + ` · ${w.n[i]} SNPs`
  })()

  return (
    <div style={{ position: 'relative' }}>
      <svg viewBox={`0 0 ${W} ${H}`} width="100%" role="img" aria-label={title || 'SNP-index plot'} onMouseMove={onMove} onMouseLeave={() => setHover(null)}
        style={{ display: 'block' }}>
        {yTicks.map((v) => (
          <g key={v} fontSize="11" fill={MUTED}>
            <line x1={L} x2={W - R} y1={py(v)} y2={py(v)} stroke={v === 0 && isDelta ? '#8aa595' : GRID} />
            <text x={L - 7} y={py(v) + 4} textAnchor="end">{v}</text>
          </g>
        ))}
        <text x={14} y={(T + H - B) / 2} fontSize="12" fill={MUTED} textAnchor="middle" transform={`rotate(-90 14 ${(T + H - B) / 2})`}>
          {isDelta ? 'Δ SNP-index' : 'SNP-index'}
        </text>
        {!single && chroms.slice(1).map((c) => <line key={c.name} x1={px(c.name, -layout.total * 0.0055)} x2={px(c.name, -layout.total * 0.0055)} y1={T} y2={H - B} stroke="#d5e3da" />)}
        {qtls.map((q) => {
          if (layout.off[q.Chromosome] === undefined) return null
          const x0 = px(q.Chromosome, q.Start), x1 = px(q.Chromosome, q.End)
          return <rect key={q.QTL_ID} x={x0} y={T} width={Math.max(x1 - x0, 2)} height={H - T - B} fill="rgba(82,183,136,0.2)" />
        })}
        {showPoints && paths.map((p) => <path key={`p${p.name}`} d={p.pts} stroke={p.color} strokeOpacity="0.42" strokeWidth="2.3" strokeLinecap="round" fill="none" />)}
        {paths.map((p) => p.ci.map((c) => (
          <g key={`ci${p.name}${c.lv}`}>
            <path d={c.up} stroke={LEVEL_COLOR[c.lv]} strokeWidth="1.1" fill="none" />
            <path d={c.lo} stroke={LEVEL_COLOR[c.lv]} strokeWidth="1.1" fill="none" />
          </g>
        )))}
        {paths.map((p) => <path key={`m${p.name}`} d={p.mean} stroke={INK} strokeWidth="1.6" fill="none" />)}
        {paths.map((p) => p.sig.map((s) => s.d && <path key={`s${p.name}${s.lv}`} d={s.d} stroke={LEVEL_COLOR[s.lv]} strokeWidth="3.2" fill="none" strokeLinecap="butt" />))}
        {qtls.map((q) => layout.off[q.Chromosome] !== undefined && (
          <path key={`k${q.QTL_ID}`} fill={INK}
            d={q.Peak_Delta_SNP_index > 0
              ? `M${px(q.Chromosome, q.Peak_Position) - 5} ${py(1) + 2}l10 0l-5 9z`
              : `M${px(q.Chromosome, q.Peak_Position) - 5} ${py(-1) - 2}l10 0l-5 -9z`} />
        ))}
        {xTicks.map((t, i) => (
          <text key={i} x={t.x} y={H - B + 16} fontSize="11" fill={MUTED} textAnchor="middle"
            style={!single && onPickChrom ? { cursor: 'pointer', textDecoration: 'underline' } : undefined}
            onClick={!single && onPickChrom ? () => onPickChrom(t.name) : undefined}>{t.label}</text>
        ))}
        <text x={(L + W - R) / 2} y={H - 3} fontSize="11.5" fill={MUTED} textAnchor="middle">{single ? `Position on ${chroms[0].name} (Mb)` : 'Chromosome (click a number to open it)'}</text>
        {hover && <line x1={hover.x} x2={hover.x} y1={T} y2={H - B} stroke={INK} strokeOpacity="0.45" strokeDasharray="3 3" />}
      </svg>
      <div style={{ minHeight: 20, fontSize: 12.5, color: MUTED }}>{tip || ' '}</div>
    </div>
  )
}

/** Null distribution at the typical read depth with the limits, and the limits by read depth. */
export function NullPlot({ nul }) {
  if (!nul) return null
  const W = 480, H = 270, L = 48, R = 12, T = 14, B = 38
  const lo = nul.edges[0], hi = nul.edges[nul.edges.length - 1]
  const mx = Math.max(...nul.counts, 1)
  const x = (v) => L + ((W - L - R) * (v - lo)) / (hi - lo)
  const yy = (v) => T + (H - T - B) * (1 - v / mx)
  return (
    <svg viewBox={`0 0 ${W} ${H}`} width="100%" role="img" aria-label="Simulated null distribution" style={{ maxWidth: W, display: 'block' }}>
      {nul.counts.map((c, i) => (
        <rect key={i} x={x(nul.edges[i]) + 0.5} y={yy(c)} width={Math.max(x(nul.edges[i + 1]) - x(nul.edges[i]) - 1, 1)} height={H - B - yy(c)} fill="#b7d9c4" />
      ))}
      {LEVELS.map((lv) => nul.limits[lv].map((v, k) => <line key={`${lv}${k}`} x1={x(v)} x2={x(v)} y1={T} y2={H - B} stroke={LEVEL_COLOR[lv]} strokeWidth="1.8" />))}
      <line x1={L} x2={W - R} y1={H - B} y2={H - B} stroke="#8aa595" />
      {[-0.5, 0, 0.5].filter((v) => v > lo && v < hi).map((v) => (
        <text key={v} x={x(v)} y={H - B + 15} fontSize="11" fill={MUTED} textAnchor="middle">{v}</text>
      ))}
      <text x={(L + W - R) / 2} y={H - 6} fontSize="11.5" fill={MUTED} textAnchor="middle">Simulated Δ SNP-index (no QTL), depth {nul.depth_high} / {nul.depth_low}</text>
      <text x={13} y={(T + H - B) / 2} fontSize="11.5" fill={MUTED} textAnchor="middle" transform={`rotate(-90 13 ${(T + H - B) / 2})`}>Replicates</text>
    </svg>
  )
}

export function DepthLimitPlot({ data }) {
  if (!data || !data.depth.length) return null
  const W = 480, H = 270, L = 48, R = 12, T = 14, B = 38
  const x0 = Math.min(...data.depth), x1 = Math.max(...data.depth, x0 + 1)
  const x = (v) => L + ((W - L - R) * (v - x0)) / (x1 - x0)
  const yy = (v) => T + (H - T - B) * (1 - (v + 1) / 2)
  const line = (arr) => data.depth.map((d, i) => `${i ? 'L' : 'M'}${x(d).toFixed(1)} ${yy(arr[i]).toFixed(1)}`).join('')
  return (
    <svg viewBox={`0 0 ${W} ${H}`} width="100%" role="img" aria-label="Limits by read depth" style={{ maxWidth: W, display: 'block' }}>
      {[-1, -0.5, 0, 0.5, 1].map((v) => (
        <g key={v} fontSize="11" fill={MUTED}>
          <line x1={L} x2={W - R} y1={yy(v)} y2={yy(v)} stroke={v === 0 ? '#8aa595' : GRID} />
          <text x={L - 6} y={yy(v) + 4} textAnchor="end">{v}</text>
        </g>
      ))}
      {LEVELS.map((lv) => (
        <g key={lv}>
          <path d={line(data[`ci${lv}u`])} stroke={LEVEL_COLOR[lv]} strokeWidth="1.8" fill="none" />
          <path d={line(data[`ci${lv}l`])} stroke={LEVEL_COLOR[lv]} strokeWidth="1.8" fill="none" />
        </g>
      ))}
      <text x={x0} y={H - B + 15} fontSize="11" fill={MUTED} textAnchor="middle">{x0}</text>
      <text x={x1} y={H - B + 15} fontSize="11" fill={MUTED} textAnchor="middle">{x1}</text>
      <text x={(L + W - R) / 2} y={H - 6} fontSize="11.5" fill={MUTED} textAnchor="middle">Read depth</text>
      <text x={13} y={(T + H - B) / 2} fontSize="11.5" fill={MUTED} textAnchor="middle" transform={`rotate(-90 13 ${(T + H - B) / 2})`}>Limit of Δ SNP-index</text>
    </svg>
  )
}
