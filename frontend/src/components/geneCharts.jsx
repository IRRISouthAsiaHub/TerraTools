/* Hand-drawn SVG charts for Gene Prioritization (the app has no chart library). */

export const CAT_COLORS = {
  association: '#2d6a4f',
  function: '#52b788',
  pathway: '#457b9d',
  ortholog: '#e9c46a',
  literature: '#b5838d',
}
export const GENETIC_COLOR = '#e76f51'
const INK = '#1a2e22'
const MUTED = '#5a7566'
const GRID = 'rgba(82,183,136,0.25)'

const clip = (s, n) => (s && s.length > n ? s.slice(0, n - 1) + '…' : s || '')
const ROW = 26

function Axis({ x0, x1, top, bottom, max = 100, step = 20 }) {
  const ticks = []
  for (let v = 0; v <= max; v += step) ticks.push(v)
  const x = (v) => x0 + ((x1 - x0) * v) / max
  return (
    <g fontSize="11" fill={MUTED}>
      {ticks.map((v) => (
        <g key={v}>
          <line x1={x(v)} x2={x(v)} y1={top} y2={bottom} stroke={GRID} />
          <text x={x(v)} y={bottom + 14} textAnchor="middle">{v}</text>
        </g>
      ))}
    </g>
  )
}

export function Legend({ items }) {
  return (
    <div style={{ display: 'flex', flexWrap: 'wrap', gap: '6px 16px', fontSize: 12.5, color: MUTED, margin: '6px 0 2px' }}>
      {items.map((it) => (
        <span key={it.label} style={{ display: 'inline-flex', alignItems: 'center', gap: 6 }}>
          <i style={{ width: 12, height: 12, borderRadius: 3, background: it.color, display: 'inline-block' }} />
          {it.label}
        </span>
      ))}
    </div>
  )
}

/** Final score (and genetic score when the run has one) of the top genes. */
export function ScoreBars({ bars, hasGenetic, onSelect }) {
  if (!bars || !bars.length) return null
  const W = 760, left = 150, right = 40, top = 8
  const per = hasGenetic ? 2 : 1
  const rowH = per * 13 + 9
  const H = top + bars.length * rowH + 28
  const x = (v) => left + ((W - left - right) * (v || 0)) / 100
  return (
    <svg viewBox={`0 0 ${W} ${H}`} width="100%" role="img" aria-label="Score of the top genes" style={{ maxWidth: W, display: 'block' }}>
      <Axis x0={left} x1={W - right} top={top} bottom={H - 28} />
      {bars.map((b, i) => {
        const y = top + i * rowH
        return (
          <g key={b.id} style={{ cursor: onSelect ? 'pointer' : undefined }} onClick={() => onSelect && onSelect(b.id)}>
            <title>{`${b.gene} (${b.id})  final ${b.final ?? 'n/a'}${hasGenetic ? `  genetic ${b.genetic ?? 'n/a'}` : ''}`}</title>
            <text x={left - 8} y={y + rowH / 2 + 1} textAnchor="end" fontSize="12" fill={INK}>{clip(b.gene, 20)}</text>
            {b.final != null
              ? <rect x={left} y={y + 2} width={Math.max(2, x(b.final) - left)} height={11} rx={3} fill={CAT_COLORS.association} />
              : <text x={left + 4} y={y + 11} fontSize="10.5" fill={MUTED}>no trait evidence available</text>}
            {b.final != null && <text x={x(b.final) + 5} y={y + 11} fontSize="10.5" fill={MUTED}>{b.final.toFixed(1)}</text>}
            {hasGenetic && b.genetic != null && (
              <>
                <rect x={left} y={y + 15} width={Math.max(2, x(b.genetic) - left)} height={11} rx={3} fill={GENETIC_COLOR} />
                <text x={x(b.genetic) + 5} y={y + 24} fontSize="10.5" fill={MUTED}>{b.genetic.toFixed(1)}</text>
              </>
            )}
          </g>
        )
      })}
    </svg>
  )
}

/** Stacked bar: how many of the final score's points each evidence category supplied. */
export function ContributionBars({ rows, weights, onSelect }) {
  if (!rows || !rows.length) return null
  const keys = Object.keys(weights)
  const W = 760, left = 150, right = 40, top = 8
  const H = top + rows.length * ROW + 28
  const x = (v) => left + ((W - left - right) * v) / 100
  return (
    <svg viewBox={`0 0 ${W} ${H}`} width="100%" role="img" aria-label="Evidence contribution to the final score" style={{ maxWidth: W, display: 'block' }}>
      <Axis x0={left} x1={W - right} top={top} bottom={H - 28} />
      {rows.map((r, i) => {
        const y = top + i * ROW
        let acc = 0
        const total = keys.reduce((s, k) => s + (r[k] || 0), 0)
        return (
          <g key={r.id} style={{ cursor: onSelect ? 'pointer' : undefined }} onClick={() => onSelect && onSelect(r.id)}>
            <text x={left - 8} y={y + 17} textAnchor="end" fontSize="12" fill={INK}>{clip(r.gene, 20)}</text>
            {keys.map((k) => {
              const v = r[k] || 0
              if (v <= 0) return null
              const rect = (
                <rect key={k} x={x(acc)} y={y + 4} width={Math.max(1.5, x(acc + v) - x(acc))} height={16} fill={CAT_COLORS[k]}>
                  <title>{`${weights[k].label}: ${v.toFixed(1)} points of ${total.toFixed(1)}`}</title>
                </rect>
              )
              acc += v
              return rect
            })}
            {total > 0 && <text x={x(total) + 5} y={y + 17} fontSize="10.5" fill={MUTED}>{total.toFixed(1)}</text>}
          </g>
        )
      })}
    </svg>
  )
}

/** Gene x pathway membership. */
export function PathwayMatrix({ data, onSelect }) {
  if (!data || !data.pathways.length || !data.genes.length) return null
  const cell = 24, left = 150, labelH = 175
  const W = left + data.pathways.length * cell + 130
  const H = labelH + data.genes.length * cell + 10
  return (
    <div style={{ overflowX: 'auto' }}>
      <svg viewBox={`0 0 ${W} ${H}`} width={W} height={H} role="img" aria-label="Genes and the pathways they belong to">
        {data.pathways.map((p, j) => {
          const cx = left + j * cell + cell / 2
          return (
            <g key={p.id}>
              <title>{p.name}</title>
              <line x1={cx} x2={cx} y1={labelH - 4} y2={H - 6} stroke={GRID} />
              <text transform={`translate(${cx + 4},${labelH - 8}) rotate(-60)`} fontSize="11" fill={INK}>{clip(p.name, 30)}</text>
            </g>
          )
        })}
        {data.genes.map((g, i) => {
          const cy = labelH + i * cell + cell / 2
          return (
            <g key={g.id} style={{ cursor: onSelect ? 'pointer' : undefined }} onClick={() => onSelect && onSelect(g.id)}>
              <line x1={left} x2={left + data.pathways.length * cell} y1={cy} y2={cy} stroke={GRID} />
              <text x={left - 8} y={cy + 4} textAnchor="end" fontSize="12" fill={INK}>{clip(g.gene, 20)}</text>
              {g.in.map((pid) => {
                const j = data.pathways.findIndex((p) => p.id === pid)
                return (
                  <circle key={pid} cx={left + j * cell + cell / 2} cy={cy} r={7} fill={CAT_COLORS.pathway}>
                    <title>{`${g.gene} — ${data.pathways[j].name}`}</title>
                  </circle>
                )
              })}
            </g>
          )
        })}
      </svg>
    </div>
  )
}
