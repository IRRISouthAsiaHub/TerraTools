import { useState } from 'react'
import { STATUSES, STATUS_COLOURS, statusOf } from '../../lib/favAllele.js'

export const SERIES_COLOURS = ['#2d6a4f', '#e9a23b', '#3f7cac', '#d1603d', '#7b5ea7', '#52b788', '#b5838d', '#264653']
export const seriesColour = (i) => SERIES_COLOURS[i % SERIES_COLOURS.length]

const W = 760
const TXT = '#3f5a4c'

/** Grouped horizontal bars: one group per trait, one bar per selected line. */
export function TraitBars({ traits, series, metric }) {
  const [hover, setHover] = useState(null)
  if (!traits.length) return <div className="fa-empty">No SNPs match the filters.</div>
  const left = 200, right = 46, barH = 13, gap = 3, groupPad = 14
  const groupH = series.length * (barH + gap) + groupPad
  const H = traits.length * groupH + 36
  const max = Math.max(1, ...series.flatMap((s) => s.values))
  const x = (v) => left + ((W - left - right) * v) / max
  const ticks = [...new Set([0, 0.25, 0.5, 0.75, 1].map((f) => Math.round(max * f)))]
  return (
    <svg viewBox={`0 0 ${W} ${H}`} className="fa-chart" role="img" aria-label={`${metric} SNP count by trait`}>
      {ticks.map((t) => (
        <g key={t}>
          <line x1={x(t)} x2={x(t)} y1={6} y2={H - 26} stroke="#dbe9e1" />
          <text x={x(t)} y={H - 10} textAnchor="middle" fontSize="11" fill={TXT}>{t}</text>
        </g>
      ))}
      {traits.map((t, ti) => {
        const y0 = 8 + ti * groupH
        return (
          <g key={t}>
            <text x={left - 10} y={y0 + (series.length * (barH + gap)) / 2 + 4} textAnchor="end" fontSize="12" fill="#1a3a2a" fontWeight="600">{t}</text>
            {series.map((s, si) => {
              const v = s.values[ti]
              const y = y0 + si * (barH + gap)
              const key = `${ti}-${si}`
              return (
                <g key={s.id} onMouseEnter={() => setHover(key)} onMouseLeave={() => setHover(null)} opacity={hover && hover !== key ? 0.55 : 1}>
                  <title>{`${s.name} · ${t}: ${v} ${metric.toLowerCase()} SNP${v === 1 ? '' : 's'}`}</title>
                  <rect x={left} y={y} width={Math.max(0, x(v) - left)} height={barH} rx="3" fill={s.colour} />
                  <text x={x(v) + 5} y={y + barH - 2.5} fontSize="10.5" fill={TXT}>{v}</text>
                </g>
              )
            })}
          </g>
        )
      })}
    </svg>
  )
}

/** Stacked horizontal bars of the status counts of the top genes of the shown line. Click a bar for the gene. */
export function GeneStack({ genes, onGene }) {
  const [hover, setHover] = useState(null)
  if (!genes.length) return <div className="fa-empty">No genes match the filters.</div>
  const left = 150, right = 40, rowH = 21
  const H = genes.length * rowH + 40
  const max = Math.max(1, ...genes.map((g) => g.total))
  const x = (v) => left + ((W - left - right) * v) / max
  return (
    <svg viewBox={`0 0 ${W} ${H}`} className="fa-chart" role="img" aria-label="Allele status by gene">
      {Array.from({ length: max + 1 }, (_, i) => i).filter((i) => max <= 8 || i % Math.ceil(max / 8) === 0).map((t) => (
        <g key={t}>
          <line x1={x(t)} x2={x(t)} y1={4} y2={H - 26} stroke="#dbe9e1" />
          <text x={x(t)} y={H - 10} textAnchor="middle" fontSize="11" fill={TXT}>{t}</text>
        </g>
      ))}
      {genes.map((g, i) => {
        let acc = 0
        const y = 6 + i * rowH
        return (
          <g key={g.gene} style={{ cursor: 'pointer' }} onClick={() => onGene(g.gene)} onMouseEnter={() => setHover(g.gene)} onMouseLeave={() => setHover(null)}>
            <title>{`${g.gene} · ${g.trait}\n${g.Favourable} favourable · ${g.Unfavourable} unfavourable · ${g.Heterozygous} heterozygous · ${g.Missing} missing`}</title>
            <rect x={0} y={y - 2} width={W} height={rowH} fill={hover === g.gene ? 'rgba(82,183,136,.12)' : 'transparent'} />
            <text x={left - 8} y={y + 12} textAnchor="end" fontSize="11.5" fill="#1a3a2a" fontWeight={hover === g.gene ? 700 : 500}>
              {g.gene.length > 22 ? `${g.gene.slice(0, 21)}…` : g.gene}
            </text>
            {['Favourable', 'Heterozygous', 'Unfavourable', 'Missing'].map((st) => {
              const v = g[st]
              if (!v) return null
              const r = <rect key={st} x={x(acc)} y={y} width={x(acc + v) - x(acc)} height={rowH - 6} fill={STATUS_COLOURS[st]} stroke="#fff" strokeWidth="0.8" />
              acc += v
              return r
            })}
          </g>
        )
      })}
    </svg>
  )
}

/** Lines (rows) × SNPs of the best-typed genes (columns), coloured by status. Click a cell for its gene. */
export function Heatmap({ rows, varieties, mode, onGene }) {
  const [hover, setHover] = useState(null)
  if (!rows.length) return <div className="fa-empty">No SNPs match the filters.</div>
  return (
    <div className="fa-heat-wrap">
      <table className="fa-heat">
        <thead>
          <tr>
            <th className="first">Line</th>
            {rows.map(({ gene, snp, n, k }) => (
              <th key={snp.i} className={hover === gene ? 'hl' : ''} onClick={() => onGene(gene)} title={gene}>
                <span>{n > 1 ? `${gene} #${k + 1}` : gene}</span>
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {varieties.map((v) => (
            <tr key={v.id}>
              <td className="row-label"><b>{v.name}</b></td>
              {rows.map(({ gene, snp }) => {
                const r = statusOf(snp, v, mode)
                return (
                  <td key={snp.i} className={hover === gene ? 'hl' : ''} onMouseEnter={() => setHover(gene)} onMouseLeave={() => setHover(null)} onClick={() => onGene(gene)}>
                    <span className="cell" style={{ background: STATUS_COLOURS[r.status] }} title={`${gene} · ${v.name}\n${r.call} (${r.status})`} />
                  </td>
                )
              })}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  )
}

export function StatusLegend() {
  return (
    <div className="fa-legend">
      {STATUSES.map((s) => <span key={s}><i style={{ background: STATUS_COLOURS[s] }} />{s}</span>)}
    </div>
  )
}
