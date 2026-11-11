import { useEffect, useMemo, useRef, useState } from 'react'
import { createPortal } from 'react-dom'
import { motion, AnimatePresence } from 'framer-motion'
import { ChevronDown, Search, X, Check } from 'lucide-react'
import { traitColour, statusOf, STATUS_COLOURS, STATUSES } from '../../lib/favAllele.js'
import { seriesColour } from './Charts.jsx'

/** A drop-down with check boxes (several choices), a search box and Select all / Clear. An empty selection means "all". */
export function MultiSelect({ label, allLabel, noun, options, value, onChange, search = false, dot, clearLabel = 'Clear (all)' }) {
  const [open, setOpen] = useState(false)
  const [q, setQ] = useState('')
  const box = useRef(null)
  useEffect(() => {
    if (!open) return undefined
    const away = (e) => { if (box.current && !box.current.contains(e.target)) setOpen(false) }
    const esc = (e) => { if (e.key === 'Escape') setOpen(false) }
    document.addEventListener('mousedown', away)
    document.addEventListener('keydown', esc)
    return () => { document.removeEventListener('mousedown', away); document.removeEventListener('keydown', esc) }
  }, [open])
  const shown = useMemo(() => options.filter((o) => !q || `${o.label} ${o.hint || ''}`.toLowerCase().includes(q.toLowerCase())), [options, q])
  const toggle = (v) => onChange(value.includes(v) ? value.filter((x) => x !== v) : [...value, v])
  const text = !value.length ? allLabel : value.length === 1 ? value[0] : `${value.length} ${noun} selected`
  return (
    <div className="field fa-multi" ref={box}>
      <label className="field-title">{label}{value.length > 0 && <b className="fa-count">{value.length}</b>}</label>
      <button type="button" className={`input fa-multi-btn ${open ? 'open' : ''}`} onClick={() => setOpen((o) => !o)} aria-expanded={open}>
        <span>{text}</span><ChevronDown size={16} />
      </button>
      <AnimatePresence>
        {open && (
          <motion.div className="fa-multi-panel" initial={{ opacity: 0, y: -6 }} animate={{ opacity: 1, y: 0 }} exit={{ opacity: 0, y: -6 }} transition={{ duration: 0.15 }}>
            {search && (
              <label className="search"><Search size={14} /><input autoFocus value={q} onChange={(e) => setQ(e.target.value)} placeholder={`Search ${noun}…`} /></label>
            )}
            <div className="fa-multi-actions">
              <button type="button" className="btn btn-secondary sm" onClick={() => onChange(options.map((o) => o.value))}>Select all</button>
              <button type="button" className="btn btn-secondary sm" onClick={() => onChange([])}>{clearLabel}</button>
            </div>
            <div className="fa-multi-list">
              {shown.map((o) => (
                <label key={o.value} className={`fa-multi-item ${value.includes(o.value) ? 'on' : ''}`}>
                  <input type="checkbox" checked={value.includes(o.value)} onChange={() => toggle(o.value)} />
                  <span className="box">{value.includes(o.value) && <Check size={12} strokeWidth={3} />}</span>
                  {dot && <i className="dot" style={{ background: dot(o.value) }} />}
                  <span className="t">{o.label}</span>
                  {o.hint && <small>{o.hint}</small>}
                </label>
              ))}
              {!shown.length && <div className="fa-empty sm">No match.</div>}
            </div>
          </motion.div>
        )}
      </AnimatePresence>
    </div>
  )
}

/** Line list with search + check boxes; the chosen ones are shown as pills (click = show on the plant, × = remove). */
export function VarietyPicker({ varieties, selected, active, onToggle, onActive, onClear, onMany, max = 12 }) {
  const [q, setQ] = useState('')
  const shown = useMemo(() => {
    const n = q.trim().toLowerCase()
    return varieties.filter((v) => !n || v.name.toLowerCase().includes(n))
  }, [varieties, q])
  const byId = useMemo(() => new Map(varieties.map((v) => [v.id, v])), [varieties])
  return (
    <>
      <div className="table-tools">
        <label className="search"><Search size={15} /><input value={q} onChange={(e) => setQ(e.target.value)} placeholder="Search lines by name…" />
          {q && <button type="button" onClick={() => setQ('')} aria-label="Clear"><X size={14} /></button>}</label>
        <small>{shown.length} of {varieties.length} lines · {selected.length} selected</small>
      </div>
      <div className="fa-pills">
        {selected.length === 0 && <span className="fa-hint">Tick at least one line below.</span>}
        {selected.map((id, i) => {
          const v = byId.get(id)
          if (!v) return null
          return (
            <span key={id} className={`fa-pill ${id === active ? 'active' : ''}`} style={{ '--c': seriesColour(i) }}>
              <button type="button" className="main" onClick={() => onActive(id)} title="Show this line on the plant">
                <i />{v.name}
              </button>
              <button type="button" className="x" onClick={() => onToggle(id)} aria-label={`Remove ${v.name}`}><X size={12} /></button>
            </span>
          )
        })}
        {selected.length > 1 && <button type="button" className="btn btn-secondary sm" onClick={onClear}>Clear selection</button>}
        {q.trim() && shown.length > 0 && shown.length <= max && (
          <button type="button" className="btn btn-secondary sm" onClick={() => onMany(shown.map((v) => v.id))}>Select these {shown.length}</button>
        )}
      </div>
      <div className="fa-vlist">
        {shown.map((v) => {
          const on = selected.includes(v.id)
          return (
            <label key={v.id} className={`fa-vitem ${on ? 'on' : ''}`}>
              <input type="checkbox" checked={on} onChange={() => onToggle(v.id)} />
              <span className="box">{on && <Check size={12} strokeWidth={3} />}</span>
              <span className="t">{v.name}</span>
            </label>
          )
        })}
        {!shown.length && <div className="fa-empty sm">No line matches “{q}”.</div>}
      </div>
    </>
  )
}

/** Details of one gene: its trait, position and the status of every SNP in every selected line. */
export function GeneDialog({ gene, varieties, mode, onClose }) {
  useEffect(() => {
    if (!gene) return undefined
    const esc = (e) => { if (e.key === 'Escape') onClose() }
    document.addEventListener('keydown', esc)
    return () => document.removeEventListener('keydown', esc)
  }, [gene, onClose])
  return createPortal(
    <AnimatePresence>
      {gene && (
        <motion.div className="fa-modal-back" initial={{ opacity: 0 }} animate={{ opacity: 1 }} exit={{ opacity: 0 }} onClick={onClose}>
          <motion.div className="fa-modal" role="dialog" aria-modal="true" aria-label={gene.name}
            initial={{ opacity: 0, y: 24, scale: 0.97 }} animate={{ opacity: 1, y: 0, scale: 1 }} exit={{ opacity: 0, y: 12 }} transition={{ duration: 0.22, ease: [0.22, 1, 0.36, 1] }}
            onClick={(e) => e.stopPropagation()}>
            <button type="button" className="fa-modal-x" onClick={onClose} aria-label="Close"><X size={18} /></button>
            <h2>{gene.name}</h2>
            <div className="fa-mono fa-sub">{gene.geneId}</div>
            <div className="fa-info">
              <div><small>Trait</small><b style={{ color: traitColour(gene.trait) }}>{gene.trait}</b></div>
              <div><small>Sub-category</small><b>{gene.sub || '—'}</b></div>
              <div><small>Chromosome · position</small><b className="fa-mono">{gene.chrom != null && gene.pos != null ? `chr${gene.chrom} : ${Number(gene.pos).toLocaleString()}` : '—'}</b></div>
              <div><small>SNPs in gene</small><b>{gene.snps.length}</b></div>
            </div>
            <div className="fa-tags">
              {STATUSES.map((s) => (
                <span key={s} style={{ '--c': STATUS_COLOURS[s] }}>{s}<b>{gene.counts[s]}</b></span>
              ))}
            </div>
            <h4>Per-line SNP status</h4>
            <div className="table-wrap" style={{ maxHeight: 320 }}>
              <table className="results-table">
                <thead>
                  <tr><th>SNP</th><th>REF</th><th>ALT</th><th>Favourable</th>{varieties.map((v, i) => <th key={v.id}><i className="fa-sw" style={{ background: seriesColour(i) }} />{v.name}</th>)}</tr>
                </thead>
                <tbody>
                  {gene.snps.map((s) => (
                    <tr key={s.i}>
                      <td className="fa-mono">{s.marker.split('_').slice(0, 2).join('_')}</td>
                      <td className="fa-mono">{s.ref || '-'}</td>
                      <td className="fa-mono">{s.alt || '-'}</td>
                      <td className="fa-mono" style={{ color: STATUS_COLOURS.Favourable, fontWeight: 700 }}>{s.fav || '-'}</td>
                      {varieties.map((v) => {
                        const r = statusOf(s, v, mode)
                        return <td key={v.id}><span className="fa-chip" style={{ '--c': STATUS_COLOURS[r.status] }} title={r.status}>{r.call}</span></td>
                      })}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </motion.div>
        </motion.div>
      )}
    </AnimatePresence>,
    document.body,
  )
}
