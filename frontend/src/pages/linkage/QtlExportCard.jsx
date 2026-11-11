import { motion } from 'framer-motion'
import { Download, FileText, FileSpreadsheet, ArrowRightLeft, Check, TriangleAlert, Waypoints } from 'lucide-react'
import { Card, StatCard } from '../../components/kit.jsx'

/*
 * Results tab: the file that goes straight into QTL Mapping — the linkage
 * map matched (by marker ID) to the genotypes that went into Linkage Mapping,
 * transposed (rows = individuals; header rows = marker, chromosome, position).
 */

const chip = (v) => {
  const k = String(v)
  const cls = k === '0' ? 'c0' : k === '1' ? 'c1' : k === '2' ? 'c2' : 'cm'
  return <span className={`code-chip ${cls}`}>{k}</span>
}

export default function QtlExportCard({ result, onDownload }) {
  const q = result.qtl_export
  if (!q) {
    return (
      <Card title="File for QTL Mapping" icon={Waypoints}>
        <div className="info-line warn">
          <TriangleAlert size={14} style={{ verticalAlign: '-2px', marginRight: 6 }} />
          The QTL Mapping file could not be created{result.qtl_export_error ? `: ${result.qtl_export_error}` : '.'} The rest of the results are unaffected.
        </div>
      </Card>
    )
  }
  const pv = q.preview
  const reasons = Object.entries(q.unmatched_input_reasons || {})

  return (
    <motion.div initial={{ opacity: 0, y: 14 }} animate={{ opacity: 1, y: 0 }} transition={{ duration: 0.45 }}>
      <Card title="File for QTL Mapping" icon={Waypoints} className="qtl-export">
        <p className="qtl-export-intro">
          The linkage map above, matched by marker ID to the genotypes you uploaded and <b>transposed</b>: one row per individual,
          with marker, chromosome and position as the first three rows. Load the CSV straight into <b>QTL Mapping → Genotype CSV</b>.
        </p>

        <div className="summary-grid">
          <StatCard className="summary-card" value={q.markers} label="Markers (matched)" />
          <StatCard className="summary-card" value={q.individuals} label="Individuals" delay={0.05} />
          <StatCard className="summary-card" value={q.chromosomes} label="Chromosomes" delay={0.1} />
          <StatCard className="summary-card" value={q.missing_pct} label="Missing calls (%)" delay={0.15} />
        </div>

        <div className="downloads" style={{ marginBottom: 14 }}>
          <motion.button className="btn btn-primary" whileHover={{ y: -2 }} whileTap={{ scale: 0.97 }} onClick={() => onDownload('qtl_input_csv')}>
            <Download size={16} /> QTL Mapping input (CSV)
          </motion.button>
          <motion.button className="btn btn-secondary" whileHover={{ y: -2 }} whileTap={{ scale: 0.97 }} onClick={() => onDownload('qtl_input_xlsx')}>
            <FileSpreadsheet size={16} /> Excel with audit sheets
          </motion.button>
        </div>

        <div className="info-line">
          <Check size={14} style={{ verticalAlign: '-2px', marginRight: 6 }} />
          {q.markers.toLocaleString()} of {q.linkage_markers.toLocaleString()} linkage-map markers matched the uploaded genotypes
          {q.unmatched_linkage ? <> — <b>{q.unmatched_linkage}</b> map marker(s) were not found in the input (see the Excel sheet <i>Unmatched_Linkage_IDs</i>)</> : ' — none missing'}.
        </div>
        {q.unmatched_input > 0 && (
          <div className="qtl-notmapped">
            <b>{q.unmatched_input.toLocaleString()}</b> input markers are not in the file because they were not mapped:
            {reasons.map(([k, v]) => <span key={k} className="badge">{k} · {v}</span>)}
          </div>
        )}
        {q.duplicates_dropped > 0 && <div className="info-line warn">{q.duplicates_dropped} duplicate marker ID(s) in the input — the first of each was used.</div>}
        {q.notes.map((n) => <div className="info-line warn" key={n}>{n}</div>)}

        <div className="plot-title" style={{ marginTop: 18 }}>
          <ArrowRightLeft size={14} style={{ marginRight: 6, verticalAlign: '-2px' }} />Preview of the file (first {pv.markers.length} markers, first {pv.rows.length} individuals)
        </div>
        <div className="table-wrap" style={{ maxHeight: 320 }}>
          <table className="results-table qtl-preview">
            <thead>
              <tr><th>Marker</th>{pv.markers.map((m) => <th key={m} title={m}>{m}</th>)}<th className="more">…</th></tr>
              <tr className="meta"><th>Chromosome</th>{pv.chromosome.map((c, i) => <td key={i}>{c}</td>)}<td className="more">…</td></tr>
              <tr className="meta"><th>Position</th>{pv.position.map((p, i) => <td key={i}>{p}</td>)}<td className="more">…</td></tr>
            </thead>
            <tbody>
              {pv.rows.map((r) => (
                <tr key={r.id}><th>{r.id}</th>{r.codes.map((c, i) => <td key={i}>{chip(c)}</td>)}<td className="more">…</td></tr>
              ))}
            </tbody>
          </table>
        </div>
        <div className="field-hint" style={{ marginTop: 8 }}>
          Codes: 0 = homozygous P1 · 1 = heterozygous · 2 = homozygous P2 · −1 = missing.
          Chromosomes are in natural order (1, 2, … 10, 11) and markers by position.
        </div>
      </Card>
    </motion.div>
  )
}

/* Browser-preview (npm run dev) sample for the Results tab. */
export const MOCK_QTL_EXPORT = {
  markers: 338, individuals: 176, chromosomes: 12, linkage_markers: 338, unmatched_linkage: 0, unmatched_input: 41,
  unmatched_input_reasons: { 'Monomorphic marker': 33, 'Excess missing data': 6, 'Unplaced (chromosome 0)': 2 },
  duplicates_dropped: 0, missing_pct: 0.29, notes: [],
  preview: {
    markers: ['chr01_1044946', 'chr01_10645259', 'chr01_10846090', 'chr01_11250486', 'chr01_1175585', 'chr01_12216652', 'chr01_12335190', 'chr01_1410467'],
    chromosome: [1, 1, 1, 1, 1, 1, 1, 1], position: [0, 4.615, 9.231, 13.846, 18.462, 23.077, 27.692, 32.308],
    rows: [['BSP2-79', [2, 0, 2, 2, 1, 0, 2, 2]], ['BSP2-87', [0, 2, 0, 0, 0, 1, 2, 0]], ['BSP2-95', [1, 1, 0, 2, -1, 0, 0, 1]], ['BSP2-8', [2, 2, 2, 0, 0, 1, 1, 2]]]
      .map(([id, codes]) => ({ id, codes })),
  },
}
