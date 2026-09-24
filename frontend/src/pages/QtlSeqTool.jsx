import { useState } from 'react'
import { Dna, SlidersHorizontal, LayoutGrid, FolderOpen, ChartArea, Star, Table2, Download, ChartLine } from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, TextField, FilePicker, Note,
  EmptyState, StatCard, ImageGallery, DownloadButtons, TableSlot, DataTable,
} from '../components/kit.jsx'
import { useBridge, useTool, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'config', label: 'Configuration', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]

const READS = [
  { key: 'parent_file_r1', label: 'Parent reads — R1 / forward (FASTQ or BAM)', ph: 'Select parent forward reads' },
  { key: 'parent_file_r2', label: 'Parent reads — R2 / reverse (optional, paired-end)', ph: 'Select parent reverse reads (leave empty for single-end)' },
  { key: 'bulk1_file_r1', label: 'Bulk 1 reads — R1 / forward (FASTQ or BAM)', ph: 'Select bulk 1 forward reads' },
  { key: 'bulk1_file_r2', label: 'Bulk 1 reads — R2 / reverse (optional, paired-end)', ph: 'Select bulk 1 reverse reads (leave empty for single-end)' },
  { key: 'bulk2_file_r1', label: 'Bulk 2 reads — R1 / forward (FASTQ or BAM)', ph: 'Select bulk 2 forward reads' },
  { key: 'bulk2_file_r2', label: 'Bulk 2 reads — R2 / reverse (optional, paired-end)', ph: 'Select bulk 2 reverse reads (leave empty for single-end)' },
]

const WIN_COLS = ['CHROM', 'POS', 'N_SNPS', 'deltaSNP', 'CI_95', 'CI_99']
const SNP_COLS = ['CHROM', 'POS', 'REF', 'PARENT_DEPTH', 'BULK1_DEPTH', 'BULK1_SNP_INDEX', 'BULK2_DEPTH', 'BULK2_SNP_INDEX', 'deltaSNP']

const delta = (c, v) =>
  c === 'deltaSNP' && typeof v === 'number'
    ? <span className={`badge ${v >= 0 ? '' : 'warn'}`}>{Number.isInteger(v) ? v : v.toFixed(3)}</span>
    : undefined

const Tbl = ({ title, rows, columns }) => (
  <div className="plot-block" style={{ marginBottom: 22 }}>
    <div className="plot-title">{title}</div>
    {rows && rows.length ? <DataTable rows={rows} columns={columns} renderCell={delta} maxHeight={340} /> : <div className="frame-wrap"><div className="empty-state" style={{ padding: 20 }}>No rows.</div></div>}
  </div>
)

export default function QtlSeqTool() {
  const { bridge } = useBridge('qtlSeqBridge', () => mockGeneric())
  const tool = useTool(bridge, 'run_qtlseq')
  const [tab, setTab] = useState('config')
  const [files, setFiles] = useState({ ref_file: '' })
  const [c, setC] = useState({ min_depth: 10, min_total_depth: 50, ref_allele_freq: 0.1, window_size: 1000, step_size: 100, threads: 2 })
  const setF = (k) => async (slot) => { const p = await tool.browse(slot); if (p) setFiles((f) => ({ ...f, [k]: p })) }
  const setP = (k) => (v) => setC((p) => ({ ...p, [k]: v }))
  tool.onResult(() => setTab('results'))

  const run = () => {
    if (!files.ref_file || !files.parent_file_r1 || !files.bulk1_file_r1 || !files.bulk2_file_r1)
      return tool.fail('Please select the reference and at least R1 reads for parent, bulk1, and bulk2.')
    tool.run({
      ref_file: files.ref_file,
      ...Object.fromEntries(READS.map((r) => [r.key, files[r.key] || ''])),
      min_depth: Number(c.min_depth), min_total_depth: Number(c.min_total_depth),
      ref_allele_freq: Number(c.ref_allele_freq),
      window_size: Number(c.window_size) * 1000, step_size: Number(c.step_size) * 1000,
      threads: Number(c.threads),
    })
  }
  const r = tool.result
  const s = (r && r.summary) || {}

  return (
    <ToolPage
      id="qs" icon={Dna} tabs={TABS} tab={tab} setTab={setTab}
      title="QTL-seq — Bulk Segregant Analysis"
      subtitle="Align reads, call SNPs and locate QTL regions with the Δ(SNP-index) method"
      run={{ label: 'Run QTL-seq Analysis', runningLabel: 'Running QTL-seq…', state: tool.state, status: tool.status, onRun: run }}
    >
      {tab === 'config' && (
        <TabPane id="config" key="config">
          <SectionHeader title="Input Files">Select the reference genome and the read files for the parent and both bulks.</SectionHeader>
          <Card title="Reference & Reads" icon={FolderOpen}>
            <div className="field">
              <label className="field-title">Reference genome (FASTA)</label>
              <FilePicker value={files.ref_file} placeholder="Select reference.fasta" onBrowse={() => setF('ref_file')('browse_ref_file')} />
            </div>
            {READS.map((rd) => (
              <div className="field" key={rd.key}>
                <label className="field-title">{rd.label}</label>
                <FilePicker value={files[rd.key]} placeholder={rd.ph} onBrowse={() => setF(rd.key)('browse_reads_file')} />
              </div>
            ))}
            <Note>
              FASTQ inputs are aligned automatically (BWA-MEM → sorted, indexed BAM). Provide R2 for paired-end alignment,
              or leave it empty for single-end reads. You can also point R1 directly at an already-sorted <code>.bam</code> file
              to skip alignment (R2 is ignored in that case).
            </Note>
          </Card>
          <Card title="Parameters" icon={SlidersHorizontal} delay={0.08}>
            <div className="grid-2" style={{ gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))' }}>
              <TextField label="Min. per-sample depth" type="number" min={1} value={c.min_depth} onChange={setP('min_depth')} />
              <TextField label="Min. total bulk depth" type="number" min={1} value={c.min_total_depth} onChange={setP('min_total_depth')} />
              <TextField label="Parent ref. allele freq. ≥" type="number" min={0} max={1} step={0.05} value={c.ref_allele_freq} onChange={setP('ref_allele_freq')} />
              <TextField label="Window size (kb)" type="number" min={1} value={c.window_size} onChange={setP('window_size')} />
              <TextField label="Step size (kb)" type="number" min={1} value={c.step_size} onChange={setP('step_size')} />
              <TextField label="Threads" type="number" min={1} value={c.threads} onChange={setP('threads')} />
            </div>
          </Card>
        </TabPane>
      )}
      {tab === 'results' && (
        <TabPane id="results" key="results">
          <SectionHeader title="Results">Δ(SNP-index) results across the genome.</SectionHeader>
          {!r ? <EmptyState icon={ChartLine}>Run the analysis to see Δ(SNP-index) results here.</EmptyState> : (
            <>
              <div className="summary-grid">
                <StatCard className="summary-card" value={s.n_snps ?? 0} label="SNPs analyzed" />
                <StatCard className="summary-card" value={s.n_windows ?? 0} label="Windows" delay={0.06} />
                <StatCard className="summary-card" value={s.n_chromosomes ?? 0} label="Chromosomes" delay={0.12} />
                <StatCard className="summary-card" value={s.n_candidate_windows ?? 0} label="Candidate QTL windows" delay={0.18} />
              </div>
              <div className="plot-title">Δ(SNP-index) Manhattan Plot</div>
              {r.plot ? <ImageGallery single images={[{ src: r.plot, label: 'Δ(SNP-index) Manhattan plot' }]} /> : <EmptyState icon={ChartArea}>No plot generated.</EmptyState>}
              <DownloadButtons onSave={tool.save} items={[
                r.plot && { icon: Download, label: 'Plot (PNG)', path: r.plot, name: 'qtlseq_plot.png' },
                r.snp_csv && { icon: Download, label: 'SNP data (CSV)', path: r.snp_csv, name: 'snp_data.csv' },
                r.window_csv && { icon: Download, label: 'Window data (CSV)', path: r.window_csv, name: 'window_data.csv' },
              ].filter(Boolean)} />
              <Note>Candidate windows are those where the smoothed Δ(SNP-index) exceeds both the 95% confidence interval and |0.3| — the suggested threshold for likely QTL regions.</Note>
              <div style={{ height: 18 }} />
              <Tbl title="Candidate QTL Windows" rows={r.candidate_windows} columns={WIN_COLS} />
              <Tbl title="Sliding Window Results" rows={r.window_table} columns={WIN_COLS} />
              <Tbl title="Top SNPs (by |Δ(SNP-index)|)" rows={r.snp_table} columns={SNP_COLS} />
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
