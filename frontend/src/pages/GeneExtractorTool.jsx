import { useState } from 'react'
import { ScanSearch, Upload, MapPin, Dna, FileSpreadsheet, Info, Download, Waypoints, ListChecks, Settings2 } from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, FilePicker, EmptyState, DataTable, TextField, SelectField, RadioChips,
} from '../components/kit.jsx'
import { useBridge, useTool, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'extraction', label: 'Data Extraction', icon: Upload },
  { id: 'markers', label: 'Regions', icon: MapPin },
  { id: 'genes', label: 'Genes in Gene Regions', icon: Dna },
]

const INTERVAL = [
  { value: 'ci', label: 'Confidence interval of the QTL (narrow)' },
  { value: 'flank', label: 'From the left to the right flanking marker (wider)' },
]

const POS_KINDS = [
  { value: 'map', label: 'Genetic map (Marker, Chromosome, bp)' },
  { value: 'vcf', label: 'No map: use a VCF' },
]

const MATCH = [
  { value: 'inside', label: 'Completely inside the region' },
  { value: 'overlap', label: 'Overlapping the region (also partly)' },
]

const dash = (_c, v) => (v === null || v === undefined || v === '' ? '—' : undefined)

const pre = { margin: '8px 0', padding: '10px 14px', borderRadius: 8, background: 'rgba(0,0,0,0.06)', fontSize: 13, overflowX: 'auto' }
const POS_EXAMPLE = ['Marker,Chromosome,bp', 'M001,chr01,19308987', 'M005,chr01,21011077'].join('\n')
const QTL_EXAMPLE = ['QTL_ID,QTL_Start,QTL_End', 'QTL_3.1,4000000,4300000', 'QTL_5.1,22500000,23100000'].join('\n')

export default function GeneExtractorTool() {
  const { bridge } = useBridge('geneExtractorBridge', () => mockGeneric())
  const tool = useTool(bridge, 'run_extraction')
  const [tab, setTab] = useState('extraction')
  const [f, setF] = useState({ marker_path: '', gff_path: '', fasta_path: '', positions_path: '', interval: 'ci', match: 'inside' })
  const [posKind, setPosKind] = useState('map')
  tool.onResult(() => setTab('markers'))

  const pick = (key, slot) => async () => { const p = await tool.browse(slot); if (p) setF((x) => ({ ...x, [key]: p })) }
  const set = (key) => (v) => setF((x) => ({ ...x, [key]: v }))
  const useQtlMapping = async () => {
    const p = await tool.browse('latest_qtl_intervals')
    if (p) setF((x) => ({ ...x, marker_path: p }))
    else tool.fail('No QTL Mapping output found. Run QTL Mapping first, or browse to the QTL results.')
  }
  const useQtlTable = async () => {
    const p = await tool.browse('latest_qtl_table')
    if (p) setF((x) => ({ ...x, marker_path: p }))
    else tool.fail('No QTL table found. Run QTL Mapping first, or browse to the QTL table (qtl_significant_*.csv).')
  }
  const run = () => {
    if (!f.marker_path) return tool.fail('Please upload the QTL results.')
    if (!f.positions_path) return tool.fail(posKind === 'vcf' ? 'Please select the VCF file.' : 'Please select the genetic map (Marker, Chromosome, bp), or choose “No map: use a VCF”.')
    if (!f.gff_path) return tool.fail('Please upload the annotation (GFF3) file.')
    tool.run({
      marker_path: f.marker_path, gff_path: f.gff_path, match: f.match, positions_path: f.positions_path, interval: f.interval,
      ...(f.fasta_path ? { fasta_path: f.fasta_path } : {}),
    })
  }
  const download = (type) => bridge && bridge.download_result(type)
  const r = tool.result

  return (
    <ToolPage
      id="gene" icon={ScanSearch} tabs={TABS} tab={tab} setTab={setTab}
      title="Gene Extractor"
      subtitle="Genes, their annotation and sequences inside the QTL of QTL Mapping"
      run={tab === 'extraction' ? { label: 'Extract Data', runningLabel: 'Extracting…', state: tool.state, status: tool.status, onRun: run } : undefined}
    >
      {tab === 'extraction' && (
        <TabPane id="extraction" key="extraction">
          <SectionHeader title="Instructions">Provide the QTL results, the genetic map (or a VCF) and the annotation, then extract.</SectionHeader>
          <div className="note" style={{ marginTop: 0, marginBottom: 20 }}>
            <p>1. Upload the QTL results (CSV, TSV or Excel) — the tool finds the columns by their <b>name</b>, in any order:</p>
            <p><b>QTL intervals</b> — the QTL intervals of QTL Mapping (click “Use QTL Mapping output”), or any table with <b>Chromosome</b>, <b>Start</b> and <b>End</b> columns
              in base pairs (QTL_ID / QTL_Start / QTL_End, Left / Right, Left_Mb / Right_Mb are all understood). When there is no Chromosome column it is read from the QTL_ID (QTL_3.1 is chromosome 3).</p>
            <pre style={pre}>{QTL_EXAMPLE}</pre>
            <p>The <b>QTL table</b> of QTL Mapping (qtl_significant_*.csv, click “Use QTL table”) also works: its positions are in cM, so the tool needs the base-pair position of its markers,
              and that is what the genetic map below gives it.</p>
            <p style={{ marginTop: 12 }}>2. Upload the <b>genetic map</b> with the base-pair position of the markers — a table with <b>Marker</b>, <b>Chromosome</b> and <b>bp</b> columns
              (the names must be the same as in the QTL table; the “Physical Map” sheet of a Linkage Mapping workbook works too). <b>If you have no genetic map, choose “No map: use a VCF”</b>{' '}
              and select the VCF: its CHROM and POS give the position and its ID the marker name (a VCF without IDs cannot be matched). The positions of the QTL’s left, peak and right markers place
              the confidence interval (cM) on the genome; the file only needs the markers of the QTL.</p>
            <pre style={pre}>{POS_EXAMPLE}</pre>
            <p style={{ marginTop: 12 }}>3. Upload the gene annotation (GFF3) — the genes in each region are read from its gene lines,
              with all of their details (name, function, biotype …). Chr01, chr1 and 1 are treated as the same chromosome</p>
            <p>4. Optionally add the reference genome FASTA to get the sequence of every gene (minus-strand genes are reverse-complemented)</p>
            <p>5. Set the options below if needed, then click “Extract Data”</p>
            <p>6. Download the gene table, the <b>gene ID list</b> or the sequences from the “Genes in Gene Regions” tab</p>
          </div>
          <div className="row">
            <Card title="QTL Results (input)" icon={FileSpreadsheet}>
              <FilePicker value={f.marker_path} placeholder="No file selected" onBrowse={pick('marker_path', 'browse_markers')} />
              <div style={{ marginTop: 10 }}>
                <button className="btn btn-secondary" onClick={useQtlMapping}><Waypoints size={15} /> Use QTL Mapping output</button>{' '}
                <button className="btn btn-secondary" onClick={useQtlTable}><Waypoints size={15} /> Use QTL table (cM)</button>
              </div>
            </Card>
            <Card title="Genetic Map (required)" icon={MapPin} delay={0.03}>
              <RadioChips options={POS_KINDS} value={posKind} onChange={(v) => { setPosKind(v); setF((x) => ({ ...x, positions_path: '' })) }} />
              <div style={{ marginTop: 10 }}>
                <FilePicker value={f.positions_path} placeholder={posKind === 'vcf' ? 'Select the VCF (.vcf or .vcf.gz)' : 'Select the genetic map'}
                  onBrowse={pick('positions_path', posKind === 'vcf' ? 'browse_vcf' : 'browse_positions')} />
              </div>
              <div className="field-hint" style={{ marginTop: 8 }}>
                {posKind === 'vcf' ? 'CHROM and POS give the position of each marker; the ID (or CHROM_POS) must be the marker name of the QTL table.'
                  : 'Marker, Chromosome and bp: places the QTL (in cM) on the genome.'}
              </div>
            </Card>
            <Card title="Annotation File (GFF3)" icon={FileSpreadsheet} delay={0.06}>
              <FilePicker value={f.gff_path} placeholder="No file selected" onBrowse={pick('gff_path', 'browse_gff')} />
            </Card>
            <Card title="Reference Genome FASTA (optional)" icon={Dna} delay={0.12}>
              <FilePicker value={f.fasta_path} placeholder="No file selected" onBrowse={pick('fasta_path', 'browse_fasta')} />
            </Card>
          </div>
          <div className="row">
            <Card title="Options" icon={Settings2} delay={0.18}>
              <SelectField label="QTL interval (QTL table in cM)" value={f.interval} onChange={set('interval')} options={INTERVAL}
                hint="Only for a QTL table in cM. The confidence interval is placed on the genome from the positions of the QTL's markers, as in the QTL Intervals file." />
              <SelectField label="A gene is taken when it lies" value={f.match} onChange={set('match')} options={MATCH}
                hint="Completely inside keeps only genes that fit in the region; overlapping also takes genes cut by its border." />
            </Card>
          </div>
        </TabPane>
      )}
      {tab === 'markers' && (
        <TabPane id="markers" key="markers">
          <SectionHeader title="Regions">The regions read from the input file and the number of genes in each.</SectionHeader>
          {!r ? <EmptyState icon={Info}>No data yet. Run extraction first.</EmptyState> : (
            <>
              {(r.warnings || []).map((w) => <div className="note" key={w} style={{ marginBottom: 12 }}><p>{w}</p></div>)}
              <div className="info-line">
                QTL intervals · {(r.region_summary || []).length} region(s) · {r.gene_ids} gene(s)
              </div>
              <DataTable rows={r.region_summary} columns={r.region_columns || []} renderCell={dash} maxHeight={480} />
            </>
          )}
        </TabPane>
      )}
      {tab === 'genes' && (
        <TabPane id="genes" key="genes">
          <SectionHeader title="Genes in Gene Regions">Genes found in each region, with the details from the annotation.</SectionHeader>
          {!r ? <EmptyState icon={Info}>No data yet. Run extraction first.</EmptyState> : (
            <>
              <DataTable rows={r.genes} columns={r.columns || []} renderCell={dash} maxHeight={480} />
              <div style={{ marginTop: 16, display: 'flex', gap: 12, flexWrap: 'wrap' }}>
                <button className="btn btn-primary" onClick={() => download('gene_list')} disabled={!r.gene_ids}><ListChecks size={16} /> Download Gene ID List ({r.gene_ids})</button>
                <button className="btn btn-secondary" onClick={() => download('genes')}><Download size={16} /> Download Genes Data</button>
                {r.fasta_file && <button className="btn btn-secondary" onClick={() => download('fasta')}><Download size={16} /> Download Gene Sequences (FASTA)</button>}
              </div>
              <div className="field-hint" style={{ marginTop: 8 }}>
                The gene ID list is a text file with one gene ID per line (each gene once). The Excel file also holds it, with the regions each gene was found in, on its Gene_ID_list sheet.
              </div>
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
