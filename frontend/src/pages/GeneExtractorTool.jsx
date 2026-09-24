import { useState } from 'react'
import { ScanSearch, Upload, MapPin, Dna, FileSpreadsheet, Info, Download } from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, FilePicker, EmptyState, DataTable,
} from '../components/kit.jsx'
import { useBridge, useTool, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'extraction', label: 'Data Extraction', icon: Upload },
  { id: 'markers', label: 'Markers with Positions', icon: MapPin },
  { id: 'genes', label: 'Genes in Gene Regions', icon: Dna },
]

const dash = (_c, v) => (v === null || v === undefined || v === '' ? '—' : undefined)

export default function GeneExtractorTool() {
  const { bridge } = useBridge('geneExtractorBridge', () => mockGeneric())
  const tool = useTool(bridge, 'run_extraction')
  const [tab, setTab] = useState('extraction')
  const [f, setF] = useState({ extracted_markers_path: '', gene_markers_path: '', rice_gene_path: '' })
  tool.onResult(() => setTab('markers'))

  const pick = (key, slot) => async () => { const p = await tool.browse(slot); if (p) setF((x) => ({ ...x, [key]: p })) }
  const run = () => {
    if (!f.extracted_markers_path || !f.gene_markers_path || !f.rice_gene_path) return tool.fail('Please upload all three required files.')
    tool.run(f)
  }
  const download = (type) => bridge && bridge.download_result(type)
  const r = tool.result

  // plain function (not a component) so the table keeps its state across renders
  const renderResults = ({ rows, type, cols, label }) => (
    !r ? <EmptyState icon={Info}>No data yet. Run extraction first.</EmptyState> : (
      <>
        <DataTable rows={rows} columns={cols} renderCell={dash} maxHeight={480} />
        <div style={{ marginTop: 16 }}>
          <button className="btn btn-secondary" onClick={() => download(type)}><Download size={16} /> {label}</button>
        </div>
      </>
    )
  )

  return (
    <ToolPage
      id="gene" icon={ScanSearch} tabs={TABS} tab={tab} setTab={setTab}
      title="Gene Marker & Gene Extractor"
      subtitle="Physical marker positions & genes within gene regions"
      run={tab === 'extraction' ? { label: 'Extract Data', runningLabel: 'Extracting…', state: tool.state, status: tool.status, onRun: run } : undefined}
    >
      {tab === 'extraction' && (
        <TabPane id="extraction" key="extraction">
          <SectionHeader title="Instructions">Provide the three Excel files, then extract.</SectionHeader>
          <div className="note" style={{ marginTop: 0, marginBottom: 20 }}>
            <p>1. Upload your extracted markers Excel file (columns: Markers, Chromosome, Position)</p>
            <p>2. Upload your gene markers Excel file (columns: Markers, Chromosome)</p>
            <p>3. Upload the rice gene list Excel file (columns: Chromosome, Start, Stop, Strand, Gene_id, Description)</p>
            <p>4. Click “Extract Data” to process the files</p>
            <p>5. Download the results from the tabs above</p>
          </div>
          <div className="row">
            <Card title="Extracted Markers File" icon={FileSpreadsheet}>
              <FilePicker value={f.extracted_markers_path} placeholder="No file selected" onBrowse={pick('extracted_markers_path', 'browse_extracted_markers')} />
            </Card>
            <Card title="Gene Markers File" icon={FileSpreadsheet} delay={0.06}>
              <FilePicker value={f.gene_markers_path} placeholder="No file selected" onBrowse={pick('gene_markers_path', 'browse_gene_markers')} />
            </Card>
            <Card title="Rice Gene List File" icon={FileSpreadsheet} delay={0.12}>
              <FilePicker value={f.rice_gene_path} placeholder="No file selected" onBrowse={pick('rice_gene_path', 'browse_rice_genes')} />
            </Card>
          </div>
        </TabPane>
      )}
      {tab === 'markers' && (
        <TabPane id="markers" key="markers">
          <SectionHeader title="Markers with Positions">Physical positions for the extracted markers.</SectionHeader>
          {renderResults({ rows: r && r.markers, type: 'markers', cols: ['Markers', 'Chromosome', 'Position'], label: 'Download Markers Data' })}
        </TabPane>
      )}
      {tab === 'genes' && (
        <TabPane id="genes" key="genes">
          <SectionHeader title="Genes in Gene Regions">Genes falling inside each gene-marker region.</SectionHeader>
          {renderResults({ rows: r && r.genes, type: 'genes', label: 'Download Genes Data',
            cols: ['Chromosome', 'Start', 'Stop', 'Strand', 'Gene_id', 'Description', 'gene_start', 'gene_stop'] })}
        </TabPane>
      )}
    </ToolPage>
  )
}
