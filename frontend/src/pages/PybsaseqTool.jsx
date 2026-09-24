import { useState } from 'react'
import { Activity, SlidersHorizontal, LayoutGrid, FileText, Settings2, Ruler, Crosshair, FileCode, Image, FolderOpen, Terminal, Info, Download } from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, TextField, SelectField, Toggle, FilePicker, SubTabs, Frame,
  EmptyState, ImageGallery, FileList, LogBox, DownloadButtons,
} from '../components/kit.jsx'
import { useBridge, useTool, mockGeneric } from '../lib/bridge.js'

const TABS = [
  { id: 'input', label: 'Input & Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]
const SUB = [
  { id: 'report', label: 'HTML Report', icon: FileCode },
  { id: 'plots', label: 'Plots', icon: Image },
  { id: 'files', label: 'All Output Files', icon: FolderOpen },
  { id: 'log', label: 'Log', icon: Terminal },
]

const DEFAULTS = {
  bulksizes: '430,385', popstruct: 'F2', pvalues: '0.01,0.01', replication: 10000,
  slidingwindow: '2000000,10000', smoothing: '51,3', chromosomeOrder: '', readLength: 100,
  region: '-1', miscParams: '20,5,1', parentRef: 'no', adjustGap: false,
}

export default function PybsaseqTool() {
  const { bridge } = useBridge('pybsaseqBridge', () => mockGeneric())
  const tool = useTool(bridge, 'run_analysis', { accept: (r) => r.output_files && r.output_files.length > 0 })
  const [tab, setTab] = useState('input')
  const [sub, setSub] = useState('report')
  const [parents, setParents] = useState('')
  const [bulks, setBulks] = useState('')
  const [c, setC] = useState(DEFAULTS)
  const set = (k) => (v) => setC((p) => ({ ...p, [k]: v }))
  tool.onResult(() => { setTab('results'); setSub('report') })

  const run = () => {
    if (!parents || !bulks) return tool.fail('Please select both a parents file and a bulks file.')
    tool.run({
      parents_file: parents, bulks_file: bulks, bulksizes: c.bulksizes, popstruct: c.popstruct, pvalues: c.pvalues,
      replication: parseInt(c.replication, 10), slidingwindow: c.slidingwindow, smoothing: c.smoothing,
      chromosome_order: c.chromosomeOrder, read_length: parseInt(c.readLength, 10), region: c.region,
      misc_params: c.miscParams, adjust_gap: c.adjustGap, parent_ref: c.parentRef,
    })
  }
  const r = tool.result
  const files = (r && r.output_files) || []
  const pngs = files.filter((f) => f.ext === 'png')

  return (
    <ToolPage
      id="bsa" icon={Activity} tabs={TABS} tab={tab} setTab={setTab}
      title="PyBSA-Seq"
      subtitle="Bulk Segregant Analysis Sequencing — G-statistic & Δ(allele frequency) based SV/QTL identification"
      run={{ label: 'Run Analysis', runningLabel: 'Running BSA-Seq…', state: tool.state, status: tool.status, onRun: run }}
    >
      {tab === 'input' && (
        <TabPane id="input" key="input">
          <SectionHeader title="Input Files">Upload the GATK4-generated parents and bulks TSV/CSV files.</SectionHeader>
          <div className="row">
            <Card title="Parents File" icon={FileText}>
              <FilePicker value={parents} placeholder="parents.tsv / parents.csv"
                onBrowse={async () => { const p = await tool.browse('browse_parents_file'); if (p) setParents(p) }} />
            </Card>
            <Card title="Bulks File" icon={FileText} delay={0.05}>
              <FilePicker value={bulks} placeholder="bulks.tsv / bulks.csv"
                onBrowse={async () => { const p = await tool.browse('browse_bulks_file'); if (p) setBulks(p) }} />
            </Card>
          </div>
          <SectionHeader title="Analysis Parameters" />
          <div className="row">
            <Card title="Population & Bulk Sizes" icon={Settings2} delay={0.1}>
              <TextField label="Bulk sizes (first,second)" value={c.bulksizes} onChange={set('bulksizes')} />
              <SelectField label="Population structure" value={c.popstruct} onChange={set('popstruct')} options={['F2', 'RIL', 'BC']} />
              <TextField label="P-value cutoffs (real_data,simulation)" value={c.pvalues} onChange={set('pvalues')} />
              <TextField label="Number of replications" type="number" min={1000} max={50000} step={1000} value={c.replication} onChange={set('replication')} />
            </Card>
            <Card title="Windows & Smoothing" icon={Ruler} delay={0.15}>
              <TextField label="Sliding window (size,incremental_step)" value={c.slidingwindow} onChange={set('slidingwindow')} />
              <TextField label="Smoothing (window_len,polyorder)" value={c.smoothing} onChange={set('smoothing')} />
              <TextField label="Chromosome order (optional)" value={c.chromosomeOrder} onChange={set('chromosomeOrder')} placeholder="e.g. 1,2,3,4,5" />
              <TextField label="Read length" type="number" step={1} value={c.readLength} onChange={set('readLength')} />
            </Card>
            <Card title="Region & Misc" icon={Crosshair} delay={0.2}>
              <TextField label="Region of interest (chrm,start,end)" value={c.region} onChange={set('region')}
                hint="Add multiple regions as repeating triplets, e.g. 1,1000001,4000000,10,20000001,40000000" />
              <TextField label="Misc parameters (GQ cutoff, min SVs, mirror index)" value={c.miscParams} onChange={set('miscParams')} />
              <TextField label="Parent reference genome" value={c.parentRef} onChange={set('parentRef')} />
              <div style={{ marginTop: 10 }}><Toggle label="Adjust gap between subplots" checked={c.adjustGap} onChange={set('adjustGap')} /></div>
            </Card>
          </div>
        </TabPane>
      )}
      {tab === 'results' && (
        <TabPane id="results" key="results">
          <SectionHeader title="Results">{r ? `${files.length} output file(s)` : 'Run the analysis to see results here.'}</SectionHeader>
          {!r ? <EmptyState icon={Activity}>Run the analysis to see results here.</EmptyState> : (
            <>
              <SubTabs tabs={SUB} value={sub} onChange={setSub} id="bsasub" />
              {sub === 'report' && (r.html_report ? (
                <>
                  <DownloadButtons onSave={tool.save} items={[{ icon: Download, label: 'Download Report', path: r.html_report, name: 'BSASeq_report.html' }]} />
                  <Frame src={r.html_report} height={720} title="BSA-Seq report" />
                </>
              ) : <EmptyState icon={Info}>No HTML report was generated.</EmptyState>)}
              {sub === 'plots' && (pngs.length ? (
                <>
                  <ImageGallery images={pngs.map((f) => ({ src: f.path, label: f.name }))} />
                </>
              ) : <EmptyState icon={Info}>No plot images found.</EmptyState>)}
              {sub === 'files' && <FileList files={files} onSave={tool.save} />}
              {sub === 'log' && <LogBox>{r.stdout || r.message || '(no output)'}</LogBox>}
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
