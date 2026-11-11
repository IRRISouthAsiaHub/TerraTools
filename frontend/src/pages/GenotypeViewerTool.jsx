import { useEffect, useState } from 'react'
import { Rows3, SlidersHorizontal, LayoutGrid, Dna, Palette, Image as ImageIcon, FileText, ScanLine, Plus } from 'lucide-react'
import {
  ToolPage, TabPane, SectionHeader, Card, Toggle, FilePicker, TextField, SelectField, RadioChips, Note,
  EmptyState, SummaryGrid, ImageGallery, DownloadButtons,
} from '../components/kit.jsx'
import { useBridge, useTool, useSignal, callJson, mockGeneric } from '../lib/bridge.js'

const isVcf = (p) => /\.vcf(\.gz)?$/i.test(p || '')

const TABS = [
  { id: 'input', label: 'Input & Settings', icon: SlidersHorizontal },
  { id: 'results', label: 'Results', icon: LayoutGrid },
]
const MODES = [
  { value: 'normal', label: 'Normal — one image per sample' },
  { value: 'compare', label: 'Compare — all samples, one image per chromosome' },
  { value: 'zoomed', label: 'Zoomed — all samples, one region' },
]
const COLORINGS = [
  { value: '2-color', label: '2-color (REF / ALT)' },
  { value: '3-color', label: '3-color (REF / HET / ALT)' },
]
const DPIS = [
  { value: 200, label: '200 dpi — screen' },
  { value: 300, label: '300 dpi — print' },
  { value: 600, label: '600 dpi — publication' },
]
// the palettes are the ones of GenoSee's color_set.json; "normal" is its bright set
const PALETTE_NOTE = { normal: 'bright', 'Red-Blue-Yellow-Grey': 'mixed: REF red, ALT blue, HET yellow, missing grey' }

const Swatches = ({ p }) => (
  <span style={{ display: 'inline-flex', gap: 3, marginRight: 8, verticalAlign: 'middle' }}>
    {[p['0|0'] || p['0'], p['0|1'], p['1|1'] || p['1']].filter(Boolean).map((c, i) => (
      <span key={i} style={{ width: 14, height: 14, borderRadius: 3, background: c, border: '1px solid rgba(0,0,0,.25)', display: 'inline-block' }} />
    ))}
  </span>
)

export default function GenotypeViewerTool() {
  const { bridge } = useBridge('genotypeViewerBridge', () => mockGeneric(['speciesAdded']))
  const tool = useTool(bridge, 'run_viewer')
  const [tab, setTab] = useState('input')
  const [file, setFile] = useState('')
  const [info, setInfo] = useState(null)               // { samples } or { error }
  const [species, setSpecies] = useState('')
  const [speciesList, setSpeciesList] = useState([])
  const [palettes, setPalettes] = useState({})
  const [mode, setMode] = useState('normal')
  const [coloring, setColoring] = useState('2-color')
  const [fill, setFill] = useState(true)
  const [names, setNames] = useState(false)
  const [guides, setGuides] = useState(true)
  const [graphic, setGraphic] = useState(true)
  const [palette, setPalette] = useState('normal')
  const [chr, setChr] = useState('')
  const [start, setStart] = useState(0)
  const [end, setEnd] = useState(0)
  const [dpi, setDpi] = useState(300)
  const [pdf, setPdf] = useState(false)
  const [samples, setSamples] = useState('')
  const [vcfNames, setVcfNames] = useState(true)
  const [legend, setLegend] = useState({ REF: '', ALT: '', HET: '' })
  const [fasta, setFasta] = useState('')
  const [newSpecies, setNewSpecies] = useState('')
  const [speciesMsg, setSpeciesMsg] = useState('')

  useEffect(() => {
    if (!bridge) return
    callJson(bridge, 'get_options').then((o) => {
      if (!o) return
      setSpeciesList(o.species || [])
      setPalettes(o.palettes || {})
      setSpecies((s) => s || (o.species || [])[0] || '')
    })
  }, [bridge])

  useSignal(bridge, 'speciesAdded', (json) => {
    let r
    try { r = JSON.parse(json) } catch { return }
    if (r.status === 'success') {
      setSpeciesList(r.species)
      setSpecies(r.added)
      setSpeciesMsg(`${r.added} added to the chromosome length database.`)
    } else setSpeciesMsg(r.message || 'The species could not be added.')
  })

  tool.onResult(() => setTab('results'))

  const browse = async () => {
    const p = await tool.browse()
    if (!p) return
    setFile(p)
    setInfo(null)
    const r = await callJson(bridge, 'get_samples', p)
    setInfo(r && r.status === 'success' ? { samples: r.samples } : { error: (r && r.message) || 'The file could not be read.' })
  }
  const browseFasta = async () => { const p = await tool.browse('browse_fasta'); if (p) setFasta(p) }
  const addSpecies = () => {
    if (!fasta || !newSpecies.trim()) return setSpeciesMsg('Choose the genome FASTA file and type the species name.')
    setSpeciesMsg('Reading the chromosome lengths…')
    bridge.add_species(fasta, newSpecies.trim())
  }

  const run = () => {
    if (!file) return tool.fail('Please select the input file.')
    if (!species) return tool.fail('Please choose the species.')
    tool.run({
      input_file: file, species, drawing_mode: mode, coloring_mode: coloring,
      fill: fill ? 'on' : 'off', display_marker_names: names ? 'on' : 'off', color_set: palette,
      chr, start, end, dpi, pdf, vcf_marker_names: vcfNames, legend_labels: legend, guide_lines: guides, graphic,
      samples: samples.split(/[,;\n]+/).map((s) => s.trim()).filter(Boolean),
    })
  }

  const r = tool.result
  const images = r ? r.images.map((i) => ({ src: i.png, label: i.name })) : []
  const pdfs = r ? r.images.filter((i) => i.pdf).map((i) => ({ icon: FileText, label: `${i.name} (PDF)`, path: i.pdf, name: `${i.name}.pdf` })) : []
  const nSel = samples.split(/[,;\n]+/).map((s) => s.trim()).filter(Boolean).length
  const nDraw = nSel || (info && info.samples ? info.samples.length : 0)

  return (
    <ToolPage
      id="gv" icon={Rows3} tabs={TABS} tab={tab} setTab={setTab}
      title="Graphical Genotype Visualisation"
      subtitle="Used for parent recovery — see which parent every stretch of every chromosome came from (GenoSee, Hashimoto 2024)"
      run={{ label: 'Draw Graphical Genotypes', runningLabel: 'Drawing…', state: tool.state, status: tool.status, onRun: run }}
    >
      {tab === 'input' && (
        <TabPane id="input" key="input">
          <SectionHeader title="Genotype Data">
            A CSV (or VCF) with the chromosome in the first column, the marker name in the second, the position (bp) in the third,
            then one column per sample. The column names do not matter.
          </SectionHeader>
          <Card title="Genotype File" icon={Dna}>
            <FilePicker value={file} onBrowse={browse} placeholder="Select the genotype file (CSV, VCF or VCF.GZ)" />
            <Note>
              Genotype calls can be simple (<code>A</code> / <code>B</code> / <code>H</code> / <code>N</code>), phased (<code>0|1</code>, <code>.|.</code> missing)
              or non-phased (<code>0/1</code>, <code>./.</code> missing). Multi-allelic sites are not supported.
              A VCF is converted with GenoSee's own VCF converter (the <code>GT</code> field of every sample).
            </Note>
            {isVcf(file) && (
              <div style={{ marginTop: 12 }}>
                <Toggle label="Name the markers CHROM_POS (untick to use the ID column of the VCF)" checked={vcfNames} onChange={setVcfNames} />
              </div>
            )}
            {info && info.error && <div className="warning-box" style={{ marginTop: 12, marginBottom: 0 }}><span>{info.error}</span></div>}
            {info && info.samples && <Note>{info.samples.length} samples: {info.samples.slice(0, 6).join(', ')}{info.samples.length > 6 ? ', …' : ''}</Note>}
          </Card>
          <div className="row">
            <Card title="Species" icon={ScanLine} delay={0.04}>
              <SelectField label="Chromosome lengths of" value={species} onChange={setSpecies} options={speciesList.length ? speciesList : ['']} />
              <div className="field-hint">Chromosome names are matched to the species' own (1, chr1 and chr01 all match chr01). Markers on other chromosomes are not drawn.</div>
            </Card>
            <Card title="Samples to draw (optional)" icon={Dna} delay={0.08}>
              <TextField label="Sample names, separated by commas" value={samples} onChange={setSamples}
                placeholder="empty = every sample" />
              {mode === 'normal' && nDraw > 20 && (
                <div className="field-hint">Normal mode writes one image per sample: {nDraw} images. Name the samples you need to draw fewer.</div>
              )}
            </Card>
          </div>

          <SectionHeader title="Drawing" />
          <Card title="Drawing mode" icon={Rows3} delay={0.05}>
            <RadioChips options={MODES} value={mode} onChange={setMode} />
            {mode === 'zoomed' && (
              <div className="row" style={{ marginTop: 14 }}>
                <TextField label="Chromosome" value={chr} onChange={setChr} placeholder="chr01" />
                <TextField label="Start (bp)" type="number" value={start} onChange={setStart} min={0} />
                <TextField label="End (bp)" type="number" value={end} onChange={setEnd} min={0} />
              </div>
            )}
          </Card>
          <div className="row">
            <Card title="Colouring" icon={Palette} delay={0.06}>
              <RadioChips options={COLORINGS} value={coloring} onChange={setColoring} />
              <div style={{ height: 12 }} />
              <Toggle label="Graphic style (rounded, shaded chromosomes with no black lines inside)" checked={graphic} onChange={setGraphic} />
              <Toggle label="Fill the area between two markers" checked={fill} onChange={setFill} />
              <Toggle label="Show marker names in the figure (off = chromosomes only, cleaner for publication)" checked={names} onChange={setNames} />
              {names && mode !== 'normal' && (
                <Toggle label="Dashed line from each marker name across the samples (off = names only)" checked={guides} onChange={setGuides} />
              )}
            </Card>
            <Card title="Colour palette" icon={Palette} delay={0.1}>
              <RadioChips value={palette} onChange={setPalette}
                options={Object.keys(palettes).length
                  ? Object.entries(palettes).map(([k, p]) => ({
                    value: k,
                    label: <span><Swatches p={p} />{k}{PALETTE_NOTE[k] ? ` — ${PALETTE_NOTE[k]}` : ''}</span>,
                  }))
                  : ['normal']} />
            </Card>
          </div>
          <Card title="Legend names (optional)" icon={Palette} delay={0.11}>
            <div className="row">
              <TextField label="REF (first colour)" value={legend.REF} onChange={(v) => setLegend({ ...legend, REF: v })} placeholder="REF  →  e.g. Parent 1" />
              <TextField label="ALT (second colour)" value={legend.ALT} onChange={(v) => setLegend({ ...legend, ALT: v })} placeholder="ALT  →  e.g. Parent 2" />
              <TextField label="HET (heterozygous)" value={legend.HET} onChange={(v) => setLegend({ ...legend, HET: v })} placeholder="HET" />
            </div>
            <Note>
              GenoSee only knows the allele types: A / <code>0|0</code> is drawn as REF and B / <code>1|1</code> as ALT. Type here what they stand for,
              for example <b>Parent 1 (IR64)</b> and <b>Parent 2 (Swarna)</b>; the words replace REF and ALT in the legend of every figure. Empty = REF / ALT / HET.
            </Note>
          </Card>
          <Card title="Image quality" icon={ImageIcon} delay={0.12}>
            <div className="row" style={{ gap: 36, alignItems: 'flex-end' }}>
              <div style={{ flex: 1, minWidth: 240 }}>
                <SelectField label="PNG resolution" value={dpi} onChange={(v) => setDpi(Number(v))} options={DPIS} />
              </div>
              <div style={{ flex: 1, minWidth: 240 }}>
                <Toggle label="Also save a vector PDF (for publication)" checked={pdf} onChange={setPdf} />
              </div>
            </div>
          </Card>

          <Card title="Add a species (optional)" icon={Plus} delay={0.14}>
            <FilePicker value={fasta} onBrowse={browseFasta} placeholder="Genome FASTA (.fa, .fasta, .gz)" />
            <div className="row" style={{ alignItems: 'flex-end', marginTop: 10 }}>
              <TextField label="Species name" value={newSpecies} onChange={setNewSpecies} placeholder="Oryza_sativa_IR64" />
              <button className="btn btn-secondary" onClick={addSpecies}><Plus size={15} /> Add species</button>
            </div>
            <Note>
              The chromosome lengths are read from the FASTA file and kept with the species list, so the species is available the next time too.
              {speciesMsg && <><br /><b>{speciesMsg}</b></>}
            </Note>
          </Card>
        </TabPane>
      )}
      {tab === 'results' && (
        <TabPane id="results" key="results">
          <SectionHeader title="Results">
            {r ? `${r.n_images} image${r.n_images === 1 ? '' : 's'} — ${r.n_samples} samples, ${r.n_markers} markers` : 'Draw the graphical genotypes to see them here.'}
          </SectionHeader>
          {!r ? <EmptyState icon={Rows3}>Draw the graphical genotypes to see them here.</EmptyState> : (
            <>
              <SummaryGrid items={[
                { label: 'Images', value: r.n_images },
                { label: 'Samples', value: r.n_samples },
                { label: 'Markers', value: r.n_markers },
                { label: 'Chromosomes', value: r.n_chromosomes },
              ]} />
              {r.n_off_species > 0 && <Note>{r.n_off_species} marker(s) lie on chromosomes that are not in the chosen species and were not drawn.</Note>}
              <ImageGallery single={images.length === 1} images={images} />
              {pdfs.length > 0 && <DownloadButtons onSave={tool.save} items={pdfs} />}
            </>
          )}
        </TabPane>
      )}
    </ToolPage>
  )
}
