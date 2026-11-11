import { useMemo } from 'react'
import { Dna, Rows3, Settings2, FileText, Grip, X } from 'lucide-react'
import { Card, RadioChips, Note, SummaryGrid, ImageGallery, DownloadButtons, EmptyState } from '../../components/kit.jsx'
import { TRAITS, VARIETIES, traitColour } from '../../lib/favAllele.js'
import { MultiSelect } from './Pickers.jsx'

const DPIS = [
  { value: 150, label: '150 dpi · quick' },
  { value: 300, label: '300 dpi' },
  { value: 600, label: '600 dpi · publication' },
]

/** The Heatmap tab: choose genes and lines, draw the balloon heatmap with the R script, show it. The state lives in the page. */
export default function HeatmapTab({
  geneOptions, geneTrait, genes, setGenes, lineOptions, lines, setLines, dpi, setDpi,
  filterGenes, filterLines, tool,
}) {
  const r = tool.result
  const addTrait = (t) => {
    const add = geneOptions.filter((o) => geneTrait.get(o.value) === t).map((o) => o.value)
    setGenes([...genes, ...add.filter((g) => !genes.includes(g))])
  }
  const traitCount = useMemo(() => {
    const m = new Map()
    geneOptions.forEach((o) => m.set(geneTrait.get(o.value), (m.get(geneTrait.get(o.value)) || 0) + 1))
    return m
  }, [geneOptions, geneTrait])
    const files = r ? r.files.map((f) => ({ icon: FileText, label: f.label, path: f.path, name: f.name })) : []

  return (
    <>
      <div className="grid-2">
        <Card title="Genes" icon={Dna}>
          <MultiSelect label="Genes in the heatmap" allLabel="Choose genes…" clearLabel="Clear" noun="genes" options={geneOptions} value={genes} onChange={setGenes} search />
          <div className="field">
            <label className="field-title">Add all the genes of a trait</label>
            <div className="fa-quick">
              {TRAITS.map((t) => (
                <button key={t} type="button" className="fa-quickbtn" style={{ '--c': traitColour(t) }} onClick={() => addTrait(t)} title={`Add the ${traitCount.get(t) || 0} genes of ${t}`}>
                  <i />{t}<small>{traitCount.get(t) || 0}</small>
                </button>
              ))}
            </div>
          </div>
          <div className="fa-quick-row">
            <button type="button" className="btn btn-secondary sm" onClick={() => setGenes(filterGenes)} disabled={!filterGenes.length}>Genes of the filters ({filterGenes.length})</button>
            <button type="button" className="btn btn-secondary sm" onClick={() => setGenes([])} disabled={!genes.length}>Clear</button>
          </div>
          {genes.length > 0 && (
            <div className="fa-chips">
              {genes.map((g) => (
                <span key={g} className="fa-chip2" style={{ '--c': traitColour(geneTrait.get(g)) }}>
                  <i />{g}
                  <button type="button" onClick={() => setGenes(genes.filter((x) => x !== g))} aria-label={`Remove ${g}`}><X size={11} /></button>
                </span>
              ))}
            </div>
          )}
        </Card>
        <Card title="Lines" icon={Rows3} delay={0.05}>
          <MultiSelect label="Lines in the heatmap" allLabel="Choose lines…" clearLabel="Clear" noun="lines" options={lineOptions} value={lines} onChange={setLines} search />
          <div className="fa-quick-row">
            <button type="button" className="btn btn-secondary sm" onClick={() => setLines(VARIETIES.map((v) => v.id))}>All lines ({VARIETIES.length})</button>
            <button type="button" className="btn btn-secondary sm" onClick={() => setLines(filterLines)} disabled={!filterLines.length}>Lines selected on the Your 4K File tab ({filterLines.length})</button>
            <button type="button" className="btn btn-secondary sm" onClick={() => setLines([])} disabled={!lines.length}>Clear</button>
          </div>
          <Note>
            {lines.length ? <><b>{lines.length}</b> of the {VARIETIES.length} lines of your file chosen.</> : <>Choose the lines of your file to compare.</>}
          </Note>
        </Card>
      </div>
      <Card title="Figure" icon={Settings2} delay={0.08}>
        <div className="field"><label className="field-title">Resolution</label><RadioChips options={DPIS} value={dpi} onChange={setDpi} /></div>
        <Note>
          The figure is drawn by the R script supplied with the Heatmap_Final folder, on the genes and lines chosen here. Each balloon is the fraction of a gene’s markers
          at which the line is homozygous for the favourable allele; a heterozygous call counts as not favourable and a cross means no call. Before drawing, markers whose
          favourable allele matches neither REF nor ALT, and markers missing in 80 % or more of the chosen lines, are removed.
        </Note>
      </Card>

      <div className="fa-hm-result">
        {!r ? (
          <EmptyState icon={Grip}>Choose the genes and lines, then press <b>Draw Balloon Heatmap</b>.</EmptyState>
        ) : (
          <>
            <SummaryGrid items={[
              { label: 'Genes', value: r.n_genes },
              { label: 'Lines', value: r.n_lines },
              { label: 'Markers used', value: `${r.n_markers_used} / ${r.n_markers_in}` },
              { label: 'Size', value: `${r.width_mm} × ${r.height_mm} mm` },
            ]} />
            {r.dropped.length > 0 && (
              <Note>
                <b>{r.dropped.length} marker{r.dropped.length === 1 ? ' was' : 's were'} removed by the quality control:</b>{' '}
                {r.dropped.slice(0, 8).map((d) => `${d.gene} ${d.marker.split('_').slice(0, 2).join('_')} (${d.reason})`).join('; ')}
                {r.dropped.length > 8 ? `; and ${r.dropped.length - 8} more (see the QC report below).` : '.'}
              </Note>
            )}
            <ImageGallery single images={[{ src: r.image, label: 'Favourable allele balloon heatmap' }]} />
            <DownloadButtons onSave={tool.save} items={files} />
          </>
        )}
      </div>
    </>
  )
}
