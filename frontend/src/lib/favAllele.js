/**
 * Favourable Allele Analysis — data model and scoring.
 *
 * The page has no data of its own. The user uploads a rice 4K haplopanel genotype file; the Python backend matches its marker
 * IDs against the built-in trait-gene markers (tools/favourable_allele/backend) and sends back only the matched markers and the
 * genotypes of the user's lines. loadDataset() takes that result; everything below then scores it in the browser.
 *
 * Dataset:
 *   snps       one record per matched SNP / indel: gene id, gene name, trait, sub-category, marker, chromosome, position,
 *              REF, ALT, favourable allele, unfavourable allele
 *   varieties  { id, name } — the lines of the user's file
 *   calls      one string per line, one character per SNP: 0 / 1 / 2 (REF/REF, REF/ALT, ALT/ALT), x = no call
 */
export let SNPS = []
export let VARIETIES = []
export let LOST_GENES = []
export let TRAITS = []
let CALLS = []
let FAV_CODE = []

export function loadDataset(data) {
  SNPS = data.snps.map((s, i) => ({ ...s, i, geneLabel: s.name || s.gene, traitLabel: s.trait && s.trait !== '-' ? s.trait : 'Unspecified' }))
  VARIETIES = data.varieties.map((v, i) => ({ ...v, i, code: v.code || '' }))
  LOST_GENES = data.lostGenes || []
  CALLS = data.calls
  FAV_CODE = SNPS.map((s) => favCode(s.ref, s.alt, s.fav))
  TRAITS = [...new Set(SNPS.map((s) => s.traitLabel))].sort()
}

export const STATUSES = ['Favourable', 'Unfavourable', 'Heterozygous', 'Missing']
export const STATUS_COLOURS = { Favourable: '#2f9e5f', Unfavourable: '#e5484d', Heterozygous: '#f2b01e', Missing: '#b4beb8' }
export const STATUS_SHORT = { Favourable: 'fav', Unfavourable: 'unfav', Heterozygous: 'het', Missing: 'miss' }

/* ───────── genotype codes ─────────
 * 0 = homozygous REF · 1 = heterozygous · 2 = homozygous ALT · -1 (or no value) = no call
 */
export const DEFAULT_MODE = 'het1'
const ALLELE_OF = { het1: { 0: 'ref', 1: 'het', 2: 'alt' } }

/* ───────── favourable allele → REF (0) or ALT (1) ───────── */
const COMP = { A: 'T', T: 'A', G: 'C', C: 'G', N: 'N' }
const revComp = (s) => (s ? s.toUpperCase().split('').reverse().map((c) => COMP[c] || 'N').join('') : null)

/** 0 when the favourable allele is REF, 1 when it is ALT, null when it matches neither (also on the opposite strand).
 *  A deletion is written "-": it counts as an allele like any other (e.g. REF "-" and favourable "-" = REF is favourable). */
export function favCode(ref, alt, fav) {
  if (!fav) return null
  if (fav === ref) return 0
  if (fav === alt) return 1
  if (fav === '-') return null
  const rr = revComp(ref), ar = revComp(alt)
  if (fav === rr) return 0
  if (fav === ar) return 1
  return null
}

/** Raw stored code of a variety at a SNP: 0 / 1 / 2, or -1 for no call. */
export const rawCall = (v, s) => { const ch = CALLS[v][s]; return ch === 'x' ? -1 : ch.charCodeAt(0) - 48 }

/**
 * Status of one variety at one SNP.
 * @returns {{ status: 'Favourable'|'Unfavourable'|'Heterozygous'|'Missing', call: string, gt: number }}
 *   call = the two nucleotides (REF/ALT written out), "NN" when there is no call.
 */
export function statusOf(snp, variety, mode = DEFAULT_MODE) {
  const s = SNPS[snp.i ?? snp]
  const gt = rawCall(variety.i ?? variety, s.i)
  const allele = ALLELE_OF[mode][gt]
  const r = s.ref || '?', a = s.alt || '?'
  if (!allele) return { status: 'Missing', call: 'NN', gt }
  if (allele === 'het') return { status: 'Heterozygous', call: r + a, gt }
  const call = allele === 'ref' ? r + r : a + a
  const fc = FAV_CODE[s.i]
  if (fc === null) return { status: 'Missing', call, gt }
  return { status: (allele === 'ref' ? 0 : 1) === fc ? 'Favourable' : 'Unfavourable', call, gt }
}

/* ───────── traits ───────── */
const TRAIT_RULES = [
  ['grain quality', '#e8741c'], ['yield', '#d6408f'], ['floral', '#3b82d6'], ['plant growth', '#6aa81a'],
  ['hybrid', '#d19a0b'], ['biotic', '#d9382f'], ['multiple stress', '#7d8aa0'], ['abiotic', '#8a62d9'],
  ['nutrient', '#12a394'], ['dsr', '#8cae17'],
]
export const traitColour = (trait) => {
  const t = String(trait || '').toLowerCase()
  for (const [k, c] of TRAIT_RULES) if (t.includes(k)) return c
  return '#8a948f'
}

/** Genes (by label) with their trait; optionally only those of the given traits. */
export function geneList(traits) {
  const m = new Map()
  SNPS.forEach((s) => { if (!traits.length || traits.includes(s.traitLabel)) m.set(s.geneLabel, s.traitLabel) })
  return [...m.entries()].sort((a, b) => a[0].localeCompare(b[0]))
}

export const snpsFor = (traits, genes) => SNPS.filter((s) => (!traits.length || traits.includes(s.traitLabel)) && (!genes.length || genes.includes(s.geneLabel)))

/* ───────── aggregation ───────── */
const empty = () => ({ Favourable: 0, Unfavourable: 0, Heterozygous: 0, Missing: 0 })

/** Counts for one variety over a SNP set. `view` ('all' or a status in lower case) keeps only that status, like the dashboard. */
export function tally(snps, variety, mode, view = 'all') {
  const counts = empty()
  const favGenes = new Set(), unfavGenes = new Set()
  const byGene = new Map()
  const rows = []
  snps.forEach((s) => {
    const r = statusOf(s, variety, mode)
    if (view !== 'all' && view !== r.status.toLowerCase()) return
    counts[r.status]++
    if (r.status === 'Favourable') favGenes.add(s.gene)
    if (r.status === 'Unfavourable') unfavGenes.add(s.gene)
    let g = byGene.get(s.geneLabel)
    if (!g) byGene.set(s.geneLabel, (g = { gene: s.geneLabel, trait: s.traitLabel, ...empty(), total: 0 }))
    g[r.status]++
    g.total++
    rows.push({ snp: s, ...r })
  })
  return { counts, favGenes: favGenes.size, unfavGenes: unfavGenes.size, byGene: [...byGene.values()], rows }
}

/** Genes of a SNP set with their SNPs and the status counts of one variety. */
export function genesOf(snps, variety, mode) {
  const m = new Map()
  snps.forEach((s) => {
    let g = m.get(s.geneLabel)
    if (!g) {
      m.set(s.geneLabel, (g = {
        name: s.geneLabel, geneId: s.gene, trait: s.traitLabel, sub: s.sub && s.sub !== '-' ? s.sub : null,
        chrom: s.chrom, pos: s.pos, region: regionForTrait(s.trait), snps: [], counts: empty(),
      }))
    }
    g.snps.push(s)
    g.counts[statusOf(s, variety, mode).status]++
  })
  return [...m.values()]
}

/* ───────── rice-plant layout ───────── */
const REGIONS = {
  panicles: [[250, 70, 40, 40], [168, 88, 36, 38], [332, 88, 36, 38], [210, 120, 34, 26], [290, 120, 34, 26], [250, 142, 30, 22]],
  upperStem: [[250, 175, 30, 34], [205, 205, 24, 30], [295, 205, 24, 30]],
  leaves: [[92, 430, 52, 14], [70, 392, 46, 18], [84, 350, 42, 24], [116, 308, 38, 30], [158, 262, 30, 38],
    [408, 430, 52, 14], [430, 392, 46, 18], [416, 350, 42, 24], [384, 308, 38, 30], [342, 262, 30, 38],
    [200, 240, 24, 40], [300, 240, 24, 40], [250, 270, 22, 44]],
  lowerStem: [[232, 320, 22, 30], [290, 332, 22, 30], [210, 360, 20, 30], [300, 366, 20, 30], [246, 392, 22, 30], [285, 410, 20, 26], [218, 420, 20, 26]],
  stem: [[250, 260, 20, 70], [220, 320, 16, 44], [280, 320, 16, 44]],
  roots: [[250, 590, 60, 34], [196, 600, 42, 28], [304, 600, 42, 28], [224, 632, 36, 24], [276, 632, 36, 24], [250, 560, 40, 22]],
}
export const REGION_LABEL = {
  panicles: 'Panicle', upperStem: 'Upper stem', leaves: 'Leaves', lowerStem: 'Lower stem', stem: 'Stem', roots: 'Roots',
}
export function regionForTrait(trait) {
  const t = String(trait || '').toLowerCase()
  if (t.includes('grain quality') || t.includes('yield')) return 'panicles'
  if (t.includes('floral') || t.includes('hybrid fertility')) return 'upperStem'
  if (t.includes('plant growth')) return 'leaves'
  if (t.includes('biotic') && !t.includes('abiotic')) return 'lowerStem'
  if (t.includes('nutrient') || t.includes('abiotic')) return 'roots'
  if (t.includes('multiple stress')) return 'lowerStem'
  return 'stem'
}

/** A dot for every gene, placed inside the region of its trait and pushed apart so that each one can be clicked. */
export function placeGenes(genes) {
  const byRegion = {}
  genes.forEach((g) => (byRegion[g.region] = byRegion[g.region] || []).push(g))
  const placed = []
  Object.entries(byRegion).forEach(([region, list]) => {
    const subs = REGIONS[region] || REGIONS.stem
    list.forEach((gene, idx) => {
      const [cx, cy, rx, ry] = subs[idx % subs.length]
      const subIdx = Math.floor(idx / subs.length)
      const subTotal = Math.ceil(list.length / subs.length)
      const angle = subIdx * 2.39996323 + (idx % subs.length) * 0.7
      const frac = subTotal > 1 ? subIdx / (subTotal - 1) : 0.5
      const r = 0.4 + 0.6 * frac
      placed.push({ gene, x: cx + Math.cos(angle) * r * rx, y: cy + Math.sin(angle) * r * ry })
    })
  })
  const MIN = 11
  for (let it = 0; it < 18; it++) {
    for (let a = 0; a < placed.length; a++) {
      for (let b = a + 1; b < placed.length; b++) {
        const dx = placed[b].x - placed[a].x, dy = placed[b].y - placed[a].y
        const d = Math.hypot(dx, dy)
        if (d < MIN && d > 1e-4) {
          const push = (MIN - d) / 2, ux = dx / d, uy = dy / d
          placed[a].x -= ux * push; placed[a].y -= uy * push
          placed[b].x += ux * push; placed[b].y += uy * push
        }
      }
    }
  }
  return placed
}

/* ───────── panicle drawing (static geometry) ───────── */
const q = (t, A, B, C) => ({ x: (1 - t) ** 2 * A.x + 2 * (1 - t) * t * B.x + t * t * C.x, y: (1 - t) ** 2 * A.y + 2 * (1 - t) * t * B.y + t * t * C.y })
function grain(x, y, dirx, diry, L, W) {
  const len = Math.hypot(dirx, diry) || 1
  const ang = (Math.atan2(diry / len, dirx / len) * 180) / Math.PI
  const d = `M0 0 C ${L * 0.14} ${-W}, ${L * 0.46} ${-W * 1.05}, ${L * 0.64} ${-W * 0.72} C ${L * 0.84} ${-W * 0.4}, ${L * 0.95} ${-W * 0.16}, ${L} 0 `
    + `C ${L * 0.95} ${W * 0.16}, ${L * 0.84} ${W * 0.4}, ${L * 0.64} ${W * 0.72} C ${L * 0.46} ${W * 1.05}, ${L * 0.14} ${W}, 0 0 Z`
  return { d, seam: `M${L * 0.1} 0 Q ${L * 0.55} 0 ${L * 0.9} 0`, transform: `translate(${x.toFixed(1)},${y.toFixed(1)}) rotate(${ang.toFixed(1)})` }
}
function panicle(baseX, baseY, reach, tipDrop) {
  const sign = reach >= 0 ? 1 : -1
  const A = { x: baseX, y: baseY }
  const B = { x: baseX + reach * 0.56, y: baseY - Math.abs(reach) * 0.64 }
  const C = { x: baseX + reach, y: baseY + tipDrop }
  const grains = []
  const start = 0.05, N = 26
  for (let i = 0; i <= N; i++) {
    const t = start + (1 - start) * (i / N)
    const p = q(t, A, B, C)
    const tf = (t - start) / (1 - start)
    const gl = 14 - 4.5 * tf, gw = 5.8 - 2.0 * tf
    const dx = sign * (0.18 + 0.1 * tf)
    grains.push(grain(p.x - sign * 2, p.y, dx, 1, gl, gw), grain(p.x + sign * 2, p.y, dx, 1, gl * 0.94, gw))
  }
  const tip = q(1, A, B, C)
  grains.push(grain(tip.x, tip.y, sign * 0.28, 1, 11, 4.4))
  return { stem: `M${A.x} ${A.y} Q${B.x.toFixed(1)} ${B.y.toFixed(1)} ${C.x.toFixed(1)} ${C.y.toFixed(1)}`, grains }
}
export const PANICLES = [panicle(238, 150, -104, 54), panicle(246, 118, -66, 62), panicle(250, 94, -26, 64),
  panicle(256, 94, 26, 64), panicle(260, 118, 66, 62), panicle(268, 150, 104, 54)]
