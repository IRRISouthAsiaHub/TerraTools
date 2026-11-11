import { memo } from 'react'
import { PANICLES } from '../../lib/favAllele.js'

/* The rice plant itself (roots, tillers, leaves, drooping panicles) never changes, so it is drawn once. */
const Plant = memo(function Plant() {
  return (
    <>
      <defs>
        <linearGradient id="faBladeLight" x1="0" y1="0" x2="1" y2="0"><stop offset="0%" stopColor="#74bd52" /><stop offset="48%" stopColor="#93d268" /><stop offset="100%" stopColor="#5da642" /></linearGradient>
        <linearGradient id="faBladeMid" x1="0" y1="0" x2="1" y2="0"><stop offset="0%" stopColor="#4f9c3a" /><stop offset="48%" stopColor="#73bf52" /><stop offset="100%" stopColor="#3f8a30" /></linearGradient>
        <linearGradient id="faBladeDark" x1="0" y1="0" x2="1" y2="0"><stop offset="0%" stopColor="#2f7a2c" /><stop offset="48%" stopColor="#479c3c" /><stop offset="100%" stopColor="#256a22" /></linearGradient>
        <linearGradient id="faSheath" x1="0" y1="1" x2="0" y2="0"><stop offset="0%" stopColor="#79bd54" /><stop offset="60%" stopColor="#8ccf63" /><stop offset="100%" stopColor="#74bd52" /></linearGradient>
        <linearGradient id="faPanStem" x1="0" y1="0" x2="1" y2="1"><stop offset="0%" stopColor="#a3d05d" /><stop offset="100%" stopColor="#7ab63f" /></linearGradient>
        <linearGradient id="faGrain" x1="0" y1="0" x2="1" y2="1"><stop offset="0%" stopColor="#fbe07a" /><stop offset="50%" stopColor="#f1c54a" /><stop offset="100%" stopColor="#d99c28" /></linearGradient>
        <linearGradient id="faRoot" x1="0" y1="0" x2="0" y2="1"><stop offset="0%" stopColor="#cbb78a" /><stop offset="45%" stopColor="#bda572" /><stop offset="100%" stopColor="#a98e57" /></linearGradient>
      </defs>

      <g stroke="url(#faRoot)" fill="none" strokeLinecap="round">
        {[['M250 540 C 238 576, 224 600, 206 636', 3.4], ['M250 540 C 262 576, 276 600, 294 636', 3.4], ['M250 540 C 244 586, 240 614, 236 652', 3.6],
          ['M250 540 C 256 586, 260 614, 264 652', 3.6], ['M250 540 C 250 586, 250 620, 250 658', 3.8], ['M250 540 C 232 572, 208 592, 182 620', 2.4],
          ['M250 540 C 268 572, 292 592, 318 620', 2.4], ['M250 540 C 226 568, 196 586, 168 606', 2], ['M250 540 C 274 568, 304 586, 332 606', 2],
          ['M250 540 C 240 570, 222 588, 200 612', 1.3], ['M250 540 C 260 570, 278 588, 300 612', 1.3], ['M236 652 C 232 660, 224 666, 214 672', 1.1],
          ['M264 652 C 268 660, 276 666, 286 672', 1.1], ['M206 636 C 200 646, 192 652, 184 660', 1], ['M294 636 C 300 646, 308 652, 316 660', 1],
          ['M250 658 C 246 666, 244 670, 240 676', 1], ['M250 658 C 254 666, 256 670, 260 676', 1], ['M182 620 C 176 628, 172 634, 166 642', 0.9],
          ['M318 620 C 324 628, 328 634, 334 642', 0.9]].map(([d, w]) => <path key={d} d={d} strokeWidth={w} />)}
      </g>

      <path d="M250 540 C 236 484, 230 430, 240 374 C 248 422, 252 486, 250 540 Z" fill="url(#faSheath)" stroke="#4f9c3a" strokeWidth="0.9" />
      <path d="M250 540 C 264 484, 270 430, 260 374 C 252 422, 248 486, 250 540 Z" fill="url(#faSheath)" stroke="#4f9c3a" strokeWidth="0.9" />
      <path d="M250 540 C 244 486, 246 426, 250 376 C 254 426, 256 486, 250 540 Z" fill="url(#faBladeMid)" opacity="0.85" />

      <g className="fa-sway-l">
        <path d="M248 536 C 226 440, 184 322, 132 196 C 110 168, 92 148, 80 136 C 100 150, 128 176, 156 212 C 200 330, 234 444, 252 538 Z" fill="url(#faBladeDark)" stroke="#1f5a1c" strokeWidth="0.7" />
        <path d="M250 532 C 228 440, 190 330, 140 206" fill="none" stroke="#7cc25a" strokeWidth="1.1" opacity="0.4" />
      </g>
      <g className="fa-sway-r">
        <path d="M252 536 C 274 440, 316 322, 368 196 C 390 168, 408 148, 420 136 C 400 150, 372 176, 344 212 C 300 330, 266 444, 248 538 Z" fill="url(#faBladeDark)" stroke="#1f5a1c" strokeWidth="0.7" />
        <path d="M250 532 C 272 440, 310 330, 360 206" fill="none" stroke="#7cc25a" strokeWidth="1.1" opacity="0.4" />
      </g>
      <g className="fa-sway-l">
        <path d="M248 538 C 224 430, 192 312, 150 200 C 130 170, 114 148, 104 132 C 122 150, 146 180, 170 216 C 208 330, 232 440, 252 539 Z" fill="url(#faBladeMid)" stroke="#256a22" strokeWidth="0.7" />
        <path d="M250 534 C 226 430, 196 320, 156 208" fill="none" stroke="#a6e07e" strokeWidth="1.2" opacity="0.5" />
      </g>
      <g className="fa-sway-r">
        <path d="M252 538 C 276 430, 308 312, 350 200 C 370 170, 386 148, 396 132 C 378 150, 354 180, 330 216 C 292 330, 268 440, 248 539 Z" fill="url(#faBladeMid)" stroke="#256a22" strokeWidth="0.7" />
        <path d="M250 534 C 274 430, 304 320, 344 208" fill="none" stroke="#a6e07e" strokeWidth="1.2" opacity="0.5" />
      </g>
      <g className="fa-sway-l">
        <path d="M249 539 C 238 430, 218 290, 198 180 C 190 148, 184 122, 180 100 C 192 126, 202 158, 212 198 C 230 320, 244 432, 251 540 Z" fill="url(#faBladeLight)" stroke="#3f8a30" strokeWidth="0.7" />
        <path d="M250 535 C 236 428, 216 298, 202 188" fill="none" stroke="#bdeb90" strokeWidth="1.2" opacity="0.55" />
      </g>
      <g className="fa-sway-r">
        <path d="M251 539 C 262 430, 282 290, 302 180 C 310 148, 316 122, 320 100 C 308 126, 298 158, 288 198 C 270 320, 256 432, 249 540 Z" fill="url(#faBladeLight)" stroke="#3f8a30" strokeWidth="0.7" />
        <path d="M250 535 C 264 428, 284 298, 298 188" fill="none" stroke="#bdeb90" strokeWidth="1.2" opacity="0.55" />
      </g>
      <path d="M250 540 C 244 430, 238 230, 242 122 C 244 92, 248 68, 250 52 C 252 68, 256 92, 258 122 C 262 230, 256 430, 250 540 Z" fill="url(#faBladeLight)" stroke="#3f8a30" strokeWidth="0.7" />
      <path d="M250 536 L250 58" stroke="#cdf0a0" strokeWidth="1.2" opacity="0.5" />

      {PANICLES.map((p, i) => (
        <g key={i}>
          <path d={p.stem} fill="none" stroke="url(#faPanStem)" strokeWidth="2" strokeLinecap="round" />
          {p.grains.map((g, j) => (
            <g key={j} transform={g.transform}>
              <path d={g.d} fill="url(#faGrain)" stroke="#b07d1c" strokeWidth="0.5" />
              <path d={g.seam} stroke="#a8730f" strokeWidth="0.4" opacity="0.5" fill="none" />
            </g>
          ))}
        </g>
      ))}
    </>
  )
})

/**
 * The plant with one dot per gene. `dots` = [{ gene, x, y, colour, dim }].
 * Hovering shows the gene in the tooltip, clicking opens its details.
 */
export default function RicePlant({ dots, hover, onHover, onSelect, svgRef }) {
  return (
    <svg ref={svgRef} viewBox="0 0 500 680" xmlns="http://www.w3.org/2000/svg" className="fa-plant" role="img" aria-label="Rice plant with one dot for each gene">
      <rect width="500" height="680" fill="#ffffff" />
      <Plant />
      <g>
        {dots.map(({ gene, x, y, colour, dim }) => {
          const on = hover === gene.name
          return (
            <g key={gene.name} className="fa-dot" style={{ opacity: dim ? 0.18 : 1 }}
              onMouseEnter={() => onHover(gene.name)} onMouseLeave={() => onHover(null)} onClick={() => onSelect(gene)}
              tabIndex={0} role="button" aria-label={`${gene.name}, ${gene.trait}`}
              onKeyDown={(e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); onSelect(gene) } }}>
              {on && <circle cx={x} cy={y} r="9.5" fill="none" stroke={colour} strokeWidth="1.6" opacity="0.55" />}
              <circle cx={x} cy={y} r={on ? 5.6 : 4} fill={colour} stroke="#ffffff" strokeWidth="1.2" />
            </g>
          )
        })}
      </g>
    </svg>
  )
}
