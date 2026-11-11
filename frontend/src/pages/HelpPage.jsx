import { useMemo, useState } from 'react'
import { BookOpen, Search } from 'lucide-react'
import chapters from '../lib/manualChapters.json'

// The manual lives in public/help (built by docs/manual/build_manual.py) and is shown in a frame.
const PAGE = 'help/index.html'

export default function HelpPage() {
  const [target, setTarget] = useState('')
  const [q, setQ] = useState('')
  const shown = useMemo(() => {
    const n = q.trim().toLowerCase()
    return n ? chapters.filter((c) => c.title.toLowerCase().includes(n)) : chapters
  }, [q])
  let last = null

  return (
    <div className="help-page">
      <aside className="help-toc">
        <div className="help-head"><BookOpen size={18} /><b>User Manual</b></div>
        <label className="help-search"><Search size={14} />
          <input value={q} onChange={(e) => setQ(e.target.value)} placeholder="Find a chapter…" />
        </label>
        <nav>
          <button className={`help-link ${target === '' ? 'active' : ''}`} onClick={() => setTarget('')}>Start of the manual</button>
          {shown.map((c) => {
            const head = c.group !== last ? <div className="help-group" key={`g-${c.group}`}>{c.group}</div> : null
            last = c.group
            return (
              <div key={c.id}>
                {head}
                <button className={`help-link ${target === c.id ? 'active' : ''}`} onClick={() => setTarget(c.id)}>
                  <span className="help-num">{c.n}</span>{c.title}
                </button>
              </div>
            )
          })}
        </nav>
      </aside>
      <iframe key={target} className="help-frame" title="TerraTools User Manual" src={target ? `${PAGE}#${target}` : PAGE} />
    </div>
  )
}
