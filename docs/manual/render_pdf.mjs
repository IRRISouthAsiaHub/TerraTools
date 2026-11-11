// Prints an HTML file to PDF with headless Edge (Chrome DevTools protocol), with a page-number footer.
// node render_pdf.mjs <input.html> <output.pdf>
import { spawn, spawnSync } from 'node:child_process'
import { writeFileSync, mkdtempSync } from 'node:fs'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { pathToFileURL } from 'node:url'

const [, , input, output] = process.argv
// headless Edge starts helper processes: end the whole tree, or the next run finds a stale instance on the port
const quit = (code) => {
  try {
    if (process.platform === 'win32') {
      // the browser may hand over to a process that is not the one we started, so find ours by its debugging port
      spawnSync('powershell', ['-NoProfile', '-Command',
        `Get-CimInstance Win32_Process -Filter "Name='msedge.exe'" | Where-Object { $_.CommandLine -match 'remote-debugging-port=${PORT}' } | ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }`])
    } else edge.kill()
  } catch { /* already gone */ }
  process.exit(code)
}
const EDGE = process.env.EDGE || 'C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe'
const PORT = 9444
const sleep = (ms) => new Promise((r) => setTimeout(r, ms))

const edge = spawn(EDGE, ['--headless=new', '--disable-gpu', `--remote-debugging-port=${PORT}`,
  `--user-data-dir=${mkdtempSync(join(tmpdir(), 'tt-pdf-'))}`, 'about:blank'], { stdio: 'ignore' })
let targets
for (let i = 0; i < 60; i++) {
  try { targets = await (await fetch(`http://127.0.0.1:${PORT}/json`)).json(); if (targets.find((t) => t.type === 'page')) break } catch { /* wait */ }
  await sleep(250)
}
const ws = new WebSocket(targets.find((t) => t.type === 'page').webSocketDebuggerUrl)
await new Promise((r) => (ws.onopen = r))
let id = 0
const waits = new Map()
const events = []
ws.onmessage = (ev) => {
  const m = JSON.parse(ev.data)
  if (m.id && waits.has(m.id)) { waits.get(m.id)(m); waits.delete(m.id) } else if (m.method) events.push(m.method)
}
const send = (method, params = {}) => new Promise((res) => { const my = ++id; waits.set(my, res); ws.send(JSON.stringify({ id: my, method, params })) })

await send('Page.enable')
await send('Page.navigate', { url: pathToFileURL(input).href })
for (let i = 0; i < 240 && !events.includes('Page.loadEventFired'); i++) await sleep(250)
await sleep(1500)   // images and fonts
const footer = `<div style="width:100%;font-size:8px;color:#5a7566;padding:0 16mm;display:flex;justify-content:space-between;font-family:Segoe UI,Arial,sans-serif">
  <span>TerraTools User Manual</span><span>Page <span class="pageNumber"></span> of <span class="totalPages"></span></span></div>`
const r = await send('Page.printToPDF', {
  printBackground: true, preferCSSPageSize: true, displayHeaderFooter: true,
  headerTemplate: '<span></span>', footerTemplate: footer,
})
if (!r.result) { console.error('printToPDF failed', JSON.stringify(r)); quit(1) }
writeFileSync(output, Buffer.from(r.result.data, 'base64'))
console.log('wrote', output)
ws.close()
quit(0)
