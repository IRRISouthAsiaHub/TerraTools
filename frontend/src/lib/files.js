import { useCallback } from 'react'
import { useBridge } from './bridge.js'
import { pushNotice } from './runState.js'

/** "Hybrid Performance" -> "Hybrid_Performance.png" */
export const imageFileName = (label, fallback = 'image') => {
  const base = String(label || fallback).replace(/[\/:*?"<>|\u0000-\u001f]+/g, '_').replace(/\s+/g, '_').replace(/^[._]+|[._]+$/g, '')
  return (base || fallback) + (/\.[A-Za-z0-9]{2,4}$/.test(base) ? '' : '.png')
}

/* Plain-browser fallback (npm run dev): there is no Qt, so use an <a download>. */
function browserSave(src, name) {
  const a = document.createElement('a')
  a.href = src
  a.download = name
  document.body.appendChild(a)
  a.click()
  a.remove()
}
const mockFiles = () => ({
  save_image: (src, name, cb) => { browserSave(src, name); cb && cb(name) },
  save_images: (json, cb) => {
    const items = JSON.parse(json)
    items.forEach((it, i) => setTimeout(() => browserSave(it.src, it.name), i * 150))
    cb && cb(JSON.stringify({ status: 'success', saved: items.length, folder: 'Downloads', failed: [] }))
  },
})

/**
 * useFiles() -> { saveImage(src, label), saveImages([{ src, label }]) }
 * Opens the system Save dialog (one image) or Choose Folder dialog (many)
 * and confirms with a toast. Works for any tool's images.
 */
export function useFiles() {
  const { bridge } = useBridge('terraFiles', mockFiles)

  const saveImage = useCallback((src, label) => {
    if (!bridge || !src) return
    const name = imageFileName(label)
    bridge.save_image(src, name, (res) => {
      if (!res) return                                   // cancelled
      if (String(res).startsWith('ERROR')) return pushNotice({ kind: 'err', text: `Could not save the image — ${String(res).slice(7)}` })
      pushNotice({ text: `Image saved: ${res}` })
    })
  }, [bridge])

  const saveImages = useCallback((items) => {
    if (!bridge || !items.length) return
    const payload = items.map((it) => ({ src: it.src, name: imageFileName(it.label) }))
    bridge.save_images(JSON.stringify(payload), (json) => {
      let r
      try { r = JSON.parse(json) } catch { r = { status: 'error', message: 'Unexpected response.' } }
      if (r.status === 'cancelled') return
      if (r.status !== 'success') return pushNotice({ kind: 'err', text: r.message || 'Could not save the images.' })
      const bad = (r.failed || []).length
      pushNotice({
        kind: bad ? 'err' : 'ok',
        text: `${r.saved} image${r.saved === 1 ? '' : 's'} saved to ${r.folder}` + (bad ? ` — ${bad} could not be saved` : ''),
      })
    })
  }, [bridge])

  return { saveImage, saveImages }
}
