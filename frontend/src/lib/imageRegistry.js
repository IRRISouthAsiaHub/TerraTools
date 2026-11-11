import { useEffect, useSyncExternalStore } from 'react'

/*
 * Which result images each tool page is currently showing, so the page can
 * offer "Download all images". Galleries register while mounted, so this
 * always matches what is on screen (the active tab).
 */

const byRoute = new Map()          // route -> Map(galleryId -> images[])
const subs = new Set()
const snapshots = new Map()        // route -> { version, list }
let version = 0

const emit = () => { version += 1; subs.forEach((f) => f()) }

export function registerImages(route, id, images) {
  if (!byRoute.has(route)) byRoute.set(route, new Map())
  byRoute.get(route).set(id, images)
  emit()
}

export function unregisterImages(route, id) {
  const m = byRoute.get(route)
  if (m && m.delete(id)) emit()
}

function snapshot(route) {
  const cached = snapshots.get(route)
  if (cached && cached.version === version) return cached.data
  const groups = byRoute.get(route) ? [...byRoute.get(route).values()] : []
  const data = { galleries: groups.length, images: groups.flat() }
  snapshots.set(route, { version, data })
  return data
}

/** { galleries, images } for a route. */
export function useRouteImages(route) {
  return useSyncExternalStore((f) => { subs.add(f); return () => subs.delete(f) }, () => snapshot(route))
}

/** Register a gallery's images for the lifetime of the component. */
export function useRegisterImages(route, id, images) {
  const key = images.map((i) => `${i.src}|${i.label}`).join('\n')
  useEffect(() => {
    if (!images.length) return undefined
    registerImages(route, id, images)
    return () => unregisterImages(route, id)
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [route, id, key])
}
