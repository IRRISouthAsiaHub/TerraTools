import { createContext, useSyncExternalStore } from 'react'

/* Which tool page is rendering a component (set by the app shell). */
export const RouteContext = createContext('')

/* Run state per tool route ('idle' | 'running' | 'success' | 'error'),
   shared with the sidebar and the "finished" notifications. */
let states = {}
const subs = new Set()

export function setRunState(route, state) {
  if (!route || states[route] === state) return
  states = { ...states, [route]: state }
  subs.forEach((f) => f())
}

export function useRunStates() {
  return useSyncExternalStore(
    (f) => { subs.add(f); return () => subs.delete(f) },
    () => states,
  )
}

/* ───────────── Small confirmations ("Saved to …") ───────────── */

let notices = []
let noticeId = 0
const noticeSubs = new Set()
const emitNotices = () => noticeSubs.forEach((f) => f())

export function dismissNotice(id) {
  notices = notices.filter((n) => n.id !== id)
  emitNotices()
}

/** pushNotice({ text, kind: 'ok' | 'err' }) — shown as a toast for a few seconds. */
export function pushNotice({ text, kind = 'ok' }) {
  const id = ++noticeId
  notices = [...notices, { id, text, kind }]
  emitNotices()
  setTimeout(() => dismissNotice(id), kind === 'err' ? 9000 : 5000)
}

export function useNotices() {
  return useSyncExternalStore((f) => { noticeSubs.add(f); return () => noticeSubs.delete(f) }, () => notices)
}
