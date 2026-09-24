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
