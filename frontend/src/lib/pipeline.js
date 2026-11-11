import { useSyncExternalStore } from 'react'

/**
 * Files handed from one tool to the next in the genomic-selection workflow:
 *   Clustering  ->  Additive Dominance  ->  genomic prediction
 * Kept in memory only: the paths point at result folders of this session.
 *
 *   clustering: { train_csv, test_csv, train_has_pheno, test_has_pheno }
 *   matrices:   { kind, train_pheno_file, train_geno_file, test_pheno_file, test_geno_file }
 */
let state = {}
const subs = new Set()

export const getPipe = () => state

export function setPipe(patch) {
  state = { ...state, ...patch }
  subs.forEach((f) => f())
}

export function usePipe() {
  return useSyncExternalStore((f) => { subs.add(f); return () => subs.delete(f) }, getPipe)
}
