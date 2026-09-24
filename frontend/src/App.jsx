import { useEffect, useState } from 'react'
import Dashboard from './pages/Dashboard.jsx'
import ClusteringTool from './pages/ClusteringTool.jsx'
import GwasTool from './pages/GwasTool.jsx'

// Each Qt window loads this same bundle with a different hash:
//   index.html            -> launcher dashboard
//   index.html#/gwas      -> GWAS window
// Tools not yet migrated keep using their original HTML pages.
const ROUTES = {
  '': Dashboard,
  gwas: GwasTool,
  clustering: ClusteringTool,
}

const parse = () => window.location.hash.replace(/^#\/?/, '')

export default function App() {
  const [route, setRoute] = useState(parse)
  useEffect(() => {
    const on = () => setRoute(parse())
    window.addEventListener('hashchange', on)
    return () => window.removeEventListener('hashchange', on)
  }, [])
  const Page = ROUTES[route] || Dashboard
  return <Page />
}
