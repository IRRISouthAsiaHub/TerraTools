import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import { viteSingleFile } from 'vite-plugin-singlefile'

// The app is loaded from file:// inside QtWebEngine. Chromium blocks
// <script type="module" src="..."> over file://, so everything (JS, CSS,
// fonts, background photo) is inlined into a single dist/index.html.
export default defineConfig({
  base: './',
  plugins: [react(), viteSingleFile()],
  build: {
    assetsInlineLimit: 100_000_000,
    chunkSizeWarningLimit: 4000,
  },
})
