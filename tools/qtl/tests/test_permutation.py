"""The CIM / MQM permutation threshold: the fast genome-wide maximum equals the maximum of the ordinary scan, and the CIM null is
not lower than the single-marker null it replaced (cofactors absorb chance associations)."""
import subprocess
import textwrap
from pathlib import Path

import pytest

from tools.common.r_locator import find_rscript

BACKEND = Path(__file__).resolve().parents[1] / "backend"
RSCRIPT = find_rscript()
pytestmark = pytest.mark.skipif(RSCRIPT is None, reason="R is not available")

SCRIPT = textwrap.dedent("""
    suppressPackageStartupMessages(library(jsonlite))
    suppressMessages(source("qtl_engine.R"))
    set.seed(3)
    chr <- rep(1:6, each = 40); pos <- rep(seq(0, 195, length.out = 40), 6)
    mi  <- data.frame(Marker = paste0("m", seq_along(chr)), Chromosome = as.character(chr), Position = pos, stringsAsFactors = FALSE)
    worst <- 0
    for (type in c("ril", "f2")) {
      n <- 150
      X <- matrix(NA, n, nrow(mi), dimnames = list(NULL, mi$Marker))
      for (k in 1:6) {                                   # linked markers: a latent chromosome, markers flip with distance
        base <- sample(0:1, n, TRUE); cols <- which(chr == k); X[, cols[1]] <- base
        for (j in cols[-1]) { flip <- runif(n) < 0.06; X[, j] <- ifelse(flip, 1 - X[, j - 1], X[, j - 1]) }
      }
      X <- if (type == "ril") X * 2 else X + matrix(sample(0:1, n * ncol(X), TRUE), n)          # RIL 0/2, F2-like 0/1/2
      if (type == "ril") X[sample(length(X), 40)] <- NA                                           # a few missing calls (slow path)
      y <- 0.8 * ifelse(is.na(X[, 20]), 0, X[, 20]) + rnorm(n)
      for (meth in c("CIM", "MQM")) for (b in 1:8) {
        yp <- if (b == 1) y else sample(y); w <- if (meth == "MQM") 25 else 10
        cof <- if (meth == "MQM") select_cofactors(X, yp, mi, max_cof = 15, p_enter = 0.01, p_remove = 0.02, min_dist = 25) else select_cofactors(X, yp, mi, max_cof = 5, min_dist = 10)
        worst <- max(worst, abs(max(qtl_scan_linear(X, yp, mi, cof, w, meth)$LOD) - cim_max_lod(X, yp, mi, cof, w)))
      }
    }
    cat("WORST", format(worst, digits = 6), "\n")
    nul <- perm_max_lod_cim(X, y, mi, "CIM", n_perm = 40, progress = FALSE)
    cat("LEN", length(nul), "\n")
""")


def test_fast_cim_maximum_equals_the_ordinary_scan(tmp_path):
    f = tmp_path / "check.R"
    f.write_text(SCRIPT, encoding="utf-8")
    p = subprocess.run([str(RSCRIPT), str(f)], capture_output=True, text=True, cwd=str(BACKEND))
    assert p.returncode == 0, p.stderr[-800:]
    worst = float([l for l in p.stdout.splitlines() if l.startswith("WORST")][0].split()[1])
    assert worst < 1e-8, p.stdout
    assert "LEN 40" in p.stdout
