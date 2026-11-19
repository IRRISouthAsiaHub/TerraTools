"""The machine-learning scans (RF / SVR / XGBoost): importance is tested by permutation, so a run on noise does not report QTL
as a matter of course, a planted QTL is found at the right place, and the fast SVR importance equals the exact loop."""
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
    set.seed(11)
    n <- 160; pos <- seq(0, 96, by = 4)                       # 2 chromosomes x 25 markers, RIL-like
    mk <- function(ch) { M <- matrix(0, n, length(pos)); for (i in 1:n) { g <- rbinom(1, 1, .5); M[i, 1] <- g
        for (j in 2:length(pos)) { if (runif(1) < 0.5 * tanh(4 / 50)) g <- 1 - g; M[i, j] <- g } }; colnames(M) <- paste0("c", ch, "_", pos); M }
    X  <- cbind(mk(1), mk(2))
    mi <- data.frame(Marker = colnames(X), Chromosome = rep(c("1", "2"), each = length(pos)), Position = rep(pos, 2), stringsAsFactors = FALSE)
    y  <- 5 + 1.6 * (X[, "c1_40"] - 0.5) + rnorm(n)             # one QTL, chromosome 1 at 40 cM

    # 1. fast SVR importance == the slow exact loop
    fit <- svm(x = X, y = y, kernel = "radial", scale = TRUE)
    set.seed(5); f <- svr_importance(fit, X, y)
    set.seed(5); mse0 <- mean((y - predict(fit, X))^2); sl <- numeric(ncol(X))
    for (j in seq_len(ncol(X))) { a <- 0; for (s in 1:3) { Xp <- X; Xp[, j] <- sample(Xp[, j]); a <- a + mean((y - predict(fit, Xp))^2) }; sl[j] <- (a / 3 - mse0) / mse0 }
    cat("SVRDIFF", format(max(abs(f - sl)), digits = 6), "\\n")

    # 2. the planted QTL is found by every ML method, at about the right place
    for (m in c("RF", "SVR", "XGBOOST")) {
      out <- capture.output(sc <- suppressWarnings(qtl_scan(X, y, mi, method = m, ml_perm = 30, ml_alpha = 0.05)))
      top <- sc[which.max(sc$LOD), ]
      cat("PEAK", m, top$Chromosome, top$Position, format(top$LOD, digits = 4), format(top$P_Value, digits = 3), "\\n")
    }

    # 3. noise: the genome-wide error rate is about alpha, not 100 %
    hits <- 0
    for (r in 1:6) {
      yn <- rnorm(n)
      out <- capture.output(sc <- suppressWarnings(qtl_scan(X, yn, mi, method = "XGBOOST", ml_perm = 40, ml_alpha = 0.05)))
      if (!"Message" %in% names(sc) && any(sc$LOD >= 1)) hits <- hits + 1
    }
    cat("NOISEHITS", hits, "\\n")
""")


@pytest.fixture(scope="module")
def output(tmp_path_factory):
    f = tmp_path_factory.mktemp("ml") / "check.R"
    f.write_text(SCRIPT, encoding="utf-8")
    p = subprocess.run([str(RSCRIPT), str(f)], capture_output=True, text=True, cwd=str(BACKEND))
    assert p.returncode == 0, p.stderr[-1000:]
    return p.stdout


def test_fast_svr_importance_equals_the_exact_loop(output):
    diff = float([l for l in output.splitlines() if l.startswith("SVRDIFF")][0].split()[1])
    assert diff < 1e-6


@pytest.mark.parametrize("method", ["RF", "SVR", "XGBOOST"])
def test_a_planted_qtl_is_found_at_the_right_place(output, method):
    line = [l for l in output.splitlines() if l.startswith("PEAK " + method)][0].split()
    chrom, pos, score, p = line[2], float(line[3]), float(line[4]), float(line[5])
    assert chrom == "1" and abs(pos - 40) <= 12
    assert score >= 1 and p <= 0.05


def test_noise_does_not_give_a_qtl_every_time(output):
    hits = int([l for l in output.splitlines() if l.startswith("NOISEHITS")][0].split()[1])
    assert hits <= 2          # 6 noise runs at alpha 0.05: far from "always" (earlier versions reported QTL in every run)
