# =============================================================================
#  LINKAGE MAP CONSTRUCTION ENGINE
#
#  Follows the published logic of two map-building programs:
#
#  JoinMap 5 (Van Ooijen)            - independence-LOD grouping with "strongest cross link"
#                                      checks, identical loci represented by one locus, order
#                                      search that minimises the sum of adjacent recombination
#                                      frequencies, multipoint recombination-frequency estimation
#                                      by forward-backward EM, probability of each genotype given
#                                      its neighbours (genotyping-error diagnostics).
#  Lep-MAP3 (Rastas 2017)            - pairwise-LOD grouping (SeparateChromosomes2, union-find),
#                                      assignment of left-over markers to the best group
#                                      (JoinSingles2All), maximum-likelihood ordering with an
#                                      HMM that carries a genotyping-error model, improvement of
#                                      the order by permuting small windows of markers
#                                      (polishFast) and distances from the multipoint solution
#                                      through a mapping function.
#
#  Population types: RIL / DH-like (two genotype classes, het calls treated as missing) and F2
#  (three classes, two independent gametes).
# =============================================================================

# -----------------------------------------------------------------------------
#  Genotype coding
#     RIL : a -> 0, b -> 1, anything else (including ab) -> NA
#     F2  : a -> 0, ab -> 1, b -> 2
# -----------------------------------------------------------------------------
lm_geno_codes <- function(df, geno_cols_idx, pop) {
  m <- as.matrix(df[, geno_cols_idx, drop = FALSE])
  codes <- matrix(NA_integer_, nrow(m), ncol(m))
  codes[which(m == "a")] <- 0L
  if (pop == "F2") {
    codes[which(m == "ab")] <- 1L
    codes[which(m == "b")]  <- 2L
  } else {
    codes[which(m == "b")]  <- 1L
  }
  rownames(codes) <- as.character(df$Marker)
  codes
}

lm_xlogy <- function(x, y) ifelse(x > 0, x * log(y), 0)

lm_map_function <- function(r, fun = "kosambi") {
  r <- pmin(pmax(r, 0), 0.4999)
  if (identical(fun, "haldane")) -50 * log(1 - 2 * r)
  else                           25 * log((1 + 2 * r) / (1 - 2 * r))
}

# recombination fraction on the scale of the HMM  ->  recombination fraction per meiosis
#   RIL (selfing): adjacent loci differ with probability R = 2r/(1+2r)  =>  r = R / (2 (1 - R))
lm_r_from_hmm <- function(R, pop) {
  R <- pmin(pmax(R, 0), 0.4999)
  if (pop == "RIL") R / (2 * (1 - R)) else R
}

# -----------------------------------------------------------------------------
#  Pairwise statistics between every row of GA and every row of GB (matrix products, so a
#  whole chromosome is done in one go).
#     R    : recombination frequency on the HMM scale (ML estimate)
#     lod  : JoinMap "independence LOD" = G2 / (2 ln 10), G2 of the genotype contingency table
#            (transformed to one degree of freedom when df > 1). Not affected by segregation
#            distortion and blind to linkage phase, so it is used for grouping.
#     n    : number of individuals scored for both markers
# -----------------------------------------------------------------------------
lm_pairwise_block <- function(GA, GB, pop, min_n = 8L) {
  ind <- function(G, v) { x <- (G == v); x[is.na(x)] <- FALSE; storage.mode(x) <- "double"; x }
  ln10 <- log(10)

  if (pop == "F2") {
    A <- lapply(0:2, function(v) ind(GA, v))
    B <- lapply(0:2, function(v) ind(GB, v))
    N <- vector("list", 9)
    for (i in 1:3) for (j in 1:3) N[[(i - 1) * 3 + j]] <- A[[i]] %*% t(B[[j]])
    n <- Reduce(`+`, N)
    # N[[ (g1) * 3 + g2 + 1 ]]  with g in 0..2
    N00 <- N[[1]]; N01 <- N[[2]]; N02 <- N[[3]]
    N10 <- N[[4]]; N11 <- N[[5]]; N12 <- N[[6]]
    N20 <- N[[7]]; N21 <- N[[8]]; N22 <- N[[9]]
    r <- matrix(0.25, nrow(n), ncol(n))
    for (it in 1:40) {
      q <- 1 - r; pc <- q^2 + r^2
      rec <- 2 * (N02 + N20) + (N01 + N10 + N12 + N21) + N11 * 2 * r^2 / pc
      r_new <- pmin(pmax(rec / (2 * n), 1e-6), 0.4999)
      r_new[!is.finite(r_new)] <- 0.25
      done <- max(abs(r_new - r)) < 1e-7
      r <- r_new
      if (done) break
    }
    # independence G2 on the 3 x 3 table
    rows <- list(N00 + N01 + N02, N10 + N11 + N12, N20 + N21 + N22)
    cols <- list(N00 + N10 + N20, N01 + N11 + N21, N02 + N12 + N22)
    G2 <- 0
    for (i in 1:3) for (j in 1:3) {
      o <- N[[(i - 1) * 3 + j]]
      e <- rows[[i]] * cols[[j]] / n
      G2 <- G2 + 2 * lm_xlogy(o, o / e)
    }
    nr <- (rows[[1]] > 0) + (rows[[2]] > 0) + (rows[[3]] > 0)
    nc <- (cols[[1]] > 0) + (cols[[2]] > 0) + (cols[[3]] > 0)
    df <- (nr - 1) * (nc - 1)
    G1 <- ifelse(df > 1, { e <- exp(-G2 / (2 * pmax(df - 1, 1))); ((4 - e) * e - 3) * (df - 1) + G2 }, G2)
    G1[df < 1] <- 0
    lod <- pmax(G1, 0) / (2 * ln10)
    Rh <- r
  } else {
    A0 <- ind(GA, 0); A1 <- ind(GA, 1); B0 <- ind(GB, 0); B1 <- ind(GB, 1)
    n11 <- A1 %*% t(B1); n00 <- A0 %*% t(B0); n10 <- A1 %*% t(B0); n01 <- A0 %*% t(B1)
    n <- n11 + n00 + n10 + n01
    Rh <- (n10 + n01) / n
    r1 <- n11 + n10; r0 <- n01 + n00; c1 <- n11 + n01; c0 <- n10 + n00
    G2 <- 2 * (lm_xlogy(n11, n11 * n / (r1 * c1)) + lm_xlogy(n10, n10 * n / (r1 * c0)) +
               lm_xlogy(n01, n01 * n / (r0 * c1)) + lm_xlogy(n00, n00 * n / (r0 * c0)))
    G2[!is.finite(G2)] <- 0
    lod <- pmax(G2, 0) / (2 * ln10)
  }
  bad <- n < min_n
  Rh[bad] <- NA_real_
  lod[bad] <- NA_real_
  list(R = Rh, lod = lod, n = n)
}

# recombination frequencies of consecutive markers only (used when a group is too large for a
# full pairwise matrix)
lm_adjacent_R <- function(codes, ord, pop, chunk = 400L) {
  R <- length(ord)
  out <- rep(NA_real_, R - 1)
  i <- 1L
  while (i < R) {
    j <- min(i + chunk - 1L, R - 1L)
    a <- ord[i:j]; b <- ord[(i + 1L):(j + 1L)]
    blk <- lm_pairwise_block(codes[a, , drop = FALSE], codes[b, , drop = FALSE], pop)
    out[i:j] <- blk$R[cbind(seq_along(a), seq_along(a))]
    i <- j + 1L
  }
  out
}

# -----------------------------------------------------------------------------
#  Order search on the recombination-frequency matrix (JoinMap: smallest sum of adjacent
#  recombination frequencies). Greedy construction + 2-opt + segment insertion (with reversal)
#  + a few random restarts (iterated local search).
# -----------------------------------------------------------------------------
lm_path_cost <- function(D, o) sum(D[cbind(o[-length(o)], o[-1])])

lm_two_opt <- function(D, o, max_pass = 100L) {
  n <- length(o)
  if (n < 4) return(o)
  for (pass in seq_len(max_pass)) {
    improved <- FALSE
    for (i in 1:(n - 1)) {
      js <- (i + 1):n
      left_new <- if (i > 1) D[o[i - 1], o[js]] else 0
      left_old <- if (i > 1) D[o[i - 1], o[i]] else 0
      has_right <- (js + 1) <= n
      o_next <- o[pmin(js + 1, n)]
      right_new <- ifelse(has_right, D[cbind(rep(o[i], length(js)), o_next)], 0)
      right_old <- ifelse(has_right, D[cbind(o[js], o_next)], 0)
      delta <- (left_new + right_new) - (left_old + right_old)
      k <- which.min(delta)
      if (delta[k] < -1e-12) {
        j <- js[k]
        o[i:j] <- rev(o[i:j])
        improved <- TRUE
      }
    }
    if (!improved) break
  }
  o
}

lm_or_opt <- function(D, o, max_seg = 3L, max_pass = 50L) {
  n <- length(o)
  if (n < 5) return(o)
  for (pass in seq_len(max_pass)) {
    improved <- FALSE
    for (L in seq_len(min(max_seg, n - 2L))) {
      i <- 1L
      while (i <= n - L + 1L) {
        idx <- i:(i + L - 1L)
        seg <- o[idx]
        rest <- o[-idx]
        m <- length(rest)
        prv <- if (i > 1) o[i - 1L] else NA_integer_
        nxt <- if (i + L - 1L < n) o[i + L] else NA_integer_
        rem_gain <- (if (!is.na(prv)) D[prv, seg[1]] else 0) +
                    (if (!is.na(nxt)) D[seg[L], nxt] else 0) -
                    (if (!is.na(prv) && !is.na(nxt)) D[prv, nxt] else 0)
        dlr <- c(0, D[cbind(rest[-m], rest[-1])], 0)
        fwd <- c(0, D[rest, seg[1]]) + c(D[seg[L], rest], 0) - dlr
        rvs <- c(0, D[rest, seg[L]]) + c(D[seg[1], rest], 0) - dlr
        bf <- which.min(fwd); br <- which.min(rvs)
        if (min(fwd[bf], rvs[br]) < rem_gain - 1e-12) {
          if (fwd[bf] <= rvs[br]) o <- append(rest, seg, after = bf - 1L)
          else                    o <- append(rest, rev(seg), after = br - 1L)
          improved <- TRUE
        }
        i <- i + 1L
      }
    }
    if (!improved) break
  }
  o
}

lm_local_search <- function(D, o) {
  repeat {
    c0 <- lm_path_cost(D, o)
    o <- lm_two_opt(D, o)
    o <- lm_or_opt(D, o)
    if (lm_path_cost(D, o) >= c0 - 1e-12) break
  }
  o
}

# cheapest insertion of markers into an existing path (cost = added sum of recombination frequencies)
lm_insert_unplaced <- function(D, base, unplaced) {
  o <- base
  for (u in unplaced) {
    m <- length(o)
    left <- c(0, D[o, u]); right <- c(D[u, o], 0)
    join <- c(0, D[cbind(o[-m], o[-1])], 0)
    cost <- left + right - join
    t <- which.min(cost)
    o <- append(o, u, after = t - 1L)
  }
  o
}

lm_nn_path <- function(D, start) {
  n <- nrow(D)
  o <- integer(n); o[1] <- start
  left <- rep(TRUE, n); left[start] <- FALSE
  cur <- start
  for (k in 2:n) {
    cand <- which(left)
    nx <- cand[which.min(D[cur, cand])]
    o[k] <- nx; left[nx] <- FALSE; cur <- nx
  }
  o
}

lm_tsp_order <- function(D, start_order, deadline = Inf, seed = 1L) {
  n <- nrow(D)
  if (n <= 2) return(seq_len(n))
  set.seed(seed)
  cands <- list(lm_local_search(D, start_order))
  # greedy constructions from the two ends of the given order, from the marker whose nearest
  # neighbour is farthest (a likely chromosome end) and from random markers
  far <- which.max(apply(D + diag(Inf, n), 1, min))
  starts <- unique(c(start_order[1], start_order[n], far, sample.int(n, min(n, 4L))))
  for (s in starts) {
    if (Sys.time() > deadline) break
    cands[[length(cands) + 1L]] <- lm_local_search(D, lm_nn_path(D, s))
  }
  costs <- vapply(cands, function(o) lm_path_cost(D, o), 0)
  best <- cands[[which.min(costs)]]; best_cost <- min(costs)
  iters <- if (n <= 60) 25L else if (n <= 300) 12L else if (n <= 1000) 5L else 2L
  if (n >= 8) for (it in seq_len(iters)) {
    if (Sys.time() > deadline) break
    cuts <- sort(sample.int(n - 1L, 3L))
    if (length(unique(cuts)) < 3) next
    a <- best[1:cuts[1]]; b <- best[(cuts[1] + 1L):cuts[2]]
    cc <- best[(cuts[2] + 1L):cuts[3]]; d <- if (cuts[3] < n) best[(cuts[3] + 1L):n] else integer(0)
    cand <- lm_local_search(D, c(a, cc, b, d))
    cst <- lm_path_cost(D, cand)
    if (cst < best_cost - 1e-12) { best <- cand; best_cost <- cst }
  }
  best
}

# -----------------------------------------------------------------------------
#  Hidden Markov model of the inheritance of the parental alleles along a chromosome
#     RIL : 2 states (allele of parent 1 / parent 2), switch probability R between neighbours
#     F2  : 4 states = (allele of the paternal gamete, allele of the maternal gamete), each
#           gamete switches with probability r; the observed genotype is the sum.
#  The emission model contains a genotyping-error rate (Lep-MAP3 style).
# -----------------------------------------------------------------------------
lm_hmm_setup <- function(G, pop, eps) {
  eps <- max(eps, 1e-4)
  n <- ncol(G)
  if (pop == "F2") {
    S <- 4L
    sumv <- c(0L, 1L, 1L, 2L)
    em <- matrix(eps / 2, 3, 4)
    for (s in 1:4) em[sumv[s] + 1L, s] <- 1 - eps
    Tf <- function(x) { t2 <- matrix(c(1 - x, x, x, 1 - x), 2); kronecker(t2, t2) }
    a <- c(0L, 0L, 1L, 1L); b <- c(0L, 1L, 0L, 1L)
    Dm <- outer(a, a, `!=`) * 1 + outer(b, b, `!=`) * 1
    pi0 <- rep(0.25, 4); denom <- 2
  } else {
    S <- 2L
    em <- matrix(c(1 - eps, eps, eps, 1 - eps), 2)
    Tf <- function(x) matrix(c(1 - x, x, x, 1 - x), 2)
    Dm <- matrix(c(0, 1, 1, 0), 2)
    pi0 <- c(0.5, 0.5); denom <- 1
  }
  E <- lapply(seq_len(nrow(G)), function(i) {
    g <- G[i, ]
    M <- em[g + 1L, , drop = FALSE]
    M[is.na(g), ] <- 1
    M
  })
  # a %*% T(x) without building T: the hot loop of every search and EM step
  step <- if (pop == "F2") {
    function(a, x) {
      p <- 1 - x; q <- x; pp <- p * p; pq <- p * q; qq <- q * q
      a %*% matrix(c(pp, pq, pq, qq,  pq, pp, qq, pq,  pq, qq, pp, pq,  qq, pq, pq, pp), 4L)
    }
  } else {
    function(a, x) { y <- 1 - x; cbind(a[, 1] * y + a[, 2] * x, a[, 1] * x + a[, 2] * y) }
  }
  list(pop = pop, S = S, n = n, E = E, T = Tf, step = step, TD = function(x) Tf(x) * Dm, pi = pi0,
       pi_mat = matrix(pi0, n, S, byrow = TRUE), denom = denom, obs = !is.na(G), em = em)
}

lm_hmm_pass <- function(md, ord, rv) {
  R <- length(ord); n <- md$n; S <- md$S
  alpha <- vector("list", R); cs <- vector("list", R); bet <- vector("list", R)
  a <- md$pi_mat * md$E[[ord[1]]]
  s <- rowSums(a); alpha[[1]] <- a / s; cs[[1]] <- s
  if (R > 1) for (m in 2:R) {
    a <- md$step(alpha[[m - 1]], rv[m - 1]) * md$E[[ord[m]]]
    s <- rowSums(a); alpha[[m]] <- a / s; cs[[m]] <- s
  }
  ll <- 0
  for (m in seq_len(R)) ll <- ll + sum(log(cs[[m]]))
  bet[[R]] <- matrix(1, n, S)
  if (R > 1) for (m in (R - 1):1)
    bet[[m]] <- md$step(md$E[[ord[m + 1]]] * bet[[m + 1]], rv[m]) / cs[[m + 1]]
  list(alpha = alpha, beta = bet, cs = cs, ll = ll)
}

# Multipoint estimation of the recombination frequencies between neighbours (EM with a
# forward-backward pass; expected number of switches per interval).
lm_hmm_em <- function(md, ord, rv0, max_iter = 40L, tol = 1e-5, deadline = Inf) {
  R <- length(ord)
  rv <- pmin(pmax(rv0, 1e-4), 0.4999)
  ps <- NULL; iters <- 0L
  if (R < 2) return(list(rv = numeric(0), ll = NA_real_, iter = 0L, pass = lm_hmm_pass(md, ord, numeric(0))))
  for (it in seq_len(max_iter)) {
    ps <- lm_hmm_pass(md, ord, rv)
    esw <- numeric(R - 1)
    for (m in 1:(R - 1)) {
      Eb <- md$E[[ord[m + 1]]] * ps$beta[[m + 1]]
      esw[m] <- sum(rowSums((ps$alpha[[m]] %*% md$TD(rv[m])) * Eb) / ps$cs[[m + 1]])
    }
    rv_new <- pmin(pmax(esw / (md$n * md$denom), 1e-6), 0.4999)
    delta <- max(abs(rv_new - rv))
    rv <- rv_new
    iters <- it
    if (delta < tol || Sys.time() > deadline) break
  }
  ps <- lm_hmm_pass(md, ord, rv)
  list(rv = rv, ll = ps$ll, iter = iters, pass = ps)
}

# Probability of each observed genotype given only the genotypes of its neighbours on the map
# (JoinMap "genotype probabilities"): values far below 1 point to genotyping errors.
lm_hmm_genotype_prob <- function(md, ord, rv, ps = NULL) {
  R <- length(ord)
  if (is.null(ps)) ps <- lm_hmm_pass(md, ord, rv)
  P <- matrix(NA_real_, R, md$n)
  for (m in seq_len(R)) {
    left <- if (m == 1) md$pi_mat else md$step(ps$alpha[[m - 1]], rv[m - 1])
    pred <- left * ps$beta[[m]]
    pred <- pred / rowSums(pred)
    p <- rowSums(pred * md$E[[ord[m]]])
    p[!md$obs[ord[m], ]] <- NA_real_
    P[m, ] <- p
  }
  P
}

# -----------------------------------------------------------------------------
#  Polishing by permuting a sliding window of k neighbouring markers and keeping the
#  arrangement with the highest multipoint likelihood (Lep-MAP3 polishFast / JoinMap "ripple").
#  Forward messages up to the window and backward messages after it are cached, so only the
#  window itself is re-evaluated.
# -----------------------------------------------------------------------------
lm_dilate <- function(flag, w) {
  R <- length(flag); out <- rep(FALSE, R)
  for (i in which(flag)) out[max(1L, i - w + 1L):min(R, i + w - 1L)] <- TRUE
  out
}

lm_polish <- function(md, ord, Rmat, k = 4L, max_sweeps = 6L, deadline = Inf, active = NULL) {
  R <- length(ord)
  k <- as.integer(min(k, R))
  touched <- rep(FALSE, R)
  if (R < 3 || k < 2) return(list(ord = ord, moves = 0L, sweeps = 0L, touched = touched))
  Rc <- Rmat; Rc[is.na(Rc)] <- 0.4999; Rc <- pmin(pmax(Rc, 1e-4), 0.4999)
  rp <- function(a, b) Rc[a, b]
  moves <- 0L; sweeps <- 0L
  if (is.null(active)) active <- rep(TRUE, R)
  for (sw in seq_len(max_sweeps)) {
    if (Sys.time() > deadline || !any(active)) break
    sweeps <- sw
    bet <- vector("list", R)
    bet[[R]] <- matrix(1, md$n, md$S)
    for (m in (R - 1):1) {
      b <- md$step(md$E[[ord[m + 1]]] * bet[[m + 1]], rp(ord[m], ord[m + 1]))
      bet[[m]] <- b / rowSums(b)
    }
    changed <- 0L
    a <- NULL
    sweep_touched <- rep(FALSE, R)
    for (p in 1:(R - k + 1)) {
      if (Sys.time() > deadline) break
      win <- ord[p:(p + k - 1L)]
      run_window <- any(active[p:(p + k - 1L)])
      right <- if (p + k <= R) ord[p + k] else NA_integer_
      gamR <- if (!is.na(right)) md$E[[right]] * bet[[p + k]] else NULL
      left <- if (p > 1) ord[p - 1L] else NA_integer_
      best_ls <- -Inf; best_seq <- seq_len(k)
      rec <- function(a_in, prev, used, seq, ls) {
        if (length(seq) == k) {
          if (!is.na(right)) {
            z <- rowSums(md$step(a_in, rp(prev, right)) * gamR)
            ls <- ls + sum(log(z))
          }
          if (ls > best_ls + 1e-7) { best_ls <<- ls; best_seq <<- seq }
          return(invisible())
        }
        for (j in which(!used)) {
          mk <- win[j]
          a2 <- if (is.na(prev)) md$pi_mat * md$E[[mk]] else md$step(a_in, rp(prev, mk)) * md$E[[mk]]
          s <- rowSums(a2)
          u2 <- used; u2[j] <- TRUE
          rec(a2 / s, mk, u2, c(seq, j), ls + sum(log(s)))
        }
      }
      if (run_window) rec(a, left, rep(FALSE, k), integer(0), 0)
      if (run_window && !identical(best_seq, seq_len(k))) {
        sweep_touched[p - 1L + which(win[best_seq] != win)] <- TRUE
        ord[p:(p + k - 1L)] <- win[best_seq]
        changed <- changed + 1L
      }
      a1 <- if (is.null(a)) md$pi_mat * md$E[[ord[p]]]
            else md$step(a, rp(ord[p - 1L], ord[p])) * md$E[[ord[p]]]
      a <- a1 / rowSums(a1)
    }
    moves <- moves + changed
    touched <- touched | sweep_touched
    if (changed == 0L) break
    active <- lm_dilate(sweep_touched, k)
  }
  list(ord = ord, moves = moves, sweeps = sweeps, touched = touched)
}

# Second polishing stage: inside a window of w markers try moving ONE marker to any other place of
# the window and reversing any block of the window (the typical mistakes of an order found from noisy
# pairwise frequencies: a marker a few positions off, a short inverted block). Same likelihood,
# same forward/backward caching as lm_polish().
lm_polish_moves <- function(md, ord, Rmat, w = 8L, max_sweeps = 4L, deadline = Inf, active = NULL) {
  R <- length(ord)
  w <- as.integer(min(w, R))
  touched <- rep(FALSE, R)
  if (w < 3L) return(list(ord = ord, moves = 0L, sweeps = 0L, touched = touched))
  if (is.null(active)) active <- rep(TRUE, R)
  Rc <- Rmat; Rc[is.na(Rc)] <- 0.4999; Rc <- pmin(pmax(Rc, 1e-4), 0.4999)
  rp <- function(a, b) Rc[a, b]
  ident <- seq_len(w)
  cands <- list()
  for (i in ident) for (j in ident) if (i != j) {
    rest <- ident[-i]
    cands[[length(cands) + 1L]] <- append(rest, i, after = j - 1L)
  }
  for (a in 1:(w - 1L)) for (b in (a + 1L):w) {
    p <- ident; p[a:b] <- rev(p[a:b]); cands[[length(cands) + 1L]] <- p
  }
  cands <- unique(cands)
  first_diff <- vapply(cands, function(p) which(p != ident)[1], 1L)
  moves <- 0L; sweeps <- 0L
  for (sw in seq_len(max_sweeps)) {
    if (Sys.time() > deadline || !any(active)) break
    sweeps <- sw
    bet <- vector("list", R)
    bet[[R]] <- matrix(1, md$n, md$S)
    for (m in (R - 1):1) {
      b <- md$step(md$E[[ord[m + 1]]] * bet[[m + 1]], rp(ord[m], ord[m + 1]))
      bet[[m]] <- b / rowSums(b)
    }
    changed <- 0L
    a <- NULL
    sweep_touched <- rep(FALSE, R)
    for (p in 1:(R - w + 1L)) {
      if (Sys.time() > deadline) break
      win <- ord[p:(p + w - 1L)]
      left <- if (p > 1) ord[p - 1L] else NA_integer_
      best <- NULL
      if (any(active[p:(p + w - 1L)])) {
        right <- if (p + w <= R) ord[p + w] else NA_integer_
        gamR <- if (!is.na(right)) md$E[[right]] * bet[[p + w]] else NULL
        # running states along the current arrangement: st[[i + 1]] = after i markers of the window
        st_a <- vector("list", w + 1L); st_ls <- numeric(w + 1L)
        st_a[[1]] <- a; st_ls[1] <- 0
        prev <- left
        for (i in ident) {
          mk <- win[i]
          a2 <- if (is.na(prev)) md$pi_mat * md$E[[mk]] else md$step(st_a[[i]], rp(prev, mk)) * md$E[[mk]]
          s <- rowSums(a2)
          st_a[[i + 1L]] <- a2 / s; st_ls[i + 1L] <- st_ls[i] + sum(log(s)); prev <- mk
        }
        fin <- function(a_in, last, ls) {
          if (is.na(right)) return(ls)
          ls + sum(log(rowSums(md$step(a_in, rp(last, right)) * gamR)))
        }
        id_ls <- fin(st_a[[w + 1L]], win[w], st_ls[w + 1L])
        best_ls <- id_ls
        for (ci in seq_along(cands)) {
          pm <- cands[[ci]]; f <- first_diff[ci]
          a_cur <- st_a[[f]]; ls <- st_ls[f]
          prev <- if (f == 1L) left else win[pm[f - 1L]]
          for (i in f:w) {
            mk <- win[pm[i]]
            a2 <- if (is.na(prev)) md$pi_mat * md$E[[mk]] else md$step(a_cur, rp(prev, mk)) * md$E[[mk]]
            s <- rowSums(a2)
            a_cur <- a2 / s; ls <- ls + sum(log(s)); prev <- mk
          }
          ls <- fin(a_cur, prev, ls)
          if (ls > best_ls + 1e-6) { best_ls <- ls; best <- pm }
        }
      }
      if (!is.null(best)) {
        sweep_touched[p - 1L + which(win[best] != win)] <- TRUE
        ord[p:(p + w - 1L)] <- win[best]
        changed <- changed + 1L
      }
      a1 <- if (is.null(a)) md$pi_mat * md$E[[ord[p]]]
            else md$step(a, rp(ord[p - 1L], ord[p])) * md$E[[ord[p]]]
      a <- a1 / rowSums(a1)
    }
    moves <- moves + changed
    touched <- touched | sweep_touched
    if (changed == 0L) break
    active <- lm_dilate(sweep_touched, w)
  }
  list(ord = ord, moves = moves, sweeps = sweeps, touched = touched)
}

# Both polishing stages in turn, each time looking only where the previous round changed something.
lm_polish_all <- function(md, ord, Rmat, k = 4L, w = 8L, deadline = Inf, rounds = 5L) {
  R <- length(ord)
  active <- rep(TRUE, R); moves <- 0L
  for (round in seq_len(rounds)) {
    if (Sys.time() > deadline) break
    p1 <- lm_polish(md, ord, Rmat, k = k, deadline = deadline, active = active)
    ord <- p1$ord
    m1 <- lm_polish_moves(md, ord, Rmat, w = w, deadline = deadline, active = active | p1$touched)
    ord <- m1$ord
    moves <- moves + p1$moves + m1$moves
    if (p1$moves + m1$moves == 0L) break
    active <- lm_dilate(p1$touched | m1$touched, max(k, w))
  }
  list(ord = ord, moves = moves)
}

# -----------------------------------------------------------------------------
#  Grouping (LOD) -- union-find on pairwise independence LOD (Lep-MAP3 SeparateChromosomes2,
#  JoinMap grouping with the independence LOD)
# -----------------------------------------------------------------------------
lm_denovo_groups <- function(codes, pop, lod_limit, block = 400L) {
  R <- nrow(codes)
  lab <- seq_len(R)
  find <- function(x) { while (lab[x] != x) { lab[x] <<- lab[lab[x]]; x <- lab[x] }; x }
  for (i0 in seq(1L, R, by = block)) {
    rows <- i0:min(R, i0 + block - 1L)
    blk <- lm_pairwise_block(codes[rows, , drop = FALSE], codes, pop)$lod
    hit <- which(!is.na(blk) & blk >= lod_limit, arr.ind = TRUE)
    if (!nrow(hit)) next
    gi <- rows[hit[, 1]]; gj <- hit[, 2]
    keep <- gj > gi
    gi <- gi[keep]; gj <- gj[keep]
    if (!length(gi)) next
    comp <- vapply(seq_len(R), find, 0L)
    pr <- unique(cbind(pmin(comp[gi], comp[gj]), pmax(comp[gi], comp[gj])))
    pr <- pr[pr[, 1] != pr[, 2], , drop = FALSE]
    for (e in seq_len(nrow(pr))) {
      a <- find(pr[e, 1]); b <- find(pr[e, 2])
      if (a != b) lab[max(a, b)] <- min(a, b)
    }
  }
  vapply(seq_len(R), find, 0L)
}

# strongest linkage of each marker to every group (via up to n_anchor anchor markers per group):
# JoinMap "strongest cross link"; Lep-MAP3 JoinSingles2All
lm_group_scan <- function(codes, chr, pop, support_lod = 5, exclude = NULL, n_anchor = 40L, block = 2000L) {
  groups <- unique(chr[!is.na(chr)])
  anchors <- lapply(groups, function(g) {
    idx <- setdiff(which(chr == g), exclude)
    if (length(idx) <= n_anchor) idx else idx[unique(round(seq(1, length(idx), length.out = n_anchor)))]
  })
  names(anchors) <- groups
  all_anc <- unlist(anchors, use.names = FALSE)
  anc_grp <- rep(groups, lengths(anchors))
  R <- nrow(codes)
  best <- matrix(NA_real_, R, length(groups), dimnames = list(NULL, groups))
  support <- matrix(0L, R, length(groups), dimnames = list(NULL, groups))
  for (i0 in seq(1L, R, by = block)) {
    rows <- i0:min(R, i0 + block - 1L)
    lod <- lm_pairwise_block(codes[rows, , drop = FALSE], codes[all_anc, , drop = FALSE], pop)$lod
    self_col <- match(rows, all_anc)                 # a marker is not evidence for its own group
    has_self <- which(!is.na(self_col))
    if (length(has_self)) lod[cbind(has_self, self_col[has_self])] <- NA_real_
    for (g in seq_along(groups)) {
      cols <- which(anc_grp == groups[g])
      sub <- lod[, cols, drop = FALSE]
      mx <- suppressWarnings(apply(sub, 1, max, na.rm = TRUE))
      mx[!is.finite(mx)] <- NA_real_
      best[rows, g] <- mx
      support[rows, g] <- rowSums(sub >= support_lod, na.rm = TRUE)
    }
  }
  list(best = best, support = support, n_anchor = setNames(lengths(anchors), groups))
}

# What the scan says about marker u (own = its chromosome label in the file).
#   A group counts as linked only when at least two of its anchor markers (one if it has one)
#   reach the LOD limit, so a single mislabelled marker cannot pull others out of their group.
lm_group_call <- function(scan, u, own, unknown, lod_limit) {
  b <- scan$best[u, ]; s <- scan$support[u, ]
  need <- pmin(2L, scan$n_anchor)
  valid <- !is.na(b) & b >= lod_limit & s >= need
  top_any <- suppressWarnings(max(b, na.rm = TRUE))
  if (unknown) {
    if (!any(valid)) return(list(action = "none", lod = top_any, target = NA_character_, own = NA_real_))
    vb <- ifelse(valid, b, -Inf)
    t <- names(which.max(vb))
    rest <- b; rest[t] <- NA
    second <- suppressWarnings(max(rest, na.rm = TRUE))
    if (!is.finite(second)) second <- -Inf
    if (b[t] - second >= 1) return(list(action = "assign", lod = b[t], target = t, own = NA_real_))
    return(list(action = "none", lod = b[t], target = t, own = NA_real_))
  }
  own_lod <- if (own %in% names(b) && !is.na(b[own])) b[[own]] else -Inf
  cand <- valid; if (own %in% names(cand)) cand[own] <- FALSE
  if (any(cand)) {
    t <- names(which.max(ifelse(cand, b, -Inf)))
    if (b[t] > own_lod) return(list(action = "move", lod = b[[t]], target = t, own = own_lod))
  }
  list(action = "ok", lod = top_any, target = NA_character_, own = own_lod)
}

# -----------------------------------------------------------------------------
#  Map of one linkage group
# -----------------------------------------------------------------------------
lm_map_group <- function(codes, start_rank, pop, cfg, label = "", phys = NULL) {
  # codes: reps x individuals (rows in the GIVEN order, start_rank = position in the given order)
  R <- nrow(codes)
  info <- list(label = label, R = R)
  deadline <- Sys.time() + cfg$max_seconds
  md <- lm_hmm_setup(codes, pop, cfg$genotyping_error)
  ord <- seq_len(R)                       # rows are already in the given order
  notes <- character(0)

  full <- R <= cfg$max_group_markers
  Rmat <- NULL; Lmat <- NULL
  if (full && R >= 2) {
    pw <- lm_pairwise_block(codes, codes, pop)
    Rmat <- pw$R; Lmat <- pw$lod
    diag(Rmat) <- 0; diag(Lmat) <- NA_real_
  }

  # --- weakly linked markers (JoinMap: ungrouped; Lep-MAP3: singles) -----------------------
  dropped <- integer(0); drop_detail <- character(0)
  if (!is.null(Lmat) && R >= 4 && cfg$min_lod > 0) {
    mx <- suppressWarnings(apply(Lmat, 1, max, na.rm = TRUE)); mx[!is.finite(mx)] <- 0
    weak <- which(mx < cfg$min_lod)
    if (length(weak) && length(weak) < R - 2) {
      dropped <- weak
      drop_detail <- sprintf("Weakly linked: strongest LOD with any marker of the group is %.1f (< %.1f)", mx[weak], cfg$min_lod)
    }
  }
  keep_idx <- setdiff(seq_len(R), dropped)

  method <- cfg$map_method
  Rk <- if (!is.null(Rmat)) Rmat[keep_idx, keep_idx, drop = FALSE] else NULL
  o_local <- seq_along(keep_idx)           # order among kept markers, in indices of keep_idx
  phys_k <- if (!is.null(phys)) phys[keep_idx] else rep(NA_real_, length(keep_idx))
  n_phys <- sum(!is.na(phys_k))
  phys_ok <- n_phys >= 5 && n_phys >= 0.7 * length(phys_k)

  # Markers whose name carries no position (gene markers, array IDs) would sit at the end of a
  # position-sorted list. Put each one where it fits the data best instead (JoinMap adds loci one
  # by one at their best position).
  if (!is.null(Rk) && phys_ok && n_phys < length(phys_k) && !identical(cfg$marker_order, "file")) {
    placed <- which(!is.na(phys_k))
    if (!is.unsorted(phys_k[placed])) {
      Dg <- Rk; Dg[is.na(Dg)] <- 0.5; Dg <- pmin(Dg, 0.5); diag(Dg) <- 0
      o_local <- lm_insert_unplaced(Dg, placed, which(is.na(phys_k)))
    }
  }

  if (method == "ml" && !is.null(Rk) && length(keep_idx) >= 4) {
    D <- Rk; D[is.na(D)] <- 0.5; diag(D) <- 0
    D <- pmin(D, 0.5)                      # R > 0.5 means phase problem, not "farther than unlinked"
    cand_tsp <- lm_tsp_order(D, o_local, deadline = deadline)
    # direction of the map: along the genome when the marker names give positions, else along the file
    if (length(cand_tsp) > 2) {
      ref <- if (phys_ok) phys_k else start_rank[keep_idx]
      rc <- suppressWarnings(cor(seq_along(cand_tsp), ref[cand_tsp], method = "spearman", use = "complete.obs"))
      if (is.finite(rc) && rc < 0) cand_tsp <- rev(cand_tsp)
    }
    mdk <- lm_hmm_setup(codes[keep_idx, , drop = FALSE], pop, cfg$genotyping_error)
    cands <- list(data = cand_tsp)
    if (length(keep_idx) <= cfg$max_given_candidates) cands$given <- o_local
    res <- lapply(cands, function(o) {
      w_eff <- if (length(o) > 700L) min(cfg$move_window, 6L) else cfg$move_window   # keep very large groups tractable
      pol <- lm_polish_all(mdk, o, Rk, k = cfg$polish_window, w = w_eff, deadline = deadline)
      rv0 <- vapply(seq_len(length(o) - 1L), function(i) { x <- Rk[pol$ord[i], pol$ord[i + 1L]]; if (is.na(x)) 0.4999 else x }, 0)
      ll <- lm_hmm_pass(mdk, pol$ord, pmin(pmax(rv0, 1e-4), 0.4999))$ll
      list(ord = pol$ord, ll = ll, moves = pol$moves, cost = lm_path_cost(D, pol$ord))
    })
    # The given order is a prior (physical positions, a previous map): it is replaced by the
    # data-driven order only when the likelihood is clearly higher (cfg$prefer_given_lod LOD units).
    pick <- "data"
    if (!is.null(res$given)) {
      gain_lod <- (res$data$ll - res$given$ll) / log(10)
      pick <- if (gain_lod >= cfg$prefer_given_lod || !phys_ok && gain_lod > 0) "data" else "given"
      info$gain_lod <- gain_lod
    }
    chosen <- res[[pick]]
    info$ll_data <- res$data$ll
    info$ll_given <- if (!is.null(res$given)) res$given$ll else NA_real_
    info$polish_moves <- chosen$moves
    info$order_source <- if (pick == "data") "data-driven order (rec. frequency search + likelihood polish)" else "given/physical order (data-driven order not clearly better)"
    o_local <- chosen$ord
    if (Sys.time() > deadline) notes <- c(notes, "time budget reached while ordering; order may not be fully optimised")
  } else if (method == "ml" && is.null(Rk)) {
    notes <- c(notes, sprintf("%d unique markers exceed the %d limit for a full pairwise search; kept the given order", R, cfg$max_group_markers))
    info$order_source <- "given order (group too large to search)"
  } else {
    info$order_source <- "given order (kept)"
    # still polish a fixed-order map? No: 'fixed' means exactly the given order.
  }

  ord_global <- keep_idx[o_local]          # indices into codes rows

  # --- poor-fit markers (large jump in map length because of one marker) --------------------
  fit_removed <- integer(0); fit_detail <- character(0)
  if (!is.null(Rmat) && cfg$fit_threshold_cm > 0 && length(ord_global) >= 5) {
    repeat {
      n_o <- length(ord_global)
      if (n_o < 5) break
      rr <- function(a, b) lm_r_from_hmm(Rmat[cbind(a, b)], pop)
      r12 <- rr(ord_global[1:(n_o - 2)], ord_global[2:(n_o - 1)])
      r23 <- rr(ord_global[2:(n_o - 1)], ord_global[3:n_o])
      r13 <- rr(ord_global[1:(n_o - 2)], ord_global[3:n_o])
      okr <- is.finite(r12) & is.finite(r23) & is.finite(r13) & r12 < 0.45 & r23 < 0.45 & r13 < 0.45
      delta <- ifelse(okr, lm_map_function(r12, cfg$mapping_function) + lm_map_function(r23, cfg$mapping_function) -
                           lm_map_function(r13, cfg$mapping_function), 0)
      w <- which.max(delta)
      if (!length(w) || delta[w] <= cfg$fit_threshold_cm) break
      if (length(fit_removed) >= max(1, floor(0.2 * R))) break
      pos <- w + 1L
      fit_removed <- c(fit_removed, ord_global[pos])
      fit_detail <- c(fit_detail, sprintf("Poor fit: placing it between its neighbours adds %.1f cM to the map (limit %.1f cM)", delta[w], cfg$fit_threshold_cm))
      ord_global <- ord_global[-pos]
    }
  }

  # --- multipoint estimation on the final order -------------------------------------------
  Rn <- length(ord_global)
  mdf <- lm_hmm_setup(codes[ord_global, , drop = FALSE], pop, cfg$genotyping_error)
  ordf <- seq_len(Rn)
  if (Rn >= 2) {
    rv0 <- if (!is.null(Rmat)) vapply(seq_len(Rn - 1L), function(i) { x <- Rmat[ord_global[i], ord_global[i + 1L]]; if (is.na(x)) 0.4999 else x }, 0)
           else lm_adjacent_R(codes, ord_global, pop)
    rv0[is.na(rv0)] <- 0.4999
    em <- lm_hmm_em(mdf, ordf, rv0, deadline = deadline + 30)
    rv <- em$rv; ll <- em$ll
    P <- lm_hmm_genotype_prob(mdf, ordf, rv, em$pass)
  } else {
    rv <- numeric(0); ll <- NA_real_
    P <- matrix(NA_real_, Rn, ncol(codes))
  }
  r_map <- lm_r_from_hmm(rv, pop)
  cm <- lm_map_function(r_map, cfg$mapping_function)
  broken <- r_map >= cfg$recombination_threshold
  cm[broken] <- 50
  pos <- c(0, cumsum(cm))

  info$ll <- ll
  info$n_breaks <- sum(broken)
  info$em_iter <- if (Rn >= 2) em$iter else 0L
  info$notes <- notes
  list(
    ord = ord_global, pos = pos, rv = rv, r_map = r_map, cm = cm, broken = broken,
    dropped = dropped, drop_detail = drop_detail,
    fit_removed = fit_removed, fit_detail = fit_detail,
    suspect = rowSums(P < cfg$suspect_p, na.rm = TRUE),
    scored = rowSums(!is.na(P)),
    info = info
  )
}

# -----------------------------------------------------------------------------
#  Whole data set
# -----------------------------------------------------------------------------
lm_build_maps <- function(df, pop, cfg, removal_log) {
  geno_cols_idx <- setdiff(3:ncol(df), which(colnames(df) == ".input_order"))
  codes <- lm_geno_codes(df, geno_cols_idx, pop)
  N <- nrow(df)
  markers <- as.character(df$Marker)
  chr <- as.character(df$Chromosome)
  chr_file <- chr                          # chromosome labels as given in the file
  in_ord <- df$.input_order
  # physical position (bp) carried by the marker name, used for orientation and as the prior order
  phys_all <- vapply(seq_len(N), function(i) parse_marker_position(markers[i], chr_file[i]), 0)
  report <- character(0)
  log_add <- function(idx, reason, detail) {
    if (!length(idx)) return(invisible())
    removal_log <<- rbind(removal_log, data.frame(
      Marker = markers[idx], Chromosome = chr[idx], Reason = reason,
      Detail = rep_len(detail, length(idx)), stringsAsFactors = FALSE))
  }
  alive <- rep(TRUE, N)

  # ---- 1. identical markers: one representative per set (JoinMap "identicals") --------------
  pat <- apply(codes, 1, function(x) paste(ifelse(is.na(x), "-", x), collapse = ""))
  grouping_mode <- cfg$grouping
  key <- if (grouping_mode == "lod") pat else paste(chr, pat, sep = "|")
  first_of <- match(key, key)
  is_rep <- first_of == seq_len(N)
  rep_of <- first_of
  n_ident <- sum(!is_rep)

  rep_idx <- which(is_rep)
  codes_r <- codes[rep_idx, , drop = FALSE]

  # ---- 2. linkage groups ---------------------------------------------------------------------
  chr_r <- chr[rep_idx]
  unknown <- is.na(chr_r) | chr_r %in% c("Unknown", "NA", "")
  scan_note <- NULL
  if (grouping_mode == "lod") {
    lab <- lm_denovo_groups(codes_r, pop, cfg$group_lod)
    sizes <- table(lab)
    ordg <- names(sizes)[order(-as.integer(sizes), as.integer(names(sizes)))]
    keepg <- ordg[as.integer(sizes[ordg]) >= cfg$min_group_size]
    newlab <- setNames(seq_along(keepg), keepg)
    chr_new <- ifelse(as.character(lab) %in% keepg, as.character(newlab[as.character(lab)]), NA_character_)
    # relate each new group to the chromosome label that most of its markers carry in the file
    lines <- vapply(seq_along(keepg), function(i) {
      mem <- which(chr_new == as.character(i))
      tb <- sort(table(chr_r[mem]), decreasing = TRUE)
      sprintf("    LG%-3d %5d unique markers   mostly labelled chromosome %s (%.0f%%)", i, length(mem), names(tb)[1], 100 * tb[1] / length(mem))
    }, "")
    report <- c(report, sprintf("  Grouping            : de novo, independence LOD >= %.1f -> %d linkage groups", cfg$group_lod, length(keepg)), lines)
    unl <- which(is.na(chr_new))
    if (length(unl)) {
      for (u in unl) {
        members <- which(rep_of == rep_idx[u]); alive[members] <- FALSE
        log_add(members, "Not linked to any group", sprintf("No linkage with LOD >= %.1f to a group of at least %d unique markers", cfg$group_lod, cfg$min_group_size))
      }
    }
    chr_r <- chr_new
    chr[rep_idx] <- chr_r
    for (i in seq_len(N)) chr[i] <- chr[rep_of[i]]
    df$Chromosome <- chr
  } else {
    do_scan <- cfg$group_check != "off" && length(unique(chr_r[!unknown])) >= 2
    if (do_scan || any(unknown)) {
      grp_ok <- chr_r; grp_ok[unknown] <- NA_character_
      scan <- lm_group_scan(codes_r, grp_ok, pop, support_lod = cfg$group_lod)
      # Second pass without the markers the first pass says belong elsewhere: a mislabelled marker
      # used as an anchor of the wrong group would otherwise make its true neighbours look misassigned.
      calls1 <- lapply(seq_len(nrow(scan$best)), function(u) lm_group_call(scan, u, chr_r[u], unknown[u], cfg$group_lod))
      bad1 <- which(vapply(calls1, function(cl) cl$action == "move", TRUE))
      if (length(bad1)) scan <- lm_group_scan(codes_r, grp_ok, pop, support_lod = cfg$group_lod, exclude = bad1)
      moved <- 0L; assigned <- 0L; unassigned <- 0L; flagged <- character(0)
      for (u in seq_len(nrow(scan$best))) {
        mem <- which(rep_of == rep_idx[u])
        cl <- lm_group_call(scan, u, chr_r[u], unknown[u], cfg$group_lod)
        if (unknown[u]) {
          if (cl$action == "assign") {
            chr[mem] <- cl$target; chr_r[u] <- cl$target; assigned <- assigned + 1L
          } else {
            alive[mem] <- FALSE; unassigned <- unassigned + 1L
            log_add(mem, "Unassigned (no chromosome)",
                    if (is.finite(cl$lod)) sprintf("No convincing group: best LOD %.1f (needs >= %.1f from at least two markers of one group, and a margin of 1 over the next group)", cl$lod, cfg$group_lod)
                    else "No chromosome in the file and no linkage to any group")
          }
        } else if (do_scan && cl$action == "move") {
          moved <- moved + 1L
          flagged <- c(flagged, sprintf("%s (chr %s -> strongest link to chr %s, LOD %.1f vs %s)", markers[rep_idx[u]], chr_r[u], cl$target, cl$lod,
                                        if (is.finite(cl$own)) sprintf("%.1f within own group", cl$own) else "none within own group"))
          if (cfg$group_check == "remove") {
            alive[mem] <- FALSE
            log_add(mem, "Linked to another chromosome",
                    sprintf("Strongest linkage is to chromosome %s (LOD %.1f), own chromosome %s (%s)", cl$target, cl$lod, chr_r[u],
                            if (is.finite(cl$own)) sprintf("LOD %.1f", cl$own) else "no linkage"))
          }
        }
      }
      report <- c(report, sprintf("  Group check         : %d marker(s) linked more strongly to another chromosome (%s)%s%s",
                                  moved, if (cfg$group_check == "remove") "removed" else "reported only",
                                  if (assigned + unassigned > 0) sprintf("; chromosome-less markers: %d assigned, %d left out", assigned, unassigned) else "",
                                  ""))
      if (length(flagged) && cfg$group_check == "report") report <- c(report, paste0("    ", head(flagged, 12)), if (length(flagged) > 12) "    ...")
      df$Chromosome <- chr
    } else {
      report <- c(report, "  Group check         : skipped (one chromosome or switched off)")
    }
  }

  # members of a removed representative are removed too; a representative that is itself
  # removed leaves its other members without a representative -> they are dropped as well
  for (i in which(!is_rep)) if (!alive[rep_of[i]]) alive[i] <- FALSE

  # ---- 3. map every linkage group ------------------------------------------------------------
  positions_list <- list(); rec_rows <- list(); fit_rows <- list(); sum_rows <- list()
  rank_out <- rep(NA_integer_, N)
  chroms <- unique(chr[alive])
  cn <- suppressWarnings(as.numeric(chroms)); chroms <- if (all(!is.na(cn))) chroms[order(cn)] else sort(chroms)
  grp_lines <- character(0)
  for (cc in chroms) {
    members <- which(alive & chr == cc)
    reps <- members[is_rep[members]]
    reps <- reps[order(in_ord[reps])]
    if (length(reps) == 0) next
    start_rank <- seq_along(reps)
    t0 <- Sys.time()
    ph <- phys_all[reps]
    miss_ph <- which(is.na(ph))
    if (length(miss_ph)) ph[miss_ph] <- vapply(reps[miss_ph], function(i) parse_marker_position(markers[i], cc), 0)
    g <- lm_map_group(codes[reps, , drop = FALSE], start_rank, pop, cfg, label = cc, phys = ph)
    # markers dropped from the group
    drop_all <- c(g$dropped, g$fit_removed)
    if (length(g$dropped)) for (k in seq_along(g$dropped)) {
      mem <- which(rep_of == reps[g$dropped[k]] & alive); alive[mem] <- FALSE
      log_add(mem, "Weakly linked", g$drop_detail[k])
    }
    if (length(g$fit_removed)) for (k in seq_along(g$fit_removed)) {
      mem <- which(rep_of == reps[g$fit_removed[k]] & alive); alive[mem] <- FALSE
      log_add(mem, "Poor fit in map", g$fit_detail[k])
    }
    rep_ord <- reps[g$ord]
    L <- length(rep_ord)
    pos_rep <- g$pos
    n_ident_grp <- 0L
    for (k in seq_len(L)) {
      mem <- which(rep_of == rep_ord[k] & alive)
      mem <- mem[order(in_ord[mem])]
      for (mm in mem) positions_list[[markers[mm]]] <- pos_rep[k]
      rank_out[mem] <- k
      fit_rows[[length(fit_rows) + 1L]] <- data.frame(
        Marker = markers[rep_ord[k]], Chromosome = cc, Position = round(pos_rep[k], 3),
        Markers_At_Position = length(mem),
        Scored = g$scored[k], Suspect_Genotypes = g$suspect[k],
        Suspect_Pct = if (g$scored[k] > 0) round(100 * g$suspect[k] / g$scored[k], 2) else NA_real_,
        stringsAsFactors = FALSE)
      n_ident_grp <- n_ident_grp + length(mem) - 1L
    }
    if (L >= 2) for (k in seq_len(L - 1L))
      rec_rows[[length(rec_rows) + 1L]] <- data.frame(
        Marker1 = markers[rep_ord[k]], Marker2 = markers[rep_ord[k + 1L]],
        RecombinationFraction = g$r_map[k], MapDistance_cM = g$cm[k],
        Chromosome = cc, stringsAsFactors = FALSE)
    len <- if (L >= 2) tail(pos_rep, 1) else 0
    sum_rows[[length(sum_rows) + 1L]] <- data.frame(
      Chromosome = cc, Markers = sum(alive & chr == cc), Unique_Positions = L,
      Length_cM = round(len, 2), Mean_Spacing_cM = if (L >= 2) round(len / (L - 1L), 2) else NA_real_,
      Largest_Gap_cM = if (L >= 2) round(max(g$cm), 2) else NA_real_,
      Linkage_Breaks = g$info$n_breaks,
      Log_Likelihood = round(g$info$ll, 2),
      Suspect_Genotypes = sum(g$suspect), Order_Source = g$info$order_source,
      stringsAsFactors = FALSE)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    grp_lines <- c(grp_lines, sprintf("    chr %-4s %4d markers  %4d positions  %7.1f cM  logL %.1f  [%s, %.1fs]%s%s",
                                      cc, sum(alive & chr == cc), L, len, g$info$ll, g$info$order_source, el,
                                      if (g$info$n_breaks > 0) sprintf("  WARNING: %d interval(s) with r >= %.2f (unlinked) counted as 50 cM", g$info$n_breaks, cfg$recombination_threshold) else "",
                                      if (length(g$info$notes)) paste0("  NOTE: ", paste(g$info$notes, collapse = "; ")) else ""))
    if (!is.null(g$info$ll_given) && is.finite(g$info$ll_given) && is.finite(g$info$ll_data))
      grp_lines <- c(grp_lines, sprintf("             likelihood: data-driven order %.1f vs given order %.1f (difference %+.1f)",
                                        g$info$ll_data, g$info$ll_given, g$info$ll_data - g$info$ll_given))
  }

  list(
    keep = alive, positions_list = positions_list,
    recombination_rows = rec_rows, removal_log = removal_log,
    marker_fit = if (length(fit_rows)) do.call(rbind, fit_rows) else data.frame(),
    chrom_summary = if (length(sum_rows)) do.call(rbind, sum_rows) else data.frame(),
    map_rank = rank_out, n_identical = n_ident,
    chromosome = chr,
    report = paste(c(
      sprintf("  Map method          : %s", switch(cfg$map_method,
              ml = "optimised order (JoinMap/Lep-MAP3 logic) + multipoint ML distances",
              fixed = "given order kept, multipoint ML distances",
              cfg$map_method)),
      sprintf("  Mapping function    : %s", cfg$mapping_function),
      sprintf("  Genotyping error    : %.3f (HMM emission model)", cfg$genotyping_error),
      sprintf("  Identical markers   : %d markers share the position of an identical marker (one representative each)", n_ident),
      report, "  Linkage groups:", grp_lines), collapse = "\n")
  )
}
