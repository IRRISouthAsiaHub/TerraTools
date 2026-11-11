"""
ocs_optimizer.py — the plan of crosses of optimal cross selection (OCS).

The candidates are the genotypes of the marker data. Every cross has a value U (its usefulness, from
ocs_usefulness.py) and every pair of candidates a coancestry C (half of the relationship matrix calculated from the
allele dosage). OCS chooses a plan of `n_crosses` crosses, each parent used in at most `max_per_parent` of them. A plan
gives every candidate a contribution  x_i = n_i / (2 n_crosses)  (n_i = crosses it is in) and is scored by

    expected gain            g = mean of U over the crosses of the plan   (mean usefulness)
    group coancestry         c = x'Cx                                     (diversity kept in the next generation)
    progeny inbreeding       f = mean of C over the mates                 (how related the two parents of a cross are)

The balance between gain and diversity is set with one "penalty degree": 0 degrees maximises gain only,
90 degrees minimises coancestry only. Gain and coancestry are rescaled to a common 0..1 range and the plan
maximises

    J = cos(angle) * gain_scaled - sin(angle) * coancestry_scaled        coancestry = c + f

The plan is searched by simulated annealing over the list of crosses. Running the optimisation over a range of
angles gives the gain / coancestry frontier, from which the angle is chosen.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

ANGLES = [0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90]


class OCSError(Exception):
    pass


# ------------------------------------------------------------
# PLAN SCORING
# ------------------------------------------------------------

class Problem:
    """`a` = genomic value of each candidate, `C` = coancestry matrix, `U` (n x n) = usefulness of every cross."""

    def __init__(self, a, C, n_crosses, max_per_parent, U):
        self.a, self.C, self.U = a, C, U
        self.n = len(a)
        self.nc = int(n_crosses)
        self.m = int(max_per_parent)
        k = max(2, math.ceil(2 * self.nc / self.m))
        top = np.argsort(-a)[:k]
        # reference points for the scaling: random mating of everybody, and the best crosses whatever the limit per parent
        vals = U[np.triu_indices(self.n, 1)]
        self.g0 = float(vals.mean())
        self.g1 = float(np.sort(vals)[-self.nc:].mean())
        self.c0 = float(C.mean())
        self.c1 = float(C[np.ix_(top, top)].mean())
        self.gscale = (self.g1 - self.g0) or 1.0
        self.cscale = (self.c1 - self.c0) or 1.0

    def score(self, plan):
        """(gain, group coancestry, progeny inbreeding) of a plan = int array (n_crosses, 2)."""
        cnt = np.bincount(plan.ravel(), minlength=self.n)
        x = cnt / (2.0 * self.nc)
        idx = np.nonzero(cnt)[0]
        g = float(self.U[plan[:, 0], plan[:, 1]].mean())
        c = float(x[idx] @ self.C[np.ix_(idx, idx)] @ x[idx])
        f = float(self.C[plan[:, 0], plan[:, 1]].mean())
        return g, c, f

    def objective(self, g, c, f, theta):
        gs = (g - self.g0) / self.gscale
        cs = (c + f - self.c0) / self.cscale
        return math.cos(theta) * gs - math.sin(theta) * cs


def _initial_plan(prob, rng):
    """The best parents, each used `max_per_parent` times, paired at random without selfing or repeats."""
    order = list(np.argsort(-prob.a))
    slots = []
    for p in order:
        slots += [p] * prob.m
        if len(slots) >= 2 * prob.nc:
            break
    slots = np.array(slots[:2 * prob.nc])
    for _ in range(500):
        rng.shuffle(slots)
        pairs = slots.reshape(-1, 2)
        keys = {tuple(sorted(p)) for p in pairs}
        if len(keys) == prob.nc and not (pairs[:, 0] == pairs[:, 1]).any():
            return pairs.copy()
    # fall back: walk through distinct pairs of the best parents
    pairs = []
    top = order[:max(2, math.ceil(2 * prob.nc / prob.m) + 1)]
    for i in range(len(top)):
        for j in range(i + 1, len(top)):
            pairs.append((top[i], top[j]))
            if len(pairs) == prob.nc:
                return np.array(pairs)
    raise OCSError("There are not enough candidates for that many different crosses.")


def optimise(prob, theta_deg, rng, start=None, iters=60000, t0=0.05, t1=0.0005):
    """Simulated annealing over the crossing plan; returns (plan, gain, coancestry, inbreeding, objective)."""
    theta = math.radians(theta_deg)
    plan = (start if start is not None else _initial_plan(prob, rng)).copy()
    n, nc, m = prob.n, prob.nc, prob.m
    C, a = prob.C, prob.a
    d = 1.0 / (2.0 * nc)

    cnt = np.bincount(plan.ravel(), minlength=n).astype(np.int64)
    x = cnt * d
    Cx = C @ x
    U = prob.U
    g = float(U[plan[:, 0], plan[:, 1]].mean())
    c = float(x @ Cx)
    f_sum = float(C[plan[:, 0], plan[:, 1]].sum())
    pairs = {tuple(sorted(p)) for p in plan.tolist()}
    cos_t, sin_t = math.cos(theta), math.sin(theta)

    def J(g_, c_, fs_):
        return cos_t * (g_ - prob.g0) / prob.gscale - sin_t * (c_ + fs_ / nc - prob.c0) / prob.cscale

    cur = J(g, c, f_sum)
    best, best_plan = cur, plan.copy()
    ratio = (t1 / t0) ** (1.0 / max(1, iters))
    temp = t0
    top = np.argsort(-a)[:max(10, min(n, 3 * math.ceil(2 * nc / m)))]
    rand_k = rng.integers(0, nc, iters)
    rand_s = rng.integers(0, 2, iters)
    rand_j = rng.integers(0, n, iters)
    rand_biased = rng.random(iters) < 0.3
    rand_top = rng.integers(0, len(top), iters)
    rand_u = rng.random(iters)

    for it in range(iters):
        k, s = rand_k[it], rand_s[it]
        j = int(top[rand_top[it]]) if rand_biased[it] else int(rand_j[it])
        i = int(plan[k, s])
        o = int(plan[k, 1 - s])
        temp *= ratio
        if j == i or j == o or cnt[j] >= m:
            continue
        key = (j, o) if j < o else (o, j)
        if key in pairs:
            continue
        # contribution moves from i to j: x_j += d, x_i -= d
        g2 = g + (U[j, o] - U[i, o]) / nc
        c2 = c + 2 * d * (Cx[j] - Cx[i]) + d * d * (C[j, j] + C[i, i] - 2 * C[i, j])
        fs2 = f_sum + C[j, o] - C[i, o]
        new = J(g2, c2, fs2)
        if new >= cur or rand_u[it] < math.exp((new - cur) / temp):
            old_key = (i, o) if i < o else (o, i)
            pairs.discard(old_key)
            pairs.add(key)
            plan[k, s] = j
            cnt[i] -= 1
            cnt[j] += 1
            Cx += d * (C[:, j] - C[:, i])
            g, c, f_sum, cur = g2, c2, fs2, new
            if cur > best:
                best, best_plan = cur, plan.copy()

    g, c, f = prob.score(best_plan)
    return best_plan, g, c, f, prob.objective(g, c, f, theta)


# ------------------------------------------------------------
# PUBLIC ENTRYPOINT
# ------------------------------------------------------------

def run_ocs(config, results_dir):
    """Marker data (allele dosage matrix, marker effects, genetic map or VCF) -> usefulness of every cross -> the plan.
    Writes ocs_crosses.csv, ocs_contributions.csv, ocs_frontier.csv, ocs_frontier.png and returns the 'ocs' part of the result."""
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    from tools.cross_selection.backend.ocs_usefulness import run_usefulness
    return run_usefulness(config, results_dir)


def plan_and_write(names, a, C, config, results_dir, U, pairinfo, gain_label="mean usefulness of the crosses", crit_label="Genomic"):
    """Optimise the plan over the penalty degrees and write the three CSV files and the frontier plot.
    U / pairinfo: usefulness of every cross (n x n) and the mean and standard deviation of every pair."""
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    nc = int(config.get("n_crosses") or 20)
    m = int(config.get("max_crosses_per_parent") or 4)
    angle = float(config.get("penalty_degree", 30))
    if nc < 1 or m < 1:
        raise OCSError("The number of crosses and the maximum crosses per parent must be at least 1.")
    if not 0 <= angle <= 90:
        raise OCSError("The penalty degree must be between 0 and 90.")
    n = len(names)
    if n * m < 2 * nc:
        raise OCSError(f"{n} candidates cannot make {nc} crosses with at most {m} crosses per parent "
                       f"(needs at least {math.ceil(2 * nc / m)} candidates).")
    if nc > n * (n - 1) // 2:
        raise OCSError(f"There are only {n * (n - 1) // 2} different crosses among {n} candidates.")

    rng = np.random.default_rng(int(config.get("seed", 1)))
    prob = Problem(a, C, nc, m, U)

    # frontier: the same problem at a range of angles, each starting from the previous plan
    angles = sorted(set(ANGLES) | {angle})
    frontier, plans = [], {}
    start = None
    for ang in angles:
        plan, g, c, f, obj = optimise(prob, ang, rng, start=start, iters=40000)
        start = plan
        plans[ang] = (plan, g, c, f, obj)
        frontier.append({"Penalty_Degree": ang, "Expected_Gain": g, "Group_Coancestry": c, "Progeny_Inbreeding": f,
                         "N_Parents": int(len(np.unique(plan))), "Objective": obj})

    # the chosen angle gets a longer search from several starts
    best = plans[angle]
    for r in range(3):
        cand = optimise(prob, angle, rng, start=(plans[angle][0] if r == 0 else None), iters=150000)
        if cand[4] > best[4]:
            best = cand
    plan, g, c, f, obj = best
    for row in frontier:
        if row["Penalty_Degree"] == angle:
            row.update({"Expected_Gain": g, "Group_Coancestry": c, "Progeny_Inbreeding": f,
                        "N_Parents": int(len(np.unique(plan))), "Objective": obj})

    # reference: truncation selection (the best parents, mated at random)
    tr_plan = _initial_plan(prob, rng)
    tg, tc, tf = prob.score(tr_plan)

    cnt = np.bincount(plan.ravel(), minlength=n)
    info = pairinfo
    crosses = pd.DataFrame({
        "Cross_ID": np.arange(1, nc + 1),
        "Parent1": [names[i] for i in plan[:, 0]],
        "Parent2": [names[i] for i in plan[:, 1]],
        "P1_Genomic_Value": a[plan[:, 0]],
        "P2_Genomic_Value": a[plan[:, 1]],
        "Cross_Mean": [info["Mean"][i, j] for i, j in plan],
        "Progeny_SD": [info["sd"][i, j] for i, j in plan],
        "Usefulness": U[plan[:, 0], plan[:, 1]],
        "Kinship_Coefficient": C[plan[:, 0], plan[:, 1]],
    }).sort_values("Usefulness", ascending=False)
    crosses["Cross_ID"] = np.arange(1, nc + 1)
    used = np.nonzero(cnt)[0]
    contrib = pd.DataFrame({"Genotype": [names[i] for i in used], crit_label + "_Value": a[used],
                            "N_Crosses": cnt[used], "Contribution": cnt[used] / (2.0 * nc)})
    contrib = contrib.sort_values(["N_Crosses", crit_label + "_Value"], ascending=False)
    front = pd.DataFrame(frontier)

    crosses.round(6).to_csv(results_dir / "ocs_crosses.csv", index=False)
    contrib.round(6).to_csv(results_dir / "ocs_contributions.csv", index=False)
    front.round(6).to_csv(results_dir / "ocs_frontier.csv", index=False)
    plot = _plot_frontier(front, angle, (tg, tc + tf), results_dir / "ocs_frontier.png",
                          ylabel=f"Expected gain ({gain_label})")

    def rows(df):
        return df.round(6).astype(object).where(df.notna(), None).to_dict("records")

    return {
        "summary": {
            "n_candidates": n, "n_crosses": nc, "max_crosses_per_parent": m, "penalty_degree": angle,
            "n_parents": int(len(used)),
            "expected_gain": round(g, 4),
            "mean_usefulness": round(prob.g0, 4), "gain_label": gain_label,
            "group_coancestry": round(c, 5), "progeny_inbreeding": round(f, 5),
            "truncation_gain": round(tg, 4), "truncation_coancestry": round(tc, 5),
            "truncation_inbreeding": round(tf, 5),
        },
        "crosses": rows(crosses), "contributions": rows(contrib), "frontier": rows(front),
        "crosses_csv": _p(results_dir / "ocs_crosses.csv"),
        "contributions_csv": _p(results_dir / "ocs_contributions.csv"),
        "frontier_csv": _p(results_dir / "ocs_frontier.csv"),
        "frontier_plot": plot,
    }


def _p(path):
    return str(Path(path).resolve()).replace("\\", "/")


def _plot_frontier(front, angle, truncation, path, ylabel="Expected gain (mean usefulness of the crosses)"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=110)
    cc = front["Group_Coancestry"] + front["Progeny_Inbreeding"]
    ax.plot(cc, front["Expected_Gain"], "-o", color="#2d6a4f", ms=4, lw=1.4, label="OCS at each penalty degree")
    for _, r in front.iterrows():
        if r["Penalty_Degree"] in (0, 30, 60, 90):
            ax.annotate(f"{int(r['Penalty_Degree'])}°", (r["Group_Coancestry"] + r["Progeny_Inbreeding"], r["Expected_Gain"]),
                        textcoords="offset points", xytext=(5, -10), fontsize=8, color="#40594c")
    sel = front[front["Penalty_Degree"] == angle].iloc[0]
    ax.plot(sel["Group_Coancestry"] + sel["Progeny_Inbreeding"], sel["Expected_Gain"], "o", ms=11, mfc="none",
            mec="#e65a14", mew=2, label=f"chosen ({angle:g}°)")
    ax.plot(truncation[1], truncation[0], "s", color="#990000", ms=7, label="truncation selection")
    ax.set_xlabel("Coancestry (group coancestry + progeny inbreeding)  →  less diversity")
    ax.set_ylabel(ylabel)
    ax.set_title("Usefulness versus diversity")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, facecolor="white")
    plt.close(fig)
    return _p(path)
