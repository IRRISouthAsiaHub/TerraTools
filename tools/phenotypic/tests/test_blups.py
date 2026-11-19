"""BLUPs module of the Phenotypic Data Analysis tool (needs R with lme4, Matrix, emmeans ...; skipped without it)."""
import csv
import random
import tempfile
import unittest
from pathlib import Path

from tools.common.r_locator import find_rscript
from tools.phenotypic.backend.phenotypic_runner import PhenotypicEngine
from tools.phenotypic.config import PhenotypicConfig


def _rcbd_file(path: Path, n_geno=40, reps=3, envs=("A", "B"), seed=3):
    rnd = random.Random(seed)
    g_eff = {f"G{i:03d}": rnd.gauss(0, 2) for i in range(n_geno)}
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Genotype", "Environment", "Rep", "Yield"])
        for e in envs:
            for r in range(1, reps + 1):
                for g, eff in g_eff.items():
                    w.writerow([g, e, r, round(20 + eff + (r - 2) + rnd.gauss(0, 1.5), 3)])


@unittest.skipIf(find_rscript() is None, "Rscript not available")
class BlupTests(unittest.TestCase):
    def test_balanced_rcbd(self):
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "trial.csv"
            _rcbd_file(f)
            cfg = PhenotypicConfig.from_dict({"input_file": str(f), "design": "rcbd", "module": "blups",
                                              "settings": {"traits": ["Yield"], "by_env": False}})
            r = PhenotypicEngine(cfg).run()
            self.assertEqual(r["status"], "success", r.get("message"))
            self.assertFalse([m for m in r.get("msgs", []) if m["type"] == "error"])
            summ = r["tables_data"]["blup_summary"]
            self.assertEqual([s["Environment"] for s in summ], ["Pooled"])
            # one BLUP per accession, ranked best first
            rows = r["tables_data"]["blups"]
            self.assertEqual(len(rows), 40)
            self.assertEqual([int(x["Rank"]) for x in rows], sorted(int(x["Rank"]) for x in rows))
            blups = [float(x["BLUP"]) for x in rows]
            self.assertEqual(blups, sorted(blups, reverse=True))
            # a balanced design: the Cullis H2 (from the PEV matrix) equals the classical entry-mean H2
            h_cullis = float(summ[0]["H² Cullis (BLUP)"])
            h_class = float(summ[0]["H² entry-mean (classical)"])
            self.assertAlmostEqual(h_cullis, h_class, delta=0.002)
            self.assertTrue(0.3 < h_cullis < 1)
            # the BLUPs workbook is offered
            e = PhenotypicEngine().export(r["outdir"], "blups")
            self.assertEqual(e["status"], "success", e.get("message"))
            self.assertTrue(Path(e["path"]).is_file())


if __name__ == "__main__":
    unittest.main()
