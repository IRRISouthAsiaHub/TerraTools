"""
BAM input: validation, allele counting (samtools or the built-in reader) and the marker catalogue.

Every selector below works on arrays of read counts (4 alleles A/C/G/T x sites) so that the same code serves both
engines. "Index allele" = the allele the SNP-index counts; "ref allele" = the other allele of the pair (named Ref in
the tables, it is not necessarily the reference-genome base).

Marker catalogue
    two parents    both parents homozygous for different alleles (confident calls, Phred >= min_parent_gq,
                   depth >= min_parent_depth). The index allele is the allele of the chosen index parent.
    one parent     the parent is homozygous (confident call); the other allele of the pair is the commonest allele
                   in the bulks that the parent does not carry. Parental polymorphism cannot be fully validated.
    no parent      the reference-genome base against the commonest other allele in the bulks; the index allele is
                   the non-reference one. No parental validation at all.
"""

from __future__ import annotations

import os
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from tools.qtl_seq.backend import deps
from tools.qtl_seq.backend.bamreader import BamCounter, BamError, read_header
from tools.qtl_seq.backend.calling import homozygous_call

TILE = 2_000_000
BATCH_LINES = 200_000
BASES = "ACGT"
_LUT = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _LUT[_b] = _i
    _LUT[_b + 32] = _i


ENGINE_LABEL = {"samtools": "samtools mpileup", "pysam": "pysam", "builtin": "built-in BAM reader"}


class BamPipelineError(Exception):
    pass


# ============================================================
# VALIDATION
# ============================================================

def has_index(path):
    p = str(path)
    return any(os.path.exists(x) for x in (p + ".bai", p + ".csi", re.sub(r"\.bam$", ".bai", p)))


def validate_bam_inputs(cfg, fasta=None):
    """Checks every BAM: is it a BAM, indexed, coordinate-sorted, on the reference of the FASTA and on the same
    reference as the other BAMs. Returns {'errors', 'warnings', 'files'}."""
    errors, warnings, files = [], [], {}
    roles = [("high_bam", "high bulk"), ("low_bam", "low bulk"), ("parent1_bam", "parent 1"), ("parent2_bam", "parent 2")]
    first = None
    for key, label in roles:
        path = cfg.get(key) or ""
        if not path:
            if key in ("high_bam", "low_bam"):
                errors.append(f"The {label} BAM is not selected.")
            continue
        if not os.path.isfile(path):
            errors.append(f"The {label} BAM was not found: {path}")
            continue
        try:
            refs, text = read_header(path)
        except BamError as e:
            errors.append(f"{label}: {e}")
            continue
        info = {"file": Path(path).name, "contigs": len(refs), "indexed": has_index(path)}
        so = re.search(r"SO:(\w+)", text)
        info["sort_order"] = so.group(1) if so else "unknown"
        files[label] = info
        if not info["indexed"]:
            warnings.append(f"{label} BAM has no index (.bai); the built-in reader does not need one, samtools will create it.")
        if so and so.group(1) != "coordinate":
            errors.append(f"{label} BAM is not coordinate-sorted (SO:{so.group(1)}). Run samtools sort.")
        names = {n: l for n, l in refs}
        if first is None:
            first = (label, names)
        elif names != first[1]:
            diff = sorted(set(names) ^ set(first[1]))[:5]
            errors.append(f"The {label} BAM and the {first[0]} BAM were not aligned to the same reference "
                          f"(contigs differ{': ' + ', '.join(diff) if diff else ' in length'}).")
        if fasta is not None:
            fa = fasta.contigs()
            common = set(fa) & set(names)
            if not common:
                errors.append(f"The {label} BAM and the reference FASTA have no contig name in common "
                              f"(BAM: {', '.join(list(names)[:3])} ... FASTA: {', '.join(list(fa)[:3])} ...).")
            else:
                bad = [n for n in common if fa[n] != names[n]]
                if bad:
                    errors.append(f"{label} BAM: contig length differs from the FASTA for {', '.join(bad[:5])} "
                                  "(a different assembly version?).")
    if not cfg.get("ref_fasta"):
        errors.append("The reference FASTA is not selected.")
    elif not os.path.isfile(cfg["ref_fasta"]):
        errors.append(f"The reference FASTA was not found: {cfg['ref_fasta']}")
    return {"errors": errors, "warnings": warnings, "files": files}


# ============================================================
# COUNTING ENGINES
# ============================================================

_CARET = re.compile(r"\^.")


def _strip_indels(s):
    out, i, n = [], 0, len(s)
    while i < n:
        c = s[i]
        if c in "+-":
            j = i + 1
            while j < n and s[j].isdigit():
                j += 1
            i = j + int(s[i + 1:j] or 0)
        else:
            out.append(c)
            i += 1
    return "".join(out)


def count_pileup_string(reads, ref_idx):
    """[A, C, G, T] read counts from the read-base column of one sample in `samtools mpileup` output."""
    if "^" in reads:
        reads = _CARET.sub("", reads)
    if "+" in reads or "-" in reads:
        reads = _strip_indels(reads)
    c = [reads.count(b) + reads.count(b.lower()) for b in BASES]
    if 0 <= ref_idx < 4:
        c[ref_idx] += reads.count(".") + reads.count(",")
    return c


def parse_mpileup_line(line, n_samples):
    """(chrom, pos, ref index, counts (n_samples, 4)) or None."""
    f = line.rstrip("\n").split("\t")
    if len(f) < 3 + 3 * n_samples:
        return None
    ref = _LUT[ord(f[2][0])] if f[2] else 4
    if ref > 3:
        return None
    return f[0], int(f[1]), int(ref), [count_pileup_string(f[4 + 3 * i], int(ref)) for i in range(n_samples)]


class Engine:
    """Iterates the sites of a set of BAMs and hands arrays to a handler:
    handler(chrom, positions (n,), reference base codes (n,), {role: counts (4, n)})."""

    def __init__(self, cfg, fasta, chroms, log, cancelled):
        self.cfg, self.fasta, self.chroms, self.log, self.cancelled = cfg, fasta, chroms, log, cancelled
        name = cfg.get("engine", "auto")
        if name not in ("samtools", "pysam", "builtin"):
            name = "samtools" if deps.where("samtools") else ("pysam" if deps.pysam_where() else "builtin")
        self.name = name
        if name == "samtools" and not deps.where("samtools"):
            raise BamPipelineError("samtools was not found on PATH or inside WSL. Choose another reader or install samtools.")
        if name == "pysam" and not deps.pysam_where():
            raise BamPipelineError("pysam was not found (neither in this Python nor in the Python of WSL). "
                                   "Choose another reader or install pysam.")

    def check(self):
        if self.cancelled and self.cancelled():
            raise InterruptedError("cancelled")

    # ---------------------------------------------------------- samtools
    def _scan_samtools(self, roles, handler, positions):
        cfg, n = self.cfg, len(roles)
        for ci, chrom in enumerate(self.chroms):
            self.check()
            pos_file = None
            if positions is not None:
                cat = positions.get(chrom)
                if cat is None or len(cat) == 0:
                    continue
                pos_file = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
                pos_file.write("".join(f"{chrom}\t{int(p)}\n" for p in cat))
                pos_file.close()
            args = ["mpileup", "-f", cfg["ref_fasta"], "-r", chrom, "-q", cfg["min_mapq"], "-Q", cfg["min_baseq"], "-B",
                    "-d", 100000]
            if pos_file is not None:
                args += ["-l", pos_file.name]
            args += [p for _r, p in roles]
            cmd = deps.command("samtools", *args)
            self.log(f"samtools mpileup {chrom} ({ci + 1}/{len(self.chroms)})")
            err = tempfile.TemporaryFile()
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, text=True, bufsize=1 << 20, **deps.run_kwargs())
            try:
                pos, ref, cnt = [], [], []

                def flush():
                    if not pos:
                        return
                    arr = np.asarray(cnt, dtype=np.int32)                 # (sites, samples, 4)
                    handler(chrom, np.asarray(pos, dtype=np.int64), np.asarray(ref, dtype=np.uint8),
                            {role: arr[:, k, :].T.copy() for k, (role, _p) in enumerate(roles)})
                    pos.clear(), ref.clear(), cnt.clear()

                for line in proc.stdout:
                    site = parse_mpileup_line(line, n)
                    if site is None:
                        continue
                    pos.append(site[1]), ref.append(site[2]), cnt.append(site[3])
                    if len(pos) >= BATCH_LINES:
                        flush()
                        self.check()
                flush()
                proc.wait()
                if proc.returncode != 0:
                    err.seek(0)
                    raise BamPipelineError("samtools mpileup failed: " + err.read().decode(errors="replace")[-1500:])
            finally:
                if proc.poll() is None:
                    proc.kill()
                err.close()
                if pos_file is not None:
                    os.unlink(pos_file.name)

    # ---------------------------------------------------------- pysam
    def _scan_pysam(self, roles, handler, positions):
        """pysam counts per site: in this process where pysam can be imported, otherwise by pysam_counts.py under WSL."""
        cfg, n = self.cfg, len(roles)
        paths = [p for _r, p in roles]
        native = deps.pysam_where() == "python"
        script = Path(__file__).with_name("pysam_counts.py")
        if not native and not script.exists():
            raise BamPipelineError("pysam_counts.py is missing from this installation, so pysam cannot be used through WSL.")

        def emit(chrom, pos, cnt):
            pos = np.asarray(pos, dtype=np.int64)
            cnt = np.asarray(cnt, dtype=np.int32).reshape(len(pos), n, 4)
            seq = self.fasta.fetch(chrom, int(pos[0]), int(pos[-1]))
            codes = _LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)] if len(seq) == pos[-1] - pos[0] + 1 \
                else np.full(pos[-1] - pos[0] + 1, 4, np.uint8)
            handler(chrom, pos, codes[pos - pos[0]], {role: cnt[:, k, :].T.copy() for k, (role, _p) in enumerate(roles)})

        for ci, chrom in enumerate(self.chroms):
            self.check()
            cat = None
            if positions is not None:
                cat = positions.get(chrom)
                if cat is None or len(cat) == 0:
                    continue
            self.log(f"pysam {chrom} ({ci + 1}/{len(self.chroms)})")
            if native:
                from tools.qtl_seq.backend.pysam_counts import iter_counts
                for pos, rows in iter_counts(paths, chrom, cfg["min_mapq"], cfg["min_baseq"],
                                             None if cat is None else [int(x) for x in cat], should_stop=self.cancelled):
                    emit(chrom, pos, rows)
                    self.check()
                continue
            pos_file = None
            args = ["--chrom", chrom, "--mapq", str(cfg["min_mapq"]), "--baseq", str(cfg["min_baseq"])]
            if cat is not None:
                pos_file = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
                pos_file.write("".join(f"{chrom}\t{int(p)}\n" for p in cat))
                pos_file.close()
                args += ["--positions", pos_file.name]
            err = tempfile.TemporaryFile()
            proc = subprocess.Popen(deps.pysam_script_command(script, *args, *paths), stdout=subprocess.PIPE, stderr=err, text=True,
                                    bufsize=1 << 20, **deps.run_kwargs())
            try:
                lines = []

                def flush():
                    if not lines:
                        return
                    a = np.array("\t".join(lines).split(), dtype=np.int64).reshape(len(lines), 1 + 4 * n)
                    emit(chrom, a[:, 0], a[:, 1:])
                    lines.clear()

                for line in proc.stdout:
                    lines.append(line.rstrip("\n"))
                    if len(lines) >= BATCH_LINES:
                        flush()
                        self.check()
                flush()
                proc.wait()
                if proc.returncode != 0:
                    err.seek(0)
                    raise BamPipelineError("pysam failed: " + err.read().decode(errors="replace")[-1500:])
            finally:
                if proc.poll() is None:
                    proc.kill()
                err.close()
                if pos_file is not None:
                    os.unlink(pos_file.name)

    # ---------------------------------------------------------- built-in
    def _scan_builtin(self, roles, handler, positions):
        cfg = self.cfg
        counters = [(role, BamCounter(path, cfg["min_mapq"], cfg["min_baseq"])) for role, path in roles]
        try:
            names = [n for n, _l in counters[0][1].refs]
            lengths = dict(counters[0][1].refs)
            for chrom in self.chroms:
                ref_id = names.index(chrom)
                cat = None if positions is None else positions.get(chrom)
                if positions is not None and (cat is None or len(cat) == 0):
                    continue
                self.log(f"Reading {chrom} ({self.chroms.index(chrom) + 1}/{len(self.chroms)})")
                for start in range(0, lengths[chrom], TILE):
                    self.check()
                    end = min(start + TILE, lengths[chrom])
                    arrays = {role: c.tile(ref_id, start, end) for role, c in counters}
                    pos = np.arange(start + 1, end + 1, dtype=np.int64)
                    if cat is not None:
                        sel = np.flatnonzero((cat >= start + 1) & (cat <= end))
                        if sel.size == 0:
                            continue
                        idx = cat[sel] - start - 1
                        pos = pos[idx]
                        arrays = {r: a[:, idx] for r, a in arrays.items()}
                    else:
                        idx = None
                    seq = self.fasta.fetch(chrom, start + 1, end)
                    ref = _LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)] if len(seq) == end - start else np.full(end - start, 4, np.uint8)
                    handler(chrom, pos, ref if idx is None else ref[idx], {r: a.astype(np.int32) for r, a in arrays.items()})
        finally:
            for _r, c in counters:
                c.close()

    def scan(self, roles, handler, positions=None):
        """roles: [(role name, bam path)]. positions: {chrom: sorted array of 1-based positions} to restrict to."""
        if self.name == "samtools":
            for _r, path in roles:
                if not has_index(path):
                    self.log(f"Indexing {Path(path).name}")
                    subprocess.run(deps.command("samtools", "index", path), check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   **deps.run_kwargs())
            self._scan_samtools(roles, handler, positions)
        elif self.name == "pysam":
            self._scan_pysam(roles, handler, positions)
        else:
            self._scan_builtin(roles, handler, positions)


# ============================================================
# SELECTORS
# ============================================================

def _top_depth(a):
    return a.argmax(axis=0), a.sum(axis=0)


class Rows:
    """Accumulates the selected sites: chrom, pos, ref allele index, index allele index and the 4 bulk counts."""

    def __init__(self):
        self.parts = []

    def add(self, chrom, pos, ref_idx, alt_idx, H, L):
        k = np.arange(len(pos))
        self.parts.append(pd.DataFrame({
            "Chromosome": chrom, "Position": pos, "Ref_i": ref_idx, "Alt_i": alt_idx,
            "Ref_High": H[ref_idx, k], "Alt_High": H[alt_idx, k], "Ref_Low": L[ref_idx, k], "Alt_Low": L[alt_idx, k]}))

    def table(self):
        if not self.parts:
            return pd.DataFrame(columns=["Chromosome", "Position", "Ref", "Alt", "Ref_High", "Alt_High", "Ref_Low", "Alt_Low"])
        t = pd.concat(self.parts, ignore_index=True)
        t.insert(2, "Ref", np.array(list(BASES))[t.pop("Ref_i").to_numpy()])
        t.insert(3, "Alt", np.array(list(BASES))[t.pop("Alt_i").to_numpy()])
        return t


class TwoParents:
    """Stage 1 of the two-parent mode: sites where the parents are confidently homozygous for different alleles."""

    def __init__(self, cfg):
        self.cfg = cfg
        self.identified = self.passed = 0
        self.cat = {}                                   # chrom -> [positions, parent1 allele, parent2 allele]

    def __call__(self, chrom, pos, ref, arr):
        c = self.cfg
        t1, d1 = _top_depth(arr["P1"])
        t2, d2 = _top_depth(arr["P2"])
        cand = np.flatnonzero((d1 >= c["min_parent_depth"]) & (d2 >= c["min_parent_depth"]) & (t1 != t2))
        if cand.size == 0:
            return
        self.identified += int(cand.size)
        _t, _d, g1 = homozygous_call(arr["P1"][:, cand], c["base_error"])
        _t, _d, g2 = homozygous_call(arr["P2"][:, cand], c["base_error"])
        ok = (g1 >= c["min_parent_gq"]) & (g2 >= c["min_parent_gq"])
        sel = cand[ok]
        self.passed += int(sel.size)
        if sel.size:
            self.cat.setdefault(chrom, []).append((pos[sel], t1[sel], t2[sel]))

    def catalogue(self):
        out = {}
        for chrom, parts in self.cat.items():
            p = np.concatenate([x[0] for x in parts])
            o = np.argsort(p)
            out[chrom] = (p[o], np.concatenate([x[1] for x in parts])[o], np.concatenate([x[2] for x in parts])[o])
        return out


class BulksAtCatalogue:
    """Stage 2 of the two-parent mode: bulk counts at the catalogue positions."""

    def __init__(self, cat, index_parent, rows):
        self.cat, self.index_parent, self.rows = cat, index_parent, rows

    def __call__(self, chrom, pos, ref, arr):
        cp, a1, a2 = self.cat[chrom]
        i = np.clip(np.searchsorted(cp, pos), 0, len(cp) - 1)
        ok = cp[i] == pos
        if not ok.any():
            return
        i, pos = i[ok], pos[ok]
        alt = (a2 if self.index_parent == "parent2" else a1)[i]
        refa = (a1 if self.index_parent == "parent2" else a2)[i]
        self.rows.add(chrom, pos, refa, alt, arr["H"][:, ok], arr["L"][:, ok])


class OneParent:
    """The parent is homozygous; the other allele is the commonest allele in the bulks that the parent lacks."""

    def __init__(self, cfg, which, rows):
        self.cfg, self.which, self.rows = cfg, which, rows
        self.identified = self.passed = 0

    def __call__(self, chrom, pos, ref, arr):
        c = self.cfg
        P = arr["P"]
        top, d = _top_depth(P)
        pooled = arr["H"] + arr["L"]
        k = np.arange(pos.size)
        other = pooled.copy()
        other[top, k] = 0
        y = other.argmax(axis=0)
        ny = other[y, k]
        share = ny / np.maximum(pooled.sum(axis=0), 1)
        cand = np.flatnonzero((d >= c["min_parent_depth"]) & (ny >= c["min_alt_reads"]) & (share >= c["min_alt_fraction"]))
        if cand.size == 0:
            return
        self.identified += int(cand.size)
        _t, _d, gq = homozygous_call(P[:, cand], c["base_error"])
        sel = cand[gq >= c["min_parent_gq"]]
        self.passed += int(sel.size)
        if sel.size == 0:
            return
        x_par, y_oth = top[sel], y[sel]
        index_is_parent_allele = self.which == c["index_parent"]
        alt = x_par if index_is_parent_allele else y_oth
        refa = y_oth if index_is_parent_allele else x_par
        self.rows.add(chrom, pos[sel], refa, alt, arr["H"][:, sel], arr["L"][:, sel])


class NoParent:
    """Reference base against the commonest other allele of the pooled bulks; the non-reference allele is the index allele."""

    def __init__(self, cfg, rows):
        self.cfg, self.rows = cfg, rows
        self.identified = self.passed = 0

    def __call__(self, chrom, pos, ref, arr):
        c = self.cfg
        pooled = arr["H"] + arr["L"]
        k = np.arange(pos.size)
        valid = ref < 4
        other = pooled.copy()
        other[np.minimum(ref, 3), k] = 0
        y = other.argmax(axis=0)
        ny = other[y, k]
        share = ny / np.maximum(pooled.sum(axis=0), 1)
        sel = np.flatnonzero(valid & (ny >= c["min_alt_reads"]) & (share >= c["min_alt_fraction"]))
        if sel.size == 0:
            return
        self.identified += int(sel.size)
        self.passed += int(sel.size)
        self.rows.add(chrom, pos[sel], ref[sel].astype(np.int64), y[sel], arr["H"][:, sel], arr["L"][:, sel])


# ============================================================
# PIPELINE
# ============================================================

def choose_chromosomes(cfg, fasta, bam_contigs, warnings):
    fa = fasta.contigs()
    common = [c for c in bam_contigs if c in fa]
    wanted = [c.strip() for c in str(cfg.get("chromosomes", "")).replace(";", ",").split(",") if c.strip()]
    if wanted:
        missing = [c for c in wanted if c not in common]
        if missing:
            raise BamPipelineError("These chromosomes are not in both the BAM files and the FASTA: " + ", ".join(missing))
        return wanted
    if len(common) > 100:
        keep = [c for c in common if fa[c] >= 1_000_000]
        warnings.append(f"The reference has {len(common)} contigs; only the {len(keep)} of at least 1 Mb are analysed. "
                        "List the chromosomes explicitly in the advanced settings to change this.")
        return keep
    return common


def build_bam_snps(cfg, fasta, log=print, cancelled=None):
    """BAM mode. Returns (snps, ledger, chrom_lengths, warnings, input info)."""
    p = cfg
    check = validate_bam_inputs(p, fasta)
    if check["errors"]:
        raise BamPipelineError("\n".join(check["errors"]))
    warnings = list(check["warnings"])
    refs, _t = read_header(p["high_bam"])
    chroms = choose_chromosomes(p, fasta, [n for n, _l in refs], warnings)
    lengths = {c: fasta.contigs()[c] for c in chroms}
    eng = Engine(p, fasta, chroms, log, cancelled)
    log(f"BAM engine: {ENGINE_LABEL[eng.name]}")
    if eng.name == "pysam":
        warnings.append("Reads were counted with pysam, which counts both mates of an overlapping read pair (samtools mpileup counts "
                        "the overlap once), so depths can be slightly higher at the overlap of paired-end reads.")
    if eng.name == "builtin":
        warnings.append("samtools was not used: the built-in BAM reader is slower (it is meant for small genomes and "
                        "tests). Install samtools for whole-genome data.")
    p1, p2 = p.get("parent1_bam") or "", p.get("parent2_bam") or ""
    ledger = {}
    rows = Rows()
    if p1 and p2:
        mode = "two_parents"
        log("Step 1: parental genotypes")
        sel = TwoParents(p)
        eng.scan([("P1", p1), ("P2", p2)], sel)
        cat = sel.catalogue()
        log(f"Step 1 done: {sel.identified:,} sites differ between the parents, {sel.passed:,} are confident homozygous calls")
        ledger.update(total_snps=sel.identified, removed_parental=sel.identified - sel.passed, passing_parental=sel.passed)
        if not sel.passed:
            raise BamPipelineError("No SNP is homozygous for different alleles in the two parents with the current depth/quality "
                                   "settings. Check that the parents are the parents of this population, or lower the "
                                   "parental depth / genotype quality.")
        log("Step 2: allele counts of the bulks at the parental SNPs")
        eng.scan([("H", p["high_bam"]), ("L", p["low_bam"])], BulksAtCatalogue(cat, p["index_parent"], rows),
                 positions={c: v[0] for c, v in cat.items()})
    elif p1 or p2:
        mode = "one_parent"
        which = "parent1" if p1 else "parent2"
        warnings.append(f"Only {which} was supplied: the other parent is assumed to carry the other allele, so parental "
                        "polymorphism cannot be completely validated.")
        sel = OneParent(p, which, rows)
        log("Scanning the parent and both bulks")
        eng.scan([("P", p1 or p2), ("H", p["high_bam"]), ("L", p["low_bam"])], sel)
        ledger.update(total_snps=sel.identified, removed_parental=sel.identified - sel.passed, passing_parental=sel.passed)
    else:
        mode = "no_parent"
        warnings.append("No parental BAM: sites are the reference base against the commonest other allele of the pooled bulks, "
                        "and the index allele is the non-reference allele. Parental polymorphism cannot be checked, so sites "
                        "that are not polymorphic between the parents and sequencing errors can remain, and the null model "
                        "assumes the non-reference allele has frequency 1/2 in the population.")
        sel = NoParent(p, rows)
        log("Scanning both bulks")
        eng.scan([("H", p["high_bam"]), ("L", p["low_bam"])], sel)
        ledger.update(total_snps=sel.identified, removed_parental=0, passing_parental=sel.passed)
    ledger["mode"] = mode
    ledger["removed_quality"] = 0                       # base/mapping quality is applied while counting reads
    table = rows.table()
    if table.empty:
        raise BamPipelineError("No informative SNP was found in the BAM files.")
    return table, ledger, lengths, warnings, {"engine": eng.name, "files": check["files"], "chromosomes": chroms}
