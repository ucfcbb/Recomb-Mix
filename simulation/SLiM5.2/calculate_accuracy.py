# This program was generated with the assistance of Claude
import sys
import argparse
import bisect
import gzip


def resolve_method(name):
    norm = "".join(ch for ch in name.lower() if ch not in "-_ ")
    return norm if norm in METHODS else None


def fields(line):
    return [f.strip() for f in line.rstrip("\r\n").split("\t") if f.strip() != ""]


def die(msg):
    sys.exit(f"ERROR: {msg}")


# ----------------------------------------------------------------------------
# Readers
# ----------------------------------------------------------------------------


def read_sites(path, chrom=None):
    sites, order, chroms, seen_data = set(), [], set(), False
    with open(path) as fh:
        for n, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            f = fields(line)
            if len(f) < 2:
                die(f"{path} line {n}: expected at least 2 tab-separated columns")
            try:
                pos = int(f[1])
            except ValueError:
                if not seen_data:          # header line written without '#'
                    continue
                die(f"{path} line {n}: position '{f[1]}' is not an integer")
            seen_data = True
            chroms.add(f[0])
            if chrom is None or f[0] == chrom:
                sites.add(pos)
                order.append(pos)
    if chrom is None and len(chroms) > 1:
        die(f"{path} contains {len(chroms)} chromosomes; select one with --chrom (the ground truth and prediction files have no chromosome column)")
    if not sites:
        die(f"no positions found in {path}" + (f" for chromosome '{chrom}'" if chrom else ""))
    return sorted(sites), order


def parse_segments(path, n, hap, vals, id2name):
    if len(vals) % 3 != 0:
        die(f"{path} line {n}: haplotype '{hap}' has {len(vals)} values after the ID; "
            "expected a multiple of 3 (start, end, label)")
    segs = []
    for i in range(0, len(vals), 3):
        try:
            start, end = int(vals[i]), int(vals[i + 1])
        except ValueError:
            die(f"{path} line {n}: non-integer position in segment {i // 3 + 1} of '{hap}'")
        label = vals[i + 2]
        if label not in id2name:
            die(f"{path} line {n}: label id '{label}' of '{hap}' is not defined in the "
                "population/ancestry list")
        segs.append((start, end, id2name[label]))
    return segs


def read_truth(path):
    id2name, haps = {}, {}
    with open(path) as fh:
        for n, line in enumerate(fh, start=1):
            if not line.strip():
                continue
            if line.startswith("#"):
                if line.startswith("#Population_label_and_ID"):
                    for tok in fields(line)[1:]:
                        name, eq, idx = tok.partition("=")
                        if not eq:
                            die(f"{path} line {n}: expected NAME=ID, got '{tok}'")
                        id2name[idx.strip()] = name.strip()
                continue
            if not id2name:
                die(f"{path}: the '#Population_label_and_ID:' header line was not found before the first data line")
            f = fields(line)
            if f[0] in haps:
                die(f"{path} line {n}: duplicate haplotype id '{f[0]}'")
            haps[f[0]] = parse_segments(path, n, f[0], f[1:], id2name)
    if not haps:
        die(f"no haplotypes found in {path}")
    return haps, id2name


def read_recombmix(path):
    id2name, positions, haps = {}, [], {}
    section = None
    with open_text(path) as fh:
        for n, line in enumerate(fh, start=1):
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                key = s[1:].strip().lower()
                if key.startswith("ancestry"):
                    section = "ancestry"
                elif key.startswith("position"):
                    section = "position"
                elif key.startswith("result"):
                    section = "result"
                else:
                    die(f"{path} line {n}: unrecognised section header '{s}'")
                continue
            f = fields(line)
            if section == "ancestry":
                if len(f) < 2:
                    die(f"{path} line {n}: expected 'population<TAB>id'")
                id2name[f[1]] = f[0]
            elif section == "position":
                try:
                    positions.append(int(f[0]))
                except ValueError:
                    die(f"{path} line {n}: position '{f[0]}' is not an integer")
            elif section == "result":
                if f[0] in haps:
                    die(f"{path} line {n}: duplicate haplotype id '{f[0]}'")
                haps[f[0]] = parse_segments(path, n, f[0], f[1:], id2name)
            else:
                die(f"{path} line {n}: data found before any '#Ancestry/#Position/#Result' header")
    if not id2name:
        die(f"{path}: no #Ancestry section found")
    if not haps:
        die(f"{path}: no #Result section / haplotypes found")
    return haps, positions


def open_text(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def norm_chrom(c):
    c = c.strip()
    return c[3:] if c.lower().startswith("chr") else c


def read_flare(path, chrom=None):
    id2name, samples = {}, None
    cur = start = last = None
    haps, positions, chroms = {}, [], set()
    idx_cache = {}

    def close(h):
        label = cur[h]
        if label is None or label == ".":
            return
        if label not in id2name:
            die(f"{path}: ancestry id '{label}' is not defined in the ##ANCESTRY header line")
        haps[hap_ids[h]].append((start[h], last[h], id2name[label]))

    with open_text(path) as fh:
        for n, line in enumerate(fh, start=1):
            if line.startswith("##"):
                if line.startswith("##ANCESTRY="):
                    body = line[len("##ANCESTRY="):].strip().strip("<>")
                    for tok in body.split(","):
                        name, eq, idx = tok.partition("=")
                        if not eq:
                            die(f"{path} line {n}: cannot parse '##ANCESTRY' entry '{tok}'")
                        id2name[idx.strip()] = name.strip()
                continue
            if line.startswith("#"):                 # #CHROM header
                cols = line.rstrip("\r\n").split("\t")
                samples = cols[9:]
                if not id2name:
                    die(f"{path}: no '##ANCESTRY=<...>' header line found (is this a FLARE output VCF?)")
                if not samples:
                    die(f"{path}: no sample columns found")
                hap_ids = [f"{sm}_{k}" for sm in samples for k in (0, 1)]
                haps = {h: [] for h in hap_ids}
                cur = [None] * len(hap_ids)
                start = [0] * len(hap_ids)
                last = [0] * len(hap_ids)
                continue
            if samples is None:
                die(f"{path} line {n}: data line before the '#CHROM' header line")
            f = line.rstrip("\r\n").split("\t")
            chroms.add(norm_chrom(f[0]))
            if chrom is not None and norm_chrom(f[0]) != norm_chrom(chrom):
                continue
            try:
                pos = int(f[1])
            except ValueError:
                die(f"{path} line {n}: POS '{f[1]}' is not an integer")
            fmt = f[8]
            if fmt not in idx_cache:
                keys = fmt.split(":")
                if "AN1" not in keys or "AN2" not in keys:
                    die(f"{path} line {n}: FORMAT '{fmt}' has no AN1/AN2 fields")
                idx_cache[fmt] = (keys.index("AN1"), keys.index("AN2"))
            i1, i2 = idx_cache[fmt]
            positions.append(pos)
            for k in range(len(samples)):
                parts = f[9 + k].split(":")
                for h, val in ((2 * k, parts[i1]), (2 * k + 1, parts[i2])):
                    if val != cur[h]:
                        close(h)
                        cur[h], start[h] = val, pos
                    last[h] = pos
    if chrom is None and len(chroms) > 1:
        die(f"{path} contains {len(chroms)} chromosomes; select one with --chrom")
    if not positions:
        die(f"{path}: no variant records" + (f" for chromosome '{chrom}'" if chrom else ""))
    for h in range(len(hap_ids)):
        close(h)
    return haps, positions


def read_msp(path, chrom=None, end_exclusive=False):
    id2name, hap_ids, ncols, chroms, haps = {}, None, 0, set(), {}
    pending = None

    def add_window(start, end, labs, n):
        for h, lab in zip(hap_ids, labs):
            lab = lab.strip()
            if lab not in id2name:
                die(f"{path} line {n}: ancestry id '{lab}' is not in the population list")
            segs = haps[h]
            if segs and segs[-1][2] == id2name[lab] and start > segs[-1][1]:
                segs[-1] = (segs[-1][0], end, segs[-1][2])        # merge equal neighbours
            else:
                segs.append((start, end, id2name[lab]))

    with open_text(path) as fh:
        for n, line in enumerate(fh, start=1):
            if not line.strip():
                continue
            if line.startswith("#"):
                if line.lower().startswith("#chm"):                  # column header
                    cols = line.rstrip("\r\n").split("\t")
                    if len(cols) < 7 or "snps" not in cols[5].lower().replace(" ", ""):
                        die(f"{path} line {n}: expected the six columns '#chm, spos, epos, "
                            f"sgpos, egpos, n snps' before the haplotype columns, got {cols[:7]}")
                    hap_ids = []
                    for c in cols[6:]:
                        sample, dot, hap = c.strip().rpartition(".")
                        if not dot or not sample:
                            die(f"{path} line {n}: haplotype column '{c}' has no '.<0|1>' suffix")
                        hap_ids.append(f"{sample}_{hap}")
                    if len(set(hap_ids)) != len(hap_ids):
                        die(f"{path} line {n}: duplicate haplotype column names")
                    haps = {h: [] for h in hap_ids}
                    ncols = len(cols)
                elif ":" in line and "=" in line and not id2name:   # population list
                    body = line.split(":", 1)[1].strip()
                    toks = [t.strip() for t in body.split("\t") if t.strip()]
                    if len(toks) == 1 and toks[0].count("=") > 1:
                        toks = body.split()
                    for tok in toks:
                        name, eq, idx = tok.rpartition("=")
                        if not eq:
                            die(f"{path} line {n}: expected NAME=ID, got '{tok}'")
                        id2name[idx.strip()] = name.strip()
                continue
            if hap_ids is None or not id2name:
                die(f"{path} line {n}: data line before the '#Subpopulation order/codes:' and '#chm ...' header lines")
            f = line.rstrip("\r\n").split("\t")
            if len(f) != ncols:
                die(f"{path} line {n}: {len(f)} columns, expected {ncols}")
            chroms.add(norm_chrom(f[0]))
            if chrom is not None and norm_chrom(f[0]) != norm_chrom(chrom):
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                die(f"{path} line {n}: spos/epos are not integers")
            if end_exclusive:
                if pending is not None and pending[1] - 1 >= pending[0]:
                    add_window(pending[0], pending[1] - 1, pending[2], pending[3])
                pending = (start, end, f[6:], n)
            else:
                add_window(start, end, f[6:], n)
    if pending is not None:                       # last line: end position is inclusive
        add_window(pending[0], pending[1], pending[2], pending[3])
    if hap_ids is None or not id2name:
        die(f"{path}: header lines not found (is this a Gnomix/RFMix .msp file?)")
    if chrom is None and len(chroms) > 1:
        die(f"{path} contains {len(chroms)} chromosomes; select one with --chrom")
    if not any(haps.values()):
        die(f"{path}: no windows" + (f" for chromosome '{chrom}'" if chrom else ""))
    return haps, None


def read_sparsepainter(path, map_order):
    names, haps, cur, prev_end = None, {}, None, None
    n_map = len(map_order)
    with open_text(path) as fh:
        for n, line in enumerate(fh, start=1):
            tok = line.split()
            if not tok or line.startswith("#"):
                continue
            if names is None:                                    # column header line
                if len(tok) < 3 or not tok[0].lower().startswith("snpidx"):
                    die(f"{path} line {n}: expected the column header line 'SNPidx_start SNPidx_end <population names...>'")
                names = tok[2:]
                continue
            if len(tok) == 1:                                    # haplotype id line
                if tok[0] in haps:
                    die(f"{path} line {n}: duplicate haplotype id '{tok[0]}'")
                cur = haps[tok[0]] = []
                prev_end = None
                continue
            if cur is None:
                die(f"{path} line {n}: range line before any haplotype id line")
            if len(tok) != 2 + len(names):
                die(f"{path} line {n}: {len(tok)} values, expected {2 + len(names)} "
                    "(start, end and one probability per population)")
            try:
                a, b = int(tok[0]), int(tok[1])
                probs = [float(x) for x in tok[2:]]
            except ValueError:
                die(f"{path} line {n}: non-numeric value")
            if not (1 <= a <= b <= n_map):
                die(f"{path} line {n}: index range {a}-{b} is outside 1..{n_map} (the number of map rows used); indices must be 1-based positions in the map file")
            best = names[max(range(len(probs)), key=probs.__getitem__)]
            if cur and cur[-1][2] == best and prev_end is not None and a == prev_end + 1:
                cur[-1] = (cur[-1][0], map_order[b - 1], best)        # merge equal neighbours
            else:
                cur.append((map_order[a - 1], map_order[b - 1], best))
            prev_end = b
    if names is None:
        die(f"{path}: SparsePainter column header line not found")
    if not haps:
        die(f"{path}: no haplotypes found")
    return haps, None


def read_gnomix(path, chrom=None):
    return read_msp(path, chrom, end_exclusive=False)


def read_rfmix(path, chrom=None):
    return read_msp(path, chrom, end_exclusive=True)


def read_prediction(method, path, chrom=None, map_order=None):
    if method == "recombmix":
        return read_recombmix(path)
    if method == "flare":
        return read_flare(path, chrom)
    if method == "gnomix":
        return read_gnomix(path, chrom)
    if method == "rfmix":
        return read_rfmix(path, chrom)
    if method == "sparsepainter":
        return read_sparsepainter(path, map_order)
    die(f"the predicted-label format of '{method}' is not available (currently supported: recombmix, flare, gnomix, rfmix, sparsepainter)")


INCLUSIVE_END = {"recombmix": True, "flare": True, "gnomix": True, "rfmix": True, "sparsepainter": True}

# ----------------------------------------------------------------------------
# Evaluation
# ----------------------------------------------------------------------------


def to_index_intervals(segs, sites, inclusive_end, what):
    out = []
    for start, end, name in sorted(segs, key=lambda x: x[0]):
        lo = bisect.bisect_left(sites, start)
        hi = bisect.bisect_right(sites, end) if inclusive_end else bisect.bisect_left(sites, end)
        if hi > lo:
            out.append((lo, hi, name))
    for a, b in zip(out, out[1:]):
        if b[0] < a[1]:
            die(f"overlapping segments in {what}: a position lies in two ranges (if a range's end equals the next range's start, the end is probably exclusive - try --truth-end exclusive)")
    return out


def count_correct(truth_iv, pred_iv):
    correct, j = {}, 0
    for lo, hi, name in truth_iv:
        while j < len(pred_iv) and pred_iv[j][1] <= lo:
            j += 1
        k = j
        while k < len(pred_iv) and pred_iv[k][0] < hi:
            plo, phi, pname = pred_iv[k]
            if pname == name:
                correct[name] = correct.get(name, 0) + min(hi, phi) - max(lo, plo)
            k += 1
    return correct


if __name__ == "__main__":
    METHODS = ["rfmix", "recombmix", "flare", "gnomix", "sparsepainter"]

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--method", required=True, help="one of: " + ", ".join(METHODS))
    ap.add_argument("--map-in", required=True)
    ap.add_argument("--truth-in", required=True)
    ap.add_argument("--pred-in", required=True)
    ap.add_argument("--out", required=True, help="output text file")
    ap.add_argument("--truth-end", choices=["inclusive", "exclusive"], default="inclusive", help="whether the end position of each ground truth range is inclusive (default) or exclusive")
    ap.add_argument("--chrom", help="chromosome to use if the map has several")
    a = ap.parse_args()

    method = resolve_method(a.method)
    if method is None:
        die(f"unknown method '{a.method}'. Available: {', '.join(METHODS)}")

    sites, map_order = read_sites(a.map_in, a.chrom)
    truth, truth_pops = read_truth(a.truth_in)
    pred, pred_positions = read_prediction(method, a.pred_in, a.chrom, map_order)
    warnings = []

    # Cross-check the method's own position list against the map, if it has one.
    if pred_positions is not None:
        pset, sset = set(pred_positions), set(sites)
        if pset != sset:
            warnings.append(f"positions in the map ({len(sset)}) and in the prediction's position list ({len(pset)}) differ: {len(sset - pset)} only in map, {len(pset - sset)} only in prediction. The map positions are used for evaluation.")

    common = [h for h in truth if h in pred]
    if not common:
        die("no haplotype IDs are shared between the ground truth and the prediction. Examples - truth: {list(truth)[:3]}, prediction: {list(pred)[:3]}")
    missing_pred = [h for h in truth if h not in pred]
    extra_pred = [h for h in pred if h not in truth]
    if missing_pred:
        warnings.append(f"{len(missing_pred)} truth haplotype(s) have no prediction and are counted as entirely incorrect (e.g. {missing_pred[:3]})")
    if extra_pred:
        warnings.append(f"{len(extra_pred)} predicted haplotype(s) are not in the truth and are ignored (e.g. {extra_pred[:3]})")

    truth_names = {seg[2] for segs in truth.values() for seg in segs}
    pred_names = {seg[2] for segs in pred.values() for seg in segs}
    if pred_names - truth_names:
        warnings.append("predicted population name(s) absent from the ground truth (always counted as incorrect): {sorted(pred_names - truth_names)}")

    rows, tot_correct, tot_eval, partial = [], 0, 0, 0
    anc_names = list(dict.fromkeys(truth_pops.values()))     # header order
    anc_ok = {x: 0 for x in anc_names}                       # correct positions per ancestry
    anc_n = {x: 0 for x in anc_names}                        # true positions per ancestry
    anc_hap_rates = {x: [] for x in anc_names}               # per-haplotype rates per ancestry
    hap_anc = {}                                             # hap -> {ancestry: accuracy}
    incl = INCLUSIVE_END.get(method, False)
    for hap, tsegs in truth.items():
        t_iv = to_index_intervals(tsegs, sites, a.truth_end == "inclusive", f"ground truth haplotype '{hap}'")
        n_eval = sum(hi - lo for lo, hi, _ in t_iv)
        if n_eval < len(sites):
            partial += 1
        if hap in pred:
            p_iv = to_index_intervals(pred[hap], sites, incl, f"prediction haplotype '{hap}'")
            ok_by = count_correct(t_iv, p_iv)
        else:
            p_iv = []
            ok_by = {}
        n_ok = sum(ok_by.values())
        n_by = {}
        for lo, hi, name in t_iv:
            n_by[name] = n_by.get(name, 0) + (hi - lo)
        hap_anc[hap] = {name: ok_by.get(name, 0) / cnt for name, cnt in n_by.items()}
        for name, cnt in n_by.items():
            anc_n[name] += cnt
            anc_ok[name] += ok_by.get(name, 0)
            anc_hap_rates[name].append(ok_by.get(name, 0) / cnt)
        rows.append((hap, n_ok, n_eval))
        tot_correct += n_ok
        tot_eval += n_eval
    if partial:
        warnings.append(f"{partial} haplotype(s) have ground truth covering fewer than all {len(sites)} map positions; uncovered positions are excluded")
    if tot_eval == 0:
        die("no positions were evaluated (do the ground truth ranges overlap the map positions?)")

    overall = tot_correct / tot_eval
    rates = [(hap, ok / n) for hap, ok, n in rows if n > 0]
    if len(rates) < len(rows):
        warnings.append(f"{len(rows) - len(rates)} haplotype(s) have no evaluable positions and are excluded from the per-haplotype list and the average accuracy")
    average = sum(r for _, r in rates) / len(rates)
    with open(a.out, "w") as out:
        out.write("#Summary\n")
        out.write(f"method\t{method}\n")
        out.write(f"map_file\t{a.map_in}\n")
        out.write(f"ground_truth_file\t{a.truth_in}\n")
        out.write(f"predicted_file\t{a.pred_in}\n")
        out.write(f"number_of_haplotypes\t{len(truth)}\n")
        out.write(f"number_of_positions\t{len(sites)}\n")
        out.write(f"number_of_populations\t{len(truth_pops)}\n")
        out.write(f"overall_accuracy\t{overall:.6f}\n")
        out.write(f"average_accuracy\t{average:.6f}\n")
        for w in warnings:
            out.write(f"# WARNING: {w}\n")
        out.write("#Ancestry\tOverall_Accuracy\tAverage_Accuracy\n")
        for x in anc_names:
            if anc_n[x] == 0:
                out.write(f"{x}\tNA\tNA\n")
            else:
                avg_x = sum(anc_hap_rates[x]) / len(anc_hap_rates[x])
                out.write(f"{x}\t{anc_ok[x] / anc_n[x]:.6f}\t{avg_x:.6f}\n")
        out.write("#Haplotype_ID\tOverall_Accuracy"
                  + "".join(f"\t{x}_Accuracy" for x in anc_names) + "\n")
        for hap, r in rates:
            cols = [f"{hap_anc[hap][x]:.6f}" if x in hap_anc[hap] else "NA" for x in anc_names]
            out.write(f"{hap}\t{r:.6f}" + "".join("\t" + c for c in cols) + "\n")

    for w in warnings:
        print(f"WARNING: {w}", file=sys.stderr)
    print(f"[{method}] overall accuracy = {overall:.6f}, average accuracy = {average:.6f} ({len(truth)} haplotypes) -> {a.out}")
    sys.exit(0)
