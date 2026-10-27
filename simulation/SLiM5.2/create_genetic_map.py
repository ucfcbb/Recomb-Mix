# This program was generated with the assistance of Claude
import sys
import argparse
import bisect


def read_header(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\r\n")
    if not header:
        sys.exit(f"ERROR: {path} is empty")
    return header


def same_chrom(a, b):
    norm = lambda c: (c[3:] if c.lower().startswith("chr") else c).lower()
    return norm(a) == norm(b)


def read_map(path):
    rows = []
    with open(path) as fh:
        header = fh.readline()
        if not header:
            sys.exit(f"ERROR: {path} is empty")
        for n, line in enumerate(fh, start=2):
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 4:
                sys.exit(f"ERROR: {path} line {n}: expected 4 tab-separated columns, got {len(f)}")
            try:
                float(f[1]); float(f[2]); float(f[3])
            except ValueError:
                sys.exit(f"ERROR: {path} line {n}: non-numeric value in position/rate/map columns")
            rows.append({"chrom": f[0].strip(), "pos": f[1].strip(), "rate": f[2].strip(), "cm": f[3].strip()})
    return rows


def fix_chrom(name, style):
    if style == "as-is":
        return name
    bare = name[3:] if name.lower().startswith("chr") else name
    return "chr" + bare if style == "add-chr" else bare


def check_map_sorted(rows):
    last = {}
    for r in rows:
        p, c = float(r["pos"]), float(r["cm"])
        prev = last.get(r["chrom"])
        if prev and (p <= prev[0] or c <= prev[1]):
            print(f"WARNING: map position or cM not strictly increasing on chromosome {r['chrom']} at position {r['pos']}", file=sys.stderr)
            return
        last[r["chrom"]] = (p, c)


def read_positions(path):
    with open(path) as fh:
        lines = [ln.strip() for ln in fh.read().splitlines()]
    while lines and not lines[-1]:
        lines.pop()
    if not lines:
        sys.exit(f"ERROR: {path} is empty")
    for n, ln in enumerate(lines, start=1):
        try:
            float(ln)
        except ValueError:
            sys.exit(f"ERROR: {path} line {n}: '{ln}' is not a valid position (blank lines and headers are not allowed)")
    vals = [float(x) for x in lines]
    for n in range(1, len(vals)):
        if vals[n] <= vals[n - 1]:
            sys.exit(f"ERROR: {path} line {n + 1}: positions must be strictly increasing (duplicate or out-of-order position '{lines[n]}')")
    return lines


def to_int(text, what):
    try:
        v = float(text)
    except ValueError:
        sys.exit(f"ERROR: {what} '{text}' is not a number")
    if v != int(v):
        sys.exit(f"ERROR: {what} '{text}' is not an integer")
    return int(v)


def map_to_positions(gmap, positions, decimals=12, rate_decimals=6, flat_step=1e-9):
    """Return one map row per entry of `positions`, in the same order.

    Step 1 - base cM of every position:
      * position found in the dictionary -> the ORIGINAL cM text of the map is copied;
      * position between two map positions -> linear interpolation
            cM = cM_next + (pos - pos_next) * (cM_next - cM_prev) / (pos_next - pos_prev);
      * position below the first / above the last map position -> the cM of the
        first / last map position.

    Step 2 - strictly increasing cM. Step 1 gives equal cM values for positions outside
      the map range and inside flat map intervals. Positions that are NOT found in the
      map are adjusted (never the ones copied from the map) so that cM increases by at
      least `flat_step` cM from row to row: a forward pass raises tied values, a backward
      pass lowers values that would otherwise reach the next map-copied value.
      Only positions that would be equal (or too close) to a neighbour are changed.

    Rate(cM/Mb) is computed from the WRITTEN rows:
        rate[i] = (cM[i+1]-cM[i]) / ((pos[i+1]-pos[i]) / 1e6); the last row gets 0.
    """
    chrom = gmap[0]["chrom"]
    genetic_maps = {}
    for r in gmap:
        genetic_maps.setdefault(to_int(r["pos"], "map position"), (float(r["cm"]), r["cm"]))
    xs = sorted(genetic_maps)

    pv = [to_int(p, "position") for p in positions]
    n = len(pv)
    base, text, anchor = [0.0] * n, [None] * n, [False] * n
    n_exact = n_interp = n_outside = 0
    for i, p in enumerate(pv):
        hit = genetic_maps.get(p)
        if hit is not None:                              # exact: copy the map's cM
            base[i], text[i], anchor[i] = hit[0], hit[1], True
            n_exact += 1
            continue
        j = bisect.bisect_right(xs, p)                   # upper_bound
        if j == len(xs):                                 # above the last map position
            base[i] = genetic_maps[xs[-1]][0]
            n_outside += 1
        elif j == 0:                                     # below the first map position
            base[i] = genetic_maps[xs[0]][0]
            n_outside += 1
        else:                                            # linear interpolation
            nx, pr = xs[j], xs[j - 1]
            gn, gp = genetic_maps[nx][0], genetic_maps[pr][0]
            base[i] = gn + (p - nx) * (gn - gp) / (nx - pr)
            n_interp += 1

    for i in range(1, n):
        if base[i] < base[i - 1] - 1e-9:
            sys.exit(f"ERROR: cM decreases between positions {pv[i - 1]} and {pv[i]}; "
                     "the input map must have non-decreasing cM.")

    # --- Make cM strictly increasing without touching map-copied values ---
    w = list(base)
    for i in range(1, n):                                # forward: lift ties
        if not anchor[i]:
            w[i] = max(w[i], w[i - 1] + flat_step)
    for i in range(n - 2, -1, -1):                       # backward: stay below next value
        if not anchor[i]:
            w[i] = min(w[i], w[i + 1] - flat_step)
    changed = [i for i in range(n) if not anchor[i] and w[i] != base[i]]
    max_change = max((abs(w[i] - base[i]) for i in changed), default=0.0)

    if n and min(w) < 0:
        i = w.index(min(w))
        sys.exit(f"ERROR: making cM strictly increasing before position {pv[i]} would need a negative cM value (too many positions below the first map-copied value).")
    for i in range(1, n):
        if not w[i] > w[i - 1]:
            sys.exit(f"ERROR: cM cannot be made strictly increasing at positions {pv[i - 1]} -> {pv[i]} (cM {w[i - 1]:.12g} -> {w[i]:.12g}). The two values are fixed by the map (equal or too close) or --flat-step is 0.")

    # --- format: map-copied rows keep their original text; others share `d` decimals ---
    for d in range(decimals, 16):
        for i in range(n):
            text[i] = text[i] if anchor[i] else f"{w[i]:.{d}f}"
        back = [float(t) for t in text]
        if all(back[i] > back[i - 1] for i in range(1, n)):
            break
    else:
        sys.exit("ERROR: cM values cannot be written strictly increasing with up to 15 decimals.")

    notes = []
    if changed:
        notes.append(f"{len(changed)} position(s) not in the map (outside the map range or in flat map intervals) were adjusted by at most {max_change:.3g} cM to keep cM strictly increasing; the {n_exact} positions found in the map keep their original cM")
    if d > decimals:
        notes.append(f"cM written with {d} decimals (instead of {decimals})")

    # --- Rate computed from the written values ---
    out = []
    for i, ptxt in enumerate(positions):
        rate = (back[i + 1] - back[i]) / ((pv[i + 1] - pv[i]) / 1e6) if i < n - 1 else 0.0
        out.append({"chrom": chrom, "pos": ptxt, "cm": text[i], "rate": f"{rate:.{rate_decimals}f}"})
    return out, n_exact, n_interp, n_outside, notes


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--map-in", required=True)
    ap.add_argument("--positions-in", required=True, help="One position (bp) per line, no header, strictly increasing")
    ap.add_argument("--map-out", required=True)
    ap.add_argument("--chrom", help="Keep only this chromosome of the input map (e.g. 21 or chr21); required if the map has several chromosomes")
    ap.add_argument("--chrom-format", choices=["short", "long"], default=None, help="Chromosome column of the output: short (21) or long (chr21). Default: as written in the input map")
    ap.add_argument("--flat-step", type=float, default=1e-9, help="Minimum cM increase between consecutive output rows, used for positions that are not in the map and would otherwise get equal cM (outside the map range, flat map intervals). Default: 1e-9 cM")
    ap.add_argument("--decimals", type=int, default=12, help="Minimum decimal places for interpolated cM values (default: 12; the rate always uses 6). More are used automatically, up to 15, if needed to keep cM strictly increasing")
    a = ap.parse_args()

    header = read_header(a.map_in)
    gmap = read_map(a.map_in)
    if a.chrom is not None:
        gmap = [r for r in gmap if same_chrom(r["chrom"], a.chrom)]
        if not gmap:
            sys.exit(f"ERROR: no rows for chromosome '{a.chrom}' in {a.map_in}")
    check_map_sorted(gmap)
    chroms = sorted({r["chrom"] for r in gmap})
    if len(chroms) > 1:
        sys.exit("ERROR: the position file has no chromosome column, but the map contains {len(chroms)} chromosomes; select one with --chrom.")
    if not gmap:
        sys.exit(f"ERROR: no data rows in {a.map_in}")

    positions = read_positions(a.positions_in)
    rows, n_exact, n_interp, n_outside, notes = map_to_positions(gmap, positions, a.decimals, flat_step=a.flat_step)
    for note in notes:
        print(f"NOTE: {note}", file=sys.stderr)
    if n_outside:
        print(f"WARNING: {n_outside} position(s) outside the map range (set to the first/last map cM)", file=sys.stderr)

    style = {None: "as-is", "short": "strip-chr", "long": "add-chr"}[a.chrom_format]
    with open(a.map_out, "w") as out:
        out.write(header + "\n")
        for r in rows:
            out.write("\t".join([fix_chrom(r["chrom"], style), r["pos"], r["rate"], r["cm"]]) + "\n")
    print(f"positions: {n_exact} found in map, {n_interp} interpolated, {n_outside} outside map range; wrote {len(rows)} map rows -> {a.map_out}")
    sys.exit(0)
    