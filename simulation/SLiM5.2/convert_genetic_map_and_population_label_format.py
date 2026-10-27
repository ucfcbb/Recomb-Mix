# This program was generated with the assistance of Claude
import sys
import argparse
from dataclasses import dataclass
from typing import List, Optional


# ----------------------------------------------------------------------------
# Method specifications
# ----------------------------------------------------------------------------
# Map columns are chosen from these source fields:
#   "chrom" : chromosome           "pos" : position (bp)
#   "rate"  : rate (cM/Mb)         "cm"  : cumulative map (cM)
# A column starting with "=" is a constant (e.g. "=." for a placeholder ID column).
# Label columns are chosen from: "sample", "label".


@dataclass
class MethodSpec:
    map_columns: List[str]
    map_header: Optional[str] = None
    map_delim: str = "\t"
    label_columns: Optional[List[str]] = None
    label_header: Optional[str] = None
    label_delim: str = "\t"
    notes: str = ""

    def __post_init__(self):
        if self.label_columns is None:
            self.label_columns = ["sample", "label"]


METHODS = {
    "rfmix": MethodSpec(
        map_columns=["chrom", "pos", "cm"],
        map_header="#(RFmix)Chromosome\tPosition\tGenetic_Map(cM)",
        label_header="#(RFmix)Sample\tPanel",
        notes="Tab-delimited; 3-column map with header; labels with header.",
    ),
    "flare": MethodSpec(
        map_columns=["chrom", "rate", "cm", "pos"],
        map_delim=" ",
        label_delim="\t",
        notes="Map: chrom, rate, cM, pos (space-delimited, no header); "
              "labels: sample, label (tab-delimited, no header).",
    ),
    "gnomix": MethodSpec(
        map_columns=["chrom", "pos", "cm"],
        map_header="#(Gnomix)Chromosome\tPosition\tGenetic_Map(cM)",
        label_header="#(Gnomix)Sample\tPanel",
        notes="Tab-delimited; 3-column map with header; labels with header.",
    ),
    "recombmix": MethodSpec(
        map_columns=["chrom", "pos", "rate", "cm"],
        map_header="#Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)",
        label_header="#Sample_ID\tAncestry_Label",
        notes="Tab-delimited; map in input column order with '#' header; labels with header.",
    ),
    "sparsepainter": MethodSpec(
        map_columns=["pos", "cm"],
        map_header="position(base)\tgenetic_distance(centiMorgan)",
        notes="Tab-delimited; 2-column map (pos, cM) with header; labels: sample, label, no header.",
    ),
}


def resolve_method(name):
    """Case-insensitive; ignores '-', '_' and spaces."""
    norm = "".join(ch for ch in name.lower() if ch not in "-_ ")
    return norm if norm in METHODS else None


# ----------------------------------------------------------------------------
# Parsing and conversion
# ----------------------------------------------------------------------------


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


def read_labels(path):
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
            if len(f) < 2:
                sys.exit(f"ERROR: {path} line {n}: expected 2 tab-separated columns, got {len(f)}")
            rows.append({"sample": f[0].strip(), "label": f[1].strip()})
    return rows


def fix_chrom(name, style):
    if style == "as-is":
        return name
    bare = name[3:] if name.lower().startswith("chr") else name
    return "chr" + bare if style == "add-chr" else bare


def write_table(path, rows, columns, delim, header):
    with open(path, "w") as out:
        if header is not None:
            out.write(header + "\n")
        for r in rows:
            vals = [c[1:] if c.startswith("=") else r[c] for c in columns]
            out.write(delim.join(vals) + "\n")


def check_map_sorted(rows):
    last = {}
    for r in rows:
        p, c = float(r["pos"]), float(r["cm"])
        prev = last.get(r["chrom"])
        if prev and (p <= prev[0] or c <= prev[1]):
            print(f"WARNING: map position or cM not strictly increasing on chromosome {r['chrom']} at position {r['pos']}", file=sys.stderr)
            return
        last[r["chrom"]] = (p, c)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--method", help="LAI method name (see --list-methods)")
    ap.add_argument("--map-in")
    ap.add_argument("--labels-in")
    ap.add_argument("--map-out")
    ap.add_argument("--labels-out")
    ap.add_argument("--chrom", help="Keep only this chromosome (as named in the input map)")
    ap.add_argument("--chrom-style", choices=["as-is", "add-chr", "strip-chr"], default="as-is", help="Rename chromosomes in the output map")
    ap.add_argument("--list-methods", action="store_true")
    a = ap.parse_args()

    if a.list_methods:
        for name, s in METHODS.items():
            print(f"{name:10s} {s.notes}")
    else:
        missing = [k for k in ("method", "map_in", "labels_in", "map_out", "labels_out") if getattr(a, k) is None]
        if missing:
            ap.error("missing required arguments: " + ", ".join("--" + m.replace("_", "-") for m in missing))

        key = resolve_method(a.method)
        if key is None:
            sys.exit(f"ERROR: unknown method '{a.method}'. Available: {', '.join(METHODS)}")
        spec = METHODS[key]

        gmap = read_map(a.map_in)
        if a.chrom is not None:
            gmap = [r for r in gmap if r["chrom"] == a.chrom]
            if not gmap:
                sys.exit(f"ERROR: no rows for chromosome '{a.chrom}' in {a.map_in}")
        check_map_sorted(gmap)
        for r in gmap:
            r["chrom"] = fix_chrom(r["chrom"], a.chrom_style)

        labels = read_labels(a.labels_in)

        write_table(a.map_out, gmap, spec.map_columns, spec.map_delim, spec.map_header)
        write_table(a.labels_out, labels, spec.label_columns, spec.label_delim, spec.label_header)
        print(f"[{key}] wrote {len(gmap)} map rows -> {a.map_out}; {len(labels)} labels -> {a.labels_out}")
    sys.exit(0)
