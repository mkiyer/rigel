"""Phase-1: bucket the 4-mode (prod/off/max/eb) matrix. For each (strand x capture) bucket, the mean
|gdna surplus| (the gDNA/RNA split error) and the mean PURE siphon (nrNONE only, where true nascent=0)
per mode. Best mode per bucket = min |gdna|. Message value vs off/prod annotated."""
import json, sys
import numpy as np

d = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "phase1_matrix.json"))
arms = d["arms"]
conds = d["conditions"]


def parse(name):
    g = "gdna300" if "gdna300" in name else "gNONE"
    ss = "ss99" if "0.99" in name else "ss50"
    nr = "nrRND" if "nrna_rnd" in name else "nrNONE"
    cap = "capON" if "capture_on" in name else "capOFF"
    return g, ss, nr, cap


def sur(arm, c, p):
    return d["data"][arm][c]["pools"][p]["surplus"]


# per-condition gdna surplus table
print("gDNA surplus (assigned-true; <0 = gDNA->RNA leak) by mode")
print(f"{'condition':44} " + " ".join(f"{a:>10}" for a in arms))
for c in sorted(conds, key=parse):
    g, ss, nr, cap = parse(c)
    print(f"{g:8}{ss:5}{nr:7}{cap:6}{'':13} " + " ".join(f"{sur(a,c,'gdna'):>10.0f}" for a in arms))

keys = {c: parse(c) for c in conds}


def bucket(pred, label):
    cs = [c for c in conds if pred(keys[c])]
    if not cs:
        return
    gd = {a: np.mean([abs(sur(a, c, "gdna")) for c in cs]) for a in arms}
    # pure siphon: nrNONE only (true nascent=0 => nrna surplus IS the siphon)
    cs_none = [c for c in cs if keys[c][2] == "nrNONE"]
    sp = {a: np.mean([sur(a, c, "nrna") for c in cs_none]) for a in arms} if cs_none else None
    best = min(arms, key=lambda a: gd[a])
    print(f"\n  {label}  (n={len(cs)}, nrNONE={len(cs_none)})")
    print(f"    |gdna| : " + "  ".join(f"{a}={gd[a]:.0f}" for a in arms) + f"   [best split: {best}]")
    if sp:
        bsp = min(arms, key=lambda a: abs(sp[a]))
        print(f"    siphon : " + "  ".join(f"{a}={sp[a]:+.0f}" for a in arms) + f"   [best siphon: {bsp}]")


print("\n=== BUCKETED (mean |gdna split err| + pure siphon) ===")
for ss in ("ss99", "ss50"):
    for cap in ("capOFF", "capON"):
        bucket(lambda k, ss=ss, cap=cap: k[1] == ss and k[3] == cap, f"{ss} {cap}")
print("\n-- rollups --")
bucket(lambda k: True, "ALL 16")
bucket(lambda k: k[1] == "ss50", "UNSTRANDED (ss50)")
bucket(lambda k: k[1] == "ss99", "STRANDED (ss99)")
