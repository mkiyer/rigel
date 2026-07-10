"""Phase-0 message-precision baseline: bucket the benchmark_ab_report JSON into OFF-vs-prod.

The question: does the current global-scalar message layer net-HELP or net-HURT vs no messages,
per (capture x strand x gDNA x nascent) bucket, now that the nascent siphon is fixed?

Pools (assigned - true fragments): gDNA / nascent(=SIPHON, synthetic-only) / mature. Ideal = 0.
- gdna_surplus < 0  => gDNA under-assigned => gDNA->RNA leak.
- nrna_surplus       => nascent siphon (on nrna_none conditions true=0, so it IS the pure siphon).
Message value = |error(off)| - |error(prod)|  (POSITIVE => messages HELP; NEGATIVE => messages HURT).

    python phase0_summarize.py phase0_baseline.json
"""
import json, sys, re
import numpy as np

d = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "phase0_baseline.json"))
arms = d["arms"]                      # ["prod", "off"]
conds = d["conditions"]


def parse(name):
    g = "gdna300" if "gdna300" in name else "gNONE"
    ss = "ss99" if "0.99" in name else "ss50"
    nr = "nrRND" if "nrna_rnd" in name else "nrNONE"
    cap = "capON" if "capture_on" in name else "capOFF"
    return g, ss, nr, cap


def pool(arm, cond, p):
    return d["data"][arm][cond]["pools"][p]["surplus"]


print(f"{'condition':52} | {'gdna surplus':>18} | {'SIPHON (nrna)':>18} | {'mature surplus':>18}")
print(f"{'':52} | {'prod':>8} {'off':>8} | {'prod':>8} {'off':>8} | {'prod':>8} {'off':>8}")
print("-" * 120)
rows = []
for c in sorted(conds, key=parse):
    g, ss, nr, cap = parse(c)
    gp, go = pool("prod", c, "gdna"), pool("off", c, "gdna")
    np_, no = pool("prod", c, "nrna"), pool("off", c, "nrna")
    mp, mo = pool("prod", c, "mrna"), pool("off", c, "mrna")
    rows.append((g, ss, nr, cap, gp, go, np_, no, mp, mo))
    tag = f"{g:8} {ss:5} {nr:7} {cap:6}"
    print(f"{tag:52} | {gp:>8.0f} {go:>8.0f} | {np_:>8.0f} {no:>8.0f} | {mp:>8.0f} {mo:>8.0f}")

R = np.array([r[4:] for r in rows], float)  # gp,go,np,no,mp,mo
keys = [r[:4] for r in rows]


def bucket(mask, label):
    if not mask.any():
        return
    sub = R[mask]
    # message value on the gDNA/RNA split (|gdna err|) and the siphon (|nrna err|): |off| - |prod|
    dg = np.abs(sub[:, 1]) - np.abs(sub[:, 0])   # |go|-|gp|
    dn = np.abs(sub[:, 3]) - np.abs(sub[:, 2])   # |no|-|np|
    print(f"  {label:34} n={mask.sum():2}  "
          f"|gdna| prod={np.abs(sub[:,0]).mean():8.0f} off={np.abs(sub[:,1]).mean():8.0f} "
          f"(msg {'HELPS' if dg.mean()>0 else 'hurts'} {dg.mean():+8.0f}) | "
          f"SIPHON prod={sub[:,2].mean():8.0f} off={sub[:,3].mean():8.0f} "
          f"(msg {'HELPS' if dn.mean()>0 else 'hurts'} {dn.mean():+8.0f})")


kg = np.array([k[0] for k in keys]); kss = np.array([k[1] for k in keys])
knr = np.array([k[2] for k in keys]); kcap = np.array([k[3] for k in keys])
print("\n=== BUCKETED message value ( |off| - |prod| ; POSITIVE => messages help ) ===")
print("-- by strand --")
bucket(kss == "ss99", "STRANDED (ss99)")
bucket(kss == "ss50", "UNSTRANDED (ss50)")
print("-- by strand x capture --")
for ss in ("ss99", "ss50"):
    for cap in ("capOFF", "capON"):
        bucket((kss == ss) & (kcap == cap), f"{ss} {cap}")
print("-- by gDNA level (gdna300 only has gDNA to relay) --")
bucket(kg == "gdna300", "gdna300")
bucket(kg == "gNONE", "gNONE (no gDNA)")
print("-- overall --")
bucket(np.ones(len(rows), bool), "ALL 16")
