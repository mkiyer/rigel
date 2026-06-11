# The strand deconvolution and its prior — a ground-up explainer

A node has **N** unspliced fragments, split by orientation into **`pos`** (matching transcript sense)
and **`neg`** (antisense). They are a mix of two populations:

- **gDNA** — genomic, unstranded. Its sense rate is exactly **½** (a gDNA fragment is equally likely
  to look pos or neg).
- **RNA** (mature-unspliced + nascent) — stranded. Its sense rate is **`ss`** (the library's RNA sense
  fraction). `ss = 1` ⇒ all RNA is pos; `ss = 0` ⇒ all RNA is neg (a reverse-stranded library);
  `ss = ½` ⇒ RNA is 50/50, indistinguishable from gDNA.

We want **`g`** = the gDNA *fraction* of the N fragments (so gDNA count = `g·N`, RNA count = `(1−g)·N`).
"Strand-specificity" in the intuitive sense is the **distance of `ss` from ½**, i.e. `|2·ss − 1|`: it is
1 at `ss=0` *or* `ss=1` (both fully specific, just opposite orientations) and 0 at `ss=½`.

## 1. How the counts inform `g` — the mixing line

If a fraction `g` of the fragments are gDNA (sense rate ½) and `1−g` are RNA (sense rate `ss`), the
node's overall expected sense rate is a straight line in `g`:

```
p(g) = ½·g + ss·(1−g) = ss + (½ − ss)·g
```

- At `g = 1` (all gDNA): `p = ½`.
- At `g = 0` (all RNA): `p = ss`.

So the node's *observed* sense fraction `s = pos/N` is read off this line. Inverting it:

```
g = (s − ss) / (½ − ss)
```

This is the whole deconvolution in one formula — **and it already shows the trouble.** The denominator
is `½ − ss`. As `ss → ½`, that denominator → 0, and `g` blows up: a tiny change in `s` swings `g`
across all of [0,1]. **At `ss = ½` exactly the formula is 0/0 — undefined.** That is the mathematical
face of "50/50 at `ss=½` is uninterpretable."

The reason is visible in the line: at `ss = ½`, `p(g) = ½ + 0·g = ½` for *every* `g`. gDNA and RNA have
the *same* sense rate, so no observed split can tell them apart. The data carry **zero information**
about `g`.

## 2. The likelihood — what the data actually say

We don't just invert a point; we score every candidate `g` by how well it predicts the observed split.
For `pos` sense out of `N`, the likelihood is Binomial:

```
L(g) ∝ p(g)^pos · (1 − p(g))^neg          # in rigel, Beta-Binomial — overdispersed; same shape
log L(g) = pos·log p(g) + neg·log(1 − p(g))
```

For a balanced node (`pos = neg = 50`, `N = 100`) this is maximised where `p(g) = ½` — i.e. at
`g = (½ − ss)/(½ − ss) = 1` for any `ss ≠ ½`. So **a perfectly balanced node always points at "all
gDNA"** (only gDNA produces a 50/50 split; RNA produces an `ss`-skewed split). What changes with `ss`
is not *where* the likelihood peaks but **how sharp** that peak is.

How sharp? The curvature of `log L` at the peak is the **Fisher information** of `g`:

```
I(g) = N · (dp/dg)² / [p(1−p)] = N · (½ − ss)² / [p(1−p)]
```

At the balanced peak (`p = ½`): `I = N · (½ − ss)² / ¼ = N · (2·ss − 1)²`. So the standard error of the
deconvolved `g` is

```
SE(g) ≈ 1 / sqrt(I) = 1 / ( 2·sqrt(N) · |½ − ss| )
```

**This is the key quantity.** The deconvolution's uncertainty is governed entirely by `|½ − ss|`
(how specific the library is) and `√N` (how many fragments). It → ∞ as `ss → ½`.

## 3. The prior — what `Beta(½, ½)` is and what it does

The likelihood alone has two problems: at `ss = ½` it's perfectly flat (no answer), and near `ss = ½`
the raw inversion `g = (s−ss)/(½−ss)` can fly outside [0,1] on noisy data. A **prior** on `g` fixes
both. Rigel uses **`Beta(½, ½)`** — the Jeffreys prior for a proportion. The posterior is

```
posterior(g) ∝ Beta(g; ½, ½) · L(g) ,   read off on a grid of g ∈ [0,1]; report the median.
```

What `Beta(½, ½)` looks like and why:

- **U-shaped.** Its density goes to **infinity at `g = 0` and `g = 1`** and dips to a minimum at
  `g = ½`. It places *more* prior weight on the extremes than on the middle.
- **Why U-shaped, not flat?** A flat prior (`Beta(1,1)`, uniform) sounds "uninformative" but isn't,
  in an invariance sense: a Binomial carries more information per observation near the extremes (where
  its variance is small) than near ½. The Jeffreys prior is the unique one that is **invariant to how
  you parameterise the problem**, and it compensates by weighting the extremes. Intuitively here it
  encodes a mild belief that *a node is usually mostly-one-thing* — predominantly gDNA **or**
  predominantly RNA — rarely a clean half-and-half mixture.
- **It is symmetric**, so its **median is exactly ½** even though its *mass* sits at the ends. Hold on
  to this distinction (§5) — it matters.

What the prior *does*, concretely:

1. **Regularises the inversion** — keeps `g` in [0,1] and tames the near-singular `1/(½−ss)` blow-up; an
   "impossible" value (e.g. `g=0` when `ss=1`, which would need RNA to produce antisense reads) gets
   likelihood 0 and is excluded cleanly.
2. **Is the entire answer when the data are silent** — at `ss = ½` the likelihood is flat, so
   `posterior = prior`. The strand module falls back to its prior belief.
3. **Vanishes when the data speak** — once `ss` is away from ½ and `N` is decent, `L(g)` is far sharper
   than the prior, and the posterior is driven by the data. With `N = 100` the prior only matters in a
   thin band around `ss = ½`.

## 4. Walk-through: `N = 100`, `pos = neg = 50`, swept over `ss`

| `ss` | RNA sense rate | line `p(g)` | peak `g` | per-fragment info `(2ss−1)²` | `SE(g)=1/(20·\|½−ss\|)` | what the module knows |
|---|---|---|---|---|---|---|
| **0.5** | 0.50 | `≡ ½` (flat) | undefined | **0.00** | **∞** | **Nothing.** posterior = prior. |
| 0.45 | 0.45 | `0.45+0.05g` | 1.0 | 0.01 | 1.00 | Almost nothing; SE spans all of [0,1]. |
| 0.2 | 0.20 | `0.20+0.30g` | 1.0 | 0.36 | 0.17 | Good; confidently *mostly* gDNA. |
| 0.01 | 0.01 | `0.01+0.49g` | 1.0 | 0.96 | 0.10 | Excellent; ≈ all gDNA. |
| 0.999 | 0.999 | `0.999−0.499g` | ≈1.0 | 0.996 | 0.10 | Excellent; ≈ all gDNA. |
| **1.0** | 1.00 | `1−0.5g` | 1.0 | **1.00** | **0.10** | **Perfect; all 100 are gDNA.** |

Read it as a story. **`ss = 1.0`:** RNA is *all* pos, so the 50 **neg** reads *cannot* be RNA — they
must be gDNA; gDNA is 50/50, so the matching 50 pos are gDNA too → **all 100 gDNA**, tight (SE 0.10).
This is your intuition, exactly. **`ss = 0.01`** (reverse-stranded) and **`ss = 0.2`** tell the same
story with a little more noise — a 50/50 split is far more gDNA-like than RNA-like, so still "mostly/all
gDNA," just wider. **`ss = 0.45`:** gDNA (½) and RNA (0.45) are barely different; SE = 1.0 means the data
can't pin `g` at all → the prior leads. **`ss = 0.5`:** the line is flat; the data say nothing.

Notice the peak is at `g = 1` for *every* `ss ≠ ½` — what moves is the **confidence** (the SE column),
collapsing from ∞ at `ss=½` to 0.10 at `ss=1`. The deconvolution doesn't change its *guess* with `ss`;
it changes how *sure* it is.

## 5. The subtlety at `ss = ½` — the point estimate lies, the variance tells the truth

At `ss = ½` the posterior is the bare `Beta(½,½)`: its **median is ½**, but its **density is bimodal**
— piled at 0 and 1. So the honest reading is *not* "this node is 50% gDNA." It is **"this node is
almost surely *pure* — either all gDNA or all RNA — and I cannot tell which."** The median ½ is an
artefact of symmetry, not a belief in a half-mixture.

This is the crux for any pipeline that *passes the deconvolution onward*: **the point estimate alone is
misleading; you must carry the variance.** At `ss = ½` the variance is maximal (SE = ∞), which is the
module's way of saying *"ignore me, defer to whatever else you have"* (the count module). At `ss` near 0
or 1 the variance is small (SE ≈ 0.10), saying *"trust me."*

## 6. Why this *is* the weight `w`

Look again at the per-fragment information column: it is exactly **`(2·ss − 1)²`** — which is the blend
weight `w` used elsewhere in calibration. That is not a coincidence:

```
posterior precision (for a balanced node) = I(g) = N · (2·ss − 1)² = N · w
```

So **`w` is the per-fragment Fisher information of the deconvolution**, and `N·w` is its precision. The
"how much do I trust strand" weight and the "how uncertain is the deconvolution" variance are *the same
object*, seen from two sides. The deconvolution precision the downstream count module needs is `N·w`:
near zero at `ss=½` (defer entirely), large at high specificity (lean on strand).

## 7. Summary

- The deconvolution is the mixing line `p(g) = ss + (½−ss)·g`, inverted; it is singular at `ss = ½`.
- The likelihood always *points* at all-gDNA for a balanced node; `|½−ss|·√N` sets how *sharp* that is.
- `Beta(½,½)` regularises the inversion, supplies the answer when the data are silent (`ss=½`), and
  fades when they aren't. It is U-shaped (mass at the ends, "nodes are usually pure") but symmetric
  (median ½).
- At `ss=½` the posterior is bimodal "pure-but-unknown-which," **not** a half-mixture — so propagate the
  **variance**, never the point estimate alone.
- The deconvolution precision is `N·(2·ss−1)²` — identical to the blend weight `w`. Uncertainty and
  weight are one quantity.
