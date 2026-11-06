<!-- SPDX-License-Identifier: LGPL-3.0-only -->
<!-- SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev> -->

# mitoleaf cheatsheet

```jennifer
import "@deckyard/mitoleaf/" as mitoleaf;
```

Errors and input rules: [mitoleaf](mitoleaf.md). `tools/check-docs.sh` checks
every signature, field and value on this page against the code.

## `mitoleaf`

| Function | Returns | Does |
|---|---|---|
| `mitoleaf.parse(doc)` | `Tree` | parses a phm XML document; motifs are pre-parsed into `ops`, not rejected |
| `mitoleaf.fromNodes(root, nodes)` | `Tree` | a `Tree` from hand-built nodes, `ops` filled in |
| `mitoleaf.haplogroups(t)` | `list of string` | every haplogroup, in document order |
| `mitoleaf.has(t, hg)` | `bool` | whether `hg` is in the tree (exact, case-sensitive) |
| `mitoleaf.node(t, hg)` | `Node` | one haplogroup's node |
| `mitoleaf.children(t, hg)` | `list of string` | direct children, in document order |
| `mitoleaf.motif(t, hg)` | `list of string` | the tokens `hg` adds to its parent, verbatim |
| `mitoleaf.path(t, hg)` | `list of string` | root to `hg`, both included |
| `mitoleaf.haplotype(t, hg)` | `lineage.MtHaplotype` | motifs accumulated from the root, back mutations applied |
| `mitoleaf.distance(t, hg, h)` | `int` | `lineage.mtDistance` from `hg`'s haplotype to `h` |
| `mitoleaf.assign(t, h)` | `Match` | nearest haplogroup to `h`; ties go to the first in document order |
| `mitoleaf.parseMotifItem(token)` | `MotifItem` | interprets one motif token |

| Constant | Type | Value | Means |
|---|---|---|---|
| `mitoleaf.WHOLE_GENOME` | `string` | `"1-16569"` | the one range of every phylogeny haplotype |

### `mitoleaf.Node`

| Field | Holds |
|---|---|
| `id as string` | haplogroup name, e.g. `"L0a1b2a"` |
| `parent as string` | parent haplogroup, `""` at the root |
| `motif as list of string` | tokens this node adds to its parent, verbatim |
| `children as list of string` | child haplogroups, in document order |

### `mitoleaf.Tree`

| Field | Holds |
|---|---|
| `root as string` | root haplogroup, `"mt-MRCA"` in the published tree |
| `nodes as map of string to Node` | haplogroup name to node, in document order |
| `ops as map of string to string` | haplogroup name to its motif pre-parsed (internal encoding); built by `parse` / `fromNodes` |

### `mitoleaf.MotifItem`

| Field | Holds |
|---|---|
| `reversion as bool` | true for a back mutation (`A263`, `-5899.1`) |
| `position as int` | rCRS position the token names |
| `insertion as int` | n for the n-th base inserted after `position`, else 0 |
| `designation as string` | normalised rCRS form of a gain; for a back mutation the token, trimmed and upper-cased |

### `mitoleaf.Match`

| Field | Holds |
|---|---|
| `haplogroup as string` | haplogroup name |
| `distance as int` | `lineage.mtDistance` from the observed haplotype |
