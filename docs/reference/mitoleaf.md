# `mitoleaf`

```jennifer
import "@deckyard/mitoleaf/" as mitoleaf;
```

Signatures, constants and struct fields are in the
[cheatsheet](cheatsheet.md). This page covers input rules and errors.

## Input

`parse` reads the phm XML as downloaded: comments and wrapper elements around
the root `<Node>` are skipped. The document holds exactly one top-level
`<Node>`. The element name matches in any case; the
attributes are read as `Id`, `id` or `ID` and `HG`, `hg` or `Hg`. A missing
`HG` is an empty motif. Motif tokens are split on spaces and stored verbatim
in `Node.motif`, and parsed once into `Tree.ops`. A malformed token is not
rejected by `parse`; it throws when `haplotype`, `distance` or `assign` folds
it.

## Building a `Tree` by hand

`Tree.ops` caches each haplogroup's parsed motif. `fromNodes(root, nodes)`
builds a `Tree` with `ops` filled from the motifs; call it again after
changing a motif, or `ops` goes stale and results follow the old motif. A
`Tree{...}` literal must name `ops` (since 0.3.0); `ops: {}` is valid, and a
haplogroup missing from `ops` has its motif parsed on every use.

`parseMotifItem` trims and upper-cases a token before reading it. A back
mutation names a reference position: `A<position>` and
`-<position>.<insertion>` take a position in [1, 16569] and an insertion index
in [1, 9999]. A back mutation with nothing to revert is ignored, not rejected,
because published trees restate some. `-<position>.<insertion>` drops the
highest-indexed inserted base still present at that position, so a partial
loss shortens an insertion run from its end.

An observed haplotype passed to `distance` or `assign` carries at most one
designation per position and insertion index: `263G 263G` or `263G 263A` would
count twice in `lineage.mtDistance` and is rejected.

Haplogroup names match exactly, case included.

## Errors

Errors are `Error{kind: "mitoleaf"}`, except an invalid motif gain, which
raises `kind: "lineage"` from `lineage.parseDifference`.

| Condition | `kind` | Thrown by |
| --- | --- | --- |
| the document is not well-formed XML | `mitoleaf` | `parse` |
| the document holds no `<Node>` | `mitoleaf` | `parse` |
| more than one top-level `<Node>` | `mitoleaf` | `parse` |
| a `<Node>` has no `Id` | `mitoleaf` | `parse` |
| a haplogroup name appears twice | `mitoleaf` | `parse` |
| unknown haplogroup | `mitoleaf` | `node`, `children`, `motif`, `path`, `haplotype`, `distance` |
| empty token | `mitoleaf` | `parseMotifItem` |
| a leading-minus token not of the form `-<position>.<insertion>` | `mitoleaf` | `parseMotifItem`, `haplotype`, `distance`, `assign` |
| a back mutation at a position outside [1, 16569] or an insertion index outside [1, 9999] | `mitoleaf` | `parseMotifItem`, `haplotype`, `distance`, `assign` |
| the observed haplotype carries two designations at one position | `mitoleaf` | `distance`, `assign` |
| a gain that is not valid rCRS nomenclature | `lineage` | `parseMotifItem`, `haplotype`, `distance`, `assign` |

`has` and `haplogroups` never throw. `assign` on a `Tree` with no nodes (which
`parse` never returns) fails with a runtime error, not a `mitoleaf` error.
