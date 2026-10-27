# `mitoleaf`

```jennifer
import "@deckyard/mitoleaf/" as mitoleaf;
```

Signatures, constants and struct fields are in the
[cheatsheet](cheatsheet.md). This page covers input rules and errors.

## Input

`parse` reads the phm XML as downloaded: comments and wrapper elements around
the root `<Node>` are skipped. The element name matches in any case; the
attributes are read as `Id`, `id` or `ID` and `HG`, `hg` or `Hg`. A missing
`HG` is an empty motif. Motif tokens are split on spaces and stored verbatim;
they are interpreted, and a malformed one rejected, only when `haplotype`,
`distance` or `assign` folds them.

`parseMotifItem` trims and upper-cases a token before reading it. A back
mutation with nothing to revert is ignored, not rejected, because published
trees restate some.

Haplogroup names match exactly, case included.

## Errors

Errors are `Error{kind: "mitoleaf"}`, except an invalid motif gain, which
raises `kind: "lineage"` from `lineage.parseDifference`.

| Condition | `kind` | Thrown by |
| --- | --- | --- |
| the document holds no `<Node>` | `mitoleaf` | `parse` |
| a `<Node>` has no `Id` | `mitoleaf` | `parse` |
| a haplogroup name appears twice | `mitoleaf` | `parse` |
| unknown haplogroup | `mitoleaf` | `node`, `children`, `motif`, `path`, `haplotype`, `distance` |
| empty token | `mitoleaf` | `parseMotifItem` |
| a leading-minus token not of the form `-<position>.<insertion>` | `mitoleaf` | `parseMotifItem`, `haplotype`, `distance`, `assign` |
| a gain that is not valid rCRS nomenclature | `lineage` | `parseMotifItem`, `haplotype`, `distance`, `assign` |

`has` and `haplogroups` never throw. `assign` on a `Tree` with no nodes (which
`parse` never returns) fails with a runtime error, not a `mitoleaf` error.
