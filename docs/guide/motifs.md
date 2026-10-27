# Motifs and back mutations

mitoLEAF publishes the tree as nested `<Node>` elements. Each carries a
haplogroup name and the motif it adds to its parent:

```xml
<Node Id="mt-MRCA" HG="73G 146C 263G 315.1c 16519Y">
  <Node Id="L0" HG="A263 1048T">
    <Node Id="L0a" HG="T146 5231A 16278T"/>
  </Node>
</Node>
```

## Token forms

The published tree (mitoLEAF 1.6.1: 6409 nodes, 18114 motif tokens) uses
five token forms: three kinds of gain and two of back mutation. Gains are rCRS
nomenclature, passed to `lineage.parseDifference`, which also accepts the
redundant `A263G` form:

| Token | Meaning |
| --- | --- |
| `263G` | substitution at 263 to G |
| `A263G` | the same, with the redundant reference base |
| `315.1c` | first base inserted after 315 (case is ignored) |
| `8281-` | deletion of 8281 |
| `16519Y` | IUPAC ambiguity code, a heteroplasmic position |

Two are back mutations, which drop what an ancestor contributed. Both
spellings are mitoLEAF's own and `lineage.parseDifference` rejects them, so the
deck reads them itself:

| Token | Meaning |
| --- | --- |
| `A263` | base first, nothing after the position: 263 reverts to the reference |
| `-5899.1` | leading minus: the base inserted at 5899.1 is gone |

`8281-` and `-8281.1` are opposite operations: minus last gains a deletion,
minus first drops an insertion. `A263` is a reversion; `A263G` is a
substitution. A back mutation with nothing to revert is ignored, because
published trees restate some.

```jennifer
mitoleaf.parseMotifItem("315.1c");   # reversion false, position 315, insertion 1
mitoleaf.parseMotifItem("A263");     # reversion true,  position 263, insertion 0
mitoleaf.parseMotifItem("-5899.1");  # reversion true,  position 5899, insertion 1
```

## Accumulation

A haplogroup's haplotype is not the concatenation of its ancestors' motifs:
positions are lost as well as gained. `haplotype` folds them from the root
down. Each token acts on one slot, its position plus insertion index, so
`-315.1` removes `315.1C` and leaves a change at 315 itself alone:

```
mt-MRCA   73G 146C 263G 315.1C 16519Y
L0        73G 146C      315.1C 16519Y 1048T        <- A263 dropped 263G
L0a       73G           315.1C 16519Y 1048T ...    <- T146 dropped 146C
```
