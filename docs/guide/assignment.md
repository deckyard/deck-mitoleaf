# Haplotypes and assignment

## A haplogroup's haplotype

```jennifer
def h as lineage.MtHaplotype init mitoleaf.haplotype($t, "L0a1");
lineage.mtDesignations($h);
```

The motifs from the root down, with back mutations applied (see
[Motifs and back mutations](motifs.md)). The haplotype's one range is
`mitoleaf.WHOLE_GENOME`, `"1-16569"`, because a motif may name any position.

An observed haplotype usually covers only the control region (HVR-I, or the
whole CR), so the ranges disagree:

```jennifer
lineage.mtRangesAgree(mitoleaf.haplotype($t, "L0a1"), $observed);   # false
```

That is correct. `lineage` keeps "do the differences match" and "were these
sequenced alike" as separate questions, and `distance` and `assign` answer only
the first.

## Distance

```jennifer
mitoleaf.distance($t, "L0a1", $observed);
```

`lineage.mtDistance` over the two difference lists: a designation only one side
carries counts 1, so a position both carry with different bases counts 2.

## Assignment

```jennifer
def m as mitoleaf.Match init mitoleaf.assign($t, $observed);
$m.haplogroup;
$m.distance;
```

The haplogroup with the smallest `distance`. One depth-first pass in document
order carries the accumulated motif state down the tree, so each haplogroup is
scored once. Ties go to the first node reached: an ancestor beats an equally
distant descendant, and otherwise the earlier node in the file wins.

`assign` is a nearest-neighbour report, not a haplogroup call; the reasons are
in [Limitations](../project/limitations.md). Check the result against the tree:

```jennifer
for (def hg in mitoleaf.path($t, $m.haplogroup)) {
    io.printf("%s %v\n", $hg, mitoleaf.motif($t, $hg));
}
```

## Performance

On the published 6409-haplogroup tree, one desktop core:

| Operation | Time |
| --- | --- |
| `parse` of the 549 KB phm XML | ~0.5 s |
| `assign`, one pass over all 6409 haplogroups | ~2 s |

Parse the tree once and reuse the `Tree`. `parse` builds a flat node list and
assembles the tree in one pass; threading the `Tree` through the recursion
deep-copied the node map at every node and took 2m14s on the same file.
`assign` avoids rebuilding each haplotype from the root, which would multiply
the work by the tree's depth.
