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
carries counts 1, so a position both carry with different bases counts 2. The
observed haplotype may carry one designation per position and insertion index;
`263G 263G` or `263G 263A` throws `kind: "mitoleaf"` from `distance` and
`assign` alike, since `lineage.mtDistance` would count the pair twice.

## Assignment

```jennifer
def m as mitoleaf.Match init mitoleaf.assign($t, $observed);
$m.haplogroup;
$m.distance;
```

The haplogroup with the smallest `distance`. One depth-first pass in document
order keeps one motif state and one running distance, adjusts both only for the
slots each node's motif touches, and undoes them on the way back up, so each
haplogroup is scored once in time proportional to its motif. Ties go to the first node reached: an ancestor beats an equally
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
| `parse` of the 549 KB phm XML, motifs pre-parsed | ~0.7 s |
| `assign`, one pass over all 6409 haplogroups | ~0.4 s |
| `has`, `path`, `node`, `children`, `motif`, per call from an importing program | ~5.5 ms |
| `haplotype`, `distance`, per call from an importing program | ~7.5 ms |

Parse the tree once and reuse the `Tree`. `parse` parses every motif token once
into `Tree.ops`, so `assign`, `haplotype` and `distance` never re-parse one; it
builds a flat node list and assembles the tree in one pass; threading the `Tree` through the recursion
deep-copied the node map at every node and took 2m14s on the same file.
`assign` avoids rebuilding each haplotype from the root, which would multiply
the work by the tree's depth, and walks with an explicit stack: as a recursion
that passed the `Tree` along it took 56 s in an importing program.

Every call into the deck that passes a `Tree` costs ~5.5 ms on this tree, before
any work: the interpreter rebuilds a value holding the deck's own structs
(`Tree`, `Node`) at each call across the module boundary. A loop calling
`distance` once per haplogroup therefore takes ~1 min; `assign` scores the
whole tree in one call.
