# Limitations

**`assign` is not a haplogroup classifier.** It reports the haplogroup at the
smallest symmetric difference. It:

- ignores which positions were sequenced, so an unobserved position and an
  absent one count the same
- weights every difference alike, where a classifier weights by how
  phylogenetically informative a position is and discounts unstable ones
- applies none of the phylogenetic reasoning, alignment rules included, of
  EMPOP's or HaploGrep's classifiers

Treat the result as a starting point and check it with `path` and `haplotype`.

**Range mismatch is unhandled by design.** Phylogeny haplotypes span the whole
genome; observations usually cover the control region only. `assign` compares
designations regardless, so a haplogroup defined by coding-region variants
cannot be told apart from its neighbours by a CR-only observation, and the
distance does not show that this happened. Check `lineage.mtRangesAgree`.

**No phylogenetic alignment.** EMPOP shifts indels to the position the mtDNA
phylogeny prefers. Neither this deck nor `lineage.align` does, so a complex
indel may be designated differently from how mitoLEAF or EMPOP write it, and
then it will not match.

**The tree moves.** mitoLEAF adds haplogroups, removes variants and respells
back mutations between releases. An assignment means something only together
with the release that produced it. The deck does not record the release; the
phm file states it in a leading comment.

**No tree is bundled.** `testdata/tree.xml` is a seven-haplogroup synthetic
fixture for the tests and the example, not the phylogeny.

**Whole-tree operations take seconds.** `parse` takes ~0.25 s and `assign` ~1.5 s
on the published tree ([figures](../guide/assignment.md#performance)): fine
interactively, too slow inside a loop over thousands of samples without
restructuring.

**Not a validated forensic tool.** Neither this deck nor `@deckyard/forensicgenetics`
has been through the validation any jurisdiction requires of software used in
casework.
