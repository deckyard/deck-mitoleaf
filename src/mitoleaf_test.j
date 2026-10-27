# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# White-box test overlay for mitoleaf.j. Run with `jennifer test mitoleaf_test.j`.
#
# Expected haplotypes are worked out by hand from testdata/tree.xml rather than
# captured from this implementation, so a regression in motif accumulation or
# back-mutation handling fails the suite instead of being blessed.
use testing;
use fs;

# testdata/tree.xml, inline: a consumer runs this overlay from its own tree,
# where testdata/ does not exist. testTreeMatchesTestdata keeps the two copies
# equal.
def const TREE as string init '<?xml version="1.0" encoding="UTF-8"?>
<!--
Synthetic mitoLEAF-shaped fixture for the test overlay and the example.
It is NOT the mitoLEAF phylogeny - four haplogroups, chosen to exercise the
parts that matter: an insertion, an IUPAC ambiguity, a deletion, and back
mutations in both spellings - base-first (A263) and leading-minus (-315.1).
Get the real tree from https://forensicgenomics.github.io/mitoLeaf/.
-->
<Node Id="mt-MRCA" HG="73G 146C 263G 315.1c 16519Y">
  <Node Id="L0" HG="A263 1048T">
    <Node Id="L0a" HG="T146 5231A 16278T">
      <Node Id="L0a1" HG="C16278 8281-"/>
    </Node>
    <Node Id="L0b" HG="16187T"/>
  </Node>
  <Node Id="L1" HG="3594T 16187T">
    <Node Id="L1a" HG="-315.1 8994A"/>
  </Node>
</Node>';

func tree() {
    return parse(TREE);
}

# --- parsing ----------------------------------------------------------------

func testParseFindsTheRootPastTheComments() {
    testing.assertEqual(tree().root, "mt-MRCA");
}

func testParseBuildsTheWholeTree() {
    def t as Tree init tree();
    testing.assertEqual(haplogroups($t), ["mt-MRCA", "L0", "L0a", "L0a1", "L0b", "L1", "L1a"]);
}

func testParentAndChildrenAreLinkedBothWays() {
    def t as Tree init tree();
    testing.assertEqual(children($t, "mt-MRCA"), ["L0", "L1"]);
    testing.assertEqual(children($t, "L0"), ["L0a", "L0b"]);
    testing.assertEqual(children($t, "L0a1"), []);
    testing.assertEqual(node($t, "L0a").parent, "L0");
    testing.assertEqual(node($t, "mt-MRCA").parent, "");
}

func testMotifIsKeptVerbatim() {
    testing.assertEqual(motif(tree(), "L0a"), ["T146", "5231A", "16278T"]);
}

func testPathRunsRootFirst() {
    testing.assertEqual(path(tree(), "L0a1"), ["mt-MRCA", "L0", "L0a", "L0a1"]);
}

func testHasAnswersForEveryDepthAndOnlyExactNames() {
    # tree.xml names mt-MRCA at the root and L0a1 as its deepest leaf; Z9 is
    # absent, and names are matched exactly, so "l0a" is not "L0a".
    def t as Tree init tree();
    testing.assertTrue(has($t, "mt-MRCA"));
    testing.assertTrue(has($t, "L0a1"));
    testing.assertTrue(not has($t, "Z9"));
    testing.assertTrue(not has($t, "l0a"));
    testing.assertTrue(not has($t, ""));
}

func missingHaplogroup() {
    node(tree(), "Z9");
}

func testUnknownHaplogroupThrows() {
    testing.assertThrows("missingHaplogroup", "mitoleaf");
}

func noNodeElement() {
    parse("<phylotree><clade Id=\"x\"/></phylotree>");
}

func testDocumentWithoutNodesThrows() {
    testing.assertThrows("noNodeElement", "mitoleaf");
}

# --- motif tokens -----------------------------------------------------------

func testGainParsesThroughLineage() {
    def m as MotifItem init parseMotifItem("315.1c");
    testing.assertTrue(not $m.reversion);
    testing.assertEqual($m.position, 315);
    testing.assertEqual($m.insertion, 1);
    testing.assertEqual($m.designation, "315.1C");
}

func testDeletionAndAmbiguityAreGains() {
    testing.assertEqual(parseMotifItem("8281-").designation, "8281d");
    testing.assertEqual(parseMotifItem("16519Y").designation, "16519Y");
}

func testBackMutationIsRecognised() {
    def m as MotifItem init parseMotifItem("A263");
    testing.assertTrue($m.reversion);
    testing.assertEqual($m.position, 263);
}

func testRedundantReferenceBaseIsNotABackMutation() {
    # "A263G" is a substitution written with the reference base; only "A263"
    # with nothing after the position is a reversion.
    def m as MotifItem init parseMotifItem("A263G");
    testing.assertTrue(not $m.reversion);
    testing.assertEqual($m.designation, "263G");
}

# --- accumulated haplotypes -------------------------------------------------

func testRootHaplotypeIsItsOwnMotif() {
    testing.assertEqual(
        lineage.mtDesignations(haplotype(tree(), "mt-MRCA")),
        ["73G", "146C", "263G", "315.1C", "16519Y"]);
}

func testBackMutationDropsAnInheritedDifference() {
    # L0 reverts 263, so 263G from the root is gone and 1048T is added.
    testing.assertEqual(
        lineage.mtDesignations(haplotype(tree(), "L0")),
        ["73G", "146C", "315.1C", "1048T", "16519Y"]);
}

func testAccumulationRunsDownTheWholeLineage() {
    # L0a reverts 146 as well, and adds 5231A and 16278T.
    testing.assertEqual(
        lineage.mtDesignations(haplotype(tree(), "L0a")),
        ["73G", "315.1C", "1048T", "5231A", "16278T", "16519Y"]);
}

func testABackMutationCanUndoAParentsGain() {
    # L0a1 reverts 16278, which L0a itself contributed, and deletes 8281.
    testing.assertEqual(
        lineage.mtDesignations(haplotype(tree(), "L0a1")),
        ["73G", "315.1C", "1048T", "5231A", "8281d", "16519Y"]);
}

func testSiblingsDoNotShareEachOthersMotifs() {
    testing.assertEqual(
        lineage.mtDesignations(haplotype(tree(), "L1")),
        ["73G", "146C", "263G", "315.1C", "3594T", "16187T", "16519Y"]);
}

func testHaplotypeCoversTheWholeGenome() {
    testing.assertEqual(haplotype(tree(), "L0a").ranges, [WHOLE_GENOME]);
}

# --- distance and assignment ------------------------------------------------

func observed(text as string) {
    return lineage.mtHaplotype("evidence", [WHOLE_GENOME], $text);
}

func testDistanceToItsOwnHaplogroupIsZero() {
    def h as lineage.MtHaplotype init observed("73G 315.1C 1048T 5231A 16278T 16519Y");
    testing.assertEqual(distance(tree(), "L0a", $h), 0);
}

func testDistanceCountsBothDirections() {
    # One designation short of L0a, and one it does not carry.
    def h as lineage.MtHaplotype init observed("73G 315.1C 1048T 5231A 9999A 16519Y");
    testing.assertEqual(distance(tree(), "L0a", $h), 2);
}

func testDifferentBaseAtOnePositionCountsTwice() {
    def h as lineage.MtHaplotype init observed("73G 315.1C 1048T 5231A 16278C 16519Y");
    testing.assertEqual(distance(tree(), "L0a", $h), 2);
}

func testAssignFindsTheExactHaplogroup() {
    def h as lineage.MtHaplotype init observed("73G 315.1C 1048T 5231A 8281d 16519Y");
    def m as Match init assign(tree(), $h);
    testing.assertEqual($m.haplogroup, "L0a1");
    testing.assertEqual($m.distance, 0);
}

func testAssignFindsTheNearestWhenNothingMatchesExactly() {
    # L0a plus one extra difference: still nearer L0a than anything else.
    def h as lineage.MtHaplotype init observed("73G 315.1C 1048T 5231A 16278T 16519Y 9999A");
    def m as Match init assign(tree(), $h);
    testing.assertEqual($m.haplogroup, "L0a");
    testing.assertEqual($m.distance, 1);
}

func testAssignAgreesWithDistance() {
    def h as lineage.MtHaplotype init observed("73G 146C 263G 315.1C 3594T 16187T 16519Y");
    def m as Match init assign(tree(), $h);
    testing.assertEqual($m.haplogroup, "L1");
    testing.assertEqual($m.distance, distance(tree(), "L1", $h));
}

# --- the leading-minus back mutation ----------------------------------------

func testInsertionLossIsRecognised() {
    def m as MotifItem init parseMotifItem("-5899.1");
    testing.assertTrue($m.reversion);
    testing.assertEqual($m.position, 5899);
    testing.assertEqual($m.insertion, 1);
}

func testInsertionLossDropsTheInsertedBase() {
    # L1a carries "-315.1", so the root's 315.1C is gone while 315 itself is
    # untouched - the two are separate slots.
    testing.assertEqual(
        lineage.mtDesignations(haplotype(tree(), "L1a")),
        ["73G", "146C", "263G", "3594T", "8994A", "16187T", "16519Y"]);
}

func testDeletionGainIsNotAnInsertionLoss() {
    # "8281-" gains a deletion; "-8281.1" would drop an insertion. Minus last
    # and minus first are opposite operations.
    testing.assertTrue(not parseMotifItem("8281-").reversion);
    testing.assertTrue(parseMotifItem("-8281.1").reversion);
}

func malformedInsertionLoss() {
    parseMotifItem("-59x9.1");
}

func testMalformedInsertionLossThrows() {
    testing.assertThrows("malformedInsertionLoss", "mitoleaf");
}

# --- the inline fixture -------------------------------------------------------

# Run from the deck root, the inline copy must equal the file the example and the
# docs use; elsewhere (an installed deck) there is no file to compare against.
func testTreeMatchesTestdata() {
    if (not fs.isFile("testdata/tree.xml")) {
        return;
    }
    testing.assertEqual(TREE + "\n", fs.readString("testdata/tree.xml"));
}
