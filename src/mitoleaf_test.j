# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

# White-box test overlay for mitoleaf.j. Run with `jennifer test mitoleaf_test.j`.
#
# Expected haplotypes are worked out by hand from testdata/tree.xml rather than
# captured from this implementation, so a regression in motif accumulation or
# back-mutation handling fails the suite instead of being blessed.
use testing;
use fs;
use strings;
use maps;

# testdata/tree.xml, inline: a consumer runs this overlay from its own tree,
# where testdata/ does not exist. testTreeMatchesTestdata keeps the two copies
# equal.
def const TREE as string init '<?xml version="1.0" encoding="UTF-8"?>
<!--
Synthetic mitoLEAF-shaped fixture for the test overlay and the example.
It is NOT the mitoLEAF phylogeny - seven haplogroups, chosen to exercise the
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

# --- review regressions -----------------------------------------------------

# haplotypeIn parses a one-off tree and returns one haplogroup's designations.
func haplotypeIn(doc as string, hg as string) {
    return lineage.mtDesignations(haplotype(parse($doc), $hg));
}

func testBackMutationWithNothingToRevertIsIgnored() {
    # A73 and -315.1 name slots nothing on the lineage set, so A keeps exactly
    # its parent's 263G; assign must not abort on them either.
    def doc as string init '<Node Id="R" HG="263G"><Node Id="A" HG="A73 -315.1"/></Node>';
    testing.assertEqual(haplotypeIn($doc, "A"), ["263G"]);
    def m as Match init assign(parse($doc), observed("263G"));
    testing.assertEqual($m.haplogroup, "R");
    testing.assertEqual($m.distance, 0);
}

func testPartialInsertionLossDropsTheHighestIndex() {
    # R inserts two bases at 5899; one loss shortens the run from its end
    # (5899.2C goes, 5899.1C stays), a second loss empties it.
    def doc as string init '<Node Id="R" HG="5899.1C 5899.2C"><Node Id="A" HG="-5899.1">' +
        '<Node Id="B" HG="-5899.1"/></Node></Node>';
    testing.assertEqual(haplotypeIn($doc, "A"), ["5899.1C"]);
    testing.assertEqual(haplotypeIn($doc, "B"), []);
    # assign follows the same rule: an observed 5899.1C is A exactly
    def m as Match init assign(parse($doc), observed("5899.1C"));
    testing.assertEqual($m.haplogroup, "A");
    testing.assertEqual($m.distance, 0);
}

func twoRoots() {
    parse('<x><Node Id="R" HG="263G"/><Node Id="S" HG="9999A"/></x>');
}

func testTwoTopLevelNodesThrow() {
    testing.assertThrows("twoRoots", "mitoleaf");
    def msg as string init "";
    try {
        twoRoots();
    } catch (e) {
        $msg = $e.message;
    }
    testing.assertTrue(strings.contains($msg, "'R'") and strings.contains($msg, "'S'"));
}

func malformedXml() {
    parse('<Node Id="R"');
}

func testMalformedXmlThrowsMitoleaf() {
    testing.assertThrows("malformedXml", "mitoleaf");
}

# rawObserved builds a haplotype without lineage.mtHaplotype, which in newer
# forensicgenetics releases rejects these inputs itself.
func rawObserved(designations as list of string) {
    def ds as list of lineage.Difference init [];
    for (def d in $designations) {
        $ds[] = lineage.parseDifference($d);
    }
    return lineage.MtHaplotype{id: "evidence", ranges: [WHOLE_GENOME], differences: $ds};
}

func assignConflict() {
    assign(tree(), rawObserved(["263G", "263A"]));
}

func distanceConflict() {
    distance(tree(), "L0", rawObserved(["263G", "263A"]));
}

func assignDuplicate() {
    assign(tree(), rawObserved(["9999A", "9999A"]));
}

func distanceDuplicate() {
    distance(tree(), "L0", rawObserved(["9999A", "9999A"]));
}

func testTwoDesignationsAtOnePositionThrowInAssignAndDistance() {
    testing.assertThrows("assignConflict", "mitoleaf");
    testing.assertThrows("distanceConflict", "mitoleaf");
    testing.assertThrows("assignDuplicate", "mitoleaf");
    testing.assertThrows("distanceDuplicate", "mitoleaf");
}

func testDistanceIsLineageMtDistanceForEveryHaplogroup() {
    def t as Tree init tree();
    def h as lineage.MtHaplotype init observed("73G 146C 315.1C 1048T 5231A 9999A 16519Y");
    for (def hg in haplogroups($t)) {
        testing.assertEqual(distance($t, $hg, $h), lineage.mtDistance(haplotype($t, $hg), $h));
    }
}

func testAssignAgreesWithDistanceEverywhere() {
    # Brute force: the minimum of `distance` over the tree in document order,
    # first minimum kept, for each haplogroup's own haplotype plus one extra
    # difference.
    def t as Tree init tree();
    for (def target in haplogroups($t)) {
        def text as string init strings.join(lineage.mtDesignations(haplotype($t, $target)), " ");
        def h as lineage.MtHaplotype init observed($text + " 9999A");
        def want as Match init Match{haplogroup: "", distance: -1};
        for (def hg in haplogroups($t)) {
            def d as int init distance($t, $hg, $h);
            if ($want.distance < 0 or $d < $want.distance) {
                $want = Match{haplogroup: $hg, distance: $d};
            }
        }
        testing.assertEqual(assign($t, $h), $want);
        testing.assertEqual($want, Match{haplogroup: $target, distance: 1});
    }
}

func reversionAtZero() {
    parseMotifItem("A0");
}

func reversionPastTheGenome() {
    parseMotifItem("A16570");
}

func reversionOverflowing() {
    parseMotifItem("A99999999999999999999999");
}

func insertionLossAtIndexZero() {
    parseMotifItem("-263.0");
}

func insertionLossAtPositionZero() {
    parseMotifItem("-0.1");
}

func insertionLossPastTheGenome() {
    parseMotifItem("-16570.1");
}

func insertionLossOverlongIndex() {
    parseMotifItem("-5899.99999999999999999999");
}

func testBackMutationOutsideTheReferenceThrows() {
    testing.assertThrows("reversionAtZero", "mitoleaf");
    testing.assertThrows("reversionPastTheGenome", "mitoleaf");
    testing.assertThrows("reversionOverflowing", "mitoleaf");
    testing.assertThrows("insertionLossAtIndexZero", "mitoleaf");
    testing.assertThrows("insertionLossAtPositionZero", "mitoleaf");
    testing.assertThrows("insertionLossPastTheGenome", "mitoleaf");
    testing.assertThrows("insertionLossOverlongIndex", "mitoleaf");
    # the ends of the range are accepted: 16569 is the last rCRS position
    testing.assertEqual(parseMotifItem("A16569").position, 16569);
    testing.assertEqual(parseMotifItem("-16569.1").position, 16569);
    testing.assertEqual(parseMotifItem("A1").position, 1);
}

# --- the parsed-motif cache --------------------------------------------------

func testParseFillsOpsForEveryHaplogroup() {
    def t as Tree init tree();
    testing.assertEqual(maps.keys($t.ops), haplogroups($t));
    # L0 is "A263 1048T": a reversion, then a gain
    testing.assertEqual(len(opsOf($t, "L0")), 2);
    # 16187T: kind G, position padded to 016187, insertion 0000, designation
    testing.assertEqual(opsOf($t, "L0b"), ["G016187.000016187T"]);
}

func testFromNodesRebuildsTheParsedTree() {
    def t as Tree init tree();
    testing.assertEqual(fromNodes($t.root, $t.nodes), $t);
}

func testTreeWithoutOpsGivesTheSameResults() {
    # an empty `ops` makes every fold parse the motif itself
    def t as Tree init tree();
    def bare as Tree init Tree{root: $t.root, nodes: $t.nodes, ops: {}};
    def h as lineage.MtHaplotype init observed("73G 315.1C 1048T 5231A 9999A 16519Y");
    for (def hg in haplogroups($t)) {
        testing.assertEqual(haplotype($bare, $hg), haplotype($t, $hg));
        testing.assertEqual(distance($bare, $hg, $h), distance($t, $hg, $h));
    }
    testing.assertEqual(assign($bare, $h), assign($t, $h));
}

func testFromNodesPicksUpAnEditedMotif() {
    # L0b's motif 16187T replaced by 9999A: its haplotype is L0's plus 9999A
    def t as Tree init tree();
    def nodes as map of string to Node init $t.nodes;
    $nodes["L0b"].motif = ["9999A"];
    testing.assertEqual(
        lineage.mtDesignations(haplotype(fromNodes($t.root, $nodes), "L0b")),
        ["73G", "146C", "315.1C", "1048T", "9999A", "16519Y"]);
}

func badTokenTree() {
    return parse('<Node Id="R" HG="263G"><Node Id="A" HG="9999A"/><Node Id="B" HG="263Q"/></Node>');
}

func badTokenHaplotype() {
    haplotype(badTokenTree(), "B");
}

func badTokenAssign() {
    assign(badTokenTree(), observed("263G"));
}

func testMalformedTokenParsesButThrowsWhenFolded() {
    # parse keeps the tree; the token's lineage error comes from the fold
    testing.assertEqual(lineage.mtDesignations(haplotype(badTokenTree(), "A")), ["263G", "9999A"]);
    testing.assertThrows("badTokenHaplotype", "lineage");
    testing.assertThrows("badTokenAssign", "lineage");
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
