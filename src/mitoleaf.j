# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

/**
 * Reader for the mitoLEAF human mitochondrial phylogeny.
 *
 * mitoLEAF publishes the haplogroup tree as nested `<Node>` elements (the "phm"
 * XML). Each node carries a haplogroup name and the motif it adds to its
 * parent:
 *
 * ```xml
 * <Node Id="mt-MRCA" HG="73G 146C 263G 315.1c 16519Y">
 *   <Node Id="L0" HG="A263 1048T 3516A">
 *     <Node Id="L0a" HG="T146 5231A C16278"/>
 *   </Node>
 * </Node>
 * ```
 *
 * A motif token is a gain in rCRS nomenclature (`73G`, `315.1c`, `8281-`,
 * `16519Y`) or a back mutation in one of mitoLEAF's two spellings: `A263`
 * (base first, nothing after the position) reverts a substitution, `-5899.1`
 * (leading minus) drops an inserted base. `lineage.parseDifference` rejects
 * both back-mutation spellings, so this module reads them itself and hands
 * only gains to `lineage`.
 *
 * `haplotype` accumulates the motifs from the root with back mutations applied
 * and returns a `lineage.MtHaplotype`, comparable with `lineage.mtMatches` and
 * `lineage.mtDistance`.
 *
 * No tree is bundled: download the phm XML from mitoLEAF and pass its text to
 * `parse`. Haplogroup definitions change between mitoLEAF releases, so record
 * the release with every result.
 * @module mitoleaf
 * @see https://forensicgenomics.github.io/mitoLeaf/ - the mitoLEAF project
 * @see https://doi.org/10.1093/nargab/lqaf079 - the mitoLEAF paper
 * @example
 * import "@deckyard/mitoleaf/" as mitoleaf;
 * import "@deckyard/forensicgenetics/lineage.j" as lineage;
 *
 * def t as mitoleaf.Tree init mitoleaf.parse(fs.readString("mitoLEAF_phm.xml"));
 * def h as lineage.MtHaplotype init mitoleaf.haplotype($t, "L0a");
 * def m as mitoleaf.Match init mitoleaf.assign($t, $observed);
 */
use xml;
use strings;
use maps;
use lists;
use convert;

import "@deckyard/forensicgenetics/lineage.j" as lineage;

/**
 * The range every phylogeny haplotype declares: the whole rCRS, because a
 * motif may name any position.
 */
export def const WHOLE_GENOME as string init "1-16569";

/**
 * One haplogroup in the tree.
 * @field id {string} haplogroup name, e.g. "L0a1b2a"
 * @field parent {string} parent haplogroup, "" at the root
 * @field motif {list of string} tokens this node adds to its parent, verbatim
 * @field children {list of string} child haplogroups, in document order
 */
export def struct Node {
    id as string,
    parent as string,
    motif as list of string,
    children as list of string
};

/**
 * A parsed phylogeny.
 * @field root {string} root haplogroup, "mt-MRCA" in the published tree
 * @field nodes {map of string to Node} haplogroup name to node, in document order
 */
export def struct Tree {
    root as string,
    nodes as map of string to Node
};

/**
 * One motif token, interpreted.
 * @field reversion {bool} true for a back mutation (`A263`, `-5899.1`)
 * @field position {int} rCRS position the token names
 * @field insertion {int} n for the n-th base inserted after `position`, else 0
 * @field designation {string} normalised rCRS form of a gain; for a back
 *   mutation the token, trimmed and upper-cased
 */
export def struct MotifItem {
    reversion as bool,
    position as int,
    insertion as int,
    designation as string
};

/**
 * The haplogroup nearest an observed haplotype.
 * @field haplogroup {string} haplogroup name
 * @field distance {int} `lineage.mtDistance` from the observed haplotype
 */
export def struct Match {
    haplogroup as string,
    distance as int
};

# --- errors -----------------------------------------------------------------

func fail(msg as string) {
    throw Error{kind: "mitoleaf", message: $msg, file: "", line: 0, col: 0};
}

# --- motif tokens -----------------------------------------------------------

func isDigit(ch as string) {
    return strings.indexOf("0123456789", $ch) >= 0;
}

# The IUPAC nucleotide alphabet, so an ambiguity code standing for a
# heteroplasmic position parses like any other base.
func isBase(ch as string) {
    return strings.indexOf("ACGTURYSWKMBDHVN", $ch) >= 0;
}

func allDigits(ch as list of string, from as int) {
    if ($from >= len($ch)) {
        return false;
    }
    def i as int init $from;
    while ($i < len($ch)) {
        if (not isDigit($ch[$i])) {
            return false;
        }
        $i = $i + 1;
    }
    return true;
}

# A reverting substitution is written base-first with nothing after the
# position - "A263", never "A263G" (which is a redundantly spelled
# substitution).
func isReversion(ch as list of string) {
    if (len($ch) < 2 or not isBase($ch[0])) {
        return false;
    }
    return allDigits($ch, 1);
}

# The loss of an inserted base is written with a leading minus: "-5899.1"
# drops whatever an ancestor inserted at 5899.1.
func insertionLoss(token as string, ch as list of string) {
    def parts as list of string init strings.split(strings.join($ch[1..], ""), ".");
    if (len($parts) != 2) {
        fail("motif token '" + $token + "': expected -<position>.<insertion>");
    }
    def pos as list of string init strings.chars($parts[0]);
    def ins as list of string init strings.chars($parts[1]);
    if (not allDigits($pos, 0) or not allDigits($ins, 0)) {
        fail("motif token '" + $token + "': expected -<position>.<insertion>");
    }
    return MotifItem{
        reversion: true,
        position: convert.toInt($parts[0]),
        insertion: convert.toInt($parts[1]),
        designation: $token
    };
}

/**
 * Interpret one motif token. Case and surrounding whitespace are ignored.
 *
 * Gains go through `lineage.parseDifference`, so every form it accepts
 * (`263G`, `A263G`, `315.1C`, `8281-`, `16093Y`) works here and `designation`
 * comes back normalised. `A263` (nothing after the position) reverts 263 to
 * the reference; `-5899.1` drops the base inserted at 5899.1. `8281-` (minus
 * last) is the opposite of `-8281.1`: it gains a deletion.
 * @param token {string} motif token as published
 * @return {MotifItem} the interpreted token
 * @throws {Error} kind "mitoleaf" on an empty token or a leading-minus token
 *   not of the form `-<position>.<insertion>`; kind "lineage" on a gain that
 *   is not valid rCRS nomenclature
 */
export func parseMotifItem(token as string) {
    def t as string init strings.upper(strings.trim($token));
    if (len($t) == 0) {
        fail("empty motif token");
    }
    def ch as list of string init strings.chars($t);
    if ($ch[0] == "-") {
        return insertionLoss($t, $ch);
    }
    if (isReversion($ch)) {
        return MotifItem{
            reversion: true,
            position: convert.toInt(strings.join($ch[1..], "")),
            insertion: 0,
            designation: $t
        };
    }
    def d as lineage.Difference init lineage.parseDifference($t);
    return MotifItem{
        reversion: false,
        position: $d.position,
        insertion: $d.insertion,
        designation: lineage.formatDifference($d)
    };
}

# motifKey identifies the reference position a token acts on, so a gain and a
# later reversion at the same place meet in the same slot.
func motifKey(position as int, insertion as int) {
    def p as string init convert.toString($position);
    def i as string init convert.toString($insertion);
    return $p + "." + $i;
}

# splitMotif splits an HG attribute on whitespace, dropping the empty pieces a
# run of spaces produces.
func splitMotif(text as string) {
    def out as list of string init [];
    for (def part in strings.split(strings.trim($text), " ")) {
        if (len(strings.trim($part)) > 0) {
            $out[] = strings.trim($part);
        }
    }
    return $out;
}

# applyMotif folds one node's motif into the accumulated state: a gain sets its
# slot, a back mutation clears it. A reversion with nothing to revert is
# ignored rather than rejected - a published tree may restate one.
func applyMotif(state as map of string to string, motif as list of string) {
    def out as map of string to string init $state;
    for (def token in $motif) {
        def m as MotifItem init parseMotifItem($token);
        def key as string init motifKey($m.position, $m.insertion);
        if ($m.reversion) {
            $out = maps.delete($out, $key);
        } else {
            $out[$key] = $m.designation;
        }
    }
    return $out;
}

# --- parsing ----------------------------------------------------------------

# attrValue reads the first attribute present from a list of spellings, so a
# document using `id` / `hg` reads the same as one using `Id` / `HG`.
func attrValue(e as xml.Value, names as list of string) {
    for (def n in $names) {
        if (xml.hasAttr($e, $n)) {
            return xml.attr($e, $n);
        }
    }
    return "";
}

# Flat is one node as read off the document, before the tree is assembled.
# Collecting these first keeps parsing linear: threading the whole `Tree`
# through the recursion would deep-copy every node's map at every node.
def struct Flat {
    id as string,
    parent as string,
    motif as list of string
};

# Pending is one element still to flatten, under the haplogroup that encloses it.
def struct Pending {
    e as xml.Value,
    parent as string
};

# flatten walks the document in preorder with an explicit stack: recursing
# and splicing each subtree's list into its parent's copies every node once
# per level above it.
func flatten(root as xml.Value) {
    def out as list of Flat init [];
    def stack as list of Pending init [Pending{e: $root, parent: ""}];
    def size as int init 1;
    while ($size > 0) {
        $size = $size - 1;
        def top as Pending init $stack[$size];
        if (xml.typeOf($top.e) != "element") {
            continue;
        }
        def under as string init $top.parent;
        if (strings.lower(xml.tag($top.e)) == "node") {
            def id as string init strings.trim(attrValue($top.e, ["Id", "id", "ID"]));
            if (len($id) == 0) {
                fail("a <Node> element carries no Id attribute");
            }
            $out[] = Flat{
                id: $id,
                parent: $under,
                motif: splitMotif(attrValue($top.e, ["HG", "hg", "Hg"]))
            };
            $under = $id;
        }
        # pushed in reverse, so children come off in document order
        def kids as list of xml.Value init xml.children($top.e);
        for (def i as int init len($kids) - 1; $i >= 0; $i = $i - 1) {
            def p as Pending init Pending{e: $kids[$i], parent: $under};
            if ($size < len($stack)) {
                $stack[$size] = $p;
            } else {
                $stack[] = $p;
            }
            $size = $size + 1;
        }
    }
    return $out;
}

/**
 * Parse a mitoLEAF phm XML document.
 *
 * Comments and wrapper elements around the root `<Node>` are skipped, so the
 * published file parses as downloaded. The element name matches in any case;
 * the attributes are read as `Id`, `id` or `ID` and `HG`, `hg` or `Hg`. Motif
 * tokens are stored verbatim and not validated here: a malformed one throws
 * later, from `haplotype`, `distance` or `assign`.
 * @param doc {string} XML document text
 * @return {Tree} the parsed phylogeny
 * @throws {Error} kind "mitoleaf" if the document holds no `<Node>`, a node
 *   has no `Id`, or a haplogroup name appears twice
 */
export func parse(doc as string) {
    def flat as list of Flat init flatten(xml.decode($doc));
    if (len($flat) == 0) {
        fail("no <Node> element found in the document");
    }
    def t as Tree init Tree{root: $flat[0].id, nodes: {}};
    for (def f in $flat) {
        if (maps.has($t.nodes, $f.id)) {
            fail("duplicate haplogroup '" + $f.id + "' in the tree");
        }
        $t.nodes[$f.id] = Node{id: $f.id, parent: $f.parent, motif: $f.motif, children: []};
    }
    # Linked in place on the one binding: routing this through a helper that
    # takes and returns a Tree would deep-copy all 6000-odd nodes per link.
    for (def f in $flat) {
        if (len($f.parent) > 0) {
            $t.nodes[$f.parent].children[] = $f.id;
        }
    }
    return $t;
}

# --- navigation -------------------------------------------------------------

func requireNode(t as Tree, hg as string) {
    if (not maps.has($t.nodes, $hg)) {
        fail("no such haplogroup in the tree: '" + $hg + "'");
    }
}

/**
 * Every haplogroup in the tree, in document order.
 * @param t {Tree} parsed tree
 * @return {list of string} haplogroup names
 */
export func haplogroups(t as Tree) {
    return maps.keys($t.nodes);
}

/**
 * Whether the tree names a haplogroup. Names match exactly, case included.
 * @param t {Tree} parsed tree
 * @param hg {string} haplogroup name
 * @return {bool} true if present
 */
export func has(t as Tree, hg as string) {
    return maps.has($t.nodes, $hg);
}

/**
 * One haplogroup's node.
 * @param t {Tree} parsed tree
 * @param hg {string} haplogroup name
 * @return {Node} the node
 * @throws {Error} kind "mitoleaf" if the tree has no such haplogroup
 */
export func node(t as Tree, hg as string) {
    requireNode($t, $hg);
    return $t.nodes[$hg];
}

/**
 * A haplogroup's direct children.
 * @param t {Tree} parsed tree
 * @param hg {string} haplogroup name
 * @return {list of string} child haplogroup names, in document order
 * @throws {Error} kind "mitoleaf" if the tree has no such haplogroup
 */
export func children(t as Tree, hg as string) {
    requireNode($t, $hg);
    return $t.nodes[$hg].children;
}

/**
 * The motif a haplogroup adds to its parent, verbatim as published.
 * @param t {Tree} parsed tree
 * @param hg {string} haplogroup name
 * @return {list of string} motif tokens
 * @throws {Error} kind "mitoleaf" if the tree has no such haplogroup
 */
export func motif(t as Tree, hg as string) {
    requireNode($t, $hg);
    return $t.nodes[$hg].motif;
}

/**
 * The lineage from the root down to a haplogroup, both ends included.
 * @param t {Tree} parsed tree
 * @param hg {string} haplogroup name
 * @return {list of string} haplogroup names, root first
 * @throws {Error} kind "mitoleaf" if the tree has no such haplogroup
 */
export func path(t as Tree, hg as string) {
    requireNode($t, $hg);
    def up as list of string init [];
    def cur as string init $hg;
    while (len($cur) > 0) {
        $up[] = $cur;
        $cur = $t.nodes[$cur].parent;
    }
    return lists.reverse($up);
}

# --- haplotypes -------------------------------------------------------------

# stateOf accumulates the motifs along a lineage into position -> designation.
func stateOf(t as Tree, hg as string) {
    def out as map of string to string init {};
    for (def id in path($t, $hg)) {
        $out = applyMotif($out, $t.nodes[$id].motif);
    }
    return $out;
}

func stateHaplotype(id as string, state as map of string to string) {
    def parts as list of string init [];
    for (def key in $state) {
        $parts[] = $state[$key];
    }
    return lineage.mtHaplotype($id, [WHOLE_GENOME], strings.join($parts, " "));
}

/**
 * A haplogroup's full haplotype: its motifs accumulated from the root, with
 * back mutations applied. A back mutation with nothing to revert is ignored.
 *
 * The haplotype's id is `hg` and its one range is `WHOLE_GENOME`. An observed
 * haplotype usually covers only the control region, so
 * `lineage.mtRangesAgree` returns false for the pair; `distance` compares
 * designations regardless.
 * @param t {Tree} parsed tree
 * @param hg {string} haplogroup name
 * @return {lineage.MtHaplotype} the accumulated haplotype
 * @throws {Error} kind "mitoleaf" if the tree has no such haplogroup or a
 *   motif on the lineage is malformed; kind "lineage" on an invalid gain
 */
export func haplotype(t as Tree, hg as string) {
    return stateHaplotype($hg, stateOf($t, $hg));
}

/**
 * `lineage.mtDistance` between a haplogroup's haplotype and an observed one:
 * a designation only one side carries counts 1, so a position both carry with
 * different bases counts 2. Ranges are not compared.
 * @param t {Tree} parsed tree
 * @param hg {string} haplogroup name
 * @param h {lineage.MtHaplotype} observed haplotype
 * @return {int} the distance
 * @throws {Error} kind "mitoleaf" if the tree has no such haplogroup or a
 *   motif on the lineage is malformed; kind "lineage" on an invalid gain
 */
export func distance(t as Tree, hg as string, h as lineage.MtHaplotype) {
    return lineage.mtDistance(haplotype($t, $hg), $h);
}

# --- assignment -------------------------------------------------------------

# observedState keys an observed haplotype the same way a motif state is keyed,
# so the two compare without rebuilding an MtHaplotype per candidate.
func observedState(h as lineage.MtHaplotype) {
    def out as map of string to string init {};
    for (def d in $h.differences) {
        $out[motifKey($d.position, $d.insertion)] = lineage.formatDifference($d);
    }
    return $out;
}

# stateDistance is `lineage.mtDistance` over two keyed states: one point for a
# position only one side carries, two for a position both carry differently.
func stateDistance(a as map of string to string, b as map of string to string) {
    def total as int init 0;
    for (def key in $a) {
        if (not maps.has($b, $key)) {
            $total = $total + 1;
        } elseif ($a[$key] != $b[$key]) {
            $total = $total + 2;
        }
    }
    for (def key in $b) {
        if (not maps.has($a, $key)) {
            $total = $total + 1;
        }
    }
    return $total;
}

# Frame is one pending node of the search in `assign`, with the state its
# parent reached.
def struct Frame {
    id as string,
    state as map of string to string
};

/**
 * The haplogroup whose haplotype is nearest an observed one, by the same
 * measure as `distance`.
 *
 * One depth-first pass in document order scores every haplogroup once (about
 * 1.5 s on the 6409-haplogroup tree). Ties go to the first node reached, so an
 * ancestor beats an equally distant descendant.
 *
 * This is a nearest-neighbour report, not a validated haplogroup call: it
 * ignores which positions were sequenced, weights every difference alike and
 * applies none of the phylogenetic rules of a classifier such as EMPOP or
 * HaploGrep. Check the result by hand with `path` and `haplotype`.
 * @param t {Tree} tree from `parse`; a tree with no nodes is a runtime error
 * @param h {lineage.MtHaplotype} observed haplotype
 * @return {Match} the nearest haplogroup and its distance
 * @throws {Error} kind "mitoleaf" or "lineage" on a malformed motif token
 */
export func assign(t as Tree, h as lineage.MtHaplotype) {
    # An explicit preorder stack, not recursion: a recursive call that passes
    # the Tree deep-copies all of it per node once the deck is imported.
    def observed as map of string to string init observedState($h);
    def best as Match init Match{haplogroup: "", distance: -1};
    def empty as map of string to string init {};
    # `size` is the live depth: slots above it are reused, never sliced off,
    # since a slice copies every pending state
    def stack as list of Frame init [Frame{id: $t.root, state: $empty}];
    def size as int init 1;
    while ($size > 0) {
        $size = $size - 1;
        def top as Frame init $stack[$size];
        def here as map of string to string init applyMotif($top.state, $t.nodes[$top.id].motif);
        def d as int init stateDistance($here, $observed);
        if ($best.distance < 0 or $d < $best.distance) {
            $best = Match{haplogroup: $top.id, distance: $d};
        }
        # pushed in reverse, so children are visited in document order
        def kids as list of string init $t.nodes[$top.id].children;
        for (def i as int init len($kids) - 1; $i >= 0; $i = $i - 1) {
            def f as Frame init Frame{id: $kids[$i], state: $here};
            if ($size < len($stack)) {
                $stack[$size] = $f;
            } else {
                $stack[] = $f;
            }
            $size = $size + 1;
        }
    }
    return $best;
}
