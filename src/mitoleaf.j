# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.26.0

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
 * A parsed phylogeny. Build one with `parse` or `fromNodes`; a hand-written
 * `Tree{...}` literal must set `ops` too (`{}` is valid but slower).
 * @field root {string} root haplogroup, "mt-MRCA" in the published tree
 * @field nodes {map of string to Node} haplogroup name to node, in document order
 * @field ops {map of string to string} haplogroup name to its motif
 *   pre-parsed, in an internal encoding; derived from `nodes`, so rebuild it
 *   with `fromNodes` after changing a motif. A haplogroup missing here has its
 *   motif parsed on every use
 */
export def struct Tree {
    root as string,
    nodes as map of string to Node,
    ops as map of string to string
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

# The last rCRS position: a back mutation names a reference position, so it
# must lie on the reference.
def const RCRS_LENGTH as int init 16569;

# reversionPosition reads the position of a back mutation: at most five digits,
# in [1, RCRS_LENGTH].
func reversionPosition(token as string, digits as string) {
    if (len($digits) > 5) {
        fail("motif token '" + $token + "': position must be in [1, {RCRS_LENGTH}]");
    }
    def p as int init convert.toInt($digits);
    if ($p < 1 or $p > RCRS_LENGTH) {
        fail("motif token '" + $token + "': position must be in [1, {RCRS_LENGTH}]");
    }
    return $p;
}

# The loss of an inserted base is written with a leading minus: "-5899.1"
# drops a base an ancestor inserted at 5899. The index has at most four digits
# and starts at 1, as in `lineage.parseDifference`.
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
    def position as int init reversionPosition($token, $parts[0]);
    if (len($ins) > 4 or convert.toInt($parts[1]) < 1) {
        fail("motif token '" + $token + "': insertion index must be in [1, 9999]");
    }
    return MotifItem{
        reversion: true,
        position: $position,
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
 * @throws {Error} kind "mitoleaf" on an empty token, a leading-minus token
 *   not of the form `-<position>.<insertion>`, a back mutation whose position
 *   is outside [1, 16569] or whose insertion index is outside [1, 9999]; kind
 *   "lineage" on a gain that is not valid rCRS nomenclature
 */
export func parseMotifItem(token as string) {
    def t as string init strings.upper(strings.trim($token));
    if (len($t) == 0) {
        fail("empty motif token");
    }
    # A gain ends in a base or a deletion mark; only a back mutation starts
    # with "-" or ends in a digit, so most tokens skip the character split.
    def last as string init strings.substring($t, len($t) - 1, len($t));
    if (not strings.startsWith($t, "-") and not isDigit($last)) {
        return gain($t);
    }
    def ch as list of string init strings.chars($t);
    if ($ch[0] == "-") {
        return insertionLoss($t, $ch);
    }
    if (isReversion($ch)) {
        return MotifItem{
            reversion: true,
            position: reversionPosition($t, strings.join($ch[1..], "")),
            insertion: 0,
            designation: $t
        };
    }
    return gain($t);
}

# gain reads a gain token through `lineage.parseDifference`.
func gain(t as string) {
    def d as lineage.Difference init lineage.parseDifference($t);
    return MotifItem{
        reversion: false,
        position: $d.position,
        insertion: $d.insertion,
        designation: lineage.formatDifference($d)
    };
}

# motifKey identifies the reference slot a token acts on, so a gain and a later
# reversion at the same place meet in the same slot. Position and insertion are
# zero-padded to the widths `lineage.parseDifference` allows (6 and 4 digits),
# so `lists.sort` on keys orders slots by position, then insertion index.
func motifKey(position as int, insertion as int) {
    def p as string init convert.toString($position);
    def i as string init convert.toString($insertion);
    return strings.repeat("0", 6 - len($p)) + $p + "." + strings.repeat("0", 4 - len($i)) + $i;
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

# An op is one motif token parsed once. A node's ops are stored space-joined
# in one string per node: `Tree.ops` then holds no struct of this module (the
# interpreter rebuilds those at every call across the module boundary) and no
# nested list for its retag check to walk. Layout: a kind letter, then for
# G/R/L the 11-character `motifKey`:
#   "G" + key + designation   a gain
#   "R" + key                 `A263`, a reverting substitution
#   "L" + key                 `-5899.1`, the loss of an inserted base
#   "!" + index               a malformed token, the index-th of the motif,
#                             re-parsed to throw where it is folded
def const KEY_WIDTH as int init 11;

# opOf encodes one token.
func opOf(token as string) {
    def m as MotifItem init parseMotifItem($token);
    def key as string init motifKey($m.position, $m.insertion);
    if (not $m.reversion) {
        return "G" + $key + $m.designation;
    }
    if ($m.insertion == 0) {
        return "R" + $key;
    }
    return "L" + $key;
}

# compileMotif encodes a node's tokens. A malformed token is kept for its
# error, which is raised where the token is folded, so `parse` still accepts a
# tree with a bad token in it.
func compileMotif(motif as list of string) {
    def out as list of string init [];
    for (def i in 0..len($motif)) {
        def op as string init "!{$i}";
        try {
            $op = opOf($motif[$i]);
        } catch (e) {
            $op = "!{$i}";
        }
        $out[] = $op;
    }
    return strings.join($out, " ");
}

# Change is the effect of one op on the accumulated state: `value` is the
# designation left in slot `key`, "" for an emptied slot. An empty `key` means
# the op changes nothing. `insertion` is set for a gained inserted base, so its
# position can be recorded in the `top` bound.
def struct Change {
    key as string,
    value as string,
    insertion as int
};

# change applies one op to the accumulated state, which maps a slot
# (`motifKey`) to its designation, "" for a slot a back mutation emptied. `top`
# maps a padded position to the highest insertion index any gain so far used
# there; it only bounds the search, so it may run ahead of `state`.
#
# A gain sets its slot. `A263` empties slot 263. `-5899.1` drops the last base
# of the insertion run at 5899 - the highest index present, whatever index the
# token names - so a partial loss shortens a run from its end. A back mutation
# with nothing to revert changes nothing: a published tree may restate one.
func change(op as string, state as map of string to string, top as map of string to int) {
    def kind as string init strings.substring($op, 0, 1);
    if ($kind == "!") {
        # throws: the token did not parse when the tree was built
        parseMotifItem(strings.substring($op, 1, len($op)));
    }
    def key as string init strings.substring($op, 1, KEY_WIDTH + 1);
    if ($kind == "G") {
        def ins as string init strings.substring($key, 7, KEY_WIDTH);
        def insertion as int init 0;
        if ($ins != "0000") {
            $insertion = convert.toInt($ins);
        }
        return Change{
            key: $key,
            value: strings.substring($op, KEY_WIDTH + 1, len($op)),
            insertion: $insertion
        };
    }
    def none as Change init Change{key: "", value: "", insertion: 0};
    if ($kind == "R") {
        if (maps.has($state, $key) and len($state[$key]) > 0) {
            return Change{key: $key, value: "", insertion: 0};
        }
        return $none;
    }
    def p as string init strings.substring($key, 0, 6);
    if (not maps.has($top, $p)) {
        return $none;
    }
    for (def i as int init $top[$p]; $i >= 1; $i = $i - 1) {
        def n as string init convert.toString($i);
        def slot as string init $p + "." + strings.repeat("0", 4 - len($n)) + $n;
        if (maps.has($state, $slot) and len($state[$slot]) > 0) {
            return Change{key: $slot, value: "", insertion: 0};
        }
    }
    return $none;
}

# opsOf returns a node's ops, from `Tree.ops` or, for a node it lacks, from
# the motif parsed now. A malformed token's op comes back as "!" + the token.
func opsOf(t as Tree, id as string) {
    def text as string init "";
    if (maps.has($t.ops, $id)) {
        $text = $t.ops[$id];
    } else {
        $text = compileMotif($t.nodes[$id].motif);
    }
    def out as list of string init [];
    if (len($text) == 0) {
        return $out;
    }
    for (def op in strings.split($text, " ")) {
        if (strings.startsWith($op, "!")) {
            def i as int init convert.toInt(strings.substring($op, 1, len($op)));
            $out[] = "!" + $t.nodes[$id].motif[$i];
        } else {
            $out[] = $op;
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
 * tokens are stored verbatim and parsed once into `Tree.ops`; a malformed one
 * is not rejected here but throws later, from `haplotype`, `distance` or
 * `assign`.
 * @param doc {string} XML document text
 * @return {Tree} the parsed phylogeny
 * @throws {Error} kind "mitoleaf" if the document is not well-formed XML,
 *   holds no `<Node>` or more than one top-level `<Node>`, a node has no `Id`,
 *   or a haplogroup name appears twice
 */
export func parse(doc as string) {
    def root as xml.Value;
    try {
        $root = xml.decode($doc);
    } catch (e) {
        fail("malformed XML: " + $e.message);
    }
    def flat as list of Flat init flatten($root);
    if (len($flat) == 0) {
        fail("no <Node> element found in the document");
    }
    def t as Tree init Tree{root: $flat[0].id, nodes: {}, ops: {}};
    for (def f in $flat) {
        if (maps.has($t.nodes, $f.id)) {
            fail("duplicate haplogroup '" + $f.id + "' in the tree");
        }
        # `assign` walks down from `root`, so a second top-level node would be
        # unreachable
        if (len($f.parent) == 0 and $f.id != $t.root) {
            fail("two top-level <Node> elements, '" + $t.root + "' and '" + $f.id +
                "': a tree has one root");
        }
        $t.nodes[$f.id] = Node{id: $f.id, parent: $f.parent, motif: $f.motif, children: []};
        $t.ops[$f.id] = compileMotif($f.motif);
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

/**
 * Build a `Tree` from nodes assembled by hand, filling `ops` from the motifs.
 *
 * Nothing is checked: `parent` and `children` must already link the nodes,
 * and `root` must name one of them, or `assign` fails with a runtime error.
 * Use it after editing a parsed tree's motifs, so `ops` matches them again.
 * @param root {string} root haplogroup
 * @param nodes {map of string to Node} haplogroup name to node
 * @return {Tree} the tree, with `ops` built
 */
export func fromNodes(root as string, nodes as map of string to Node) {
    def t as Tree init Tree{root: $root, nodes: $nodes, ops: {}};
    for (def id in $nodes) {
        $t.ops[$id] = compileMotif($nodes[$id].motif);
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

# stateOf accumulates the motifs along a lineage into slot -> designation, ""
# for an emptied slot (see `change`).
func stateOf(t as Tree, hg as string) {
    def state as map of string to string init {};
    def top as map of string to int init {};
    for (def id in path($t, $hg)) {
        for (def op in opsOf($t, $id)) {
            def c as Change init change($op, $state, $top);
            if (len($c.key) == 0) {
                continue;
            }
            $state[$c.key] = $c.value;
            if ($c.insertion > 0) {
                $top = raiseTop($top, $c);
            }
        }
    }
    return $state;
}

# raiseTop records a gained insertion index in the per-position bound `change`
# searches down from.
func raiseTop(top as map of string to int, c as Change) {
    def p as string init strings.substring($c.key, 0, 6);
    if (not maps.has($top, $p) or $top[$p] < $c.insertion) {
        $top[$p] = $c.insertion;
    }
    return $top;
}

# stateHaplotype builds the MtHaplotype directly from slots sorted by key, which
# is `lineage`'s order (position, then insertion index). `lineage.mtHaplotype`
# would sort through `lists.sortBy`, and a func value in the call tree makes the
# interpreter deep-copy the caller's Tree per call.
func stateHaplotype(id as string, state as map of string to string) {
    def diffs as list of lineage.Difference init [];
    for (def key in lists.sort(maps.keys($state))) {
        if (len($state[$key]) > 0) {
            $diffs[] = lineage.parseDifference($state[$key]);
        }
    }
    return lineage.MtHaplotype{id: $id, ranges: [WHOLE_GENOME], differences: $diffs};
}

/**
 * A haplogroup's full haplotype: its motifs accumulated from the root, with
 * back mutations applied. A back mutation with nothing to revert is ignored;
 * `-<position>.<insertion>` drops the highest-indexed inserted base still
 * present at that position.
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

# --- distance ---------------------------------------------------------------

# observedState keys an observed haplotype the same way a motif state is keyed,
# so the two compare without rebuilding an MtHaplotype per candidate. Two
# designations in one slot (263G 263G, 263G 263A) would count twice in
# `lineage.mtDistance` but once here, so they are rejected.
func observedState(h as lineage.MtHaplotype) {
    def out as map of string to string init {};
    for (def d in $h.differences) {
        def key as string init motifKey($d.position, $d.insertion);
        def des as string init lineage.formatDifference($d);
        if (maps.has($out, $key)) {
            fail("haplotype '" + $h.id + "': '" + $out[$key] + "' and '" + $des +
                "' name the same position");
        }
        $out[$key] = $des;
    }
    return $out;
}

# slotDistance scores one slot: 0 if both sides agree (both empty included), 1
# if only one side carries a designation, 2 if both carry different ones.
func slotDistance(a as string, b as string) {
    if ($a == $b) {
        return 0;
    }
    if (len($a) == 0 or len($b) == 0) {
        return 1;
    }
    return 2;
}

# observedAt reads one slot of an observed state, "" if the slot is empty.
func observedAt(observed as map of string to string, key as string) {
    if (maps.has($observed, $key)) {
        return $observed[$key];
    }
    return "";
}

# stateDistance is `lineage.mtDistance` over two keyed states.
func stateDistance(state as map of string to string, observed as map of string to string) {
    def total as int init 0;
    for (def key in $state) {
        $total = $total + slotDistance($state[$key], observedAt($observed, $key));
    }
    for (def key in $observed) {
        if (not maps.has($state, $key)) {
            $total = $total + 1;
        }
    }
    return $total;
}

/**
 * `lineage.mtDistance` between a haplogroup's haplotype and an observed one:
 * a designation only one side carries counts 1, so a position both carry with
 * different bases counts 2. Ranges are not compared.
 * @param t {Tree} parsed tree
 * @param hg {string} haplogroup name
 * @param h {lineage.MtHaplotype} observed haplotype
 * @return {int} the distance
 * @throws {Error} kind "mitoleaf" if the tree has no such haplogroup, a motif
 *   on the lineage is malformed, or `h` carries two designations at one
 *   position (`263G 263G`, `263G 263A`); kind "lineage" on an invalid gain
 */
export func distance(t as Tree, hg as string, h as lineage.MtHaplotype) {
    def state as map of string to string init stateOf($t, $hg);
    return stateDistance($state, observedState($h));
}

# --- assignment -------------------------------------------------------------

/**
 * The haplogroup whose haplotype is nearest an observed one, by the same
 * measure as `distance`.
 *
 * One depth-first pass in document order scores every haplogroup once. It
 * keeps one motif state and one running distance, updates both only for the
 * slots a node's motif touches and undoes them on the way back up, reading
 * the motifs pre-parsed in `Tree.ops` (about 0.4 s on a 6409-haplogroup
 * tree). Ties go to the first node reached, so an ancestor beats an equally
 * distant descendant.
 *
 * This is a nearest-neighbour report, not a validated haplogroup call: it
 * ignores which positions were sequenced, weights every difference alike and
 * applies none of the phylogenetic rules of a classifier such as EMPOP or
 * HaploGrep. Check the result by hand with `path` and `haplotype`.
 * @param t {Tree} tree from `parse` or `fromNodes`; a tree with no nodes, or
 *   whose `root` is not among its nodes, is a runtime error
 * @param h {lineage.MtHaplotype} observed haplotype
 * @return {Match} the nearest haplogroup and its distance
 * @throws {Error} kind "mitoleaf" if `h` carries two designations at one
 *   position; kind "mitoleaf" or "lineage" on a malformed motif token
 */
export func assign(t as Tree, h as lineage.MtHaplotype) {
    # An explicit preorder stack, not recursion: a recursive call that passes
    # the Tree deep-copies all of it per node once the deck is imported.
    def observed as map of string to string init observedState($h);
    def state as map of string to string init {};
    def top as map of string to int init {};
    # undo log: a slot and the designation it held before ("" for empty)
    def undoKey as list of string init [];
    def undoValue as list of string init [];
    def undoSize as int init 0;
    # Pending frames in parallel scalar lists (a struct frame would be copied
    # out per pop). An enter frame has mark -1; an exit frame restores the undo
    # log to `mark` and the distance to `dist`. `size` is the live depth: slots
    # above it are reused, never sliced off.
    def frameId as list of string init [$t.root];
    def frameMark as list of int init [-1];
    def frameDist as list of int init [0];
    def size as int init 1;
    def d as int init len($observed);
    def best as Match init Match{haplogroup: "", distance: -1};
    while ($size > 0) {
        $size = $size - 1;
        def id as string init $frameId[$size];
        if ($frameMark[$size] >= 0) {
            while ($undoSize > $frameMark[$size]) {
                $undoSize = $undoSize - 1;
                $state[$undoKey[$undoSize]] = $undoValue[$undoSize];
            }
            $d = $frameDist[$size];
            continue;
        }
        def mark as int init $undoSize;
        def before as int init $d;
        for (def op in opsOf($t, $id)) {
            def c as Change init change($op, $state, $top);
            if (len($c.key) == 0) {
                continue;
            }
            def old as string init "";
            if (maps.has($state, $c.key)) {
                $old = $state[$c.key];
            }
            if ($old == $c.value) {
                continue;
            }
            def o as string init observedAt($observed, $c.key);
            $d = $d - slotDistance($old, $o) + slotDistance($c.value, $o);
            if ($undoSize < len($undoKey)) {
                $undoKey[$undoSize] = $c.key;
                $undoValue[$undoSize] = $old;
            } else {
                $undoKey[] = $c.key;
                $undoValue[] = $old;
            }
            $undoSize = $undoSize + 1;
            $state[$c.key] = $c.value;
            if ($c.insertion > 0) {
                $top = raiseTop($top, $c);
            }
        }
        if ($best.distance < 0 or $d < $best.distance) {
            $best = Match{haplogroup: $id, distance: $d};
        }
        # the exit frame goes under the children, so it runs after them
        def kids as list of string init $t.nodes[$id].children;
        for (def i as int init len($kids); $i >= 0; $i = $i - 1) {
            def fid as string init $id;
            def fmark as int init $mark;
            if ($i < len($kids)) {
                $fid = $kids[$i];
                $fmark = -1;
            }
            if ($size < len($frameId)) {
                $frameId[$size] = $fid;
                $frameMark[$size] = $fmark;
                $frameDist[$size] = $before;
            } else {
                $frameId[] = $fid;
                $frameMark[] = $fmark;
                $frameDist[] = $before;
            }
            $size = $size + 1;
        }
    }
    return $best;
}
