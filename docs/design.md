<!-- SPDX-License-Identifier: LGPL-3.0-only -->
<!-- SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev> -->

# Design

## Scope

A reader for the mitoLEAF human mitochondrial phylogeny, the "phm" XML of
nested `<Node>` elements. `parse` builds a haplogroup `Tree`; `path`,
`motif`, `children` and `node` walk it; `haplotype` accumulates a
haplogroup's motifs from the root with back mutations applied; `distance`
compares a haplogroup with an observed haplotype; and `assign` reports the
nearest haplogroups. It is pure computation on a string the caller has read:
no capability, no file or network access.

Not in scope:

- **A classifier.** `assign` is a nearest-neighbour report by symmetric
  difference: no range awareness, no weighting by how informative a position
  is, none of EMPOP's or HaploGrep's phylogenetic reasoning.
  [Limitations](project/limitations.md) spells this out.
- **Phylogenetic alignment** of indels, which neither this deck nor its
  dependency does.
- **Bundled tree data.** mitoLEAF's data is MPL-2.0 and not redistributed;
  `testdata/tree.xml` is a seven-haplogroup synthetic fixture.
- **Forensic validation.** The deck has not been through the validation any
  jurisdiction requires for casework.

## Decisions

- **One nomenclature implementation.** Haplotypes are
  `lineage.MtHaplotype` values from `@deckyard/forensicgenetics`, and gains
  are parsed with its `lineage.parseDifference`, so rCRS notation has one
  implementation. Only `lineage.j` is imported.
- **Back mutations parsed here.** `lineage.parseDifference` rejects both
  back-mutation spellings, so `parseMotifItem` classifies a token itself: a
  leading `-` (`-5899.1`) is the loss of an inserted base; a base and digits
  only (`A263`) is a reverting substitution; anything else is a gain
  (`263G`, `A263G`, `315.1C`, `8281-`, IUPAC codes). `8281-` is a deletion
  gained, `-8281.1` an insertion lost, and `A263G` a redundantly spelled
  gain, not a reversion. Positions are 1 to 16569, insertion indexes 1 to
  9999.
- **Reverting an insertion removes the highest index present** at that
  position; a reversion with nothing to revert does nothing.
- **Lenient structure, strict identity.** `Id` / `HG` attributes and the
  tag name are read case-insensitively; comments and wrapper elements are
  skipped. The first `<Node>` is the root. A duplicate haplogroup or a second
  top-level node (which `assign` could never reach) is an error.
- **Motifs compiled once, errors deferred.** `parse` compiles every node's
  motif into `Tree.ops`, so `haplotype`, `distance` and `assign` never parse
  a token again. A malformed token does not fail `parse`; it throws its
  original error when a call folds that node, so one bad motif deep in a
  release does not make the whole tree unreadable.
- **`Tree.ops` holds plain strings.** One space-joined string of operations
  per haplogroup, rather than a list or a struct of the deck's own, because
  the interpreter checks values crossing the module boundary and nested
  lists made every call measurably slower. `fromNodes` is the supported
  constructor and derives `ops`; a hand-built `Tree` without them falls back
  to compiling motifs on the fly.
- **`distance` agrees with `lineage.mtDistance`.** It compares states keyed
  by position directly instead of building a haplotype, and refuses an
  observation with two designations at one position, the rule of
  forensicgenetics' `mtHaplotype`. Two tests compare it with `mtDistance` and
  with `assign` over every haplogroup, since forensicgenetics does not know
  about this deck and nothing upstream guards the contract.
- **`assign` in one pass.** A preorder walk with an explicit stack applies
  each node's changes to one shared state, keeps a running distance, and
  rolls the changes back on the way out. Ties go to the node reached first,
  so ancestors win.
- **Explicit stacks, never recursion over the `Tree`.** Passing the tree
  into a recursive helper copied it per call: `parse` took 2 min 14 s and
  `assign` 56 s on the 6409-node tree before both moved to explicit stacks.
- **Ranges are not compared.** Phylogeny haplotypes span the whole genome
  (`WHOLE_GENOME`); observations usually cover the control region. Checking
  that is the caller's, with `lineage.mtRangesAgree`.

## Errors

Failures of the deck are `Error{kind: "mitoleaf"}`: malformed XML, no
`<Node>`, a duplicate or second top-level haplogroup, an unknown haplogroup,
a malformed back mutation, two designations at one position in an
observation. Kind `"lineage"` passes through unchanged from the dependency
for an invalid gain, and kind `"limit"` from `xml.decode` for a document past
its nesting or node caps. `assign` on a tree without nodes fails with a runtime
error.

## Security and limits

- The input is a string the caller has read; the deck performs no I/O.
- XML is decoded by `xml.decode`, under its nesting and node caps; a
  reached cap passes on as kind `"limit"`, any other decode failure is
  reported as malformed XML (kind `"mitoleaf"`).
- No deck-level size cap: the published tree is about 550 KB and 6409
  haplogroups.

## Testing and validation

- 43 tests on the synthetic seven-haplogroup tree, which exercises an
  insertion, an IUPAC ambiguity, a deletion and both back-mutation
  spellings. The tree is inlined in the tests so an installed deck runs them
  without `testdata/`; a test keeps it byte-equal to `testdata/tree.xml`.
- Agreement with `lineage.mtDistance` and between `assign` and `distance` is
  tested over every haplogroup.
- Timings were taken on the published 6409-haplogroup tree: `parse` about
  0.7 s, `assign` about 0.4 s, and every call passing the `Tree` about 5.5 ms
  ([figures](guide/assignment.md#performance)).
- The tests need `vendor/` populated (`tools/vendor-deps.sh`, network),
  which fetches the released forensicgenetics through the registry.
- Not validated against EMPOP or HaploGrep assignments.

## Known limitations and open questions

- Every call from an importing program that passes a `Tree` costs about
  5.5 ms on the published tree before doing any work, because the
  interpreter rebuilds values holding a module's own structs at the
  boundary. A loop calling `distance` per haplogroup takes about a minute;
  `assign` scores the whole tree in one call. The fix belongs to the
  interpreter (`jennifer`'s M29.14).
- Code reachable with the `Tree` in its call chain must not use a callback
  builtin such as `lists.sortBy`: it made every `haplotype` and `distance`
  call copy the whole tree. The deck sorts plain keys instead
  (`jennifer`'s M30.8 would lift this).
- `Tree.ops` is derived from the motifs and not checked against them;
  editing `Node.motif` without `fromNodes` leaves results on the old motif.
- The deck does not record which mitoLEAF release a tree came from; the phm
  file's leading comment states it.

## Dependencies on the interpreter

- Floor `>=0.26.0`, the default binary (no `jennifer-tiny` engine declared).
- `xml.decode`; `@deckyard/forensicgenetics` `>=0.1.0` through `[decks]`.
- Performance depends on `jennifer`'s M29.14 (cheaper calls into another
  module) and M30.8 (borrow through named callbacks).
