# SPDX-License-Identifier: LGPL-3.0-only
# SPDX-FileCopyrightText: Copyright (C) 2026 mplx <mplx@jennifer-lang.dev>
# pragma-jennifer-version: >=0.25.0

# Walk a mitoLEAF phylogeny: lineages, accumulated haplotypes, back mutations,
# and a nearest-haplogroup lookup.
#
#   jennifer run examples/haplogroups.j
#
# It runs on testdata/tree.xml - a seven-haplogroup fixture, not the real
# phylogeny. Point it at a download from https://forensicgenomics.github.io/
# mitoLeaf/ to walk the actual 6409-haplogroup tree.

use io;
use fs;

import "../src/mitoleaf.j" as mitoleaf;
import "@deckyard/forensicgenetics/lineage.j" as lineage;

def t as mitoleaf.Tree init mitoleaf.parse(fs.readString("testdata/tree.xml"));

io.printf("root %s, %d haplogroups\n\n", $t.root, len(mitoleaf.haplogroups($t)));

# --- a lineage and what each step contributes -------------------------------

io.printf("lineage of L0a1:\n");
for (def hg in mitoleaf.path($t, "L0a1")) {
    io.printf("  %s|pad=8 %v\n", $hg, mitoleaf.motif($t, $hg));
}

# --- the accumulated haplotype ----------------------------------------------

# Each step's motif is folded into the one above it. A token like "A263" is a
# back mutation: it removes 263G rather than adding anything, so the root's
# 263G is absent from L0's haplotype and below.
io.printf("\naccumulated haplotypes down that lineage:\n");
for (def hg in mitoleaf.path($t, "L0a1")) {
    def h as lineage.MtHaplotype init mitoleaf.haplotype($t, $hg);
    io.printf("  %s|pad=8 %v\n", $hg, lineage.mtDesignations($h));
}

# --- nearest haplogroup -----------------------------------------------------

# An observed haplotype that is L0a plus one extra difference.
def observed as lineage.MtHaplotype init lineage.mtHaplotype(
    "evidence",
    [mitoleaf.WHOLE_GENOME],
    "73G 315.1C 1048T 5231A 16278T 16519Y 9999A");

def m as mitoleaf.Match init mitoleaf.assign($t, $observed);
io.printf("\nnearest haplogroup: %s (distance %d)\n", $m.haplogroup, $m.distance);
io.printf("its lineage: %v\n", mitoleaf.path($t, $m.haplogroup));

# The distance is symmetric-difference over designations, so it counts what
# the observation lacks as well as what it carries beyond the haplogroup.
io.printf("\ndistance to each haplogroup:\n");
for (def hg in mitoleaf.haplogroups($t)) {
    io.printf("  %s|pad=8 %d\n", $hg, mitoleaf.distance($t, $hg, $observed));
}
