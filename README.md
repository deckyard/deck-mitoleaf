# @deckyard/mitoleaf

[![test](https://github.com/deckyard/deck-mitoleaf/actions/workflows/test.yml/badge.svg)](https://github.com/deckyard/deck-mitoleaf/actions/workflows/test.yml)
[![docs](https://github.com/deckyard/deck-mitoleaf/actions/workflows/docs.yml/badge.svg)](https://github.com/deckyard/deck-mitoleaf/actions/workflows/docs.yml)
[![release](https://img.shields.io/github/v/release/deckyard/deck-mitoleaf)](https://github.com/deckyard/deck-mitoleaf/releases/latest)

Reads the [mitoLEAF](https://forensicgenomics.github.io/mitoLeaf/) human
mitochondrial phylogeny (the maintained successor to PhyloTree), turns a
haplogroup into its defining haplotype with back mutations applied, and finds
the haplogroup nearest an observed haplotype. Haplotypes are
`lineage.MtHaplotype` values from
[`@deckyard/forensicgenetics`](https://github.com/deckyard/deck-forensicgenetics).
Pure Jennifer.

Manual: **<https://deckyard.github.io/deck-mitoleaf/>**

## Install

Requires Jennifer 0.26.0 or newer.

```sh
jvc add "@deckyard/mitoleaf" "^0.2.0"
jvc install
```

jvc vendors the deck to `vendor/deckyard/mitoleaf/` and pulls in
`@deckyard/forensicgenetics` with it. Without jvc, copy the `src/` of both
decks (from a checkout or a `deckyard-mitoleaf-<version>.tar.gz` release
archive) by hand:

```sh
mkdir -p vendor/deckyard/mitoleaf vendor/deckyard/forensicgenetics
cp path/to/deck-mitoleaf/src/*.j         vendor/deckyard/mitoleaf/
cp path/to/deck-forensicgenetics/src/*.j vendor/deckyard/forensicgenetics/
rm -f vendor/deckyard/*/*_test.j
```

## Example

No tree is bundled; download the phm XML from mitoLEAF.

```jennifer
import "@deckyard/mitoleaf/" as mitoleaf;
import "@deckyard/forensicgenetics/lineage.j" as lineage;

def t as mitoleaf.Tree init mitoleaf.parse(fs.readString("mitoLEAF_phm.xml"));

mitoleaf.path($t, "L0a1");                               # root first
mitoleaf.motif($t, "L0a1");                              # this node's own tokens
lineage.mtDesignations(mitoleaf.haplotype($t, "L0a1"));  # accumulated haplotype

def m as mitoleaf.Match init mitoleaf.assign($t, $observed);
```

`assign` is a nearest-neighbour report, not a haplogroup classifier; see
[Limitations](https://deckyard.github.io/deck-mitoleaf/project/limitations.html).

## Tests

From the deck root (the example reads `testdata/tree.xml`):

```sh
jennifer test src/mitoleaf_test.j
jennifer run examples/haplogroups.j
```

## Citing mitoLEAF

> Huber N, Hurmer N, Dür A, Parson W. mitoLEAF: mitochondrial DNA Lineage,
> Evolution, Annotation Framework. *NAR Genomics and Bioinformatics* 7(2),
> 2025. <https://doi.org/10.1093/nargab/lqaf079>

mitoLEAF data is MPL-2.0; none is redistributed here.

## Licence

LGPL-3.0-only, see [LICENSE](LICENSE). Copyright (C) 2026 mplx &lt;mplx@jennifer-lang.dev&gt;.
