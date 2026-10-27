# @deckyard/mitoleaf

Reads the [mitoLEAF](https://forensicgenomics.github.io/mitoLeaf/) human
mitochondrial phylogeny, the maintained successor to PhyloTree (whose updates
stopped in 2016), and turns a haplogroup into the haplotype that defines it.

This deck is part of [Deckyard](https://github.com/deckyard), community decks for [Jennifer](https://jennifer-lang.dev/).

```jennifer
import "@deckyard/mitoleaf/" as mitoleaf;
import "@deckyard/forensicgenetics/lineage.j" as lineage;

def t as mitoleaf.Tree init mitoleaf.parse(fs.readString("mitoLEAF_phm.xml"));
lineage.mtDesignations(mitoleaf.haplotype($t, "L0a1"));
```

Haplotypes are `lineage.MtHaplotype` values from
[`@deckyard/forensicgenetics`](https://deckyard.github.io/deck-forensicgenetics/),
and motif gains are parsed by `lineage.parseDifference`, so rCRS nomenclature
has one implementation and a phylogeny haplotype compares against observed
data with the same `lineage.mtMatches` and `lineage.mtDistance`. The
dependency runs one way: `lineage` knows nothing about mitoLEAF.

The deck bundles no phylogeny, is not a haplogroup classifier and applies none
of EMPOP's phylogenetic alignment rules; see [Limitations](project/limitations.md).

## Citing mitoLEAF

The phylogeny is mitoLEAF's work:

> Huber N, Hurmer N, Dür A, Parson W. mitoLEAF: mitochondrial DNA Lineage,
> Evolution, Annotation Framework. *NAR Genomics and Bioinformatics* 7(2),
> 2025. <https://doi.org/10.1093/nargab/lqaf079>

mitoLEAF data is MPL-2.0 licensed; none is redistributed in this deck.
