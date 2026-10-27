# Installing the deck

## Requirements

Jennifer 0.25.0 or newer. Every source file declares the floor, and the
interpreter refuses to load the deck on an older build:

```jennifer
# pragma-jennifer-version: >=0.25.0
```

One deck dependency,
[`@deckyard/forensicgenetics`](https://github.com/deckyard/deck-forensicgenetics),
for `lineage.MtHaplotype` and `lineage.parseDifference`.

## With jvc

From a project with a `deck.toml`:

```sh
jvc add "@deckyard/mitoleaf" "^0.1.0"
jvc install
```

`add` writes the requirement; `install` resolves the graph, verifies each
download against its `sha256` and vendors it:

```
your-project/
  app.j
  deck.toml
  camcorder.lock
  vendor/
    deckyard/
      mitoleaf/
        mitoleaf.j
      forensicgenetics/          <- resolved transitively
        forensicgenetics.j
        lineage.j
        kinship.j
        pedigree.j
```

`@deckyard/forensicgenetics` is installed because this deck's manifest requires
it, and its constraint is unified with anything else in your graph that needs
it. Only `src/` is vendored, without the `*_test.j` overlays. `camcorder.lock`
pins the resolved set, so `git clone` + `jvc install` reproduces it exactly;
`jvc update` moves it forward. `jvc check` verifies the running interpreter
satisfies the deck's `[engines]`.

The imports then resolve with no path:

```jennifer
import "@deckyard/mitoleaf/" as mitoleaf;
import "@deckyard/forensicgenetics/lineage.j" as lineage;
```

The trailing `/` on the first import expands to the package-named entry file,
`mitoleaf/mitoleaf.j`. A sibling module is named outright.

### From a git remote

If the deck is not on the registry you resolve against, point the requirement
at its repository. `[decks]` says which version, `[sources]` where from:

```sh
jvc source "@deckyard/mitoleaf" https://github.com/deckyard/deck-mitoleaf.git
jvc install
```

```toml
[decks]
"@deckyard/mitoleaf" = "^0.1.0"

[sources]
"@deckyard/mitoleaf" = "https://github.com/deckyard/deck-mitoleaf.git"
```

Versions are the repository's SemVer tags (this deck tags bare, `0.1.0`; jvc
also accepts a `v` prefix), and the tag must agree with the manifest it
carries. `camcorder.lock` pins the commit, so moving a tag later cannot change
what a lockfile installs. This needs `git` on `PATH` and an interpreter with
the `exec` capability.

## By hand into a vendor tree

The `@scope/package` resolver reads a plain directory, so a hand-populated tree
works without jvc or a network. Take `src/` from a checkout or from a release
archive (`deckyard-mitoleaf-<version>.tar.gz`):

```sh
mkdir -p vendor/deckyard/mitoleaf vendor/deckyard/forensicgenetics
cp path/to/deck-mitoleaf/src/*.j         vendor/deckyard/mitoleaf/
cp path/to/deck-forensicgenetics/src/*.j vendor/deckyard/forensicgenetics/
rm -f vendor/deckyard/*/*_test.j
```

The installed tree is flat: the files from `src/` sit at the package root. You
own the dependency's version; there is no lockfile and nothing checks that the
two decks agree.

The resolver takes the vendor root from `--vendor`, then `JENNIFER_VENDOR`,
then the nearest `vendor/` directory above the program being run.

## By relative path

For work on a checkout, import the source directly:

```jennifer
import "./src/mitoleaf.j" as mitoleaf;
```

The bundled example does this. Struct identity is keyed by the file's canonical
path, so a `Tree` built through either import is the same type.

## Getting the tree

Download the phm XML from
[mitoLEAF](https://forensicgenomics.github.io/mitoLeaf/) (at the time of
writing `docs/data/mitoLEAF_phm.xml`, 549 KB, 6409 haplogroups) and pass its
text to `parse`. The file's leading comments are a changelog of variants added,
removed and respelled; record which release you used with any result.

## Checking the install

From the deck root:

```sh
jennifer test src/mitoleaf_test.j
jennifer run examples/haplogroups.j
```

Expect 29 tests and `0 failed`. The example reads `testdata/tree.xml` relative
to the working directory, so run it from the deck root; the tests carry their
own copy of the tree and run from anywhere, including an installed deck.
