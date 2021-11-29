![](logo.svg)

# Project description

This is an application to analyze base pairing patterns of DNA/RNA 3D
structures to find and classify tetrads and quadruplexes. ElTetrado
assigns tetrads to one of the ONZ classes (O, N, Z) alongside with the
directionality of the tetrad (+/-) determined by the bonds between bases
and their non-canonical interactions. The interactions follow
Leontis/Westhof classification (Leontis and Westhof, 2001). Watson-Crick
(W) edge of first base in the tetrad structure exposed to the Hoogsteen
(H) edge of the next nucleobase from the same tetrad sets the tetrad
directionality, clockwise (+) or anticlockwise (-). For more details,
please refer to Zok *et al.* (2020) and Popenda *et al.* (2020)

# Installation

Please run:

    pip install eltetrado

If you have both Python 2 and Python 3 installed, you need to explicitly
call `pip3`:

    pip3 install eltetrado

# Dependencies

The project is written in Python 3.6+ and requires
[mmcif-pdbx](https://github.com/Electrostatics/mmcif_pdbx),
[orjson](https://github.com/ijl/orjson) and [NumPy](https://numpy.org/).

ElTetrado depends on DSSR (Lu *et al.*, 2015) in terms of detection of
base pairing and stacking. The binary `x3dna-dssr` can be
[downloaded](http://forum.x3dna.org/site-announcements/download-instructions/)
and it needs to be put in `$PATH` i.e. typically in `/usr/bin`.
Alternatively, one can pre-process the 3D data with `x3dna-dssr --json`
and provide the path to a JSON result as an input to ElTetrado (see
Usage section below).

Visualization is created by `R` 3.6+ script which uses
[R4RNA](https://www.e-rna.org/r-chie/) (Lai *et al.*, 2012) library. The
dependency will be automatically installed if not present.

# Usage

ElTetrado is a command line application, which requires to be provided
with:

-   either `--dssr-json` and the path to JSON generated with
    `x3dna-dssr --json` (fast, but quadruplex parameters like `rise` or
    `twist` are not calculated)
-   or `--pdb` and the path to PDB or PDBx/mmCIF file (slow, because the
    execution time is a sum of ElTetrado and DSSR times)
-   or both `--pdb` and `--dssr-json` (recommended, all analyses are
    made and only ElTetrado is executed while DSSR results are read from
    the file)

By default, ElTetrado outputs textual results on the standard output. A
JSON version of the output can be obtained with `--output` switch
followed by a path where the file is supposed to be created.

ElTetrado prepares visualization of the whole structure and of each
N4-helices, quadruplexes and tetrads. This can be supplemented with
canonical base pairs visualization when `--complete-2d` is set. All
color settings are located in the first several lines of the `quadraw.R`
file, you can easily change them without knowledge of R language. If you
want ElTetrado to not visualize anything, pass `--no-image` switch to
it.

    usage: eltetrado [-h] [--pdb PDB] [--dssr-json DSSR_JSON] [--output OUTPUT]
                     [--stacking-mismatch STACKING_MISMATCH] [--strict]
                     [--no-reorder] [--complete-2d] [--no-image] [--version]

    optional arguments:
      -h, --help            show this help message and exit
      --pdb PDB             path to input PDB or PDBx/mmCIF file
      --dssr-json DSSR_JSON
                            path to input JSON file generated with `x3dna-dssr
                            --json`
      --output OUTPUT       (optional) path for output JSON file
      --stacking-mismatch STACKING_MISMATCH
                            a perfect tetrad stacking covers 4 nucleotides; this
                            option can be used with value 1 or 2 to allow this
                            number of nucleotides to be non-stacked with otherwise
                            well aligned tetrad [default=2]
      --strict              nucleotides in tetrad are found when linked only by
                            cWH pairing
      --no-reorder          chains of bi- and tetramolecular quadruplexes are
                            reordered to be able to have them classified; when
                            this is set, chains will be processed in original
                            order and bi-/tetramolecular quadruplexes will not be
                            classified
      --complete-2d         when set, the visualization will also show canonical
                            base pairs to provide context for the quadruplex
      --no-image            when set, the visualization will not be created at all
      --version             show program's version number and exit

# Chains reorder

ElTetrado keeps a global and unique 5’-3’ index for every nucleotide
which is independent from residue numbers. For example, if a structure
has chain M with 60 nucleotides and chain N with 15 nucleotides, then
ElTetrado will keep index between 0 and 74 which uniquely identifies
every nucleotide. Initially, ElTetrado assigns this indices according to
the order of chains in the input file. Therefore, if M preceded N then
nucleotides in M will be indexed from 0 to 59 and in N from 60 to 74.
Otherwise, nucleotides in N will be indexed from 0 to 14 and in M from
15 to 74.

When `--no-reorder` is present, this initial assignment is used.
Otherwise, ElTetrado exhaustively checks all permutations of chains’
orders. Every permutation check induces recalculation of the global and
unique 5’-3’ index and in effect it changes ONZ classification of
tetrads.

ElTetrado keeps a table of tetrad classification scores according to
these rules:

-   Type preference: `O` > `N` > `Z`
-   Direction preference: `+` > `-`

The table keeps low values for preferred classes i.e. `O+` is 0, `O-` is
1 and so on up to `Z-` with score 5. For every permutation of chain
orders, ElTetrado computes sum of scores for tetrads classification
induced by 5’-3’ indexing. We select permutation with the minimum value.

# Examples

## 1MY9: Solution structure of a K+ cation stabilized dimeric RNA quadruplex containing two G:G(:A):G:G(:A) hexads, G:G:G:G tetrads and UUUU loops

![](1my9.png)

    $ curl ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/my/1my9.cif.gz | gzip -d > 1my9.cif
    $ ./eltetrado --pdb 1my9.cif --output 1my9.json

    Chain order: A, B
    n4-helix with 4 tetrads
      Op+ VIII 1a quadruplex with 2 tetrads
        A.G2 A.G5 A.G11 A.G14 cWH-cWH-cWH-cWH O+ VIIIa planarity=0.5 ions_channel= ions_outside={}
          direction=parallel rise=4.12 twist=31.08
        A.G1 A.G4 A.G10 A.G13 cWH-cWH-cWH-cWH O+ VIIIa planarity=0.63 ions_channel= ions_outside={}

        Tracts:
          A.G2, A.G1
          A.G5, A.G4
          A.G11, A.G10
          A.G14, A.G13

        Loops:
          propeller- A.A3
          propeller- A.U6, A.U7, A.U8, A.U9
          propeller- A.A12

      Op+ VIII 1a quadruplex with 2 tetrads
        B.G15 B.G18 B.G24 B.G27 cWH-cWH-cWH-cWH O+ VIIIa planarity=0.26 ions_channel= ions_outside={}
          direction=parallel rise=4.21 twist=34.69
        B.G16 B.G19 B.G25 B.G28 cWH-cWH-cWH-cWH O+ VIIIa planarity=0.31 ions_channel= ions_outside={}

        Tracts:
          B.G15, B.G16
          B.G18, B.G19
          B.G24, B.G25
          B.G27, B.G28

        Loops:
          propeller- B.A17
          propeller- B.U20, B.U21, B.U22, B.U23
          propeller- B.A26

    GGAGGUUUUGGAGG-GGAGGUUUUGGAGG
    ([.)]....([.)]-([.)]....([.)]
    ([.([....)].)]-([.([....)].)]

<details>
<summary>
Click to see the output JSON
</summary>

``` json
{
  "metals": [],
  "nucleotides": [
    {
      "index": 1,
      "model": 1,
      "chain": "A",
      "number": 1,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.G1",
      "shortName": "G",
      "chi": -114.063,
      "glycosidicBond": "anti"
    },
    {
      "index": 2,
      "model": 1,
      "chain": "A",
      "number": 2,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.G2",
      "shortName": "G",
      "chi": -136.905,
      "glycosidicBond": "anti"
    },
    {
      "index": 3,
      "model": 1,
      "chain": "A",
      "number": 3,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.A3",
      "shortName": "A",
      "chi": -53.884,
      "glycosidicBond": "syn"
    },
    {
      "index": 4,
      "model": 1,
      "chain": "A",
      "number": 4,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.G4",
      "shortName": "G",
      "chi": 167.763,
      "glycosidicBond": "anti"
    },
    {
      "index": 5,
      "model": 1,
      "chain": "A",
      "number": 5,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.G5",
      "shortName": "G",
      "chi": -98.198,
      "glycosidicBond": "anti"
    },
    {
      "index": 6,
      "model": 1,
      "chain": "A",
      "number": 6,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.U6",
      "shortName": "U",
      "chi": -150.069,
      "glycosidicBond": "anti"
    },
    {
      "index": 7,
      "model": 1,
      "chain": "A",
      "number": 7,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.U7",
      "shortName": "U",
      "chi": -130.523,
      "glycosidicBond": "anti"
    },
    {
      "index": 8,
      "model": 1,
      "chain": "A",
      "number": 8,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.U8",
      "shortName": "U",
      "chi": -158.504,
      "glycosidicBond": "anti"
    },
    {
      "index": 9,
      "model": 1,
      "chain": "A",
      "number": 9,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.U9",
      "shortName": "U",
      "chi": -149.743,
      "glycosidicBond": "anti"
    },
    {
      "index": 10,
      "model": 1,
      "chain": "A",
      "number": 10,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.G10",
      "shortName": "G",
      "chi": -113.245,
      "glycosidicBond": "anti"
    },
    {
      "index": 11,
      "model": 1,
      "chain": "A",
      "number": 11,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.G11",
      "shortName": "G",
      "chi": -138.466,
      "glycosidicBond": "anti"
    },
    {
      "index": 12,
      "model": 1,
      "chain": "A",
      "number": 12,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.A12",
      "shortName": "A",
      "chi": -70.627,
      "glycosidicBond": "syn"
    },
    {
      "index": 13,
      "model": 1,
      "chain": "A",
      "number": 13,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.G13",
      "shortName": "G",
      "chi": 150.585,
      "glycosidicBond": "anti"
    },
    {
      "index": 14,
      "model": 1,
      "chain": "A",
      "number": 14,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.G14",
      "shortName": "G",
      "chi": -158.594,
      "glycosidicBond": "anti"
    },
    {
      "index": 15,
      "model": 1,
      "chain": "B",
      "number": 15,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.G15",
      "shortName": "G",
      "chi": -114.303,
      "glycosidicBond": "anti"
    },
    {
      "index": 16,
      "model": 1,
      "chain": "B",
      "number": 16,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.G16",
      "shortName": "G",
      "chi": -134.388,
      "glycosidicBond": "anti"
    },
    {
      "index": 17,
      "model": 1,
      "chain": "B",
      "number": 17,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.A17",
      "shortName": "A",
      "chi": -52.702,
      "glycosidicBond": "syn"
    },
    {
      "index": 18,
      "model": 1,
      "chain": "B",
      "number": 18,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.G18",
      "shortName": "G",
      "chi": 156.482,
      "glycosidicBond": "anti"
    },
    {
      "index": 19,
      "model": 1,
      "chain": "B",
      "number": 19,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.G19",
      "shortName": "G",
      "chi": -100.234,
      "glycosidicBond": "anti"
    },
    {
      "index": 20,
      "model": 1,
      "chain": "B",
      "number": 20,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.U20",
      "shortName": "U",
      "chi": -146.822,
      "glycosidicBond": "anti"
    },
    {
      "index": 21,
      "model": 1,
      "chain": "B",
      "number": 21,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.U21",
      "shortName": "U",
      "chi": -144.09,
      "glycosidicBond": "anti"
    },
    {
      "index": 22,
      "model": 1,
      "chain": "B",
      "number": 22,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.U22",
      "shortName": "U",
      "chi": -160.945,
      "glycosidicBond": "anti"
    },
    {
      "index": 23,
      "model": 1,
      "chain": "B",
      "number": 23,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.U23",
      "shortName": "U",
      "chi": -144.171,
      "glycosidicBond": "anti"
    },
    {
      "index": 24,
      "model": 1,
      "chain": "B",
      "number": 24,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.G24",
      "shortName": "G",
      "chi": -121.507,
      "glycosidicBond": "anti"
    },
    {
      "index": 25,
      "model": 1,
      "chain": "B",
      "number": 25,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.G25",
      "shortName": "G",
      "chi": -132.739,
      "glycosidicBond": "anti"
    },
    {
      "index": 26,
      "model": 1,
      "chain": "B",
      "number": 26,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.A26",
      "shortName": "A",
      "chi": -72.432,
      "glycosidicBond": "syn"
    },
    {
      "index": 27,
      "model": 1,
      "chain": "B",
      "number": 27,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.G27",
      "shortName": "G",
      "chi": -178.717,
      "glycosidicBond": "anti"
    },
    {
      "index": 28,
      "model": 1,
      "chain": "B",
      "number": 28,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.G28",
      "shortName": "G",
      "chi": -130.902,
      "glycosidicBond": "anti"
    }
  ],
  "basePairs": [
    {
      "nt1": "A.G1",
      "nt2": "A.A3",
      "lw": "tSH",
      "stericity": "trans",
      "edge5": "Sugar",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.A3",
      "nt2": "A.G1",
      "lw": "tHS",
      "stericity": "trans",
      "edge5": "Hoogsteen",
      "edge3": "Sugar"
    },
    {
      "nt1": "A.G1",
      "nt2": "A.G4",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.G4",
      "nt2": "A.G1",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "A.G1",
      "nt2": "A.G13",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "A.G13",
      "nt2": "A.G1",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.G2",
      "nt2": "A.G5",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.G5",
      "nt2": "A.G2",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "A.G2",
      "nt2": "A.G11",
      "lw": "tWW",
      "stericity": "trans",
      "edge5": "Watson-Crick",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "A.G11",
      "nt2": "A.G2",
      "lw": "tWW",
      "stericity": "trans",
      "edge5": "Watson-Crick",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "A.G2",
      "nt2": "A.G14",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "A.G14",
      "nt2": "A.G2",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.G4",
      "nt2": "A.G10",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.G10",
      "nt2": "A.G4",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "A.G4",
      "nt2": "B.A17",
      "lw": "tHH",
      "stericity": "trans",
      "edge5": "Hoogsteen",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.A17",
      "nt2": "A.G4",
      "lw": "tHH",
      "stericity": "trans",
      "edge5": "Hoogsteen",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.G5",
      "nt2": "A.G11",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.G11",
      "nt2": "A.G5",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "A.G10",
      "nt2": "A.A12",
      "lw": "tSH",
      "stericity": "trans",
      "edge5": "Sugar",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.A12",
      "nt2": "A.G10",
      "lw": "tHS",
      "stericity": "trans",
      "edge5": "Hoogsteen",
      "edge3": "Sugar"
    },
    {
      "nt1": "A.G10",
      "nt2": "A.G13",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.G13",
      "nt2": "A.G10",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "A.G11",
      "nt2": "A.G14",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.G14",
      "nt2": "A.G11",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "A.G13",
      "nt2": "B.G27",
      "lw": "cHH",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.G27",
      "nt2": "A.G13",
      "lw": "cHH",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.G15",
      "nt2": "B.A17",
      "lw": "tSH",
      "stericity": "trans",
      "edge5": "Sugar",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.A17",
      "nt2": "B.G15",
      "lw": "tHS",
      "stericity": "trans",
      "edge5": "Hoogsteen",
      "edge3": "Sugar"
    },
    {
      "nt1": "B.G15",
      "nt2": "B.G18",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.G18",
      "nt2": "B.G15",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "B.G15",
      "nt2": "B.G27",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "B.G27",
      "nt2": "B.G15",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.G16",
      "nt2": "B.G19",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.G19",
      "nt2": "B.G16",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "B.G16",
      "nt2": "B.G28",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "B.G28",
      "nt2": "B.G16",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.G18",
      "nt2": "B.G24",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.G24",
      "nt2": "B.G18",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "B.G19",
      "nt2": "B.G25",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.G25",
      "nt2": "B.G19",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "B.G24",
      "nt2": "B.A26",
      "lw": "tSH",
      "stericity": "trans",
      "edge5": "Sugar",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.A26",
      "nt2": "B.G24",
      "lw": "tHS",
      "stericity": "trans",
      "edge5": "Hoogsteen",
      "edge3": "Sugar"
    },
    {
      "nt1": "B.G24",
      "nt2": "B.G27",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.G27",
      "nt2": "B.G24",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "B.G25",
      "nt2": "B.G28",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.G28",
      "nt2": "B.G25",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    }
  ],
  "helices": [
    {
      "quadruplexes": [
        {
          "tetrads": [
            {
              "id": "A.G2-A.G5-A.G11-A.G14",
              "nt1": "A.G2",
              "nt2": "A.G5",
              "nt3": "A.G11",
              "nt4": "A.G14",
              "onz": "O+",
              "gbaClassification": "VIIIa",
              "planarityDeviation": 0.49745508591228504,
              "ionsChannel": [],
              "ionsOutside": []
            },
            {
              "id": "A.G1-A.G4-A.G10-A.G13",
              "nt1": "A.G1",
              "nt2": "A.G4",
              "nt3": "A.G10",
              "nt4": "A.G13",
              "onz": "O+",
              "gbaClassification": "VIIIa",
              "planarityDeviation": 0.6331162511735106,
              "ionsChannel": [],
              "ionsOutside": []
            }
          ],
          "onzm": "Op+",
          "loopClassification": "1a",
          "gbaClassification": [
            "VIII"
          ],
          "tracts": [
            [
              "A.G2",
              "A.G1"
            ],
            [
              "A.G5",
              "A.G4"
            ],
            [
              "A.G11",
              "A.G10"
            ],
            [
              "A.G14",
              "A.G13"
            ]
          ],
          "loops": [
            {
              "type": "propeller-",
              "nucleotides": [
                "A.A3"
              ]
            },
            {
              "type": "propeller-",
              "nucleotides": [
                "A.U6",
                "A.U7",
                "A.U8",
                "A.U9"
              ]
            },
            {
              "type": "propeller-",
              "nucleotides": [
                "A.A12"
              ]
            }
          ]
        },
        {
          "tetrads": [
            {
              "id": "B.G15-B.G18-B.G24-B.G27",
              "nt1": "B.G15",
              "nt2": "B.G18",
              "nt3": "B.G24",
              "nt4": "B.G27",
              "onz": "O+",
              "gbaClassification": "VIIIa",
              "planarityDeviation": 0.26169770251188723,
              "ionsChannel": [],
              "ionsOutside": []
            },
            {
              "id": "B.G16-B.G19-B.G25-B.G28",
              "nt1": "B.G16",
              "nt2": "B.G19",
              "nt3": "B.G25",
              "nt4": "B.G28",
              "onz": "O+",
              "gbaClassification": "VIIIa",
              "planarityDeviation": 0.31002671255877334,
              "ionsChannel": [],
              "ionsOutside": []
            }
          ],
          "onzm": "Op+",
          "loopClassification": "1a",
          "gbaClassification": [
            "VIII"
          ],
          "tracts": [
            [
              "B.G15",
              "B.G16"
            ],
            [
              "B.G18",
              "B.G19"
            ],
            [
              "B.G24",
              "B.G25"
            ],
            [
              "B.G27",
              "B.G28"
            ]
          ],
          "loops": [
            {
              "type": "propeller-",
              "nucleotides": [
                "B.A17"
              ]
            },
            {
              "type": "propeller-",
              "nucleotides": [
                "B.U20",
                "B.U21",
                "B.U22",
                "B.U23"
              ]
            },
            {
              "type": "propeller-",
              "nucleotides": [
                "B.A26"
              ]
            }
          ]
        }
      ],
      "tetradPairs": [
        {
          "tetrad1": "A.G2-A.G5-A.G11-A.G14",
          "tetrad2": "A.G1-A.G4-A.G10-A.G13",
          "direction": "parallel",
          "rise": 4.121721951290505,
          "twist": 31.077291799731125
        },
        {
          "tetrad1": "A.G1-A.G4-A.G10-A.G13",
          "tetrad2": "B.G15-B.G18-B.G24-B.G27",
          "direction": "parallel",
          "rise": 3.184510372490565,
          "twist": 9.27463267295774
        },
        {
          "tetrad1": "B.G15-B.G18-B.G24-B.G27",
          "tetrad2": "B.G16-B.G19-B.G25-B.G28",
          "direction": "parallel",
          "rise": 4.211160606798321,
          "twist": 34.69391705015754
        }
      ]
    }
  ],
  "dotBracket": {
    "sequence": "GGAGGUUUUGGAGG-GGAGGUUUUGGAGG",
    "line1": "([.)]....([.)]-([.)]....([.)]",
    "line2": "([.([....)].)]-([.([....)].)]"
  }
}
```

</details>

## 4RJ1: Structural variations and solvent structure of UGGGGU quadruplexes stabilized by Sr2+ ions

![](4rj1.png)

    $ curl https://www.ebi.ac.uk/pdbe/static/entry/download/4rj1-assembly-1.cif.gz | gzip -d > 4rj1-1.cif
    $ ./eltetrado --pdb 4rj1-1.cif --output 4rj1-1.json

    Chain order: A, AB, AA, AC, B, BC, BA, BB
    n4-helix with 10 tetrads
      Op* VIII n/a quadruplex with 5 tetrads
        A.U1006 AC.U1006 AA.U1006 AB.U1006 cWH-cWH-cWH-cWH O- VIIIa planarity=1.06 ions_channel=NA ions_outside={A.U1006: 'SR', AA.U1006: 'SR', AB.U1006: 'SR', AC.U1006: 'SR'}
          direction=parallel rise=3.37 twist=39.96
        A.G1005 AB.G1005 AA.G1005 AC.G1005 cWH-cWH-cWH-cWH O+ VIIIa planarity=0.8 ions_channel=SR ions_outside={}
          direction=parallel rise=3.31 twist=25.9
        A.G1004 AB.G1004 AA.G1004 AC.G1004 cWH-cWH-cWH-cWH O+ VIIIa planarity=0.41 ions_channel= ions_outside={}
          direction=parallel rise=3.34 twist=35.81
        A.G1003 AB.G1003 AA.G1003 AC.G1003 cWH-cWH-cWH-cWH O+ VIIIa planarity=0.55 ions_channel=SR ions_outside={}
          direction=parallel rise=3.29 twist=27.12
        A.G1002 AB.G1002 AA.G1002 AC.G1002 cWH-cWH-cWH-cWH O+ VIIIa planarity=0.54 ions_channel= ions_outside={}

        Tracts:
          A.U1006, A.G1005, A.G1004, A.G1003, A.G1002
          AC.U1006, AC.G1005, AC.G1004, AC.G1003, AC.G1002
          AA.U1006, AA.G1005, AA.G1004, AA.G1003, AA.G1002
          AB.U1006, AB.G1005, AB.G1004, AB.G1003, AB.G1002

      Op* VIII n/a quadruplex with 5 tetrads
        B.G2002 BC.G2002 BA.G2002 BB.G2002 cWH-cWH-cWH-cWH O+ VIIIa planarity=0.67 ions_channel= ions_outside={}
          direction=parallel rise=3.37 twist=27.41
        B.G2003 BC.G2003 BA.G2003 BB.G2003 cWH-cWH-cWH-cWH O+ VIIIa planarity=0.58 ions_channel=SR ions_outside={}
          direction=parallel rise=3.32 twist=35.04
        B.G2004 BC.G2004 BA.G2004 BB.G2004 cWH-cWH-cWH-cWH O+ VIIIa planarity=0.23 ions_channel=SR ions_outside={}
          direction=parallel rise=3.27 twist=25.15
        B.G2005 BC.G2005 BA.G2005 BB.G2005 cWH-cWH-cWH-cWH O+ VIIIa planarity=0.78 ions_channel=NA ions_outside={}
          direction=parallel rise=7.14 twist=43.41
        B.U2006 BB.U2006 BA.U2006 BC.U2006 cWH-cWH-cWH-cWH O- VIIIa planarity=1.58 ions_channel=NA,NA ions_outside={}

        Tracts:
          B.G2002, B.G2003, B.G2004, B.G2005, B.U2006
          BC.G2002, BC.G2003, BC.G2004, BC.G2005, BC.U2006
          BA.G2002, BA.G2003, BA.G2004, BA.G2005, BA.U2006
          BB.G2002, BB.G2003, BB.G2004, BB.G2005, BB.U2006

    UGGGGU-UGGGGU-UGGGGU-UGGGGU-UGGGGU-UGGGGU-UGGGGU-UGGGGU
    .([{<A-.)]}>A-.([{<a-.)]}>a-.([{<A-.)]}>A-.([{<a-.)]}>a
    .([{<A-.([{<a-.)]}>A-.)]}>a-.([{<A-.([{<a-.)]}>A-.)]}>a

<details>
<summary>
Click to see the output JSON
</summary>

``` json
{
  "metals": [
    {
      "symbol": "Sr",
      "count": 8
    },
    {
      "symbol": "Na",
      "count": 4
    },
    {
      "symbol": "Ca",
      "count": 12
    }
  ],
  "nucleotides": [
    {
      "index": 1,
      "model": 1,
      "chain": "A",
      "number": 1001,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.U1001",
      "shortName": "U",
      "chi": -141.927,
      "glycosidicBond": "anti"
    },
    {
      "index": 2,
      "model": 1,
      "chain": "A",
      "number": 1002,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.G1002",
      "shortName": "G",
      "chi": -165.93,
      "glycosidicBond": "anti"
    },
    {
      "index": 3,
      "model": 1,
      "chain": "A",
      "number": 1003,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.G1003",
      "shortName": "G",
      "chi": -121.565,
      "glycosidicBond": "anti"
    },
    {
      "index": 4,
      "model": 1,
      "chain": "A",
      "number": 1004,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.G1004",
      "shortName": "G",
      "chi": -156.01,
      "glycosidicBond": "anti"
    },
    {
      "index": 5,
      "model": 1,
      "chain": "A",
      "number": 1005,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.G1005",
      "shortName": "G",
      "chi": -148.101,
      "glycosidicBond": "anti"
    },
    {
      "index": 6,
      "model": 1,
      "chain": "A",
      "number": 1006,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "A.U1006",
      "shortName": "U",
      "chi": -137.28,
      "glycosidicBond": "anti"
    },
    {
      "index": 13,
      "model": 1,
      "chain": "AA",
      "number": 1001,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AA.U1001",
      "shortName": "U",
      "chi": -141.927,
      "glycosidicBond": "anti"
    },
    {
      "index": 14,
      "model": 1,
      "chain": "AA",
      "number": 1002,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AA.G1002",
      "shortName": "G",
      "chi": -165.93,
      "glycosidicBond": "anti"
    },
    {
      "index": 15,
      "model": 1,
      "chain": "AA",
      "number": 1003,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AA.G1003",
      "shortName": "G",
      "chi": -121.565,
      "glycosidicBond": "anti"
    },
    {
      "index": 16,
      "model": 1,
      "chain": "AA",
      "number": 1004,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AA.G1004",
      "shortName": "G",
      "chi": -156.01,
      "glycosidicBond": "anti"
    },
    {
      "index": 17,
      "model": 1,
      "chain": "AA",
      "number": 1005,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AA.G1005",
      "shortName": "G",
      "chi": -148.101,
      "glycosidicBond": "anti"
    },
    {
      "index": 18,
      "model": 1,
      "chain": "AA",
      "number": 1006,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AA.U1006",
      "shortName": "U",
      "chi": -137.28,
      "glycosidicBond": "anti"
    },
    {
      "index": 7,
      "model": 1,
      "chain": "AB",
      "number": 1001,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AB.U1001",
      "shortName": "U",
      "chi": -141.927,
      "glycosidicBond": "anti"
    },
    {
      "index": 8,
      "model": 1,
      "chain": "AB",
      "number": 1002,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AB.G1002",
      "shortName": "G",
      "chi": -165.93,
      "glycosidicBond": "anti"
    },
    {
      "index": 9,
      "model": 1,
      "chain": "AB",
      "number": 1003,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AB.G1003",
      "shortName": "G",
      "chi": -121.565,
      "glycosidicBond": "anti"
    },
    {
      "index": 10,
      "model": 1,
      "chain": "AB",
      "number": 1004,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AB.G1004",
      "shortName": "G",
      "chi": -156.01,
      "glycosidicBond": "anti"
    },
    {
      "index": 11,
      "model": 1,
      "chain": "AB",
      "number": 1005,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AB.G1005",
      "shortName": "G",
      "chi": -148.101,
      "glycosidicBond": "anti"
    },
    {
      "index": 12,
      "model": 1,
      "chain": "AB",
      "number": 1006,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AB.U1006",
      "shortName": "U",
      "chi": -137.28,
      "glycosidicBond": "anti"
    },
    {
      "index": 19,
      "model": 1,
      "chain": "AC",
      "number": 1001,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AC.U1001",
      "shortName": "U",
      "chi": -141.927,
      "glycosidicBond": "anti"
    },
    {
      "index": 20,
      "model": 1,
      "chain": "AC",
      "number": 1002,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AC.G1002",
      "shortName": "G",
      "chi": -165.93,
      "glycosidicBond": "anti"
    },
    {
      "index": 21,
      "model": 1,
      "chain": "AC",
      "number": 1003,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AC.G1003",
      "shortName": "G",
      "chi": -121.565,
      "glycosidicBond": "anti"
    },
    {
      "index": 22,
      "model": 1,
      "chain": "AC",
      "number": 1004,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AC.G1004",
      "shortName": "G",
      "chi": -156.01,
      "glycosidicBond": "anti"
    },
    {
      "index": 23,
      "model": 1,
      "chain": "AC",
      "number": 1005,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AC.G1005",
      "shortName": "G",
      "chi": -148.101,
      "glycosidicBond": "anti"
    },
    {
      "index": 24,
      "model": 1,
      "chain": "AC",
      "number": 1006,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "AC.U1006",
      "shortName": "U",
      "chi": -137.28,
      "glycosidicBond": "anti"
    },
    {
      "index": 25,
      "model": 1,
      "chain": "B",
      "number": 2001,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.U2001",
      "shortName": "U",
      "chi": -146.462,
      "glycosidicBond": "anti"
    },
    {
      "index": 26,
      "model": 1,
      "chain": "B",
      "number": 2002,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.G2002",
      "shortName": "G",
      "chi": -170.797,
      "glycosidicBond": "anti"
    },
    {
      "index": 27,
      "model": 1,
      "chain": "B",
      "number": 2003,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.G2003",
      "shortName": "G",
      "chi": -117.687,
      "glycosidicBond": "anti"
    },
    {
      "index": 28,
      "model": 1,
      "chain": "B",
      "number": 2004,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.G2004",
      "shortName": "G",
      "chi": -153.886,
      "glycosidicBond": "anti"
    },
    {
      "index": 29,
      "model": 1,
      "chain": "B",
      "number": 2005,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.G2005",
      "shortName": "G",
      "chi": -148.852,
      "glycosidicBond": "anti"
    },
    {
      "index": 30,
      "model": 1,
      "chain": "B",
      "number": 2006,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "B.U2006",
      "shortName": "U",
      "chi": -159.437,
      "glycosidicBond": "anti"
    },
    {
      "index": 37,
      "model": 1,
      "chain": "BA",
      "number": 2001,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BA.U2001",
      "shortName": "U",
      "chi": -146.462,
      "glycosidicBond": "anti"
    },
    {
      "index": 38,
      "model": 1,
      "chain": "BA",
      "number": 2002,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BA.G2002",
      "shortName": "G",
      "chi": -170.797,
      "glycosidicBond": "anti"
    },
    {
      "index": 39,
      "model": 1,
      "chain": "BA",
      "number": 2003,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BA.G2003",
      "shortName": "G",
      "chi": -117.687,
      "glycosidicBond": "anti"
    },
    {
      "index": 40,
      "model": 1,
      "chain": "BA",
      "number": 2004,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BA.G2004",
      "shortName": "G",
      "chi": -153.886,
      "glycosidicBond": "anti"
    },
    {
      "index": 41,
      "model": 1,
      "chain": "BA",
      "number": 2005,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BA.G2005",
      "shortName": "G",
      "chi": -148.852,
      "glycosidicBond": "anti"
    },
    {
      "index": 42,
      "model": 1,
      "chain": "BA",
      "number": 2006,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BA.U2006",
      "shortName": "U",
      "chi": -159.437,
      "glycosidicBond": "anti"
    },
    {
      "index": 43,
      "model": 1,
      "chain": "BB",
      "number": 2001,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BB.U2001",
      "shortName": "U",
      "chi": -146.462,
      "glycosidicBond": "anti"
    },
    {
      "index": 44,
      "model": 1,
      "chain": "BB",
      "number": 2002,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BB.G2002",
      "shortName": "G",
      "chi": -170.797,
      "glycosidicBond": "anti"
    },
    {
      "index": 45,
      "model": 1,
      "chain": "BB",
      "number": 2003,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BB.G2003",
      "shortName": "G",
      "chi": -117.687,
      "glycosidicBond": "anti"
    },
    {
      "index": 46,
      "model": 1,
      "chain": "BB",
      "number": 2004,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BB.G2004",
      "shortName": "G",
      "chi": -153.886,
      "glycosidicBond": "anti"
    },
    {
      "index": 47,
      "model": 1,
      "chain": "BB",
      "number": 2005,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BB.G2005",
      "shortName": "G",
      "chi": -148.852,
      "glycosidicBond": "anti"
    },
    {
      "index": 48,
      "model": 1,
      "chain": "BB",
      "number": 2006,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BB.U2006",
      "shortName": "U",
      "chi": -159.437,
      "glycosidicBond": "anti"
    },
    {
      "index": 31,
      "model": 1,
      "chain": "BC",
      "number": 2001,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BC.U2001",
      "shortName": "U",
      "chi": -146.462,
      "glycosidicBond": "anti"
    },
    {
      "index": 32,
      "model": 1,
      "chain": "BC",
      "number": 2002,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BC.G2002",
      "shortName": "G",
      "chi": -170.797,
      "glycosidicBond": "anti"
    },
    {
      "index": 33,
      "model": 1,
      "chain": "BC",
      "number": 2003,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BC.G2003",
      "shortName": "G",
      "chi": -117.687,
      "glycosidicBond": "anti"
    },
    {
      "index": 34,
      "model": 1,
      "chain": "BC",
      "number": 2004,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BC.G2004",
      "shortName": "G",
      "chi": -153.886,
      "glycosidicBond": "anti"
    },
    {
      "index": 35,
      "model": 1,
      "chain": "BC",
      "number": 2005,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BC.G2005",
      "shortName": "G",
      "chi": -148.852,
      "glycosidicBond": "anti"
    },
    {
      "index": 36,
      "model": 1,
      "chain": "BC",
      "number": 2006,
      "icode": " ",
      "molecule": "RNA",
      "fullName": "BC.U2006",
      "shortName": "U",
      "chi": -159.437,
      "glycosidicBond": "anti"
    }
  ],
  "basePairs": [
    {
      "nt1": "A.U1001",
      "nt2": "B.G2003",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "B.G2003",
      "nt2": "A.U1001",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "A.G1002",
      "nt2": "AB.G1002",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "AB.G1002",
      "nt2": "A.G1002",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "A.G1002",
      "nt2": "AC.G1002",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AC.G1002",
      "nt2": "A.G1002",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.G1003",
      "nt2": "AB.G1003",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "AB.G1003",
      "nt2": "A.G1003",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "A.G1003",
      "nt2": "AC.G1003",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AC.G1003",
      "nt2": "A.G1003",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.G1003",
      "nt2": "B.U2001",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "B.U2001",
      "nt2": "A.G1003",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "A.G1004",
      "nt2": "AB.G1004",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "AB.G1004",
      "nt2": "A.G1004",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "A.G1004",
      "nt2": "AC.G1004",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AC.G1004",
      "nt2": "A.G1004",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.G1005",
      "nt2": "AB.G1005",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "AB.G1005",
      "nt2": "A.G1005",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "A.G1005",
      "nt2": "AC.G1005",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AC.G1005",
      "nt2": "A.G1005",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.U1006",
      "nt2": "AB.U1006",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AB.U1006",
      "nt2": "A.U1006",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "A.U1006",
      "nt2": "AC.U1006",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "AC.U1006",
      "nt2": "A.U1006",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AA.U1001",
      "nt2": "BA.G2003",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "BA.G2003",
      "nt2": "AA.U1001",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "AA.G1002",
      "nt2": "AB.G1002",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AB.G1002",
      "nt2": "AA.G1002",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "AA.G1002",
      "nt2": "AC.G1002",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "AC.G1002",
      "nt2": "AA.G1002",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AA.G1003",
      "nt2": "AB.G1003",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AB.G1003",
      "nt2": "AA.G1003",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "AA.G1003",
      "nt2": "AC.G1003",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "AC.G1003",
      "nt2": "AA.G1003",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AA.G1003",
      "nt2": "BA.U2001",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "BA.U2001",
      "nt2": "AA.G1003",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "AA.G1004",
      "nt2": "AB.G1004",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AB.G1004",
      "nt2": "AA.G1004",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "AA.G1004",
      "nt2": "AC.G1004",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "AC.G1004",
      "nt2": "AA.G1004",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AA.G1005",
      "nt2": "AB.G1005",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AB.G1005",
      "nt2": "AA.G1005",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "AA.G1005",
      "nt2": "AC.G1005",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "AC.G1005",
      "nt2": "AA.G1005",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AA.U1006",
      "nt2": "AB.U1006",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "AB.U1006",
      "nt2": "AA.U1006",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AA.U1006",
      "nt2": "AC.U1006",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "AC.U1006",
      "nt2": "AA.U1006",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "AB.U1001",
      "nt2": "BB.G2003",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "BB.G2003",
      "nt2": "AB.U1001",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "AB.G1003",
      "nt2": "BB.U2001",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "BB.U2001",
      "nt2": "AB.G1003",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "AC.U1001",
      "nt2": "BC.G2003",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "BC.G2003",
      "nt2": "AC.U1001",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "AC.G1003",
      "nt2": "BC.U2001",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "BC.U2001",
      "nt2": "AC.G1003",
      "lw": "cSS",
      "stericity": "cis",
      "edge5": "Sugar",
      "edge3": "Sugar"
    },
    {
      "nt1": "B.G2002",
      "nt2": "BB.G2002",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "BB.G2002",
      "nt2": "B.G2002",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.G2002",
      "nt2": "BC.G2002",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BC.G2002",
      "nt2": "B.G2002",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "B.G2003",
      "nt2": "BB.G2003",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "BB.G2003",
      "nt2": "B.G2003",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.G2003",
      "nt2": "BC.G2003",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BC.G2003",
      "nt2": "B.G2003",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "B.G2004",
      "nt2": "BB.G2004",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "BB.G2004",
      "nt2": "B.G2004",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.G2004",
      "nt2": "BC.G2004",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BC.G2004",
      "nt2": "B.G2004",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "B.G2005",
      "nt2": "BB.G2005",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "BB.G2005",
      "nt2": "B.G2005",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "B.G2005",
      "nt2": "BC.G2005",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BC.G2005",
      "nt2": "B.G2005",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "B.U2006",
      "nt2": "BB.U2006",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BB.U2006",
      "nt2": "B.U2006",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "B.U2006",
      "nt2": "BC.U2006",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "BC.U2006",
      "nt2": "B.U2006",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BA.G2002",
      "nt2": "BB.G2002",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BB.G2002",
      "nt2": "BA.G2002",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "BA.G2002",
      "nt2": "BC.G2002",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "BC.G2002",
      "nt2": "BA.G2002",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BA.G2003",
      "nt2": "BB.G2003",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BB.G2003",
      "nt2": "BA.G2003",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "BA.G2003",
      "nt2": "BC.G2003",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "BC.G2003",
      "nt2": "BA.G2003",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BA.G2004",
      "nt2": "BB.G2004",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BB.G2004",
      "nt2": "BA.G2004",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "BA.G2004",
      "nt2": "BC.G2004",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "BC.G2004",
      "nt2": "BA.G2004",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BA.G2005",
      "nt2": "BB.G2005",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BB.G2005",
      "nt2": "BA.G2005",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "BA.G2005",
      "nt2": "BC.G2005",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "BC.G2005",
      "nt2": "BA.G2005",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BA.U2006",
      "nt2": "BB.U2006",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    },
    {
      "nt1": "BB.U2006",
      "nt2": "BA.U2006",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BA.U2006",
      "nt2": "BC.U2006",
      "lw": "cWH",
      "stericity": "cis",
      "edge5": "Watson-Crick",
      "edge3": "Hoogsteen"
    },
    {
      "nt1": "BC.U2006",
      "nt2": "BA.U2006",
      "lw": "cHW",
      "stericity": "cis",
      "edge5": "Hoogsteen",
      "edge3": "Watson-Crick"
    }
  ],
  "helices": [
    {
      "quadruplexes": [
        {
          "tetrads": [
            {
              "id": "A.U1006-AC.U1006-AA.U1006-AB.U1006",
              "nt1": "A.U1006",
              "nt2": "AC.U1006",
              "nt3": "AA.U1006",
              "nt4": "AB.U1006",
              "onz": "O-",
              "gbaClassification": "VIIIa",
              "planarityDeviation": 1.061,
              "ionsChannel": [
                "Na"
              ],
              "ionsOutside": [
                {
                  "nt": "A.U1006",
                  "ion": "Sr"
                },
                {
                  "nt": "AA.U1006",
                  "ion": "Sr"
                },
                {
                  "nt": "AB.U1006",
                  "ion": "Sr"
                },
                {
                  "nt": "AC.U1006",
                  "ion": "Sr"
                }
              ]
            },
            {
              "id": "A.G1005-AB.G1005-AA.G1005-AC.G1005",
              "nt1": "A.G1005",
              "nt2": "AB.G1005",
              "nt3": "AA.G1005",
              "nt4": "AC.G1005",
              "onz": "O+",
              "gbaClassification": "VIIIa",
              "planarityDeviation": 0.7999999999999972,
              "ionsChannel": [
                "Sr"
              ],
              "ionsOutside": []
            },
            {
              "id": "A.G1004-AB.G1004-AA.G1004-AC.G1004",
              "nt1": "A.G1004",
              "nt2": "AB.G1004",
              "nt3": "AA.G1004",
              "nt4": "AC.G1004",
              "onz": "O+",
              "gbaClassification": "VIIIa",
              "planarityDeviation": 0.4059999999999988,
              "ionsChannel": [],
              "ionsOutside": []
            },
            {
              "id": "A.G1003-AB.G1003-AA.G1003-AC.G1003",
              "nt1": "A.G1003",
              "nt2": "AB.G1003",
              "nt3": "AA.G1003",
              "nt4": "AC.G1003",
              "onz": "O+",
              "gbaClassification": "VIIIa",
              "planarityDeviation": 0.5549999999999997,
              "ionsChannel": [
                "Sr"
              ],
              "ionsOutside": []
            },
            {
              "id": "A.G1002-AB.G1002-AA.G1002-AC.G1002",
              "nt1": "A.G1002",
              "nt2": "AB.G1002",
              "nt3": "AA.G1002",
              "nt4": "AC.G1002",
              "onz": "O+",
              "gbaClassification": "VIIIa",
              "planarityDeviation": 0.541999999999998,
              "ionsChannel": [],
              "ionsOutside": []
            }
          ],
          "onzm": "Op*",
          "loopClassification": "n/a",
          "gbaClassification": [
            "VIII"
          ],
          "tracts": [
            [
              "A.U1006",
              "A.G1005",
              "A.G1004",
              "A.G1003",
              "A.G1002"
            ],
            [
              "AC.U1006",
              "AC.G1005",
              "AC.G1004",
              "AC.G1003",
              "AC.G1002"
            ],
            [
              "AA.U1006",
              "AA.G1005",
              "AA.G1004",
              "AA.G1003",
              "AA.G1002"
            ],
            [
              "AB.U1006",
              "AB.G1005",
              "AB.G1004",
              "AB.G1003",
              "AB.G1002"
            ]
          ],
          "loops": []
        },
        {
          "tetrads": [
            {
              "id": "B.G2002-BC.G2002-BA.G2002-BB.G2002",
              "nt1": "B.G2002",
              "nt2": "BC.G2002",
              "nt3": "BA.G2002",
              "nt4": "BB.G2002",
              "onz": "O+",
              "gbaClassification": "VIIIa",
              "planarityDeviation": 0.6730000000000018,
              "ionsChannel": [],
              "ionsOutside": []
            },
            {
              "id": "B.G2003-BC.G2003-BA.G2003-BB.G2003",
              "nt1": "B.G2003",
              "nt2": "BC.G2003",
              "nt3": "BA.G2003",
              "nt4": "BB.G2003",
              "onz": "O+",
              "gbaClassification": "VIIIa",
              "planarityDeviation": 0.5769999999999982,
              "ionsChannel": [
                "Sr"
              ],
              "ionsOutside": []
            },
            {
              "id": "B.G2004-BC.G2004-BA.G2004-BB.G2004",
              "nt1": "B.G2004",
              "nt2": "BC.G2004",
              "nt3": "BA.G2004",
              "nt4": "BB.G2004",
              "onz": "O+",
              "gbaClassification": "VIIIa",
              "planarityDeviation": 0.2289999999999992,
              "ionsChannel": [
                "Sr"
              ],
              "ionsOutside": []
            },
            {
              "id": "B.G2005-BC.G2005-BA.G2005-BB.G2005",
              "nt1": "B.G2005",
              "nt2": "BC.G2005",
              "nt3": "BA.G2005",
              "nt4": "BB.G2005",
              "onz": "O+",
              "gbaClassification": "VIIIa",
              "planarityDeviation": 0.7810000000000006,
              "ionsChannel": [
                "Na"
              ],
              "ionsOutside": []
            },
            {
              "id": "B.U2006-BB.U2006-BA.U2006-BC.U2006",
              "nt1": "B.U2006",
              "nt2": "BB.U2006",
              "nt3": "BA.U2006",
              "nt4": "BC.U2006",
              "onz": "O-",
              "gbaClassification": "VIIIa",
              "planarityDeviation": 1.5840000000000005,
              "ionsChannel": [
                "Na",
                "Na"
              ],
              "ionsOutside": []
            }
          ],
          "onzm": "Op*",
          "loopClassification": "n/a",
          "gbaClassification": [
            "VIII"
          ],
          "tracts": [
            [
              "B.G2002",
              "B.G2003",
              "B.G2004",
              "B.G2005",
              "B.U2006"
            ],
            [
              "BC.G2002",
              "BC.G2003",
              "BC.G2004",
              "BC.G2005",
              "BC.U2006"
            ],
            [
              "BA.G2002",
              "BA.G2003",
              "BA.G2004",
              "BA.G2005",
              "BA.U2006"
            ],
            [
              "BB.G2002",
              "BB.G2003",
              "BB.G2004",
              "BB.G2005",
              "BB.U2006"
            ]
          ],
          "loops": []
        }
      ],
      "tetradPairs": [
        {
          "tetrad1": "A.U1006-AC.U1006-AA.U1006-AB.U1006",
          "tetrad2": "A.G1005-AB.G1005-AA.G1005-AC.G1005",
          "direction": "parallel",
          "rise": 3.366499999999995,
          "twist": 39.962531742191736
        },
        {
          "tetrad1": "A.G1005-AB.G1005-AA.G1005-AC.G1005",
          "tetrad2": "A.G1004-AB.G1004-AA.G1004-AC.G1004",
          "direction": "parallel",
          "rise": 3.308,
          "twist": 25.89614444631925
        },
        {
          "tetrad1": "A.G1004-AB.G1004-AA.G1004-AC.G1004",
          "tetrad2": "A.G1003-AB.G1003-AA.G1003-AC.G1003",
          "direction": "parallel",
          "rise": 3.3394999999999904,
          "twist": 35.81115298630443
        },
        {
          "tetrad1": "A.G1003-AB.G1003-AA.G1003-AC.G1003",
          "tetrad2": "A.G1002-AB.G1002-AA.G1002-AC.G1002",
          "direction": "parallel",
          "rise": 3.2865000000000073,
          "twist": 27.11515971986807
        },
        {
          "tetrad1": "A.G1002-AB.G1002-AA.G1002-AC.G1002",
          "tetrad2": "B.G2002-BC.G2002-BA.G2002-BB.G2002",
          "direction": "parallel",
          "rise": 3.3694999999999986,
          "twist": 28.993180312675587
        },
        {
          "tetrad1": "B.G2002-BC.G2002-BA.G2002-BB.G2002",
          "tetrad2": "B.G2003-BC.G2003-BA.G2003-BB.G2003",
          "direction": "parallel",
          "rise": 3.371000000000002,
          "twist": 27.410084968596852
        },
        {
          "tetrad1": "B.G2003-BC.G2003-BA.G2003-BB.G2003",
          "tetrad2": "B.G2004-BC.G2004-BA.G2004-BB.G2004",
          "direction": "parallel",
          "rise": 3.3180000000000014,
          "twist": 35.04072146975963
        },
        {
          "tetrad1": "B.G2004-BC.G2004-BA.G2004-BB.G2004",
          "tetrad2": "B.G2005-BC.G2005-BA.G2005-BB.G2005",
          "direction": "parallel",
          "rise": 3.2689999999999984,
          "twist": 25.149997949938147
        },
        {
          "tetrad1": "B.G2005-BC.G2005-BA.G2005-BB.G2005",
          "tetrad2": "B.U2006-BB.U2006-BA.U2006-BC.U2006",
          "direction": "parallel",
          "rise": 7.140500000000001,
          "twist": 43.40609492262336
        }
      ]
    }
  ],
  "dotBracket": {
    "sequence": "UGGGGU-UGGGGU-UGGGGU-UGGGGU-UGGGGU-UGGGGU-UGGGGU-UGGGGU",
    "line1": ".([{<A-.)]}>A-.([{<a-.)]}>a-.([{<A-.)]}>A-.([{<a-.)]}>a",
    "line2": ".([{<A-.([{<a-.)]}>A-.)]}>a-.([{<A-.([{<a-.)]}>A-.)]}>a"
  }
}
```

</details>

# Funding

This research was supported by the National Science Centre, Poland
\[2016/23/B/ST6/03931, 2019/35/B/ST6/03074\] and Mloda Kadra project
\[09/91/SBAD/0684\] from Poznan University of Technology, and carried
out in the European Centre for Bioinformatics and Genomics (Poland). The
authors also acknowledge partial support by the statutory funds of
Poznan University of Technology, Polish Ministry of Science and Higher
Education, and the Institute of Bioorganic Chemistry, PAS within
intramural financing program.

# Bibliography

<div id="refs" class="references csl-bib-body">

1.  Topology-based classification of tetrads and quadruplex
    structures. M. Popenda, J. Miskiewicz, J. Sarzynska, T. Zok, M.
    Szachniuk. *Bioinformatics*. 2020. 36(4):1129–1134.
    doi:[10.1093/bioinformatics/btz738](https://doi.org/10.1093/bioinformatics/btz738)

2.  ElTetrado: A tool for identification and classification of tetrads
    and quadruplexes. T. Zok, M. Popenda, M. Szachniuk. *BMC
    Bioinformatics*. 2020. 21(1):40.
    doi:[10.1186/s12859-020-3385-1](https://doi.org/10.1186/s12859-020-3385-1)

3.  DSSR: An integrated software tool for dissecting the spatial
    structure of RNA. X.-J. Lu, H.J. Bussemaker, W.K. Olson. *Nucleic
    Acids Research*. 2015. 43(21):e142.
    doi:[10.1093/nar/gkv716](https://doi.org/10.1093/nar/gkv716)

4.  R-chie : A web server and R package for visualizing RNA secondary
    structures. D. Lai, J.R. Proctor, J.Y.A. Zhu, I.M. Meyer. *Nucleic
    Acids Research*. 2012. 40(12):e95.
    doi:[10.1093/nar/gks241](https://doi.org/10.1093/nar/gks241)

5.  Geometric nomenclature and classification of RNA base pairs. N.B.
    Leontis, E. Westhof. *RNA*. 2001. 7(4):499–512.
    doi:[10.1017/S1355838201002515](https://doi.org/10.1017/S1355838201002515)

</div>
