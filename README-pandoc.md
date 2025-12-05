![](logo.svg)

# Project description

This is an application to analyze base pairing patterns of DNA/RNA 3D structures to find and classify tetrads and quadruplexes. ElTetrado assigns tetrads to one of the ONZ classes (O, N, Z) alongside with the directionality of the tetrad (+/-) determined by the bonds between bases and their non-canonical interactions. The interactions follow Leontis/Westhof classification [@Leontis2001]. Watson-Crick (W) edge of first base in the tetrad structure exposed to the Hoogsteen (H) edge of the next nucleobase from the same tetrad sets the tetrad directionality, clockwise (+) or anticlockwise (-). For more details, please refer to @Zok2020 and @Popenda2020

# Installation

This project uses [Poetry](https://python-poetry.org/) for dependency management.

To install the package, run:

```bash
poetry install
```

# Dependencies

The project is written in Python 3.12+ and requires [mmcif](https://pypi.org/project/mmcif/), [orjson](https://github.com/ijl/orjson), [NumPy](https://numpy.org/) and [rnapolis](https://github.com/tzok/rnapolis-py).

Visualization is created by `R` 3.6+ script which uses [R4RNA](https://www.e-rna.org/r-chie/) [@Lai2012] library. The dependency will be automatically installed if not present.

Base pairs and stacking interactions are identified by [RNApolis](https://github.com/tzok/rnapolis-py).

# Usage

ElTetrado is a command line application, which requires to be provided with `--input` and a path to a PDB or PDBx/mmCIF file.

By default, ElTetrado outputs textual results on the standard output. A JSON version of the output can be obtained with `--output` switch followed by a path where the file is supposed to be created.

ElTetrado prepares visualization of the whole structure and of each N4-helices, quadruplexes and tetrads. This can be supplemented with canonical base pairs visualization when `--complete-2d` is set. All color settings are located in the first several lines of the `quadraw.R` file, you can easily change them without knowledge of R language. If you want ElTetrado to not visualize anything, pass `--no-image` switch to it.

```
!include assets/help.txt
```

# Chains reorder

ElTetrado keeps a global and unique 5'-3' index for every nucleotide which is independent from residue numbers. For example, if a structure has chain M with 60 nucleotides and chain N with 15 nucleotides, then ElTetrado will keep index between 0 and 74 which uniquely identifies every nucleotide. Initially, ElTetrado assigns this indices according to the order of chains in the input file. Therefore, if M preceded N then nucleotides in M will be indexed from 0 to 59 and in N from 60 to 74. Otherwise, nucleotides in N will be indexed from 0 to 14 and in M from 15 to 74.

When `--no-reorder` is present, this initial assignment is used. Otherwise, ElTetrado exhaustively checks all permutations of chains' orders. Every permutation check induces recalculation of the global and unique 5'-3' index and in effect it changes ONZ classification of tetrads.

ElTetrado keeps a table of tetrad classification scores according to these rules:

- Type preference: `O` > `N` > `Z`
- Direction preference: `+` > `-`

The table keeps low values for preferred classes i.e. `O+` is 0, `O-` is 1 and so on up to `Z-` with score 5. For every permutation of chain orders, ElTetrado computes sum of scores for tetrads classification induced by 5'-3' indexing. We select permutation with the minimum value.

# Examples

## 2HY9: Human telomere DNA quadruplex structure in K+ solution hybrid-1 form

![](2hy9.png)

```
$ curl ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/my/2hy9.cif.gz | gzip -d > 2hy9.cif
$ ./eltetrado --input 2hy9.cif --output 2hy9.json
```

```
!include assets/2hy9.txt
```

<details>
<summary>Click to see the output JSON</summary>

```json
!include assets/2hy9.json
```

</details>

## 4RJ1: Structural variations and solvent structure of UGGGGU quadruplexes stabilized by Sr2+ ions

![](4rj1.png)

```
$ curl https://www.ebi.ac.uk/pdbe/static/entry/download/4rj1-assembly-1.cif.gz | gzip -d > 4rj1-1.cif
$ ./eltetrado --input 4rj1-1.cif --output 4rj1-1.json
```

```
!include assets/4rj1-1.txt
```

<details>
<summary>Click to see the output JSON</summary>

```json
!include assets/4rj1-1.json
```

</details>
# Funding

This research was supported by the National Science Centre, Poland [2016/23/B/ST6/03931, 2019/35/B/ST6/03074] and Mloda Kadra project [09/91/SBAD/0684] from Poznan University of Technology, and carried out in the European Centre for Bioinformatics and Genomics (Poland). The authors also acknowledge partial support by the statutory funds of Poznan University of Technology, Polish Ministry of Science and Higher Education, and the Institute of Bioorganic Chemistry, PAS within intramural financing program.

# Bibliography
