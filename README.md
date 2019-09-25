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
please refer to Zok, Popenda and Szachniuk (2019) and Popenda *et al.*
(2019)

# Dependencies

The project is written in Python 3.6+ and requires
[Biopython](https://biopython.org/) and [NumPy](https://numpy.org/) in
to run. These can be installed with the following command:

    pip install -r requirements.txt

ElTetrado depends on DSSR (Lu, Bussemaker and Olson, 2015) in terms of
detection of base pairing and stacking. The binary `x3dna-dssr` can be
[downloaded](http://forum.x3dna.org/site-announcements/download-instructions/)
and put in the same directory as `eltetrado` which will execute it
during analysis. Alternatively, one can preprocess the 3D data with
`x3dna-dssr --json` and provide the JSON result as an input to ElTetrado
(see Usage section below).

Visualization is created by `R` 3.6+ script which uses
[R4RNA](https://www.e-rna.org/r-chie/) (Lai *et al.*, 2012) library. The
dependency will be automatically installed if not present.

# Usage

ElTetrado is a command line application, which requires to be provided
with:

  - either `--dssr-json` and the path to JSON generated with `x3dna-dssr
    --json` (fast, but quadruplex parameters like `rise` or `twist` are
    not calculated)
  - or `--pdb` and the path to PDB or PDBx/mmCIF file (slow, because the
    execution time is a sum of ElTetrado and DSSR times)
  - or both `--pdb` and `--dssr-json` (recommended, all analyses are
    made and only ElTetrado is executed while DSSR results are read from
    the file)

By default, ElTetrado outputs textual results on the standard output. A
JSON version of the output can be obtained with `--output` switch
followed by a path where the file is supposed to be created.

ElTetrado prepares visualization of the whole structure and of each
N4-helices, quadruplexes and tetrads. This can be supplemented with
canonical base pairs visualization when `--complete-2d` is set. All
color settings are located in the first several lines of the `quadraw.R`
file, you can easily change them without knowledge of R langauge. If you
want ElTetrado to not visualize anything, pass `--no-image` switch to
it.

    usage: eltetrado [-h] [--pdb PDB] [--dssr-json DSSR_JSON] [--output OUTPUT]
                     [--stacking-mismatch STACKING_MISMATCH]
                     [--relaxed-stem-definition] [--strict] [--no-reorder]
                     [--complete-2d] [--no-image] [--version]
    
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
      --relaxed-stem-definition
                            when set, two sequentially close tetrades will be
                            considered a stem regardless of their stacking
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

# Examples

## 1MY9: Solution structure of a K+ cation stabilized dimeric RNA quadruplex containing two G:G(:A):G:G(:A) hexads, G:G:G:G tetrads and UUUU loops

![](1my9.png)

    $ curl ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/my/1my9.cif.gz | gzip -d > 1my9.cif
    
    $ ./eltetrado --pdb 1my9.cif
    
    Chain order: A, B
    n4-helix with 4 tetrads
      Op quadruplex with 2 tetrads
        A.G1 A.G4 A.G10 A.G13 cWH-cWH-cWH-cWH O+ planarity=0.63
          direction=parallel rise=4.12 twist=28.66
        A.G2 A.G5 A.G11 A.G14 cWH-cWH-cWH-cWH O+ planarity=0.5
      Op quadruplex with 2 tetrads
        B.G15 B.G18 B.G24 B.G27 cWH-cWH-cWH-cWH O+ planarity=0.26
          direction=parallel rise=4.21 twist=33.45
        B.G16 B.G19 B.G25 B.G28 cWH-cWH-cWH-cWH O+ planarity=0.31
    
    GGAGGUUUUGGAGG-GGAGGUUUUGGAGG
    ([.)]....([.)]-([.)]....([.)]
    ([.([....)].)]-([.([....)].)]
    
    Plot: 1my9-str.pdf
    
    Plot: 1my9-h1.pdf
    
    Plot: 1my9-h1-q1.pdf
    
    Plot: 1my9-h1-q1-t1.pdf
    
    Plot: 1my9-h1-q1-t2.pdf
    
    Plot: 1my9-h1-q2.pdf
    
    Plot: 1my9-h1-q2-t1.pdf
    
    Plot: 1my9-h1-q2-t2.pdf

## 4RJ1: Structural variations and solvent structure of UGGGGU quadruplexes stabilized by Sr2+ ions

![](4rj1.png)

    $ curl https://www.ebi.ac.uk/pdbe/static/entry/download/4rj1-assembly-1.cif.gz | gzip -d > 4rj1-1.cif
    
    $ ./eltetrado --pdb 4rj1-1.cif
    
    Chain order: A, AB, AA, AC, B, BC, BA, BB
    n4-helix with 9 tetrads
      Op quadruplex with 5 tetrads
        A.G1002 AB.G1002 AA.G1002 AC.G1002 cWH-cWH-cWH-cWH O+ planarity=0.54
          direction=parallel rise=3.29 twist=27.12
        A.G1003 AB.G1003 AA.G1003 AC.G1003 cWH-cWH-cWH-cWH O+ planarity=0.56
          direction=parallel rise=3.34 twist=35.81
        A.G1004 AB.G1004 AA.G1004 AC.G1004 cWH-cWH-cWH-cWH O+ planarity=0.41
          direction=parallel rise=3.31 twist=25.9
        A.G1005 AB.G1005 AA.G1005 AC.G1005 cWH-cWH-cWH-cWH O+ planarity=0.8
          direction=parallel rise=3.37 twist=39.96
        A.U1006 AC.U1006 AA.U1006 AB.U1006 cWH-cWH-cWH-cWH O- planarity=1.06
      Op quadruplex with 4 tetrads
        B.G2002 BC.G2002 BA.G2002 BB.G2002 cWH-cWH-cWH-cWH O+ planarity=0.67
          direction=parallel rise=3.37 twist=27.41
        B.G2003 BC.G2003 BA.G2003 BB.G2003 cWH-cWH-cWH-cWH O+ planarity=0.58
          direction=parallel rise=3.32 twist=35.04
        B.G2004 BC.G2004 BA.G2004 BB.G2004 cWH-cWH-cWH-cWH O+ planarity=0.23
          direction=parallel rise=3.27 twist=25.15
        B.G2005 BC.G2005 BA.G2005 BB.G2005 cWH-cWH-cWH-cWH O+ planarity=0.78
    single tetrad without stacking
      single tetrad
        B.U2006 BB.U2006 BA.U2006 BC.U2006 cWH-cWH-cWH-cWH O- planarity=1.58
    
    UGGGGU-UGGGGU-UGGGGU-UGGGGU-UGGGGU-UGGGGU-UGGGGU-UGGGGU
    .([{<A-.)]}>A-.([{<a-.)]}>a-.([{<A-.)]}>A-.([{<a-.)]}>a
    .([{<A-.([{<a-.)]}>A-.)]}>a-.([{<A-.([{<a-.)]}>A-.)]}>a
    
    Plot: 4rj1-1-str.pdf
    
    Plot: 4rj1-1-h1.pdf
    
    Plot: 4rj1-1-h1-q1.pdf
    
    Plot: 4rj1-1-h1-q1-t1.pdf
    
    Plot: 4rj1-1-h1-q1-t2.pdf
    
    Plot: 4rj1-1-h1-q1-t3.pdf
    
    Plot: 4rj1-1-h1-q1-t4.pdf
    
    Plot: 4rj1-1-h1-q1-t5.pdf
    
    Plot: 4rj1-1-h1-q2.pdf
    
    Plot: 4rj1-1-h1-q2-t1.pdf
    
    Plot: 4rj1-1-h1-q2-t2.pdf
    
    Plot: 4rj1-1-h1-q2-t3.pdf
    
    Plot: 4rj1-1-h1-q2-t4.pdf
    
    Plot: 4rj1-1-h2.pdf
    
    Plot: 4rj1-1-h2-q1.pdf
    
    Plot: 4rj1-1-h2-q1-t1.pdf

# Bibliography

<div id="refs" class="references">

<div id="ref-Lai2012">

Lai, D., Proctor, J. R., Zhu, J. Y. A. and Meyer, I. M. (2012) ‘R-CHIE:
A web server and R package for visualizing RNA secondary structures’,
*Nucleic Acids Research*, 40(12), p. e95. doi:
[10.1093/nar/gks241](https://doi.org/10.1093/nar/gks241).

</div>

<div id="ref-Leontis2001">

Leontis, N. B. and Westhof, E. (2001) ‘Geometric nomenclature and
classification of RNA base pairs’, *RNA*, 7(4), pp. 499–512. doi:
[10.1017/s1355838201002515](https://doi.org/10.1017/s1355838201002515).

</div>

<div id="ref-Lu2015">

Lu, X. J., Bussemaker, H. J. and Olson, W. K. (2015) ‘DSSR: An
integrated software tool for dissecting the spatial structure of RNA’,
*Nucleic Acids Research*, 43(21), p. gkv716. doi:
[10.1093/nar/gkv716](https://doi.org/10.1093/nar/gkv716).

</div>

<div id="ref-Popenda2019">

Popenda, M., Miskiewicz, J., Sarzynska, J., Zok, T. and Szachniuk, M.
(2019) ‘Topology-based classification of tetrads and quadruplex
structures’, *Bioinformatics*.

</div>

<div id="ref-Zok2019">

Zok, T., Popenda, M. and Szachniuk, M. (2019) ‘ElTetrado: a tool for
identification and ONZ classification of tetrads and quadruplexes’, *in
submission*.

</div>

</div>
