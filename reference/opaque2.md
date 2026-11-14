# Maize Opaque-2 Genomic DNA Sequence

The genomic nucleotide sequence of the maize (*Zea mays*) *opaque
endosperm 2* gene, NCBI Gene ID 542375, also known as *o2* or
*Opaque-2*. This sequence provides an example of observed categorical
data for fitting a Markov chain.

## Usage

``` r
opaque2
```

## Format

A character vector of length 3,158. Each element is one nucleotide,
encoded as `"A"`, `"C"`, `"G"`, or `"T"`, in genomic order. Counts are
798 A, 755 C, 765 G, and 840 T. There are no missing or ambiguous bases.

## Source

National Center for Biotechnology Information (NCBI), [Gene
542375](https://www.ncbi.nlm.nih.gov/gene/542375) and [RefSeq
NC_050102.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_050102.1). NCBI
places no restrictions on use or distribution of its molecular data; see
its [data usage policies and
disclaimers](https://www.ncbi.nlm.nih.gov/home/about/policies/),
including the qualification concerning rights claimed by original
submitters.

## Details

The source snapshot identifies the gene as `LOC542375` on chromosome 7,
plus strand, of assembly `Zm-B73-REFERENCE-NAM-5.0` (`GCF_902167145.1`),
NCBI Zea mays Annotation Release 103 (2020-08-09). The sequence is the
inclusive interval 10,951,530–10,954,687 on RefSeq `NC_050102.1`.
Element \\k\\ corresponds to the one-based chromosome coordinate
\\10951529 + k\\. This is a genomic gene interval, including introns,
rather than a spliced transcript or a coding-sequence-only record.

The sequence was taken from the `LOC542375.csv` file prepared by
`exe/LOC542375.R` in the MarkovExplorer project. The CSV nucleotide
column was verified against the accompanying `gene.fna`; its sequential
index column is represented by the vector positions. No bases were
removed, reordered, or recoded. The original download date was not
recorded; identifiers and coordinates here describe the archived
snapshot.

Positions along DNA serve as the sequence index in modelling examples. A
first-order, time-homogeneous Markov model is an illustrative
assumption, not a verified biological property of these data. This is
one sequence, not a collection of independent trajectories for
[`test_const()`](test_const.md).

## See also

[`fit()`](fit.md), [`test_prob()`](test_prob.md)

## Examples

``` r
data(opaque2)
head(opaque2)
#> [1] "A" "A" "T" "C" "C" "A"
table(opaque2)
#> opaque2
#>   A   C   G   T 
#> 798 755 765 840 
fit(opaque2)
#> Markov Chain
#>   Number of states: 4 
#> 
#>   Initial distribution (a):
#>     A T C G 
#>     1 0 0 0 
#> 
#>   Transition matrix (P):
#>               A         T         C         G
#>     A 0.2446675 0.2609787 0.2308657 0.2634881
#>     T 0.2321429 0.2857143 0.2214286 0.2607143
#>     C 0.2860927 0.2860927 0.2264901 0.2013245
#>     G 0.2496732 0.2300654 0.2797386 0.2405229
```
