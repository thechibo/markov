# Maize Opaque-2 Genomic DNA Sequence

The genomic nucleotide sequence of the maize (*Zea mays*) *opaque
endosperm 2* gene, NCBI Gene ID 542375, also known as *o2* or
*Opaque-2*.

## Usage

``` r
opaque2
```

## Format

A character vector of length 3,158. Each element is one nucleotide,
encoded as `"A"`, `"C"`, `"G"`, or `"T"`, in genomic order. Counts are
798 A, 755 C, 765 G, and 840 T. There are no missing or ambiguous bases.

## Source

Reference Sequence (RefSeq) database, National Center for Biotechnology
Information (NCBI), U.S. National Library of Medicine, National
Institutes of Health. See the [RefSeq chromosome 7 record
(NC_050102.1)](https://www.ncbi.nlm.nih.gov/nuccore/NC_050102.1) and the
[NCBI Gene record (542375)](https://www.ncbi.nlm.nih.gov/gene/542375).

## Details

The gene is annotated as `LOC542375` on chromosome 7 of the maize
assembly `Zm-B73-REFERENCE-NAM-5.0` (`GCF_902167145.1`) in NCBI Zea mays
Annotation Release 103 (2020-08-09). The sequence is the plus strand of
the inclusive interval 10,951,530–10,954,687 on RefSeq `NC_050102.1`.
Element \\k\\ corresponds to the one-based chromosome coordinate
\\10951529 + k\\. The genomic interval includes introns.

## Retrieving the sequence

The dataset can be retrieved as FASTA from NCBI RefSeq via EFetch using
the accession, coordinates, and strand specified above.

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

if (FALSE) { # \dontrun{
refseq_url <- paste0(
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
  "db=nuccore&id=NC_050102.1&strand=1",
  "&seq_start=10951530&seq_stop=10954687&rettype=fasta&retmode=text"
)
fasta <- readLines(refseq_url)
sequence <- paste0(fasta[-1L], collapse = "")
bases <- strsplit(sequence, "", fixed = TRUE)[[1L]]
stopifnot(identical(bases, opaque2))
} # }
```
