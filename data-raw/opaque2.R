# Run from the package root: Rscript data-raw/opaque2.R
# Download the plus-strand gene interval from NCBI RefSeq via EFetch.
refseq_url <- paste0(
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
  "db=nuccore&id=NC_050102.1&strand=1",
  "&seq_start=10951530&seq_stop=10954687&rettype=fasta&retmode=text"
)
fasta <- readLines(refseq_url)
sequence <- paste0(fasta[-1L], collapse = "")
opaque2 <- strsplit(sequence, "", fixed = TRUE)[[1L]]
stopifnot(
  startsWith(fasta[1L], ">NC_050102.1:10951530-10954687 "),
  length(opaque2) == 3158L,
  all(opaque2 %in% c("A", "C", "G", "T"))
)
save(opaque2, file = file.path("data", "opaque2.rda"),
     compress = "xz", version = 2)
