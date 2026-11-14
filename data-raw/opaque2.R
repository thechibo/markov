# Run from the package root. No downloads or additional packages are required.
# Source snapshots were copied from MarkovExplorer/data; exe/LOC542375.R
# created the CSV by splitting the first sequence in gene.fna into bases.
# The original download date was not recorded.
source_dir <- file.path("data-raw", "opaque2")
csv <- read.csv(file.path(source_dir, "LOC542375.csv"),
                stringsAsFactors = FALSE)
fasta <- readLines(file.path(source_dir, "gene.fna"), warn = FALSE)
stopifnot(
  identical(names(csv), c("i", "nucl")),
  identical(csv$i, seq_len(3158L)),
  sum(startsWith(fasta, ">")) == 1L,
  identical(fasta[1L], paste0(
    ">NC_050102.1:10951530-10954687 LOC542375 ",
    "[organism=Zea mays] [GeneID=542375] [chromosome=7]"))
)
opaque2 <- csv$nucl
stopifnot(
  identical(opaque2, strsplit(paste(fasta[-1L], collapse = ""), "",
                             fixed = TRUE)[[1L]]),
  !anyNA(opaque2),
  all(opaque2 %in% c("A", "C", "G", "T"))
)
dir.create("data", showWarnings = FALSE)
save(opaque2, file = file.path("data", "opaque2.rda"),
     compress = "xz", version = 2)
