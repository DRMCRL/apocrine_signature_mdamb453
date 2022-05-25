library(memes)
library(universalmotif)
library(tidyverse)
library(plyranges)
library(extraChIPs)
library(rtracklayer)
library(magrittr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(scales)

options(meme_bin = "/opt/meme/bin/")
check_meme_install()

hg19 <- BSgenome.Hsapiens.UCSC.hg19
sq19 <- seqinfo(hg19)

## Load gene annotations
all_gr <- here::here("data", "annotations", "all_gr.rds") %>%
  read_rds()

## Load the RNA-Seq data
counts <- here::here("data", "rnaseq", "counts.out.gz") %>%
  read_tsv(skip = 1) %>%
  dplyr::select(Geneid, ends_with("bam")) %>%
  dplyr::rename(gene_id = Geneid)
detected <- all_gr$gene %>%
  as_tibble() %>%
  distinct(gene_id, gene_name) %>%
  left_join(counts) %>%
  pivot_longer(
    cols = ends_with("bam"), names_to = "sample", values_to = "counts"
  ) %>%
  mutate(detected = counts > 0) %>%
  group_by(gene_id, gene_name) %>%
  summarise(detected = mean(detected) > 0.25, .groups = "drop") %>%
  dplyr::filter(detected) %>%
  dplyr::select(starts_with("gene"))

## Load H3K27ac features
features <- here::here("data", "h3k27ac") %>%
  list.files(full.names = TRUE, pattern = "bed$") %>%
  sapply(import.bed, seqinfo = sq19) %>%
  lapply(granges) %>%
  setNames(basename(names(.))) %>%
  setNames(str_remove_all(names(.), "s.bed")) %>%
  GRangesList() %>%
  unlist() %>%
  names_to_column("feature") %>%
  sort()

## Load the database
meme_db <- here::here("data", "external", "JASPAR2022_CORE_Hsapiens.meme") %>%
  read_meme() %>%
  to_df() %>%
  as_tibble() %>%
  mutate(
    across(everything(), vctrs::vec_proxy),
    altname = str_remove(altname, "MA[0-9]+\\.[0-9]*\\."),
    organism = "Homo sapiens"
  )
nrow(meme_db) # [1] 727

## Reduce down to only the expressed genes
meme_db_detected <- meme_db %>%
  mutate(
    gene_name = str_split(altname, pattern = "::", n = 2)
  ) %>%
  dplyr::filter(
    vapply(
      gene_name,
      function(x) all(x %in% detected$gene_name),
      logical(1)
    )
  )
nrow(meme_db_detected) # [1] 465

## Check similarities
db_cors <- meme_db_detected %>%
  to_list() %>%
  compare_motifs()
db_cors[is.na(db_cors)] <- 0
motif_cluster <- (1 - db_cors) %>%
  as.dist() %>%
  hclust() %>%
  cutree(h = 0.05)


## Load the ranges
dht_peaks <- here::here("data", "peaks") %>%
  list.files(recursive = TRUE, pattern = "oracle", full.names = TRUE) %>%
  sapply(read_rds, simplify = FALSE) %>%
  lapply(function(x) x[["DHT"]]) %>%
  lapply(setNames, nm = c()) %>%
  setNames(str_extract_all(names(.), "AR|FOXA1|GATA3|TFAP2B"))
targets <- names(dht_peaks)
## Shift from GRCh37 to hg19
seqinfo(dht_peaks$AR) <- sq19
seqinfo(dht_peaks$FOXA1) <- sq19
seqinfo(dht_peaks$GATA3) <- sq19
seqinfo(dht_peaks$TFAP2B) <- sq19
dht_consensus <- dht_peaks %>%
  lapply(granges) %>%
  GRangesList() %>%
  unlist() %>%
  reduce() %>%
  mutate(
    AR = overlapsAny(., dht_peaks$AR),
    FOXA1 = overlapsAny(., dht_peaks$FOXA1),
    GATA3 = overlapsAny(., dht_peaks$GATA3),
    TFAP2B = overlapsAny(., dht_peaks$TFAP2B),
    promoter = bestOverlap(., features, var = "feature", missing = "None") == "promoter",
    enhancer = bestOverlap(., features, var = "feature", missing = "None") == "enhancer",
    H3K27ac = promoter | enhancer
  )


## Get the sequences
## These should realistically be centred by each TF and the width set at 100-200
all4_by_h3k27ac <- dht_consensus %>%
  filter(!!!syms(targets)) %>%
  resize(width = 200, fix = 'center') %>%
  split(.$H3K27ac) %>%
  setNames(c("no_h3k27ac", "h3k27ac")) %>%
  get_sequence(hg19) %>%
  as.list()
all4_by_h3k27ac$no_targets <- features %>%
  filter_by_non_overlaps(dht_consensus) %>%
  resize(width = 200, fix = 'center') %>%
  get_sequence(hg19) %>%
  as("BStringSet")
all4_by_h3k27ac <- BStringSetList(all4_by_h3k27ac)
sapply(all4_by_h3k27ac, length)

## Motif Discovery using DREME is DEPRECATED
dreme_all4 <- all4_by_h3k27ac %>%
  .[str_ends(names(.), "h3k27ac")] %>%
  runDreme(control = "no_h3k27ac")
## Try STREME instead
streme_all4 <- runStreme(
  input = all4_by_h3k27ac$h3k27ac,
  control = all4_by_h3k27ac$no_h3k27ac
)


## Compare to the database
tomtom_all4 <- streme_all4 %>%
  to_list() %>%
  runTomTom(
    database = to_list(meme_db_detected),
    # evalue = TRUE,
    thresh = 0.1
  )

## RunAME. Do we need to command line, or a better install?
ame_all4_by_h3k27ac <- all4_by_h3k27ac %>%
  runAme(
    control = "no_targets",
    database = to_list(meme_db_detected),
  )

ame_all4_by_h3k27ac %>%
  bind_rows(.id = "group") %>%
  dplyr::filter(motif_alt_id %in% targets) %>%
  mutate(
    group = str_replace_all(group, "_", " ") %>%
      str_to_title() %>%
      str_replace_all("H3k", "H3K")
  ) %>%
  plot_ame_heatmap(
    id = motif_id,
    value = tp_percent,
    group = group
  )  +
  geom_label(aes(label = percent(tp_percent/100, accuracy = 0.1))) +
  facet_grid(.~motif_alt_id, scales = "free", space = "free") +
  scale_x_discrete(expand = expansion(c(0,0))) +
  scale_y_discrete(expand = expansion(c(0,0))) +
  labs(x = c(), y = c(), fill = "TP (%)") +
  ggtitle("Estimated % True Positive Sites (AME)") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

## It might be worth plotting the position of matches for target motifs in the merged
## peaks, against the summits themselves
