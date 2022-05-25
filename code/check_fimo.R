# FIMO was run on all sites containing all 4 targets, which also overlap H3K27ac
# using: fimo  --oc output/meme/all4_fimo  data/external/all4.meme  output/all4_h3k27ac.fa

library(memes)
library(rtracklayer)
library(BiocParallel)
library(universalmotif)
library(TFBSTools)
register(MulticoreParam(8))

fimo_matches <- import.gff(
  here::here("output", "meme", "all4_fimo", "fimo.gff")
)
fimo_matches %>%
  as_tibble() %>%
  dplyr::select(-source, -type, -phase) %>%
  unnest(Alias) %>%
  mutate(
    seq_range = str_extract(range, "chr[0-9XY]*:[0-9]+-[0-9]+")
  ) %>%
  group_by(Alias) %>%
  summarise(
    n = length(unique(seq_range)),
    .groups = "drop"
  ) %>%
  mutate(`%` = n / length(
    filter(dht_consensus, AR, FOXA1, GATA3, TFAP2B, H3K27ac)
  ))

fimo_matches %>%
  as_tibble() %>%
  dplyr::select(-source, -type, -phase) %>%
  unnest(Alias) %>%
  mutate(
    seq_range = str_extract(range, "chr[0-9XY]*:[0-9]+-[0-9]+")
  ) %>%
  dplyr::filter(str_detect(Alias, "GATA3")) %>%
  group_by(seq_range) %>%
  summarise(
    AR = sum(str_detect(Alias, "GATA3")),
    best = min(pvalue),
    .groups= "drop"
  ) %>%
  arrange(desc(best))

all4_motifs <- read_meme("data/external/all4.meme")
all4_h3k <- getSeq(
  hg19,
  consenus_hg19 %>%
    filter(AR, FOXA1, GATA3, TFAP2B, H3K27ac) %>%
    resize(width = 500, fix = 'center') %>%
    setNames(as.character(.))
)
no_target_h3k <- getSeq(
  hg19,
  features_hg19 %>%
    filter_by_non_overlaps(consenus_hg19) %>%
    resize(width = 500, fix = 'center') %>%
    setNames(as.character(.))
)



min_score <- '80%'
checks <- bplapply(
  all4_motifs,
  function(x) {
    tibble(
      name = x@name,
      altname = x@altname,
      range = names(all4_h3k),
      n = vapply(
        all4_h3k,
        function(seq) {
          fwd <- countPWM(reverseComplement(x@motif), seq, min.score = min_score)
          rev <- countPWM(x@motif, seq, min.score = min_score)
          fwd + rev
        },
        integer(1)
      ),
      any = n > 0
    )
  }
) %>%
  bind_rows()

control_checks <- bplapply(
  all4_motifs,
  function(x) {
    tibble(
      name = x@name,
      altname = x@altname,
      range = names(no_target_h3k),
      n = vapply(
        no_target_h3k,
        function(seq) {
          fwd <- countPWM(reverseComplement(x@motif), seq, min.score = min_score)
          rev <- countPWM(x@motif, seq, min.score = min_score)
          fwd + rev
        },
        integer(1)
      ),
      any = n > 0
    )
  }
) %>%
  bind_rows()

## THis is much better!!!
## Just scan for presence/absence then run Fisher's Exact Test
bind_rows(
  mutate(checks, dataset = "test"),
  mutate(control_checks, dataset = "control")
) %>%
  group_by(name, altname, dataset) %>%
  summarise(
    detected = sum(any),
    not_detected = sum(!any)
  ) %>%
  split(f = .$name) %>%
  lapply(
    function(x) {
      ft <- fisher.test(cbind(x$not_detected, x$detected))
      mutate(distinct(x, name, altname), or = ft$estimate, p = ft$p.value)
    }
  ) %>%
  bind_rows() %>%
  mutate(fdr = p.adjust(p, "fdr"))

## Try pushing up to a few more motifs
db <- read_meme("data/external/JASPAR2022_CORE.meme")[1:50]
db_test <- bplapply(
  db,
  function(x) {
    tibble(
      name = x@name,
      # altname = x@altname, HOCOMO desn't have any altname field
      range = names(all4_h3k),
      n = vapply(
        all4_h3k,
        function(seq) {
          fwd <- countPWM(reverseComplement(x@motif), seq, min.score = min_score)
          rev <- countPWM(x@motif, seq, min.score = min_score)
          fwd + rev
        },
        integer(1)
      ),
      any = n > 0
    )
  }
) %>%
  bind_rows()
db_tbl <- db %>%
  lapply(
    function(x) {
      tibble(
        name = x@name,
        altname = x@altname
      )
    }
  ) %>%
  bind_rows()

db_control <- bplapply(
  db,
  function(x) {
    tibble(
      name = x@name,
      # altname = x@altname,
      range = names(no_target_h3k),
      n = vapply(
        no_target_h3k,
        function(seq) {
          fwd <- countPWM(reverseComplement(x@motif), seq, min.score = min_score)
          rev <- countPWM(x@motif, seq, min.score = min_score)
          fwd + rev
        },
        integer(1)
      ),
      any = n > 0
    )
  }
) %>%
  bind_rows()

bind_rows(
  mutate(db_test, dataset = "test"),
  mutate(db_control, dataset = "control")
) %>%
  group_by(name, dataset) %>%
  summarise(
    detected = sum(any),
    not_detected = sum(!any)
  ) %>%
  split(f = .$name) %>%
  lapply(
    function(x) {
      ft <- fisher.test(cbind(x$not_detected, x$detected))
      mutate(
        distinct(x, name),
        n_test = dplyr::filter(x, dataset == "test")$detected,
        n_control = dplyr::filter(x, dataset == "control")$ detected,
        prop_test = dplyr::filter(x, dataset == "test") %>%
          mutate(prop = detected / (detected + not_detected)) %>%
          pull("prop"),
        prop_control = dplyr::filter(x, dataset == "control") %>%
          mutate(prop = detected / (detected + not_detected)) %>%
          pull("prop") ,
        or = ft$estimate, p = ft$p.value)
    }
  ) %>%
  bind_rows() %>%
  mutate(fdr = p.adjust(p, "fdr")) %>%
  arrange(desc(or)) %>%
  left_join(db_tbl)

bind_rows(
  mutate(db_test, dataset = "test"),
  mutate(db_control, dataset = "control")
) %>%
  group_by(name) %>%
  summarise(
    glm = list(glm(n ~ dataset, family = poisson()))
  ) %>%
  mutate(
    summary = lapply(glm, broom::tidy)
  ) %>%
  unnest(summary) %>%
  dplyr::filter(term == "datasettest") %>%
  mutate(
    rank = - sign(estimate) * log10(p.value),
    p_adj = p.adjust(p.value, "bonf")
  ) %>%
  dplyr::select(-std.error) %>%
  arrange(desc(rank)) %>%
  left_join(db_tbl)


## Check similarities
db <- read_meme("data/external/JASPAR2022_CORE.meme")
cor <- db %>%
  compare_motifs()
cor[is.na(cor)] <- 0
%>%
  set_rownames(
    vapply(seq_along(db), function(x) db[[x]]@altname, character(1))
  ) %>%
  set_colnames(rownames(.))

corrplot(cor, type = "lower", diag = FALSE, order = "hclust")

## Now convert to a distance matrix & square
cor %>%
  raise_to_power(5) %>%
  pheatmap::pheatmap(
    # cutree_rows = 6,
    # cutree_cols = 6
  )

## We know that the lower corner could be merged
motif_to_merge <- c("Foxq1", "FOXF2", "FOXD1", "SOX9", "SRY")
db %>%
  .[vapply(., function(x) x@altname %in% motif_to_merge, logical(1))] %>%
  merge_motifs(normalise.scores = TRUE, min.position.ic = 1) %>%
  view_motifs(names.pos = "right")


cl <- (1 - cor) %>%
  # raise_to_power(10) %>%
  as.dist() %>%
  # cluster::diana() %>%
  # cluster::agnes() %>%
  hclust(method = "ward.D2") %>%
  # plot()
  # abline(h = 0.05)
  cutree(h = 0.05) %>%
  split(f = as.numeric(.)) %>%
  .[vapply(., length, integer(1)) > 1] %>%
  lapply(names)


db %>%
  .[vapply(., function(x) x@name %in% cl[[138]], logical(1))] %>%
  merge_motifs(normalise.scores = TRUE) %>%
  view_motifs(names.pos = "right")

db %>%
  .[vapply(., function(x) !x@name %in% unlist(cl), logical(1))] %>%
  length()
unlist(cl) %>% length()
