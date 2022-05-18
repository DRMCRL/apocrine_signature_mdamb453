## Let's use the dht_consensus peaks to check the distance between summits as a measure of direct overlap
## Summits are estimated in a treatment-specific manner, so finding a target-specific summit
## i.e. treatment agnostic may be the first step
##
summits <- here::here("data", "macs2", targets) %>%
  lapply(list.files, pattern = "summits", full.names = TRUE) %>%
  setNames(targets) %>%
  lapply(
    function(x) {
      lapply(x, import.bed) %>%
        setNames(
          str_remove_all(basename(x), "_merged_summits.bed")
        )
    }
  )

merged_summits <- summits %>%
  lapply(GRangesList) %>%
  lapply(unlist) %>%
  lapply(sort) %>%
  lapply(names_to_column, var = "treat") %>%
  lapply(plyranges::select, all_of(c("treat", "name"))) %>%
  lapply(reduce, min.gapwidth = 100) %>%
  lapply(resize, width = 1, fix = 'center')

summits_sd <- dht_consensus %>%
  mutate(n_targets = AR + FOXA1 + GATA3 + TFAP2B) %>%
  filter(n_targets == 4) %>%
  as_tibble() %>%
  dplyr::select(range, n_targets) %>%
  mutate(
    ar_summit = findOverlaps(GRanges(range), merged_summits$AR) %>%
        as_tibble() %>%
        chop(subjectHits) %>%
        pull("subjectHits")
  ) %>%
  unnest(ar_summit) %>%
  mutate(
    ar_summit = merged_summits$AR[ar_summit] %>%
      start()
  ) %>%
  chop(ar_summit) %>%
  mutate(
    foxa1_summit = findOverlaps(GRanges(range), merged_summits$FOXA1) %>%
      as_tibble() %>%
      chop(subjectHits) %>%
      pull("subjectHits")
  ) %>%
  unnest(foxa1_summit) %>%
  mutate(
    foxa1_summit = merged_summits$FOXA1[foxa1_summit] %>%
      start()
  ) %>%
  chop(foxa1_summit)  %>%
  mutate(
    gata3_summit = findOverlaps(GRanges(range), merged_summits$GATA3) %>%
      as_tibble() %>%
      chop(subjectHits) %>%
      pull("subjectHits")
  ) %>%
  unnest(gata3_summit) %>%
  mutate(
    gata3_summit = merged_summits$GATA3[gata3_summit] %>%
      start()
  ) %>%
  chop(gata3_summit) %>%
  mutate(
    tfap2b_summit = findOverlaps(GRanges(range), merged_summits$TFAP2B) %>%
      as_tibble() %>%
      chop(subjectHits) %>%
      pull("subjectHits")
  ) %>%
  unnest(tfap2b_summit) %>%
  mutate(
    tfap2b_summit = merged_summits$TFAP2B[tfap2b_summit] %>%
      start()
  ) %>%
  chop(tfap2b_summit) %>%
  mutate(
    all_single = (
      vapply(ar_summit, length, integer(1)) == 1
    ) & (
      vapply(foxa1_summit, length, integer(1)) == 1
    ) & (
      vapply(gata3_summit, length, integer(1)) == 1
    ) & (
      vapply(tfap2b_summit, length, integer(1)) == 1
    )
  ) %>%
  dplyr::filter(all_single) %>%
  unnest(everything()) %>%
  pivot_longer(cols = ends_with("summit"), names_to = "target", values_to = "summit") %>%
  group_by(range) %>%
  summarise(sd = sd(summit), .groups = "drop") %>%
  arrange(sd) %>%
  mutate(rank = seq_along(sd))

summits_sd %>%
  ggplot(aes(sd, rank)) +
  geom_line()

## See if the proximity to each other predicts H3K27ac
summits_sd %>%
  colToRanges("range") %>%
  join_overlap_inner(dht_consensus) %>%
  mutate(mn = cummean(H3K27ac)) %>%
  as_tibble() %>%
  ggplot(aes(rank, mn)) +
  geom_line()
## No

## See if it predicts being mapped to a DE gene or an Apocrine Gene

