capt_history.t %<>% 
  group_by(Mark) %>%
  mutate(index = seq(n())) %>%
  mutate(
    first_len_w = which(!is.na(len)) %>% first()
  , last_len_w  = which(!is.na(len)) %>% last()
  , first_len   = len[!is.na(len)] %>% first()
  , last_len    = len[!is.na(len)] %>% last()
    ) %>%
  mutate(
    len = ifelse(is.na(len) & (index < first_len_w), first_len, len)
  , len = ifelse(is.na(len) & (index > last_len_w), last_len, len)
  ) %>% group_by(Mark, Year) %>%
  mutate(
    len = ifelse(is.na(len), mean(len, na.rm = T), len)
  ) %>% ungroup() %>%
  group_by(Mark) %>%
  mutate(
    len = ifelse(is.na(len), mean(len, na.rm = T), len)
  ) %>% ungroup() %>%
  dplyr::select(
    -c(first_len_w, first_len, last_len_w, last_len)
  )
