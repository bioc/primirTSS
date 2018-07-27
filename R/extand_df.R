#' @importFrom stringr str_detect
#' @importFrom stringr str_split
#' @importFrom dplyr bind_cols

extand_df <- function(df, colnm, sep = ":") {
  info <- df[[colnm]]
  index <- str_detect(info, sep)
  if (any(index)) {
    double_tss_index <- which(index)
    split_res <- str_split(info[double_tss_index], sep)
    rep_times <- lapply(split_res, length) %>% unlist()
    new_add_seq <- rep(double_tss_index, rep_times)
    new_add_candidate_tss <- unlist(split_res)

    new_seq <- c(seq(length(info))[-double_tss_index], new_add_seq)
    new_info <- c(info[-double_tss_index], new_add_candidate_tss)
  } else {
    new_seq <- seq(length(info))
    new_info <- info
  }
  new_df <- bind_cols(df[new_seq, ], new_info = new_info)
  cmd <- sprintf("new_df <- select(new_df, -%s)", colnm)
  eval(parse(text = cmd))
  new_df
}
