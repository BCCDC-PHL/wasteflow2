#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Usage: parse_pafs_influenza_A.R <flu_db> <paf> <sample> <out_dir> <score_thresh>")

Flu_A_DB <- args[1]
paf_file <- args[2]
samp <- args[3]
out_dir <- args[4]
score_thresh <- as.numeric(args[5])

#-----------------------------------------
# Parameters
#-----------------------------------------
CHUNK_SIZE <- 2e6  # number of lines per read, adjust to memory (2e6 = ~500 MB typical)
cat(sprintf("Reading PAF in chunks of %s lines...\n", format(CHUNK_SIZE, big.mark = ",")))

#-----------------------------------------
# Read Flu info database
#-----------------------------------------
flu_info_dt <- fread(Flu_A_DB, sep = "\t", header = TRUE)
flu_info_dt[, accession := sub("\\.[0-9]+$", "", accession)]

#-----------------------------------------
# Function to process one chunk
#-----------------------------------------
process_chunk <- function(chunk_dt, flu_info_dt) {
  merge_dt <- merge(chunk_dt, flu_info_dt, by.x = "tname", by.y = "accession", all.x = TRUE)
  if (nrow(merge_dt) == 0) return(NULL)

  assignment_dt <- merge_dt %>%
    group_by(qname, tname, serotype, segment, strand) %>%
    summarize(
      tot_read_length = sum(qlength, na.rm = TRUE),
      tot_align = sum(align_length, na.rm = TRUE),
      tot_match = sum(num_matches, na.rm = TRUE),
      ANI = tot_match / tot_align,
      AF = tot_align / tot_read_length,
      align_score = ANI * AF,
      .groups = "drop"
    )
  return(assignment_dt)
}

#-----------------------------------------
# Stream the PAF file
#-----------------------------------------
con <- file(paf_file, "r")
header_cols <- c(
  "qname","qlength","qstart","qend","strand",
  "tname","tlength","tstart","tend",
  "num_matches","align_length","mapq"
)

all_assignments <- list()
chunk_index <- 0
lines_buffer <- character(0)

repeat {
  # Read next chunk
  chunk_lines <- readLines(con, n = CHUNK_SIZE)
  if (length(chunk_lines) == 0) break

  chunk_index <- chunk_index + 1
  cat(sprintf("\n📦 Processing chunk %d (%s lines)...\n", chunk_index, format(length(chunk_lines), big.mark = ",")))

  # Combine with leftover lines from previous iteration (if same qname)
  chunk_lines <- c(lines_buffer, chunk_lines)

  # Extract last qname to avoid splitting group
  last_fields <- strsplit(chunk_lines[length(chunk_lines)], "\t")[[1]]
  last_qname <- last_fields[1]

  # Read some lookahead lines until qname changes
  repeat {
    next_line <- readLines(con, n = 1)
    if (length(next_line) == 0) break
    if (strsplit(next_line, "\t")[[1]][1] != last_qname) {
      # different read name → save this for next chunk
      lines_buffer <- next_line
      break
    } else {
      # same qname, keep appending
      chunk_lines <- c(chunk_lines, next_line)
    }
  }

  # Parse this chunk to data.table
  paf_dt <- fread(
    text = chunk_lines,
    sep = "\t",
    header = FALSE,
    fill = TRUE,
    select = 1:12,
    col.names = header_cols,
    showProgress = FALSE
  )
  paf_dt[, tname := sub("\\.[0-9]+$", "", tname)]

  # Process and collect result
  chunk_result <- process_chunk(paf_dt, flu_info_dt)
  if (!is.null(chunk_result)) all_assignments[[length(all_assignments) + 1]] <- chunk_result
  rm(paf_dt); gc()
}

close(con)

#-----------------------------------------
# Combine all chunk results
#-----------------------------------------
assignment_dt <- bind_rows(all_assignments)
rm(all_assignments); gc()

cat(sprintf("\n✅ Total combined alignments: %s\n", format(nrow(assignment_dt), big.mark = ",")))

#-----------------------------------------
# Summarize and assign
#-----------------------------------------
sum_dt <- assignment_dt %>%
  ungroup() %>%
  group_by(qname, serotype, segment) %>%
  summarize(
    n = n(),
    top_score = max(align_score),
    avg_score = mean(align_score),
    .groups = "drop"
  ) %>%
  arrange(qname, desc(top_score)) %>%
  group_by(qname) %>%
  slice_max(top_score, n = 2, with_ties = TRUE) %>%
  mutate(
    read_assignment = case_when(
      (max(top_score - 0.003)) >= min(top_score) |
        n_distinct(serotype) == 1 ~ first(serotype),
      TRUE ~ "ambiguous"
    )
  ) %>%
  filter(max(top_score) >= score_thresh)

#-----------------------------------------
# Write outputs
#-----------------------------------------
out_summary <- file.path(out_dir, sprintf("%s_read_summary.tsv", samp))
write.table(sum_dt, file = out_summary, quote = FALSE, row.names = FALSE, sep = "\t")
cat(sprintf("✅ Summary written to %s\n", out_summary))

# Plot assignment counts
assign_plot <- sum_dt %>%
  slice_max(top_score, with_ties = FALSE) %>%
  ggplot(aes(x = read_assignment)) +
  geom_bar(fill = "steelblue") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Read assignment", y = "Count", title = paste0("Serotype assignment for ", samp))

out_plot <- file.path(out_dir, sprintf("%s_read_serotype_assignment.pdf", samp))
ggsave(assign_plot, file = out_plot, width = 7, height = 5)
cat(sprintf("📊 Plot saved to %s\n", out_plot))

# Per-serotype read lists
sum_dt %>%
  ungroup() %>%
  group_by(qname) %>%
  slice_max(top_score, with_ties = FALSE) %>%
  group_by(read_assignment) %>%
  select(read_assignment, qname) %>%
  group_walk(~ write.table(
    .x,
    file = file.path(out_dir, sprintf("%s_%s.txt", samp, .y$read_assignment)),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  ))

cat("\n🎯 Processing complete.\n")