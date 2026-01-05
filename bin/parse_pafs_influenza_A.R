#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  ### >>> CHANGED: prefer data.table + ggplot2 for speed (no tidyverse)
  library(data.table)
  library(ggplot2)
  ### <<< CHANGED
})

# ------------------------
# Argument parsing (--threads N support)
# ------------------------
raw_args <- commandArgs(trailingOnly = TRUE)

# default threads
threads <- 1L

# parse --threads flag (Option B)
if (length(raw_args) >= 2 && raw_args[1] == "--threads") {
  if (length(raw_args) < 3) stop("Usage: --threads N then 5 positional args")
  threads_val <- as.integer(raw_args[2])
  if (is.na(threads_val) || threads_val <= 0) stop("Invalid threads value: ", raw_args[2])
  threads <- threads_val
  # remaining args
  args <- raw_args[-(1:2)]
} else {
  args <- raw_args
}

if (length(args) != 5) {
  stop("Usage: parse_pafs_influenza_A_fast_safe.R [--threads N] <flu_db> <paf> <sample> <out_dir> <score_thresh>")
}

Flu_A_DB <- args[1]
paf_file <- args[2]
samp <- args[3]
out_dir <- args[4]
score_thresh <- as.numeric(args[5])
if (is.na(score_thresh)) stop("score_thresh must be numeric")

cat("Using threads =", threads, "\n")

# -----------------------------------------
# Parameters
# -----------------------------------------
### >>> CHANGED: tuned chunk size for large PAFs; adjust if memory constrained
CHUNK_SIZE <- 2e6  # number of lines per read chunk (tweak if memory issues)
cat(sprintf("Reading PAF in chunks of %s lines (SAFE mode: parse all columns)...\n", format(CHUNK_SIZE, big.mark = ",")))
### <<< CHANGED

# -----------------------------------------
# Read Flu info DB (data.table)
# -----------------------------------------
flu_info_dt <- fread(Flu_A_DB, sep = "\t", header = TRUE, nThread = threads)
flu_info_dt[, accession := sub("\\.[0-9]+$", "", accession)]
setkey(flu_info_dt, accession)

# -----------------------------------------
# Helper: parse chunk lines into data.table (SAFE: parse all columns)
# We will then keep only the first 12 mandatory PAF fields by position.
# -----------------------------------------
header_cols12 <- c(
  "qname","qlength","qstart","qend","strand",
  "tname","tlength","tstart","tend",
  "num_matches","align_length","mapq"
)

parse_chunk_safe <- function(lines, threads) {
  # fread will parse variable number of columns because optional tags exist
  # Use fill=TRUE to accept ragged rows
  dt <- tryCatch(
    fread(text = lines, sep = "\t", header = FALSE, fill = TRUE, showProgress = FALSE, nThread = threads),
    error = function(e) {
      warning("fread() failed on chunk — retrying with fill=Inf")
      fread(text = lines, sep = "\t", header = FALSE, fill = Inf, showProgress = FALSE, nThread = threads)
    }
  )

  # Ensure we have at least 12 columns
  if (ncol(dt) < 12) {
    # pad with NA columns
    for (i in (ncol(dt)+1):12) dt[, (paste0("V", i)) := NA]
  }

  # rename first 12 columns to header_cols12
  setnames(dt, old = names(dt)[1:12], new = header_cols12, skip_absent = TRUE)

  # Coerce necessary columns to appropriate types (numeric where expected)
  # Some columns might be read as character because of optional tags presence; coerce safely
  numcols <- c("qlength","qstart","qend","tlength","tstart","tend","num_matches","align_length","mapq")
  for (c in numcols) {
    if (c %in% names(dt)) dt[, (c) := as.numeric(get(c))]
  }

  # Normalize tname by stripping accession version (like original)
  dt[, tname := sub("\\.[0-9]+$", "", tname)]

  return(dt)
}

# -----------------------------------------
# Function: process one chunk -> aggregated assignments (data.table)
# -----------------------------------------
process_chunk_dt <- function(paf_dt, flu_info_dt) {
  # Left join like original (all.x = TRUE behaviour)
  # paf_dt may have tname values not in flu_info_dt; we want them preserved with NA serotype
  # We'll perform a fast data.table merge (setkey required)
  # ensure both have keys
  setkey(paf_dt, tname)
  # Use merge with all.x = TRUE to replicate original merge(..., all.x = TRUE)
  merged <- merge(paf_dt, flu_info_dt, by.x = "tname", by.y = "accession", all.x = TRUE, sort = FALSE)

  if (nrow(merged) == 0L) return(NULL)

  # compute aggregated statistics matching original dplyr summarize
  # group by qname, tname, serotype, segment, strand
  # Note: merged may not have 'segment' column in flu_info_dt; that's okay (will be NA)
  agg <- merged[,
    .(
      tot_read_length = sum(as.numeric(qlength), na.rm = TRUE),
      tot_align = sum(as.numeric(align_length), na.rm = TRUE),
      tot_match = sum(as.numeric(num_matches), na.rm = TRUE)
    ),
    by = .(qname, tname, serotype, segment, strand)
  ]

  # derived columns: ANI, AF, align_score (handle division by zero)
  agg[, ANI := fifelse(tot_align > 0, tot_match / tot_align, 0)]
  agg[, AF := fifelse(tot_read_length > 0, tot_align / tot_read_length, 0)]
  agg[, align_score := ANI * AF]

  return(agg[])
}

# -----------------------------------------
# Stream the PAF file (chunked read) - SAFE parsing of all columns
# -----------------------------------------
con <- file(paf_file, "r")
all_assignments <- list()
chunk_index <- 0L
lines_buffer <- character(0)

repeat {
  chunk_lines <- readLines(con, n = CHUNK_SIZE)
  if (length(chunk_lines) == 0) break

  chunk_index <- chunk_index + 1L
  cat(sprintf("\n📦 Processing chunk %d (%s lines)...\n", chunk_index, format(length(chunk_lines), big.mark = ",")))

  # prepend leftover lines from prior iteration (to keep qname groups intact)
  chunk_lines <- c(lines_buffer, chunk_lines)

  # Determine last qname to avoid splitting its group across chunks
  last_fields <- strsplit(chunk_lines[length(chunk_lines)], "\t")[[1]]
  last_qname <- last_fields[1]

  # read lookahead lines until qname changes (or EOF)
  repeat {
    next_line <- readLines(con, n = 1)
    if (length(next_line) == 0) {
      lines_buffer <- character(0)
      break
    }
    if (strsplit(next_line, "\t")[[1]][1] != last_qname) {
      lines_buffer <- next_line
      break
    } else {
      chunk_lines <- c(chunk_lines, next_line)
    }
  }

  # parse chunk safely (all columns)
  paf_dt <- parse_chunk_safe(chunk_lines, threads = threads)

  # process chunk -> aggregated assignments
  chunk_result <- process_chunk_dt(paf_dt, flu_info_dt)
  if (!is.null(chunk_result)) {
    all_assignments[[length(all_assignments) + 1L]] <- chunk_result
  }

  rm(paf_dt); gc()
}

close(con)

# -----------------------------------------
# Combine chunk results
# -----------------------------------------
if (length(all_assignments) == 0L) {
  stop("No assignments produced from PAF (empty result).")
}

### >>> CHANGED: use rbindlist for speed
assignment_dt <- rbindlist(all_assignments, use.names = TRUE, fill = TRUE)
rm(all_assignments); gc()
### <<< CHANGED

cat(sprintf("\n✅ Total combined alignments: %s\n", format(nrow(assignment_dt), big.mark = ",")))

# -----------------------------------------
# Summarize + top-2 + assignment logic (preserve exactly original)
# Implemented with data.table for speed (dense ranks for ties)
# -----------------------------------------
### >>> CHANGED: data.table equivalent of the original dplyr pipeline
# Step 1: summarise per qname, serotype, segment
sum_dt_dt <- assignment_dt[,
  .(
    n = .N,
    top_score = max(align_score, na.rm = TRUE),
    avg_score = mean(align_score, na.rm = TRUE)
  ),
  by = .(qname, serotype, segment)
]

# Step 2: order and pick up to top 2 per qname (including ties)
setorder(sum_dt_dt, qname, -top_score)
# dense rank per qname ensures ties share rank
sum_dt_dt[, rnk := frank(-top_score, ties.method = "dense"), by = qname]
top2_dt <- sum_dt_dt[rnk <= 2]
top2_dt[, rnk := NULL]

# Step 3: compute read_assignment per qname (exact logic)
assign_dt <- top2_dt[, {
  mx <- max(top_score - 0.003, na.rm = TRUE)
  mn <- min(top_score, na.rm = TRUE)
  n_serotypes <- uniqueN(serotype)
  ra <- if ((mx >= mn) | (n_serotypes == 1)) serotype[1L] else "ambiguous"
  # replicate original row expansion: return rows for top-two entries but include read_assignment value
  .(serotype = serotype, segment = segment, n = n, top_score = top_score, avg_score = avg_score, read_assignment = ra)
}, by = qname]

# Step 4: filter by score threshold (qnames whose max(top_score) >= score_thresh)
qname_max <- assign_dt[, .(qname_max_top = max(top_score, na.rm = TRUE)), by = qname]
keep_qnames <- qname_max[qname_max_top >= score_thresh, qname]
sum_dt_final <- assign_dt[qname %in% keep_qnames]

# Order to match original output ordering
setorder(sum_dt_final, qname, -top_score)
### <<< CHANGED

# -----------------------------------------
# Write summary TSV (fast fwrite + nThread)
# -----------------------------------------
out_summary <- file.path(out_dir, sprintf("%s_read_summary.tsv", samp))
fwrite(sum_dt_final, file = out_summary, sep = "\t", quote = FALSE, nThread = threads)
cat(sprintf("✅ Summary written to %s\n", out_summary))

# -----------------------------------------
# Plot assignment counts (keep original ggplot layout)
# -----------------------------------------
# replicate: slice_max(top_score, with_ties = FALSE) -> one row per qname with highest top_score
plot_dt <- unique(sum_dt_final[, .(qname, read_assignment, top_score)], by = c("qname", "read_assignment", "top_score"))
plot_dt <- plot_dt[order(qname, -top_score)]
plot_dt <- plot_dt[, .SD[1L], by = qname]

assign_plot <- ggplot(plot_dt, aes(x = read_assignment)) +
  geom_bar() +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Read assignment", y = "Count", title = paste0("Serotype assignment for ", samp))

out_plot <- file.path(out_dir, sprintf("%s_read_serotype_assignment.pdf", samp))
ggsave(out_plot, plot = assign_plot, width = 7, height = 5)
cat(sprintf("📊 Plot saved to %s\n", out_plot))

# -----------------------------------------
# Per-serotype read lists (fast, preserves logic)
# Keep top-per-qname (slice_max top_score, with_ties = FALSE) then write per serotype.
# -----------------------------------------
### >>> CHANGED: use data.table fast selection + fwrite in parallel
# top-per-qname: pick single row per qname with highest top_score
top_per_qname <- sum_dt_final[order(qname, -top_score), .SD[1L], by = qname]

# drop NA/empty assignments (same safeguard)
top_per_qname <- top_per_qname[!is.na(read_assignment) & read_assignment != ""]

serotypes <- unique(top_per_qname$read_assignment)
for (s in serotypes) {
  if (is.na(s) || s == "") {
    message("⚠️ Skipping invalid serotype group: ", s)
    next
  }
  outfile <- file.path(out_dir, sprintf("%s_%s.txt", samp, s))
  fwrite(top_per_qname[read_assignment == s, .(qname)], file = outfile, sep = "\t", quote = FALSE, col.names = FALSE, nThread = threads)
  message("📝 Wrote read list for serotype: ", s)
}
### <<< CHANGED

cat("\n🎯 Processing complete.\n")