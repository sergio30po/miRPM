counts <- matrix(
  c(
    120, 135, 110, 125, 850, 900, 780, 920,
    300, 320, 310, 305, 295, 315, 300, 325,
    700, 750, 680, 720, 130, 120, 150, 140,
    0,   2,   0,   1,   3,   0,   2,   1,
    50,  60,  55,  45, 200, 210, 190, 220,
    20, 200,  15, 180,  25, 210,  18, 190
  ),
  nrow = 6,
  byrow = TRUE,
  dimnames = list(
    c("miR-1", "miR-2", "miR-3", "miR-4", "miR-5", "miR-6"),
    paste0("S", 1:8)
  )
)

storage.mode(counts) <- "integer"

metadata <- data.frame(
  Sample = paste0("S", 1:8),
  Condition = c(
    rep("Control", 4),
    rep("Disease", 4)
  ),
  stringsAsFactors = FALSE
)

stopifnot(
  identical(colnames(counts), metadata$Sample),
  all(counts >= 0),
  all(is.finite(counts))
)

miRPM_example <- list(
  counts = counts,
  metadata = metadata
)

dir.create(
  "data",
  recursive = TRUE,
  showWarnings = FALSE
)

save(
  miRPM_example,
  file = file.path("data", "miRPM_example.rda"),
  compress = "xz",
  version = 2
)

rm(counts, metadata, miRPM_example)
