# Run from fhs_application/ with the dbGaP files in data_raw/.
library(igraph)
library(ape)

# Read the requested columns from both consent groups.
paths <- list.files("data_raw", recursive = TRUE, full.names = TRUE)
read_pair <- function(prefix, fields) {
  parts <- lapply(c("c1", "c2"), function(consent) {
    path <- paths[startsWith(basename(paths), paste0(prefix, consent, ".")) &
                    grepl("[.]txt[.]gz$", basename(paths))]
    con <- gzfile(path, "rt")
    lines <- readLines(con, n = 100, warn = FALSE)
    close(con)
    lines <- gsub('"', '', sub("^#\\s*", "", lines, perl = TRUE), fixed = TRUE)
    header <- which(vapply(strsplit(lines, "\t", fixed = TRUE),
                           function(x) "shareid" %in% x, logical(1)))
    columns <- strsplit(lines[header], "\t", fixed = TRUE)[[1]]
    con <- gzfile(path, "rt")
    on.exit(close(con))
    data <- read.delim(con, skip = header, header = FALSE, col.names = columns,
      colClasses = ifelse(columns %in% fields, "character", "NULL"),
      quote = '"', comment.char = "", fill = FALSE, check.names = FALSE,
      na.strings = c("", ".", "NA"), stringsAsFactors = FALSE)
    data <- data[, fields, drop = FALSE]
    data$shareid <- trimws(data$shareid)
    data$consent <- consent
    attr(data, "source_file") <- basename(path)
    data
  })
  data <- do.call(rbind, parts)
  rownames(data) <- NULL
  attr(data, "source_files") <- vapply(parts, function(x) attr(x, "source_file"), character(1))
  data
}

positive <- function(x) !is.na(x) & is.finite(x) & x > 0

e4_clinic <- read_pair(
  "phs000007.v35.pht000033.v13.p16.",
  c("shareid", "idtype", "D401", "D402", "D093",
    "D082", "D085", "D088", "D192", "D288")
)
e4_echo <- read_pair(
  "phs000007.v35.pht000300.v12.p16.",
  c("shareid", "IVSD4", "LVPD4", "LVDD4")
)
e5_clinic <- read_pair(
  "phs000007.v35.pht000034.v12.p16.",
  c("shareid", "E025")
)
e5_echo <- read_pair(
  "phs000007.v35.pht000301.v12.p16.",
  c("shareid", "idtype", "IVSD5", "LVPD5", "LVDD5")
)
demographics <- read_pair(
  "phs000007.v35.pht003099.v10.p16.",
  c("shareid", "age4", "sex", "date4")
)
relationships <- read_pair(
  "phs000153.v9.pht000836.v8.p8.",
  c("shareid", "sharealterid", "ALTERTYPE", "SPELLBEGIN", "SPELLEND")
)

source_files <- data.frame(
  table = c("E4 clinic", "E4 echo", "E5 clinic", "E5 echo",
            "Dates/demographics", "Social network"),
  c1_file = c(
    attr(e4_clinic, "source_files")[[1]],
    attr(e4_echo, "source_files")[[1]],
    attr(e5_clinic, "source_files")[[1]],
    attr(e5_echo, "source_files")[[1]],
    attr(demographics, "source_files")[[1]],
    attr(relationships, "source_files")[[1]]
  ),
  c2_file = c(
    attr(e4_clinic, "source_files")[[2]],
    attr(e4_echo, "source_files")[[2]],
    attr(e5_clinic, "source_files")[[2]],
    attr(e5_echo, "source_files")[[2]],
    attr(demographics, "source_files")[[2]],
    attr(relationships, "source_files")[[2]]
  ),
  stringsAsFactors = FALSE
)

names(e4_clinic)[names(e4_clinic) == "consent"] <- "consent_e4_clinic"
names(e4_echo)[names(e4_echo) == "consent"] <- "consent_e4_echo"
names(e5_clinic)[names(e5_clinic) == "consent"] <- "consent_e5_clinic"
names(e5_echo)[names(e5_echo) == "consent"] <- "consent_e5_echo"
names(e5_echo)[names(e5_echo) == "idtype"] <- "idtype5"
names(demographics)[names(demographics) == "consent"] <- "consent_demo"

d <- e4_clinic
for (z in list(e4_echo, e5_clinic, e5_echo, demographics))
  d <- merge(d, z, by = "shareid", all.x = TRUE, sort = FALSE)

consent_fields <- c("consent_e4_echo", "consent_e5_clinic",
                    "consent_e5_echo", "consent_demo")
nonnumeric <- c("shareid", "consent_e4_clinic", consent_fields)
for (field in setdiff(names(d), nonnumeric))
  d[[field]] <- as.numeric(d[[field]])

# Baseline covariates and height-indexed LVM at E4 and E5.
d$height4_m <- d$D402 * 0.0254
d$weight4_kg <- d$D401 * 0.45359237
d$bmi4 <- d$weight4_kg / d$height4_m^2
d$obesity4 <- as.integer(d$bmi4 >= 30)
d$female <- ifelse(d$sex == 2, 1L, ifelse(d$sex == 1, 0L, NA_integer_))
d$smoking4 <- d$D093
alcohol_complete <- !is.na(d$D082) & !is.na(d$D085) & !is.na(d$D088)
d$alcohol_any4 <- NA_integer_
d$alcohol_any4[alcohol_complete] <- as.integer(
  d$D082[alcohol_complete] > 0 |
    d$D085[alcohol_complete] > 0 |
    d$D088[alcohol_complete] > 0
)
d$sbp4 <- (d$D192 + d$D288) / 2
d$sbp_difference4 <- d$D192 - d$D288

# E4/E5 manuals: ASE/Devereux form; all dimensions are already in cm.
d$lvm4_g <- 0.8 * 1.04 * (
  (d$LVDD4 + d$IVSD4 + d$LVPD4)^3 - d$LVDD4^3
) + 0.6
d$lvm_height4 <- d$lvm4_g / d$height4_m
d$height5_m <- d$E025 * 0.0254
d$lvm5_g <- 0.8 * 1.04 * (
  (d$LVDD5 + d$IVSD5 + d$LVPD5)^3 - d$LVDD5^3
) + 0.6
d$lvm_height5 <- d$lvm5_g / d$height5_m


# Apply the cohort restrictions.
offspring <- d$idtype == 1 & d$idtype5 == 1
offspring[is.na(offspring)] <- FALSE
paired_measurements <- offspring
for (field in c("D401", "D402", "IVSD4", "LVPD4", "LVDD4",
                "E025", "IVSD5", "LVPD5", "LVDD5"))
  paired_measurements <- paired_measurements & positive(d[[field]])
adult <- paired_measurements & !is.na(d$age4) & d$age4 >= 21
baseline_complete <- adult & d$sex %in% 1:2 &
  alcohol_complete & positive(d$D192) & positive(d$D288) &
  d$D093 %in% 0:1
lvm_valid <- is.finite(d$lvm4_g) & d$lvm4_g > 0 &
  is.finite(d$lvm_height4) & d$lvm_height4 > 0 &
  is.finite(d$lvm5_g) & d$lvm5_g > 0 &
  is.finite(d$lvm_height5) & d$lvm_height5 > 0
eligible <- baseline_complete & lvm_valid
eligible[is.na(eligible)] <- FALSE

count_stage <- function(stage, keep) {
  data.frame(
    stage = stage,
    n = sum(keep),
    obese = sum(keep & d$obesity4 == 1, na.rm = TRUE),
    nonobese = sum(keep & d$obesity4 == 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
flow <- do.call(rbind, list(
  count_stage("E4 clinic Offspring with E5 Offspring echo", offspring),
  count_stage("Positive paired echo components, both heights, E4 weight", paired_measurements),
  count_stage("Plus E4 age >=21", adult),
  count_stage("Plus E4 sex, alcohol, both physician SBP, D093", baseline_complete),
  count_stage("Plus positive finite E4/E5 LVM and height-indexed LVM", eligible)
))

cohort <- data.frame(
  id = as.character(d$shareid[eligible]),
  consent_group = d$consent_e4_clinic[eligible],
  obesity4 = d$obesity4[eligible],
  lvm_height5 = d$lvm_height5[eligible],
  age4 = d$age4[eligible],
  female = d$female[eligible],
  smoking4 = d$smoking4[eligible],
  alcohol_any4 = d$alcohol_any4[eligible],
  sbp4 = d$sbp4[eligible],
  lvm_height4 = d$lvm_height4[eligible],
  height4_m = d$height4_m[eligible],
  bmi4 = d$bmi4[eligible],
  weight4_kg = d$weight4_kg[eligible],
  D192 = d$D192[eligible],
  D288 = d$D288[eligible],
  sbp_difference4 = d$sbp_difference4[eligible],
  lvm4_g = d$lvm4_g[eligible],
  height5_m = d$height5_m[eligible],
  lvm5_g = d$lvm5_g[eligible],
  date4 = d$date4[eligible],
  stringsAsFactors = FALSE
)
cohort <- cohort[order(cohort$id), , drop = FALSE]
rownames(cohort) <- NULL


# E4-active ties within the analysis cohort.
relationships$shareid <- trimws(as.character(relationships$shareid))
relationships$sharealterid <- trimws(as.character(relationships$sharealterid))
relationships$ALTERTYPE <- trimws(as.character(relationships$ALTERTYPE))
relationships$SPELLBEGIN <- as.numeric(relationships$SPELLBEGIN)
relationships$SPELLEND <- as.numeric(relationships$SPELLEND)

active <- !is.na(relationships$shareid) & nzchar(relationships$shareid) &
  !is.na(relationships$sharealterid) & nzchar(relationships$sharealterid) &
  relationships$shareid != relationships$sharealterid &
  !is.na(relationships$ALTERTYPE) & nzchar(relationships$ALTERTYPE) &
  is.finite(relationships$SPELLBEGIN) & relationships$SPELLBEGIN <= 252 &
  is.finite(relationships$SPELLEND) & relationships$SPELLEND > 192
active_edges <- relationships[active, c("shareid", "sharealterid"), drop = FALSE]
active_edges <- active_edges[
  active_edges$shareid %in% cohort$id &
    active_edges$sharealterid %in% cohort$id,
  ,
  drop = FALSE
]
edges <- unique(data.frame(
  from = pmin(active_edges$shareid, active_edges$sharealterid),
  to = pmax(active_edges$shareid, active_edges$sharealterid),
  stringsAsFactors = FALSE
))

graph <- igraph::graph_from_data_frame(
  edges,
  directed = FALSE,
  vertices = data.frame(name = cohort$id, stringsAsFactors = FALSE)
)
graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)

degree <- igraph::degree(graph)
components <- igraph::components(graph)
graph_summary <- data.frame(
  n_vertices = igraph::vcount(graph),
  n_edges = igraph::ecount(graph),
  n_isolates = sum(degree == 0),
  n_components = components$no,
  largest_component = max(components$csize),
  mean_degree = mean(degree),
  median_degree = stats::median(degree),
  max_degree = max(degree),
  spellbegin_max = 252,
  spellend_min_exclusive = 192,
  stringsAsFactors = FALSE
)

adjacency <- as.matrix(igraph::as_adjacency_matrix(graph, sparse = FALSE))
adjacency <- adjacency[cohort$id, cohort$id]
storage.mode(adjacency) <- "double"
diag(adjacency) <- 0

moran_variables <- c(
  "obesity4", "age4", "female", "smoking4",
  "alcohol_any4", "sbp4", "lvm_height4", "lvm_height5"
)
moran_report <- do.call(rbind, lapply(moran_variables, function(variable) {
  result <- ape::Moran.I(
    cohort[[variable]],
    adjacency,
    alternative = "two.sided"
  )
  data.frame(
    variable = variable,
    observed_I = unname(result$observed),
    expected_I = unname(result$expected),
    sd = unname(result$sd),
    p_value = unname(result$p.value),
    n_vertices = igraph::vcount(graph),
    n_edges = igraph::ecount(graph),
    n_isolates = sum(degree == 0),
    weight_input = "binary symmetric cohort-induced adjacency",
    test = "ape::Moran.I two-sided analytic p-value",
    stringsAsFactors = FALSE
  )
}))

make_table1 <- function(data) {
  groups <- list(
    Overall = seq_len(nrow(data)),
    `No E4 obesity` = which(data$obesity4 == 0),
    `E4 obesity` = which(data$obesity4 == 1)
  )
  continuous <- c(
    age4 = "Age, years",
    height4_m = "E4 height, m",
    bmi4 = "E4 BMI, kg/m^2",
    sbp4 = "E4 mean physician SBP, mm Hg",
    lvm_height4 = "E4 height-indexed LVM, g/m",
    lvm_height5 = "E5 height-indexed LVM, g/m"
  )
  binary <- c(
    female = "Female sex",
    smoking4 = "Past-year regular smoking",
    alcohol_any4 = "Any current alcohol use"
  )
  labels <- c(continuous, binary)
  long <- do.call(rbind, lapply(names(groups), function(group_name) {
    z <- data[groups[[group_name]], , drop = FALSE]
    do.call(rbind, lapply(names(labels), function(variable) {
      x <- z[[variable]]
      is_continuous <- variable %in% names(continuous)
      data.frame(
        group = group_name,
        variable = variable,
        characteristic = unname(labels[[variable]]),
        n = length(x),
        missing = sum(is.na(x)),
        mean = if (is_continuous) mean(x) else NA_real_,
        sd = if (is_continuous) stats::sd(x) else NA_real_,
        median = if (is_continuous) stats::median(x) else NA_real_,
        q25 = if (is_continuous) unname(stats::quantile(x, 0.25)) else NA_real_,
        q75 = if (is_continuous) unname(stats::quantile(x, 0.75)) else NA_real_,
        count = if (!is_continuous) sum(x == 1) else NA_integer_,
        percent = if (!is_continuous) 100 * mean(x == 1) else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  }))
  wide <- data.frame(
    Characteristic = c("N", unname(continuous), paste0(unname(binary), ", n (%)")),
    stringsAsFactors = FALSE
  )
  for (group_name in names(groups)) {
    subset <- long[long$group == group_name, , drop = FALSE]
    values <- ifelse(
      subset$variable %in% names(continuous),
      sprintf("%.2f (%.2f)", subset$mean, subset$sd),
      sprintf("%d (%.1f%%)", subset$count, subset$percent)
    )
    wide[[group_name]] <- c(as.character(length(groups[[group_name]])), values)
  }
  list(wide = wide, long = long)
}
table1 <- make_table1(cohort)

distance <- igraph::distances(graph, mode = "all")
distance <- distance[cohort$id, cohort$id, drop = FALSE]
diag(distance) <- 0


data_dir <- "data"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
analysis_data <- list(
  cohort = cohort,
  graph = graph,
  distance = distance,
  cohort_flow = flow,
  graph_summary = graph_summary,
  table1 = table1$wide,
  table1_long = table1$long,
  moran_report = moran_report,
  source_files = source_files,
  specification = list(
    exposure = "E4 BMI >= 30",
    outcome = "E5 height-indexed LVM",
    covariates = c("age4", "female", "smoking4", "alcohol_any4", "sbp4", "lvm_height4"),
    network_window = "SPELLBEGIN <= 252 and SPELLEND > 192"
  )
)
out <- file.path(data_dir, "fhs_analysis_data.rds")
saveRDS(analysis_data, out, compress = TRUE)

cat("\n--- E4->E5 cohort flow ---\n")
print(flow, row.names = FALSE)
cat("\n--- Table 1 ---\n")
print(table1$wide, row.names = FALSE)
cat("\n--- E4 graph summary ---\n")
print(graph_summary, row.names = FALSE)
cat("\n--- Moran's I (ape analytic p-values) ---\n")
print(moran_report[, c("variable", "observed_I", "expected_I", "sd", "p_value")],
      row.names = FALSE)

cat("\nSaved local analysis object: ", out, "\n", sep = "")
