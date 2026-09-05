# =============================================================================
# viewer.R -- R-native builder of the portable, reference-free dataset viewer.
#
# Every arm reads the SAME CSS / JS / HTML skeleton from resources/static/viewer/ and
# splices the same payload into it, so `viewer/index.html` is identical down to the byte
# no matter which arm wrote it -- which is exactly what the parity comparator asserts
# for .html files.
#
# Reads only delivered files (tables, JSONL, figures); never re-decides anything.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

.viewer_script_dir <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      value <- frame$ofile
      if (is.null(value)) "" else as.character(value)[1]
    },
    ""
  )
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files)) {
    return(dirname(normalizePath(frame_files[length(frame_files)])))
  }
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep("^--file=", args, value = TRUE)
  if (length(match)) {
    return(dirname(normalizePath(sub("^--file=", "", match[1]))))
  }
  normalizePath(getwd())
}
.VIEWER_SELF <- .viewer_script_dir()
.VIEWER_STATIC <- file.path(
  dirname(.VIEWER_SELF), "resources", "static", "viewer"
)

.viewer_read_static <- function(name) {
  path <- file.path(.VIEWER_STATIC, name)
  if (!file.exists(path)) {
    stop(sprintf("missing shared viewer asset: %s", path), call. = FALSE)
  }
  readChar(path, file.info(path)$size, useBytes = TRUE)
}

#' Replace every occurrence of `token` with `value`, both taken literally.
#'
#' gsub() interprets backslashes in its replacement even under fixed = TRUE, and the
#' payload is full of them; splitting sidesteps the whole escape question.
.viewer_splice <- function(text, token, value) {
  parts <- strsplit(text, token, fixed = TRUE)[[1]]
  paste(parts, collapse = value)
}

#' Escape the five HTML entities, ampersand first.
.viewer_html_escape <- function(text) {
  text <- gsub("&", "&amp;", text, fixed = TRUE)
  text <- gsub("<", "&lt;", text, fixed = TRUE)
  text <- gsub(">", "&gt;", text, fixed = TRUE)
  text <- gsub("\"", "&quot;", text, fixed = TRUE)
  gsub("'", "&#x27;", text, fixed = TRUE)
}

#' JSON string escaping: quote, backslash and control characters; everything else raw
#' UTF-8. Kept local so the viewer stays sourceable on its own.
.viewer_json_escape <- function(text) {
  text <- enc2utf8(as.character(text))
  text <- gsub("\\", "\\\\", text, fixed = TRUE)
  text <- gsub("\"", "\\\"", text, fixed = TRUE)
  text <- gsub("\n", "\\n", text, fixed = TRUE)
  text <- gsub("\r", "\\r", text, fixed = TRUE)
  text <- gsub("\t", "\\t", text, fixed = TRUE)
  if (grepl("[\x01-\x1f]", text, useBytes = TRUE)) {
    for (code in c(1:7, 11L, 14:31)) {
      text <- gsub(intToUtf8(code), sprintf("\\u%04x", code), text, fixed = TRUE)
    }
    text <- gsub("\b", "\\b", text, fixed = TRUE)
    text <- gsub("\f", "\\f", text, fixed = TRUE)
  }
  paste0("\"", text, "\"")
}

#' An integer written with thousands separators.
.viewer_thousands <- function(n) {
  format(as.integer(n), big.mark = ",", scientific = FALSE, trim = TRUE)
}

#' A delivered CSV as an array of JSON objects, one per row, keys in column order,
#' every value the exact string the file carries.
.viewer_table_json <- function(path) {
  table <- fread(
    path,
    colClasses = "character", na.strings = NULL, encoding = "UTF-8",
    keepLeadingZeros = TRUE
  )
  if (!nrow(table)) {
    return("[]")
  }
  columns <- names(table)
  keys <- vapply(columns, .viewer_json_escape, "")
  rows <- vapply(
    seq_len(nrow(table)),
    function(i) {
      values <- vapply(
        columns,
        function(column) .viewer_json_escape(table[[column]][i]),
        ""
      )
      paste0("{", paste(paste0(keys, ":", values), collapse = ","), "}")
    },
    ""
  )
  paste0("[", paste(rows, collapse = ","), "]")
}

# The five figure panels the viewer offers, in display order; absent files are skipped.
.VIEWER_FIGURES <- list(
  list(id = "identity-markers", label = "Identity markers",
       file = "dotplot_celltype_markers.png"),
  list(id = "identity-markers-celltype", label = "Identity markers by cell type",
       file = "dotplot_celltype_markers_by_celltype.png"),
  list(id = "identity-markers-support", label = "Publication support",
       file = "dotplot_celltype_markers_publication_support.png"),
  list(id = "cell-types", label = "Cell types", file = "umap_celltype_labeled.png"),
  list(id = "clusters", label = "Clusters", file = "umap_cluster.png")
)

.VIEWER_DOWNLOAD_TABLES <- c(
  "cluster_summary.csv", "marker_evidence.csv", "cluster_evidence.jsonl",
  "markers_all_by_cluster.csv", "markers_significant_by_cluster.csv"
)

.viewer_copy_files <- function(result_dir, viewer_dir) {
  assets <- file.path(viewer_dir, "assets", "figures")
  downloads <- file.path(viewer_dir, "downloads")
  dir.create(assets, recursive = TRUE)
  dir.create(downloads, recursive = TRUE)
  copied <- character(0)
  for (figure in .VIEWER_FIGURES) {
    source_path <- file.path(result_dir, "figures", figure$file)
    if (file.exists(source_path)) {
      file.copy(source_path, file.path(assets, figure$file), copy.date = TRUE)
      copied <- c(copied, file.path("assets", "figures", figure$file))
    }
  }
  for (name in .VIEWER_DOWNLOAD_TABLES) {
    source_path <- file.path(result_dir, name)
    if (!file.exists(source_path)) next
    file.copy(source_path, file.path(downloads, name), copy.date = TRUE)
    copied <- c(copied, file.path("downloads", name))
  }
  pdfs <- sort(list.files(file.path(result_dir, "figures"), pattern = "\\.pdf$"))
  for (name in pdfs) {
    file.copy(
      file.path(result_dir, "figures", name),
      file.path(downloads, name),
      copy.date = TRUE
    )
    copied <- c(copied, file.path("downloads", name))
  }
  copied
}

.viewer_figure_html <- function(viewer_dir) {
  available <- Filter(
    function(figure) {
      file.exists(file.path(viewer_dir, "assets", "figures", figure$file))
    },
    .VIEWER_FIGURES
  )
  tabs <- vapply(
    seq_along(available),
    function(index) {
      figure <- available[[index]]
      sprintf(
        '<button class="tab figure-tab%s" data-figure="%s">%s</button>',
        if (index == 1L) " active" else "",
        figure$id,
        .viewer_html_escape(figure$label)
      )
    },
    ""
  )
  images <- vapply(
    seq_along(available),
    function(index) {
      figure <- available[[index]]
      sprintf(
        '<img id="%s" class="%s" src="assets/figures/%s" alt="%s">',
        figure$id,
        if (index == 1L) "active" else "",
        .viewer_html_escape(figure$file),
        .viewer_html_escape(figure$label)
      )
    },
    ""
  )
  list(tabs = paste(tabs, collapse = ""), images = paste(images, collapse = ""))
}

.viewer_build_html <- function(tag, result_dir, viewer_dir) {
  clusters <- fread(
    file.path(result_dir, "cluster_summary.csv"),
    colClasses = "character", na.strings = NULL, encoding = "UTF-8"
  )
  total_cells <- sum(as.integer(clusters$n_cells))
  resolutions <- sort(
    unique(clusters$resolution_status[nzchar(clusters$resolution_status)]),
    method = "radix"
  )
  options <- paste0(
    '<option value="">All resolution states</option>',
    paste(
      vapply(
        resolutions,
        function(value) {
          escaped <- .viewer_html_escape(value)
          sprintf('<option value="%s">%s</option>', escaped, escaped)
        },
        ""
      ),
      collapse = ""
    )
  )
  figures <- .viewer_figure_html(viewer_dir)

  # The payload: the two tables re-serialized cell for cell, and the audit sidecar's
  # own lines spliced in verbatim -- parsing and re-serializing a line this package
  # wrote returns that line, so embedding the bytes equals re-serializing them.
  evidence_lines <- readLines(
    file.path(result_dir, "cluster_evidence.jsonl"),
    encoding = "UTF-8", warn = FALSE
  )
  evidence_lines <- evidence_lines[nzchar(trimws(evidence_lines))]
  data <- paste0(
    '{"clusters":', .viewer_table_json(file.path(result_dir, "cluster_summary.csv")),
    ',"markers":', .viewer_table_json(file.path(result_dir, "marker_evidence.csv")),
    ',"evidence":[', paste(evidence_lines, collapse = ","), "]}"
  )
  data <- gsub("</", "<\\\\/", data, fixed = TRUE)
  data <- gsub("<!--", "<\\\\!--", data, fixed = TRUE)

  page <- .viewer_read_static("index_template.html")
  replacements <- list(
    "__SCMA_TAG__" = .viewer_html_escape(tag),
    "__SCMA_CSS__" = .viewer_read_static("viewer.css"),
    "__SCMA_JS__" = .viewer_read_static("viewer.js"),
    "__SCMA_TOTAL_CELLS__" = .viewer_thousands(total_cells),
    "__SCMA_N_CLUSTERS__" = as.character(nrow(clusters)),
    "__SCMA_N_MARKERS__" = .viewer_thousands(
      nrow(fread(
        file.path(result_dir, "marker_evidence.csv"),
        colClasses = "character", na.strings = NULL, encoding = "UTF-8"
      ))
    ),
    "__SCMA_RESOLUTION_OPTIONS__" = options,
    "__SCMA_FIGURE_TABS__" = figures$tabs,
    "__SCMA_FIGURE_IMAGES__" = figures$images,
    "__SCMA_DATA__" = data
  )
  for (token in names(replacements)) {
    page <- .viewer_splice(page, token, replacements[[token]])
  }
  page
}

.viewer_write_zip <- function(viewer_dir, output) {
  if (!requireNamespace("zip", quietly = TRUE)) {
    stop("the 'zip' package is required to write the results archive", call. = FALSE)
  }
  members <- sort(list.files(viewer_dir, recursive = TRUE, all.files = TRUE))
  if (file.exists(output)) unlink(output)
  # mirror mode stores each member under its path relative to root -- cherry-pick
  # would flatten them to basenames, which is a different archive, not a different
  # spelling of this one.
  zip::zip(
    zipfile = output, files = members, root = viewer_dir,
    mode = "mirror", include_directories = FALSE
  )
  invisible(output)
}

#' Build the viewer directory, its manifest and the delivered ZIP for one dataset --
#' the counterpart of `viewer.build_dataset_viewer`.
scma_build_dataset_viewer <- function(tag, outdir) {
  result_dir <- normalizePath(outdir)
  viewer_dir <- file.path(result_dir, "viewer")
  if (dir.exists(viewer_dir)) {
    unlink(viewer_dir, recursive = TRUE)
  }
  dir.create(viewer_dir, recursive = TRUE)
  copied <- .viewer_copy_files(result_dir, viewer_dir)

  html_path <- file.path(viewer_dir, "index.html")
  connection <- file(html_path, open = "wb")
  writeChar(
    enc2utf8(.viewer_build_html(tag, result_dir, viewer_dir)),
    connection,
    eos = NULL, useBytes = TRUE
  )
  close(connection)
  copied <- c(copied, "index.html")

  manifest <- list(
    schema_version = "scmarkeragent-dataset-viewer-v1",
    dataset = tag,
    reference_free = TRUE,
    generated_at_utc = strftime(
      as.POSIXlt(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%OS6+00:00"
    ),
    entrypoint = "index.html",
    files = as.list(sort(copied, method = "radix"))
  )
  manifest_path <- file.path(viewer_dir, "viewer_manifest.json")
  manifest_json <- jsonlite::toJSON(
    manifest,
    auto_unbox = TRUE, pretty = 2, null = "null"
  )
  manifest_connection <- file(manifest_path, open = "wb")
  writeChar(
    enc2utf8(paste0(as.character(manifest_json), "\n")),
    manifest_connection,
    eos = NULL, useBytes = TRUE
  )
  close(manifest_connection)

  zip_path <- file.path(result_dir, sprintf("%s_results.zip", tag))
  .viewer_write_zip(viewer_dir, zip_path)
  list(
    viewer = normalizePath(html_path),
    viewer_zip = normalizePath(zip_path)
  )
}
