#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})
.ub_norm <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("\\s*\\([^)]*\\)\\s*$", "", x)
  x <- gsub("[\u2010-\u2015\u2212]", "-", x)
  x <- gsub("\\s+", " ", x)
  x <- gsub(" ?- ?", "-", x)
  trimws(x)
}
uberon_term_blocks <- function(lines) {
  starts <- which(lines == "[Term]")
  if (!length(starts)) {
    return(list())
  }
  ends <- c(starts[-1] - 1L, length(lines))
  Map(function(first, last) lines[first:last], starts, ends)
}
parse_uberon_term <- function(block) {
  id_line <- block[startsWith(block, "id: UBERON:")]
  if (!length(id_line) || any(block == "is_obsolete: true")) {
    return(NULL)
  }
  ontology_id <- sub("^id: ", "", id_line[1])
  name_line <- block[startsWith(block, "name: ")]
  name <- if (length(name_line)) sub("^name: ", "", name_line[1]) else NA_character_
  edge_lines <- block[
    startsWith(block, "is_a:") |
      startsWith(block, "relationship: part_of")
  ]
  parents <- unique(unlist(regmatches(
    edge_lines, regexpr("UBERON:\\d+", edge_lines)
  )))
  synonym_lines <- block[startsWith(block, "synonym: ")]
  synonyms <- vapply(synonym_lines, function(line) {
    match <- regexec('^synonym: "([^"]*)" *(EXACT).*$', line)
    groups <- regmatches(line, match)[[1]]
    if (length(groups) >= 3L) groups[2] else NA_character_
  }, "")
  list(
    id = ontology_id,
    name = name,
    parents = parents[nzchar(parents)],
    synonyms = unique(synonyms[!is.na(synonyms) & nzchar(synonyms)])
  )
}
unambiguous_name_map <- function(terms) {
  primary <- list()
  synonyms <- list()
  for (term in terms) {
    if (!is.na(term$name)) {
      key <- .ub_norm(term$name)
      primary[[key]] <- c(primary[[key]], term$id)
    }
    for (synonym in term$synonyms) {
      key <- .ub_norm(synonym)
      synonyms[[key]] <- c(synonyms[[key]], term$id)
    }
  }
  unique_id <- function(values) {
    values <- unique(values)
    if (length(values) == 1L) values[1] else NA_character_
  }
  primary_values <- vapply(primary, unique_id, "")
  primary_values <- primary_values[!is.na(primary_values)]
  synonym_values <- vapply(synonyms, unique_id, "")
  synonym_values <- synonym_values[!is.na(synonym_values)]
  synonym_values <- synonym_values[
    !(names(synonym_values) %in% names(primary_values))
  ]
  c(primary_values, synonym_values)
}
invert_uberon_parents <- function(parents) {
  children <- list()
  for (child in names(parents)) {
    for (parent in unique(parents[[child]])) {
      children[[parent]] <- c(children[[parent]], child)
    }
  }
  lapply(children, unique)
}
uberon_load <- function(obo) {
  terms <- Filter(Negate(is.null), lapply(
    uberon_term_blocks(readLines(obo, warn = FALSE)),
    parse_uberon_term
  ))
  ids <- vapply(terms, function(term) term$id, "")
  names <- vapply(terms, function(term) term$name, "")
  parents <- setNames(lapply(terms, function(term) term$parents), ids)
  list(
    id2name = setNames(names, ids),
    parents = parents,
    children = invert_uberon_parents(parents),
    namemap = unambiguous_name_map(terms)
  )
}
uberon_name <- function(ub, id) {
  v <- if (!is.na(id)) ub$id2name[[id]] else NA
  if (is.null(v)) v <- NA
  if (is.na(v)) id else v
}
uberon_resolve <- function(ub, name) {
  if (is.null(name) || is.na(name) || !nzchar(name)) {
    return(NA_character_)
  }
  s <- trimws(as.character(name))
  if (grepl("^UBERON:", toupper(s))) {
    return(toupper(s))
  }
  hit <- ub$namemap[[.ub_norm(s)]]
  if (is.null(hit)) NA_character_ else hit
}
uberon_descendants <- function(ub, id) {
  if (is.null(id) || is.na(id)) {
    return(character(0))
  }
  seen <- character(0)
  frontier <- ub$children[[id]]
  if (is.null(frontier)) frontier <- character(0)
  while (length(frontier)) {
    nxt <- character(0)
    for (f in frontier) {
      if (!(f %in% seen)) {
        seen <- c(seen, f)
        children <- ub$children[[f]]
        if (!is.null(children)) nxt <- c(nxt, children)
      }
    }
    frontier <- unique(nxt)
  }
  unique(seen)
}
uberon_closure_ids <- function(ub, root) {
  rid <- uberon_resolve(ub, root)
  if (is.na(rid)) {
    return(character(0))
  }
  unique(c(rid, uberon_descendants(ub, rid)))
}
uberon_closure_names <- function(ub, root) {
  ids <- uberon_closure_ids(ub, root)
  if (!length(ids)) {
    return(character(0))
  }
  unique(.ub_norm(vapply(ids, function(i) uberon_name(ub, i), "")))
}
.scma_one_tissue_closure <- function(ub, tissue, tissue_root) {
  tr <- tolower(trimws(if (!nzchar(tissue_root)) "self" else tissue_root))
  ids <- character(0)
  nms <- .ub_norm(tissue)
  rid <- uberon_resolve(ub, tissue)
  if (!is.na(rid)) ids <- rid
  if (tr %in% c("exact", "none", "off")) {
    return(list(ids = ids, names = nms))
  }
  root <- if (tr %in% c("self", "")) tissue else tissue_root
  cids <- uberon_closure_ids(ub, root)
  if (length(cids)) {
    ids <- unique(c(ids, cids))
    nms <- unique(c(nms, uberon_closure_names(ub, root)))
  }
  list(ids = ids, names = nms)
}
scma_tissue_list <- function(tissue) {
  values <- trimws(as.character(tissue))
  values <- values[!is.na(values) & nzchar(values)]
  unique(values)
}
scma_parse_tissue_arg <- function(value) {
  text <- if (is.null(value) || is.na(value)) "" else as.character(value)[1]
  if (!grepl("|", text, fixed = TRUE)) {
    return(text)
  }
  names <- scma_tissue_list(strsplit(text, "|", fixed = TRUE)[[1]])
  if (!length(names)) stop("--tissue lists no tissue name", call. = FALSE)
  names
}
scma_tissue_closures <- function(ub, tissue, tissue_root) {
  names <- scma_tissue_list(tissue)
  out <- lapply(names, function(name) .scma_one_tissue_closure(ub, name, tissue_root))
  names(out) <- names
  out
}
scma_tissue_closure <- function(ub, tissue, tissue_root) {
  ids <- character(0)
  nms <- character(0)
  for (closure in scma_tissue_closures(ub, tissue, tissue_root)) {
    ids <- unique(c(ids, closure$ids))
    nms <- unique(c(nms, closure$names))
  }
  list(ids = ids, names = nms)
}
