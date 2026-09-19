# Driven by test_borrowed_sources.py with the same synthetic sources and frozen results.
args <- commandArgs(TRUE)
root <- args[[1]]
fixture <- args[[2]]
e <- new.env(parent = globalenv())
for (name in c('config.R', 'uberon_ontology.R', 'evidence_gate.R', 'marker_sources.R',
               'ortho_map.R', 'annotator_pool.R', 'borrowed_context.R', 'reporting.R')) {
  source(file.path(root, 'rflow', name), local = e)
}
e$DB_SOURCES <- file.path(fixture, 'sources.csv')
e$ORTHO_DIR <- fixture
e$CROSS_SPECIES <- FALSE
e$uberon_load <- function(...) NULL
e$scma_tissue_closure <- function(ub, tissue, tissue_root) list(names = tolower(tissue))
fr <- jsonlite::fromJSON(file.path(fixture, 'fr.json'), simplifyVector = FALSE)
fr$dm$de <- data.table::rbindlist(fr$dm$de)
fr$scoring$marker_specificity <- unlist(fr$scoring$marker_specificity)
after <- e$scma_build_marker_evidence('fixture', fr)
native <- e$scma_source_context(NULL, e$DB_SOURCES, 'Human', 'liver', 'Normal')
before <- e$scma_build_marker_evidence('fixture', fr, sources = native)
cols <- setdiff(names(after), c('pmid', 'pmcid', 'source_sentence'))
stopifnot(identical(before[, ..cols], after[, ..cols]))
stopifnot(after[candidate_annotation == 'borrowed', pmcid] == 'PMC2')
stopifnot(after[candidate_annotation == 'donated', pmcid] == 'PMC3')
stopifnot(after[candidate_annotation == 'unborrowed', pmcid] == 'N/A')
extra <- e$scma_source_types_across_tissues(e$DB_SOURCES, 'Human', 'donated', 'Normal',
  donor_species_by_name = list(donated = c('Mouse', 'Rat')), ortho_dir = fixture)
recs <- e$scma_source_records(extra, 'donated', 'TARGET')
stopifnot(setequal(vapply(recs, function(r) r$pmcid, ''), c('PMC3', 'PMC4')))
for (gene in c('A', 'B', 'Unmapped')) stopifnot(!length(e$scma_source_records(extra, 'donated', gene)))
events <- list(a = list(borrowed = list(list(cell_type = 'donated', donor_species = 'Mouse'))),
               b = list(borrowed = list(list(cell_type = 'donated', donor_species = 'Rat'))))
stopifnot(setequal(e$scma_borrowed_donor_species(events)$donated, c('Mouse', 'Rat')))
data.table::fwrite(after, file.path(fixture, 'r_table.csv'), quote = TRUE)
cat('R borrowed-source regressions passed\n')
