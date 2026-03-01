# ============================================================================
# LLM API Functions
# ============================================================================
# GPT API calls, prompt construction, and response parsing

suppressPackageStartupMessages({
  if (!requireNamespace("writexl", quietly = TRUE)) {
    install.packages("writexl", repos = "https://cloud.r-project.org")
  }
  library(writexl)
})

#' Call GPT API with dual-channel fault tolerance (up to 10 retries each)
#' 
#' @param prompt Input prompt string
#' @param max_retries Maximum number of retries per channel
#' @return Text content returned by the API
call_gpt_api <- function(prompt, max_retries = 10) {
  
  get_provider <- function(kind = c("primary", "secondary")) {
    kind <- match.arg(kind)
    
    if (kind == "primary") {
      list(
        name = "primary",
        url = Sys.getenv("PRIMARY_API_URL", unset = "https://api.openai.com/v1/chat/completions"),
        model = Sys.getenv("PRIMARY_MODEL", unset = "gpt-4o"),
        key = Sys.getenv("PRIMARY_API_KEY", unset = "")
      )
    } else {
      list(
        name = "secondary",
        url = Sys.getenv("SECONDARY_API_URL", unset = ""),
        model = Sys.getenv("SECONDARY_MODEL", unset = ""),
        key = Sys.getenv("SECONDARY_API_KEY", unset = "")
      )
    }
  }
  
  provider_primary <- get_provider("primary")
  provider_secondary <- get_provider("secondary")
  
  if (provider_primary$key == "" && provider_secondary$key == "") {
    stop("No API key configured. Please set PRIMARY_API_KEY or SECONDARY_API_KEY.")
  }
  
  request_timeout_secs <- suppressWarnings(
    as.integer(Sys.getenv("API_REQUEST_TIMEOUT", unset = "3600"))
  )
  if (is.na(request_timeout_secs) || request_timeout_secs < 1) {
    request_timeout_secs <- 3600
  }
  
  # Generic API caller with retry logic
  call_api_with_retry <- function(provider, max_retries) {
    request_body <- list(
      model = provider$model,
      messages = list(list(role = "user", content = prompt)),
      temperature = 0
    )
    
    headers <- c(
      "Authorization" = paste("Bearer", provider$key),
      "Content-Type" = "application/json"
    )
    
    for (attempt in 1:max_retries) {
      response <- tryCatch({
        httr::POST(
          url = provider$url,
          httr::add_headers(.headers = headers),
          body = request_body,
          encode = "json",
          httr::timeout(request_timeout_secs)
        )
      }, error = function(e) {
        if (attempt < max_retries) {
          message(sprintf("  [%s] Request failed (attempt %d/%d): %s",
                         provider$name, attempt, max_retries, e$message))
        }
        NULL
      })
      
      if (!is.null(response) && httr::status_code(response) == 200) {
        response_content <- httr::content(response, "parsed")
        result <- response_content$choices[[1]]$message$content
        
        if (is.null(result) || nchar(trimws(result)) == 0) {
          if (attempt < max_retries) {
            message(sprintf("  [%s] Empty response received (attempt %d/%d)",
                           provider$name, attempt, max_retries))
          }
          next
        }
        
        return(result)
      } else if (!is.null(response)) {
        if (attempt < max_retries) {
          message(sprintf("  [%s] HTTP status %d (attempt %d/%d)",
                         provider$name, httr::status_code(response), attempt, max_retries))
        }
      }
      
      if (attempt < max_retries) {
        Sys.sleep(0.5)
      }
    }
    
    return(NULL)
  }
  
  # Try primary channel (10 retries)
  result <- call_api_with_retry(provider_primary, max_retries)
  if (!is.null(result)) {
    return(result)
  }
  
  # Fall back to secondary channel (10 retries)
  if (provider_secondary$key != "" && 
      provider_secondary$url != "" && 
      provider_secondary$model != "") {
    
    message("  Switching to secondary channel...")
    
    result <- call_api_with_retry(provider_secondary, max_retries)
    if (!is.null(result)) {
      return(result)
    }
  }
  
  stop(sprintf("API call failed: both primary and secondary channels exhausted (%d retries each)", max_retries))
}


#' Build a structured JSON annotation prompt for a single cluster
#' 
#' @param cluster_id Cluster identifier
#' @param n_cluster_markers Number of significant cluster markers
#' @param candidates_data Candidate data.frame with marker overlap information
#' @return Prompt string in JSON format
build_annotation_prompt <- function(cluster_id, n_cluster_markers, candidates_data) {
  
  # Use all candidates without hypergeometric pre-filtering
  sig_candidates <- candidates_data
  
  # Build candidate list
  candidates_list <- lapply(seq_len(nrow(sig_candidates)), function(i) {
    cand <- sig_candidates[i, ]
    
    # Parse overlap genes from comma-separated strings
    overlap_pos <- if (nchar(cand$overlap_genes_pos) > 0) {
      strsplit(cand$overlap_genes_pos, ", ")[[1]]
    } else {
      character(0)
    }
    
    overlap_neg <- if (nchar(cand$overlap_genes_neg) > 0) {
      strsplit(cand$overlap_genes_neg, ", ")[[1]]
    } else {
      character(0)
    }
    
    list(
      rank = cand$rank,
      celltype = cand$celltype,
      overlap_positive_markers = overlap_pos,
      overlap_negative_markers = overlap_neg
    )
  })
  
  # Build structured payload
  payload <- list(
    instruction = "Select the most biologically appropriate cell type annotation for this cluster from the ranked candidates.",
    
    rules = c(
      "If none of the candidates are biologically reasonable, return 'Unknown' as selected_celltype",
      "Candidates with negative markers overlapping cluster markers must be rejected due to biological conflict",
      "Select the candidate whose overlapping positive markers most coherently represent that cell type's lineage-defining features and functional identity in vivo",
      "Prefer candidates where overlapping markers match known cell-surface phenotypes used in flow cytometry and transcriptional signatures validated in single-cell studies",
      "When multiple candidates are indistinguishable based on the same marker evidence, prefer the most specific valid annotation consistent with that evidence",
      "The selected_celltype must exactly match one candidate celltype name or be 'Unknown'",
      "key_markers_validated must contain only genes from the selected candidate's overlap_positive_markers",
      "In the reasoning field, explain both why the selected cell type is most appropriate AND why other highly-ranked candidates were rejected"
    ),
    
    inputs = list(
      cluster_id = as.character(cluster_id),
      n_significant_cluster_markers = n_cluster_markers,
      candidates = candidates_list
    ),
    
    output_format = list(
      format = "minified JSON",
      schema = list(
        selected_celltype = "string (must be from candidates or 'Unknown')",
        confidence = "string (must be 'high', 'medium', or 'low')",
        reasoning = "string (explain why selected cell type is chosen AND why others were rejected)",
        key_markers_validated = "array of strings (2-5 key gene symbols)"
      )
    )
  )
  
  # Generate JSON prompt
  prompt <- paste0(
    "You are a cell biology and single-cell transcriptomics expert.\n\n",
    "CRITICAL INSTRUCTIONS:\n",
    "1. Reply with ONLY a single line of valid JSON - NO other text before or after\n",
    "2. Do NOT use markdown code blocks (no ```)\n",
    "3. Do NOT add explanations or comments\n",
    "4. The JSON must be minified (no line breaks inside the JSON)\n",
    "5. Start your response directly with { and end with }\n\n",
    "Required JSON schema:\n",
    '{"selected_celltype": "string", "confidence": "high|medium|low", "reasoning": "string", "key_markers_validated": ["gene1", "gene2"]}\n\n',
    "Task input:\n",
    jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = FALSE),
    "\n\nYour response (pure JSON only):"
  )
  
  return(prompt)
}


#' Parse and validate the LLM JSON response
#' 
#' @param response_text Raw text returned by the API
#' @param candidate_celltypes Valid candidate cell type names (for validation)
#' @return Parsed list with validated fields
parse_llm_response <- function(response_text, candidate_celltypes) {
  
  response_text <- trimws(response_text)
  
  # Strip markdown code fences if present
  if (grepl("^```", response_text)) {
    code_block_pattern <- "^```(?:json)?\\s*\\n?(.+?)```\\s*$"
    if (grepl(code_block_pattern, response_text, perl = TRUE)) {
      response_text <- gsub(code_block_pattern, "\\1", response_text, perl = TRUE)
      response_text <- trimws(response_text)
    }
  }
  
  # Parse JSON
  response <- tryCatch({
    jsonlite::fromJSON(response_text, simplifyVector = FALSE)
  }, error = function(e) {
    stop(paste0("JSON parsing failed: ", e$message, 
                ". Response preview: ", substr(response_text, 1, 200)))
  })
  
  if (is.null(response) || !is.list(response)) {
    stop("Parsed response is not a valid list object")
  }
  
  # Validate required fields
  required_fields <- c("selected_celltype", "confidence", "reasoning", "key_markers_validated")
  missing_fields <- setdiff(required_fields, names(response))
  if (length(missing_fields) > 0) {
    available_fields <- paste(names(response), collapse = ", ")
    stop(paste0("Missing required fields: ", paste(missing_fields, collapse = ", "),
                ". Available fields: ", available_fields))
  }
  
  # Validate selected_celltype
  valid_types <- c(candidate_celltypes, "Unknown")
  if (!response$selected_celltype %in% valid_types) {
    warning(sprintf(
      "LLM returned invalid celltype '%s', not in candidates. Setting to 'Unknown'.",
      response$selected_celltype
    ))
    response$selected_celltype <- "Unknown"
  }
  
  # Validate confidence
  if (!response$confidence %in% c("high", "medium", "low")) {
    warning(sprintf(
      "Invalid confidence '%s', setting to 'low'",
      response$confidence
    ))
    response$confidence <- "low"
  }
  
  # Ensure key_markers is a character vector
  if (!is.null(response$key_markers_validated)) {
    response$key_markers_validated <- unlist(response$key_markers_validated)
  } else {
    response$key_markers_validated <- character(0)
  }
  
  return(response)
}

cat("✓ LLM API functions loaded\n")
