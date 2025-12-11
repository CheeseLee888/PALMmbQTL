read_table_with_id <- function(path, id_col = NULL) {
  # path: file path
  # id_col:
  #   - NULL  -> treat the *first column* as ID (A-type table),
  #              rename it to IID.
  #   - name  -> treat the column with this name as ID (B-type table),
  #              it must exist; then rename it to IID.
  
  # fread returns a data.frame if data.table = FALSE
  df <- data.table::fread(
    path,
    data.table   = FALSE,
    check.names  = FALSE
  )
  
  ################## Set id_col=IID for now if user does not specify (keep case B only) ##################
  if (is.null(id_col)) {
    ## Case A: user does NOT specify an ID column
    ## -> use the first column as ID and rename it to IID
    if (ncol(df) < 1) {
      stop(sprintf("File '%s' has no columns.", path))
    }
    
    # Force first column to be character to avoid 001 -> 1 problems
    df[[1]] <- as.character(df[[1]])
    colnames(df)[1] <- "IID"
    
  } else {
    ## Case B: user specifies the ID column name
    ## -> find that column; if missing, throw an error
    if (!id_col %in% colnames(df)) {
      stop(sprintf(
        "ID column '%s' not found in file '%s'.",
        id_col, path
      ))
    }
    
    # Ensure ID column is character
    df[[id_col]] <- as.character(df[[id_col]])
    
    # Rename specified ID column to IID
    colnames(df)[colnames(df) == id_col] <- "IID"
  }
  
  # Move IID to the first column (for readability)
  df <- df[, c("IID", setdiff(colnames(df), "IID"))]
  
  return(df)
}

merge_abd_cov <- function(abd, cov) {
  # abd and cov must both contain an IID column
  if (!("IID" %in% names(abd))) {
    stop("abd does not contain column 'IID'.")
  }
  if (!("IID" %in% names(cov))) {
    stop("cov does not contain column 'IID'.")
  }
  
  merged <- merge(
    cov,
    abd,
    by   = "IID",
    sort = FALSE
  )
  
  merged <- merged[, c("IID", setdiff(names(merged), "IID"))]
  return(merged)
}
