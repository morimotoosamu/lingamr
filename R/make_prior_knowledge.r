#' Create a prior knowledge matrix
#'
#' @param n_variables Number of variables
#' @param exogenous_variables Exogenous variables (1-based index or variable name, NULL allowed)
#'   The specified variables are assumed not to be influenced by any other variable
#' @param sink_variables Sink variables (1-based index or variable name, NULL allowed)
#'   The specified variables are assumed not to influence any other variable
#' @param paths Variable pairs that have a directed path (NULL allowed)
#'   Of the form list(c(from, to), ...). Specified by index or variable name
#' @param no_paths Variable pairs that have no directed path (NULL allowed)
#'   Of the form list(c(from, to), ...). Specified by index or variable name
#' @param labels Vector of variable names (NULL allowed)
#'   Required when specifying by variable name. Pass e.g. colnames() of a data.frame
#' @return Prior knowledge matrix (n_variables x n_variables)
#'   -1: unknown, 0: no path, 1: path exists
#' @export
#' @examples
#' # Specify by index
#' pk <- make_prior_knowledge(6, exogenous_variables = c(4))
#'
#' # Specify by variable name
#' pk <- make_prior_knowledge(6,
#'   exogenous_variables = "x3",
#'   sink_variables = c("x1", "x4"),
#'   paths = list(c("x3", "x0"), c("x3", "x2")),
#'   no_paths = list(c("x5", "x2")),
#'   labels = c("x0", "x1", "x2", "x3", "x4", "x5")
#' )
make_prior_knowledge <- function(n_variables,
                                 exogenous_variables = NULL,
                                 sink_variables = NULL,
                                 paths = NULL,
                                 no_paths = NULL,
                                 labels = NULL) {
  # --- Basic validation ---
  if (!is.numeric(n_variables) || n_variables < 2) {
    stop("n_variables must be an integer >= 2.")
  }
  n_variables <- as.integer(n_variables)

  # --- labels validation ---
  if (!is.null(labels)) {
    if (length(labels) != n_variables) {
      stop(sprintf(
        "Length of 'labels' (%d) must equal n_variables (%d).",
        length(labels), n_variables
      ))
    }
    if (anyDuplicated(labels) > 0) {
      stop("'labels' must not contain duplicate names.")
    }
  }

  # --- Helper to convert variable names to indices ---
  resolve_variable <- function(x, arg_name) {
    if (is.character(x)) {
      if (is.null(labels)) {
        stop(sprintf(
          "'%s' contains variable name(s), but 'labels' is not specified.\n  Provide labels or use integer indices.",
          arg_name
        ))
      }
      pos <- match(x, labels)
      if (any(is.na(pos))) {
        bad <- x[is.na(pos)]
        stop(sprintf(
          "Variable name(s) %s in '%s' not found in labels.\n  Available: %s",
          paste0("'", bad, "'", collapse = ", "),
          arg_name,
          paste(labels, collapse = ", ")
        ))
      }
      return(pos)
    } else if (is.numeric(x)) {
      x <- as.integer(x)
      if (any(x < 1) || any(x > n_variables)) {
        stop(sprintf("'%s' indices must be between 1 and %d.", arg_name, n_variables))
      }
      return(x)
    } else {
      stop(sprintf("'%s' must be integer indices or character variable names.", arg_name))
    }
  }

  resolve_pairs <- function(pairs, arg_name) {
    if (is.null(pairs)) {
      return(NULL)
    }
    if (!is.list(pairs)) stop(sprintf("'%s' must be a list of length-2 vectors.", arg_name))

    resolved <- lapply(seq_along(pairs), function(i) {
      p <- pairs[[i]]
      if (length(p) != 2) {
        stop(sprintf("Each element of '%s' must be a length-2 vector (from, to).", arg_name))
      }
      # Resolve each element individually even if types are mixed within a pair
      from_idx <- resolve_variable(p[1], sprintf("%s[[%d]][1]", arg_name, i))
      to_idx <- resolve_variable(p[2], sprintf("%s[[%d]][2]", arg_name, i))
      c(from_idx, to_idx)
    })
    return(resolved)
  }

  # --- Resolve arguments ---
  if (!is.null(exogenous_variables)) {
    exogenous_variables <- resolve_variable(exogenous_variables, "exogenous_variables")
  }
  if (!is.null(sink_variables)) {
    sink_variables <- resolve_variable(sink_variables, "sink_variables")
  }
  paths <- resolve_pairs(paths, "paths")
  no_paths <- resolve_pairs(no_paths, "no_paths")

  # --- Initialize the prior knowledge matrix (all -1 = unknown) ---
  prior_knowledge <- matrix(-1L, nrow = n_variables, ncol = n_variables)

  if (!is.null(labels)) {
    rownames(prior_knowledge) <- labels
    colnames(prior_knowledge) <- labels
  }

  # --- Set no_paths: (from, to) -> prior_knowledge[to, from] = 0 ---
  if (!is.null(no_paths)) {
    for (pair in no_paths) {
      prior_knowledge[pair[2], pair[1]] <- 0L
    }
  }

  # --- Set paths: (from, to) -> prior_knowledge[to, from] = 1 ---
  if (!is.null(paths)) {
    for (pair in paths) {
      prior_knowledge[pair[2], pair[1]] <- 1L
    }
  }

  # --- sink_variables: do not influence other variables -> set entire column to 0 ---
  if (!is.null(sink_variables)) {
    for (var in sink_variables) {
      prior_knowledge[, var] <- 0L
    }
  }

  # --- exogenous_variables: not influenced by other variables -> set entire row to 0 ---
  if (!is.null(exogenous_variables)) {
    for (var in exogenous_variables) {
      prior_knowledge[var, ] <- 0L
    }
  }

  # --- Diagonal elements are -1 ---
  diag(prior_knowledge) <- -1L

  return(prior_knowledge)
}
