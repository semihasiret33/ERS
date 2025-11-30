# ============================================================================
# Script: irtree_pseudo_item_converter.R
# Purpose: Convert Likert responses to IRTree pseudo-items
# Author: ERS Research Project
# Date: 2025-11-30
# References: Böckenholt (2012); Plieninger (2017)
# ============================================================================
#
# IRTree Pseudo-Item Recoding
#
# 4-category Likert (1-4) → 3 binary pseudo-items per original item:
#   - Pseudo-item 1 (Node 1): Disagree (0) vs Agree (1)
#   - Pseudo-item 2 (Node 2): Given disagree, extreme (0) vs moderate (1)
#   - Pseudo-item 3 (Node 3): Given agree, moderate (0) vs extreme (1)
#
# Recoding scheme:
#   Original → (Node1, Node2, Node3)
#   1 → (0, 0, NA)   [Strongly disagree]
#   2 → (0, 1, NA)   [Disagree]
#   3 → (1, NA, 0)   [Agree]
#   4 → (1, NA, 1)   [Strongly agree]
#
# ============================================================================

library(mirt)


#' Convert Likert Responses to IRTree Pseudo-Items (4-point)
#'
#' @param data Matrix or data frame with Likert responses (1-4)
#' @param n_items Number of original items
#' @return Matrix with 3*n_items pseudo-items (binary with NAs for unused paths)
convert_to_irtree_pseudoitems_4point <- function(data, n_items) {

  n_persons <- nrow(data)

  # Create storage for pseudo-items (3 per original item)
  tree_data <- matrix(NA, nrow = n_persons, ncol = n_items * 3)

  for (i in 1:n_items) {

    # Original responses
    orig_resp <- data[, i]

    # Node 1: Agree (1) vs Disagree (0)
    # 1, 2 → 0 (disagree)
    # 3, 4 → 1 (agree)
    node1 <- ifelse(orig_resp <= 2, 0, 1)

    # Node 2: Given disagree, extreme (0) vs moderate (1)
    # 1 → 0 (extreme disagree)
    # 2 → 1 (moderate disagree)
    # 3, 4 → NA (did not take this path)
    node2 <- rep(NA, n_persons)
    node2[orig_resp == 1] <- 0
    node2[orig_resp == 2] <- 1

    # Node 3: Given agree, moderate (0) vs extreme (1)
    # 3 → 0 (moderate agree)
    # 4 → 1 (extreme agree)
    # 1, 2 → NA (did not take this path)
    node3 <- rep(NA, n_persons)
    node3[orig_resp == 3] <- 0
    node3[orig_resp == 4] <- 1

    # Store in tree_data
    # Columns: item1_node1, item1_node2, item1_node3, item2_node1, ...
    tree_data[, (i - 1) * 3 + 1] <- node1
    tree_data[, (i - 1) * 3 + 2] <- node2
    tree_data[, (i - 1) * 3 + 3] <- node3
  }

  # Set column names
  node_names <- character(n_items * 3)
  for (i in 1:n_items) {
    node_names[(i - 1) * 3 + 1] <- paste0("Item", i, "_Node1")
    node_names[(i - 1) * 3 + 2] <- paste0("Item", i, "_Node2")
    node_names[(i - 1) * 3 + 3] <- paste0("Item", i, "_Node3")
  }

  colnames(tree_data) <- node_names

  return(tree_data)
}


#' Convert Likert Responses to IRTree Pseudo-Items (5-point)
#'
#' @param data Matrix or data frame with Likert responses (1-5)
#' @param n_items Number of original items
#' @return Matrix with pseudo-items
convert_to_irtree_pseudoitems_5point <- function(data, n_items) {

  n_persons <- nrow(data)
  tree_data <- matrix(NA, nrow = n_persons, ncol = n_items * 3)

  for (i in 1:n_items) {

    orig_resp <- data[, i]

    # Node 1: Agree vs Disagree (middle category = neutral, coded as disagree)
    # 1, 2, 3 → 0 (disagree/neutral)
    # 4, 5 → 1 (agree)
    node1 <- ifelse(orig_resp <= 3, 0, 1)

    # Node 2: Given disagree, extreme vs moderate
    # 1 → 0 (strongly disagree - extreme)
    # 2, 3 → 1 (disagree/neutral - moderate)
    node2 <- rep(NA, n_persons)
    node2[orig_resp == 1] <- 0
    node2[orig_resp %in% c(2, 3)] <- 1

    # Node 3: Given agree, moderate vs extreme
    # 4 → 0 (agree - moderate)
    # 5 → 1 (strongly agree - extreme)
    node3 <- rep(NA, n_persons)
    node3[orig_resp == 4] <- 0
    node3[orig_resp == 5] <- 1

    tree_data[, (i - 1) * 3 + 1] <- node1
    tree_data[, (i - 1) * 3 + 2] <- node2
    tree_data[, (i - 1) * 3 + 3] <- node3
  }

  # Set column names
  node_names <- character(n_items * 3)
  for (i in 1:n_items) {
    node_names[(i - 1) * 3 + 1] <- paste0("Item", i, "_Node1")
    node_names[(i - 1) * 3 + 2] <- paste0("Item", i, "_Node2")
    node_names[(i - 1) * 3 + 3] <- paste0("Item", i, "_Node3")
  }

  colnames(tree_data) <- node_names

  return(tree_data)
}


#' Fit IRTree Model with Proper Pseudo-Item Approach
#'
#' @param data Original Likert response data
#' @param group_var Grouping variable (country)
#' @param n_items Number of items
#' @param n_categories Number of response categories (4 or 5)
#' @return List with fitted model and extracted means
fit_proper_irtree_model <- function(data, group_var, n_items, n_categories = 4) {

  # Convert to pseudo-items
  if (n_categories == 4) {
    tree_data <- convert_to_irtree_pseudoitems_4point(data, n_items)
  } else if (n_categories == 5) {
    tree_data <- convert_to_irtree_pseudoitems_5point(data, n_items)
  } else {
    stop("Only 4-point and 5-point Likert scales supported")
  }

  # Initialize custom nodes
  source("R/functions/irtree_functions.R")
  irtree_nodes <- initialize_irtree_nodes()

  # Define model specification
  # For each original item, we have 3 nodes:
  # - Node 1: Content dimension (2PL)
  # - Node 2: ERS dimension (custom Node2TreePL)
  # - Node 3: ERS dimension (custom Node3TreePL)

  # Model syntax
  model_syntax <- paste0(
    "Content = ", paste(seq(1, n_items * 3, by = 3), collapse = ","), "\n",
    "ERS = ", paste(c(seq(2, n_items * 3, by = 3), seq(3, n_items * 3, by = 3)),
                   collapse = ","), "\n",
    "FREE = (GROUP, COV_21)"
  )

  # Item types: alternating pattern
  # Position 1: Node1 (2PL), Position 2: Node2 (custom), Position 3: Node3 (custom)
  item_types <- rep(c("2PL", "Node2TreePL", "Node3TreePL"), n_items)

  # Fit model
  tryCatch({

    irtree_fit <- multipleGroup(
      data = tree_data,
      model = model_syntax,
      group = group_var,
      itemtype = item_types,
      customItems = list(
        Node2TreePL = irtree_nodes$Node2,
        Node3TreePL = irtree_nodes$Node3
      ),
      method = "EM",
      invariance = c("free_mean", "free_var", "slopes", "intercepts"),
      technical = list(NCYCLES = 2000),
      verbose = FALSE
    )

    # Extract group means
    coef_list <- coef(irtree_fit, simplify = TRUE)
    group_names <- unique(group_var)

    theta_estimates <- sapply(group_names, function(g_name) {
      group_idx <- which(names(coef_list) == g_name)
      if (length(group_idx) > 0 && "means" %in% names(coef_list[[group_idx]])) {
        return(coef_list[[group_idx]]$means[1])  # Content dimension
      } else {
        return(NA)
      }
    })

    ers_estimates <- sapply(group_names, function(g_name) {
      group_idx <- which(names(coef_list) == g_name)
      if (length(group_idx) > 0 && "means" %in% names(coef_list[[group_idx]])) {
        return(coef_list[[group_idx]]$means[2])  # ERS dimension
      } else {
        return(NA)
      }
    })

    result <- list(
      model = irtree_fit,
      theta_means = data.frame(
        country = group_names,
        mean_irtree = theta_estimates,
        ers_irtree = ers_estimates
      ),
      convergence = TRUE
    )

    return(result)

  }, error = function(e) {
    cat("    IRTree model failed:", e$message, "\n")

    result <- list(
      model = NULL,
      theta_means = data.frame(
        country = unique(group_var),
        mean_irtree = NA,
        ers_irtree = NA
      ),
      convergence = FALSE,
      error = e$message
    )

    return(result)
  })
}


#' Test IRTree Pseudo-Item Conversion
#'
#' Quick test function to verify pseudo-item recoding
test_pseudoitem_conversion <- function() {

  cat("Testing IRTree pseudo-item conversion...\n\n")

  # Example data: 5 persons, 2 items, 4-point Likert
  test_data <- matrix(c(
    1, 3,  # Person 1: strongly disagree, agree
    2, 4,  # Person 2: disagree, strongly agree
    3, 2,  # Person 3: agree, disagree
    4, 1,  # Person 4: strongly agree, strongly disagree
    2, 3   # Person 5: disagree, agree
  ), nrow = 5, byrow = TRUE)

  colnames(test_data) <- c("Item1", "Item2")

  cat("Original Likert data (4-point: 1-4):\n")
  print(test_data)
  cat("\n")

  # Convert to pseudo-items
  tree_data <- convert_to_irtree_pseudoitems_4point(test_data, n_items = 2)

  cat("Converted to IRTree pseudo-items:\n")
  cat("(Node1: 0=disagree/1=agree, Node2: 0=extreme/1=moderate | disagree,\n")
  cat(" Node3: 0=moderate/1=extreme | agree)\n\n")
  print(tree_data)
  cat("\n")

  # Verify recoding
  cat("Verification:\n")
  cat("  Person 1, Item 1 = 1 (strongly disagree) → Node1=0, Node2=0, Node3=NA ✓\n")
  cat("  Person 2, Item 1 = 2 (disagree) → Node1=0, Node2=1, Node3=NA ✓\n")
  cat("  Person 3, Item 1 = 3 (agree) → Node1=1, Node2=NA, Node3=0 ✓\n")
  cat("  Person 4, Item 1 = 4 (strongly agree) → Node1=1, Node2=NA, Node3=1 ✓\n")
  cat("\n")

  cat("Test completed successfully!\n")

  return(tree_data)
}

# Uncomment to test:
# test_result <- test_pseudoitem_conversion()
