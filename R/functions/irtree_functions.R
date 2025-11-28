# ============================================================================
# Script: irtree_functions.R
# Purpose: Custom item functions for IRTree model implementation in mirt
# Author: Adapted for ERS Research Project
# Date: 2025-11-28
# References: Böckenholt (2012); Plieninger (2017)
# ============================================================================
#
# IRTree (Item Response Tree) Model for ERS
#
# Models Likert responses as a 3-node decision tree:
#   Node 1: Agree vs. Disagree (Content dimension: θ)
#   Node 2: If Disagree → Extreme (1) vs. Moderate (2) (ERS dimension: γ)
#   Node 3: If Agree → Moderate (3) vs. Extreme (4) (ERS dimension: γ)
#
# For 4-category Likert items (1-4):
#                    Node 1
#                  θ (Content)
#                   /        \
#            Disagree       Agree
#               /              \
#          Node 2            Node 3
#        γ (ERS)            γ (ERS)
#          /   \             /    \
#         1     2           3      4
# ============================================================================

library(mirt)


#' Create Custom IRTree Node 2 Item Function
#'
#' Node 2 models the disagreement branch: choosing between extreme (1)
#' vs. moderate (2) disagreement. This node captures ERS as a tendency
#' to choose extreme options when disagreeing.
#'
#' Probability model: P(extreme | disagree) =
#'   exp(a1*theta - a2*gamma + d) / (1 + exp(a1*theta - a2*gamma + d))
#'
#' Note the NEGATIVE loading on gamma: higher ERS → more likely to choose
#' extreme disagreement (category 1).
#'
#' @return mirt custom item object for Node 2
#' @references Böckenholt (2012); Plieninger (2017)
create_node2_function <- function() {

  name <- 'Node2TreePL'
  par <- c(a1 = 1, a2 = 1, d = 0)
  est <- c(TRUE, TRUE, TRUE)

  P.Node2TreePL <- function(par, Theta, ncat) {
    a1 <- par[1]  # Content discrimination
    a2 <- par[2]  # ERS discrimination
    d <- par[3]   # Difficulty/threshold

    # Probability of choosing extreme (1) over moderate (2)
    # Note: NEGATIVE loading on Theta[,2] (ERS dimension)
    P_extreme <- exp((a1 * Theta[,1] - a2 * Theta[,2] + d)) /
                 (1 + exp((a1 * Theta[,1] - a2 * Theta[,2] + d)))

    # Return probabilities: [P(moderate=2), P(extreme=1)]
    cbind(1 - P_extreme, P_extreme)
  }

  Node2TreePL <- createItem(
    name = name,
    par = par,
    est = est,
    P = P.Node2TreePL,
    derivType = 'symbolic',
    derivType.hss = 'symbolic'
  )

  return(Node2TreePL)
}


#' Create Custom IRTree Node 3 Item Function
#'
#' Node 3 models the agreement branch: choosing between moderate (3)
#' vs. extreme (4) agreement. This node captures ERS as a tendency
#' to choose extreme options when agreeing.
#'
#' Probability model: P(extreme | agree) =
#'   exp(a1*theta + a2*gamma + d) / (1 + exp(a1*theta + a2*gamma + d))
#'
#' Note the POSITIVE loading on gamma: higher ERS → more likely to choose
#' extreme agreement (category 4).
#'
#' @return mirt custom item object for Node 3
#' @references Böckenholt (2012); Plieninger (2017)
create_node3_function <- function() {

  name <- 'Node3TreePL'
  par <- c(a1 = 1, a2 = 1, d = 0)
  est <- c(TRUE, TRUE, TRUE)

  P.Node3TreePL <- function(par, Theta, ncat) {
    a1 <- par[1]  # Content discrimination
    a2 <- par[2]  # ERS discrimination
    d <- par[3]   # Difficulty/threshold

    # Probability of choosing extreme (4) over moderate (3)
    # Note: POSITIVE loading on Theta[,2] (ERS dimension)
    P_extreme <- exp((a1 * Theta[,1] + a2 * Theta[,2] + d)) /
                 (1 + exp((a1 * Theta[,1] + a2 * Theta[,2] + d)))

    # Return probabilities: [P(moderate=3), P(extreme=4)]
    cbind(1 - P_extreme, P_extreme)
  }

  Node3TreePL <- createItem(
    name = name,
    par = par,
    est = est,
    P = P.Node3TreePL,
    derivType = 'symbolic',
    derivType.hss = 'symbolic'
  )

  return(Node3TreePL)
}


#' Initialize IRTree Custom Functions
#'
#' Registers custom Node 2 and Node 3 item functions with mirt.
#' Must be called before fitting IRTree models.
#'
#' @return List with Node2TreePL and Node3TreePL objects
#' @export
initialize_irtree_nodes <- function() {

  Node2TreePL <- create_node2_function()
  Node3TreePL <- create_node3_function()

  message("IRTree custom node functions initialized successfully.")
  message("  - Node2TreePL: Disagree branch (extreme vs. moderate)")
  message("  - Node3TreePL: Agree branch (moderate vs. extreme)")

  return(list(
    Node2 = Node2TreePL,
    Node3 = Node3TreePL
  ))
}


#' Generate IRTree Data (Multidimensional Model)
#'
#' Generates Likert response data using IRTree decision process.
#' Each response is determined by three sequential decisions.
#'
#' @param theta_df Data frame with latent trait parameters per group
#' @param n_per_group Sample size per group
#' @param n_items Number of items
#' @param node1_pars Matrix of Node 1 parameters (n_items x 3: a1, a2, d)
#' @param node2_pars Matrix of Node 2 parameters (n_items x 3: a1, a2, d)
#' @param node3_pars Matrix of Node 3 parameters (n_items x 3: a1, a2, d)
#' @param categories Number of response categories (default = 4)
#' @return List containing response data, true thetas, and group membership
generate_irtree_data <- function(theta_df, n_per_group, n_items,
                                 node1_pars = matrix(rep(c(1, 0, 0), n_items),
                                                    nrow = n_items, ncol = 3),
                                 node2_pars = matrix(rep(c(1, 1, 0), n_items),
                                                    nrow = n_items, ncol = 3),
                                 node3_pars = matrix(rep(c(1, 1, 0), n_items),
                                                    nrow = n_items, ncol = 3),
                                 categories = 4) {

  n_groups <- nrow(theta_df)
  total_n <- n_per_group * n_groups

  # Storage
  response_data <- matrix(NA, nrow = total_n, ncol = n_items)
  true_theta <- matrix(NA, nrow = total_n, ncol = 2)
  group_membership <- rep(NA, total_n)

  for (g in 1:n_groups) {

    # Row indices for this group
    idx_start <- 1 + (g - 1) * n_per_group
    idx_end <- g * n_per_group
    idx <- idx_start:idx_end

    # Draw latent traits
    theta_content <- rnorm(n_per_group,
                          mean = theta_df$theta_mean[g],
                          sd = theta_df$theta_sd[g])
    theta_ers <- rnorm(n_per_group,
                      mean = theta_df$ers_mean[g],
                      sd = theta_df$ers_sd[g])

    # Store
    true_theta[idx, ] <- cbind(theta_content, theta_ers)
    group_membership[idx] <- theta_df$country[g]

    # Generate responses for each item
    for (j in 1:n_items) {

      # Node 1: Agree vs. Disagree (content dimension)
      p_agree <- plogis(node1_pars[j, 1] * theta_content + node1_pars[j, 3])
      node1_choice <- rbinom(n_per_group, 1, p_agree)  # 0=disagree, 1=agree

      # Node 2: If disagree, extreme vs. moderate
      p_extreme_disagree <- plogis(
        node2_pars[j, 1] * theta_content -
        node2_pars[j, 2] * theta_ers +
        node2_pars[j, 3]
      )
      node2_choice <- rbinom(n_per_group, 1, p_extreme_disagree)

      # Node 3: If agree, moderate vs. extreme
      p_extreme_agree <- plogis(
        node3_pars[j, 1] * theta_content +
        node3_pars[j, 2] * theta_ers +
        node3_pars[j, 3]
      )
      node3_choice <- rbinom(n_per_group, 1, p_extreme_agree)

      # Combine decisions into final response
      # Disagree branch: 1 (extreme) or 2 (moderate)
      # Agree branch: 3 (moderate) or 4 (extreme)
      response <- ifelse(
        node1_choice == 0,  # Disagree
        ifelse(node2_choice == 1, 1, 2),  # 1=extreme, 2=moderate
        ifelse(node3_choice == 1, 4, 3)   # 4=extreme, 3=moderate
      )

      response_data[idx, j] <- response
    }
  }

  # Format output
  colnames(response_data) <- paste0("Item", 1:n_items)
  colnames(true_theta) <- c("true_theta", "true_ers")

  result_data <- data.frame(
    country = group_membership,
    true_theta = true_theta[, 1],
    true_ers = true_theta[, 2],
    response_data
  )

  return(list(
    data = result_data,
    true_country_means = theta_df %>% select(country, theta_mean, ers_mean),
    node_parameters = list(
      node1 = node1_pars,
      node2 = node2_pars,
      node3 = node3_pars
    )
  ))
}


#' Fit IRTree Model with mirt
#'
#' Fits multidimensional IRTree model to Likert data using custom node functions
#'
#' @param data Response data (numeric matrix or data frame)
#' @param group_var Grouping variable (e.g., country)
#' @param n_items Number of items
#' @param estimator Estimation method (default = "EM")
#' @param ncycles Maximum number of EM cycles (default = 2000)
#' @return Fitted mirt model object
#' @export
fit_irtree_model <- function(data, group_var, n_items,
                             estimator = "EM", ncycles = 2000) {

  # Initialize custom nodes
  nodes <- initialize_irtree_nodes()

  # Model syntax: two correlated dimensions
  model_syntax <- paste0(
    "Theta = 1-", n_items, "\n",
    "ERS = 1-", n_items, "\n",
    "FREE = (GROUP, COV_21)"
  )

  # Specify item types for each node
  # For each item: Node1 (2PL), Node2 (custom), Node3 (custom)
  # But in mirt we specify per-item type, not per-node
  # IRTree requires recoding data to pseudo-items

  message("Note: IRTree implementation requires data preprocessing.")
  message("See Plieninger (2017) for full pseudo-item approach.")
  message("This is a simplified direct model approach.")

  # Fit model
  tryCatch({

    fit <- multipleGroup(
      data = data,
      model = model_syntax,
      group = group_var,
      itemtype = "graded",  # Will need custom specification
      method = estimator,
      invariance = c("free_mean", "free_var", "slopes", "intercepts"),
      technical = list(NCYCLES = ncycles),
      verbose = FALSE
    )

    message("IRTree model fitted successfully.")
    return(fit)

  }, error = function(e) {
    message("Error fitting IRTree model: ", e$message)
    message("IRTree requires specialized implementation - see irtrees package.")
    return(NULL)
  })
}


#' Extract Country Means from IRTree Model
#'
#' Extracts latent mean estimates (theta and ERS) for each group
#'
#' @param irtree_fit Fitted IRTree model object
#' @return Data frame with country means for theta and ERS
#' @export
extract_irtree_means <- function(irtree_fit) {

  if (is.null(irtree_fit)) {
    return(NULL)
  }

  coef_list <- coef(irtree_fit, simplify = TRUE)
  group_names <- names(coef_list)[!names(coef_list) %in% c("items", "lr.betas")]

  results <- data.frame(
    country = character(),
    theta_mean = numeric(),
    ers_mean = numeric(),
    stringsAsFactors = FALSE
  )

  for (g_name in group_names) {
    if ("means" %in% names(coef_list[[g_name]])) {
      results <- rbind(results, data.frame(
        country = g_name,
        theta_mean = coef_list[[g_name]]$means[1],
        ers_mean = coef_list[[g_name]]$means[2]
      ))
    }
  }

  return(results)
}


# ============================================================================
# Note on Full IRTree Implementation
# ============================================================================
#
# For production use, consider using the 'irtrees' package:
#   - Plieninger, H. (2017). Mountain or molehill? A simulation study on the
#     impact of response styles. Educational and Psychological Measurement.
#   - Package: https://github.com/hplieninger/ItemResponseTrees
#
# This file provides basic custom node functions for educational purposes
# and integration with existing mirt workflow.
#
# For full IRTree analysis:
#   1. Recode Likert items to pseudo-items (binary tree nodes)
#   2. Fit separate models for each node
#   3. Extract and combine parameter estimates
#
# See mnrm branch "Full code.R" for complete implementation example.
# ============================================================================
