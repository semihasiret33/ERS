# ============================================================================
# Script: ers_indices.R
# Purpose: Functions to calculate various ERS (Extreme Response Style) indices
# Author: Adapted for ERS Research Project
# Date: 2025-11-28
# ============================================================================

#' Calculate Greenleaf ERS Index
#'
#' Computes the Greenleaf Extreme Response Style index for each respondent
#' based on their standard deviation across Likert items. Higher values
#' indicate more extreme responding (using endpoints of scale).
#'
#' @param data Data frame containing response data
#' @param items Character vector of item column names to use for calculation
#' @param na_rm Logical, should NA values be removed? Default is TRUE
#' @return Numeric vector of ERS index values (one per respondent)
#' @references Greenleaf, E. A. (1992). Improving rating scale measures by
#'   controlling for extreme response style. Journal of Marketing Research,
#'   29(2), 176-185.
#' @examples
#' \dontrun{
#' ers_scores <- calculate_greenleaf_ers(pisa_data, c("ST001Q01", "ST001Q02"))
#' }
calculate_greenleaf_ers <- function(data, items, na_rm = TRUE) {

  # Validate inputs
  if (!all(items %in% names(data))) {
    stop("Not all specified items found in data")
  }

  # Select item columns
  item_data <- data[, items, drop = FALSE]

  # Calculate standard deviation across items for each respondent
  ers_index <- apply(item_data, 1, sd, na.rm = na_rm)

  return(ers_index)
}


#' Calculate Alternative ERS Index (Count-based)
#'
#' Counts the proportion of extreme responses (1 or max category) for each
#' respondent. This is an alternative to the Greenleaf SD-based index.
#'
#' @param data Data frame containing response data
#' @param items Character vector of item column names
#' @param max_category Integer, maximum category value (e.g., 4 for 1-4 scale)
#' @param min_category Integer, minimum category value (default = 1)
#' @param na_rm Logical, should NA values be removed? Default is TRUE
#' @return Numeric vector of proportion of extreme responses per respondent
#' @examples
#' \dontrun{
#' ers_count <- calculate_ers_count(data, items = c("Q1", "Q2"), max_category = 4)
#' }
calculate_ers_count <- function(data, items, max_category,
                                min_category = 1, na_rm = TRUE) {

  # Select item columns
  item_data <- data[, items, drop = FALSE]

  # Count extreme responses (min or max category)
  is_extreme <- (item_data == min_category | item_data == max_category)

  # Calculate proportion of extreme responses per person
  if (na_rm) {
    ers_count <- rowSums(is_extreme, na.rm = TRUE) /
      rowSums(!is.na(item_data))
  } else {
    ers_count <- rowMeans(is_extreme)
  }

  return(ers_count)
}


#' Calculate Midpoint Response Style Index
#'
#' Calculates the proportion of midpoint responses for each respondent.
#' Useful as a complementary index to ERS.
#'
#' @param data Data frame containing response data
#' @param items Character vector of item column names
#' @param max_category Integer, maximum category value
#' @param min_category Integer, minimum category value (default = 1)
#' @param na_rm Logical, should NA values be removed? Default is TRUE
#' @return Numeric vector of proportion of midpoint responses per respondent
#' @examples
#' \dontrun{
#' mrs <- calculate_midpoint_rs(data, items = c("Q1", "Q2"), max_category = 4)
#' }
calculate_midpoint_rs <- function(data, items, max_category,
                                  min_category = 1, na_rm = TRUE) {

  # Calculate midpoint(s)
  range_vals <- min_category:max_category
  n_cats <- length(range_vals)

  if (n_cats %% 2 == 0) {
    # Even number of categories - two middle categories
    midpoints <- range_vals[c(n_cats/2, n_cats/2 + 1)]
  } else {
    # Odd number of categories - one middle category
    midpoints <- range_vals[ceiling(n_cats/2)]
  }

  # Select item columns
  item_data <- data[, items, drop = FALSE]

  # Count midpoint responses
  is_midpoint <- sapply(midpoints, function(mp) item_data == mp)
  if (length(midpoints) > 1) {
    is_midpoint <- Reduce("|",
                          lapply(1:length(midpoints), function(i) is_midpoint[, seq(i, ncol(is_midpoint), by = length(midpoints))]))
  }

  # Calculate proportion
  if (na_rm) {
    mrs <- rowSums(is_midpoint, na.rm = TRUE) / rowSums(!is.na(item_data))
  } else {
    mrs <- rowMeans(is_midpoint)
  }

  return(mrs)
}


#' Add ERS Index to Dataset
#'
#' Convenience function to add Greenleaf ERS index as a new column to a dataset
#'
#' @param data Data frame
#' @param items Character vector of item column names
#' @param index_name Name for the new ERS index column (default = "ers_index")
#' @param na_rm Logical, should NA values be removed? Default is TRUE
#' @return Data frame with added ERS index column
#' @examples
#' \dontrun{
#' pisa_data <- add_ers_index(pisa_data, items = paste0("MATHQ", 1:9))
#' }
add_ers_index <- function(data, items, index_name = "ers_index", na_rm = TRUE) {

  data[[index_name]] <- calculate_greenleaf_ers(data, items, na_rm)

  return(data)
}
