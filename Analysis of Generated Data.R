GPCM_dat_1 <- list(test_GPCM_average_10_1,
                   test_GPCM_average_10_2,
                   test_GPCM_average_10_3,
                   test_GPCM_average_10_4,
                   test_GPCM_average_10_5)

GPCM_dat_2 <- list(test_GPCM_difficult_10_1,
                   test_GPCM_difficult_10_2,
                   test_GPCM_difficult_10_3,
                   test_GPCM_difficult_10_4,
                   test_GPCM_difficult_10_5)

GPCM_dat_3 <- list(test_GPCM_average_20_1,
                   test_GPCM_average_20_2,
                   test_GPCM_average_20_3,
                   test_GPCM_average_20_4,
                   test_GPCM_average_20_5)

GPCM_dat_4 <- list(test_GPCM_difficult_20_1,
                   test_GPCM_difficult_20_2,
                   test_GPCM_difficult_20_3,
                   test_GPCM_difficult_20_4,
                   test_GPCM_difficult_20_5)

MNRM_dat <- list(test_MNRM_average_10_ERS_min1,
                 test_MNRM_difficult_10_ERS_min1,
                 test_MNRM_average_20_ERS_min1,
                 test_MNRM_difficult_20_ERS_min1,
                 test_MNRM_average_10,
                 test_MNRM_difficult_10,
                 test_MNRM_average_20,
                 test_MNRM_difficult_20,
                 test_MNRM_average_10_ERS_1,
                 test_MNRM_difficult_10_ERS_1,
                 test_MNRM_average_20_ERS_1,
                 test_MNRM_difficult_20_ERS_1)

Tree_dat <- list(test_Tree_average_10_ERS_min1,
                 test_Tree_difficult_10_ERS_min1,
                 test_Tree_average_20_ERS_min1,
                 test_Tree_difficult_20_ERS_min1,
                 test_Tree_average_10,
                 test_Tree_difficult_10,
                 test_Tree_average_20,
                 test_Tree_difficult_20,
                 test_Tree_average_10_ERS_1,
                 test_Tree_difficult_10_ERS_1,
                 test_Tree_average_20_ERS_1,
                 test_Tree_difficult_20_ERS_1)

reps <- 500

NA_matrix_GPCM_1 <- matrix(FALSE, ncol = 5, nrow = (reps/5))
for(i in 1:5){
  for(j in 1:(reps/5)){
    NA_matrix_GPCM_1[j, i] <- any(is.na(GPCM_dat_1[[i]][[j]]))
  }
}

sum(NA_matrix_GPCM_1)
which(NA_matrix_GPCM_1[,1])
which(NA_matrix_GPCM_1[,2])
which(NA_matrix_GPCM_1[,3])
which(NA_matrix_GPCM_1[,4])
which(NA_matrix_GPCM_1[,5])

NA_matrix_GPCM_2 <- matrix(FALSE, ncol = 5, nrow = (reps/5))
for(i in 1:5){
  for(j in 1:(reps/5)){
    NA_matrix_GPCM_2[j, i] <- any(is.na(GPCM_dat_2[[i]][[j]]))
  }
}

sum(NA_matrix_GPCM_2)
which(NA_matrix_GPCM_2[,1])
which(NA_matrix_GPCM_2[,2])
which(NA_matrix_GPCM_2[,3])
which(NA_matrix_GPCM_2[,4])
which(NA_matrix_GPCM_2[,5])

NA_matrix_GPCM_3 <- matrix(FALSE, ncol = 5, nrow = (reps/5))
for(i in 1:5){
  for(j in 1:(reps/5)){
    NA_matrix_GPCM_3[j, i] <- any(is.na(GPCM_dat_3[[i]][[j]]))
  }
}

sum(NA_matrix_GPCM_3)
which(NA_matrix_GPCM_3[,1])
which(NA_matrix_GPCM_3[,2])
which(NA_matrix_GPCM_3[,3])
which(NA_matrix_GPCM_3[,4])
which(NA_matrix_GPCM_3[,5])

NA_matrix_GPCM_4 <- matrix(FALSE, ncol = 5, nrow = (reps/5))
for(i in 1:5){
  for(j in 1:(reps/5)){
    NA_matrix_GPCM_4[j, i] <- any(is.na(GPCM_dat_4[[i]][[j]]))
  }
}

sum(NA_matrix_GPCM_4)
which(NA_matrix_GPCM_4[,1])
which(NA_matrix_GPCM_4[,2])
which(NA_matrix_GPCM_4[,3])
which(NA_matrix_GPCM_4[,4])
which(NA_matrix_GPCM_4[,5])

NA_matrix_MNRM <- matrix(FALSE, ncol = 12, nrow = reps)
for(i in 1:12){
  for(j in 1:reps){
    NA_matrix_MNRM[j, i] <- any(is.na(MNRM_dat[[i]][[j]]))
  }
}

sum(NA_matrix_MNRM)

NA_matrix_Tree <- matrix(FALSE, ncol = 12, nrow = reps)
for(i in 1:12){
  for(j in 1:reps){
    NA_matrix_Tree[j, i] <- any(is.na(Tree_dat[[i]][[j]]))
  }
}

sum(NA_matrix_Tree)

GPCM_dat_1 <- list(c(test_GPCM_average_10_1,
                     test_GPCM_average_10_2,
                     test_GPCM_average_10_3,
                     test_GPCM_average_10_4,
                     test_GPCM_average_10_5))

GPCM_dat_2 <- list(c(test_GPCM_difficult_10_1,
                     test_GPCM_difficult_10_2,
                     test_GPCM_difficult_10_3,
                     test_GPCM_difficult_10_4,
                     test_GPCM_difficult_10_5))

GPCM_dat_3 <- list(c(test_GPCM_average_20_1,
                     test_GPCM_average_20_2,
                     test_GPCM_average_20_3,
                     test_GPCM_average_20_4,
                     test_GPCM_average_20_5))

GPCM_dat_4 <- list(c(test_GPCM_difficult_20_1,
                     test_GPCM_difficult_20_2,
                     test_GPCM_difficult_20_3,
                     test_GPCM_difficult_20_4,
                     test_GPCM_difficult_20_5))

GPCM_dat <- list(GPCM_dat_1,
                 GPCM_dat_2,
                 GPCM_dat_3,
                 GPCM_dat_4)


control_table <- matrix(NA, nrow = 4, ncol = 6)
control_fit_table <- matrix(0, nrow = 4, ncol = 15)
MNRM_results <- Tree_results <- GPCM_results <- matrix(NA, nrow = reps, ncol = 2)
for(j in 1:4){
  for(i in 1:reps){
    MNRM_results[i, 1] <- GPCM_dat[[j]][[1]][[i]][["MNRM"]][["Group2pars"]][[length(GPCM_dat[[j]][[1]][[i]][["MNRM"]][["Group2pars"]])]][1]
    MNRM_results[i, 2] <- GPCM_dat[[j]][[1]][[i]][["MNRM"]][["Group2pars"]][[length(GPCM_dat[[j]][[1]][[i]][["MNRM"]][["Group2pars"]])]][3] - 1
    
    Tree_results[i, 1] <- GPCM_dat[[j]][[1]][[i]][["IRTree"]][["Group2pars"]][[length(GPCM_dat[[j]][[1]][[i]][["IRTree"]][["Group2pars"]])]][1]
    Tree_results[i, 2] <- GPCM_dat[[j]][[1]][[i]][["IRTree"]][["Group2pars"]][[length(GPCM_dat[[j]][[1]][[i]][["IRTree"]][["Group2pars"]])]][3] - 1
    
    GPCM_results[i, 1] <- GPCM_dat[[j]][[1]][[i]][["GPCM"]][["Group2pars"]][[length(GPCM_dat[[j]][[1]][[i]][["GPCM"]][["Group2pars"]])]][1]
    GPCM_results[i, 2] <- GPCM_dat[[j]][[1]][[i]][["GPCM"]][["Group2pars"]][[length(GPCM_dat[[j]][[1]][[i]][["GPCM"]][["Group2pars"]])]][2] - 1
    
    control_fit_table[j, which.min(c(GPCM_dat[[j]][[1]][[i]][["MNRM"]][["Fit"]][["AIC"]],
                                     GPCM_dat[[j]][[1]][[i]][["IRTree"]][["Fit"]][["AIC"]],
                                     GPCM_dat[[j]][[1]][[i]][["GPCM"]][["Fit"]][["AIC"]]))] <-
      control_fit_table[j, which.min(c(GPCM_dat[[j]][[1]][[i]][["MNRM"]][["Fit"]][["AIC"]],
                                       GPCM_dat[[j]][[1]][[i]][["IRTree"]][["Fit"]][["AIC"]],
                                       GPCM_dat[[j]][[1]][[i]][["GPCM"]][["Fit"]][["AIC"]]))] + 1
    
    control_fit_table[j, (3 + which.min(c(GPCM_dat[[j]][[1]][[i]][["MNRM"]][["Fit"]][["BIC"]],
                                          GPCM_dat[[j]][[1]][[i]][["IRTree"]][["Fit"]][["BIC"]],
                                          GPCM_dat[[j]][[1]][[i]][["GPCM"]][["Fit"]][["BIC"]])))] <-
      control_fit_table[j, (3 + which.min(c(GPCM_dat[[j]][[1]][[i]][["MNRM"]][["Fit"]][["BIC"]],
                                            GPCM_dat[[j]][[1]][[i]][["IRTree"]][["Fit"]][["BIC"]],
                                            GPCM_dat[[j]][[1]][[i]][["GPCM"]][["Fit"]][["BIC"]])))] + 1
    
    control_fit_table[j, (6 + which.min(c(GPCM_dat[[j]][[1]][[i]][["MNRM"]][["Fit"]][["SABIC"]],
                                          GPCM_dat[[j]][[1]][[i]][["IRTree"]][["Fit"]][["SABIC"]],
                                          GPCM_dat[[j]][[1]][[i]][["GPCM"]][["Fit"]][["SABIC"]])))] <-
      control_fit_table[j, (6 + which.min(c(GPCM_dat[[j]][[1]][[i]][["MNRM"]][["Fit"]][["SABIC"]],
                                            GPCM_dat[[j]][[1]][[i]][["IRTree"]][["Fit"]][["SABIC"]],
                                            GPCM_dat[[j]][[1]][[i]][["GPCM"]][["Fit"]][["SABIC"]])))] + 1
    
    control_fit_table[j, (9 + which.min(c(GPCM_dat[[j]][[1]][[i]][["MNRM"]][["Fit"]][["HQ"]],
                                          GPCM_dat[[j]][[1]][[i]][["IRTree"]][["Fit"]][["HQ"]],
                                          GPCM_dat[[j]][[1]][[i]][["GPCM"]][["Fit"]][["HQ"]])))] <-
      control_fit_table[j, (9 + which.min(c(GPCM_dat[[j]][[1]][[i]][["MNRM"]][["Fit"]][["HQ"]],
                                            GPCM_dat[[j]][[1]][[i]][["IRTree"]][["Fit"]][["HQ"]],
                                            GPCM_dat[[j]][[1]][[i]][["GPCM"]][["Fit"]][["HQ"]])))] + 1
    
    control_fit_table[j, (12 + which.max(c(GPCM_dat[[j]][[1]][[i]][["MNRM"]][["Fit"]][["logLik"]],
                                           GPCM_dat[[j]][[1]][[i]][["IRTree"]][["Fit"]][["logLik"]],
                                           GPCM_dat[[j]][[1]][[i]][["GPCM"]][["Fit"]][["logLik"]])))] <-
      control_fit_table[j, (12 + which.max(c(GPCM_dat[[j]][[1]][[i]][["MNRM"]][["Fit"]][["logLik"]],
                                             GPCM_dat[[j]][[1]][[i]][["IRTree"]][["Fit"]][["logLik"]],
                                             GPCM_dat[[j]][[1]][[i]][["GPCM"]][["Fit"]][["logLik"]])))] + 1
    
  }
  control_table[j,1] <- mean(MNRM_results[,1])
  control_table[j,2] <- mean(MNRM_results[,2])
  
  control_table[j,3] <- mean(Tree_results[,1])
  control_table[j,4] <- mean(Tree_results[,2])
  
  control_table[j,5] <- mean(GPCM_results[,1])
  control_table[j,6] <- mean(GPCM_results[,2])
}
write.csv(control_table, file = "Controltable.csv")
write.csv(control_fit_table, file = "Control_fit_table.csv")

condition_table <- matrix(NA, nrow = 24, ncol = 6)
MNRM_fit_table <- matrix(0, nrow = 12, ncol = 15)
MNRM_results <- Tree_results <- GPCM_results <- matrix(NA, nrow = reps, ncol = 2)
for(j in 1:12){
  for(i in 1:reps){
    MNRM_results[i, 1] <- MNRM_dat[[j]][[i]][["MNRM"]][["Group2pars"]][[length(MNRM_dat[[j]][[i]][["MNRM"]][["Group2pars"]])]][1]
    MNRM_results[i, 2] <- MNRM_dat[[j]][[i]][["MNRM"]][["Group2pars"]][[length(MNRM_dat[[j]][[i]][["MNRM"]][["Group2pars"]])]][3] - 1
    
    Tree_results[i, 1] <- MNRM_dat[[j]][[i]][["IRTree"]][["Group2pars"]][[length(MNRM_dat[[j]][[i]][["IRTree"]][["Group2pars"]])]][1]
    Tree_results[i, 2] <- MNRM_dat[[j]][[i]][["IRTree"]][["Group2pars"]][[length(MNRM_dat[[j]][[i]][["IRTree"]][["Group2pars"]])]][3] - 1
    
    GPCM_results[i, 1] <- MNRM_dat[[j]][[i]][["GPCM"]][["Group2pars"]][[length(MNRM_dat[[j]][[i]][["GPCM"]][["Group2pars"]])]][1]
    GPCM_results[i, 2] <- MNRM_dat[[j]][[i]][["GPCM"]][["Group2pars"]][[length(MNRM_dat[[j]][[i]][["GPCM"]][["Group2pars"]])]][2] - 1
    
    MNRM_fit_table[j, which.min(c(MNRM_dat[[j]][[i]][["MNRM"]][["Fit"]][["AIC"]],
                                  MNRM_dat[[j]][[i]][["IRTree"]][["Fit"]][["AIC"]],
                                  MNRM_dat[[j]][[i]][["GPCM"]][["Fit"]][["AIC"]]))] <-
      MNRM_fit_table[j, which.min(c(MNRM_dat[[j]][[i]][["MNRM"]][["Fit"]][["AIC"]],
                                    MNRM_dat[[j]][[i]][["IRTree"]][["Fit"]][["AIC"]],
                                    MNRM_dat[[j]][[i]][["GPCM"]][["Fit"]][["AIC"]]))] + 1
    
    MNRM_fit_table[j, (3 + which.min(c(MNRM_dat[[j]][[i]][["MNRM"]][["Fit"]][["BIC"]],
                                       MNRM_dat[[j]][[i]][["IRTree"]][["Fit"]][["BIC"]],
                                       MNRM_dat[[j]][[i]][["GPCM"]][["Fit"]][["BIC"]])))] <-
      MNRM_fit_table[j, (3 + which.min(c(MNRM_dat[[j]][[i]][["MNRM"]][["Fit"]][["BIC"]],
                                         MNRM_dat[[j]][[i]][["IRTree"]][["Fit"]][["BIC"]],
                                         MNRM_dat[[j]][[i]][["GPCM"]][["Fit"]][["BIC"]])))] + 1
    
    MNRM_fit_table[j, (6 + which.min(c(MNRM_dat[[j]][[i]][["MNRM"]][["Fit"]][["SABIC"]],
                                       MNRM_dat[[j]][[i]][["IRTree"]][["Fit"]][["SABIC"]],
                                       MNRM_dat[[j]][[i]][["GPCM"]][["Fit"]][["SABIC"]])))] <-
      MNRM_fit_table[j, (6 + which.min(c(MNRM_dat[[j]][[i]][["MNRM"]][["Fit"]][["SABIC"]],
                                         MNRM_dat[[j]][[i]][["IRTree"]][["Fit"]][["SABIC"]],
                                         MNRM_dat[[j]][[i]][["GPCM"]][["Fit"]][["SABIC"]])))] + 1
    
    MNRM_fit_table[j, (9 + which.min(c(MNRM_dat[[j]][[i]][["MNRM"]][["Fit"]][["HQ"]],
                                       MNRM_dat[[j]][[i]][["IRTree"]][["Fit"]][["HQ"]],
                                       MNRM_dat[[j]][[i]][["GPCM"]][["Fit"]][["HQ"]])))] <-
      MNRM_fit_table[j, (9 + which.min(c(MNRM_dat[[j]][[i]][["MNRM"]][["Fit"]][["HQ"]],
                                         MNRM_dat[[j]][[i]][["IRTree"]][["Fit"]][["HQ"]],
                                         MNRM_dat[[j]][[i]][["GPCM"]][["Fit"]][["HQ"]])))] + 1
    
    MNRM_fit_table[j, (12 + which.max(c(MNRM_dat[[j]][[i]][["MNRM"]][["Fit"]][["logLik"]],
                                        MNRM_dat[[j]][[i]][["IRTree"]][["Fit"]][["logLik"]],
                                        MNRM_dat[[j]][[i]][["GPCM"]][["Fit"]][["logLik"]])))] <-
      MNRM_fit_table[j, (12 + which.max(c(MNRM_dat[[j]][[i]][["MNRM"]][["Fit"]][["logLik"]],
                                          MNRM_dat[[j]][[i]][["IRTree"]][["Fit"]][["logLik"]],
                                          MNRM_dat[[j]][[i]][["GPCM"]][["Fit"]][["logLik"]])))] + 1
  }
  condition_table[j,1] <- mean(MNRM_results[,1])
  condition_table[j,2] <- mean(MNRM_results[,2])
  
  condition_table[j,3] <- mean(Tree_results[,1])
  condition_table[j,4] <- mean(Tree_results[,2])
  
  condition_table[j,5] <- mean(GPCM_results[,1])
  condition_table[j,6] <- mean(GPCM_results[,2])
}

write.csv(MNRM_fit_table, file = "MNRM_fit_table.csv")

Tree_fit_table <- matrix(0, nrow = 12, ncol = 15)
for(j in 1:12){
  for(i in 1:reps){
    MNRM_results[i, 1] <- Tree_dat[[j]][[i]][["MNRM"]][["Group2pars"]][[length(Tree_dat[[j]][[i]][["MNRM"]][["Group2pars"]])]][1]
    MNRM_results[i, 2] <- Tree_dat[[j]][[i]][["MNRM"]][["Group2pars"]][[length(Tree_dat[[j]][[i]][["MNRM"]][["Group2pars"]])]][3] - 1
    
    Tree_results[i, 1] <- Tree_dat[[j]][[i]][["IRTree"]][["Group2pars"]][[length(Tree_dat[[j]][[i]][["IRTree"]][["Group2pars"]])]][1]
    Tree_results[i, 2] <- Tree_dat[[j]][[i]][["IRTree"]][["Group2pars"]][[length(Tree_dat[[j]][[i]][["IRTree"]][["Group2pars"]])]][3] - 1
    
    GPCM_results[i, 1] <- Tree_dat[[j]][[i]][["GPCM"]][["Group2pars"]][[length(Tree_dat[[j]][[i]][["GPCM"]][["Group2pars"]])]][1]
    GPCM_results[i, 2] <- Tree_dat[[j]][[i]][["GPCM"]][["Group2pars"]][[length(Tree_dat[[j]][[i]][["GPCM"]][["Group2pars"]])]][2] - 1
    
    
    Tree_fit_table[j, which.min(c(Tree_dat[[j]][[i]][["MNRM"]][["Fit"]][["AIC"]],
                                  Tree_dat[[j]][[i]][["IRTree"]][["Fit"]][["AIC"]],
                                  Tree_dat[[j]][[i]][["GPCM"]][["Fit"]][["AIC"]]))] <-
      Tree_fit_table[j, which.min(c(Tree_dat[[j]][[i]][["MNRM"]][["Fit"]][["AIC"]],
                                    Tree_dat[[j]][[i]][["IRTree"]][["Fit"]][["AIC"]],
                                    Tree_dat[[j]][[i]][["GPCM"]][["Fit"]][["AIC"]]))] + 1
    
    Tree_fit_table[j, (3 + which.min(c(Tree_dat[[j]][[i]][["MNRM"]][["Fit"]][["BIC"]],
                                       Tree_dat[[j]][[i]][["IRTree"]][["Fit"]][["BIC"]],
                                       Tree_dat[[j]][[i]][["GPCM"]][["Fit"]][["BIC"]])))] <-
      Tree_fit_table[j, (3 + which.min(c(Tree_dat[[j]][[i]][["MNRM"]][["Fit"]][["BIC"]],
                                         Tree_dat[[j]][[i]][["IRTree"]][["Fit"]][["BIC"]],
                                         Tree_dat[[j]][[i]][["GPCM"]][["Fit"]][["BIC"]])))] + 1
    
    Tree_fit_table[j, (6 + which.min(c(Tree_dat[[j]][[i]][["MNRM"]][["Fit"]][["SABIC"]],
                                       Tree_dat[[j]][[i]][["IRTree"]][["Fit"]][["SABIC"]],
                                       Tree_dat[[j]][[i]][["GPCM"]][["Fit"]][["SABIC"]])))] <-
      Tree_fit_table[j, (6 + which.min(c(Tree_dat[[j]][[i]][["MNRM"]][["Fit"]][["SABIC"]],
                                         Tree_dat[[j]][[i]][["IRTree"]][["Fit"]][["SABIC"]],
                                         Tree_dat[[j]][[i]][["GPCM"]][["Fit"]][["SABIC"]])))] + 1
    
    Tree_fit_table[j, (9 + which.min(c(Tree_dat[[j]][[i]][["MNRM"]][["Fit"]][["HQ"]],
                                       Tree_dat[[j]][[i]][["IRTree"]][["Fit"]][["HQ"]],
                                       Tree_dat[[j]][[i]][["GPCM"]][["Fit"]][["HQ"]])))] <-
      Tree_fit_table[j, (9 + which.min(c(Tree_dat[[j]][[i]][["MNRM"]][["Fit"]][["HQ"]],
                                         Tree_dat[[j]][[i]][["IRTree"]][["Fit"]][["HQ"]],
                                         Tree_dat[[j]][[i]][["GPCM"]][["Fit"]][["HQ"]])))] + 1
    
    Tree_fit_table[j, (12 + which.max(c(Tree_dat[[j]][[i]][["MNRM"]][["Fit"]][["logLik"]],
                                        Tree_dat[[j]][[i]][["IRTree"]][["Fit"]][["logLik"]],
                                        Tree_dat[[j]][[i]][["GPCM"]][["Fit"]][["logLik"]])))] <-
      Tree_fit_table[j, (12 + which.max(c(Tree_dat[[j]][[i]][["MNRM"]][["Fit"]][["logLik"]],
                                          Tree_dat[[j]][[i]][["IRTree"]][["Fit"]][["logLik"]],
                                          Tree_dat[[j]][[i]][["GPCM"]][["Fit"]][["logLik"]])))] + 1
  }
  condition_table[j + 12,1] <- mean(MNRM_results[,1])
  condition_table[j + 12,2] <- mean(MNRM_results[,2])
  
  condition_table[j + 12,3] <- mean(Tree_results[,1])
  condition_table[j + 12,4] <- mean(Tree_results[,2])
  
  condition_table[j + 12,5] <- mean(GPCM_results[,1])
  condition_table[j + 12,6] <- mean(GPCM_results[,2])
}

write.csv(Tree_fit_table, file = "Tree_fit_table.csv")
write.csv(condition_table, file = "Conditiontable.csv")

itempars_MNRM <- vector(mode = "list", length = 6)
for(i in 1:12){
  itempars_MNRM[[i]] <- matrix(NA, nrow = length(MNRM_dat[[i]][[1]][["MNRM"]][["Group2pars"]]) - 1, ncol = 6)
  for(j in 1:reps){
    for(k in 1:(length(MNRM_dat[[i]][[j]][["MNRM"]][["Group2pars"]]) - 1)){
      if(j == 1){
        itempars_MNRM[[i]][k,] <- MNRM_dat[[i]][[j]][["MNRM"]][["Group1pars"]][[k]][c(1,2,11,12,13,14)]
      } else {
        itempars_MNRM[[i]][k,] <- itempars_MNRM[[i]][k,] + MNRM_dat[[i]][[j]][["MNRM"]][["Group1pars"]][[k]][c(1,2,11,12,13,14)]
      }
      
    }
  }
  itempars_MNRM[[i]] <- itempars_MNRM[[i]] / reps
}


###############################################################################
#Quick check
################################################################################
node1 <- node2 <- node3 <- matrix(NA, ncol = 3, nrow = 1)
node1_track <- node2_track <- node3_track <- matrix(NA, ncol = 3, nrow = 500 * 10)

for(i in 1:reps){
  for(k in 1:10){
    node1_track[((i-1) * 10 + k),] <- test_Tree_difficult_20[[i]][["IRTree"]][["Group1pars"]][[((k-1) * 3 + 1)]][1:3]
    node2_track[((i-1) * 10 + k),] <- test_Tree_difficult_20[[i]][["IRTree"]][["Group1pars"]][[((k-1) * 3 + 2)]]
    node3_track[((i-1) * 10 + k),] <- test_Tree_difficult_20[[i]][["IRTree"]][["Group1pars"]][[((k-1) * 3 + 3)]]
  }
}
node1 <- colMeans(node1_track)
node2 <- colMeans(node2_track)
node3 <- colMeans(node3_track)

node1
node2
node3
