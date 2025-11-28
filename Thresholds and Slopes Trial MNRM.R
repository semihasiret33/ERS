if(!require(mirt)){install.packages("mirt")}
if(!require(ItemResponseTrees)){install.packages("ItemResponseTrees")}
if(!require(tidyverse)){install.packages("tidyverse")}
if(!require(mixRasch)){install.packages("mixRasch")}
if(!require(future)){install.packages("future")}
if(!require(listenv)){install.packages("listenv")}
if(!require(simsalapar)){install.packages("simsalapar")}
if(!require(irtrees)){install.packages("irtrees")}
if(!require(reshape)){install.packages("reshape")}
if(!require(microbenchmark)){install.packages("microbenchmark")}
if(!require(future.apply)){install.packages("future.apply")}
if(!require(TAM)){install.packages("TAM")}
if(!require(lavaan)){install.packages("lavaan")}
if(!require(PCMRS)){install.packages("PCMRS")}
if(!require(shiny)){install.packages("shiny")}

generate <- function(gen_model = "MNRM", Nitems = 10, alpha = c(1,1), 
                     intercept = c(-2, -0.667, -0.667, -2), N = 100, 
                     categories = 4, groups = 1, theta_df,
                     estimator = "EM", node1pars = c(1, 0, 0), 
                     node2pars = c(-1, 1, 1.5), node3pars = c(1, 1, 1.5)){
  
  #Determine which groups should have their latent means and variances freely
  #estimated. For now, every group except group 1.
  Freegroups <- 2:groups
  
  #Generate matrices to store observed answers and group membership in
  obs_mat <- matrix(rep(NA, N*groups*Nitems), ncol = Nitems)
  obs_mat <- cbind(obs_mat, rep(NA, N*groups))
  
  #Confirmatory MNRM by Falk and Cai (2016) data generation. Equal to PCM if 
  #alpha is 1 and s_gen = c(0,1,2,...).
  if(gen_model == "MNRM-Falk"){
    
    #Generate s matrix
    s_gen <- matrix(0:(categories - 1), nrow = 1)
    
    #Generate ERS scoring matrix
    ers_mat <- matrix(rep(0, categories), nrow = 1)
    ers_mat[,1] <- ers_mat[,length(ers_mat)] <- 1
    
    #Combine ERS and ORS matrices
    s_gen <- rbind(s_gen, ers_mat)
    
    #initialize category probability, observed answer, newalpha and theta 
    #matrices
    newalpha <- matrix(rep(NA, categories * nrow(s_gen)), ncol = categories)
    cat_prob <- matrix(rep(NA, categories*N), ncol = categories)
    
    theta_store <- matrix(NA, ncol = 2, nrow = N*groups)
    for(g in 1:groups){
      
      #Draw participant ability from normal distribution with means and standard
      #deviations as specified in theta_df
      theta <- rnorm(N, mean = theta_df$theta_mean[g], 
                     sd = theta_df$theta_sd[g])
      
      #Draw participant ERS score from normal distribution with means and sd's
      #as specified in theta_df
      theta_ersmrs <- rnorm(N, mean = theta_df$ers_mean[g], 
                            sd = theta_df$ers_sd[g])
      
      #Combine into single matrix
      theta <- cbind(theta, theta_ersmrs)
      
      
      for(j in 1:Nitems){
        #Calculate the new item slope parameters based on the a priori s matrices
        newalpha <- alpha * s_gen
        
        #Calculate category probabilities
        numerator <- exp(apply(theta, 1, "%*%", newalpha) + intercept)
        denominator <- colSums(numerator)
        cat_prob <- t(numerator)/denominator
        
        #Convert probabilities into observed answers
        cat_prob_list <- split(cat_prob, row(cat_prob))
        obs_ans <- sapply(cat_prob_list, 
                          function(x)sample(1:categories, 1, replace = TRUE, prob = x),
                          simplify = TRUE)
        
        #Assign observed answers to observed answers matrix
        obs_mat[(1 + (g-1) * N):(N * g), j] <- obs_ans
      }
      #Store group membership in obs_mat matrix
      obs_mat[(1 + (g-1) * N):(N * g), Nitems + 1] <- g
      
      theta_store[(1 + (g-1)*N):(N*g),] <- theta
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #Data generation under IRTree multidimensional model 
    #(B?ckenholt & Meiser, 2017)
  } else if (gen_model == "IRTree-M"){
    for(g in 1:groups){
      
      #Draw substantive and ERS scores from normal distribution with mean and
      #sd specified in theta_df
      theta <- rnorm(N, mean = theta_df$theta_mean[g],
                     sd = theta_df$theta_sd[g])
      theta_e <- rnorm(N,mean = theta_df$ers_mean[g], 
                       sd = theta_df$ers_sd[g])
      
      #Calculates 1 to N, then N+1 to N * 2, etc. Used to indicate which rows
      #of obs_mat to write answers to
      groupmembers <- (1 + (g - 1) * N):(N * g)
      
      #4 category items now, so no middle option
      for(j in 1:Nitems){
        
        #IRTree category probability calculation
        #Step 1: Choose between agree or disagree
        score <- exp(node1pars[1] * theta - node1pars[3])/
          (1+exp(node1pars[1] * theta - node1pars[3]))
        score <- matrix(c(score, 1 - score), nrow = N)
        
        #Step 2a: Choose between extreme or non-extreme disagreeing
        score2 <- exp(node2pars[2] * theta_e + node2pars[1] * theta - node2pars[3])/
          (1+exp(node2pars[2] * theta_e + node2pars[1] * theta - node2pars[3]))
        score2 <- matrix(c(score2, 1-score2), nrow = N)
        
        #Step 2b: Choose between extreme or non-extreme agreeing
        score3 <- exp(node3pars[2] * theta_e + node3pars[1] * theta - node3pars[3])/
          (1+exp(node3pars[2] * theta_e + node3pars[1] * theta - node3pars[3]))
        score3 <- matrix(c(score3, 1-score3), nrow = N)
        
        #Calculate probabilites of participant answering a category for every
        #participant
        cat_prob_tree <- matrix(NA, nrow = N, ncol = categories)
        cat_prob_tree[,1] <- score[,2] * score2[,1]
        cat_prob_tree[,2] <- score[,2] * score2[,2]
        cat_prob_tree[,3] <- score[,1] * score3[,2]
        cat_prob_tree[,4] <- score[,1] * score3[,1]
        
        #Convert probabilities to actual answers
        cat_prob_list <- split(cat_prob_tree, row(cat_prob_tree))
        obs_mat[groupmembers,j] <- sapply(cat_prob_list, function(x)sample(1:4, 1, 
                                                                           replace = TRUE, prob = x),
                                          simplify = TRUE)
      }
      #Add group membership to obs_mat
      obs_mat[groupmembers,(Nitems + 1)] <- g
    }
    
  } else {
    stop("Data generating model is not correctly specified. Please specify one 
    of the following models for gen_model argument:
         1. MNRM-Falk
         2. IRTree-M")
  }
  
  items_g1 <- obs_mat[1:N,1:Nitems]
  items_g2 <- obs_mat[(N+1):(N*2), 1:Nitems]
  proportions <- rep(NA, 8)
  proportions[1] <- sum(items_g1 == 1)/(length(items_g1))
  proportions[2] <- sum(items_g1 == 2)/(length(items_g1))
  proportions[3] <- sum(items_g1 == 3)/(length(items_g1))
  proportions[4] <- sum(items_g1 == 4)/(length(items_g1))
  
  proportions[5] <- sum(items_g2 == 1)/(length(items_g2))
  proportions[6] <- sum(items_g2 == 2)/(length(items_g2))
  proportions[7] <- sum(items_g2 == 3)/(length(items_g2))
  proportions[8] <- sum(items_g2 == 4)/(length(items_g2))
  
  theta_item_cor_g1 <- ERS_extreme_cor_g1 <- theta_item_cor_g2 <- 
    ERS_extreme_cor_g2 <- rep(NA, Nitems)
  for(i in 1: Nitems){
    theta_item_cor_g1[i] <- cor(theta_store[1:N,1], items_g1[,i])
    ERS_extreme_cor_g1[i] <- cor(theta_store[1:N,2], (items_g1[,i] == 1| items_g1[,i] == 4))
    
    theta_item_cor_g2[i] <- cor(theta_store[(N+1):(N*2), 1], items_g2[,i])
    ERS_extreme_cor_g2[i] <- cor(theta_store[(N+1):(N*2), 2], (items_g2[,i] == 1| items_g2[,i] == 4))
  }
  
  #Make sure items have names for mirt
  itemnames <- paste0("item", 1:Nitems)
  colnames(obs_mat) <- c(itemnames, "group")
  
  GPCM_mod <- multipleGroup(obs_mat[,1:Nitems],
                            group = factor(obs_mat[, (Nitems + 1)]),
                            model = paste0("Theta = 1-", Nitems),
                            itemtype = "gpcm",
                            method = estimator,
                            invariance = c("free_mean", "free_var", "slopes",
                                           "intercepts"),
                            technical = list(NCYCLES = 2000))
  
  #Extract relevant outcomes
  outcome <- rbind(GPCM_mod@ParObjects[["pars"]][[1]]@ParObjects[["pars"]][[Nitems + 1]]@par,
                    GPCM_mod@ParObjects[["pars"]][[2]]@ParObjects[["pars"]][[Nitems + 1]]@par)
  colnames(outcome) <- GPCM_mod@ParObjects[["pars"]][[1]]@ParObjects[["pars"]][[Nitems + 1]]@parnames
  
  return(list(category_dist = proportions,
              theta_cors_g1 = theta_item_cor_g1,
              theta_cors_g2 = theta_item_cor_g2,
              ers_cors_g1 = ERS_extreme_cor_g1,
              ers_cors_g2 = ERS_extreme_cor_g2,
              outcome))
  
}

theta_df <- data.frame(theta_mean = c(0, 0), theta_sd = c(1,1),
                       ers_mean = c(0, 1), ers_sd = c(1,1))

grid <- expand.grid(shift = c("0", "0.5", "1", "1.5"),
                    alpha = c("1, 1", "1.5, 1.5", "2, 2", "1, 1.5", "1.5, 1"),
                    treshold = c("-1, 0, 1", "-1.5, 0, 1.5", "-1.25, 0, 1.25"),
                    stringsAsFactors = F)

outcomes <- vector(mode = "list", length = 60)

for(i in 1:nrow(grid)){
  
  tresholds <- as.numeric(unlist(strsplit(grid[i,3], split = ",")))
  alpha <- as.numeric(unlist(strsplit(grid[i,2], split = ",")))
  treshold_shift <- as.numeric(unlist(strsplit(grid[i,1], split = ",")))
  
  intercept <- c(0, tresholds[3], tresholds[3], 0)
  intercept[2] <- intercept[2] - treshold_shift
  intercept[3] <- intercept[3] - 2*treshold_shift
  intercept[4] <- intercept[4] - 3*treshold_shift
  
  intercept <- intercept * alpha[1]
  
  outcomes[[i]] <- generate(gen_model = "MNRM-Falk", N = 50000,
                            categories = 4, groups = 2,
                            theta_df = theta_df, Nitems = 10,
                            intercept = intercept, alpha = alpha)
}

table <- matrix(NA, nrow = 60, ncol = 14)

for(i in 1:60){
  table[i, 1:8] <- outcomes[[i]][["category_dist"]]
  table[i, 9] <- mean(outcomes[[i]][["theta_cors_g1"]])
  table[i, 10] <- mean(outcomes[[i]][["theta_cors_g2"]])
  table[i, 11] <- mean(outcomes[[i]][["ers_cors_g1"]])
  table[i, 12] <- mean(outcomes[[i]][["ers_cors_g2"]])
  table[i, 13] <- outcomes[[i]][[6]][2,1]
  table[i, 14] <- (outcomes[[i]][[6]][2,2] - 1)
}

table <- table[,c(1:9, 11, 13, 14)]
write.csv(table, file = "table.csv")

generate(gen_model = "MNRM-Falk", N = 500000,
         categories = 4, groups = 2,
         theta_df = theta_df, Nitems = 10,
         intercept = c(0, 1.5, 1.5, 0), alpha = c(1.5, 1.5))

colnames(MNRM_data) <- c(paste0("item", 1:10), "group")

GPCM_mod <- multipleGroup(MNRM_data[,1:10],
                          group = factor(MNRM_data[, 11]),
                          model = paste0("Theta = 1-", 10),
                          itemtype = "gpcm",
                          method = "EM",
                          invariance = c("free_mean", "free_var", "slopes",
                                         "intercepts"),
                          technical = list(NCYCLES = 2000))

GPCM_mod@ParObjects[["pars"]][[1]]@ParObjects[["pars"]][[1]]@par

tree_dat <- matrix(NA, nrow = 500000*2, ncol = 10 * (4 - 1))
for(i in 1:10){
  tree_dat[MNRM_data[,i] <= 2, (i-1)*3 +1] <- 0 
  tree_dat[MNRM_data[,i] >= 3, (i-1)*3 +1] <- 1
  
  tree_dat[MNRM_data[,i] == 2, (i-1)*3 + 2] <- 0 
  tree_dat[MNRM_data[,i] == 1, (i-1)*3 + 2] <- 1
  
  tree_dat[MNRM_data[,i] == 3, (i-1)*3 + 3] <- 0 
  tree_dat[MNRM_data[,i] == 4, (i-1)*3 + 3] <- 1 
}

items <- sort(rep(1:10, 3))
items <- paste0("item", items)
nodes <- 1:3
nodes <- rep(paste0("node", nodes), 10)
colnames(tree_dat) <- paste0(items, nodes)

ersloadings <- rep(NA, (10 * 2))
node1alpha <- node2alpha <- node3alpha <- rep(NA, 10)
for(i in 1:10){
  ersloadings[(i-1)*2+1] <- (i-1)*3+2
  ersloadings[(i-1)*2+2] <- (i-1)*3+3
  node1alpha[i] <- (i-1)*3+1
  node2alpha[i] <- (i-1)*3+2
  node3alpha[i] <- (i-1)*3+3
}
#Triple number of items for tree model
treeNitems <- 10 * 3

model <- paste0("Theta = 1-", treeNitems, "
                ", "ERS = ", paste0(ersloadings, collapse = ","), "
                CONSTRAINB = (1-", treeNitems, ", a1), (1-", treeNitems, ", a2), (1-",
                treeNitems, ", d)", "
                ", "CONSTRAIN = ", paste0("(",node2alpha, ",", node3alpha, ", a1)", collapse = ","), ",", "
                     ", paste0("(",node2alpha, ",", node3alpha, ", a2)", collapse = ","), "
                    ")

model2 <- paste0("FREE = (GROUP, COV_21)", "
                     ", "FREE [",2, "] = (GROUP, COV_11)", "
        ", "FREE [",2, "] = (GROUP, COV_22)", "
        ", "FREE [",2, "] = (GROUP, MEAN_1)", "
        ", "FREE [",2, "] = (GROUP, MEAN_2)",
                 collapse = "")

model <- paste0(model, model2, collapse = "")

Tree_mod <- multipleGroup(tree_dat,
                          model = model,
                          method = "EM",
                          itemtype = rep(c("2PL", "Node2TreePL", "Node3TreePL"), 10),
                          group = factor(MNRM_data[, 11]),
                          technical = list(NCYCLES = 2000),
                          customItems=list(Node2TreePL = Node2TreePL,
                                           Node3TreePL = Node3TreePL))

Tree_mod@ParObjects[["pars"]][[1]]@ParObjects[["pars"]][[1]]@par
Tree_mod@ParObjects[["pars"]][[1]]@ParObjects[["pars"]][[2]]@par
Tree_mod@ParObjects[["pars"]][[1]]@ParObjects[["pars"]][[3]]@par
