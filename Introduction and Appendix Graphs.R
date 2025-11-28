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
library(RColorBrewer)

old.par <- par()

#Create customitems
name <- 'Node2TreePL'
par <- c(a1 = 1, a2 = 1, d = 0)
est <- c(TRUE, TRUE, TRUE)
P.Node2TreePL <- function(par,Theta, ncat){
  a1 <- par[1]
  a2 <- par[2]
  d <- par[3]
  P1 <- exp((a1 * Theta[,1] - a2 * Theta[,2] + d)) /
    (1 + exp((a1 * Theta[,1] - a2 * Theta[,2] + d)))
  cbind(1-P1, P1)
}
Node2TreePL <- createItem(name, par=par, est=est, P=P.Node2TreePL, 
                          derivType = 'symbolic', derivType.hss = 'symbolic')

name <- 'Node3TreePL'
par <- c(a1 = 1, a2 = 1, d = 0)
est <- c(TRUE, TRUE, TRUE)
P.Node3TreePL <- function(par,Theta, ncat){
  a1 <- par[1]
  a2 <- par[2]
  d <- par[3]
  P1 <- exp((a1 * Theta[,1] + a2 * Theta[,2] + d)) /
    (1 + exp((a1 * Theta[,1] + a2 * Theta[,2] + d)))
  cbind(1-P1, P1)
}
Node3TreePL <- createItem(name, par=par, est=est, P=P.Node3TreePL,
                          derivType = 'symbolic', derivType.hss = 'symbolic')

rnorm(1)
generate_graph <- function(gen_model = "MNRM-Falk", Nitems = 10, alpha = c(1,1), 
                           thresholds = c(-1, 0, 1), thetavalues = c(0, 0.5, 1, -0.5, -1),
                           N = 100, categories = 4, groups = 1, theta_df,
                           estimator = "EM", xaxis = "ERS",
                           node1pars = matrix(rep(c(1, 0, 0), 10),
                                              nrow = 10, ncol = 3), 
                           node2pars = matrix(rep(c(1, 0, 0), 10),
                                              nrow = 10, ncol = 3), 
                           node3pars = matrix(rep(c(1, 0, 0), 10),
                                              nrow = 10, ncol = 3),
                           est_GPCM = FALSE,
                           extra_lines = TRUE,
                           seperate_tree_plot = FALSE){
  
  #Save seed that was used
  seed <- .Random.seed
  
  #Generate matrices to store observed answers + matrix to
  #store true trait values
  obs_mat <- matrix(rep(NA, N*Nitems), ncol = Nitems)
  theta_store <- matrix(NA, nrow = N, ncol = 2)
  
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
    
    #Draw participant ability from normal distribution with means and standard
    #deviations as specified in theta_df
    theta <- rnorm(N, mean = theta_df$theta_mean, 
                   sd = theta_df$theta_sd)
    
    #Draw participant ERS score from normal distribution with means and sd's
    #as specified in theta_df
    theta_ersmrs <- rnorm(N, mean = theta_df$ers_mean, 
                          sd = theta_df$ers_sd)
    
    #Combine into single matrix
    theta <- cbind(theta, theta_ersmrs)
    
    #Set item difficulty spread. 0 For graph function
    difficulty <- rep(0, length = Nitems)
    
    for(j in 1:Nitems){
      
      #Convert item thresholds into intercepts for mirt
      item_threshold <- thresholds - difficulty[j]
      int_2 <- -item_threshold[1] * alpha[1]
      int_3 <- int_2 - item_threshold[2] * alpha[1]
      int_4 <- int_3 - item_threshold[3] * alpha[1]
      intercept <- c(0, int_2, int_3, int_4)
      
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
      obs_mat[, j] <- obs_ans
    }
    
    #Store theta in theta_store
    theta_store <- theta
    
    #Data generation under IRTree multidimensional model 
    #(B?ckenholt & Meiser, 2017)
  } else if (gen_model == "IRTree-M"){
    
    #Draw substantive and ERS scores from normal distribution with mean and
    #sd specified in theta_df
    theta <- rnorm(N, mean = theta_df$theta_mean,
                   sd = theta_df$theta_sd)
    theta_e <- rnorm(N,mean = theta_df$ers_mean, 
                     sd = theta_df$ers_sd)
    
    #Store theta in theta_store
    theta_store <- cbind(theta, theta_e)
    
    #4 category items now, so no middle option
    for(j in 1:Nitems){
      
      
      #IRTree category probability calculation
      #Step 1: Choose between agree or disagree. 1 = agree, 0 = disagree
      score <- exp(node1pars[j,1] * theta + node1pars[j,3])/
        (1+exp(node1pars[j,1] * theta + node1pars[j,3]))
      score <- matrix(c(score, 1 - score), nrow = N)
      
      #Step 2a: Choose between extreme or non-extreme disagreeing.
      #1 = non-extreme disagree, 0 = extreme disagree
      score2 <- exp(node2pars[j,1] * theta - node2pars[j,2] * theta_e + node2pars[j,3])/
        (1+exp(node2pars[j,1] * theta - node2pars[j,2] * theta_e + node2pars[j,3]))
      score2 <- matrix(c(score2, 1-score2), nrow = N)
      
      #Step 2b: Choose between extreme or non-extreme agreeing
      #1 = non-extreme agree, 0 = extreme agree
      score3 <- exp(node3pars[j,2] * theta_e + node3pars[j,1] * theta + node3pars[j,3])/
        (1+exp(node3pars[j,2] * theta_e + node3pars[j,1] * theta + node3pars[j,3]))
      score3 <- matrix(c(score3, 1-score3), nrow = N)
      
      #Calculate probabilites of participant answering a category for every
      #participant
      cat_prob_tree <- matrix(NA, nrow = N, ncol = categories)
      cat_prob_tree[,1] <- score[,2] * score2[,2]
      cat_prob_tree[,2] <- score[,2] * score2[,1]
      cat_prob_tree[,3] <- score[,1] * score3[,2]
      cat_prob_tree[,4] <- score[,1] * score3[,1]
      
      #Convert probabilities to actual answers
      cat_prob_list <- split(cat_prob_tree, row(cat_prob_tree))
      obs_mat[,j] <- sapply(cat_prob_list, 
                            function(x)sample(1:4, 1, replace = TRUE, prob = x),
                            simplify = TRUE)
    }
    
  } else {
    stop("Data generating model is not correctly specified. Please specify one 
    of the following models for gen_model argument:
         1. MNRM-Falk
         2. IRTree-M")
  }
  
  #Set names of true trait matrix
  colnames(theta_store) <- c("Theta", "ERS")
  
  #MNRM evaluation
  #Generate s matrix
  s_gen <- matrix(0:(categories - 1), nrow = 1)
  
  #Generate ERS scoring matrix
  ers_mat <- matrix(rep(0, categories), nrow = 1)
  ers_mat[,1] <- ers_mat[,length(ers_mat)] <- 1
  
  #Combine ORS and ERS matrix
  s_gen <- rbind(s_gen, ers_mat)
  
  #Generate s_mirt matrix to evaluate model
  s_mirt <- vector(mode = "list", length = Nitems)
  items <- 1:Nitems
  for(j in 1:Nitems){
    s_mirt[[j]] <- t(s_gen)
    
    #Check if all items actually have 4 categories. If not, shorten the s_gen
    #matrix
    if(length(table(obs_mat[,j])) < 4){
      s_mirt[[j]] <- t(s_gen[,1:3])
      s_mirt[[j]][3,2] <- 1
      items <- items[-j]
    }
  }
  
  #Make sure items have names for mirt
  itemnames <- paste0("item", 1:Nitems)
  colnames(obs_mat) <- itemnames
  
  #Generate mirt model syntax. All dimensions load on all items. In addition,
  #all items are set to be completely equal across groups. Covariance between
  #theta and ERS is freely estimated. Covariances and means of latent dimensions
  #in all groups except group 1 are freely estimated.
  model <- paste0("Theta = 1-", Nitems, "
                ", "ERS = 1-", Nitems, "
                ", "FREE = (GROUP, COV_21)", "
                ")
  
  tryCatch({
    #Run mirt evaluation
    MNRM_mod <- mirt(obs_mat[,1:Nitems],
                     model = model,
                     itemtype = "gpcm",
                     method = estimator,
                     gpcm_mats = s_mirt,
                     technical = list(NCYCLES = 2000))
    
    #Gather itempars and take mean for graphs
    itempars_MNRM <- matrix(NA, nrow = Nitems, ncol = 14)
    for(i in 1:Nitems){
      itempars_MNRM[i,] <- MNRM_mod@ParObjects[["pars"]][[i]]@par
    }
    mean_item_MNRM <- colMeans(itempars_MNRM)
    newalpha <- mean_item_MNRM[1:2] * s_gen
    intercept <- mean_item_MNRM[11:14]
    
    #Get trait estimates from MNRM
    theta_est <<- fscores(MNRM_mod)
    
    
    
    
  }, 
  error = function(e){
    MNRM_small <<- NA
    MNRM_est_traits <<- NA
  }
  )
  
  #Convert observed answers to IRTree format
  tree_dat <- matrix(NA, nrow = N, ncol = Nitems * (categories - 1))
  start_val_d1 <- start_val_d2 <- rep(NA, Nitems)
  for(i in 1:Nitems){
    tree_dat[obs_mat[,i] <= 2, (i-1)*3 +1] <- 0 
    tree_dat[obs_mat[,i] >= 3, (i-1)*3 +1] <- 1
    
    tree_dat[obs_mat[,i] == 2, (i-1)*3 + 2] <- 1 
    tree_dat[obs_mat[,i] == 1, (i-1)*3 + 2] <- 0
    
    tree_dat[obs_mat[,i] == 3, (i-1)*3 + 3] <- 0 
    tree_dat[obs_mat[,i] == 4, (i-1)*3 + 3] <- 1 
  }
  
  #Determine start value of intercept based on observed proportion correct
  start_val_d1 <- log(table(tree_dat[,(ncol(tree_dat)/2 + 2)])[1]/
                        table(tree_dat[, (ncol(tree_dat)/2 + 2)])[2])
  start_val_d2 <- -log(table(tree_dat[,(ncol(tree_dat)/2 + 3)])[1]/
                         table(tree_dat[, (ncol(tree_dat)/2 + 3)])[2])
  
  #Evaluate IRTree model using mirt
  #Assign itemnodes names for MIRT
  items <- sort(rep(1:Nitems, 3))
  items <- paste0("item", items)
  nodes <- 1:3
  nodes <- rep(paste0("node", nodes), Nitems)
  colnames(tree_dat) <- paste0(items, nodes)
  
  #Determine which item nodes ERS loads on (second and third nodes). In addition,
  #determine which node alphas are constrained to -1 or 1 respectively.
  ersloadings <- rep(NA, (Nitems * 2))
  node1alpha <- node2alpha <- node3alpha <- rep(NA, Nitems)
  for(i in 1:Nitems){
    ersloadings[(i-1)*2+1] <- (i-1)*3+2
    ersloadings[(i-1)*2+2] <- (i-1)*3+3
    node1alpha[i] <- (i-1)*3+1
    node2alpha[i] <- (i-1)*3+2
    node3alpha[i] <- (i-1)*3+3
  }
  
  tree_dat_adj <- tree_dat
  
  itemtypes <- rep(c("2PL", "Node2TreePL", "Node3TreePL"), Nitems)
  
  for(i in 1:ncol(tree_dat)){
    if(length(table(tree_dat[,i])) < 2){
      tree_dat_adj <- tree_dat_adj[,-i]
      
      ers_shift <- which(ersloadings == i)
      ersloadings <- ersloadings[-ers_shift]
      ersloadings[ers_shift:length(ersloadings)] <- 
        ersloadings[ers_shift:length(ersloadings)] - 1
      
      itemtypes <- itemtypes[-i]
      
      constraints_shift <- ceiling(i/3)
      node2alpha <- node2alpha[-constraints_shift]
      node3alpha <- node3alpha[-constraints_shift]
      node2alpha[constraints_shift : length(node2alpha)] <-
        node2alpha[constraints_shift : length(node2alpha)] - 1
      node3alpha[constraints_shift : length(node3alpha)] <- 
        node3alpha[constraints_shift : length(node3alpha)] - 1
      
    }
  }
  #Triple number of items for tree model
  treeNitems <- ncol(tree_dat_adj)
  
  #Generate mirt syntax. Theta loads on all item nodes, ERS loads on item nodes
  #2 and 3. All items are constrained equal across groups. a1 parameters
  #are set to -1 and 1 for node 2 and node 3 respectively. Covariance between
  #ERS and theta is freely estimated. Latent means and covariances are freely
  #estimated in every group except group 1.
  
  #Constraints
  model <- paste0("Theta = 1-", treeNitems, "
                ", "ERS = ", paste0(ersloadings, collapse = ","), "
                ", "CONSTRAIN = ", paste0("(",node2alpha, ",", node3alpha, ", a1)", collapse = ","), ",", "
                     ", paste0("(",node2alpha, ",", node3alpha, ", a2)", collapse = ","), "
                     ", "START = (", paste0(node2alpha, collapse = ","), ", d, ", start_val_d1, "), (",
                  paste0(node3alpha, collapse = ","), ", d, ", start_val_d2, ") ", "
                    ", " FREE = (GROUP, COV_21)")
  
  
  tryCatch({
    #Evaluate model using MIRT 
    Tree_mod <- mirt(tree_dat_adj,
                     model = model,
                     method = estimator,
                     itemtype = itemtypes,
                     technical = list(NCYCLES = 2000),
                     customItems=list(Node2TreePL = Node2TreePL,
                                      Node3TreePL = Node3TreePL))
    
    itempars_Tree <- matrix(NA, nrow = treeNitems, ncol = 3)
    for(i in 1:treeNitems){
      itempars_Tree[i,] <- Tree_mod@ParObjects[["pars"]][[i]]@par[1:3]
    }
    
    node1pars <- colMeans(itempars_Tree[c(1,4,7,10,13,16,19,22,25,28),])
    node2pars <- colMeans(itempars_Tree[c(2,5,8,11,14,17,20,23,26,29),])
    node3pars <- colMeans(itempars_Tree[c(3,6,9,12,15,18,21,24,27,30),])
    
    theta_est_tree <- fscores(Tree_mod)
  }, 
  error = function(e){
    tryCatch({
      Tree_mod <- mirt(tree_dat_adj,
                       model = model,
                       method = "QMCEM",
                       itemtype = itemtypes,
                       technical = list(NCYCLES = 2000),
                       customItems=list(Node2TreePL = Node2TreePL,
                                        Node3TreePL = Node3TreePL))
      
      itempars_Tree <- matrix(NA, nrow = treeNitems, ncol = 3)
      for(i in 1:treeNitems){
        itempars_Tree[i,] <- Tree_mod@ParObjects[["pars"]][[i]]@par[1:3]
      }
      
      node1pars <- colMeans(itempars_Tree[c(1,4,7,10,13,16,19,22,25,28),])
      node2pars <- colMeans(itempars_Tree[c(2,5,8,11,14,17,20,23,26,29),])
      node3pars <- colMeans(itempars_Tree[c(3,6,9,12,15,18,21,24,27,30),])
      
      theta_est_tree <<- fscores(Tree_mod)
    }, error = function(e){
      Tree_small <<- NA
      Tree_est_traits <<- NA
    }
    )
  }
  )
  
  tryCatch({
    #GPCM estimation. Theta loads on all items.
    GPCM_mod <- mirt(obs_mat[,1:Nitems],
                     model = paste0("Theta = 1-", Nitems),
                     itemtype = "gpcm",
                     method = estimator,
                     technical = list(NCYCLES = 2000))
    
    itempars_GPCM <- matrix(NA, nrow = Nitems, ncol = 9)
    for(i in 1:Nitems){
      itempars_GPCM[i,] <- GPCM_mod@ParObjects[["pars"]][[i]]@par
    }
    mean_item_GPCM <- colMeans(itempars_GPCM)
    
    gpcm_a1 <- mean_item_GPCM[1]
    gpcm_a1 <- gpcm_a1 * c(0,1,2,3)
    gpcm_ds <- mean_item_GPCM[6:9]
    
    GPCM_est_traits <- fscores(GPCM_mod)
  }, 
  error = function(e){
    GPCM_small <<- NA
    GPCM_est_traits <<- NA}
  )
  
  if(xaxis == "ERS"){
    constant <- c(1, 3, 5, 7, 9, 11, 13, 15, 17)
    changing <- c(2, 4, 6, 8, 10, 12, 14, 16, 18)
    label <- "ERS score"
    label2 <- "Theta"
  } else {
    constant <- c(2, 4, 6, 8, 10, 12, 14, 16, 18)
    changing <- c(1, 3, 5, 7, 9, 11, 13, 15, 17)
    label <- "Theta score"
    label2 <- "ERS"
  }
  
  theta <- matrix(NA, nrow = 3000, ncol = 18)
  theta[,constant] <- c(rep(thetavalues[1], 3000),
                        rep(thetavalues[2], 3000),
                        rep(thetavalues[3], 3000), 
                        rep(thetavalues[4], 3000),
                        rep(thetavalues[5], 3000), 
                        rep(thetavalues[6], 3000),
                        rep(thetavalues[7], 3000), 
                        rep(thetavalues[8], 3000),
                        rep(thetavalues[9], 3000))
  theta[,changing] <- seq(from = -3, to = 3, length = 3000)
  
  #Initialize category probability matrices
  cat_probs_MNRM <- cat_probs_tree <- cat_probs_gpcm <-
    matrix(NA, nrow = 3000, ncol = 36)
  obs_mat_MNRM <- obs_mat_tree <- matrix(NA, nrow = 3000, ncol = 9)
  
  for(i in 1:9){
    #Specify which columns of theta matrix are substantive and which are ERS
    thetacols <- ((i-1)*2 +1)
    erscols <- ((i-1)*2+2)
    
    #MNRM category probability calculation
    numerator <- exp(apply(theta[,thetacols:erscols], 1, "%*%", 
                           newalpha) + intercept)
    denominator <- colSums(numerator)
    cat_probs_MNRM[,((i-1)*4 +1):((i-1)*4 +4)] <- cat_prob_MNRM <- t(numerator)/denominator
    
    #GPCM category probability calucaltion
    if(est_GPCM == TRUE){
      numerator <- exp(apply(matrix(theta[,thetacols], ncol = 1), 1, "%*%", 
                             matrix(gpcm_a1, nrow = 1)) + gpcm_ds)
      denominator <- colSums(numerator)
      cat_probs_gpcm[,((i-1)*4 +1):((i-1)*4 +4)] <- t(numerator)/denominator
    }
    
    #IRTree category probability calculation
    score <- exp(node1pars[1] * theta[,thetacols] + node1pars[3])/
      (1 + exp(node1pars[1] * theta[,thetacols] + node1pars[3]))
    score <- matrix(c(score, 1 - score), nrow = 3000)
    
    score2 <- exp(node2pars[1] * theta[,thetacols] - node2pars[2] * theta[,erscols] + node2pars[3])/
      (1 + exp(node2pars[1] * theta[,thetacols] - node2pars[2] * theta[,erscols] + node2pars[3]))
    score2 <- matrix(c(score2, 1-score2), nrow = 3000)
    
    score3 <- exp(node3pars[1] * theta[,thetacols] + node3pars[2] * theta[,erscols] + node3pars[3])/
      (1+ exp(node3pars[1] * theta[,thetacols] + node3pars[2] * theta[,erscols] + node3pars[3]))
    score3 <- matrix(c(score3, 1-score3), nrow = 3000)
    
    cat_prob_tree <- matrix(NA, nrow = 3000, ncol = 4)
    cat_prob_tree[,1] <- score[,2] * score2[,2]
    cat_prob_tree[,2] <- score[,2] * score2[,1]
    cat_prob_tree[,3] <- score[,1] * score3[,2]
    cat_prob_tree[,4] <- score[,1] * score3[,1]
    
    cat_probs_tree[,((i-1)*4 +1):((i-1)*4 +4)] <- cat_prob_tree
    
    #Make white canvas to draw lines on
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    plot(x = theta[,changing[1]], y = rep(0, 3000), col = "white", ylim = c(0,1),
         ylab = "Probability",
         xlab = label)
    #MNRM lines
    lines(x = theta[,changing[1]], y = cat_probs_MNRM[,(i-1) * 4 + 1], 
          col = "red", lwd = 3)
    lines(x = theta[,changing[1]], y = cat_probs_MNRM[,(i-1) * 4 + 2], 
          col = "orange", lwd = 3)
    lines(x = theta[,changing[1]], y = cat_probs_MNRM[,(i-1) * 4 + 3], 
          col = "black", lwd = 3)
    lines(x = theta[,changing[1]], y = cat_probs_MNRM[,(i-1) * 4 + 4], 
          col = "green", lwd = 3)
    
    legend("topright", 
           legend = c("Category 1", "Category 2", "Category 3", "Category 4"), 
           col = c("red", "orange", "black", "green"), 
           lty = "solid", lwd = 3, inset = c(-0.20, 0))
    
    
    if(seperate_tree_plot == TRUE){
      plot(x = theta[,changing[1]], y = rep(0, 3000), col = "white", ylim = c(0,1),
           ylab = "Probability",
           xlab = label)
    }
    #IRTree lines
    lines(x = theta[,changing[1]], y = cat_probs_tree[,(i-1) * 4 + 1], 
          col = "red", lty = "dotted", lwd = 3)
    lines(x = theta[,changing[1]], y = cat_probs_tree[,(i-1) * 4 + 2], 
          col = "orange", lty = "dotted", lwd = 3)
    lines(x = theta[,changing[1]], y = cat_probs_tree[,(i-1) * 4 + 3],
          col = "black", lty = "dotted", lwd = 3)
    lines(x = theta[,changing[1]], y = cat_probs_tree[,(i-1) * 4 + 4], 
          col = "green", lty = "dotted", lwd = 3)
    
    legend("topright", 
           legend = c("Category 1", "Category 2", "Category 3", "Category 4"), 
           col = c("red", "orange", "black", "green"), 
           lty = "dotted", lwd = 3, inset = c(-0.20, 0))
    
    par(old.par)
    
  }
  
  #Calculate agreement probabilites under both MNRM and IRTree
  endorse_prob <- rowSums(cat_probs_MNRM[,c(3,4)])
  endorse_prob_05 <- rowSums(cat_probs_MNRM[,c(7,8)])
  endorse_prob_1 <- rowSums(cat_probs_MNRM[,c(11,12)])
  endorse_prob_15 <- rowSums(cat_probs_MNRM[,c(15,16)])
  endorse_prob_2 <- rowSums(cat_probs_MNRM[,c(19,20)])
  endorse_prob_min05 <- rowSums(cat_probs_MNRM[,c(23,24)])
  endorse_prob_min1 <- rowSums(cat_probs_MNRM[,c(27,28)])
  endorse_prob_min15 <- rowSums(cat_probs_MNRM[,c(31,32)])
  endorse_prob_min2 <- rowSums(cat_probs_MNRM[,c(35,36)])
  
  endorse_prob_tree <- rowSums(cat_probs_tree[,c(3,4)])
  endorse_prob_tree_05 <- rowSums(cat_probs_tree[,c(7,8)])
  endorse_prob_tree_1 <- rowSums(cat_probs_tree[,c(11,12)])
  endorse_prob_tree_15 <- rowSums(cat_probs_tree[,c(15,16)])
  endorse_prob_tree_2 <- rowSums(cat_probs_tree[,c(19,20)])
  endorse_prob_tree_min05 <- rowSums(cat_probs_tree[,c(23,24)])
  endorse_prob_tree_min1 <- rowSums(cat_probs_tree[,c(27,28)])
  endorse_prob_tree_min15 <- rowSums(cat_probs_tree[,c(31,32)])
  endorse_prob_tree_min2 <- rowSums(cat_probs_tree[,c(35,36)])
  
  #Calculate conditional probability of answering 4 if extreme response for 
  #MNRM and IRTree
  cond_cat_prob_14_MNRM <- cat_probs_MNRM[,4]/(cat_probs_MNRM[,1] + 
                                                 cat_probs_MNRM[,4])
  cond_cat_prob_14_MNRM_05 <- cat_probs_MNRM[,8]/(cat_probs_MNRM[,5] + 
                                                    cat_probs_MNRM[,8])
  cond_cat_prob_14_MNRM_1 <- cat_probs_MNRM[,12]/(cat_probs_MNRM[,9] + 
                                                    cat_probs_MNRM[,12])
  cond_cat_prob_14_MNRM_15 <- cat_probs_MNRM[,16]/(cat_probs_MNRM[,13] + 
                                                     cat_probs_MNRM[,16])
  cond_cat_prob_14_MNRM_2 <- cat_probs_MNRM[,20]/(cat_probs_MNRM[,17] + 
                                                    cat_probs_MNRM[,20])
  cond_cat_prob_14_MNRM_min05 <- cat_probs_MNRM[,24]/(cat_probs_MNRM[,21] + 
                                                        cat_probs_MNRM[,24])
  cond_cat_prob_14_MNRM_min1 <- cat_probs_MNRM[,28]/(cat_probs_MNRM[,25] + 
                                                       cat_probs_MNRM[,28])
  cond_cat_prob_14_MNRM_min15 <- cat_probs_MNRM[,32]/(cat_probs_MNRM[,29] + 
                                                        cat_probs_MNRM[,32])
  cond_cat_prob_14_MNRM_min2 <- cat_probs_MNRM[,36]/(cat_probs_MNRM[,33] + 
                                                       cat_probs_MNRM[,36])
  
  cond_cat_prob_14_tree <- cat_probs_tree[,4]/(cat_probs_tree[,1] +
                                                 cat_probs_tree[,4])
  cond_cat_prob_14_tree_05 <- cat_probs_tree[,8]/(cat_probs_tree[,5] +
                                                    cat_probs_tree[,8])
  cond_cat_prob_14_tree_1 <- cat_probs_tree[,12]/(cat_probs_tree[,9] +
                                                    cat_probs_tree[,12])
  cond_cat_prob_14_tree_15 <- cat_probs_tree[,16]/(cat_probs_tree[,13] +
                                                     cat_probs_tree[,16])
  cond_cat_prob_14_tree_2 <- cat_probs_tree[,20]/(cat_probs_tree[,17] +
                                                    cat_probs_tree[,20])
  cond_cat_prob_14_tree_min05 <- cat_probs_tree[,24]/(cat_probs_tree[,21] +
                                                        cat_probs_tree[,24])
  cond_cat_prob_14_tree_min1 <- cat_probs_tree[,28]/(cat_probs_tree[,25] +
                                                       cat_probs_tree[,28])
  cond_cat_prob_14_tree_min15 <- cat_probs_tree[,32]/(cat_probs_tree[,29] +
                                                        cat_probs_tree[,32])
  cond_cat_prob_14_tree_min2 <- cat_probs_tree[,36]/(cat_probs_tree[,33] +
                                                       cat_probs_tree[,36])
  
  #Draw plots of agreement probability
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  
  plot(x = theta[,changing[1]], y = endorse_prob, type = "l", 
       xlab = 
         label,
       ylab = "Probability of agreement", ylim = c(0,1), lwd = 3)
  #MNRM lines
  lines(x = theta[,changing[1]], y = endorse_prob_05, col = "cyan", lwd = 3)
  lines(x = theta[,changing[1]], y = endorse_prob_1, col = "blue", lwd = 3)
  
  if(extra_lines == TRUE){
    lines(x = theta[,changing[1]], y = endorse_prob_15, col = "wheat", lwd = 3)
    lines(x = theta[,changing[1]], y = endorse_prob_2, col = "lightgreen", lwd = 3)
  }
  
  lines(x = theta[,changing[1]], y = endorse_prob_min05, col = "purple", lwd = 3)
  lines(x = theta[,changing[1]], y = endorse_prob_min1, col = "darkorange", lwd = 3)
  
  if(extra_lines == TRUE){
    lines(x = theta[,changing[1]], y = endorse_prob_min15, col = "red", lwd = 3)
    lines(x = theta[,changing[1]], y = endorse_prob_min2, col = "darkred", lwd = 3)
  }
  
  #IRTree lines
  lines(x = theta[,changing[1]], y = endorse_prob_tree, col = "black", 
        lwd = 3, lty = "dotted")
  lines(x = theta[,changing[1]], y = endorse_prob_tree_05, col = "cyan", 
        lwd = 3, lty = "dotted")
  lines(x = theta[,changing[1]], y = endorse_prob_tree_1, col = "blue", 
        lwd = 3, lty = "dotted")
  
  if(extra_lines == TRUE){
    lines(x = theta[,changing[1]], y = endorse_prob_tree_15, col = "wheat", 
          lwd = 3, lty = "dotted")
    lines(x = theta[,changing[1]], y = endorse_prob_tree_2, col = "lightgreen", 
          lwd = 3, lty = "dotted")
  }
  
  lines(x = theta[,changing[1]], y = endorse_prob_tree_min05, col = "purple", 
        lwd = 3, lty = "dotted")
  lines(x = theta[,changing[1]], y = endorse_prob_tree_min1, col = "darkorange", 
        lwd = 3, lty = "dotted")
  
  if(extra_lines == TRUE){
    lines(x = theta[,changing[1]], y = endorse_prob_tree_min15, col = "red", 
          lwd = 3, lty = "dotted")
    lines(x = theta[,changing[1]], y = endorse_prob_tree_min2, col = "darkred", 
          lwd = 3, lty = "dotted")
  }
  
  legend("topright",
         legend = c(expression(paste(theta[1], "= 1")),
                    expression(paste(theta[1], "= 0.5")),
                    expression(paste(theta[1], "= 0")),
                    expression(paste(theta[1], "= - 0.5")),
                    expression(paste(theta[1], "= - 1"))
         ),
         col = c("blue", "cyan", "black", "purple", "darkorange"),
         lwd = 5,
         inset = c(-0.18, 0))
  
  par(old.par)
  
  #Probability of 4 given extreme response
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  
  plot(x = theta[,changing[1]], y = cond_cat_prob_14_MNRM, type = "l", 
       ylab = "Probability",
       xlab =  label, ylim = c(0, 1), lwd = 3)
  lines(x = theta[,changing[1]], y = cond_cat_prob_14_MNRM_05, 
        col = "cyan", lwd = 3)
  lines(x = theta[,changing[1]], y = cond_cat_prob_14_MNRM_1,
        col = "blue", lwd = 3)
  
  if(extra_lines == TRUE){
    lines(x = theta[,changing[1]], y = cond_cat_prob_14_MNRM_15,
          col = "wheat", lwd = 3)
    lines(x = theta[,changing[1]], y = cond_cat_prob_14_MNRM_2,
          col = "lightgreen", lwd = 3)
  }
  lines(x = theta[,changing[1]], y = cond_cat_prob_14_MNRM_min05,
        col = "purple", lwd = 3)
  lines(x = theta[,changing[1]], y = cond_cat_prob_14_MNRM_min1, 
        col = "darkorange", lwd = 3)
  
  if(extra_lines == TRUE){
    lines(x = theta[,changing[1]], y = cond_cat_prob_14_MNRM_min15,
          col = "red", lwd = 3)
    lines(x = theta[,changing[1]], y = cond_cat_prob_14_MNRM_min2, 
          col = "darkred", lwd = 3)
  }
  
  
  lines(x = theta[,changing[1]], y = cond_cat_prob_14_tree, col = "black",
        lwd = 3, lty = "dotted")
  lines(x = theta[,changing[1]], y = cond_cat_prob_14_tree_05, col = "cyan",
        lwd = 3, lty = "dotted")
  lines(x = theta[,changing[1]], y = cond_cat_prob_14_tree_1, col = "blue",
        lwd = 3, lty = "dotted")
  
  if(extra_lines == TRUE){
    lines(x = theta[,changing[1]], y = cond_cat_prob_14_tree_15, col = "wheat",
          lwd = 3, lty = "dotted")
    lines(x = theta[,changing[1]], y = cond_cat_prob_14_tree_2, col = "lightgreen",
          lwd = 3, lty = "dotted")
  }
  
  lines(x = theta[,changing[1]], y = cond_cat_prob_14_tree_min05, col = "purple",
        lwd = 3, lty = "dotted")
  lines(x = theta[,changing[1]], y = cond_cat_prob_14_tree_min1, col = "darkorange",
        lwd = 3, lty = "dotted")
  
  if(extra_lines == TRUE){
    lines(x = theta[,changing[1]], y = cond_cat_prob_14_tree_min15, col = "red",
          lwd = 3, lty = "dotted")
    lines(x = theta[,changing[1]], y = cond_cat_prob_14_tree_min2, col = "darkred",
          lwd = 3, lty = "dotted")
  }
  
  legend("topright",
         legend = c(expression(paste(theta[1], "= 1")),
                    expression(paste(theta[1], "= 0.5")),
                    expression(paste(theta[1], "= 0")),
                    expression(paste(theta[1], "= - 0.5")),
                    expression(paste(theta[1], "= - 1"))
         ),
         col = c("blue", "cyan", "black", "purple", "darkorange"),
         lwd = 5,
         inset = c(-0.18, 0))
  
  par(old.par)
}

#debug(generate_graph)

theta_df <- data.frame(theta_mean = 0, theta_sd = 1,
                       ers_mean = 0, ers_sd = 1)

set.seed(112)
#Order of plots: First, category probability plots. These are made in the
#order of thetavalues (which refer to ability when xaxis = ERS (default), and to ERS when
#xaxis = theta). So first plot is theta = 0, second is theta = 0.5, etc.
#After these 1 plot of agreement probability and 1 plot of extreme response 
#probability
generate_graph(gen_model = "MNRM-Falk", N = 50000,
               categories = 4, groups = 2,
               theta_df = theta_df, Nitems = 10,
               alpha = c(1.5, 1.5), thresholds = c(-1, 0, 1),
               thetavalues = c(0, 0.5, 1, 1.5, 2, -0.5, -1, -1.5, -2))

#Same plot but now with theta on the x-axis instead of ERS. Thetavalues now 
#refers to ERS values
generate_graph(gen_model = "MNRM-Falk", N = 50000,
               categories = 4, groups = 2,
               theta_df = theta_df, Nitems = 10,
               alpha = c(1.5, 1.5), thresholds = c(-1, 0, 1),
               thetavalues = c(0, 0.5, 1, 1.5, 2, -0.5, -1, -1.5, -2),
               xaxis = "Theta",
               extra_lines = FALSE,
               seperate_tree_plot = TRUE)

generate_graph(gen_model = "MNRM-Falk", N = 50000,
               categories = 4, groups = 2,
               theta_df = theta_df, Nitems = 10,
               alpha = c(1.5, 1.5), thresholds = c(-1, 0, 1),
               thetavalues = c(0, 0.5, 1, 1.5, 2, -0.5, -1, -1.5, -2),
               extra_lines = FALSE)

set.seed(113)
generate_graph(gen_model = "MNRM-Falk", N = 50000,
               categories = 4, groups = 2,
               theta_df = theta_df, Nitems = 10,
               alpha = c(1.5, 1.5), thresholds = c(0, 1, 2),
               thetavalues = c(0, 0.5, 1, 1.5, 2, -0.5, -1, -1.5, -2))


#IRTree graphs
#Note: Load the IRTree parameters (Treenodes.RData) before running this part
set.seed(114)
generate_graph(gen_model = "IRTree-M", N = 50000,
               categories = 4, groups = 2,
               theta_df = theta_df, Nitems = 10,
               thetavalues = c(0, 0.5, 1, 1.5, 2, -0.5, -1, -1.5, -2),
               node1pars = tree_node1,
               node2pars = tree_node2, 
               node3pars = tree_node3,
               extra_lines = FALSE)

set.seed(115)
generate_graph(gen_model = "IRTree-M", N = 50000,
               categories = 4, groups = 2,
               theta_df = theta_df, Nitems = 10,
               thetavalues = c(0, 0.5, 1, 1.5, 2, -0.5, -1, -1.5, -2),
               node1pars = tree_node1_dif,
               node2pars = tree_node2_dif, 
               node3pars = tree_node3_dif,
               extra_lines = FALSE)


#Combine all data into a single list
#Note: Load the simulation data (Final Data.RData) into R before running this
GPCM_dat_3 <- c(test_GPCM_average_20_1,
                test_GPCM_average_20_2,
                test_GPCM_average_20_3,
                test_GPCM_average_20_4,
                test_GPCM_average_20_5)

GPCM_dat_4 <- c(test_GPCM_difficult_20_1,
                test_GPCM_difficult_20_2,
                test_GPCM_difficult_20_3,
                test_GPCM_difficult_20_4,
                test_GPCM_difficult_20_5)

GPCM_dat <- list(GPCM_dat_3,
                 GPCM_dat_4)


MNRM_dat <- list(test_MNRM_average_20_ERS_min1,
                 test_MNRM_difficult_20_ERS_min1,
                 test_MNRM_average_20,
                 test_MNRM_difficult_20,
                 test_MNRM_average_20_ERS_1,
                 test_MNRM_difficult_20_ERS_1)

Tree_dat <- list(test_Tree_average_20_ERS_min1,
                 test_Tree_difficult_20_ERS_min1,
                 test_Tree_average_20,
                 test_Tree_difficult_20,
                 test_Tree_average_20_ERS_1,
                 test_Tree_difficult_20_ERS_1)

all_data <- list(GPCM_dat,
                 MNRM_dat,
                 Tree_dat)

heatmaps <- function(true_values = TRUE, generating_model_values = FALSE,
                     plot_heatmaps = TRUE, condition = "GPCM",
                     plot_unidimensional_theta = TRUE,
                     plot_unidimensional_ERS = TRUE,
                     y_axis = c(-1.5, 1.5),
                     y_axis_names = c("GPCM substantive trait bias", 
                                      "IRTree substantive trait bias",
                                      "MNRM substantive trait bias"),
                     plot_GPCM_only = FALSE,
                     split_plots = FALSE,
                     plot_unidimensional_theta_ERS = FALSE){
  
  #Initialize matrices for true and estimated traits
  true_traits <- Tree_est <- MNRM_est <- 
    matrix(NA, nrow = 250000, ncol = 2)
  GPCM_est <- matrix(NA, nrow = 250000, ncol = 1)
  
  #initialize matrices to store average bias per bin in the unidimensional plots
  MNRM_theta_bin_1 <- Tree_theta_bin_1 <- GPCM_theta_bin_1 <- 
    MNRM_theta_bin_2 <- Tree_theta_bin_2 <- GPCM_theta_bin_2 <- 
    matrix(NA, nrow = 1, ncol = 20)
  
  GPCM_colours_1 <- GPCM_colours_2 <- 
    MNRM_colours_1 <- MNRM_colours_2 <- 
    Tree_colours_1 <- Tree_colours_2 <-
    rep(NA, 20)
  
  if(condition == "GPCM"){
    k <- 1
  } else if(condition == "MNRM"){
    k <- 2
  } else if(condition == "IRTree"){
    k <- 3
  }
  
  #Extract true, MNRM, IRTree and GPCM theta and ERS values from data
  for(j in 1:length(all_data[[k]])){
    for(i in 1:500){
      true_traits[((i-1) * 500 + 1): (i * 500),] <- 
        all_data[[k]][[j]][[i]][["True_traits"]][501:1000,]
      MNRM_est[((i-1) * 500 + 1): (i * 500),] <- 
        all_data[[k]][[j]][[i]][["MNRM_traits"]][501:1000,]
      Tree_est[((i-1) * 500 + 1): (i * 500),] <- 
        all_data[[k]][[j]][[i]][["Tree_traits"]][501:1000,]
      GPCM_est[((i-1) * 500 + 1): (i * 500),] <- 
        all_data[[k]][[j]][[i]][["GPCM_traits"]][501:1000,]
    }
    
    
    if(true_values == TRUE){
      #Change theta and ERS values into difference between true and estimated
      #values if true values = TRUE. Otherwise, change theta and ERS values into
      #difference between generating model estimates and other models estimated
      #values depending on which model is generating
      MNRM_est <- MNRM_est - true_traits
      Tree_est<- Tree_est - true_traits
      GPCM_est <- GPCM_est - true_traits[,1]
    } else if(generating_model_values == TRUE){
      if(k == 1){
        MNRM_est[,1] <- MNRM_est[,1] - GPCM_est
        Tree_est[,1] <- Tree_est[,1] - GPCM_est
        GPCM_est <- GPCM_est - GPCM_est
      } else if(k == 2){
        Tree_est <- Tree_est - MNRM_est
        GPCM_est <- GPCM_est - MNRM_est[,1]
        MNRM_est <- MNRM_est - MNRM_est
      } else if(k == 3){
        MNRM_est <- MNRM_est - Tree_est
        GPCM_est <- GPCM_est - Tree_est[,1]
        Tree_est <- Tree_est - Tree_est
      }
    }
    
    #Create dataframe of true values and bins true values into bins of width 0.2 from -2 to 2
    true_bins <- data.frame(True_theta = true_traits[,1], True_ERS = true_traits[,2],
                            theta_bin = cut(true_traits[,1], seq(from = -2, to = 2, by = 0.2)),
                            ERS_bin = cut(true_traits[,2], seq(from = -2, to = 2, by = 0.2)))
    
    #Append bins based on true values to difference between true and estimated values
    Tree_df <- data.frame(Tree_theta_bias = Tree_est[,1], Tree_ERS_bias = Tree_est[,2],
                          theta_bin = true_bins[,3], ERS_bin = true_bins[,4])
    GPCM_df <- data.frame(GPCM_theta_bias = GPCM_est[,1], theta_bin = true_bins[,3],
                          ERS_bin = true_bins[,4])
    MNRM_df <- data.frame(MNRM_theta_bias = MNRM_est[,1], MNRM_ERS_bias = MNRM_est[,2],
                          theta_bin = true_bins[,3], ERS_bin = true_bins[,4])
    
    if(plot_unidimensional_theta == TRUE | 
       plot_unidimensional_theta_ERS == TRUE){
      if(plot_unidimensional_theta == TRUE){
        matrix_col <- 1
      } else if (plot_unidimensional_theta_ERS == TRUE){
        matrix_col <- 2
      }
      #Reset graphical parameters
      par(old.par)
      
      #Make 1-dimensional graphs based on just the theta and just the ERS bin
      #Remove all theta values above 2 or below -2
      MNRM_df_plot <- MNRM_df[is.na(MNRM_df[,3]) == FALSE,]
      Tree_df_plot <- Tree_df[is.na(MNRM_df[,3]) == FALSE,]
      GPCM_df_plot <- GPCM_df[is.na(MNRM_df[,3]) == FALSE,]
      bins <- names(table(true_bins[,3]))
      
      #Only do this in conditions without threshold shift
      if(j == 1 | j == 3 | j == 5){
        #Calculate average bias per theta bin and set the colours of the bars
        #to blue if they are negative or red if they are positive
        for(i in 1:20){
          MNRM_theta_bin_1[,i] <- 
            mean(MNRM_df_plot[MNRM_df_plot[,3] == bins[i], matrix_col])
          Tree_theta_bin_1[,i] <- 
            mean(Tree_df_plot[Tree_df_plot[,3] == bins[i], matrix_col])
          if(matrix_col != 2){
            GPCM_theta_bin_1[,i] <- 
              mean(GPCM_df_plot[GPCM_df_plot[,2] == bins[i], matrix_col])
            
            if(GPCM_theta_bin_1[,i] < 0){
              GPCM_colours_1[i] <- "blue"
            } else{
              GPCM_colours_1[i] <- "red"
            }
          }
          if(MNRM_theta_bin_1[,i] < 0){
            MNRM_colours_1[i] <- "blue"
          } else{
            MNRM_colours_1[i] <- "red"
          }
          
          if(Tree_theta_bin_1[,i] < 0){
            Tree_colours_1[i] <- "blue"
          } else{
            Tree_colours_1[i] <- "red"
          }
        }
        #Conditions with threshold shift
      } else if (j == 2| j == 4| j == 6){
        for(i in 1:20){
          MNRM_theta_bin_2[,i] <- 
            mean(MNRM_df_plot[MNRM_df_plot[,3] == bins[i], matrix_col])
          Tree_theta_bin_2[,i] <- 
            mean(Tree_df_plot[Tree_df_plot[,3] == bins[i], matrix_col])
          
          if(matrix_col != 2){
            GPCM_theta_bin_2[,i] <- 
              mean(GPCM_df_plot[GPCM_df_plot[,2] == bins[i], matrix_col])
            
            if(GPCM_theta_bin_2[,i] < 0){
              GPCM_colours_2[i] <- "blue"
            } else{
              GPCM_colours_2[i] <- "red"
            }
          }
          
          if(MNRM_theta_bin_2[,i] < 0){
            MNRM_colours_2[i] <- "blue"
          } else{
            MNRM_colours_2[i] <- "red"
          }
          
          if(Tree_theta_bin_2[,i] < 0){
            Tree_colours_2[i] <- "blue"
          } else{
            Tree_colours_2[i] <- "red"
          }
        }
        
        if(plot_GPCM_only == TRUE){
          split.screen(c(1,2))
          screen(1)
          plot(x = 1:20,
               y = as.numeric(GPCM_theta_bin_1), 
               type = "h",
               lwd = 3,
               ylab = "GPCM substantive trait bias",
               xlab = "True substantive trait",
               xaxt = "n",
               ylim = y_axis,
               col = GPCM_colours_1)
          axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))
          
          screen(2)
          plot(x = 1:20,
               y = as.numeric(GPCM_theta_bin_2), 
               type = "h",
               lwd = 3,
               ylab = "GPCM substantive trait bias",
               xlab = "True substantive trait",
               xaxt = "n",
               ylim = y_axis,
               col = GPCM_colours_2)
          axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))
          
          close.screen(all.screens = TRUE)
          
          split.screen(c(1,2))
          screen(1)
          plot(x = 1:20,
               y = as.numeric(MNRM_theta_bin_1), 
               type = "h",
               lwd = 3,
               ylab = "MNRM substantive trait bias",
               xlab = "True substantive trait",
               xaxt = "n",
               ylim = y_axis,
               col = MNRM_colours_1)
          axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))
          
          screen(2)
          plot(x = 1:20,
               y = as.numeric(MNRM_theta_bin_2), 
               type = "h",
               lwd = 3,
               ylab = "MNRM substantive trait bias",
               xlab = "True substantive trait",
               xaxt = "n",
               ylim = y_axis,
               col = MNRM_colours_2)
          axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))
          
          close.screen(all.screens = TRUE)
          
          split.screen(c(1,2))
          screen(1)
          plot(x = 1:20,
               y = as.numeric(Tree_theta_bin_1), 
               type = "h",
               lwd = 3,
               ylab = "IRTree substantive trait bias",
               xlab = "True substantive trait",
               xaxt = "n",
               ylim = y_axis,
               col = Tree_colours_1)
          axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))
          
          screen(2)
          plot(x = 1:20,
               y = as.numeric(Tree_theta_bin_2), 
               type = "h",
               lwd = 3,
               ylab = "IRTree substantive trait bias",
               xlab = "True substantive trait",
               xaxt = "n",
               ylim = y_axis,
               col = Tree_colours_2)
          axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))
          
          close.screen(all.screens = TRUE)
        }
        
        #Generate a blank plot, suppress x axis and add custom one
        plot(x = 1:20,
             y = as.numeric(MNRM_theta_bin_1), 
             type = "n",
             lwd = 3,
             ylab = "Substantive trait bias",
             xlab = "True substantive trait",
             xaxt = "n",
             ylim = y_axis)
        
        axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))
        
        if(matrix_col != 2){
          #Generate GPCM lines
          if(!(k == 1 & generating_model_values == TRUE)){
            lines(x = 1:20,
                  y = as.numeric(GPCM_theta_bin_1),
                  type = "l",
                  lwd = 3,
                  lty = "solid",
                  col = "black",
                  ylim = y_axis)
            
            lines(x = 1:20,
                  y = as.numeric(GPCM_theta_bin_2),
                  type = "l",
                  lwd = 3,
                  lty = "dotted",
                  col = "black",
                  ylim = y_axis)
          }
        }
        
        #Generate MNRM lines
        if(!(k == 2 & generating_model_values == TRUE)){
          lines(x = 1:20,
                as.numeric(MNRM_theta_bin_1),
                type = "l", 
                lwd = 3,
                lty = "solid",
                ylim = y_axis,
                col = "red")
          
          lines(x = 1:20,
                y = as.numeric(MNRM_theta_bin_2),
                type = "l", 
                lwd = 3,
                lty = "dotted", 
                col = "red",
                ylim = y_axis)
        }
        
        #Generate IRTree lines
        if(!(k == 3 & generating_model_values == TRUE)){
          lines(x = 1:20,
                y = as.numeric(Tree_theta_bin_1),
                type = "l",
                lwd = 3,
                lty = "solid",
                ylim = y_axis,
                col = "blue")
          
          lines(x = 1:20,
                y = as.numeric(Tree_theta_bin_2),
                type = "l",
                lwd = 3,
                lty = "dotted",
                col = "blue",
                ylim = y_axis)
        }
        
        if(split_plots == TRUE & !(k == 3 & generating_model_values == TRUE)){
          split.screen(c(1,2))
          screen(1)
          plot(x = 1:20,
               y = as.numeric(GPCM_theta_bin_1), 
               type = "l",
               lwd = 3,
               ylab = "Substantive trait bias",
               xlab = "True substantive trait",
               xaxt = "n",
               ylim = y_axis,
               col = "black")
          
          lines(x = 1:20,
                y = as.numeric(Tree_theta_bin_1),
                type = "l",
                lwd = 3,
                lty = "solid",
                ylim = y_axis,
                col = "blue")
          
          axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))
          legend("topright", legend = c("GPCM", "IRTree"), 
                 col = c("black", "blue"), lty = "solid", lwd = 3)
          
          screen(2)
          plot(x = 1:20,
               y = as.numeric(GPCM_theta_bin_2), 
               type = "l",
               lwd = 3,
               ylab = "Substantive trait bias",
               xlab = "True substantive trait",
               xaxt = "n",
               ylim = y_axis,
               col = "black",
               lty = "solid")
          
          lines(x = 1:20,
                y = as.numeric(Tree_theta_bin_2),
                type = "l",
                lwd = 3,
                lty = "solid",
                ylim = y_axis,
                col = "blue")
          
          axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))
          legend("topright", legend = c("GPCM", "IRTree"), 
                 col = c("black", "blue"), lty = "solid", lwd = 3)
          
          close.screen(all.screens = TRUE)
        }
        
        if(split_plots == TRUE & !(k == 2 & generating_model_values == TRUE)){
          split.screen(c(1,2))
          screen(1)
          plot(x = 1:20,
               y = as.numeric(GPCM_theta_bin_1), 
               type = "l",
               lwd = 3,
               ylab = "Substantive trait bias",
               xlab = "True substantive trait",
               xaxt = "n",
               ylim = y_axis,
               col = "black")
          
          lines(x = 1:20,
                y = as.numeric(MNRM_theta_bin_1),
                type = "l",
                lwd = 3,
                lty = "solid",
                ylim = y_axis,
                col = "red")
          
          axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))
          legend("topright", legend = c("GPCM", "MNRM"), 
                 col = c("black", "red"), lty = "solid", lwd = 3)
          
          screen(2)
          plot(x = 1:20,
               y = as.numeric(GPCM_theta_bin_2), 
               type = "l",
               lwd = 3,
               ylab = "Substantive trait bias",
               xlab = "True substantive trait",
               xaxt = "n",
               ylim = y_axis,
               col = "black",
               lty = "solid")
          
          lines(x = 1:20,
                y = as.numeric(MNRM_theta_bin_2),
                type = "l",
                lwd = 3,
                lty = "solid",
                ylim = y_axis,
                col = "red")
          
          axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))
          legend("topright", legend = c("GPCM", "MNRM"), 
                 col = c("black", "red"), lty = "solid", lwd = 3)
          
          close.screen(all.screens = TRUE)
        }
      }
    }
    
    
    #Only generate heatmaps and ERS unidimensional plots when not in the GPCM
    #condition
    if(k != 1 & plot_unidimensional_ERS == TRUE){
      
      #Reset graphical parameters
      par(old.par)
      
      #Same as above, but this time put the true ERS on the x-axis
      #Remove all ERS values above 2 or below -2
      MNRM_df_plot <- MNRM_df[is.na(MNRM_df[,4]) == FALSE,]
      Tree_df_plot <- Tree_df[is.na(MNRM_df[,4]) == FALSE,]
      GPCM_df_plot <- GPCM_df[is.na(MNRM_df[,4]) == FALSE,]
      bins <- names(table(true_bins[,4]))
      
      #Initialize bins and plot colours
      MNRM_ERS_bin <- Tree_ERS_bin <- GPCM_ERS_bin <- 
        matrix(NA, nrow = 1, ncol = 20)
      MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
      for(i in 1:20){
        MNRM_ERS_bin[,i] <- 
          mean(MNRM_df_plot[MNRM_df_plot[,4] == bins[i], 1])
        Tree_ERS_bin[,i] <- 
          mean(Tree_df_plot[Tree_df_plot[,4] == bins[i], 1])
        GPCM_ERS_bin[,i] <- 
          mean(GPCM_df_plot[GPCM_df_plot[,3] == bins[i], 1])
        
        if(Tree_ERS_bin[,i] < 0){
          Tree_colours[i] <- "blue"
        } else{
          Tree_colours[i] <- "red"
        }
        
        if(GPCM_ERS_bin[,i] < 0){
          GPCM_colours[i] <- "blue"
        } else{
          GPCM_colours[i] <- "red"
        }
        
        if(MNRM_ERS_bin[,i] < 0){
          MNRM_colours[i] <- "blue"
        } else{
          MNRM_colours[i] <- "red"
        }
      }
      
      
      split.screen(c(1,3))
      
      screen(1)
      #make plots, again with custom x axis
      plot(as.numeric(GPCM_ERS_bin), type = "l",lwd = 3,
           col = GPCM_colours,
           ylab = y_axis_names[1],
           xlab = "True ERS trait",
           xaxt = "n",
           ylim = y_axis)
      axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))
      
      screen(2)
      plot(as.numeric(Tree_ERS_bin), type = "l",lwd = 3,
           col = Tree_colours,
           ylab = y_axis_names[2],
           xlab = "True ERS trait",
           xaxt = "n",
           ylim = y_axis)
      axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))
      
      screen(3)
      plot(as.numeric(MNRM_ERS_bin), type = "l",lwd = 3,
           col = MNRM_colours,
           ylab = y_axis_names[3],
           xlab = "True ERS trait",
           xaxt = "n",
           ylim = y_axis)
      axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))
    }
    close.screen(all.screens = TRUE)
    
    if(k != 1 & plot_heatmaps == TRUE){ 
      #Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
      #Consider if you want to remove these or bundle them into the highest/lowest bins
      Tree_df <- Tree_df[is.na(Tree_df[,3]) == FALSE & 
                           is.na(Tree_df[,4]) == FALSE,]
      MNRM_df <- MNRM_df[is.na(MNRM_df[,3]) == FALSE & 
                           is.na(MNRM_df[,4]) == FALSE,]
      GPCM_df <- GPCM_df[is.na(GPCM_df[,2]) == FALSE & 
                           is.na(GPCM_df[,3]) == FALSE,]
      
      bins <- names(table(true_bins[,3]))
      
      #Obtain the average difference for each of the 400 bins that were created
      #Initialzie matrices and set row and colom names
      binned_tree <- binned_MNRM <- binned_GPCM <- 
        data.frame(matrix(NA, nrow = 20, ncol = 20))
      rownames(binned_MNRM) <- rownames(binned_GPCM) <- 
        rownames(binned_tree) <- bins
      colnames(binned_MNRM) <- colnames(binned_GPCM) <- 
        colnames(binned_tree) <- bins
      
      #Calculate the mean bias for every bin
      for(i in 1:20){
        for(j in 1:20){
          binned_MNRM[i,j] <- mean(MNRM_df[
            MNRM_df[,3] == bins[j] & MNRM_df[,4] == bins[i], 1])
          binned_tree[i,j] <- mean(Tree_df[
            Tree_df[,3] == bins[j] & Tree_df[,4] == bins[i], 1])
          binned_GPCM[i,j] <- mean(GPCM_df[
            GPCM_df[,2] == bins[j] & GPCM_df[,3] == bins[i], 1])
        }
      }
      
      #Set the heatmap colours
      colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))
      
      my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)
      
      #Adjust margins so that axis labels are more readable
      par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))
      
      my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)
      
      #Plot the actual heatmaps
      heatmap(as.matrix(binned_GPCM), Rowv = NA, Colv = NA, scale = "none", 
              col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
              RowSideColors = my_pallete2)
      heatmap(as.matrix(binned_tree), Rowv = NA, Colv = NA, scale = "none",
              col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
              RowSideColors = my_pallete2)
      heatmap(as.matrix(binned_MNRM), Rowv = NA, Colv = NA, scale = "none",
              col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
              RowSideColors = my_pallete2)
    }
    
    
  }
  
}


#debug(heatmaps)
heatmaps(true_values = TRUE, generating_model_values = FALSE,
         plot_heatmaps = FALSE, condition = "GPCM",
         plot_unidimensional_theta = TRUE,
         plot_unidimensional_ERS = FALSE,
         y_axis = c(-1, 1),
         plot_GPCM_only = TRUE)

#Stupid rstudio limit of 100 plots make you do this in steps
heatmaps(true_values = TRUE, generating_model_values = FALSE,
         plot_heatmaps = FALSE, condition = "MNRM",
         plot_unidimensional_theta = TRUE,
         plot_unidimensional_ERS = FALSE,
         y_axis = c(-1, 1))

heatmaps(true_values = TRUE, generating_model_values = FALSE,
         plot_heatmaps = FALSE, condition = "MNRM",
         plot_unidimensional_theta = FALSE,
         plot_unidimensional_theta_ERS = TRUE,
         plot_unidimensional_ERS = FALSE,
         y_axis = c(-1, 1))

heatmaps(true_values = FALSE, generating_model_values = TRUE,
         plot_heatmaps = FALSE, condition = "MNRM",
         plot_unidimensional_theta = TRUE,
         plot_unidimensional_ERS = FALSE,
         y_axis = c(-1, 1),
         split_plots = TRUE)

heatmaps(true_values = TRUE, generating_model_values = FALSE,
         plot_heatmaps = FALSE, condition = "MNRM",
         plot_unidimensional_theta = FALSE,
         plot_unidimensional_ERS = TRUE,
         y_axis = c(-1, 1))

heatmaps(true_values = TRUE, generating_model_values = FALSE,
         plot_heatmaps = TRUE, condition = "MNRM",
         plot_unidimensional_theta = FALSE,
         plot_unidimensional_ERS = FALSE,
         y_axis = c(-1, 1))

#IRTree plots
heatmaps(true_values = TRUE, generating_model_values = FALSE,
         plot_heatmaps = FALSE, condition = "IRTree",
         plot_unidimensional_theta = TRUE,
         plot_unidimensional_ERS = FALSE,
         y_axis = c(-1, 1))

heatmaps(true_values = FALSE, generating_model_values = TRUE,
         plot_heatmaps = FALSE, condition = "IRTree",
         plot_unidimensional_theta = TRUE,
         plot_unidimensional_ERS = FALSE,
         y_axis = c(-1, 1),
         split_plots = TRUE)

heatmaps(true_values = TRUE, generating_model_values = FALSE,
         plot_heatmaps = FALSE, condition = "IRTree",
         plot_unidimensional_theta = FALSE,
         plot_unidimensional_ERS = TRUE,
         y_axis = c(-1, 1))

heatmaps(true_values = TRUE, generating_model_values = FALSE,
         plot_heatmaps = TRUE, condition = "IRTree",
         plot_unidimensional_theta = FALSE,
         plot_unidimensional_ERS = FALSE,
         y_axis = c(-1, 1))





################################################################################
#Everything below here not included in paper, just experimentation
################################################################################






############################################################################
#All possible plots on six pages (6 plots per page)
############################################################################

#-1 ERS condition MNRM no shift
#Only group 2 included here!
#Bins based on MNRM data, threshold -1, 0, 1 with +1 ERS mean for group 2
MNRM_true_traits_min1ERS_0 <- MNRM_Tree_est_min1ERS_0 <- MNRM_MNRM_est_min1ERS_0 <- 
  matrix(NA, nrow = 250000, ncol = 2)
MNRM_GPCM_est_min1ERS_0 <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  MNRM_true_traits_min1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_min1[[i]][["True_traits"]][501:1000,]
  MNRM_MNRM_est_min1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_min1[[i]][["MNRM_traits"]][501:1000,]
  MNRM_Tree_est_min1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_min1[[i]][["Tree_traits"]][501:1000,]
  MNRM_GPCM_est_min1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_min1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
MNRM_MNRM_est_min1ERS_0 <- MNRM_MNRM_est_min1ERS_0 - MNRM_true_traits_min1ERS_0
MNRM_Tree_est_min1ERS_0 <- MNRM_Tree_est_min1ERS_0 - MNRM_true_traits_min1ERS_0
MNRM_GPCM_est_min1ERS_0 <- MNRM_GPCM_est_min1ERS_0 - MNRM_true_traits_min1ERS_0[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
MNRM_true_bins <- data.frame(True_theta = MNRM_true_traits_min1ERS_0[,1], True_ERS = MNRM_true_traits_min1ERS_0[,2],
                             theta_bin = cut(MNRM_true_traits_min1ERS_0[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(MNRM_true_traits_min1ERS_0[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
MNRM_Tree_df <- data.frame(Tree_theta_bias = MNRM_Tree_est_min1ERS_0[,1], Tree_ERS_bias = MNRM_Tree_est_min1ERS_0[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])
MNRM_GPCM_df <- data.frame(GPCM_theta_bias = MNRM_GPCM_est_min1ERS_0[,1], theta_bin = MNRM_true_bins[,3],
                           ERS_bin = MNRM_true_bins[,4])
MNRM_MNRM_df <- data.frame(MNRM_theta_bias = MNRM_MNRM_est_min1ERS_0[,1], MNRM_ERS_bias = MNRM_MNRM_est_min1ERS_0[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])

par(old.par)

#Make 1-dimensional graphs based on just the theta and just the ERS bin
MNRM_MNRM_df_plot <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
MNRM_Tree_df_plot <- MNRM_Tree_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
MNRM_GPCM_df_plot <- MNRM_GPCM_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
bins <- names(table(MNRM_true_bins[,3]))

MNRM_theta_bin <- Tree_theta_bin <- GPCM_theta_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_theta_bin[,i] <- 
    mean(MNRM_MNRM_df_plot[MNRM_MNRM_df_plot[,3] == bins[i], 1])
  Tree_theta_bin[,i] <- 
    mean(MNRM_Tree_df_plot[MNRM_Tree_df_plot[,3] == bins[i], 1])
  GPCM_theta_bin[,i] <- 
    mean(MNRM_GPCM_df_plot[MNRM_GPCM_df_plot[,2] == bins[i], 1])
  
  if(Tree_theta_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_theta_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_theta_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_theta_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_theta_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_theta_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Make 1-dimensional graphs based on just the theta and just the ERS bin
MNRM_MNRM_df_plot <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_Tree_df_plot <- MNRM_Tree_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_GPCM_df_plot <- MNRM_GPCM_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
bins <- names(table(MNRM_true_bins[,4]))

MNRM_ERS_bin <- Tree_ERS_bin <- GPCM_ERS_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_ERS_bin[,i] <- 
    mean(MNRM_MNRM_df_plot[MNRM_MNRM_df_plot[,4] == bins[i], 1])
  Tree_ERS_bin[,i] <- 
    mean(MNRM_Tree_df_plot[MNRM_Tree_df_plot[,4] == bins[i], 1])
  GPCM_ERS_bin[,i] <- 
    mean(MNRM_GPCM_df_plot[MNRM_GPCM_df_plot[,3] == bins[i], 1])
  
  if(Tree_ERS_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_ERS_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_ERS_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_ERS_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_ERS_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_ERS_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))


#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
MNRM_Tree_df <- MNRM_Tree_df[is.na(MNRM_Tree_df[,3]) == FALSE & 
                               is.na(MNRM_Tree_df[,4]) == FALSE,]
MNRM_MNRM_df <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,3]) == FALSE & 
                               is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_GPCM_df <- MNRM_GPCM_df[is.na(MNRM_GPCM_df[,2]) == FALSE & 
                               is.na(MNRM_GPCM_df[,3]) == FALSE,]

bins <- names(table(MNRM_true_bins[,3]))

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
MNRM_binned_tree_min1ERS_0 <- MNRM_binned_MNRM_min1ERS_0 <- MNRM_binned_GPCM_min1ERS_0 <- 
  data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(MNRM_binned_MNRM_min1ERS_0) <- rownames(MNRM_binned_GPCM_min1ERS_0) <- 
  rownames(MNRM_binned_tree_min1ERS_0) <- bins
colnames(MNRM_binned_MNRM_min1ERS_0) <- colnames(MNRM_binned_GPCM_min1ERS_0) <- 
  colnames(MNRM_binned_tree_min1ERS_0) <- bins

for(i in 1:20){
  for(j in 1:20){
    MNRM_binned_MNRM_min1ERS_0[i,j] <- mean(MNRM_MNRM_df[MNRM_MNRM_df[,3] == bins[j] &
                                                           MNRM_MNRM_df[,4] == bins[i], 1])
    MNRM_binned_tree_min1ERS_0[i,j] <- mean(MNRM_Tree_df[MNRM_Tree_df[,3] == bins[j] &
                                                           MNRM_Tree_df[,4] == bins[i], 1])
    MNRM_binned_GPCM_min1ERS_0[i,j] <- mean(MNRM_GPCM_df[MNRM_GPCM_df[,2] == bins[j] &
                                                           MNRM_GPCM_df[,3] == bins[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(MNRM_binned_GPCM_min1ERS_0), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_tree_min1ERS_0), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_MNRM_min1ERS_0), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)










#-1 ERS condition MNRM with shift
#Only group 2 included here!
#Bins based on MNRM data, threshold -1, 0, 1 with +1 ERS mean for group 2
#Intialize matrices
MNRM_true_traits_min1ERS_1 <- MNRM_Tree_est_min1ERS_1 <- MNRM_MNRM_est_min1ERS_1 <- 
  matrix(NA, nrow = 250000, ncol = 2)
MNRM_GPCM_est_min1ERS_1 <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  MNRM_true_traits_min1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_min1[[i]][["True_traits"]][501:1000,]
  MNRM_MNRM_est_min1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_min1[[i]][["MNRM_traits"]][501:1000,]
  MNRM_Tree_est_min1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_min1[[i]][["Tree_traits"]][501:1000,]
  MNRM_GPCM_est_min1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_min1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
MNRM_MNRM_est_min1ERS_1 <- MNRM_MNRM_est_min1ERS_1 - MNRM_true_traits_min1ERS_1
MNRM_Tree_est_min1ERS_1 <- MNRM_Tree_est_min1ERS_1 - MNRM_true_traits_min1ERS_1
MNRM_GPCM_est_min1ERS_1 <- MNRM_GPCM_est_min1ERS_1 - MNRM_true_traits_min1ERS_1[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
MNRM_true_bins <- data.frame(True_theta = MNRM_true_traits_min1ERS_1[,1], True_ERS = MNRM_true_traits_min1ERS_1[,2],
                             theta_bin = cut(MNRM_true_traits_min1ERS_1[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(MNRM_true_traits_min1ERS_1[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
MNRM_Tree_df <- data.frame(Tree_theta_bias = MNRM_Tree_est_min1ERS_1[,1], Tree_ERS_bias = MNRM_Tree_est_min1ERS_1[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])
MNRM_GPCM_df <- data.frame(GPCM_theta_bias = MNRM_GPCM_est_min1ERS_1[,1], theta_bin = MNRM_true_bins[,3],
                           ERS_bin = MNRM_true_bins[,4])
MNRM_MNRM_df <- data.frame(MNRM_theta_bias = MNRM_MNRM_est_min1ERS_1[,1], MNRM_ERS_bias = MNRM_MNRM_est_min1ERS_1[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])

par(old.par)

#Make 1-dimensional graphs based on just the theta and just the ERS bin
MNRM_MNRM_df_plot <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
MNRM_Tree_df_plot <- MNRM_Tree_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
MNRM_GPCM_df_plot <- MNRM_GPCM_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
bins <- names(table(MNRM_true_bins[,3]))

MNRM_theta_bin <- Tree_theta_bin <- GPCM_theta_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_theta_bin[,i] <- 
    mean(MNRM_MNRM_df_plot[MNRM_MNRM_df_plot[,3] == bins[i], 1])
  Tree_theta_bin[,i] <- 
    mean(MNRM_Tree_df_plot[MNRM_Tree_df_plot[,3] == bins[i], 1])
  GPCM_theta_bin[,i] <- 
    mean(MNRM_GPCM_df_plot[MNRM_GPCM_df_plot[,2] == bins[i], 1])
  
  if(Tree_theta_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_theta_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_theta_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_theta_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_theta_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_theta_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Make 1-dimensional graphs based on just the theta and just the ERS bin
MNRM_MNRM_df_plot <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_Tree_df_plot <- MNRM_Tree_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_GPCM_df_plot <- MNRM_GPCM_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
bins <- names(table(MNRM_true_bins[,4]))

MNRM_ERS_bin <- Tree_ERS_bin <- GPCM_ERS_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_ERS_bin[,i] <- 
    mean(MNRM_MNRM_df_plot[MNRM_MNRM_df_plot[,4] == bins[i], 1])
  Tree_ERS_bin[,i] <- 
    mean(MNRM_Tree_df_plot[MNRM_Tree_df_plot[,4] == bins[i], 1])
  GPCM_ERS_bin[,i] <- 
    mean(MNRM_GPCM_df_plot[MNRM_GPCM_df_plot[,3] == bins[i], 1])
  
  if(Tree_ERS_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_ERS_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_ERS_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_ERS_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_ERS_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_ERS_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
MNRM_Tree_df <- MNRM_Tree_df[is.na(MNRM_Tree_df[,3]) == FALSE & 
                               is.na(MNRM_Tree_df[,4]) == FALSE,]
MNRM_MNRM_df <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,3]) == FALSE & 
                               is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_GPCM_df <- MNRM_GPCM_df[is.na(MNRM_GPCM_df[,2]) == FALSE & 
                               is.na(MNRM_GPCM_df[,3]) == FALSE,]

bins <- names(table(MNRM_true_bins[,3]))

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
MNRM_binned_tree_min1ERS_1 <- MNRM_binned_MNRM_min1ERS_1 <- MNRM_binned_GPCM_min1ERS_1 <- 
  data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(MNRM_binned_MNRM_min1ERS_1) <- rownames(MNRM_binned_GPCM_min1ERS_1) <- 
  rownames(MNRM_binned_tree_min1ERS_1) <- bins
colnames(MNRM_binned_MNRM_min1ERS_1) <- colnames(MNRM_binned_GPCM_min1ERS_1) <- 
  colnames(MNRM_binned_tree_min1ERS_1) <- bins

for(i in 1:20){
  for(j in 1:20){
    MNRM_binned_MNRM_min1ERS_1[i,j] <- mean(MNRM_MNRM_df[MNRM_MNRM_df[,3] == bins[j] &
                                                           MNRM_MNRM_df[,4] == bins[i], 1])
    MNRM_binned_tree_min1ERS_1[i,j] <- mean(MNRM_Tree_df[MNRM_Tree_df[,3] == bins[j] &
                                                           MNRM_Tree_df[,4] == bins[i], 1])
    MNRM_binned_GPCM_min1ERS_1[i,j] <- mean(MNRM_GPCM_df[MNRM_GPCM_df[,2] == bins[j] &
                                                           MNRM_GPCM_df[,3] == bins[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(MNRM_binned_GPCM_min1ERS_1), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_tree_min1ERS_1), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_MNRM_min1ERS_1), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)










#0 ERS condition MNRM with no shift
#Only group 2 included here!
#Bins based on MNRM data, threshold -1, 0, 1 with +1 ERS mean for group 2
MNRM_true_traits_0ERS_0 <- MNRM_Tree_est_0ERS_0 <- MNRM_MNRM_est_0ERS_0 <- 
  matrix(NA, nrow = 250000, ncol = 2)
MNRM_GPCM_est_0ERS_0 <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  MNRM_true_traits_0ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20[[i]][["True_traits"]][501:1000,]
  MNRM_MNRM_est_0ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20[[i]][["MNRM_traits"]][501:1000,]
  MNRM_Tree_est_0ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20[[i]][["Tree_traits"]][501:1000,]
  MNRM_GPCM_est_0ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
MNRM_MNRM_est_0ERS_0 <- MNRM_MNRM_est_0ERS_0 - MNRM_true_traits_0ERS_0
MNRM_Tree_est_0ERS_0 <- MNRM_Tree_est_0ERS_0 - MNRM_true_traits_0ERS_0
MNRM_GPCM_est_0ERS_0 <- MNRM_GPCM_est_0ERS_0 - MNRM_true_traits_0ERS_0[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
MNRM_true_bins <- data.frame(True_theta = MNRM_true_traits_0ERS_0[,1], True_ERS = MNRM_true_traits_0ERS_0[,2],
                             theta_bin = cut(MNRM_true_traits_0ERS_0[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(MNRM_true_traits_0ERS_0[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
MNRM_Tree_df <- data.frame(Tree_theta_bias = MNRM_Tree_est_0ERS_0[,1], Tree_ERS_bias = MNRM_Tree_est_0ERS_0[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])
MNRM_GPCM_df <- data.frame(GPCM_theta_bias = MNRM_GPCM_est_0ERS_0[,1], theta_bin = MNRM_true_bins[,3],
                           ERS_bin = MNRM_true_bins[,4])
MNRM_MNRM_df <- data.frame(MNRM_theta_bias = MNRM_MNRM_est_0ERS_0[,1], MNRM_ERS_bias = MNRM_MNRM_est_0ERS_0[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])


par(old.par)

#Make 1-dimensional graphs based on just the theta and just the ERS bin
MNRM_MNRM_df_plot <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
MNRM_Tree_df_plot <- MNRM_Tree_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
MNRM_GPCM_df_plot <- MNRM_GPCM_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
bins <- names(table(MNRM_true_bins[,3]))

MNRM_theta_bin <- Tree_theta_bin <- GPCM_theta_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_theta_bin[,i] <- 
    mean(MNRM_MNRM_df_plot[MNRM_MNRM_df_plot[,3] == bins[i], 1])
  Tree_theta_bin[,i] <- 
    mean(MNRM_Tree_df_plot[MNRM_Tree_df_plot[,3] == bins[i], 1])
  GPCM_theta_bin[,i] <- 
    mean(MNRM_GPCM_df_plot[MNRM_GPCM_df_plot[,2] == bins[i], 1])
  
  if(Tree_theta_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_theta_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_theta_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_theta_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_theta_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_theta_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Make 1-dimensional graphs based on just the theta and just the ERS bin
MNRM_MNRM_df_plot <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_Tree_df_plot <- MNRM_Tree_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_GPCM_df_plot <- MNRM_GPCM_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
bins <- names(table(MNRM_true_bins[,4]))

MNRM_ERS_bin <- Tree_ERS_bin <- GPCM_ERS_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_ERS_bin[,i] <- 
    mean(MNRM_MNRM_df_plot[MNRM_MNRM_df_plot[,4] == bins[i], 1])
  Tree_ERS_bin[,i] <- 
    mean(MNRM_Tree_df_plot[MNRM_Tree_df_plot[,4] == bins[i], 1])
  GPCM_ERS_bin[,i] <- 
    mean(MNRM_GPCM_df_plot[MNRM_GPCM_df_plot[,3] == bins[i], 1])
  
  if(Tree_ERS_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_ERS_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_ERS_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_ERS_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_ERS_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_ERS_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))


#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
MNRM_Tree_df <- MNRM_Tree_df[is.na(MNRM_Tree_df[,3]) == FALSE & 
                               is.na(MNRM_Tree_df[,4]) == FALSE,]
MNRM_MNRM_df <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,3]) == FALSE & 
                               is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_GPCM_df <- MNRM_GPCM_df[is.na(MNRM_GPCM_df[,2]) == FALSE & 
                               is.na(MNRM_GPCM_df[,3]) == FALSE,]

bins <- names(table(MNRM_true_bins[,3]))

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
MNRM_binned_tree_0ERS_0 <- MNRM_binned_MNRM_0ERS_0 <- MNRM_binned_GPCM_0ERS_0 <- 
  data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(MNRM_binned_MNRM_0ERS_0) <- rownames(MNRM_binned_GPCM_0ERS_0) <- 
  rownames(MNRM_binned_tree_0ERS_0) <- bins
colnames(MNRM_binned_MNRM_0ERS_0) <- colnames(MNRM_binned_GPCM_0ERS_0) <- 
  colnames(MNRM_binned_tree_0ERS_0) <- bins

for(i in 1:20){
  for(j in 1:20){
    MNRM_binned_MNRM_0ERS_0[i,j] <- mean(MNRM_MNRM_df[MNRM_MNRM_df[,3] == bins[j] &
                                                        MNRM_MNRM_df[,4] == bins[i], 1])
    MNRM_binned_tree_0ERS_0[i,j] <- mean(MNRM_Tree_df[MNRM_Tree_df[,3] == bins[j] &
                                                        MNRM_Tree_df[,4] == bins[i], 1])
    MNRM_binned_GPCM_0ERS_0[i,j] <- mean(MNRM_GPCM_df[MNRM_GPCM_df[,2] == bins[j] &
                                                        MNRM_GPCM_df[,3] == bins[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(MNRM_binned_GPCM_0ERS_0), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_tree_0ERS_0), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_MNRM_0ERS_0), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)










#0 ERS condition MNRM with shift
#Only group 2 included here!
#Bins based on MNRM data, threshold -1, 0, 1 with +1 ERS mean for group 2
#Intialize matrices
MNRM_true_traits_0ERS_1 <- MNRM_Tree_est_0ERS_1 <- MNRM_MNRM_est_0ERS_1 <- 
  matrix(NA, nrow = 250000, ncol = 2)
MNRM_GPCM_est_0ERS_1 <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  MNRM_true_traits_0ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20[[i]][["True_traits"]][501:1000,]
  MNRM_MNRM_est_0ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20[[i]][["MNRM_traits"]][501:1000,]
  MNRM_Tree_est_0ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20[[i]][["Tree_traits"]][501:1000,]
  MNRM_GPCM_est_0ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
MNRM_MNRM_est_0ERS_1 <- MNRM_MNRM_est_0ERS_1 - MNRM_true_traits_0ERS_1
MNRM_Tree_est_0ERS_1 <- MNRM_Tree_est_0ERS_1 - MNRM_true_traits_0ERS_1
MNRM_GPCM_est_0ERS_1 <- MNRM_GPCM_est_0ERS_1 - MNRM_true_traits_0ERS_1[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
MNRM_true_bins <- data.frame(True_theta = MNRM_true_traits_0ERS_1[,1], True_ERS = MNRM_true_traits_0ERS_1[,2],
                             theta_bin = cut(MNRM_true_traits_0ERS_1[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(MNRM_true_traits_0ERS_1[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
MNRM_Tree_df <- data.frame(Tree_theta_bias = MNRM_Tree_est_0ERS_1[,1], Tree_ERS_bias = MNRM_Tree_est_0ERS_1[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])
MNRM_GPCM_df <- data.frame(GPCM_theta_bias = MNRM_GPCM_est_0ERS_1[,1], theta_bin = MNRM_true_bins[,3],
                           ERS_bin = MNRM_true_bins[,4])
MNRM_MNRM_df <- data.frame(MNRM_theta_bias = MNRM_MNRM_est_0ERS_1[,1], MNRM_ERS_bias = MNRM_MNRM_est_0ERS_1[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])



par(old.par)

#Make 1-dimensional graphs based on just the theta and just the ERS bin
MNRM_MNRM_df_plot <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
MNRM_Tree_df_plot <- MNRM_Tree_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
MNRM_GPCM_df_plot <- MNRM_GPCM_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
bins <- names(table(MNRM_true_bins[,3]))

MNRM_theta_bin <- Tree_theta_bin <- GPCM_theta_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_theta_bin[,i] <- 
    mean(MNRM_MNRM_df_plot[MNRM_MNRM_df_plot[,3] == bins[i], 1])
  Tree_theta_bin[,i] <- 
    mean(MNRM_Tree_df_plot[MNRM_Tree_df_plot[,3] == bins[i], 1])
  GPCM_theta_bin[,i] <- 
    mean(MNRM_GPCM_df_plot[MNRM_GPCM_df_plot[,2] == bins[i], 1])
  
  if(Tree_theta_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_theta_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_theta_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_theta_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_theta_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_theta_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Make 1-dimensional graphs based on just the theta and just the ERS bin
MNRM_MNRM_df_plot <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_Tree_df_plot <- MNRM_Tree_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_GPCM_df_plot <- MNRM_GPCM_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
bins <- names(table(MNRM_true_bins[,4]))

MNRM_ERS_bin <- Tree_ERS_bin <- GPCM_ERS_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_ERS_bin[,i] <- 
    mean(MNRM_MNRM_df_plot[MNRM_MNRM_df_plot[,4] == bins[i], 1])
  Tree_ERS_bin[,i] <- 
    mean(MNRM_Tree_df_plot[MNRM_Tree_df_plot[,4] == bins[i], 1])
  GPCM_ERS_bin[,i] <- 
    mean(MNRM_GPCM_df_plot[MNRM_GPCM_df_plot[,3] == bins[i], 1])
  
  if(Tree_ERS_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_ERS_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_ERS_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_ERS_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_ERS_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_ERS_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
MNRM_Tree_df <- MNRM_Tree_df[is.na(MNRM_Tree_df[,3]) == FALSE & 
                               is.na(MNRM_Tree_df[,4]) == FALSE,]
MNRM_MNRM_df <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,3]) == FALSE & 
                               is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_GPCM_df <- MNRM_GPCM_df[is.na(MNRM_GPCM_df[,2]) == FALSE & 
                               is.na(MNRM_GPCM_df[,3]) == FALSE,]

bins <- names(table(MNRM_true_bins[,3]))

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
MNRM_binned_tree_0ERS_1 <- MNRM_binned_MNRM_0ERS_1 <- MNRM_binned_GPCM_0ERS_1 <- 
  data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(MNRM_binned_MNRM_0ERS_1) <- rownames(MNRM_binned_GPCM_0ERS_1) <- 
  rownames(MNRM_binned_tree_0ERS_1) <- bins
colnames(MNRM_binned_MNRM_0ERS_1) <- colnames(MNRM_binned_GPCM_0ERS_1) <- 
  colnames(MNRM_binned_tree_0ERS_1) <- bins

for(i in 1:20){
  for(j in 1:20){
    MNRM_binned_MNRM_0ERS_1[i,j] <- mean(MNRM_MNRM_df[MNRM_MNRM_df[,3] == bins[j] &
                                                        MNRM_MNRM_df[,4] == bins[i], 1])
    MNRM_binned_tree_0ERS_1[i,j] <- mean(MNRM_Tree_df[MNRM_Tree_df[,3] == bins[j] &
                                                        MNRM_Tree_df[,4] == bins[i], 1])
    MNRM_binned_GPCM_0ERS_1[i,j] <- mean(MNRM_GPCM_df[MNRM_GPCM_df[,2] == bins[j] &
                                                        MNRM_GPCM_df[,3] == bins[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(MNRM_binned_GPCM_0ERS_1), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_tree_0ERS_1), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_MNRM_0ERS_1), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)











#########################################################################
#Only group 2 included here!
#Bins based on MNRM data, threshold -1, 0, 1 with +1 ERS mean for group 2
#Intialize matrices
MNRM_true_traits_1ERS_0 <- MNRM_Tree_est_1ERS_0 <- MNRM_MNRM_est_1ERS_0 <- 
  matrix(NA, nrow = 250000, ncol = 2)
MNRM_GPCM_est_1ERS_0 <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  MNRM_true_traits_1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_1[[i]][["True_traits"]][501:1000,]
  MNRM_MNRM_est_1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_1[[i]][["MNRM_traits"]][501:1000,]
  MNRM_Tree_est_1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_1[[i]][["Tree_traits"]][501:1000,]
  MNRM_GPCM_est_1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
MNRM_MNRM_est_1ERS_0 <- MNRM_MNRM_est_1ERS_0 - MNRM_true_traits_1ERS_0
MNRM_Tree_est_1ERS_0 <- MNRM_Tree_est_1ERS_0 - MNRM_true_traits_1ERS_0
MNRM_GPCM_est_1ERS_0 <- MNRM_GPCM_est_1ERS_0 - MNRM_true_traits_1ERS_0[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
MNRM_true_bins <- data.frame(True_theta = MNRM_true_traits_1ERS_0[,1], True_ERS = MNRM_true_traits_1ERS_0[,2],
                             theta_bin = cut(MNRM_true_traits_1ERS_0[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(MNRM_true_traits_1ERS_0[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
MNRM_Tree_df <- data.frame(Tree_theta_bias = MNRM_Tree_est_1ERS_0[,1], Tree_ERS_bias = MNRM_Tree_est_1ERS_0[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])
MNRM_GPCM_df <- data.frame(GPCM_theta_bias = MNRM_GPCM_est_1ERS_0[,1], theta_bin = MNRM_true_bins[,3],
                           ERS_bin = MNRM_true_bins[,4])
MNRM_MNRM_df <- data.frame(MNRM_theta_bias = MNRM_MNRM_est_1ERS_0[,1], MNRM_ERS_bias = MNRM_MNRM_est_1ERS_0[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])



par(old.par)

#Make 1-dimensional graphs based on just the theta and just the ERS bin
MNRM_MNRM_df_plot <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
MNRM_Tree_df_plot <- MNRM_Tree_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
MNRM_GPCM_df_plot <- MNRM_GPCM_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
bins <- names(table(MNRM_true_bins[,3]))

MNRM_theta_bin <- Tree_theta_bin <- GPCM_theta_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_theta_bin[,i] <- 
    mean(MNRM_MNRM_df_plot[MNRM_MNRM_df_plot[,3] == bins[i], 1])
  Tree_theta_bin[,i] <- 
    mean(MNRM_Tree_df_plot[MNRM_Tree_df_plot[,3] == bins[i], 1])
  GPCM_theta_bin[,i] <- 
    mean(MNRM_GPCM_df_plot[MNRM_GPCM_df_plot[,2] == bins[i], 1])
  
  if(Tree_theta_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_theta_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_theta_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_theta_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_theta_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_theta_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Make 1-dimensional graphs based on just the theta and just the ERS bin
MNRM_MNRM_df_plot <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_Tree_df_plot <- MNRM_Tree_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_GPCM_df_plot <- MNRM_GPCM_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
bins <- names(table(MNRM_true_bins[,4]))

MNRM_ERS_bin <- Tree_ERS_bin <- GPCM_ERS_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_ERS_bin[,i] <- 
    mean(MNRM_MNRM_df_plot[MNRM_MNRM_df_plot[,4] == bins[i], 1])
  Tree_ERS_bin[,i] <- 
    mean(MNRM_Tree_df_plot[MNRM_Tree_df_plot[,4] == bins[i], 1])
  GPCM_ERS_bin[,i] <- 
    mean(MNRM_GPCM_df_plot[MNRM_GPCM_df_plot[,3] == bins[i], 1])
  
  if(Tree_ERS_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_ERS_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_ERS_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_ERS_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_ERS_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_ERS_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
MNRM_Tree_df <- MNRM_Tree_df[is.na(MNRM_Tree_df[,3]) == FALSE & 
                               is.na(MNRM_Tree_df[,4]) == FALSE,]
MNRM_MNRM_df <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,3]) == FALSE & 
                               is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_GPCM_df <- MNRM_GPCM_df[is.na(MNRM_GPCM_df[,2]) == FALSE & 
                               is.na(MNRM_GPCM_df[,3]) == FALSE,]

bins <- names(table(MNRM_true_bins[,3]))

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
MNRM_binned_tree_1ERS_0 <- MNRM_binned_MNRM_1ERS_0 <- MNRM_binned_GPCM_1ERS_0 <- 
  data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(MNRM_binned_MNRM_1ERS_0) <- rownames(MNRM_binned_GPCM_1ERS_0) <- 
  rownames(MNRM_binned_tree_1ERS_0) <- bins
colnames(MNRM_binned_MNRM_1ERS_0) <- colnames(MNRM_binned_GPCM_1ERS_0) <- 
  colnames(MNRM_binned_tree_1ERS_0) <- bins

for(i in 1:20){
  for(j in 1:20){
    MNRM_binned_MNRM_1ERS_0[i,j] <- mean(MNRM_MNRM_df[MNRM_MNRM_df[,3] == bins[j] &
                                                        MNRM_MNRM_df[,4] == bins[i], 1])
    MNRM_binned_tree_1ERS_0[i,j] <- mean(MNRM_Tree_df[MNRM_Tree_df[,3] == bins[j] &
                                                        MNRM_Tree_df[,4] == bins[i], 1])
    MNRM_binned_GPCM_1ERS_0[i,j] <- mean(MNRM_GPCM_df[MNRM_GPCM_df[,2] == bins[j] &
                                                        MNRM_GPCM_df[,3] == bins[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(MNRM_binned_GPCM_1ERS_0), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_tree_1ERS_0), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_MNRM_1ERS_0), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)















##########################################################################
#Bins based on MNRM data, threshold -1, 0, 1 with +1 ERS mean for group 2
#Intialize matrices
MNRM_true_traits_1ERS_1 <- MNRM_Tree_est_1ERS_1 <- MNRM_MNRM_est_1ERS_1 <- 
  matrix(NA, nrow = 250000, ncol = 2)
MNRM_GPCM_est_1ERS_1 <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  MNRM_true_traits_1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_1[[i]][["True_traits"]][501:1000,]
  MNRM_MNRM_est_1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_1[[i]][["MNRM_traits"]][501:1000,]
  MNRM_Tree_est_1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_1[[i]][["Tree_traits"]][501:1000,]
  MNRM_GPCM_est_1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
MNRM_MNRM_est_1ERS_1 <- MNRM_MNRM_est_1ERS_1 - MNRM_true_traits_1ERS_1
MNRM_Tree_est_1ERS_1 <- MNRM_Tree_est_1ERS_1 - MNRM_true_traits_1ERS_1
MNRM_GPCM_est_1ERS_1 <- MNRM_GPCM_est_1ERS_1 - MNRM_true_traits_1ERS_1[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
MNRM_true_bins <- data.frame(True_theta = MNRM_true_traits_1ERS_1[,1], True_ERS = MNRM_true_traits_1ERS_1[,2],
                             theta_bin = cut(MNRM_true_traits_1ERS_1[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(MNRM_true_traits_1ERS_1[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
MNRM_Tree_df <- data.frame(Tree_theta_bias = MNRM_Tree_est_1ERS_1[,1], Tree_ERS_bias = MNRM_Tree_est_1ERS_1[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])
MNRM_GPCM_df <- data.frame(GPCM_theta_bias = MNRM_GPCM_est_1ERS_1[,1], theta_bin = MNRM_true_bins[,3],
                           ERS_bin = MNRM_true_bins[,4])
MNRM_MNRM_df <- data.frame(MNRM_theta_bias = MNRM_MNRM_est_1ERS_1[,1], MNRM_ERS_bias = MNRM_MNRM_est_1ERS_1[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])



par(old.par)

#Make 1-dimensional graphs based on just the theta and just the ERS bin
MNRM_MNRM_df_plot <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
MNRM_Tree_df_plot <- MNRM_Tree_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
MNRM_GPCM_df_plot <- MNRM_GPCM_df[is.na(MNRM_MNRM_df[,3]) == FALSE,]
bins <- names(table(MNRM_true_bins[,3]))

MNRM_theta_bin <- Tree_theta_bin <- GPCM_theta_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_theta_bin[,i] <- 
    mean(MNRM_MNRM_df_plot[MNRM_MNRM_df_plot[,3] == bins[i], 1])
  Tree_theta_bin[,i] <- 
    mean(MNRM_Tree_df_plot[MNRM_Tree_df_plot[,3] == bins[i], 1])
  GPCM_theta_bin[,i] <- 
    mean(MNRM_GPCM_df_plot[MNRM_GPCM_df_plot[,2] == bins[i], 1])
  
  if(Tree_theta_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_theta_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_theta_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_theta_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_theta_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_theta_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Make 1-dimensional graphs based on just the theta and just the ERS bin
MNRM_MNRM_df_plot <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_Tree_df_plot <- MNRM_Tree_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_GPCM_df_plot <- MNRM_GPCM_df[is.na(MNRM_MNRM_df[,4]) == FALSE,]
bins <- names(table(MNRM_true_bins[,4]))

MNRM_ERS_bin <- Tree_ERS_bin <- GPCM_ERS_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_ERS_bin[,i] <- 
    mean(MNRM_MNRM_df_plot[MNRM_MNRM_df_plot[,4] == bins[i], 1])
  Tree_ERS_bin[,i] <- 
    mean(MNRM_Tree_df_plot[MNRM_Tree_df_plot[,4] == bins[i], 1])
  GPCM_ERS_bin[,i] <- 
    mean(MNRM_GPCM_df_plot[MNRM_GPCM_df_plot[,3] == bins[i], 1])
  
  if(Tree_ERS_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_ERS_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_ERS_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_ERS_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_ERS_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_ERS_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
MNRM_Tree_df <- MNRM_Tree_df[is.na(MNRM_Tree_df[,3]) == FALSE & 
                               is.na(MNRM_Tree_df[,4]) == FALSE,]
MNRM_MNRM_df <- MNRM_MNRM_df[is.na(MNRM_MNRM_df[,3]) == FALSE & 
                               is.na(MNRM_MNRM_df[,4]) == FALSE,]
MNRM_GPCM_df <- MNRM_GPCM_df[is.na(MNRM_GPCM_df[,2]) == FALSE & 
                               is.na(MNRM_GPCM_df[,3]) == FALSE,]

bins <- names(table(MNRM_true_bins[,3]))

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
MNRM_binned_tree_1ERS_1 <- MNRM_binned_MNRM_1ERS_1 <- MNRM_binned_GPCM_1ERS_1 <- 
  data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(MNRM_binned_MNRM_1ERS_1) <- rownames(MNRM_binned_GPCM_1ERS_1) <- 
  rownames(MNRM_binned_tree_1ERS_1) <- bins
colnames(MNRM_binned_MNRM_1ERS_1) <- colnames(MNRM_binned_GPCM_1ERS_1) <- 
  colnames(MNRM_binned_tree_1ERS_1) <- bins

for(i in 1:20){
  for(j in 1:20){
    MNRM_binned_MNRM_1ERS_1[i,j] <- mean(MNRM_MNRM_df[MNRM_MNRM_df[,3] == bins[j] &
                                                        MNRM_MNRM_df[,4] == bins[i], 1])
    MNRM_binned_tree_1ERS_1[i,j] <- mean(MNRM_Tree_df[MNRM_Tree_df[,3] == bins[j] &
                                                        MNRM_Tree_df[,4] == bins[i], 1])
    MNRM_binned_GPCM_1ERS_1[i,j] <- mean(MNRM_GPCM_df[MNRM_GPCM_df[,2] == bins[j] &
                                                        MNRM_GPCM_df[,3] == bins[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(MNRM_binned_GPCM_1ERS_1), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_tree_1ERS_1), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_MNRM_1ERS_1), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)






##################################################################
#IRTree as generating model
#Only group 2 included
#-1 ERS, no shift
#Intialize matrices
Tree_true_traits <- Tree_Tree_est <- Tree_MNRM_est <- 
  matrix(NA, nrow = 250000, ncol = 2)
Tree_GPCM_est <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  Tree_true_traits[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_min1[[i]][["True_traits"]][501:1000,]
  Tree_MNRM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_min1[[i]][["MNRM_traits"]][501:1000,]
  Tree_Tree_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_min1[[i]][["Tree_traits"]][501:1000,]
  Tree_GPCM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_min1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
Tree_MNRM_est <- Tree_MNRM_est - Tree_true_traits
Tree_Tree_est <- Tree_Tree_est - Tree_true_traits
Tree_GPCM_est <- Tree_GPCM_est - Tree_true_traits[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
Tree_true_bins <- data.frame(True_theta = Tree_true_traits[,1], True_ERS = Tree_true_traits[,2],
                             theta_bin = cut(Tree_true_traits[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(Tree_true_traits[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
Tree_Tree_df <- data.frame(Tree_theta_bias = Tree_Tree_est[,1], Tree_ERS_bias = Tree_Tree_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])
Tree_GPCM_df <- data.frame(GPCM_theta_bias = Tree_GPCM_est[,1], theta_bin = Tree_true_bins[,3],
                           ERS_bin = Tree_true_bins[,4])
Tree_MNRM_df <- data.frame(MNRM_theta_bias = Tree_MNRM_est[,1], MNRM_ERS_bias = Tree_MNRM_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])



par(old.par)

#Make 1-dimensional graphs based on just the theta and just the ERS bin
Tree_MNRM_df_plot <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
Tree_Tree_df_plot <- Tree_Tree_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
Tree_GPCM_df_plot <- Tree_GPCM_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
bins <- names(table(Tree_true_bins[,3]))

MNRM_theta_bin <- Tree_theta_bin <- GPCM_theta_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_theta_bin[,i] <- 
    mean(Tree_MNRM_df_plot[Tree_MNRM_df_plot[,3] == bins[i], 1])
  Tree_theta_bin[,i] <- 
    mean(Tree_Tree_df_plot[Tree_Tree_df_plot[,3] == bins[i], 1])
  GPCM_theta_bin[,i] <- 
    mean(Tree_GPCM_df_plot[Tree_GPCM_df_plot[,2] == bins[i], 1])
  
  if(Tree_theta_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_theta_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_theta_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_theta_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_theta_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_theta_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Make 1-dimensional graphs based on just the theta and just the ERS bin
Tree_MNRM_df_plot <- Tree_MNRM_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_Tree_df_plot <- Tree_Tree_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df_plot <- Tree_GPCM_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
bins <- names(table(MNRM_true_bins[,4]))

MNRM_ERS_bin <- Tree_ERS_bin <- GPCM_ERS_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_ERS_bin[,i] <- 
    mean(Tree_MNRM_df_plot[Tree_MNRM_df_plot[,4] == bins[i], 1])
  Tree_ERS_bin[,i] <- 
    mean(Tree_Tree_df_plot[Tree_Tree_df_plot[,4] == bins[i], 1])
  GPCM_ERS_bin[,i] <- 
    mean(Tree_GPCM_df_plot[Tree_GPCM_df_plot[,3] == bins[i], 1])
  
  if(Tree_ERS_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_ERS_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_ERS_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_ERS_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_ERS_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_ERS_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
Tree_Tree_df <- Tree_Tree_df[is.na(Tree_Tree_df[,3]) == FALSE & is.na(Tree_Tree_df[,4]) == FALSE,]
Tree_MNRM_df <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE & is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df <- Tree_GPCM_df[is.na(Tree_GPCM_df[,2]) == FALSE & is.na(Tree_GPCM_df[,3]) == FALSE,]

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
Tree_binned_tree <- Tree_binned_MNRM <- Tree_binned_GPCM <- data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(Tree_binned_MNRM) <- rownames(Tree_binned_GPCM) <-
  rownames(Tree_binned_tree) <- 
  names(table(Tree_true_bins[,4]))
colnames(Tree_binned_MNRM) <- colnames(Tree_binned_GPCM) <- 
  colnames(Tree_binned_tree) <- names(table(Tree_true_bins[,3]))

for(i in 1:20){
  for(j in 1:20){
    Tree_binned_MNRM[i,j] <- mean(Tree_MNRM_df[Tree_MNRM_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_MNRM_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_tree[i,j] <- mean(Tree_Tree_df[Tree_Tree_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_Tree_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_GPCM[i,j] <- mean(Tree_GPCM_df[Tree_GPCM_df[,2] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_GPCM_df[,3] == names(table(Tree_true_bins[,4]))[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(Tree_binned_GPCM), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_tree), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_MNRM), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)








#Only group 2 included
#-1 ERS, with shift


#Intialize matrices
Tree_true_traits <- Tree_Tree_est <- Tree_MNRM_est <- 
  matrix(NA, nrow = 250000, ncol = 2)
Tree_GPCM_est <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  Tree_true_traits[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_min1[[i]][["True_traits"]][501:1000,]
  Tree_MNRM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_min1[[i]][["MNRM_traits"]][501:1000,]
  Tree_Tree_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_min1[[i]][["Tree_traits"]][501:1000,]
  Tree_GPCM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_min1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
Tree_MNRM_est <- Tree_MNRM_est - Tree_true_traits
Tree_Tree_est <- Tree_Tree_est - Tree_true_traits
Tree_GPCM_est <- Tree_GPCM_est - Tree_true_traits[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
Tree_true_bins <- data.frame(True_theta = Tree_true_traits[,1], True_ERS = Tree_true_traits[,2],
                             theta_bin = cut(Tree_true_traits[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(Tree_true_traits[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
Tree_Tree_df <- data.frame(Tree_theta_bias = Tree_Tree_est[,1], Tree_ERS_bias = Tree_Tree_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])
Tree_GPCM_df <- data.frame(GPCM_theta_bias = Tree_GPCM_est[,1], theta_bin = Tree_true_bins[,3],
                           ERS_bin = Tree_true_bins[,4])
Tree_MNRM_df <- data.frame(MNRM_theta_bias = Tree_MNRM_est[,1], MNRM_ERS_bias = Tree_MNRM_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])



par(old.par)

#Make 1-dimensional graphs based on just the theta and just the ERS bin
Tree_MNRM_df_plot <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
Tree_Tree_df_plot <- Tree_Tree_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
Tree_GPCM_df_plot <- Tree_GPCM_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
bins <- names(table(Tree_true_bins[,3]))

MNRM_theta_bin <- Tree_theta_bin <- GPCM_theta_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_theta_bin[,i] <- 
    mean(Tree_MNRM_df_plot[Tree_MNRM_df_plot[,3] == bins[i], 1])
  Tree_theta_bin[,i] <- 
    mean(Tree_Tree_df_plot[Tree_Tree_df_plot[,3] == bins[i], 1])
  GPCM_theta_bin[,i] <- 
    mean(Tree_GPCM_df_plot[Tree_GPCM_df_plot[,2] == bins[i], 1])
  
  if(Tree_theta_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_theta_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_theta_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_theta_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_theta_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_theta_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Make 1-dimensional graphs based on just the theta and just the ERS bin
Tree_MNRM_df_plot <- Tree_MNRM_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_Tree_df_plot <- Tree_Tree_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df_plot <- Tree_GPCM_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
bins <- names(table(MNRM_true_bins[,4]))

MNRM_ERS_bin <- Tree_ERS_bin <- GPCM_ERS_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_ERS_bin[,i] <- 
    mean(Tree_MNRM_df_plot[Tree_MNRM_df_plot[,4] == bins[i], 1])
  Tree_ERS_bin[,i] <- 
    mean(Tree_Tree_df_plot[Tree_Tree_df_plot[,4] == bins[i], 1])
  GPCM_ERS_bin[,i] <- 
    mean(Tree_GPCM_df_plot[Tree_GPCM_df_plot[,3] == bins[i], 1])
  
  if(Tree_ERS_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_ERS_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_ERS_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_ERS_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_ERS_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_ERS_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))




#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
Tree_Tree_df <- Tree_Tree_df[is.na(Tree_Tree_df[,3]) == FALSE & is.na(Tree_Tree_df[,4]) == FALSE,]
Tree_MNRM_df <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE & is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df <- Tree_GPCM_df[is.na(Tree_GPCM_df[,2]) == FALSE & is.na(Tree_GPCM_df[,3]) == FALSE,]

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
Tree_binned_tree <- Tree_binned_MNRM <- Tree_binned_GPCM <- data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(Tree_binned_MNRM) <- rownames(Tree_binned_GPCM) <-
  rownames(Tree_binned_tree) <- 
  names(table(Tree_true_bins[,4]))
colnames(Tree_binned_MNRM) <- colnames(Tree_binned_GPCM) <- 
  colnames(Tree_binned_tree) <- names(table(Tree_true_bins[,3]))

for(i in 1:20){
  for(j in 1:20){
    Tree_binned_MNRM[i,j] <- mean(Tree_MNRM_df[Tree_MNRM_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_MNRM_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_tree[i,j] <- mean(Tree_Tree_df[Tree_Tree_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_Tree_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_GPCM[i,j] <- mean(Tree_GPCM_df[Tree_GPCM_df[,2] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_GPCM_df[,3] == names(table(Tree_true_bins[,4]))[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(Tree_binned_GPCM), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_tree), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_MNRM), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)









#IRTree as generating model
#Only group 2 included
#0 ERS, no shift
#Intialize matrices
Tree_true_traits <- Tree_Tree_est <- Tree_MNRM_est <- 
  matrix(NA, nrow = 250000, ncol = 2)
Tree_GPCM_est <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  Tree_true_traits[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20[[i]][["True_traits"]][501:1000,]
  Tree_MNRM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20[[i]][["MNRM_traits"]][501:1000,]
  Tree_Tree_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20[[i]][["Tree_traits"]][501:1000,]
  Tree_GPCM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
Tree_MNRM_est <- Tree_MNRM_est - Tree_true_traits
Tree_Tree_est <- Tree_Tree_est - Tree_true_traits
Tree_GPCM_est <- Tree_GPCM_est - Tree_true_traits[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
Tree_true_bins <- data.frame(True_theta = Tree_true_traits[,1], True_ERS = Tree_true_traits[,2],
                             theta_bin = cut(Tree_true_traits[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(Tree_true_traits[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
Tree_Tree_df <- data.frame(Tree_theta_bias = Tree_Tree_est[,1], Tree_ERS_bias = Tree_Tree_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])
Tree_GPCM_df <- data.frame(GPCM_theta_bias = Tree_GPCM_est[,1], theta_bin = Tree_true_bins[,3],
                           ERS_bin = Tree_true_bins[,4])
Tree_MNRM_df <- data.frame(MNRM_theta_bias = Tree_MNRM_est[,1], MNRM_ERS_bias = Tree_MNRM_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])



par(old.par)

#Make 1-dimensional graphs based on just the theta and just the ERS bin
Tree_MNRM_df_plot <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
Tree_Tree_df_plot <- Tree_Tree_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
Tree_GPCM_df_plot <- Tree_GPCM_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
bins <- names(table(Tree_true_bins[,3]))

MNRM_theta_bin <- Tree_theta_bin <- GPCM_theta_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_theta_bin[,i] <- 
    mean(Tree_MNRM_df_plot[Tree_MNRM_df_plot[,3] == bins[i], 1])
  Tree_theta_bin[,i] <- 
    mean(Tree_Tree_df_plot[Tree_Tree_df_plot[,3] == bins[i], 1])
  GPCM_theta_bin[,i] <- 
    mean(Tree_GPCM_df_plot[Tree_GPCM_df_plot[,2] == bins[i], 1])
  
  if(Tree_theta_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_theta_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_theta_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_theta_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_theta_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_theta_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Make 1-dimensional graphs based on just the theta and just the ERS bin
Tree_MNRM_df_plot <- Tree_MNRM_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_Tree_df_plot <- Tree_Tree_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df_plot <- Tree_GPCM_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
bins <- names(table(MNRM_true_bins[,4]))

MNRM_ERS_bin <- Tree_ERS_bin <- GPCM_ERS_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_ERS_bin[,i] <- 
    mean(Tree_MNRM_df_plot[Tree_MNRM_df_plot[,4] == bins[i], 1])
  Tree_ERS_bin[,i] <- 
    mean(Tree_Tree_df_plot[Tree_Tree_df_plot[,4] == bins[i], 1])
  GPCM_ERS_bin[,i] <- 
    mean(Tree_GPCM_df_plot[Tree_GPCM_df_plot[,3] == bins[i], 1])
  
  if(Tree_ERS_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_ERS_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_ERS_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_ERS_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_ERS_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_ERS_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))




#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
Tree_Tree_df <- Tree_Tree_df[is.na(Tree_Tree_df[,3]) == FALSE & is.na(Tree_Tree_df[,4]) == FALSE,]
Tree_MNRM_df <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE & is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df <- Tree_GPCM_df[is.na(Tree_GPCM_df[,2]) == FALSE & is.na(Tree_GPCM_df[,3]) == FALSE,]

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
Tree_binned_tree <- Tree_binned_MNRM <- Tree_binned_GPCM <- data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(Tree_binned_MNRM) <- rownames(Tree_binned_GPCM) <-
  rownames(Tree_binned_tree) <- 
  names(table(Tree_true_bins[,4]))
colnames(Tree_binned_MNRM) <- colnames(Tree_binned_GPCM) <- 
  colnames(Tree_binned_tree) <- names(table(Tree_true_bins[,3]))

for(i in 1:20){
  for(j in 1:20){
    Tree_binned_MNRM[i,j] <- mean(Tree_MNRM_df[Tree_MNRM_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_MNRM_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_tree[i,j] <- mean(Tree_Tree_df[Tree_Tree_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_Tree_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_GPCM[i,j] <- mean(Tree_GPCM_df[Tree_GPCM_df[,2] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_GPCM_df[,3] == names(table(Tree_true_bins[,4]))[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(Tree_binned_GPCM), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_tree), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_MNRM), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)








#Only group 2 included
#-1 ERS, with shift


#Intialize matrices
Tree_true_traits <- Tree_Tree_est <- Tree_MNRM_est <- 
  matrix(NA, nrow = 250000, ncol = 2)
Tree_GPCM_est <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  Tree_true_traits[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20[[i]][["True_traits"]][501:1000,]
  Tree_MNRM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20[[i]][["MNRM_traits"]][501:1000,]
  Tree_Tree_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20[[i]][["Tree_traits"]][501:1000,]
  Tree_GPCM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
Tree_MNRM_est <- Tree_MNRM_est - Tree_true_traits
Tree_Tree_est <- Tree_Tree_est - Tree_true_traits
Tree_GPCM_est <- Tree_GPCM_est - Tree_true_traits[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
Tree_true_bins <- data.frame(True_theta = Tree_true_traits[,1], True_ERS = Tree_true_traits[,2],
                             theta_bin = cut(Tree_true_traits[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(Tree_true_traits[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
Tree_Tree_df <- data.frame(Tree_theta_bias = Tree_Tree_est[,1], Tree_ERS_bias = Tree_Tree_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])
Tree_GPCM_df <- data.frame(GPCM_theta_bias = Tree_GPCM_est[,1], theta_bin = Tree_true_bins[,3],
                           ERS_bin = Tree_true_bins[,4])
Tree_MNRM_df <- data.frame(MNRM_theta_bias = Tree_MNRM_est[,1], MNRM_ERS_bias = Tree_MNRM_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])


par(old.par)

#Make 1-dimensional graphs based on just the theta and just the ERS bin
Tree_MNRM_df_plot <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
Tree_Tree_df_plot <- Tree_Tree_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
Tree_GPCM_df_plot <- Tree_GPCM_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
bins <- names(table(Tree_true_bins[,3]))

MNRM_theta_bin <- Tree_theta_bin <- GPCM_theta_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_theta_bin[,i] <- 
    mean(Tree_MNRM_df_plot[Tree_MNRM_df_plot[,3] == bins[i], 1])
  Tree_theta_bin[,i] <- 
    mean(Tree_Tree_df_plot[Tree_Tree_df_plot[,3] == bins[i], 1])
  GPCM_theta_bin[,i] <- 
    mean(Tree_GPCM_df_plot[Tree_GPCM_df_plot[,2] == bins[i], 1])
  
  if(Tree_theta_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_theta_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_theta_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_theta_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_theta_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_theta_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Make 1-dimensional graphs based on just the theta and just the ERS bin
Tree_MNRM_df_plot <- Tree_MNRM_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_Tree_df_plot <- Tree_Tree_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df_plot <- Tree_GPCM_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
bins <- names(table(MNRM_true_bins[,4]))

MNRM_ERS_bin <- Tree_ERS_bin <- GPCM_ERS_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_ERS_bin[,i] <- 
    mean(Tree_MNRM_df_plot[Tree_MNRM_df_plot[,4] == bins[i], 1])
  Tree_ERS_bin[,i] <- 
    mean(Tree_Tree_df_plot[Tree_Tree_df_plot[,4] == bins[i], 1])
  GPCM_ERS_bin[,i] <- 
    mean(Tree_GPCM_df_plot[Tree_GPCM_df_plot[,3] == bins[i], 1])
  
  if(Tree_ERS_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_ERS_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_ERS_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_ERS_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_ERS_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_ERS_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))




#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
Tree_Tree_df <- Tree_Tree_df[is.na(Tree_Tree_df[,3]) == FALSE & is.na(Tree_Tree_df[,4]) == FALSE,]
Tree_MNRM_df <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE & is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df <- Tree_GPCM_df[is.na(Tree_GPCM_df[,2]) == FALSE & is.na(Tree_GPCM_df[,3]) == FALSE,]

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
Tree_binned_tree <- Tree_binned_MNRM <- Tree_binned_GPCM <- data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(Tree_binned_MNRM) <- rownames(Tree_binned_GPCM) <-
  rownames(Tree_binned_tree) <- 
  names(table(Tree_true_bins[,4]))
colnames(Tree_binned_MNRM) <- colnames(Tree_binned_GPCM) <- 
  colnames(Tree_binned_tree) <- names(table(Tree_true_bins[,3]))

for(i in 1:20){
  for(j in 1:20){
    Tree_binned_MNRM[i,j] <- mean(Tree_MNRM_df[Tree_MNRM_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_MNRM_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_tree[i,j] <- mean(Tree_Tree_df[Tree_Tree_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_Tree_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_GPCM[i,j] <- mean(Tree_GPCM_df[Tree_GPCM_df[,2] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_GPCM_df[,3] == names(table(Tree_true_bins[,4]))[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(Tree_binned_GPCM), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_tree), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_MNRM), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)












#Only group 2 included
#+1 ERS, no shift
#Intialize matrices
Tree_true_traits <- Tree_Tree_est <- Tree_MNRM_est <- 
  matrix(NA, nrow = 250000, ncol = 2)
Tree_GPCM_est <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  Tree_true_traits[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_1[[i]][["True_traits"]][501:1000,]
  Tree_MNRM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_1[[i]][["MNRM_traits"]][501:1000,]
  Tree_Tree_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_1[[i]][["Tree_traits"]][501:1000,]
  Tree_GPCM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
Tree_MNRM_est <- Tree_MNRM_est - Tree_true_traits
Tree_Tree_est <- Tree_Tree_est - Tree_true_traits
Tree_GPCM_est <- Tree_GPCM_est - Tree_true_traits[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
Tree_true_bins <- data.frame(True_theta = Tree_true_traits[,1], True_ERS = Tree_true_traits[,2],
                             theta_bin = cut(Tree_true_traits[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(Tree_true_traits[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
Tree_Tree_df <- data.frame(Tree_theta_bias = Tree_Tree_est[,1], Tree_ERS_bias = Tree_Tree_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])
Tree_GPCM_df <- data.frame(GPCM_theta_bias = Tree_GPCM_est[,1], theta_bin = Tree_true_bins[,3],
                           ERS_bin = Tree_true_bins[,4])
Tree_MNRM_df <- data.frame(MNRM_theta_bias = Tree_MNRM_est[,1], MNRM_ERS_bias = Tree_MNRM_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])



par(old.par)

#Make 1-dimensional graphs based on just the theta and just the ERS bin
Tree_MNRM_df_plot <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
Tree_Tree_df_plot <- Tree_Tree_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
Tree_GPCM_df_plot <- Tree_GPCM_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
bins <- names(table(Tree_true_bins[,3]))

MNRM_theta_bin <- Tree_theta_bin <- GPCM_theta_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_theta_bin[,i] <- 
    mean(Tree_MNRM_df_plot[Tree_MNRM_df_plot[,3] == bins[i], 1])
  Tree_theta_bin[,i] <- 
    mean(Tree_Tree_df_plot[Tree_Tree_df_plot[,3] == bins[i], 1])
  GPCM_theta_bin[,i] <- 
    mean(Tree_GPCM_df_plot[Tree_GPCM_df_plot[,2] == bins[i], 1])
  
  if(Tree_theta_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_theta_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_theta_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_theta_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_theta_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_theta_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Make 1-dimensional graphs based on just the theta and just the ERS bin
Tree_MNRM_df_plot <- Tree_MNRM_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_Tree_df_plot <- Tree_Tree_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df_plot <- Tree_GPCM_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
bins <- names(table(MNRM_true_bins[,4]))

MNRM_ERS_bin <- Tree_ERS_bin <- GPCM_ERS_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_ERS_bin[,i] <- 
    mean(Tree_MNRM_df_plot[Tree_MNRM_df_plot[,4] == bins[i], 1])
  Tree_ERS_bin[,i] <- 
    mean(Tree_Tree_df_plot[Tree_Tree_df_plot[,4] == bins[i], 1])
  GPCM_ERS_bin[,i] <- 
    mean(Tree_GPCM_df_plot[Tree_GPCM_df_plot[,3] == bins[i], 1])
  
  if(Tree_ERS_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_ERS_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_ERS_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_ERS_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_ERS_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_ERS_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
Tree_Tree_df <- Tree_Tree_df[is.na(Tree_Tree_df[,3]) == FALSE & is.na(Tree_Tree_df[,4]) == FALSE,]
Tree_MNRM_df <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE & is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df <- Tree_GPCM_df[is.na(Tree_GPCM_df[,2]) == FALSE & is.na(Tree_GPCM_df[,3]) == FALSE,]



#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
Tree_binned_tree <- Tree_binned_MNRM <- Tree_binned_GPCM <- data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(Tree_binned_MNRM) <- rownames(Tree_binned_GPCM) <-
  rownames(Tree_binned_tree) <- 
  names(table(Tree_true_bins[,4]))
colnames(Tree_binned_MNRM) <- colnames(Tree_binned_GPCM) <- 
  colnames(Tree_binned_tree) <- names(table(Tree_true_bins[,3]))

for(i in 1:20){
  for(j in 1:20){
    Tree_binned_MNRM[i,j] <- mean(Tree_MNRM_df[Tree_MNRM_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_MNRM_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_tree[i,j] <- mean(Tree_Tree_df[Tree_Tree_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_Tree_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_GPCM[i,j] <- mean(Tree_GPCM_df[Tree_GPCM_df[,2] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_GPCM_df[,3] == names(table(Tree_true_bins[,4]))[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(Tree_binned_GPCM), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_tree), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_MNRM), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)








#Only group 2 included
#+1 ERS, with shift


#Intialize matrices
Tree_true_traits <- Tree_Tree_est <- Tree_MNRM_est <- 
  matrix(NA, nrow = 250000, ncol = 2)
Tree_GPCM_est <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  Tree_true_traits[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_1[[i]][["True_traits"]][501:1000,]
  Tree_MNRM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_1[[i]][["MNRM_traits"]][501:1000,]
  Tree_Tree_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_1[[i]][["Tree_traits"]][501:1000,]
  Tree_GPCM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
Tree_MNRM_est <- Tree_MNRM_est - Tree_true_traits
Tree_Tree_est <- Tree_Tree_est - Tree_true_traits
Tree_GPCM_est <- Tree_GPCM_est - Tree_true_traits[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
Tree_true_bins <- data.frame(True_theta = Tree_true_traits[,1], True_ERS = Tree_true_traits[,2],
                             theta_bin = cut(Tree_true_traits[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(Tree_true_traits[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
Tree_Tree_df <- data.frame(Tree_theta_bias = Tree_Tree_est[,1], Tree_ERS_bias = Tree_Tree_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])
Tree_GPCM_df <- data.frame(GPCM_theta_bias = Tree_GPCM_est[,1], theta_bin = Tree_true_bins[,3],
                           ERS_bin = Tree_true_bins[,4])
Tree_MNRM_df <- data.frame(MNRM_theta_bias = Tree_MNRM_est[,1], MNRM_ERS_bias = Tree_MNRM_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])


par(old.par)

#Make 1-dimensional graphs based on just the theta and just the ERS bin
Tree_MNRM_df_plot <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
Tree_Tree_df_plot <- Tree_Tree_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
Tree_GPCM_df_plot <- Tree_GPCM_df[is.na(Tree_MNRM_df[,3]) == FALSE,]
bins <- names(table(Tree_true_bins[,3]))

MNRM_theta_bin <- Tree_theta_bin <- GPCM_theta_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_theta_bin[,i] <- 
    mean(Tree_MNRM_df_plot[Tree_MNRM_df_plot[,3] == bins[i], 1])
  Tree_theta_bin[,i] <- 
    mean(Tree_Tree_df_plot[Tree_Tree_df_plot[,3] == bins[i], 1])
  GPCM_theta_bin[,i] <- 
    mean(Tree_GPCM_df_plot[Tree_GPCM_df_plot[,2] == bins[i], 1])
  
  if(Tree_theta_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_theta_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_theta_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_theta_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_theta_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_theta_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True substantive trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Make 1-dimensional graphs based on just the theta and just the ERS bin
Tree_MNRM_df_plot <- Tree_MNRM_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_Tree_df_plot <- Tree_Tree_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df_plot <- Tree_GPCM_df[is.na(Tree_MNRM_df[,4]) == FALSE,]
bins <- names(table(MNRM_true_bins[,4]))

MNRM_ERS_bin <- Tree_ERS_bin <- GPCM_ERS_bin <- matrix(NA, nrow = 1, ncol = 20)
MNRM_colours <- Tree_colours <- GPCM_colours <- rep(NA, 20)
for(i in 1:20){
  MNRM_ERS_bin[,i] <- 
    mean(Tree_MNRM_df_plot[Tree_MNRM_df_plot[,4] == bins[i], 1])
  Tree_ERS_bin[,i] <- 
    mean(Tree_Tree_df_plot[Tree_Tree_df_plot[,4] == bins[i], 1])
  GPCM_ERS_bin[,i] <- 
    mean(Tree_GPCM_df_plot[Tree_GPCM_df_plot[,3] == bins[i], 1])
  
  if(Tree_ERS_bin[,i] < 0){
    Tree_colours[i] <- "blue"
  } else{
    Tree_colours[i] <- "red"
  }
  
  if(GPCM_ERS_bin[,i] < 0){
    GPCM_colours[i] <- "blue"
  } else{
    GPCM_colours[i] <- "red"
  }
  
  if(MNRM_ERS_bin[,i] < 0){
    MNRM_colours[i] <- "blue"
  } else{
    MNRM_colours[i] <- "red"
  }
}

plot(as.numeric(GPCM_ERS_bin), type = "h",lwd = 3,
     col = GPCM_colours,
     ylab = y_axis_names[1],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(Tree_ERS_bin), type = "h",lwd = 3,
     col = Tree_colours,
     ylab = y_axis_names[2],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))

plot(as.numeric(MNRM_ERS_bin), type = "h",lwd = 3,
     col = MNRM_colours,
     ylab = y_axis_names[3],
     xlab = "True ERS trait",
     xaxt = "n",
     ylim = y_axis)
axis(1, at = c(1, 10.5, 20), labels = c(-2, 0, 2))



#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
Tree_Tree_df <- Tree_Tree_df[is.na(Tree_Tree_df[,3]) == FALSE & is.na(Tree_Tree_df[,4]) == FALSE,]
Tree_MNRM_df <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE & is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df <- Tree_GPCM_df[is.na(Tree_GPCM_df[,2]) == FALSE & is.na(Tree_GPCM_df[,3]) == FALSE,]

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
Tree_binned_tree <- Tree_binned_MNRM <- Tree_binned_GPCM <- data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(Tree_binned_MNRM) <- rownames(Tree_binned_GPCM) <-
  rownames(Tree_binned_tree) <- 
  names(table(Tree_true_bins[,4]))
colnames(Tree_binned_MNRM) <- colnames(Tree_binned_GPCM) <- 
  colnames(Tree_binned_tree) <- names(table(Tree_true_bins[,3]))

for(i in 1:20){
  for(j in 1:20){
    Tree_binned_MNRM[i,j] <- mean(Tree_MNRM_df[Tree_MNRM_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_MNRM_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_tree[i,j] <- mean(Tree_Tree_df[Tree_Tree_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_Tree_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_GPCM[i,j] <- mean(Tree_GPCM_df[Tree_GPCM_df[,2] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_GPCM_df[,3] == names(table(Tree_true_bins[,4]))[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(Tree_binned_GPCM), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_tree), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_MNRM), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)







































################################################################################
#True model estimated values as baseline instead of true values
#All possible plots on six pages (6 plots per page)
############################################################################

#-1 ERS condition MNRM no shift
#Only group 2 included here!
#Bins based on MNRM data, threshold -1, 0, 1 with +1 ERS mean for group 2
MNRM_true_traits_min1ERS_0 <- MNRM_Tree_est_min1ERS_0 <- MNRM_MNRM_est_min1ERS_0 <- 
  matrix(NA, nrow = 250000, ncol = 2)
MNRM_GPCM_est_min1ERS_0 <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  MNRM_true_traits_min1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_min1[[i]][["True_traits"]][501:1000,]
  MNRM_MNRM_est_min1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_min1[[i]][["MNRM_traits"]][501:1000,]
  MNRM_Tree_est_min1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_min1[[i]][["Tree_traits"]][501:1000,]
  MNRM_GPCM_est_min1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_min1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
MNRM_Tree_est_min1ERS_0 <- MNRM_Tree_est_min1ERS_0 - MNRM_MNRM_est_min1ERS_0
MNRM_GPCM_est_min1ERS_0 <- MNRM_GPCM_est_min1ERS_0 - MNRM_MNRM_est_min1ERS_0[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
MNRM_true_bins <- data.frame(True_theta = MNRM_true_traits_min1ERS_0[,1], True_ERS = MNRM_true_traits_min1ERS_0[,2],
                             theta_bin = cut(MNRM_true_traits_min1ERS_0[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(MNRM_true_traits_min1ERS_0[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
MNRM_Tree_df <- data.frame(Tree_theta_bias = MNRM_Tree_est_min1ERS_0[,1], Tree_ERS_bias = MNRM_Tree_est_min1ERS_0[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])
MNRM_GPCM_df <- data.frame(GPCM_theta_bias = MNRM_GPCM_est_min1ERS_0[,1], theta_bin = MNRM_true_bins[,3],
                           ERS_bin = MNRM_true_bins[,4])

#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
MNRM_Tree_df <- MNRM_Tree_df[is.na(MNRM_Tree_df[,3]) == FALSE & 
                               is.na(MNRM_Tree_df[,4]) == FALSE,]
MNRM_GPCM_df <- MNRM_GPCM_df[is.na(MNRM_GPCM_df[,2]) == FALSE & 
                               is.na(MNRM_GPCM_df[,3]) == FALSE,]

bins <- names(table(MNRM_true_bins[,3]))

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
MNRM_binned_tree_min1ERS_0 <- MNRM_binned_GPCM_min1ERS_0 <- 
  data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(MNRM_binned_GPCM_min1ERS_0) <- rownames(MNRM_binned_tree_min1ERS_0) <-
  bins
colnames(MNRM_binned_GPCM_min1ERS_0) <- colnames(MNRM_binned_tree_min1ERS_0) <- 
  bins

for(i in 1:20){
  for(j in 1:20){
    MNRM_binned_tree_min1ERS_0[i,j] <- mean(MNRM_Tree_df[MNRM_Tree_df[,3] == bins[j] &
                                                           MNRM_Tree_df[,4] == bins[i], 1])
    MNRM_binned_GPCM_min1ERS_0[i,j] <- mean(MNRM_GPCM_df[MNRM_GPCM_df[,2] == bins[j] &
                                                           MNRM_GPCM_df[,3] == bins[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(MNRM_binned_GPCM_min1ERS_0), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_tree_min1ERS_0), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)










#-1 ERS condition MNRM with shift
#Only group 2 included here!
#Bins based on MNRM data, threshold -1, 0, 1 with +1 ERS mean for group 2
#Intialize matrices
MNRM_true_traits_min1ERS_1 <- MNRM_Tree_est_min1ERS_1 <- MNRM_MNRM_est_min1ERS_1 <- 
  matrix(NA, nrow = 250000, ncol = 2)
MNRM_GPCM_est_min1ERS_1 <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  MNRM_true_traits_min1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_min1[[i]][["True_traits"]][501:1000,]
  MNRM_MNRM_est_min1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_min1[[i]][["MNRM_traits"]][501:1000,]
  MNRM_Tree_est_min1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_min1[[i]][["Tree_traits"]][501:1000,]
  MNRM_GPCM_est_min1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_min1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
MNRM_Tree_est_min1ERS_1 <- MNRM_Tree_est_min1ERS_1 - MNRM_MNRM_est_min1ERS_1
MNRM_GPCM_est_min1ERS_1 <- MNRM_GPCM_est_min1ERS_1 - MNRM_MNRM_est_min1ERS_1[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
MNRM_true_bins <- data.frame(True_theta = MNRM_true_traits_min1ERS_1[,1], True_ERS = MNRM_true_traits_min1ERS_1[,2],
                             theta_bin = cut(MNRM_true_traits_min1ERS_1[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(MNRM_true_traits_min1ERS_1[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
MNRM_Tree_df <- data.frame(Tree_theta_bias = MNRM_Tree_est_min1ERS_1[,1], Tree_ERS_bias = MNRM_Tree_est_min1ERS_1[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])
MNRM_GPCM_df <- data.frame(GPCM_theta_bias = MNRM_GPCM_est_min1ERS_1[,1], theta_bin = MNRM_true_bins[,3],
                           ERS_bin = MNRM_true_bins[,4])

#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
MNRM_Tree_df <- MNRM_Tree_df[is.na(MNRM_Tree_df[,3]) == FALSE & 
                               is.na(MNRM_Tree_df[,4]) == FALSE,]
MNRM_GPCM_df <- MNRM_GPCM_df[is.na(MNRM_GPCM_df[,2]) == FALSE & 
                               is.na(MNRM_GPCM_df[,3]) == FALSE,]

bins <- names(table(MNRM_true_bins[,3]))

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
MNRM_binned_tree_min1ERS_1 <- MNRM_binned_GPCM_min1ERS_1 <- 
  data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(MNRM_binned_GPCM_min1ERS_1) <- 
  rownames(MNRM_binned_tree_min1ERS_1) <- bins
colnames(MNRM_binned_GPCM_min1ERS_1) <- 
  colnames(MNRM_binned_tree_min1ERS_1) <- bins

for(i in 1:20){
  for(j in 1:20){
    MNRM_binned_tree_min1ERS_1[i,j] <- mean(MNRM_Tree_df[MNRM_Tree_df[,3] == bins[j] &
                                                           MNRM_Tree_df[,4] == bins[i], 1])
    MNRM_binned_GPCM_min1ERS_1[i,j] <- mean(MNRM_GPCM_df[MNRM_GPCM_df[,2] == bins[j] &
                                                           MNRM_GPCM_df[,3] == bins[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(MNRM_binned_GPCM_min1ERS_1), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_tree_min1ERS_1), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)










#0 ERS condition MNRM with no shift
#Only group 2 included here!
#Bins based on MNRM data, threshold -1, 0, 1 with +1 ERS mean for group 2
MNRM_true_traits_0ERS_0 <- MNRM_Tree_est_0ERS_0 <- MNRM_MNRM_est_0ERS_0 <- 
  matrix(NA, nrow = 250000, ncol = 2)
MNRM_GPCM_est_0ERS_0 <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  MNRM_true_traits_0ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20[[i]][["True_traits"]][501:1000,]
  MNRM_MNRM_est_0ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20[[i]][["MNRM_traits"]][501:1000,]
  MNRM_Tree_est_0ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20[[i]][["Tree_traits"]][501:1000,]
  MNRM_GPCM_est_0ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
MNRM_Tree_est_0ERS_0 <- MNRM_Tree_est_0ERS_0 - MNRM_MNRM_est_0ERS_0
MNRM_GPCM_est_0ERS_0 <- MNRM_GPCM_est_0ERS_0 - MNRM_MNRM_est_0ERS_0[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
MNRM_true_bins <- data.frame(True_theta = MNRM_true_traits_0ERS_0[,1], True_ERS = MNRM_true_traits_0ERS_0[,2],
                             theta_bin = cut(MNRM_true_traits_0ERS_0[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(MNRM_true_traits_0ERS_0[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
MNRM_Tree_df <- data.frame(Tree_theta_bias = MNRM_Tree_est_0ERS_0[,1], Tree_ERS_bias = MNRM_Tree_est_0ERS_0[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])
MNRM_GPCM_df <- data.frame(GPCM_theta_bias = MNRM_GPCM_est_0ERS_0[,1], theta_bin = MNRM_true_bins[,3],
                           ERS_bin = MNRM_true_bins[,4])

#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
MNRM_Tree_df <- MNRM_Tree_df[is.na(MNRM_Tree_df[,3]) == FALSE & 
                               is.na(MNRM_Tree_df[,4]) == FALSE,]
MNRM_GPCM_df <- MNRM_GPCM_df[is.na(MNRM_GPCM_df[,2]) == FALSE & 
                               is.na(MNRM_GPCM_df[,3]) == FALSE,]

bins <- names(table(MNRM_true_bins[,3]))

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
MNRM_binned_tree_0ERS_0 <- MNRM_binned_MNRM_0ERS_0 <- MNRM_binned_GPCM_0ERS_0 <- 
  data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(MNRM_binned_GPCM_0ERS_0) <- 
  rownames(MNRM_binned_tree_0ERS_0) <- bins
colnames(MNRM_binned_GPCM_0ERS_0) <- 
  colnames(MNRM_binned_tree_0ERS_0) <- bins

for(i in 1:20){
  for(j in 1:20){
    MNRM_binned_tree_0ERS_0[i,j] <- mean(MNRM_Tree_df[MNRM_Tree_df[,3] == bins[j] &
                                                        MNRM_Tree_df[,4] == bins[i], 1])
    MNRM_binned_GPCM_0ERS_0[i,j] <- mean(MNRM_GPCM_df[MNRM_GPCM_df[,2] == bins[j] &
                                                        MNRM_GPCM_df[,3] == bins[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(MNRM_binned_GPCM_0ERS_0), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_tree_0ERS_0), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)









#0 ERS condition MNRM with shift
#Only group 2 included here!
#Bins based on MNRM data, threshold -1, 0, 1 with +1 ERS mean for group 2
#Intialize matrices
MNRM_true_traits_0ERS_1 <- MNRM_Tree_est_0ERS_1 <- MNRM_MNRM_est_0ERS_1 <- 
  matrix(NA, nrow = 250000, ncol = 2)
MNRM_GPCM_est_0ERS_1 <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  MNRM_true_traits_0ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20[[i]][["True_traits"]][501:1000,]
  MNRM_MNRM_est_0ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20[[i]][["MNRM_traits"]][501:1000,]
  MNRM_Tree_est_0ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20[[i]][["Tree_traits"]][501:1000,]
  MNRM_GPCM_est_0ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
MNRM_Tree_est_0ERS_1 <- MNRM_Tree_est_0ERS_1 - MNRM_MNRM_est_0ERS_1
MNRM_GPCM_est_0ERS_1 <- MNRM_GPCM_est_0ERS_1 - MNRM_MNRM_est_0ERS_1[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
MNRM_true_bins <- data.frame(True_theta = MNRM_true_traits_0ERS_1[,1], True_ERS = MNRM_true_traits_0ERS_1[,2],
                             theta_bin = cut(MNRM_true_traits_0ERS_1[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(MNRM_true_traits_0ERS_1[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
MNRM_Tree_df <- data.frame(Tree_theta_bias = MNRM_Tree_est_0ERS_1[,1], Tree_ERS_bias = MNRM_Tree_est_0ERS_1[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])
MNRM_GPCM_df <- data.frame(GPCM_theta_bias = MNRM_GPCM_est_0ERS_1[,1], theta_bin = MNRM_true_bins[,3],
                           ERS_bin = MNRM_true_bins[,4])

#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
MNRM_Tree_df <- MNRM_Tree_df[is.na(MNRM_Tree_df[,3]) == FALSE & 
                               is.na(MNRM_Tree_df[,4]) == FALSE,]
MNRM_GPCM_df <- MNRM_GPCM_df[is.na(MNRM_GPCM_df[,2]) == FALSE & 
                               is.na(MNRM_GPCM_df[,3]) == FALSE,]

bins <- names(table(MNRM_true_bins[,3]))

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
MNRM_binned_tree_0ERS_1 <- MNRM_binned_MNRM_0ERS_1 <- MNRM_binned_GPCM_0ERS_1 <- 
  data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(MNRM_binned_GPCM_0ERS_1) <- 
  rownames(MNRM_binned_tree_0ERS_1) <- bins
colnames(MNRM_binned_GPCM_0ERS_1) <- 
  colnames(MNRM_binned_tree_0ERS_1) <- bins

for(i in 1:20){
  for(j in 1:20){
    MNRM_binned_tree_0ERS_1[i,j] <- mean(MNRM_Tree_df[MNRM_Tree_df[,3] == bins[j] &
                                                        MNRM_Tree_df[,4] == bins[i], 1])
    MNRM_binned_GPCM_0ERS_1[i,j] <- mean(MNRM_GPCM_df[MNRM_GPCM_df[,2] == bins[j] &
                                                        MNRM_GPCM_df[,3] == bins[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(MNRM_binned_GPCM_0ERS_1), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_tree_0ERS_1), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)











#########################################################################
#Only group 2 included here!
#Bins based on MNRM data, threshold -1, 0, 1 with +1 ERS mean for group 2
#Intialize matrices
MNRM_true_traits_1ERS_0 <- MNRM_Tree_est_1ERS_0 <- MNRM_MNRM_est_1ERS_0 <- 
  matrix(NA, nrow = 250000, ncol = 2)
MNRM_GPCM_est_1ERS_0 <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  MNRM_true_traits_1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_1[[i]][["True_traits"]][501:1000,]
  MNRM_MNRM_est_1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_1[[i]][["MNRM_traits"]][501:1000,]
  MNRM_Tree_est_1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_1[[i]][["Tree_traits"]][501:1000,]
  MNRM_GPCM_est_1ERS_0[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_average_20_ERS_1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
MNRM_Tree_est_1ERS_0 <- MNRM_Tree_est_1ERS_0 - MNRM_MNRM_est_1ERS_0
MNRM_GPCM_est_1ERS_0 <- MNRM_GPCM_est_1ERS_0 - MNRM_MNRM_est_1ERS_0[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
MNRM_true_bins <- data.frame(True_theta = MNRM_true_traits_1ERS_0[,1], True_ERS = MNRM_true_traits_1ERS_0[,2],
                             theta_bin = cut(MNRM_true_traits_1ERS_0[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(MNRM_true_traits_1ERS_0[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
MNRM_Tree_df <- data.frame(Tree_theta_bias = MNRM_Tree_est_1ERS_0[,1], Tree_ERS_bias = MNRM_Tree_est_1ERS_0[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])
MNRM_GPCM_df <- data.frame(GPCM_theta_bias = MNRM_GPCM_est_1ERS_0[,1], theta_bin = MNRM_true_bins[,3],
                           ERS_bin = MNRM_true_bins[,4])

#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
MNRM_Tree_df <- MNRM_Tree_df[is.na(MNRM_Tree_df[,3]) == FALSE & 
                               is.na(MNRM_Tree_df[,4]) == FALSE,]
MNRM_GPCM_df <- MNRM_GPCM_df[is.na(MNRM_GPCM_df[,2]) == FALSE & 
                               is.na(MNRM_GPCM_df[,3]) == FALSE,]

bins <- names(table(MNRM_true_bins[,3]))

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
MNRM_binned_tree_1ERS_0 <- MNRM_binned_MNRM_1ERS_0 <- MNRM_binned_GPCM_1ERS_0 <- 
  data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(MNRM_binned_GPCM_1ERS_0) <- 
  rownames(MNRM_binned_tree_1ERS_0) <- bins
colnames(MNRM_binned_GPCM_1ERS_0) <- 
  colnames(MNRM_binned_tree_1ERS_0) <- bins

for(i in 1:20){
  for(j in 1:20){
    MNRM_binned_tree_1ERS_0[i,j] <- mean(MNRM_Tree_df[MNRM_Tree_df[,3] == bins[j] &
                                                        MNRM_Tree_df[,4] == bins[i], 1])
    MNRM_binned_GPCM_1ERS_0[i,j] <- mean(MNRM_GPCM_df[MNRM_GPCM_df[,2] == bins[j] &
                                                        MNRM_GPCM_df[,3] == bins[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(MNRM_binned_GPCM_1ERS_0), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_tree_1ERS_0), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)















##########################################################################
#Bins based on MNRM data, threshold -1, 0, 1 with +1 ERS mean for group 2
#Intialize matrices
MNRM_true_traits_1ERS_1 <- MNRM_Tree_est_1ERS_1 <- MNRM_MNRM_est_1ERS_1 <- 
  matrix(NA, nrow = 250000, ncol = 2)
MNRM_GPCM_est_1ERS_1 <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  MNRM_true_traits_1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_1[[i]][["True_traits"]][501:1000,]
  MNRM_MNRM_est_1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_1[[i]][["MNRM_traits"]][501:1000,]
  MNRM_Tree_est_1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_1[[i]][["Tree_traits"]][501:1000,]
  MNRM_GPCM_est_1ERS_1[((i-1) * 500 + 1): (i * 500),] <- 
    test_MNRM_difficult_20_ERS_1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
MNRM_Tree_est_1ERS_1 <- MNRM_Tree_est_1ERS_1 - MNRM_MNRM_est_1ERS_1
MNRM_GPCM_est_1ERS_1 <- MNRM_GPCM_est_1ERS_1 - MNRM_MNRM_est_1ERS_1[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
MNRM_true_bins <- data.frame(True_theta = MNRM_true_traits_1ERS_1[,1], True_ERS = MNRM_true_traits_1ERS_1[,2],
                             theta_bin = cut(MNRM_true_traits_1ERS_1[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(MNRM_true_traits_1ERS_1[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
MNRM_Tree_df <- data.frame(Tree_theta_bias = MNRM_Tree_est_1ERS_1[,1], Tree_ERS_bias = MNRM_Tree_est_1ERS_1[,2],
                           theta_bin = MNRM_true_bins[,3], ERS_bin = MNRM_true_bins[,4])
MNRM_GPCM_df <- data.frame(GPCM_theta_bias = MNRM_GPCM_est_1ERS_1[,1], theta_bin = MNRM_true_bins[,3],
                           ERS_bin = MNRM_true_bins[,4])

#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
MNRM_Tree_df <- MNRM_Tree_df[is.na(MNRM_Tree_df[,3]) == FALSE & 
                               is.na(MNRM_Tree_df[,4]) == FALSE,]
MNRM_GPCM_df <- MNRM_GPCM_df[is.na(MNRM_GPCM_df[,2]) == FALSE & 
                               is.na(MNRM_GPCM_df[,3]) == FALSE,]

bins <- names(table(MNRM_true_bins[,3]))

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
MNRM_binned_tree_1ERS_1 <- MNRM_binned_MNRM_1ERS_1 <- MNRM_binned_GPCM_1ERS_1 <- 
  data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(MNRM_binned_GPCM_1ERS_1) <- 
  rownames(MNRM_binned_tree_1ERS_1) <- bins
colnames(MNRM_binned_GPCM_1ERS_1) <- 
  colnames(MNRM_binned_tree_1ERS_1) <- bins

for(i in 1:20){
  for(j in 1:20){
    MNRM_binned_tree_1ERS_1[i,j] <- mean(MNRM_Tree_df[MNRM_Tree_df[,3] == bins[j] &
                                                        MNRM_Tree_df[,4] == bins[i], 1])
    MNRM_binned_GPCM_1ERS_1[i,j] <- mean(MNRM_GPCM_df[MNRM_GPCM_df[,2] == bins[j] &
                                                        MNRM_GPCM_df[,3] == bins[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(MNRM_binned_GPCM_1ERS_1), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(MNRM_binned_tree_1ERS_1), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)









##################################################################
#IRTree as generating model
#Only group 2 included
#-1 ERS, no shift
#Intialize matrices
Tree_true_traits <- Tree_Tree_est <- Tree_MNRM_est <- 
  matrix(NA, nrow = 250000, ncol = 2)
Tree_GPCM_est <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  Tree_true_traits[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_min1[[i]][["True_traits"]][501:1000,]
  Tree_MNRM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_min1[[i]][["MNRM_traits"]][501:1000,]
  Tree_Tree_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_min1[[i]][["Tree_traits"]][501:1000,]
  Tree_GPCM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_min1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
Tree_MNRM_est <- Tree_MNRM_est - Tree_Tree_est
Tree_GPCM_est <- Tree_GPCM_est - Tree_Tree_est[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
Tree_true_bins <- data.frame(True_theta = Tree_true_traits[,1], True_ERS = Tree_true_traits[,2],
                             theta_bin = cut(Tree_true_traits[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(Tree_true_traits[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
Tree_GPCM_df <- data.frame(GPCM_theta_bias = Tree_GPCM_est[,1], theta_bin = Tree_true_bins[,3],
                           ERS_bin = Tree_true_bins[,4])
Tree_MNRM_df <- data.frame(MNRM_theta_bias = Tree_MNRM_est[,1], MNRM_ERS_bias = Tree_MNRM_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])

#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
Tree_MNRM_df <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE &
                               is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df <- Tree_GPCM_df[is.na(Tree_GPCM_df[,2]) == FALSE &
                               is.na(Tree_GPCM_df[,3]) == FALSE,]

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
Tree_binned_tree <- Tree_binned_MNRM <- Tree_binned_GPCM <- data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(Tree_binned_GPCM) <-
  rownames(Tree_binned_MNRM) <- names(table(Tree_true_bins[,4]))
colnames(Tree_binned_GPCM) <- 
  colnames(Tree_binned_MNRM) <- names(table(Tree_true_bins[,3]))

for(i in 1:20){
  for(j in 1:20){
    Tree_binned_MNRM[i,j] <- mean(Tree_MNRM_df[Tree_MNRM_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_MNRM_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_GPCM[i,j] <- mean(Tree_GPCM_df[Tree_GPCM_df[,2] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_GPCM_df[,3] == names(table(Tree_true_bins[,4]))[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(Tree_binned_GPCM), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_MNRM), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)








#Only group 2 included
#-1 ERS, with shift


#Intialize matrices
Tree_true_traits <- Tree_Tree_est <- Tree_MNRM_est <- 
  matrix(NA, nrow = 250000, ncol = 2)
Tree_GPCM_est <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  Tree_true_traits[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_min1[[i]][["True_traits"]][501:1000,]
  Tree_MNRM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_min1[[i]][["MNRM_traits"]][501:1000,]
  Tree_Tree_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_min1[[i]][["Tree_traits"]][501:1000,]
  Tree_GPCM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_min1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
Tree_MNRM_est <- Tree_MNRM_est - Tree_Tree_est
Tree_GPCM_est <- Tree_GPCM_est - Tree_Tree_est[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
Tree_true_bins <- data.frame(True_theta = Tree_true_traits[,1], True_ERS = Tree_true_traits[,2],
                             theta_bin = cut(Tree_true_traits[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(Tree_true_traits[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
Tree_GPCM_df <- data.frame(GPCM_theta_bias = Tree_GPCM_est[,1], theta_bin = Tree_true_bins[,3],
                           ERS_bin = Tree_true_bins[,4])
Tree_MNRM_df <- data.frame(MNRM_theta_bias = Tree_MNRM_est[,1], MNRM_ERS_bias = Tree_MNRM_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])

#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
Tree_MNRM_df <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE & is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df <- Tree_GPCM_df[is.na(Tree_GPCM_df[,2]) == FALSE & is.na(Tree_GPCM_df[,3]) == FALSE,]

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
Tree_binned_tree <- Tree_binned_MNRM <- Tree_binned_GPCM <- data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(Tree_binned_MNRM) <- rownames(Tree_binned_GPCM) <-
  names(table(Tree_true_bins[,4]))
colnames(Tree_binned_MNRM) <- colnames(Tree_binned_GPCM) <- 
  names(table(Tree_true_bins[,3]))

for(i in 1:20){
  for(j in 1:20){
    Tree_binned_MNRM[i,j] <- mean(Tree_MNRM_df[Tree_MNRM_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_MNRM_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_GPCM[i,j] <- mean(Tree_GPCM_df[Tree_GPCM_df[,2] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_GPCM_df[,3] == names(table(Tree_true_bins[,4]))[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(Tree_binned_GPCM), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_MNRM), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)









#IRTree as generating model
#Only group 2 included
#0 ERS, no shift
#Intialize matrices
Tree_true_traits <- Tree_Tree_est <- Tree_MNRM_est <- 
  matrix(NA, nrow = 250000, ncol = 2)
Tree_GPCM_est <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  Tree_true_traits[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20[[i]][["True_traits"]][501:1000,]
  Tree_MNRM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20[[i]][["MNRM_traits"]][501:1000,]
  Tree_Tree_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20[[i]][["Tree_traits"]][501:1000,]
  Tree_GPCM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
Tree_MNRM_est <- Tree_MNRM_est - Tree_Tree_est
Tree_GPCM_est <- Tree_GPCM_est - Tree_Tree_est[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
Tree_true_bins <- data.frame(True_theta = Tree_true_traits[,1], True_ERS = Tree_true_traits[,2],
                             theta_bin = cut(Tree_true_traits[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(Tree_true_traits[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
Tree_GPCM_df <- data.frame(GPCM_theta_bias = Tree_GPCM_est[,1], theta_bin = Tree_true_bins[,3],
                           ERS_bin = Tree_true_bins[,4])
Tree_MNRM_df <- data.frame(MNRM_theta_bias = Tree_MNRM_est[,1], MNRM_ERS_bias = Tree_MNRM_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])

#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
Tree_MNRM_df <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE & is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df <- Tree_GPCM_df[is.na(Tree_GPCM_df[,2]) == FALSE & is.na(Tree_GPCM_df[,3]) == FALSE,]

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
Tree_binned_tree <- Tree_binned_MNRM <- Tree_binned_GPCM <- data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(Tree_binned_MNRM) <- rownames(Tree_binned_GPCM) <- 
  names(table(Tree_true_bins[,4]))
colnames(Tree_binned_MNRM) <- colnames(Tree_binned_GPCM) <-
  names(table(Tree_true_bins[,3]))

for(i in 1:20){
  for(j in 1:20){
    Tree_binned_MNRM[i,j] <- mean(Tree_MNRM_df[Tree_MNRM_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_MNRM_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_GPCM[i,j] <- mean(Tree_GPCM_df[Tree_GPCM_df[,2] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_GPCM_df[,3] == names(table(Tree_true_bins[,4]))[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(Tree_binned_GPCM), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_MNRM), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)








#Only group 2 included
#-1 ERS, with shift


#Intialize matrices
Tree_true_traits <- Tree_Tree_est <- Tree_MNRM_est <- 
  matrix(NA, nrow = 250000, ncol = 2)
Tree_GPCM_est <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  Tree_true_traits[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20[[i]][["True_traits"]][501:1000,]
  Tree_MNRM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20[[i]][["MNRM_traits"]][501:1000,]
  Tree_Tree_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20[[i]][["Tree_traits"]][501:1000,]
  Tree_GPCM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
Tree_MNRM_est <- Tree_MNRM_est - Tree_Tree_est
Tree_GPCM_est <- Tree_GPCM_est - Tree_Tree_est[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
Tree_true_bins <- data.frame(True_theta = Tree_true_traits[,1], True_ERS = Tree_true_traits[,2],
                             theta_bin = cut(Tree_true_traits[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(Tree_true_traits[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
Tree_GPCM_df <- data.frame(GPCM_theta_bias = Tree_GPCM_est[,1], theta_bin = Tree_true_bins[,3],
                           ERS_bin = Tree_true_bins[,4])
Tree_MNRM_df <- data.frame(MNRM_theta_bias = Tree_MNRM_est[,1], MNRM_ERS_bias = Tree_MNRM_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])

#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
Tree_MNRM_df <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE & is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df <- Tree_GPCM_df[is.na(Tree_GPCM_df[,2]) == FALSE & is.na(Tree_GPCM_df[,3]) == FALSE,]

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
Tree_binned_tree <- Tree_binned_MNRM <- Tree_binned_GPCM <- data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(Tree_binned_MNRM) <- rownames(Tree_binned_GPCM) <-
  names(table(Tree_true_bins[,4]))
colnames(Tree_binned_MNRM) <- colnames(Tree_binned_GPCM) <- 
  names(table(Tree_true_bins[,3]))

for(i in 1:20){
  for(j in 1:20){
    Tree_binned_MNRM[i,j] <- mean(Tree_MNRM_df[Tree_MNRM_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_MNRM_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_GPCM[i,j] <- mean(Tree_GPCM_df[Tree_GPCM_df[,2] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_GPCM_df[,3] == names(table(Tree_true_bins[,4]))[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(Tree_binned_GPCM), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_MNRM), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)












#IRTree as generating model
#Only group 2 included
#+1 ERS, no shift
#Intialize matrices
Tree_true_traits <- Tree_Tree_est <- Tree_MNRM_est <- 
  matrix(NA, nrow = 250000, ncol = 2)
Tree_GPCM_est <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  Tree_true_traits[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_1[[i]][["True_traits"]][501:1000,]
  Tree_MNRM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_1[[i]][["MNRM_traits"]][501:1000,]
  Tree_Tree_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_1[[i]][["Tree_traits"]][501:1000,]
  Tree_GPCM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_average_20_ERS_1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
Tree_MNRM_est <- Tree_MNRM_est - Tree_Tree_est
Tree_GPCM_est <- Tree_GPCM_est - Tree_Tree_est[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
Tree_true_bins <- data.frame(True_theta = Tree_true_traits[,1], True_ERS = Tree_true_traits[,2],
                             theta_bin = cut(Tree_true_traits[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(Tree_true_traits[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
Tree_GPCM_df <- data.frame(GPCM_theta_bias = Tree_GPCM_est[,1], theta_bin = Tree_true_bins[,3],
                           ERS_bin = Tree_true_bins[,4])
Tree_MNRM_df <- data.frame(MNRM_theta_bias = Tree_MNRM_est[,1], MNRM_ERS_bias = Tree_MNRM_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])

#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
Tree_MNRM_df <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE & is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df <- Tree_GPCM_df[is.na(Tree_GPCM_df[,2]) == FALSE & is.na(Tree_GPCM_df[,3]) == FALSE,]

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
Tree_binned_tree <- Tree_binned_MNRM <- Tree_binned_GPCM <- data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(Tree_binned_MNRM) <- rownames(Tree_binned_GPCM) <-
  names(table(Tree_true_bins[,4]))
colnames(Tree_binned_MNRM) <- colnames(Tree_binned_GPCM) <-
  names(table(Tree_true_bins[,3]))

for(i in 1:20){
  for(j in 1:20){
    Tree_binned_MNRM[i,j] <- mean(Tree_MNRM_df[Tree_MNRM_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_MNRM_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_GPCM[i,j] <- mean(Tree_GPCM_df[Tree_GPCM_df[,2] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_GPCM_df[,3] == names(table(Tree_true_bins[,4]))[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(Tree_binned_GPCM), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_MNRM), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)








#Only group 2 included
#+1 ERS, with shift


#Intialize matrices
Tree_true_traits <- Tree_Tree_est <- Tree_MNRM_est <- 
  matrix(NA, nrow = 250000, ncol = 2)
Tree_GPCM_est <- matrix(NA, nrow = 250000, ncol = 1)

#Extract true, MNRM, IRTree and GPCM theta and ERS values from data
for(i in 1:500){
  Tree_true_traits[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_1[[i]][["True_traits"]][501:1000,]
  Tree_MNRM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_1[[i]][["MNRM_traits"]][501:1000,]
  Tree_Tree_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_1[[i]][["Tree_traits"]][501:1000,]
  Tree_GPCM_est[((i-1) * 500 + 1): (i * 500),] <- 
    test_Tree_difficult_20_ERS_1[[i]][["GPCM_traits"]][501:1000,]
}
#Change theta and ERS values into difference between true and estimated values
Tree_MNRM_est <- Tree_MNRM_est - Tree_Tree_est
Tree_GPCM_est <- Tree_GPCM_est - Tree_Tree_est[,1]

#Create dataframe of true values and bins true values fall into of width 0.2 from -2 to 2
Tree_true_bins <- data.frame(True_theta = Tree_true_traits[,1], True_ERS = Tree_true_traits[,2],
                             theta_bin = cut(Tree_true_traits[,1], seq(from = -2, to = 2, by = 0.2)),
                             ERS_bin = cut(Tree_true_traits[,2], seq(from = -2, to = 2, by = 0.2)))

#Append bins based on true values to IRTree difference between true and estimated values
Tree_GPCM_df <- data.frame(GPCM_theta_bias = Tree_GPCM_est[,1], theta_bin = Tree_true_bins[,3],
                           ERS_bin = Tree_true_bins[,4])
Tree_MNRM_df <- data.frame(MNRM_theta_bias = Tree_MNRM_est[,1], MNRM_ERS_bias = Tree_MNRM_est[,2],
                           theta_bin = Tree_true_bins[,3], ERS_bin = Tree_true_bins[,4])

#Remove all values where one of the bins is NA (value lower than -2 or higher than +2)
#Consider if you want to remove these or bundle them into the highest/lowest bins
Tree_MNRM_df <- Tree_MNRM_df[is.na(Tree_MNRM_df[,3]) == FALSE & is.na(Tree_MNRM_df[,4]) == FALSE,]
Tree_GPCM_df <- Tree_GPCM_df[is.na(Tree_GPCM_df[,2]) == FALSE & is.na(Tree_GPCM_df[,3]) == FALSE,]

#Obtain the average difference and sum of squared deviation for each of the 400
#bins that were created
Tree_binned_tree <- Tree_binned_MNRM <- Tree_binned_GPCM <- data.frame(matrix(NA, nrow = 20, ncol = 20))
rownames(Tree_binned_MNRM) <- rownames(Tree_binned_GPCM) <-
  names(table(Tree_true_bins[,4]))
colnames(Tree_binned_MNRM) <- colnames(Tree_binned_GPCM) <-
  names(table(Tree_true_bins[,3]))

for(i in 1:20){
  for(j in 1:20){
    Tree_binned_MNRM[i,j] <- mean(Tree_MNRM_df[Tree_MNRM_df[,3] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_MNRM_df[,4] == names(table(Tree_true_bins[,4]))[i], 1])
    Tree_binned_GPCM[i,j] <- mean(Tree_GPCM_df[Tree_GPCM_df[,2] == names(table(Tree_true_bins[,3]))[j] &
                                                 Tree_GPCM_df[,3] == names(table(Tree_true_bins[,4]))[i], 1])
  }
}

colors = c(seq(-2, -0.2,length=20),seq(-0.2,0.2,length=4),seq(0.2,2,length=20))

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 43)

par(mar = c(6.5, 6.5, 0.5, 0.5), mgp = c(3, 0.7, 1))

my_pallete2 <- colorRampPalette(c("blue", "white", "red"))(n = 20)

heatmap(as.matrix(Tree_binned_GPCM), Rowv = NA, Colv = NA, scale = "none", 
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
heatmap(as.matrix(Tree_binned_MNRM), Rowv = NA, Colv = NA, scale = "none",
        col = my_palette, breaks = colors, xlab = "Theta", ylab = " ERS",
        RowSideColors = my_pallete2)
