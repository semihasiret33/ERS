if(!require(tidyverse)){install.packages("tidyverse")}
if(!require(irtrees)){install.packages("irtrees")}
if(!require(reshape)){install.packages("reshape")}
if(!require(foreach)){install.packages("foreach")}
if(!require(doParallel)){install.packages("doParallel")}
if(!require(installr)){install.packages("installr")}
if(!require(mirt)){install.packages("mirt")}


#Read in the PISA data
mydata <- read.table("INT_STU12_DEC03.txt", header = FALSE, sep = "'")

a <- separate(mydata, col = 1, into = c("COUNTRY", "SUBNATION", "STRATNUM", 
                                        "OECD", "NC", "SCHOOLID", "STUDENT ID",
                                        "INT GRAD", "NAT STUD PROG", "BIRTHMONTH",
                                        "BIRTHYEAR", "GENDER", "Unused questions",
                                        "MATHQ1", "MATHQ2", "MATHQ3", "MATHQ4",
                                        "MATHQ5", "MATHQ6", "MATHQ7", "MATHQ8",
                                        "MATHQ9", "More Unused questions"), 
              sep = c(3, 10, 17, 18, 24, 31, 36, 38, 40, 42, 46, 47, 169,
                      170, 171, 172, 173, 174, 175, 176, 177, 178))

mydata <- data.frame(COUNTRY = a[,1],
                     SUBNATION = a[,2],
                     STRATNUM = a[,3],
                     OECD = a[,4],
                     NC = a[,5],
                     SCHOOLID = a[,6],
                     STUDENT_ID = a[,7],
                     INT_GRAD = a[,8],
                     NAT_STUD_PROG = a[,9],
                     BIRTHMONTH = a[,10],
                     BIRTHYEAR = a[,11],
                     GENDER = a[,12],
                     MATHQ1 = as.numeric(a[,14]),
                     MATHQ2 = as.numeric(a[,15]),
                     MATHQ3 = as.numeric(a[,16]),
                     MATHQ4 = as.numeric(a[,17]),
                     MATHQ5 = as.numeric(a[,18]),
                     MATHQ6 = as.numeric(a[,19]),
                     MATHQ7 = as.numeric(a[,20]),
                     MATHQ8 = as.numeric(a[,21]),
                     MATHQ9 = as.numeric(a[,22]))
rm(a)

#Make sure all missing values (7, 8 and 9) are converted to NA
mydata$MATHQ1[mydata$MATHQ1 > 4] <- NA
mydata$MATHQ2[mydata$MATHQ2 > 4] <- NA
mydata$MATHQ3[mydata$MATHQ3 > 4] <- NA
mydata$MATHQ4[mydata$MATHQ4 > 4] <- NA
mydata$MATHQ5[mydata$MATHQ5 > 4] <- NA
mydata$MATHQ6[mydata$MATHQ6 > 4] <- NA
mydata$MATHQ7[mydata$MATHQ7 > 4] <- NA
mydata$MATHQ8[mydata$MATHQ8 > 4] <- NA
mydata$MATHQ9[mydata$MATHQ9 > 4] <- NA

#Reverse code the data to make sure higher scores reflect higher latent trait
#values
mydata$MATHQ1 <- 5 - mydata$MATHQ1
mydata$MATHQ2 <- 5 - mydata$MATHQ2 
mydata$MATHQ3 <- 5 - mydata$MATHQ3 
mydata$MATHQ4 <- 5 - mydata$MATHQ4 
mydata$MATHQ5 <- 5 - mydata$MATHQ5 
mydata$MATHQ6 <- 5 - mydata$MATHQ6 
mydata$MATHQ7 <- 5 - mydata$MATHQ7 
mydata$MATHQ8 <- 5 - mydata$MATHQ8 
mydata$MATHQ9 <- 5 - mydata$MATHQ9 


#Correlation matrix
cor(mydata[,13:21], use = "pairwise.complete.obs")

#Calculate overall test score with participants scoring all NA removed
mydata_no_NA <- mydata[rowSums(is.na(mydata[,13:21])) != 
                                 ncol(mydata[,13:21]), ]

mean(rowMeans(mydata_no_NA[,13:21], na.rm = TRUE))

#Take the subset of countries for the example. Country 13 is Costa Rica,
#Country 42 is Malaysia
mydata_subset <- mydata[mydata[,1] == names(table(mydata[,1]))[13] | 
                          mydata[,1] == names(table(mydata[,1]))[42], ]

mydata_subset <- mydata_subset[rowSums(is.na(mydata_subset[,13:21])) != 
                                 ncol(mydata_subset[,13:21]), ]

Nitems <- 9

model <- paste0("Theta = 1-", Nitems)

#Run the GPCM model. Forward SE is used as not all SE types work with the 
#IRTree and we wanted consistency on how the SE's are calculated between
#models
GPCM_test_mod_13_42 <- multipleGroup(mydata_subset[, 13:21],
                                     model = model,
                                     group = mydata_subset[, 1],
                                     itemtype = "gpcm",
                                     invariance = c("free_mean",
                                                    "free_var",
                                                    "slopes",
                                                    "intercepts"),
                                     method = "EM",
                                     technical = list(NCYCLES = 2000),
                                     SE = TRUE,
                                     SE.type = "forward")


#Prepare to run MNRM model
model <- paste0("Theta = 1-", Nitems, "
                ", "ERS = 1-", Nitems,"
                FREE = (GROUP, COV_21)")


#Generate s matrix
s_gen <- matrix(0:(4 - 1), nrow = 1)
#Generate ERS scoring matrix
ers_mat <- matrix(c(1,0,0,1), nrow = 1)

s_gen <- rbind(s_gen, ers_mat)

s_mirt <- vector(mode = "list", length = Nitems)

#Generate item names and s_mirt matrix to evaluate model
for(j in 1:Nitems){
  s_mirt[[j]] <- t(s_gen)
}

#Run MNRM model, again with forward SE's for consistency with IRTree
MNRM_test_mod_13_42 <- multipleGroup(mydata_subset[, 13:21],
                               model = model,
                               group = mydata_subset[, 1],
                               itemtype = "gpcm",
                               invariance = c("free_mean",
                                              "free_var",
                                              "slopes",
                                              "intercepts"),
                               method = "EM",
                               gpcm_mats = s_mirt,
                               technical = list(NCYCLES = 2000),
                               SE = TRUE,
                               SE.type = "forward")

#Calculate mean slopes and item thresholds for the MNRM
itempars <- coef(MNRM_test_mod_13_42)[[2]]

itempar_mat <- matrix(NA, nrow = 9, ncol = 6)
for(i in 1:9){
  itempar_mat[i, 1] <- itempars[[i]][1]
  itempar_mat[i, 2] <- itempars[[i]][4]
  itempar_mat[i, 3] <- itempars[[i]][31]
  itempar_mat[i, 4] <- itempars[[i]][34]
  itempar_mat[i, 5] <- itempars[[i]][37]
  itempar_mat[i, 6] <- itempars[[i]][40]
}

itempar_mat[,3] <- (itempar_mat[,3] - itempar_mat[,4])/1.235291 
itempar_mat[,4] <- (itempar_mat[,4] - itempar_mat[,5])/1.235291 
itempar_mat[,5] <- (itempar_mat[,5] - itempar_mat[,6])/1.235291 

colMeans(itempar_mat)



#Prepare to run IRTree model
name <- "Node2TreePL"
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
                          derivType = "symbolic", derivType.hss = "symbolic")

#Custom item function for Node 3
name <- "Node3TreePL"
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
                          derivType = "symbolic", derivType.hss = "symbolic")



#Manually create tree nodes
tree_subdat <- matrix(NA, nrow = nrow(mydata_subset), ncol = (3*Nitems))
for(i in 1:Nitems){
  tree_subdat[mydata_subset[,i + 12] <= 2, (i-1)*3 + 1] <- 0 
  tree_subdat[mydata_subset[,i + 12] >= 3, (i-1)*3 + 1] <- 1
  
  tree_subdat[mydata_subset[,i + 12] == 2, (i-1)*3 + 2] <- 1 
  tree_subdat[mydata_subset[,i + 12] == 1, (i-1)*3 + 2] <- 0
  
  tree_subdat[mydata_subset[,i + 12] == 3, (i-1)*3 + 3] <- 0 
  tree_subdat[mydata_subset[,i + 12] == 4, (i-1)*3 + 3] <- 1 
}

items <- sort(rep(1:Nitems, 3))
items <- paste0("item", items)
nodes <- 1:3
nodes <- rep(paste0("node", nodes), Nitems)
colnames(tree_subdat) <- paste0(items, nodes)

treeNitems <- 3*Nitems

ersloadings <- 1:(2 * Nitems)
for(i in 1:Nitems){
  ersloadings[(i-1)*2+1] <- (i-1)*3+2
  ersloadings[(i-1)*2+2] <- (i-1)*3+3
}

node1alpha <- node2alpha <- node3alpha <- 1:Nitems
for(i in 1:Nitems){
  node1alpha[i] <- (i-1)*3+1
  node2alpha[i] <- (i-1)*3+2
  node3alpha[i] <- (i-1)*3+3
}

model <- paste0("Theta = 1-", treeNitems, "
                ", "ERS = ", paste0(ersloadings, collapse = ","), "
                ", "CONSTRAIN = ", paste0("(",node2alpha, ",", node3alpha, ", a1)", collapse = ","), ",", "
                     ", paste0("(",node2alpha, ",", node3alpha, ", a2)", collapse = ","), "
                     ")

model2 <- paste0("FREE = (GROUP, COV_21)",
                 collapse = "")

model <- paste0(model, model2, collapse = "")


#Run the IRTree model with forward SE's. This is the SE type that is fastest
#and seems to work consistently with the Tree
tree_submod_constr_13_42 <- multipleGroup(tree_subdat,
                                    model = model,
                                    group = mydata_subset[,1],
                                    itemtype = rep(c("2PL", 
                                                     "Node2TreePL", 
                                                     "Node3TreePL"), Nitems),
                                    method = "EM",
                                    accelerate = "Ramsay",
                                    invariance = c("free_mean",
                                                   "free_var",
                                                   "slopes",
                                                   "intercepts"),
                                    technical = list(NCYCLES = 2000),
                                    customItems = list(Node2TreePL = Node2TreePL,
                                                       Node3TreePL = Node3TreePL),
                                    SE = TRUE,
                                    SE.type = "forward")
