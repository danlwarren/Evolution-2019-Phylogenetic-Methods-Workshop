require(treeplyr)
tree <- geiger::sim.bdtree(stop="taxa", n=50, seed=1)
dat <- data.frame(X1 = rnorm(50), X2 = rnorm(50), X3 = rnorm(50),
                  taxa=c(tree$tip.label[1:40], paste("s", 51:60, sep="")), 
                  D1=sample(c("Hello", "World"), 50, replace=TRUE), D2=rbinom(50,1,0.5), 
                  D3=0.3+rbinom(50, 4, c(0.1,0.1,0.5,0.3)), XNA1 = c(rep(NA, 10), rnorm(40)))

write.tree(tree, "~/Documents/teaching/Evolution-2019-Phylogenetic-Methods-Workshop/4 - Luke Harmon/Data/sampleTree.phy")
write.csv(dat, "~/Documents/teaching/Evolution-2019-Phylogenetic-Methods-Workshop/4 - Luke Harmon/Data/sampleData.csv")

td <- make.treedata(tree, dat)


summary(td)

filter(td, X1 + X2 > 0 & D1 == "Hello")

