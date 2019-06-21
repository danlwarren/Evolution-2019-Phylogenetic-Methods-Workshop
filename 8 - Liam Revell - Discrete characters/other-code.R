library(phytools)








set.seed(999)
sim.tree<-pbtree(n=90,scale=1)
plotTree(sim.tree,ftype="off")
nodelabels()

sim.tree<-paintSubTree(sim.tree,93,"b","a",0.5)
plot(sim.tree,setNames(c("blue","red"),c("a","b")),
	ftype="off")
Q<-setNames(list(
	0.5*matrix(c(-2,1,1,1,-2,1,1,1,-2),3,3,
	dimnames=list(LETTERS[1:3],LETTERS[1:3])),
	2.0*matrix(c(-2,1,1,1,-2,1,1,1,-2),3,3,
	dimnames=list(LETTERS[1:3],LETTERS[1:3]))),
	letters[1:2])
x<-sim.multiMk(sim.tree,Q)

write.tree(sim.tree,"simulated-tree.tre")
write.csv(x,"simulated-data.csv")

fitmultiMk(sim.tree,x,model="SYM")

## pure-birth phylogeny
tree<-pbtree(n=200,scale=1)

## 4-state 'ER' ordered
Q<-matrix(c(-1,1,0,0,0,0,0,
    1,-2,1,0,0,0,0,
    0,1,-2,1,0,0,0,
    0,0,1,-2,1,0,0,
    0,0,0,1,-2,1,0,
    0,0,0,0,1,-2,1,
    0,0,0,0,0,1,-1),7,7,byrow=TRUE)
rownames(Q)<-colnames(Q)<-c("A","A+B","B","B+C","C","C+D","D")
y<-sim.Mk(tree,Q)

write.csv(file="polymorphic-data.csv",y)
write.tree(file="polymorphic-tree.phy",tree)


