<h2>Some methods for the comparative analysis of discrete character data</h2>


```r
library(phytools)
packageVersion("phytools") ## installed from GitHub 0.6-99 is fine
```

```
## [1] '0.7.0'
```

<h3><i>Method 1</i>: Testing for heterogeneous rates of discrete character 
evolution on phylogenies.</h3>

In this first of three short modules, we'll see how to fit a extended 
M<i>k</i> model in which the rate of character evolution varies as a 
function of 'regimes' mapped onto the tree.

We can do this first using some simulated data that can be downloaded here:

1. <a href="simulated-tree.tre">simulated-tree.tre</a>
2. <a href="simulated-data.csv">simulated-data.csv</a>


```r
sim.tree<-read.tree("simulated-tree.tre")
sim.data<-read.csv("simulated-data.csv",row.names=1)
```

Plot our tree:


```r
plotTree(sim.tree,ftype="off")
nodelabels(cex=0.6)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

Now let's imagine that we had some <i>a priori</i> reason, independent of our
character data, to suspect that clade descended from node '93' had evolved 
under a different evolution. Let's map this hypothesis on the tree, jointly
with our trait data at the tips:


```r
hypothesis<-paintSubTree(sim.tree,93,
	"regime 2","regime 1")
cols<-setNames(c("blue","red"),c("regime 1","regime 2"))
plot(hypothesis,cols,ftype="off")
trait<-setNames(sim.data[,1],rownames(sim.data))
tiplabels(pie=to.matrix(trait[hypothesis$tip.label],
	levels(trait)),cex=0.3,
	piecol=c("white","grey","black"))
legend("topleft",pch=21,pt.cex=2,
	pt.bg=c("white","grey","black"),
	legend=levels(trait),bty="n")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

Let's fit two models to the data as follows:


```r
singleRateModel<-fitMk(hypothesis,trait,
	model="ER")
singleRateModel
```

```
## Object of class "fitMk".
## 
## Fitted (or set) value of Q:
##           A         B         C
## A -1.431319  0.715659  0.715659
## B  0.715659 -1.431319  0.715659
## C  0.715659  0.715659 -1.431319
## 
## Fitted (or set) value of pi:
##         A         B         C 
## 0.3333333 0.3333333 0.3333333 
## 
## Log-likelihood: -64.975367 
## 
## Optimization method used was "nlminb"
```

```r
multiRateModel<-fitmultiMk(hypothesis,trait,
	model="ER")
multiRateModel
```

```
## Object of class "fitmultiMk".
## 
## Fitted value of Q[regime 1]:
##           A         B         C
## A -0.689337  0.344668  0.344668
## B  0.344668 -0.689337  0.344668
## C  0.344668  0.344668 -0.689337
## 
## Fitted value of Q[regime 2]:
##           A         B         C
## A -3.560624  1.780312  1.780312
## B  1.780312 -3.560624  1.780312
## C  1.780312  1.780312 -3.560624
## 
## Fitted (or set) value of pi:
##         A         B         C 
## 0.3333333 0.3333333 0.3333333 
## 
## Log-likelihood: -60.517811 
## 
## Optimization method used was "nlminb"
```

These models can be compared directly:


```r
data.frame(model=c("single-rate","multi-rate"),
	logLik=c(logLik(singleRateModel),
	logLik(multiRateModel)),
	k=c(attr(AIC(singleRateModel),"df"),
	attr(AIC(multiRateModel),"df")),
	AIC=c(AIC(singleRateModel),
	AIC(multiRateModel)))
```

```
##         model    logLik k      AIC
## 1 single-rate -64.97537 1 131.9507
## 2  multi-rate -60.51781 2 125.0356
```

This tells us that the multi-rate model is much better supported than the
single rate model.

Let's proceed to plot our best-fitting model as follows:


```r
layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
plot(hypothesis,cols,ftype="off",
	direction="downwards")
tiplabels(pie=to.matrix(trait[hypothesis$tip.label],
	levels(trait)),cex=0.3,
	piecol=c("white","grey","black"))
legend("topleft",pch=21,pt.cex=2,
	pt.bg=c("white","grey","black"),
	legend=levels(trait),bty="n")
obj<-multiRateModel
obj$rates<-round(obj$rates[1],2)
obj$regimes<-NULL
class(obj)<-"fitMk"
plot(obj,show.zeros=FALSE,
	mar=rep(2.1,4),show.zeros=FALSE,
	tol=1e-3,cex.traits=0.9,
	cex.rates=0.6)
mtext(text="a) Regime 1 (blue)",
    adj=0,line=-1,cex=0.9)

obj<-multiRateModel
obj$rates<-round(obj$rates[2],2)
obj$regimes<-NULL
class(obj)<-"fitMk"
plot(obj,show.zeros=FALSE,
	mar=rep(2.1,4),show.zeros=FALSE,
	tol=1e-3,cex.traits=0.9,
	cex.rates=0.6)
mtext(text="b) Regime 2 (red)",
    adj=0,line=-1,cex=0.9)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

Now we can try the same thing with a real data set consisting of a phylogeny
of <i>Anolis</i> lizards and  the number of vertebrae in the tail:

1. <a href="VERT.CSV">VERT.CSV</a>
2. <a href="ANOLIS.PHY">ANOLIS.PHY</a>

Start by reading the data from file:


```r
X<-read.csv("VERT.CSV",row.names=1)
vert<-factor(setNames(X[,1],rownames(X)),
	levels=min(X[,1]):max(X[,1]))
```

Now we can create the model that we want to fit. For this step, I propose
an ordered reversible model - in which caudal vertebrae are gained & lost 
sequentially - although possibly with different rates in either direction:


```r
k<-length(levels(vert))
ordered<-matrix(0,k,k,dimnames=list(levels(vert),
	levels(vert)))
for(i in 1:k){
	if(i<k) ordered[i,i+1]<-1
	if(i>1) ordered[i,i-1]<-2
}
ordered[i,i-1]<-2
```

This is the design matrix of our model:


```r
ordered
```

```
##    34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55
## 34  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
## 35  2  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
## 36  0  2  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
## 37  0  0  2  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
## 38  0  0  0  2  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
## 39  0  0  0  0  2  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
## 40  0  0  0  0  0  2  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
## 41  0  0  0  0  0  0  2  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
## 42  0  0  0  0  0  0  0  2  0  1  0  0  0  0  0  0  0  0  0  0  0  0
## 43  0  0  0  0  0  0  0  0  2  0  1  0  0  0  0  0  0  0  0  0  0  0
## 44  0  0  0  0  0  0  0  0  0  2  0  1  0  0  0  0  0  0  0  0  0  0
## 45  0  0  0  0  0  0  0  0  0  0  2  0  1  0  0  0  0  0  0  0  0  0
## 46  0  0  0  0  0  0  0  0  0  0  0  2  0  1  0  0  0  0  0  0  0  0
## 47  0  0  0  0  0  0  0  0  0  0  0  0  2  0  1  0  0  0  0  0  0  0
## 48  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0  1  0  0  0  0  0  0
## 49  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0  1  0  0  0  0  0
## 50  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0  1  0  0  0  0
## 51  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0  1  0  0  0
## 52  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0  1  0  0
## 53  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0  1  0
## 54  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0  1
## 55  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  2  0
```

Now let's convert our vector of vertebrae to a matrix. We do this because
there are certain levels of our trait that are not present in the input
data vector:


```r
vert<-to.matrix(vert,levels(vert))
head(vert)
```

```
##             34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54
## ahli         0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
## allogus      0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0
## rubribarbus  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0
## imias        0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0
## sagrei       0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0
## bremeri      0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0
##             55
## ahli         0
## allogus      0
## rubribarbus  0
## imias        0
## sagrei       0
## bremeri      0
```

Our next step is to propose our hypothesis for rate variation on the tree.
In this case, to keep it simple, we propose that mainland & Caribbean
anoles have different rates. I have saved this hypothesis (mapped on the 
phylogeny) to our tree file:


```r
anolis.tree<-read.simmap("ANOLIS.PHY",format="phylip")
plot(anolis.tree,ftype="off",type="fan",
	colors=setNames(c("blue","brown"),
	c("I","M")))
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

Unfortunately, the taxa in our tree & our data matrix don't match exactly.
We can identify mismatches using <code>geiger::name.check</code> & fix 
them as follows:


```r
library(geiger)
chk<-name.check(anolis.tree,vert)
anolis.tree<-drop.tip.simmap(anolis.tree,
	chk$tree_not_data)
vert<-vert[anolis.tree$tip.label,]
name.check(anolis.tree,vert)
```

```
## [1] "OK"
```

Now we should be ready to fit our two models:


```r
fit.single<-fitMk(anolis.tree,vert,model=ordered)
fit.multi<-fitmultiMk(anolis.tree,vert,model=ordered)
```

Here's a summary of our fitted models:


```r
data.frame(model=c("single-rate","multi-rate"),
	logLik=c(logLik(fit.single),
	logLik(fit.multi)),
	k=c(attr(AIC(fit.single),"df"),
	attr(AIC(fit.multi),"df")),
	AIC=c(AIC(fit.single),
	AIC(fit.multi)))
```

```
##         model    logLik k      AIC
## 1 single-rate -280.4681 2 564.9362
## 2  multi-rate -279.8423 4 567.6846
```

Here, the two-rate model <i>is not</i> justified. Let's nonetheless
graph our fitted model as we did in the previous part of the exercise:


```r
par(mfrow=c(2,1))
obj<-fit.multi
obj$rates<-round(obj$rates[3:4],2)
obj$regimes<-NULL
class(obj)<-"fitMk"
plot(obj,show.zeros=F,mar=rep(2.1,4),show.zeros=F,
	tol=1e-3,cex.traits=0.8,cex.rates=0.4)
mtext(text="a) Mainland caudal vertebra number",
    adj=0,line=-1,cex=0.9)
obj<-fit.multi
obj$rates<-round(obj$rates[1:2],2)
obj$regimes<-NULL
class(obj)<-"fitMk"
plot(obj,show.zeros=F,mar=rep(2.1,4),show.zeros=F,
	tol=1e-3,cex.traits=0.8,cex.rates=0.4)
mtext(text="b) Caribbean caudal vertebra number",
    adj=0,line=-1,cex=0.9)
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)

<h3><i>Method 2</i>: Fitting discrete character evolution to data with 
intraspecific polymorphism.</h3>

For the second module we can use a new function called 
<code>fitpolyMk</code>. This function handles discrete character data
with intraspecific polymorphism in a fairly obvious way - that is merely as 
an intermediate state between the two or more observed character values.

The data for this part of the exercise can be downloaded here:

1. <a href="polymorphic-data.csv">polymorphic-data.csv</a>
2. <a href="polymorphic-tree.phy">polymorphic-tree.phy</a>


```r
poly.tree<-read.tree("polymorphic-tree.phy")
poly.data<-read.csv(file="polymorphic-data.csv",row.names=1)
y<-setNames(poly.data[,1],rownames(poly.data))
```

These are what our data should look like:


```r
y
```

```
##  t70  t71  t47 t169 t170  t67  t68  t62  t59  t60 t173 t174  t19  t20  t21 
##  C+D  C+D  C+D  C+D  C+D  C+D  B+C  C+D    C    C    B    B    B  B+C    B 
##  t15  t50 t109 t123 t124 t104 t105  t52  t53 t167 t168 t141  t74  t66  t48 
##    B  B+C  B+C  B+C  B+C  B+C  B+C  B+C  B+C    C    C    C    C  C+D    C 
## t181 t182   t2  t39  t40   t4  t27 t191 t192  t28 t145 t146  t77   t8 t193 
##  C+D  C+D    C    C    C    D    C  B+C    C    D  C+D  C+D  C+D    D    D 
## t194  t96  t45  t86  t87  t13   t6 t125 t183 t184  t16  t17 t135 t136 t147 
##    D    D    D    D    D  C+D  C+D  C+D    C    C    C    C    C    C    C 
## t148  t83  t88  t89  t78  t79  t24 t187 t188  t36  t29  t99 t189 t190 t115 
##    C    C  B+C    C    C    C    C    C    C    C    D  C+D    C    C    C 
## t119 t120 t121 t122  t46  t49 t165 t166 t112   t5 t106 t126 t142 t143  t98 
##    C    C    C    C    C    C  C+D  C+D  C+D  B+C  B+C  B+C  B+C  B+C    C 
## t179 t180 t144 t113 t114  t18  t30 t155 t156  t11 t128 t129 t127  t64  t65 
##  B+C  B+C  B+C  B+C  B+C  B+C    C    C    C  B+C    C  C+D    C  C+D  C+D 
## t153 t154 t185 t186  t80  t81  t33  t43 t151 t152  t82  t97 t132 t133  t61 
##    C    C    C    C    B    B    C    D  C+D  C+D  C+D    D    D    D    D 
##   t3 t102 t103  t35  t69 t149 t150   t7  t92  t93 t110 t111 t197 t198 t161 
##  A+B  B+C  B+C  B+C  A+B    B    B  C+D  C+D  C+D  C+D  C+D    C    C    C 
## t157 t158  t12   t9  t41  t42  t25  t37  t38  t26  t10  t90  t94  t95  t54 
##    C    C  B+C    C    C    C    C    C    C    C  B+C  C+D  C+D  C+D    C 
##  t55  t91 t171 t172 t162 t163  t51  t32 t177 t178  t63  t84  t85  t44  t56 
##    C    D    C    C  C+D  C+D  C+D    C    C    C  B+C    C    C    C    C 
## t199 t200 t107 t108 t130 t131  t31 t137 t138  t22 t139 t140 t118  t72  t73 
##    C    C    C    C    C    C    C  C+D    C    C  C+D  C+D  C+D  C+D    C 
##  t76 t159 t160  t14  t75 t134 t195 t196 t164  t34 t100 t101 t116 t117  t23 
##  C+D  C+D  C+D    D  C+D  C+D  C+D  C+D  C+D    D    D    D  C+D    C  C+D 
## t175 t176  t57  t58   t1 
##  C+D  C+D  B+C  B+C  C+D 
## Levels: A+B B B+C C C+D D
```

Now let's plot them:


```r
plotTree(poly.tree,ftype="off",lwd=1,type="fan")
X<-strsplit(setNames(as.character(y),names(y)),"+",fixed=TRUE)
pies<-matrix(0,Ntip(poly.tree),4,dimnames=list(poly.tree$tip.label,
	c("A","B","C","D")))
for(i in 1:Ntip(poly.tree)) 
	pies[poly.tree$tip.label[i],X[[poly.tree$tip.label[i]]]]<-
		rep(1/length(X[[poly.tree$tip.label[i]]]),
		length(X[[poly.tree$tip.label[i]]]))
tiplabels(pie=pies,piecol=c("black","yellow","red","blue"),
	cex=0.35)
legend(x="topleft",legend=c("A","B","C","D"),pt.cex=2,pch=21,
	pt.bg=c("black","yellow","red","blue"))
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)

Our next step is just to fit our different models. We can plot each of
these as we go:


```r
par(mfrow=c(2,2),fg="blue")
er.ordered<-fitpolyMk(poly.tree,y,model="ER",ordered=TRUE)
```

```
## 
## This is the design matrix of the fitted model. Does it make sense?
## 
##     A A+B B B+C C C+D D
## A   0   1 0   0 0   0 0
## A+B 1   0 1   0 0   0 0
## B   0   1 0   1 0   0 0
## B+C 0   0 1   0 1   0 0
## C   0   0 0   1 0   1 0
## C+D 0   0 0   0 1   0 1
## D   0   0 0   0 0   1 0
```

```r
plot(er.ordered,lwd=2,mar=c(1.1,1.1,2.1,1.1))
mtext("a) fitted ER \'ordered\' model",adj=0,
	line=-0.5,col="black")
er.unordered<-fitpolyMk(poly.tree,y,model="ER")
```

```
## 
## This is the design matrix of the fitted model. Does it make sense?
## 
##         A B C D A+B A+C A+D B+C B+D C+D A+B+C A+B+D A+C+D B+C+D A+B+C+D
## A       0 0 0 0   1   1   1   0   0   0     0     0     0     0       0
## B       0 0 0 0   1   0   0   1   1   0     0     0     0     0       0
## C       0 0 0 0   0   1   0   1   0   1     0     0     0     0       0
## D       0 0 0 0   0   0   1   0   1   1     0     0     0     0       0
## A+B     1 1 0 0   0   0   0   0   0   0     1     1     0     0       0
## A+C     1 0 1 0   0   0   0   0   0   0     1     0     1     0       0
## A+D     1 0 0 1   0   0   0   0   0   0     0     1     1     0       0
## B+C     0 1 1 0   0   0   0   0   0   0     1     0     0     1       0
## B+D     0 1 0 1   0   0   0   0   0   0     0     1     0     1       0
## C+D     0 0 1 1   0   0   0   0   0   0     0     0     1     1       0
## A+B+C   0 0 0 0   1   1   0   1   0   0     0     0     0     0       1
## A+B+D   0 0 0 0   1   0   1   0   1   0     0     0     0     0       1
## A+C+D   0 0 0 0   0   1   1   0   0   1     0     0     0     0       1
## B+C+D   0 0 0 0   0   0   0   1   1   1     0     0     0     0       1
## A+B+C+D 0 0 0 0   0   0   0   0   0   0     1     1     1     1       0
```

```r
plot(er.unordered,lwd=2,mar=c(1.1,1.1,2.1,1.1))
mtext("b) fitted ER \'unordered\' model",adj=0,
	line=-0.5,col="black")
transient.ordered<-fitpolyMk(poly.tree,y,model="transient",
	ordered=TRUE)
```

```
## 
## This is the design matrix of the fitted model. Does it make sense?
## 
##     A A+B B B+C C C+D D
## A   0   2 0   0 0   0 0
## A+B 1   0 1   0 0   0 0
## B   0   2 0   2 0   0 0
## B+C 0   0 1   0 1   0 0
## C   0   0 0   2 0   2 0
## C+D 0   0 0   0 1   0 1
## D   0   0 0   0 0   2 0
```

```r
plot(transient.ordered,lwd=2,mar=c(1.1,1.1,2.1,1.1))
mtext("c) fitted transient \'ordered\' model",adj=0,
	line=-0.5,col="black")
transient.unordered<-fitpolyMk(poly.tree,y,
	model="transient")
```

```
## 
## This is the design matrix of the fitted model. Does it make sense?
## 
##         A B C D A+B A+C A+D B+C B+D C+D A+B+C A+B+D A+C+D B+C+D A+B+C+D
## A       0 0 0 0   2   2   2   0   0   0     0     0     0     0       0
## B       0 0 0 0   2   0   0   2   2   0     0     0     0     0       0
## C       0 0 0 0   0   2   0   2   0   2     0     0     0     0       0
## D       0 0 0 0   0   0   2   0   2   2     0     0     0     0       0
## A+B     1 1 0 0   0   0   0   0   0   0     2     2     0     0       0
## A+C     1 0 1 0   0   0   0   0   0   0     2     0     2     0       0
## A+D     1 0 0 1   0   0   0   0   0   0     0     2     2     0       0
## B+C     0 1 1 0   0   0   0   0   0   0     2     0     0     2       0
## B+D     0 1 0 1   0   0   0   0   0   0     0     2     0     2       0
## C+D     0 0 1 1   0   0   0   0   0   0     0     0     2     2       0
## A+B+C   0 0 0 0   1   1   0   1   0   0     0     0     0     0       2
## A+B+D   0 0 0 0   1   0   1   0   1   0     0     0     0     0       2
## A+C+D   0 0 0 0   0   1   1   0   0   1     0     0     0     0       2
## B+C+D   0 0 0 0   0   0   0   1   1   1     0     0     0     0       2
## A+B+C+D 0 0 0 0   0   0   0   0   0   0     1     1     1     1       0
```

```r
plot(transient.unordered,lwd=2,mar=c(1.1,1.1,2.1,1.1))
mtext("d) fitted ER \'unordered\' model",adj=0,
	line=-0.5,col="black")
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20-1.png)

```r
par(fg="black")
```

Compare our different models:


```r
data.frame(transition_model=c("ER","ER","transient","transient"),
	ordered=c("no","yes","no","yes"),
	logLik=c(logLik(er.unordered),logLik(er.ordered),
	logLik(transient.unordered),logLik(transient.ordered)),
	k=c(attr(AIC(er.unordered),"df"),attr(AIC(er.ordered),"df"),
	attr(AIC(transient.unordered),"df"),attr(AIC(transient.ordered),"df")),
	AIC=c(AIC(er.unordered),AIC(er.ordered),
	AIC(transient.unordered),AIC(transient.ordered)))
```

```
##   transition_model ordered    logLik k      AIC
## 1               ER      no -217.7703 1 437.5407
## 2               ER     yes -180.4346 1 362.8692
## 3        transient      no -211.1883 2 426.3765
## 4        transient     yes -180.3397 2 364.6795
```

This tells us that our equal-rates ordered model is the best fit to our
data.

<h3><i>Method 3</i>: Modeling discrete character evolution using Bayesian 
MCMC in R.</h3>

In this final short module we can use the <i>phytools</i> Bayesian MCMC function 
for the extended M<i>k</i> model which is called <code>mcmcMk</code>. Note that
this is also possible to do using the <i>diversitree</i> package.

In this case we'll use the following tree & data:

1. <a href="elopomorph.tre">elopomorph.tre</a>
2. <a href="elopomorph.csv">elopomorph.csv</a>


```r
X<-read.csv("elopomorph.csv",row.names=1)
feed.mode<-setNames(X[,1],rownames(X))
feed.mode
```

```
##                 Albula_vulpes             Anguilla_anguilla 
##                       suction                       suction 
##              Anguilla_bicolor             Anguilla_japonica 
##                       suction                       suction 
##             Anguilla_rostrata                Ariosoma_anago 
##                       suction                       suction 
##           Ariosoma_balearicum           Ariosoma_shiroanago 
##                       suction                       suction 
##        Bathyuroconger_vicinus   Brachysomophis_crocodilinus 
##                          bite                          bite 
##              Conger_japonicus              Conger_myriaster 
##                       suction                       suction 
##               Conger_verreaxi                Conger_wilsoni 
##                       suction                       suction 
##        Congresox_talabonoides            Cynoponticus_ferox 
##                       suction                          bite 
##            Dysomma_anguillare                  Elops_saurus 
##                          bite                       suction 
##         Facciolella_gilbertii          Gavialiceps_taeniola 
##                          bite                          bite 
##         Gnathophis_longicauda          Gorgasia_taiwanensis 
##                       suction                       suction 
##         Gymnothorax_castaneus   Gymnothorax_flavimarginatus 
##                          bite                          bite 
##            Gymnothorax_kidako           Gymnothorax_moringa 
##                          bite                          bite 
## Gymnothorax_pseudothyrsoideus       Gymnothorax_reticularis 
##                          bite                          bite 
##            Heteroconger_hassi          Ichthyapus_ophioneus 
##                       suction                          bite 
##      Kaupichthys_hyoproroides          Kaupichthys_nuchalis 
##                          bite                          bite 
##          Megalops_cyprinoides             Moringua_edwardsi 
##                       suction                          bite 
##             Moringua_javanica              Muraenesox_bagio 
##                          bite                          bite 
##           Muraenesox_cinereus          Myrichthys_breviceps 
##                          bite                       suction 
##          Myrichthys_maculosus         Myrichthys_magnificus 
##                       suction                       suction 
##                Myrophis_vafer        Nemichthys_scolopaceus 
##                          bite                          bite 
##          Nettastoma_melanurum        Ophichthus_serpentinus 
##                          bite                       suction 
##          Ophichthus_zophochir        Oxyconger_leptognathus 
##                       suction                          bite 
## Parabathymyrus_macrophthalmus           Paraconger_notialis 
##                       suction                       suction 
##      Pisodonophis_cancrivorus          Poeciloconger_kapala 
##                          bite                       suction 
##         Rhinomuraena_quaesita          Rhynchoconger_flavus 
##                          bite                       suction 
##        Saurenchelys_fierasfer      Scolecenchelys_breviceps 
##                          bite                       suction 
##            Scuticaria_tigrina             Serrivomer_beanii 
##                          bite                          bite 
##             Serrivomer_sector        Simenchelys_parasitica 
##                          bite                       suction 
##            Uroconger_lepturus      Uropterygius_micropterus 
##                       suction                          bite 
##          Venefica_proboscidea 
##                          bite 
## Levels: bite suction
```

```r
eel.tree<-read.tree("elopomorph.tre")
eel.tree
```

```
## 
## Phylogenetic tree with 61 tips and 60 internal nodes.
## 
## Tip labels:
## 	Moringua_edwardsi, Kaupichthys_nuchalis, Gorgasia_taiwanensis, Heteroconger_hassi, Venefica_proboscidea, Anguilla_rostrata, ...
## 
## Rooted; includes branch lengths.
```

```r
mcmc<-mcmcMk(eel.tree,feed.mode,model="ARD",
	prior.rate=100,prop.var=0.001,ngen=10000)
```

```
## Running MCMC....
## gen 	[suction,bite] 	[bite,suction] 	logLik	accept
## 1	0.001	0.001	-59.2835
## 101	0.0181	0.0101	-37.6057	0.16
## 201	0.0171	0.0239	-38.0419	0.16
## 301	0.0069	0.0138	-37.5138	0.07
## 401	0.0088	0.0193	-38.338	0.13
## 501	0.0191	0.0146	-37.0677	0.15
## 601	0.0109	0.0143	-37.2342	0.07
## 701	0.0245	0.0151	-37.5074	0.18
## 801	0.0112	0.0115	-37.38	0.11
## 901	0.0218	0.0221	-37.3623	0.15
## 1001	0.0138	0.0172	-37.2442	0.19
## 1101	0.0092	0.0099	-37.7779	0.2
## 1201	0.0105	0.0075	-38.6373	0.19
## 1301	0.0024	0.011	-38.0122	0.17
## 1401	0.0083	0.014	-37.4075	0.24
## 1501	0.0065	0.0054	-40.853	0.15
## 1601	0.0043	0.013	-37.7603	0.21
## 1701	0.0142	0.0107	-37.4317	0.09
## 1801	0.0105	0.0037	-40.185	0.1
## 1901	0.0161	0.0131	-37.0964	0.17
## 2001	0.0039	0.0136	-37.8955	0.17
## 2101	0.0198	0.0067	-38.6437	0.11
## 2201	0.0174	0.0223	-37.6736	0.16
## 2301	0.032	0.0275	-37.9694	0.21
## 2401	0.0067	0.0072	-39.1681	0.19
## 2501	0.0027	0.011	-37.9605	0.06
## 2601	0.0118	0.0194	-37.84	0.1
## 2701	0.0208	0.0159	-37.0899	0.1
## 2801	0.0017	0.0129	-38.2357	0.02
## 2901	0.0012	0.0128	-38.3733	0.23
## 3001	0.0244	0.0163	-37.3674	0.18
## 3101	0.0068	0.0187	-38.5906	0.12
## 3201	0.009	0.0097	-37.8575	0.2
## 3301	0.0122	0.0149	-37.173	0.16
## 3401	0.0128	0.016	-37.205	0.19
## 3501	0.0058	0.0138	-37.6411	0.22
## 3601	0.016	0.0093	-37.7009	0.1
## 3701	0.0038	0.0061	-40.2143	0.04
## 3801	0.0176	0.0136	-37.0801	0.14
## 3901	0.0146	0.0166	-37.1257	0.17
## 4001	0.0125	0.0115	-37.3371	0.22
## 4101	0.0157	0.0088	-37.8189	0.16
## 4201	0.0024	0.0113	-37.9892	0.15
## 4301	0.0056	0.0074	-39.0447	0.15
## 4401	0.0454	0.0405	-39.1804	0.07
## 4501	0.0035	0.0118	-37.8197	0.21
## 4601	0.0218	0.0133	-37.4008	0.18
## 4701	0.0196	0.0148	-37.0818	0.24
## 4801	0.022	0.0182	-37.1255	0.13
## 4901	0.0217	0.0158	-37.1533	0.15
## 5001	0.0065	0.0148	-37.6869	0.17
## 5101	0.0134	0.0065	-38.5657	0.13
## 5201	0.0162	0.0152	-37.0135	0.1
## 5301	0.0081	0.0234	-39.967	0.14
## 5401	0.0088	0.0129	-37.3505	0.07
## 5501	0.0033	0.0109	-37.8887	0.11
## 5601	0.0123	0.0187	-37.6293	0.19
## 5701	0.0124	0.0096	-37.7373	0.15
## 5801	0.0226	0.0151	-37.2741	0.17
## 5901	0.013	0.014	-37.1175	0.12
## 6001	0.0213	0.0209	-37.254	0.16
## 6101	0.0068	0.0142	-37.5737	0.12
## 6201	0.0121	0.0114	-37.358	0.12
## 6301	0.0248	0.0216	-37.3459	0.1
## 6401	0.0068	0.009	-38.1268	0.18
## 6501	0.0164	0.021	-37.5607	0.07
## 6601	0.0045	0.0091	-38.1804	0.17
## 6701	0.0015	0.0117	-38.1886	0.16
## 6801	0.014	0.024	-38.6214	0.19
## 6901	0.0129	0.0089	-37.8887	0.11
## 7001	0.0105	0.0145	-37.2697	0.2
## 7101	0.0024	0.0142	-38.3108	0.2
## 7201	0.0073	0.0107	-37.6249	0.1
## 7301	0.0104	0.0079	-38.5012	0.13
## 7401	0.0194	0.0154	-37.0436	0.1
## 7501	0.0121	0.0144	-37.1659	0.15
## 7601	0.0084	0.0103	-37.6833	0.1
## 7701	0.0168	0.0063	-38.4474	0.14
## 7801	0.0135	0.0155	-37.1247	0.18
## 7901	0.0135	0.0179	-37.3525	0.16
## 8001	0.0058	0.0109	-37.6685	0.2
## 8101	0.0195	0.0137	-37.1539	0.07
## 8201	0.0182	0.022	-37.5502	0.15
## 8301	0.0083	0.0122	-37.4027	0.13
## 8401	0.0094	0.0093	-37.9664	0.14
## 8501	0.0064	0.0206	-39.3892	0.24
## 8601	0.0123	0.0234	-38.8129	0.25
## 8701	0.0084	0.0175	-37.9384	0.22
## 8801	0.0342	0.0245	-38.0972	0.14
## 8901	0.0176	0.0087	-37.8985	0.12
## 9001	0.0111	0.0015	-40.6695	0.15
## 9101	0.0088	0.007	-39.1233	0.14
## 9201	0.0049	0.0108	-37.747	0.23
## 9301	0.009	0.0127	-37.343	0.16
## 9401	0.0026	0.0087	-38.4827	0.19
## 9501	0.0221	0.0215	-37.3027	0.19
## 9601	0.0097	0.0137	-37.293	0.18
## 9701	0.0022	0.0092	-38.3521	0.18
## 9801	0.0042	0.0163	-38.4478	0.28
## 9901	0.0068	0.0063	-39.8621	0.12
## Done.
```

Now let's plot our results:


```r
par(mfrow=c(2,2))
mar<-c(4.1,4.1,2.1,1.1)

plot.new()
par(mar=mar)
plot.window(xlim=c(0,1),ylim=c(0,1))
library(png)
download.file(
	"http://www.phytools.org/evol2019/Enchelycore_schismatorhynchus.png",
	"eel-picture.png",mode="wb")
img<-readPNG(source="eel-picture.png")
```

```
## Warning in readPNG(source = "eel-picture.png"): libpng warning: iCCP: known
## incorrect sRGB profile
```

```
## Warning in readPNG(source = "eel-picture.png"): libpng warning: iCCP: cHRM
## chunk does not match sRGB
```

```r
rasterImage(img,0,0,1,1)
mtext(text="a) a biting eel",adj=0,line=0,cex=1)

plotTree(eel.tree,fsize=0.5,ftype="i",
	ylim=c(-8,Ntip(eel.tree)),mar=mar)
FMODE<-to.matrix(feed.mode,levels(feed.mode))
par(fg="transparent")
tiplabels(pie=FMODE[eel.tree$tip.label,],piecol=c("red","blue"),cex=0.3)
par(fg="black")
par(cex=0.8)
add.simmap.legend(colors=setNames(c("red","blue"),c(" bite"," suction")),
	vertical=FALSE,fsize=0.6,prompt=FALSE,x=2,y=-8)
mtext(text="b) phylogeny of eels",
	adj=0,line=0,cex=1)

par(mar=mar)
plot(mcmc,main="")
mtext(text="c) likelihood profile from MCMC",adj=0,
	line=1,cex=1)

d<-density(mcmc)
```

```
## Assuming 20% burn-in as no burn-in was specified....
```

```r
plot(d,main="")
mtext(text=expression(paste(
	"d) estimated posterior density for ",Q[ij])),
	adj=0,line=1,cex=1)
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23-1.png)

That's it.

