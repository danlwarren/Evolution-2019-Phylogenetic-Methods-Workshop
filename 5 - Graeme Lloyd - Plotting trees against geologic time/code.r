################################################################################
#                                                                              #
#      PLOTTING TIME-SCALED TREES IN R AND GENERATING DISCRETE CHARACTER       #
#                              PHYLOMORPHOSPACES                               #
#                                                                              #
################################################################################

################################################################################
#                                                                              #
# This script was written by Graeme T. Lloyd for the Evolution workshop held   #
# in Providence, Rhode Island (21/6/19)                                        #
#                                                                              #
################################################################################

################################################################################
#                                                                              #
#                           HOW TO USE THIS SCRIPT:                            #
#                                                                              #
# All lines below are either comments (beginning with the hash symbol) or R    #
# code. Both can be copied and pasted directly into the R console, but the     #
# former will not be executed.                                                 #
#                                                                              #
# To use this script (which assumes no prior R knowledge!) you must:           #
#                                                                              #
# 1. (If you have not already) install R from: http://cran.r-project.org/      #
# 2. Open R.                                                                   #
#                                                                              #
# Then:                                                                        #
#                                                                              #
# 3. Work through each line below, reading the explanations beginning with '#' #
#    and copying and pasting the rest straight into R. (You may have to hit    #
#    enter to execute the line.) Note that many operations in R will provide   #
#    no immediate feedback on their effect, but persevere as the next line     #
#    will usually explain things.                                              #
#                                                                              #
################################################################################

################################################################################
#                                                                              #
# SECTION I - INSTALLING PACKAGES                                              #
#                                                                              #
# (NB: You can skip this once you have installed packages - it only needs to   #
# be done once.)                                                               #
#                                                                              #
################################################################################

# We begin by making a "vector" (think a single row) of R packages that will
# need to be installed. This is long primarily because of things like
# hierarchical dependencies between packages, but also because we will be using
# functions ("doing" code) from across a variety of packages:
LloydSectionPackages <- c("ade4", "ape", "devtools", "gdata", "paleotree",
  "phytools", "msm", "strap", "rgl")

# Next we will install the packages and their dependencies. Note that this
# will typically take a while to run. (Also, a top tip for Mac users installing
# packages from the menu system instead: the dependencies there is a checkbox
# that for some reason defaults to unchecked - make sure you check it before
# installation.):
install.packages(LloydSectionPackages, dependencies = TRUE)

# Note that the above installs packages from CRAN, but sometimes it is better
# to install from github when packages are available from there as this will
# more typically get you the absolute latest version. We can still do this
# pretty easily from inside R using the devtools package we just installed.
# We will do this for Claddis with:
library(devtools)
install_github("graemetlloyd/Claddis")

# You should now have all the packages you need to run this script and you can
# ignore the above block in future.

################################################################################
#                                                                              #
# SECTION II - LOADING PACKAGES                                                #
#                                                                              #
# (NB: You will need to do this at the beginning of each R "session" as this   #
# effectively "loads" the functions you will need into memorey ready for       #
# access. Without doing this only R's "base" functions are accessible.)        #
#                                                                              #
################################################################################

# Load the various packages into memory:
library(ade4)
library(paleotree)
library(phytools)
library(strap)
library(gdata)
library(Claddis)

# For another top tip you can access the help file for a package or function
# by typing its' name directly after a question mark:
?Claddis

# If you scroll to the bottom you will also find an "index" link that will take
# you to a list of all the functions in that package. This is a good way to
# learn all the things you can use the package for (assuming the author(s)
# have done a good job with their documentation that is!).

################################################################################
#                                                                              #
# SECTION III - GETTING DATA INTO R AND PLOTTING TIME-SCALED TREES             #
#                                                                              #
# (NB: The former is often the hardest thing about using R.)                   #
#                                                                              #
################################################################################

# Before we can use any of the functions we now have in memory we need some
# data to apply them to. As this is just a workshop we will use an example data
# set chosen not for its' interest, but because it is small. This will help
# with keeping analysis run times fast, but also make it easier to process what
# we are looking at when we "print" variables to the screen (R console).
#
# Specifically, here we will use a dinosaur data set from Cullen et al (2013):
browseURL(url =
  "http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0058853")

# This is an open access paper so you should be able to read the full version
# inside your browser (where the above line of code should send you). You don't
# need to know anything about dinosaurs for this workshop though!
#
# For the practical exercises below we will specifically need to import three
# different pieces of data into R's memory. These are:
#
# 1. A discrete character-taxon matrix (summarising the morphological variation
#    in our dinosaurs).
# 2. A phylogenetic hypothesis for our dinosaur species.
# 3. The geologic ages of our dinosaur species.
#
# Note that any other data we might need we can generate from just these three
# starting points using R functions along the way.
#
# The first piece of data is accessible from my web site here:
browseURL("http://www.graemetlloyd.com/nexus/Cullen_etal_2013a.nex")

# If you like you can inspect the file for its basic structure (#NEXUS format)
# then download it to your personal computer, e.g., by right-clocking on it and
# choosing "Save (As)" or similar (e.g., File > Save (As)).
#
# We will also want a tree to play with. For now we will use just a single tree
# (the first most parsimonious tree from a larger set). (However, it is worth
# noting that good practice is to incorporate phylogenetic uncertainty by
# applying analyses across a sample of trees.) This is also available
# from my web site:
browseURL("http://www.graemetlloyd.com/firstmpt/Cullen_etal_2013a.tre")

# Again, briefly inspect the structure. This tree is in "Newick" format. Don't
# worry if you have your own tree in the alternative #NEXUS format. That can
# still be read into R. For now though, just download this file to your hard
# drive in the same way as the matrix.
#
# The final piece of data we will need for this example is the ages of our taxa.
# One thing we will use this for is to time-scale our tree, but even if you have
# a tree that is already time-scaled the ages file is important to have as
# it can be used for other purposes as we will find out later. Again, this is
# available from my web site (from an older workshop) at:
browseURL("http://www.graemetlloyd.com/teaching/SVP2017/Cullenages.txt")

# Again, take a while to inspect the data format. Although this is not a
# standard format in the same way as the tree and matrix it is the way a lot of
# R packages like age data to be in (including paleotree and strap). Again,
# download this data to your hard drive.
#
# Now what is often the trickiest part with R: getting the data in. We will
# start with the matrix. Run the line below and navigate to the folder you
# downloaded your data to and select the matrix (.nex) file:
nexus.data <- ReadMorphNexus(File = file.choose())

# In the above line we have stored (<-) the data in a variable (nexus.data).
# This data is stored in an R format known as a "list", which is one way to
# store different types of data together. The components of this list have
# names that can be seen by typing:
names(nexus.data)

# We can view individual components of this list by using the dollar symbol and
# a name from above, e.g.:
nexus.data$Topper

# This is simply the header text from the file (if there is any) and any step
# matrices used (if any). NB: Both will default to NULL if they are missing.
#
# Now try:
nexus.data$Matrix_1

# This shows the character-taxon data along with a lot of other information.
# In fact this is another list:
names(nexus.data$Matrix_1)

# (Note also that this is Matrix_1 because it is possible to have different
# matrix "blocks". These could be, for example, morphology, sequence data,
# different genes, continuous characters etc. Here there is only one as the
# characters are all discrete and morphological, but if there were more they
# would be numbered Matrix_2, Matrix_3 etc.)
#
# We won't go through all of these data types, but broadly they are:
#
# 1. BlockName - Fairly self explanatory.
# 2. Datatype - Ditto.
# 3. Matrix - The actual character matrix with taxa as row names and characters
#    as columns (these would represent loci if sequence data).
# 4. Ordering - Whether characters should be treated as "ord"(ered), i.e.,
#    Wagner optimised, or "unord"(ered), i.e., Fitch optimised. Other possible
#    values here are "step_A", "step_B" etc. which would correspond to step
#    matrices (which must be specified in the "Topper" part (above) or
#    "cont"(inuous) characters. Note that currently Claddis cannot deal with
#    morphometric or Dollo characters, although I may add these in future.
# 5. Weights - The weights to be applied to each character (column).
# 6. MinVals - The "minimum" value for each character. NB: This is only
#    really relevant if using Wagner optimised (ordered) characters.
# 7. MaxVals - The "maximum" value for each character. Same rules as above.
# 8. Characters - This is information for writing out the data. E.g., if the
#    input data are loci (coded as A, C, G, and T) then they will appear in
#    the Matrix as 0, 1, 2, and 3 respectively, but the data will still be
#    written to a file (e.g., using WriteMorphNexus) correctly using the
#    data stored here.
#
# Let's just take a quick look at the matrix itself:
nexus.data$Matrix_1$Matrix

# You should see a mixture of codings here. "0", "1" and "2" correspond to
# specific morphological codings and NA to missing data. (Note that Claddis
# can also handle "inapplicable" characters and when present they are coded
# as empty strings, i.e., "".) There is another type of coding too:
nexus.data$Matrix_1$Matrix["Sinornithomimus_dongi", 22]

# This is a polymorphism (where two or more morphologies are observable in
# the same taxon, in this case Sinotnithomimus). These are very common in
# morphology (and also exist in some sequence data). Claddis also allows
# for a different kind of ambiguous coding known as uncertainties. In
# #NEXUS format these are coded with multiple states inside curly braces,
# e.g., {01} and in Claddis are shown with a slash instead of an ampersand,
# i.e., "0/1".
#
# Polymorphisms and uncertainties tend to complicate analysis of matrices,
# but for now you just need to know that they exist.
#
# Now we can look at ordering:
nexus.data$Matrix_1$Ordering

# In this case everything is unordered (indicated by "unord"). In other
# words even if a character has three states (0, 1, and 2) then this tells
# Claddis that the transitions 0 -> 2 and 2 -> 0 are allowed, i.e., there is
# no requirement to pass through state 1 to get there.
#
# Now let's look at weights:
nexus.data$Matrix_1$Weights

# You should see that at present everything is weighted one. Because all our
# characters are unordered (or binary) this means they are equally weighted.
# We can again edit our weights if we so chose. For example, lets say we wish to
# exclude character 6 as we are unhappy with it for some reason. We could delete
# it, but that would get a little messy as we would have to remove it from the
# matrix, the ordering list, the weights list etc. A much simpler way that has
# the same mathematical outcome is to weight it zero.
nexus.data$Matrix_1$Weights[6] <- 0

# Note that if you really do need to delete characters or taxa then you should
# do so using the MatrixPruner function. We will not use it today, but you can
# see the help file for it with:
?MatrixPruner

# We can see our weighting change has worked with:
nexus.data$Matrix_1$Weights[6]

# Now let's look at the maximum values for all of our characters:
nexus.data$Matrix_1$MaxVals

# We can see some twos in that list (implying multistate characters), but we
# need to see the minimum values to be sure:
nexus.data$Matrix_1$MinVals

# In this example they are all zeroes, but note that they don't have to be. (A
# binary character could be coded with states 2 and 3 for example.) Note also
# that we shouldn't alter these values - they exist to help the disparity
# functions such as GED and MORD to calculate maximum possible distances.
#
# That covers the basic structure of the data matrix. Note that if you prefer to
# build your matrices in something like Excel you can still import them into R
# (after saving it as a plain-text file) by using the MakeMorphMatrix function.
# If this is your preference then it is worth looking at the help file for it:
?MakeMorphMatrix

# Now lets move onto the tree. We can read that in in a similar way to the
# matrix, but this time we will use a function from the ape package (read.tree).
# Again, navigate to the folder you downloaded your tree file (.tre) to and load
# it in with:
tree.data <- read.tree(file = file.choose())

# (If your tree data is in #NEXUS format then you can substitue read.tree for
# read.nexus.) This time we have stored the input as "tree.data", and we can get
# a brief summary with:
tree.data

# You should see some key information here. There are 15 tips and 12 internal
# nodes. This should immediately flag up that there are polytomies as for a
# fully bifurcating tree there should be N - 1 nodes (14 in this case). There is
# also a partial list of our tip names (labels) and two more key pieces of
# information: 1) our tree is rooted, and 2) there are no branch lengths.
# Rooting is important for palaeontological trees as we need a "starting point".
# You can root an unrooted tree with the ape command "root". Branch lengths are
# also important as we ultimately want a tree where branches reflect time.
# Collectively this summary indicate our tree is not ready for analysis so we
# will have to perform some additional steps first (more below). For now though
# lets look a bit deeper into how ape stores a tree. Like our data matrix the
# data are actually a list as they mix different types of data (numbers and
# text). We can see that again with:
names(tree.data)

# There should only be three items, but note that the more complex the tree the
# more items there may be. But our basics are listed and we can look at them
# separately. First "edge":
tree.data$edge

# This is a matrix that you can think of as having columns of "from" and "to"
# values that describe each branch. The numbers here are somewhat confusing
# though as they are ape specific and don't necessarily correspond to tip or
# node labels you might see in other software. The general rule for labels is
# as follows: tip are numbered from 1 to N (here 1 to 15), the root is numbered
# N + 1 (16 here), and the remaining nodes are numbered such that each internal
# branch goes from a low number to a high number (e.g., 19 to 20). We can make
# things much clearer with a simple plot:
plot(tree.data, show.tip.label = FALSE)
tiplabels()
nodelabels()

# Here we have plotted the tree (without our actual taxon names) and overlaid the
# tip and node numbers. Now you should be able to match up branches from $edge
# to the plot. We can also now see that our tree has two polytomies
# (multifurcations), at nodes 18 and 27. This is less obvious from $edge, but is
# noticeable as we will see 18 and 27 appear three times in the first ("from")
# column.
#
# Next up let's look at our tip labels:
tree.data$tip.label

# Unsurprisingly these are our taxon names. Note that here the genus and species
# name is separated by an underscore. I *strongly* recommend you use this format
# with your own taxon names as any other characters can cause different kinds of
# problems. It's worth stopping at this point to check these names match those
# in our data matrix as this is *by far* the most common problem I encounter
# from users. Here we can use a base R function called "setdiff" that takes two
# lists of names or numbers and returns which are in the first, but not the
# latter. As an example we can use the letters A-D and the letters A-C:
setdiff(c("A", "B", "C", "D"), c("A", "B", "C"))

# This should tell us that D is in the first list, but not the latter. If we add
# D to the latter (fixing the problem) then we get:
setdiff(c("A", "B", "C", "D"), c("A", "B", "C", "D"))

# This returns the rather odd message "character(0)". This is R's way of telling
# you that the answer is an empty vector of length zero and type "character". In
# other words there is an empty text field because all the items in the first
# list of names are also in the second list. Importantly this is an asymmetric
# question so in reality we will want to run it both ways (i.e., change the
# order we give things to setdiff). Now let's give our actual data to it. We
# already know how to get the taxon names from our tree (tree.data$tip.label),
# but what about our matrix? The answer is rownames(nexus.data$Matrix_1$Matrix).
# Note that rownames(nexus.data) will not work as we first have to identify the
# matrix part of our data. So now we can ask whether our taxon names match by
# running setdiff both ways:
setdiff(tree.data$tip.label, rownames(nexus.data$Matrix_1$Matrix))
setdiff(rownames(nexus.data$Matrix_1$Matrix), tree.data$tip.label)

# Phew. We get character(0) both times so we know everything matches up. For the
# sake of the exercise though lets break things to see how the answer changes.
# Lets add a typo to our tree, specifically we will change the name of our first
# taxon (the outgroup) Allosaurus:
tree.data$tip.label[1] <- "Alosaurus_fragilis"

# And re-run our setdiff pair:
setdiff(tree.data$tip.label, rownames(nexus.data$Matrix_1$Matrix))
setdiff(rownames(nexus.data$Matrix_1$Matrix), tree.data$tip.label)

# Now we should see that we get *two* wrong answers. An even set of names
# normally indicates a typo and you should visually be able to match them up
# into pairs. If you get an odd number or the pairs don't match it is likely
# that you simply have more taxa in your tree or matrix than your matrix or
# tree, e.g., because you deleted a taxon from one but not the other. There
# are different ways to deal with this, but remember that caution note from
# above that you should be careful about deleting taxa from the matrix *inside*
# R (use MatrixPruner to be safe). For now let's fix our mistake so we can move
# on:
tree.data$tip.label[1] <- "Allosaurus_fragilis"

# One final note on the tip labels is that the order reflects the tip numbers
# from $edge. I.e., the first entry is node 1, the second is node 2 and so on.
# The final object in our tree.data is the node number:
tree.data$Nnode

# This is just the number of (internal) nodes in our tree.
#
# At this point we should deal with one of the issues with our tree,
# specifically the polytomies. (Note that it is important that we deal with
# polytomies *before* time-scaling as otherwise we will end up with zero-length
# branches that will cause problems later.) As a general rule you should not use
# a consensus tree for analyses (as it will often be a suboptimal tree), but in
# this case we are using an optimal (most parsimonious) tree. There are several
# ways we might want to resolve such a tree, but for now we are just going to
# resolve it randomly using ape's multi2di function:
tree.data <- multi2di(tree.data)

# Note that here we have overwritten tree.data which is a bad habit to get into,
# but we can confirm that there are now 14 internal nodes:
tree.data$Nnode

# We can also visually check by plotting the tree:
plot(tree.data)

# Note that this has also changed our node numbers as there are additional
# branches.
#
# The next step is to time-scale our tree, but first we need to import our age
# data. This time we will use the base R function read.table:
ages.data <- read.table(file = file.choose(), row.names = 1, sep = "\t",
  header = TRUE)

# To get an idea what format the data is in let's look at it:
ages.data

# This time it is not a list, but a simple table with two columns (FAD = First
# Appearance Datum, LAD = Last Appearance Datum) and a row for each taxon.
# Before we time-scale the tree we should repeat our taxon name check from
# before to ensure all our names match up with both our tree and matrix:
setdiff(rownames(ages.data), rownames(nexus.data$Matrix_1$Matrix))
setdiff(rownames(ages.data), tree.data$tip.label)
setdiff(tree.data$tip.label, rownames(ages.data))
setdiff(rownames(nexus.data$Matrix_1$Matrix), rownames(ages.data))

# You should see lots of "character(0)" which is what we want to see.
#
# We're now ready to produce a time-scaled tree. Note that this is essential for
# pretty much any form of ancestral state estimation, i.e., what we will
# ultimately be doing to build phylomorphospaces below. In R we have a couple of
# useful options here: the timePaleoPhy function in the paleotree package and the
# DatePhylo function in the strap package. You could also date your tree
# "outside" R and import it with the branch-lengths alredy included, but note
# that you will also need to state the "root age" (see below) to use the plotting
# function. We will create our first time tree with:
time.tree1 <- timePaleoPhy(tree.data, ages.data, type = "basic")

# Here we are providing our tree, our ages for our tips, and a choice of method
# (basic). This method is the "traditional" way in which palaeontologists have
# time-scaled their trees: each node is considered to be as old as its oldest
# descendant. This is the way trees have often been conceptualised to work,
# but in reality these are absolute minimum estimates (assuming no ancestor-
# descendant pairs exist). This causes a major problem for us that can be
# shown with a quick visualisation:
plot(time.tree1)

# What should be immediately obvious is that there are now lots of polytomies.
# This isn't quite right though, as we already know there are no polytomies in
# our data. We can see what is going on by checking the names for our tree:
names(time.tree1)

# You should see that there are now two more items here - edge.length and
# root.time. Let's look at edge.length as this will reveal the issue:
time.tree1$edge.length

# You should see a list of numbers. These are the branch lengths for our tree
# and correspond (in order) to the rows in $edge. The apparent polytomies are
# caused by lots of zero-length branches (or ZLBs), i.e., when we plot these
# branches the x-values for each end are the same. Not only is this unsatisfying
# as a visualisation of our tree this makes many mathematical operations
# impossible. E.g., ancestral state estimation cannot be calculated as we
# would run into a divide by zero problem. Just to be clear this is not caused
# by a lack of precision in our fossil ages. Even if every tip had a unique age
# this problem would still occur. In fact, at a minimum, half of the branches in
# a fully bifurcating tree time-scaled this way will be zero. We can see how
# many are in this case by doing a simple calculation:
sum(time.tree1$edge.length == 0) / length(time.tree1$edge.length) * 100

# This should confirm that c. 60% of our branches have zero-length (NB: the exact
# answer may vary because of the random bifurcating option earlier.) Before we
# move on to methods that eliminate ZLBs one final point worth making is that
# ZLBs can be legitimate outcomes in some scenarios as they would indicate an
# ancestral relationship (see the cal3 function in paleotree). But here we can
# be fairly confident that this is not the case.
#
# It's worth pointing out here that it is better to get rid of ZLBs using a
# probablistic approach (e.g., the fossilised birth-death model), but that
# would involve a whole workshop on its' own. So for now we will stick with
# the available (and fast), albeit somewhat arbitrary approaches available
# inside R. As a concession to practicallity we can at least try multiple
# options here to see how sensitive the outcome is to this choice. Here we will
# try just two of the options available to us - the "equal" and "mbl" aproaches:
time.tree2 <- timePaleoPhy(tree.data, ages.data, type = "equal", vartime = 2)
time.tree3 <- timePaleoPhy(tree.data, ages.data, type = "mbl", vartime = 2)

# The first of these (equal) works by first time-scaling the tree using the
# "basic" method then extending each ZLB by allowing it an equal share of the
# time on the first preceding branch of positive length. This is an iterative
# approach that starts with the most deeply nested ZLB and continues until all
# ZLBs are removed. Note that we run into a problem at the root of the tree
# where there is no preceding branch. This is what the "vartime" value is used
# for. This is the value in millions of years that is assigned to the root
# "branch", i.e., the imaginary branch leading from the root backwards in time.
# There is no perfect way to chose this value - one reason why such an
# algorithm can be considered arbitrary. We can now check edge.length to see if
# we have gotten rid of our ZLBs:
time.tree2$edge.length

# You should see all numbers are positive. But we can double check with another
# simple calculation:
all(time.tree2$edge.length > 0)

# For our "mbl" tree we have used a minimum-branch length that we force to be
# positive and so avoid ZLBs. Here vartime represents the minimum length, and
# again there is no perfect way to chose this. Again we will go with 2 million
# years for this value and can see that it has also removed any ZLBs:
time.tree3$edge.length
all(time.tree3$edge.length > 0)

# We are now ready to visualise our time-scaled trees properly by plotting
# them against a geologic timescale. To do this we will use the geoscalePhylo
# function from the strap package:
geoscalePhylo(ladderize(time.tree2), cex.age = 0.6, cex.ts = 0.8, cex.tip = 1,
  x.lim = c(166, 23.03))

# Here we are handing the function one of our time-scaled trees then a bunch
# of other values that mostly decide font sizes (cex = R shorthand for character
# expansion). It is worth playing around with these to see how they affect the
# plot. The other option (x.lim) sets the maximum and minimum values (i.e.,
# lim(its)) of the x-axis. These are also often worth setting manually as they
# will typically want to extend outside of the plotting range of the tree so
# that the tip names are visible. I.e., here we extend into the Cenozoic even
# though non (non-avian) dinosaurs do.
#
# We are also seeing the "full" geologic time-scale here, with geologic
# Periods (Cretaceous, Jurassic etc.), Epochs (Upper Jurassic, Lower Cretaceous
# etc.), and Stages (Campanian, Maastrichtian etc.). We can be more selective
# about this with the "units" argument, for example getting just the stages with:
geoscalePhylo(ladderize(time.tree2), cex.age = 0.6, cex.ts = 0.8, cex.tip = 1,
  x.lim = c(166, 23.03), units = "Age")

# Or just Periods and Epochs with:
geoscalePhylo(ladderize(time.tree2), cex.age = 0.6, cex.ts = 0.8, cex.tip = 1,
  x.lim = c(166, 23.03), units = c("Period", "Epoch"))

# We can also add some background "stripes" to help reading up and down the
# plot with:
geoscalePhylo(ladderize(time.tree2), cex.age = 0.6, cex.ts = 0.8, cex.tip = 1,
  x.lim = c(166, 23.03), boxes = "Age")

# As mentioned earlier, it is worth importing age data even if your tree is
# already time-scaled as we can also give these ages to the plot so the ranges
# of the tips can be plotted as thick black lines:
geoscalePhylo(ladderize(time.tree2), cex.age = 0.6, cex.ts = 0.8, cex.tip = 1,
  x.lim = c(166, 23.03), boxes = "Age", ages = ages.data)

# There are many other options with the function you can explore yourself with:
?geoscalePhylo

# However, one important point to visit before we move on is that other item
# that appears in the names of our time-scaled trees:
names(time.tree2)

# I.e., the "root.time" value. This is, as you might have guessed, the age
# of the root of our tree. This is important for palaeontological trees, but is
# not yet incorporated in standard (Newick or #NEXUS) tree formats. Thus if you
# import a tree time-scaled elsewhere into R you will need to manually set a
# root.time value before you can do things like plot it against stratigraphy.
# Let's check our root.time values for our three time-scaled trees:
time.tree1$root.time
time.tree2$root.time
time.tree3$root.time

# Note that our basic tree has the youngest age (as we might expect) and the
# other two are exactly 2 million years older (due to the vartime value we used).
#
# We now have the three data items we need in a format we can use them so we're
# almost done for this part. However, let's do one more thing that is helpful to
# learn: how to extract graphics from R in vector format. Specifically, we are
# going to make a PDF file with three pages that correspond to our three trees
# plotted against stratigraphy. We can do this very simply with:
pdf("Timetrees.pdf", width = 10, height = 8)
geoscalePhylo(ladderize(time.tree1), cex.age = 0.6, cex.ts = 0.8, cex.tip = 1,
  x.lim = c(166, 23.03), boxes = "Age", ages = ages.data)
geoscalePhylo(ladderize(time.tree2), cex.age = 0.6, cex.ts = 0.8, cex.tip = 1,
  x.lim = c(166, 23.03), boxes = "Age", ages = ages.data)
geoscalePhylo(ladderize(time.tree3), cex.age = 0.6, cex.ts = 0.8, cex.tip = 1,
  x.lim = c(166, 23.03), boxes = "Age", ages = ages.data)
dev.off()

# The first line (pdf) indicates we are writing to a PDF, the name in quotes is
# the file name we want to create (NB: this will overwrite without warning if
# you run it again), the width is the width of the plots (in inches), and the
# height is the height of the plot in inches. We then run our three plot lines
# (note that you won't see them anymore as they are going to the PDF file and
# not the screen). Then to indicate we are done writing to the PDF we close our
# graphics device with dev.off().
#
# You might wonder where this file is because we didn't specify a specific path
# when we created it. The answer is wherever this line points you:
getwd()

# (This is short for "get working directory".) You can change this to a
# different folder using the function setwd() or by using the menu system within
# R (look for the set or change working directory option - it is different
# depending on your OS).
#
# Before we finish try changing and re-running the pdf creation lines from
# above. In particular, try changing width and height for the pdf (e.g., "width
# = 20, height = 5") and different numbers for the various cex values. This will
# help you work out which values you might prefer when dealing with your own
# data, although note that trees of different sizes might be better suited to
# different values (i.e., if you use your own data later you might want to change
# these again for that specific case, e.g., when I do this for very large trees
# I tend to make the PDF size enormous).

################################################################################
#                                                                              #
# SECTION IV - DISCRETE CHARACTER PHYLOMORPHOSPACES                            #
#                                                                              #
################################################################################

# As we shift our focus to phylomorphospaces we are also going to shift our
# focus to the Claddis package (the name is a contraction of cladistic disparity
# as the package makes use of the same cladistic matrix formatting seen in
# #NEXUS files and the "disparity" part refers to a palaeontological shorthand
# for the general measurement of morphological diversity). We already used one
# function from this package - ReadMorphNexus. This gets the required discrete
# character data into R, but we need to produce a lot more data before we can
# build a phylomorphospace. Primarily these are:
#
# 1. Some form of ordination - this will give us a set of coordinates that
#    represent where our taxa fall in the "morphospace".
# 2. Some form of ancestral state estimation - this will set the hypothetical
#    morphology of our ancestors (branching points) in the tree.
#
# Critically we can perform these two steps in either order, but the results
# are typically very different (as we saw in the lecture part). We can
# generate the pre-ordination ancestral state estimate data with:
preOASE1 <- MorphMatrix2PCoA(nexus.data, Tree = time.tree2,
  EstimateAllNodes = FALSE, EstimateTipValues = FALSE)
preOASE2 <- MorphMatrix2PCoA(nexus.data, Tree = time.tree2,
  EstimateAllNodes = FALSE, EstimateTipValues = TRUE)
preOASE3 <- MorphMatrix2PCoA(nexus.data, Tree = time.tree2,
  EstimateAllNodes = TRUE, EstimateTipValues = FALSE)
preOASE4 <- MorphMatrix2PCoA(nexus.data, Tree = time.tree2,
  EstimateAllNodes = TRUE, EstimateTipValues = TRUE)

# These should take a little while to run (the bottleneck is the ancestral
# state estimates), but this function (MorphMatrix2PCoA) is a fast
# "wrapper" for a bunch of functions inside Claddis. Specifically, it is
# working through the following steps:
#
# 1. Estimate discrete ancestral state for each character using the
#    supplied tree (here we will just use time.tree2 from before). It will
#    not perform this step if no tree is provided (this is what happens below
#    for the post-ordination approach). This uses the Claddis function
#    AncStateEstMatrix.
# 2. Create a pairwise distance matrix. This uses the Claddis function
#    MorphDistMatrix.
# 3. Checks this matrix for any missing (incalculable) values and trims
#    taxa/nodes until the matrix is complete. This uses the Claddis function
#    TrimMorphDistMatrix.
# 4. Performs an ordination using principal coordinates. This uses the
#    ape function pcoa.
#
# There are a *lot* of options here we are skipping over, see the help
# file for more information:
?MorphMatrix2PCoA

# This will also refer you to the help files for all the other
# functions called along the way, but for now we will stick to the
# defaults beyond what we specified above.
#
# The post-ordination ancestral states are a little more fiddly. First
# we will perform the ordination using the same function as above, but
# this time without the tree:
postOASE <- MorphMatrix2PCoA(nexus.data)

# We can now check if any taxa have been removed (because they lead to
# incomplete distances):
postOASE$RemovedTaxa

# You should see Pelecanimimus polyodon had to be removed. Note that
# this is a common occurrence with palaeontological data, but should
# be less of an issue with sequence data. Because this taxon is being
# pruned we also need to prune the taxon from our tree before we
# produce ancestral state estimates. We can do that with ape's drop.tip
# function:
time.tree2.trimmed <- drop.tip(time.tree2, postOASE$RemovedTaxa)

# Note that whenever you prune taxa from a time-scaled tree in R you
# need to be careful to also check this hasn't altered the "root.time"
# we encountered earlier. You can do this using either the fixRootTime
# function in paleotree, but there is a similar function in Claddis
# called CorrectRootTime that will do the same thing and uses the full
# and trimmed trees:
time.tree2.trimmed <- CorrectRootTime(time.tree2, time.tree2.trimmed)

# Note that in this case there was no need to do this (the root.time
# is still 159). However, it is worth adopting good habits like this
# as things could go badly wrong in future if this step is omitted.
#
# We can now store this tree in the correct spot in our postOASE
# variable with:
postOASE$Tree <- time.tree2.trimmed

# (This step means the plotting functions we will use below will
# automatically produce a phylomorphospace. Without this the $Tree variable
# will be NULL and the plot will default to a regular morphospace.)
#
# Now we need to generate *continuous* ancestral state estimates using
# our ordination axes, which you can see with:
postOASE$vectors

# The next line will generate these for each axis, but this gets a
# little more complex than what we have done so far, using more
# "advanced" R functions (apply, lapply, do.call). All you need to
# know is that this is estimating ancestral states for each axis:
PostOrdinationACEs <- do.call(cbind, lapply(apply(postOASE$vectors, 2,
  list), function(x) ace(x[[1]], phy = time.tree2.trimmed)$ace))

# We can now incorporate these into our postOASE variable in the
# same manner as the preOASE approaches do:
postOASE$vectors <- rbind(postOASE$vectors, PostOrdinationACEs)

# We can now visually confirm that our choices lead to different
# phylomorphospaces by plotting them out using the MorphospacePlot
# function:
par(mfrow = c(2, 3))
MorphospacePlot(preOASE1, plot_taxon_names = TRUE)
title("preOASE1")
MorphospacePlot(preOASE2, plot_taxon_names = TRUE)
title("preOASE2")
MorphospacePlot(postOASE, plot_taxon_names = TRUE)
title("postOASE")
MorphospacePlot(preOASE3, plot_taxon_names = TRUE)
title("preOASE3")
MorphospacePlot(preOASE4, plot_taxon_names = TRUE)
title("preOASE4")

# Typically you should be cautious about overinterpreting data from
# what is really a multidimensional space that you are only seeing
# two axes of here. However, we can attempt to visualise a third
# by using the size of the points being plotted to capture the
# absolute value of a third (z-)axis and the colour (black/white) to
# capture the sign (positive/negative) with:
par(mfrow = c(2, 3))
MorphospacePlot(preOASE1, x_axis = 1, y_axis = 2, z_axis = 3,
  plot_taxon_names = TRUE)
title("preOASE1")
MorphospacePlot(preOASE2, x_axis = 1, y_axis = 2, z_axis = 3,
  plot_taxon_names = TRUE)
title("preOASE2")
MorphospacePlot(postOASE, x_axis = 1, y_axis = 2, z_axis = 3,
  plot_taxon_names = TRUE)
title("postOASE")
MorphospacePlot(preOASE3, x_axis = 1, y_axis = 2, z_axis = 3,
  plot_taxon_names = TRUE)
title("preOASE3")
MorphospacePlot(preOASE4, x_axis = 1, y_axis = 2, z_axis = 3,
  plot_taxon_names = TRUE)
title("preOASE4")

# Alternatively, for a single data set (here preOASE1) you can
# simply produce multiple bivariate plots with the
# MultiMorphospacePlot function:
par(mfrow = c(1, 1))
MultiMorphospacePlot(preOASE1, N_axes = 4, plot_taxon_names = TRUE)

# Finally, if we want to go really fancy we can produce a three-
# dimensional "chronophylomorphospace" where two axes come from our
# morphospace and the third is simply geologic time. This can be done
# with the ChronoPhyloMorphospacePlot function. Note for this to work
# you need to have successfully installed the rgl package which can be
# a bit finicky. However, hopefully you managed this OK and so you can
# generate the plot for your preferred phylomorphospace (here I will
# use the postOASE one) with:
ChronoPhyloMorphospacePlot(postOASE)

# Note that this is an "interactive" plot and you should be able to
# "grab" it with your cursor and rotate it to get a better sense of what
# it shows.
#
# You have now reached the end of the script! Hopefully you were able
# to confirm for yourself that our view of evolution from a
# phylomorphospace can dramatically vary without changing the data,
# only the methods we use to produce it. I discuss which of these might be
# best (hint: preOASE1) in my 2018 paper in Palaeontology.
