# Author: dlabes
#-----------------------------------------------------------------------------
# a bunch of functions to get the design charcteristics
known.designs<-function()
{
# Nomenclature of replicate designs: treatments x sequences x periods
# Note: the df for the replicate designs are those without carry-over
# Chen, Chow and Liu used models with carry-over, i.e. one df lower
# n is the total number in case of cross-over design, 
# n is number per group in case of parallel group design
# df2 = degrees of freedom for robust analysis (aka Senn's basic estimator)
# bk is the so-called design constant (here in terms of total n, 
# also without carry over)
# in case of 2x2x4 design Chen, Chow and Liu used bk=1.1 in a model 
# with carry over
des <- ("
no  design    df      df2    steps  bk
 0  parallel 2*(n-1) 2*(n-1)  1      2
 1  2x2      n-2     n-2      2      2
 2  3x3      2*n-3   n-3      3      2
 3  4x4      3*n-5   n-4      4      2
 4  2x2x3    2*n-3   n-2      2      1.5
 5  2x2x4    3*n-4   n-2      2      1
 6  2x4x4    3*n-4   n-4      4      1
 9  2x3x3    2*n-3   n-3      3      1.5
10  2x4x2    n-2     n-2      4      8    
")
# no. 10 is Balaam's design, a mixture of crossover and parallel group.
# no. 9 is f.i. the partial replicate design TRR/RTR/RRT
#
# eventually it would be better to have steps=6 in case of 3x3 (6 seq. design)
# TODO: the df2 have to be checked, especially for Balaam's design! 
	
  des2 <- textConnection(des)
  designs <- read.table(des2, header=TRUE, sep="", strip.white=TRUE, as.is=TRUE)           
  close(des2)   # without this close() warnings are generated
	
  # nicer names for nicer output of design
  designs$name[designs$no==0] <- "2 parallel groups"
  designs$name[designs$no %in% c(1,2,3)] <- 
		  paste(designs$design[designs$no %in% c(1,2,3)],"crossover")
  designs$name[designs$no %in% c(4,5,6)] <- 
		  paste(designs$design[designs$no %in% c(4,5,6)],"replicate crossover")
  designs$name[designs$no==10] <- "Balaam's design (2x4x2)" 
  designs$name[designs$no==9]  <- "partial replicate design (2x3x3)" 
  #degrees of freedom as expression: not possible, 
  #expression in data.frame not allowed 
  return(designs)
	
}
#--- return no of design ---
# design: a character string describing the design
.design.no<-function(design)
{
  #take the first word if more then one
  desi <- unlist(strsplit(tolower(design)," "))[1]
	
  des <- known.designs()
  i   <- match(desi, des$design)
  if (!is.na(i)) return(des$no[i]) else return(NA)
}
#--- return design constant ---
.design.bk<-function(design.no)
{
  des <- known.designs()
  return(des$bk[des$no==design.no])
}
#--- return design degrees of freedom ---
# as expression
.design.df<-function(design.no, nvar="n")
{
  des <- known.designs()
  df  <- (des$df[des$no==design.no])
  #replace with calling variable name
  df  <- gsub('n',nvar,df)
  e   <- parse(text=df, srcfile=NULL)
  return(as.expression(e))
}
#--- return step for sample size search ---
.design.step<-function(design.no)
{
  des <- known.designs()
  return (des$steps[des$no==design.no])
}
#--- return step for sample size search ---
.design.name<-function(design.no)
{
  des <- known.designs()
  return (des$name[des$no==design.no])
}
#--- return all properties as dataframe ---
# or as list?
.design.props <- function(design.no)
{
	des <- known.designs()
	return (des[des$no==design.no,])
}	

