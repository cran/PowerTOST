# Author: dlabes
#-----------------------------------------------------------------------------
# a bunch of functions to get the design charcteristics
known.designs <- function()
{
# Nomenclature of replicate designs: treatments x sequences x periods
# Note: the df for the replicate designs are those without carry-over.
# Chen, Chow and Liu used models with carry-over, i.e. 1 df lower.
# n is the total number of subjects in case of cross-over designs 
# and is number of subjects per group in case of parallel group design.
# df2 = degrees of freedom for robust analysis (aka Senn's basic estimator).
# bk is the so-called design constant, here in terms of total n, 
# also without carry over.
# In case of 2x2x4 design Chen, Chow and Liu used bk=1.1 in a model 
# with carry over.
des <- ("
no  design    df      df2    steps  bk    bkni
 0  parallel 2*(n-1) 2*(n-1)  1      2    NA
 1  2x2      n-2     n-2      2      2    0.5
 1  2x2x2    n-2     n-2      2      2    0.5
 2  3x3      2*n-4   n-3      3      2    NA
 3  3x6x3    2*n-4   n-6      6      2    NA
 4  4x4      3*n-6   n-4      4      2    NA
 5  2x2x3    2*n-3   n-2      2      1.5  0.375       # 3/8
 6  2x2x4    3*n-4   n-2      2      1    0.25        # 1/4
 7  2x4x4    3*n-4   n-4      4      1    0.0625      # 1/16
 9  2x3x3    2*n-3   n-3      3      1.5  0.1666667   # 1/6
10  2x4x2    n-2     n-2      4      8    0.5
")
# no. 9 is f.i. the partial replicate design TRR/RTR/RRT
# no. 10 is Balaam's design, a mixture of crossover and parallel group.
#
# eventually it would be better to have steps=6 in case of 3x3 (6 seq. design)
# Jan 2011: 3x6x3 introduced, df for 3x3 corrected (former 2*n-3)
# also df for 4x4 corrected (former 3*n-5)
  
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
#-----------------------------------------------------------------------------
#--- return no of design ---
# design: a character string describing the design
.design.no <- function(design)
{
  #take the first word if more then one f.i. in "parallel group"
  desi <- unlist(strsplit(tolower(design)," "))[1]
	
  des <- known.designs()
  i   <- match(desi, des$design)
  if (!is.na(i)) return(des$no[i]) else return(NA)
}
#--- return all properties as dataframe ---
# or as list?
.design.props <- function(design.no)
{
	des <- known.designs()
  des <- des[des$no==design.no,]
  # in case of 2x2 the alias 2x2x2 causes 2 entries
  des <- des[1,]
  if (is.na(des$no)) stop("Design ",design.no," not defined!")
	return (des)
}	

