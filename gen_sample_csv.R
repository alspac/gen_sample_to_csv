### Sam Neaves
### Sept 2019
### Convert one snp in .gen and .sample file to a csv

library(dplyr)
library(data.table)

### This should be renamed. This function is for making the three repeats for index
my_function <- function(incol)
{
    n <- 3
    out <- rep(1:length(incol), each=n)
    return(out)
}

best_guess_to_allele <- function(bestguess,A1,A2)
{
    if(bestguess == "1,0,0")
    {
        return(paste(A1,A1,collapse=NULL,sep=""))
    }
    else if(bestguess=="0,1,0")
    {
        return(paste(A1,A2,collapse=NULL,sep=""))
    }
    else if(bestguess=="0,0,1")
    {
        return(paste(A2,A2, collapse=NULL, sep=""))
    }
    else
    {
        return("error")
    }
}




###To read the sample file:

data_samples <- read.table("data.sample", header=TRUE)

###To read the gene file:

data_gene <- read.table("subsetted.gen")

####Fix that r has read "T" as "Truth":


data_gene$V5[1] <- "T"


###Label each group of 3 with an index
transposed_label <- cbind(t(data_gene)[c(-1,-2,-3,-4,-5)],my_function(t(data_gene)[c(-1,-2,-3,-4,-5)]))

colnames(transposed_label) <- c("val","index")

transposed_label$val <- as.factor(transposed_label$val)


###colapse each group of 3
new <- as.data.frame(transposed_label,stringsAsFactors = FALSE) %>% group_by(as.numeric(index)) %>% summarise(val = paste(val, collapse =","))

###Replace collapse with AA AT TT etc
new_col <- mapply(best_guess_to_allele,new$val,t(data_gene)[4],t(data_gene)[5])

new_table <- cbind(data_samples$ID_1[-1],new_col)

colnames(new_table) <- c("id","genotype")

new_table2 <- setDT(as.data.frame(new_table), keep.rownames = "best_guess")[]

new_table3 <- new_table2[,c(2,1,3)]


###to write the gene file as a transposed csv

write.csv(new_table3, file = "rs641738.csv")
