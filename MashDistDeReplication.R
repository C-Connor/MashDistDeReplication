if (!("optparse" %in% installed.packages())){
  install.packages("optparse",dependencies = TRUE, repos = "http://cran.uk.r-project.org")
}
library("optparse")

if(!("ape" %in% installed.packages())){
  install.packages("ape", dependencies = TRUE, repos = "http://cran.uk.r-project.org")
}
library("ape")

option_list <-list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Mashdist input file", dest="input"),
  make_option(c("-o", "--out"), type="character", default="MashDeRep", help="Output path and prefix [default= %default]", dest="output"),
  make_option(c("-t","--threshold"), type = "numeric", default = 0, help = "Maximum mash distance for genomes to be considered duplicates [less than or equal].",dest = "threshold")
  )

arg_parser = OptionParser(option_list=option_list)
cl_args = parse_args(arg_parser)

#check user input
if (is.null(cl_args$input)){
  print_help(arg_parser)
  stop("Input file must be provided", call.=FALSE)
}

#remove spaces from output
cl_args$output <- gsub(" ","",cl_args$output)

#make output directory if it doesn't exist
if (!dir.exists(dirname(cl_args$output))){
  dir.create(dirname(cl_args$output))
}

#read in mash distance file
print('Reading input file...')
mashin <- read.table(cl_args$input, row.names = 1)
colnames(mashin) <- row.names(mashin)
mashdf <- as.data.frame(mashin)

#function for checking if pairing already exists
checkfun<-function(inlist,inpair){
  for(il in 1:length(inlist)){
    if(any(is.element(inlist[[il]],inpair))){
      inlist[[il]] <- union(inlist[[il]],inpair)
      return(inlist)
      }}
    inlist[[length(inlist)+1]] <- inpair
    return(inlist)
}

#analyse mash distance matrix
print('Analysing mash distance matrix...')
colindex <- c()
count <- 1
replist <- list()
for(y in 1:nrow(mashdf)){
  colindex <- append(colindex,count)
  count <- count +1
  for(x in colindex){
    if((mashdf[y,x] <= cl_args$threshold) && (rownames(mashdf)[y] != colnames(mashdf)[x])){
      print(c('Identical pair: ',rownames(mashdf)[y],colnames(mashdf)[x]))
      pair <- c(rownames(mashdf)[y],colnames(mashdf)[x])
      if(length(replist)==0){
        replist[[1]] <- pair
      }else{
        replist <- checkfun(replist,pair)
      }
    }
  }
}

#make graph
print('Writing graph file...')
pdf(paste0(cl_args$output,"_mashdist_tree.pdf"), width = 100, height = 100)
plot(as.phylo(hclust(as.dist(mashdf))),cex=0.1)
dev.off()

#write out list of duplicates
print('Writing output file...')
for(out in replist){
  write(out[-1],file = paste0(cl_args$output,"_duplicate_genomes.txt"),ncolumns = 1, append = T)
}
if (length(replist)==0){
  print('No identical genomes identified, no output file will be produced.')
}
print('Complete')