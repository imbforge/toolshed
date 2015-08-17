##################################
##
## Assess sequencing depth saturation by sampling
## reads from the counts tables
## --
## Who:  Sergi Sayols
## When: 18-feb-2015
## --
## Input:
##   folder=<input folder>
##   pattern=<analyze only files matching this pattern (valid R regular expression)>
##   out=<output file name>
##   pre=<prefix to be removed from sample name (for plotting)>
##   suf=<suffix to be removed from sample name (for plotting)>
##   steps=<number of subsampling steps>
##   simulations=<number of simulations at every step>
##   cores=<number of cores to use>
## --
## Run it from bash: double escape special characters
##   $ Rscript saturation.R folder=./test pattern="\\\\.tsv$" pre="Sample_imb_richly_2014_06_\\\\d+_" suf="\\\\_readcounts.tsv"
## Run it from bsub: double-double escape special characters
##   $ echo "Rscript saturation.R folder=./test pattern=\\"\\\\\\\\.tsv\\$\\" pre=\\"Sample_imb_richly_2014_06_\\\\\\\\d+_\\" suf=\\"\\\\\\\\_readcounts.tsv\\"" | bsub
##
########################################
library(ggplot2)
library(parallel)

##
## Parse input parms
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {

	if(length(i <- grep(string,args,fixed=T)) == 1) 
		return(do.call(convert,list(gsub(string,"",args[i]))))
    
	if(!is.null(default)) default else do.call(convert,list(NA))
}

args    <- commandArgs(T)
FOLDER  <- parseArgs(args,"folder=","./")       # folder containing the bam files
PATTERN <- parseArgs(args,"pattern=","\\_readcounts.tsv$") # files to be analyzed. Default: all bams
OUT     <- parseArgs(args,"out=","sat.pdf")    # output filename
PRE     <- parseArgs(args,"pre=","")            # pattern to remove from the file name
SUF     <- parseArgs(args,"suf=","\\_readcounts.tsv$")     # pattern to remove from the file name
STEPS   <- parseArgs(args,"steps=",10,"as.numeric")	# number of subsampling steps
SIM     <- parseArgs(args,"simulations=",10,"as.numeric")	# number of simulations at every step
CORES   <- parseArgs(args,"cores=",1,"as.numeric") # number of cores to use

print(args)
if(length(args) == 0 | args[1] == "-h" | args[1] == "--help")
	stop(paste("Rscript saturation.R [arguments|-h|--help]\n",
			   "  [folder=./]         : input folder\n",
			   "  [pattern=\"\\\\_readcounts.tsv$\"] : analyze only files matching this pattern (valid R regular expression)\n",
			   "  [out=sat.pdf]    : output filename\n",
			   "  [pre=\"\"]          : prefix to be removed from sample name (for plotting)\n",
			   "  [suf=\"\\\\_readcounts.tsv$\"]     : suffix to be removed from sample name (for plotting)\n",
			   "  [steps=10]          : number of subsampling steps",
			   "  [simulations=10]    : number of simulations at every step",
			   "  [cores=1]           : number of cores to use"))
if(is.na(STEPS))         stop("steps has to be an integer number")
if(is.na(SIM))           stop("simulations has to be an integer number")
if(is.na(CORES))         stop("cores has to be an integer number")
if(!file.exists(FOLDER)) stop(paste("Dir",FOLDER,"does NOT exist"))
files <- list.files(path=FOLDER,pattern=PATTERN)
files <- files[grep(PATTERN,files)]
samples <- gsub(PRE,"",gsub(SUF,"",files))
if(length(files) == 0) stop(paste("Dir",FOLDER,"does not contain files matching the pattern. Nothing to do"))

##
## read counts
##
counts <- lapply(paste0(FOLDER,"/",files),read.delim,head=F,row.names=1)
counts <- do.call(cbind,counts)
colnames(counts) <- samples

##
## Subsample counts STEPS times (10%..100% of the reads) and calculate the number
##   of detected genes per sample.
## How to calculate the number of features detected at every step:
##    -calculate SIM multinomially distributed random counts vectors with 'counts' probabilities
##    -count how many detected genes (counts > 0) at eevry vector
##    -average the number of detected genes over the SIM vectors
## res is a list of length #samples containing a data frame with x,y coordinates (depth,detected genes)
res <- apply(counts,2,function(x) {
	l <- mclapply(seq(.1,1,length=STEPS),function(i) {
		depth    <- round(i * sum(x))
		counts   <- rmultinom(SIM,size=depth,prob=x)
		detected <- round(mean(apply(counts,2,function(x) sum(x > 0))))
		c(genes=detected,depth=depth)
	},mc.cores=CORES)
	d <- data.frame(do.call(rbind,l))
	d$diff <- c(0,diff(d$genes))
	d
})

##
## plot the results
##
pdf(OUT)

df <- data.frame(depth =unlist(lapply(res,function(x) x$depth),use.names=F) / 10^6,
				 genes =unlist(lapply(res,function(x) x$genes),use.names=F),
				 diff  =unlist(lapply(res,function(x) x$diff) ,use.names=F),
				 sample=factor(rep(samples,each=STEPS)))
p <- ggplot(df,aes(depth,genes,color=sample)) +
	geom_point(size=3,alpha=.5) +
	geom_line(size=2,alpha=.5) +
#	geom_smooth(se=F) +
	scale_x_continuous(name="Sequencing depth (million reads)") +
	scale_y_continuous(name="Number of detected features") +
	theme_bw()

print(p)

dev.off()
