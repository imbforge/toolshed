########################################
##
## Translate the chromosome names in a BW file
## --
## Who:  Sergi Sayols
## When: 23-jan-2015
## --
## Input:
##   in=<input file>
##   out=<output file>
##   table=<translation table>
##
########################################

##
## Parse input parms
##
parseArgs <- function(args, string, default=NULL, convert="as.character") {

	if(length(i <- grep(string, args, fixed=T)) == 1) 
		return(do.call(convert, list(gsub(string, "", args[i]))))
    
	if(!is.null(default)) default else do.call(convert, list(NA))
}

args  <- commandArgs(T)
IN    <- parseArgs(args, "in=")
OUT   <- parseArgs(args, "out=")
TABLE <- parseArgs(args, "table=")

if(length(args) == 0 | args[1] == "-h" | args[1] == "--help")
	stop(paste("Rscript BWtranslator.R [arguments|-h|--help]\n", 
			   "  [in=in.bw]       : the input file in BigWig format\n", 
			   "  [out=out.bw]     : the output file in BigWig format\n", 
			   "  [table=table.txt]: <TAB> separated table with the translations"))
if(!file.exists(IN))    stop(paste("File", IN, "does NOT exist"))
if(!file.exists(TABLE)) stop(paste("File", TABLE, "does NOT exist"))

##
## PROGRAM
##
# load libraries
library(rtracklayer)

# Read in the input file
in.bw <- import.bw(BigWigFile(IN))
tab   <- read.delim(TABLE, comment.char="#")
seqlevels(in.bw) <- Reduce(function(x, i) gsub(paste0("^", tab[i, 1], "$"), tab[i, 2], x), 1:nrow(tab), init=seqlevels(in.bw))
export.bw(in.bw, BigWigFile(OUT))
