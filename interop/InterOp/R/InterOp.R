readInterOpFiles <- function(path = "./") {

	f <- c("ExtractionMetricsOut.bin",
		   "QMetricsOut.bin",
		   "ErrorMetricsOut.bin",
		   "TileMetricsOut.bin",
		   "CorrectedIntMetricsOut.bin",
		   "ControlMetricsOut.bin",
		   "ImageMetricsOut.bin")

	# progress bar
	progress <- function(msg,current) { 
		cat(paste(c("\r[",rep("=",current),rep(" ",100-current),"]",round(current),"% ",msg,rep(" ",25)),collapse="")) 
		flush.console() 
		Sys.sleep(1)
	}

	# readExtractionMetrics
	progress(paste("reading",f[1]),100*1 / (length(f) + 1))
	extraction_metrics <- readExtractionMetrics(paste(path,f[1],sep="/"))
	extraction_metrics$datetime <- as.POSIXlt(extraction_metrics$datetime / 10000000,origin="0001-01-01")

	# readQualityMetrics
	progress(paste("reading",f[2]),100*2 / (length(f) + 1))
	quality_metrics <- readQualityMetrics(paste(path,f[2],sep="/"))

	# readErrorMetrics
	progress(paste("reading",f[3]),100*3 / (length(f) + 1))
	error_metrics <- readErrorMetrics(paste(path,f[3],sep="/"))

	# readTileMetrics
	progress(paste("reading",f[4]),100*4 / (length(f) + 1))
	tile_metrics <- readTileMetrics(paste(path,f[4],sep="/"))

	# readCorrectedIntMetrics
	progress(paste("reading",f[5]),100*5 / (length(f) + 1))
	corrected_int_metrics <- readCorrectedIntMetrics(paste(path,f[5],sep="/"))

	# readControlMetrics
	progress(paste("reading",f[6]),100*6 / (length(f) + 1))
	control_metrics <- readControlMetrics(paste(path,f[6],sep="/"))

	# readImageMetrics
	progress(paste("reading",f[7]),100*7 / (length(f) + 1))
	image_metrics <- readImageMetrics(paste(path,f[7],sep="/"))

	# output object
	progress("done",100); cat("\n"); flush.console() 
	structure(list("extraction_metrics"    = extraction_metrics,
				   "quality_metrics"       = quality_metrics,
				   "error_metrics"         = error_metrics,
				   "tile_metrics"          = tile_metrics, 
				   "corrected_int_metrics" = corrected_int_metrics, 
				   "control_metrics"       = control_metrics, 
				   "image_metrics"         = image_metrics),
			  class="InterOp")

}

#########################
##
## plot min and max contrasts for the ACGT channels
## uses the ggplot2 package
##
#########################
plotImageContrasts <- function(x) UseMethod("plotImageContrasts",x)	# generic method: function dispatcher
plotImageContrasts.InterOp <- function(iop) {

	if(require("ggplot2") && require("reshape")) {
		# melting to a ggplot suitable dataframe
		df <- melt.data.frame(data=iop$image_metrics,id.vars=c("lane","channelid"),measure.vars=c("mincont","maxcont"))
		df$channelid <- sapply(as.character(df$channelid),function(x) chartr("0123","ACGT",x))

		p <- ggplot(data=df,aes(x=value)) + 
				geom_density(aes(fill=factor(channelid),y=..scaled..),alpha=.25) + 
				coord_trans(x="log10") +
				theme_bw()

		p + facet_grid(variable ~ lane)
	}
}
