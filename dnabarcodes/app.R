###############################################################################
# Shiny app to calculate Hamming distances for a set of DNA barcode sequences.
## All barcodes should be of the same length.
## --
## Who:    Martin Oti
## What:   dnabarcodes/app.R
## Which:  1.0
## When:   2018-06-12
## --
###############################################################################

library(shiny)
library(ggplot2)
library(ggseqlogo)
library(pheatmap)
##library(DNABarcodes)        # not currently used
options(stringsAsFactors = FALSE)

## General setup

# Base colors for sequence logo
nextseq_basecolors <- make_col_scheme(chars = c('G', 'A', 'C', 'T'), 
                                      cols  = c('black', 'orange', 'red', 'green'))


####################
## Browser UI code
####################

ui <- fluidPage(
   
   # Application title
   headerPanel("DNA Barcodes Hamming Distances"),
   
   sidebarLayout(
     # Sidebar with file input panel and file contents display table
     sidebarPanel(
        downloadButton("downloadExampleFile", "Download Example File"),
        fileInput("barcodesfile",
                  "Upload barcodes file (format: id <tab> sequence)",
                  multiple = FALSE,
                  accept = c("text/tsv",
                             "text/tab-separated-values,text/plain",
                             ".txt")),
        sliderInput("minimumDistance",
                  "Minimum required distance",
                  min = 0, max = 9,
                  value = 3),
        plotOutput("seqlogoPlot"),
        tableOutput("barcodesTable")
      ),
      
      # Main panel with Hamming distance plots
      mainPanel(
         plotOutput("histogramPlot"),
         plotOutput("heatmapPlot"),
         tableOutput("conflictsTable")
      )
   )
)



####################
## Server logic
####################

server <- function(input, output) {
  
  ## Read in barcodes file and calculate Hamming distances
  
  datasetInput <- reactive({
    
    req(input$barcodesfile)
    
    tryCatch(
      {
        # Read input barcodes file 
        # Should be tab-separated file with header and two columns: ID & sequence
        df <- read.delim(input$barcodesfile$datapath)
        
        # Calculate Hamming distances between all barcode sequence pairs
        seq_pairs <- combn(df[,2], 2)
        hd_vec <- apply(seq_pairs, 2, function(x){
          x1 <- unlist(strsplit(x[1], ""))
          x2 <- unlist(strsplit(x[2], ""))
          sum(x1 != x2)
        })
        # Also get all barcode ID pairs
        id_pairs <- combn(df[,1], 2)
      },
      error = function(e) {
        # Return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    
    # Collect all relevant datasets into single results list
    results <- list("df" = df, "hd_vec" = hd_vec, "seq_pairs" = seq_pairs, "id_pairs" = id_pairs)
  })
  
  
  ## Handle downloading of example file
  
  output$downloadExampleFile <- downloadHandler(
    filename = "sample_barcodes.txt",
    content = function(tmpfilepath) {
      file.copy(from = "sample_barcodes.txt", to = tmpfilepath)
    },
    contentType = "text/tsv"
  )
  
  
  ## Create sequence logo
  
  output$seqlogoPlot <- renderPlot({
    
    res <- datasetInput()
    
    ggseqlogo(res$df[,2], 
              col_scheme = nextseq_basecolors,
              method = 'prob',
              font = 'helvetica_bold') +
      ggtitle("Barcodes sequence logo") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))
  })
  
  
  ## Create histogram plot
  
  output$histogramPlot <- renderPlot({
    
    res <- datasetInput()
    
    # Draw the histogram of Hamming distances
    hist(res$hd_vec, 
         breaks = seq(min(res$hd_vec)-0.5, max(res$hd_vec)+0.5, by = 1),
         col = 'darkgray', 
         xlab = 'Hamming distance', 
         main = 'Histogram of Hamming distances')
  })
  
  
  ## Create heatmap plot
  
  output$heatmapPlot <- renderPlot({
     
     res <- datasetInput()

     # Create barcode-vs-barcode matrix and convert Hamming distances vector into matrix
     hd_mat <- matrix(0, ncol = nrow(res$df), nrow = nrow(res$df))
     colnames(hd_mat) <- res$df[,2]                    # colnames = barcode sequences
     rownames(hd_mat) <- res$df[,2]                    # rownames = barcode sequences
     hd_mat[t(res$seq_pairs)] <- res$hd_vec            # fills upper triangle of matrix
     hd_mat[lower.tri(hd_mat)] <- t(hd_mat)[lower.tri(hd_mat)]  # copy to lower triangle
     # Rename columns and rows to sample names for heatmap plot
     colnames(hd_mat) <- res$df[,1]
     rownames(hd_mat) <- res$df[,1]
     
     # Draw the histogram of Hamming distances
     pheatmap(hd_mat, display_numbers = TRUE, number_format = "%d", number_color = "black")
  })
  
  
  ## Create table of conflicting barcodes
  
  output$conflictsTable <- renderTable({
    
    res <- datasetInput()
    
    # Identify barcode pairs with distance < minimum required distance, put in data.frame
    is_bad_pair <- res$hd_vec < input$minimumDistance
    bad_pairs_df <- data.frame("ID" = t(res$id_pairs)[is_bad_pair,,drop=FALSE], 
                               "Sequence" = t(res$seq_pairs)[is_bad_pair,,drop=FALSE], 
                               "Distance"  = res$hd_vec[is_bad_pair])
    
    bad_pairs_df[order(bad_pairs_df$Distance, decreasing = FALSE),]
  }, caption = "<b><span style='color:#000000'>Potentially conflicting barcodes</span></b>",
     caption.placement = "top"
  )  # end renderTable
  
  
  ## Create HTML table with barcodes data
  
  output$barcodesTable <- renderTable({ datasetInput()$df })
}


# Run the application 
shinyApp(ui = ui, server = server)

