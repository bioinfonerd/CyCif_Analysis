#Tile: "Tutorial"
#author: "Nathan T. Johnson"
#date: "11/7/2019"
#Purpose: Provide basic analysis infrastructure to analyze CyCif Data

# Load in Library -----------------------------------------------------------------
setwd('C:/Users/Nathan/Dropbox/@Dana Farber/CyCif/git/CyCifAnalysis/vignettes')


#load all required libraries
list.of.packages <- c("tidyverse","stringr","tibble","readr",'BiocManager',"dplyr","tximport","tximportData",
                      "ggplot2","RColorBrewer","Hmisc","devtools","PTXQC","gridExtra","cowplot")
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

#load in Ajit's cell type caller
if( !require(devtools) ) install.packages("devtools")
devtools::install_github("labsyspharm/IMAAP")
library(IMAAP)

#Another possibility is Artem's cell type caller, but sticking with Ajit's for now
# #load in Artem's cell type caller (Naive Bayes)
# if( !require(devtools) ) install.packages("devtools")
# devtools::install_github( "labsyspharm/naivestates" )
# library( naivestates )

#install in local directory the package specific functions
#devtools::install()

#install documentation
devtools::document()

#until build reflects functions loaded
devtools::load_all()


# Location Variables, Data, and Data Cleaning ------------------------------------------------

#where is your data
location <- '../data/test_run_5channels_10cycles_Data_Files'

#where are images being stored
output_location <- '../results/'

#load in all data as tibble objects (function may want to be amended if alot of data)
mydata <- load_all_csv_files(location)

#clean up data column names, replace sample name from number to name, and add what type of data (Cell, Nucleus) to improve reading and be able to work with Ajit's cell type calling (IMAAP)
mydata <- clean_all_column_cleaner(mydata)

#cell type definition
cell_type_definition <- read.table(file = "../data/cell_type_definition.csv", sep=',', header = T, stringsAsFactors=FALSE)

# Basic Plots QC ------------------------------------------------

###########################################
#plot individual sample on simple x,y plot#
###########################################

  #possible markers to select
  markers <- names(mydata[[1]])
  markers

  #With natural log normalization
  data <- mydata[[1]] %>% select("X_position","Y_position","Ecadherin")

  plot <- ggplot(data, aes(X_position, Y_position, color = Ecadherin)) +
    geom_point(size = 1) +
    theme_minimal()

  #display plot (can take a while if alot of cells) - would recommend saving plot and then viewing as will render faster
  #print(plot)

  #save plot (can take a while if alot of cells)
  dir.create(file.path(output_location, "xy_images"), showWarnings = FALSE) #make directory if not found
  ggsave(file=paste(output_location,"/xy_images/test_xy.png",sep=""), plot=plot, width=10, height=8)
  dev.off() #clear plotting space

###########################################################################
# plot density plot for 1 or several images for 1 or several distributions#
###########################################################################
# probability density is the probability per unit on the x-axis


# plot density plot for 1 images for 1 marker

  #possible markers to select
  markers <- names(mydata[[1]])
  markers

  #organize data for plotting 1 image
  marker <- "Ecadherin" #selected marker
  data <- mydata[[1]] %>% select("ImageId",marker)

  #plot
  plot <- ggplot(data, aes_string(x=marker)) +
    geom_density(color="darkblue", fill="lightblue") +
    xlim(5,7)

  #display plot
  plot

  #save plot
  dir.create(file.path(output_location, "density_plots"), showWarnings = FALSE) #make directory if not found
  ggsave(file=paste(output_location,"/density_plots/test_1_image_1_marker_density.png",sep=""), plot=plot, width=10, height=8)
  dev.off() #clear plotting space

# plot density plot for 1 image for multiple markers

  #possible markers to select
  markers <- names(mydata[[1]])
  markers

  #select markers for plotting 1 image with multiple markers
  marker <- c("Ecadherin","Actin555","FOXP3","CD45")
  data <- mydata[[1]] %>% select("ImageId",marker)

  #reorganize data for plotting
  data <- data %>% gather(key="Marker",value="Expression",-ImageId)

  #plot
  plot <- ggplot(data, aes(x=Expression,fill=Marker)) +
    geom_density(alpha=0.4)

  #display plot
  plot

  #save plot
  dir.create(file.path(output_location, "density_plots"), showWarnings = FALSE) #make directory if not found
  ggsave(file=paste(output_location,"/density_plots/test_1_image_multiple_markers_density.png",sep=""), plot=plot, width=10, height=8)
  dev.off() #clear plotting space

# plot density plot for multiple images for 1 marker

  #possible markers to select
  markers <- names(mydata[[1]])
  markers

  #select markers for plotting 1 image with multiple markers (will take a while depending on number of images) *needs improvement when # of cells increases
  #marker selected *Marker must be present for each image
  marker <- "Actin555"
  for(i in 1:length(mydata)){
    if(i == 1){
      data <- mydata[[i]] %>% select("ImageId",marker)
    }
    else{
      tmp <- mydata[[i]] %>% select("ImageId",marker)
      data <- rbind(data,tmp)
    }
  }

  #list of plot objects
  plot_list <- list() #list to store plots
  for(i in 1:length(unique(data$ImageId))){
    tmp <- data %>% filter(ImageId == unique(data$ImageId)[i])
    plot_list[[i]] <- ggplot(tmp, aes_string(x=marker)) +
      geom_density(color="darkblue", fill="lightblue") + ggtitle(unique(data$ImageId)[i])
  }
  #plot
  n <- length(plot_list)
  nCol <- floor(sqrt(n))
  plot <- do.call("grid.arrange", c(plot_list, ncol=nCol))

  #save plot
  dir.create(file.path(output_location, "density_plots"), showWarnings = FALSE) #make directory if not found
  ggsave(file=paste(output_location,"/density_plots/test_multiple_images_1_marker_density.png",sep=""), plot=plot, width=10, height=8)
  dev.off() #clear plotting space

# plot density plot for multiple images for multiple markers

  #possible markers to select
  markers <- names(mydata[[1]])
  markers

  #select markers for plotting 1 image with multiple markers (will take a while depending on number of images) *needs improvement when # of cells increases
  #marker selected *Marker must be present for each image
  marker <- c("Ecadherin","Actin555","FOXP3","CD45")
  for(i in 1:length(mydata)){
    if(i == 1){
      data <- mydata[[i]] %>% select("ImageId",marker)
    }
    else{
      tmp <- mydata[[i]] %>% select("ImageId",marker)
      data <- rbind(data,tmp)
    }
  }

  #reorganize data for plotting
  data <- data %>% gather(key="Marker",value="Expression",-ImageId)

  #list of plot objects
  plot_list <- list() #list to store plots
  for(i in 1:length(unique(data$ImageId))){
    tmp <- data %>% filter(ImageId == unique(data$ImageId)[i])

    plot_list[[i]] <- ggplot(tmp, aes(x=Expression,fill=Marker)) +
      geom_density(alpha=0.4) + ggtitle(unique(data$ImageId)[i])
  }

  #plot
  n <- length(plot_list)
  nCol <- floor(sqrt(n))
  plot <- do.call("grid.arrange", c(plot_list, ncol=nCol))

  #save plot
  dir.create(file.path(output_location, "density_plots"), showWarnings = FALSE) #make directory if not found
  ggsave(file=paste(output_location,"/density_plots/test_multiple_images__multiple_markers_density.png",sep=""), plot=plot, width=10, height=8)
  dev.off() #clear plotting space

# Cell State Calling (Ajit)---------------------------------------

  #store output
  celltype_output<-list()

#careful running all of this on your entire dataset as it will take a while (example dataset only first image (i = 1) works without bugs.)
for(i in 1:length(mydata)){
  #identify what marker names to use for analysis
  #ideally only include your research/cell-type markers, so exclude DNA and NA channels (but depending on analysis that may not be ideal)
  markers <- names(mydata[[i]])
  #Display the column names
  markers
  #Remove markers with DNA in them
  markers <- markers[!grepl("DNA",markers)]
  #Remove markers with NA + a number in them
  markers <- markers[!grepl('NA[0-9]',markers)]
  #verify removed
  markers
  #what columns are marker columns (for this example have two marker sets, so run it slightly different)
  #MAKE SURE TO SELECT MARKER COLUMNS that have cell type calling present
  #SELECT MARKERS just for cell type to run faster
  markers <- markers[3:38]
  #verify selecting correct columns
  markers
  #select for those markers
  IMAAP_data <- mydata[[i]] %>% select(markers)
  #remove any rows with inifinte
  IMAAP_data <- IMAAP_data[is.finite(rowSums(IMAAP_data)),]
  #remove any rows with NA present
  IMAAP_data <- IMAAP_data %>% drop_na()
  #
  IMAAP_data <- exp(IMAAP_data)-1
  #keep Cell Id
  CellId <- as.integer(row.names(IMAAP_data))
  #annotate cells by cell type
  cell_type_definition <- read.table(file = "../data/cell_type_definition.csv", sep=',', header = T, stringsAsFactors=FALSE)
  #NATHAN's FUNCTION
  annotate_the_cells <-imaap_annotate_cells_nj (data = IMAAP_data, cheat_sheet = cell_type_definition, SD = 3, sample_name = unique(mydata[[i]]$ImageId))
  #Save Cell Type Call assign cell type as column
  output <- annotate_the_cells[[3]]
  output$CellId <- CellId
  output$ImageId <- unique(mydata[[i]]$ImageId)
  celltype_output[[i]] <- output
}

#save results
  dir.create(file.path(output_location, "Cell_Type_Calling"), showWarnings = FALSE) #make directory if not found
  for(i in 1:length(celltype_output)){
    write.csv(celltype_output[[i]],file=paste(output_location,"Cell_Type_Calling/",unique(celltype_output[[i]]$ImageId),"_cell_type_calling.csv",sep=""),sep = ",",
              row.names = FALSE)
  }

# Statistics on Cell Type Calling -----------------------------------------

#table for number of cells called
celltype_output[[1]] %>% group_by( Level.4) %>% count()
