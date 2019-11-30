
#Same function as IMAAP::imaap_annotate_cells except Nathan added sample name to file output

imaap_annotate_cells_nj <- function (data, cheat_sheet, SD,sample_name)
{
  s = imaap_model(data, SD)
  marker_prob = data.frame(s[1])
  colnames(marker_prob) = colnames(data)
  cell_prob <- imaap_prob(marker_prob, cheat_sheet)
  cell_annotation <- imaap_annotation(cell_prob, cheat_sheet)
  modified_cell_annotation <- imaap_marker_drop(cell_prob,
                                                cell_annotation, cheat_sheet)
  peaks_low = data.frame(s[2])
  peaks_high = data.frame(s[3])
  my_plots = list()
  my_plots <- imaap_model_plot(data = data.frame(s[4]), peaks_low,
                               peaks_high, SD)
  dir.create("Fitted distribution plots", showWarnings = FALSE)
  pdf(paste("./Fitted distribution plots/",sample_name,"_distibution plot",
            ".pdf",sep=""), width = 9.5, height = 12)
  grid.arrange(grobs = my_plots, ncol = 4)
  dev.off()
  row.names(cell_annotation) = cell_annotation$cell_id
  cell_annotation = cell_annotation[, -c(1, 2)]
  row.names(modified_cell_annotation) = modified_cell_annotation$cell_id
  modified_cell_annotation = modified_cell_annotation[, -c(1,
                                                           2)]
  cell_collapsed <- imaap_one_call_per_cell(cell_calling = modified_cell_annotation)
  cell_collapsed[, 1] <- as.character(cell_collapsed[, 1])
  cell_calling <- list(cell_annotation, modified_cell_annotation,
                       cell_collapsed)
  return(cell_calling)
}
