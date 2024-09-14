#########################################
# read_density.R
#
# This script reads in the RNA-seq and Ribo-seq data.
#


require(tidyverse)
require(data.table)

read_density <- function(
    data_file = NULL,
    category = "RNA-Ribo"
    ) {

    # read in the data
    data <- fread(data_file, header = TRUE, check.names = FALSE)

    # merge the data
    data <- merge(data1, data2, by = gene_column, all = TRUE)

    # fill in the missing values
    data <- data %>%
        mutate_at(vars(contains("density")), funs(replace(., is.na(.), 0)))

    return(data)

}












