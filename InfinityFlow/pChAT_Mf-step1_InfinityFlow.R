### Section 0: Preparation --------------------------------------------------

# clear the environment
rm(list = ls())
graphics.off()
gc()

# Section 1: Setting up your input data---------------------------------------------------------------

## Load requred packages
suppressMessages(suppressWarnings(library(infinityFlow)))

# set your work diretory
workdir <- "/your_work_dir/data/Peritoneal_ChAT/InfinityFlow"
setwd(workdir)
set.seed(123)

# Set up input raw FACS files
input_dir <- file.path(workdir, "raw")
list.files(input_dir)

# Annotation of Infinity antibody targets and isotypes
annot <- read.csv(file.path(workdir,'auxilliaries/LegendScreen_annotation.csv'))
# View(annot)
targets <- annot$Infinity_target
names(targets) <- annot$file
head(targets)
isotypes <- annot$Infinity_isotype
names(isotypes) <- annot$file
head(isotypes)

# specifying the Backbone and Infinity antibodies
# backbone_specification <- select_backbone_and_exploratory_markers(list.files(input_dir, pattern = ".fcs", full.names = TRUE))
# write.csv(backbone_specification,
#     file = file.path(workdir, "auxilliaries/backbone_selection_file.csv"), ,
#     row.names = FALSE
# )
backbone_selection_file <- file.path(workdir, "auxilliaries/backbone_selection_file.csv")
# View( read.csv(backbone_selection_file))

# Section 2: Running the Infinity Flow computational pipeline---------------------------------------------------------------

# set prediction events per file
prediction_events_downsampling <- 2500

imputed_data <- infinity_flow(
    path_to_fcs = input_dir,
    path_to_output = file.path(workdir, prediction_events_downsampling, "out"),
    path_to_intermediary_results = file.path(workdir, prediction_events_downsampling, "tmp"),
    backbone_selection_file = backbone_selection_file,
    annotation = targets,
    isotype = isotypes,
    # input_events_downsampling = 1000, # default:Inf
    prediction_events_downsampling = prediction_events_downsampling,
    #extra_args_export = list(FCS_export = c("split", "concatenated", "csv")),
    extra_args_export = list(FCS_export = c("none")),
    extra_args_correct_background = list(FCS_export = c("split","concatenated")),
    your_random_seed = 123,
    verbose = TRUE,
    cores = 32
)
gc()
