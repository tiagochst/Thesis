
> mae
A MultiAssayExperiment object of 2 listed
 experiments with user-defined names and respective classes.
 Containing an ExperimentList class object of length 2:
 [1] DNA methylation: RangedSummarizedExperiment with 135331 rows and 866 columns
 [2] Gene expression: RangedSummarizedExperiment with 57035 rows and 866 columns
Features:
 experiments() - obtain the ExperimentList instance
 colData() - the primary/phenotype DataFrame
 sampleMap() - the sample availability DataFrame
 `$`, `[`, `[[` - extract colData columns, subset, or experiment
 *Format() - convert ExperimentList into a long or wide DataFrame
 assays() - convert ExperimentList to a list of rectangular matrices
> table(mae$definition)
         Metastatic Primary solid Tumor Solid Tissue Normal 
                  5                 778                  83 
