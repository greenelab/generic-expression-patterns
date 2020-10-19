# recount2 contains over 70,000 uniformly processed human RNA-seq samples
# spanning TCGA, SRA and GTEx.
# This module includes a few functions to download recount2 SRA data.

library('recount')

set.seed(123)

# This function downloads the input project_id's expression data in
# `download_dir`, and saves its corresponding transposed data counts data
# as `raw_template_filename`.
get_recount2_template_experiment <- function(project_id,
                                             download_dir,
                                             raw_template_filename) {

  # Save current working directory
  #original_wd = getwd()

  # Change working directory to `download_dir
  #setwd(download_dir)

  # Get data associated with project ids
  # Download the RSE (RangedSummarizedExperiment) object.
  # The RSE object for the counts summarized at the gene level using the
  # Gencode v25 (GRCh38.p7, CHR) annotation as provided by Gencode.
  if (!file.exists(file.path(project_id, 'rse_gene.Rdata'))) {
    download_study(project_id)
  }
  load(file.path(project_id, 'rse_gene.Rdata'), verbose = TRUE)

  # Counts are the number of reads that map that each gene
  rse_gene_scaled <- scale_counts(rse_gene)
  rse_gene_counts <- assays(rse_gene_scaled)$counts
  data_counts <- t(rse_gene_counts)

  ## Save counts matrix to file
  write.table(
    data_counts,
    raw_template_filename,
    sep = '\t',
    row.names = TRUE,
    col.names = NA
  )

  # Set working directory back
  #setwd(original_wd)
}


get_recount2_sra_subset <- function(template_project_id,
                                    num_experiments,
                                    local_dir,
                                    base_dir) {
  # This function downloads 50 random experiments, including the
  # the experiment associated with the input project_id.
  # These experiments makeup a compendium of gene expression experiments
  # that are saved in the local_dir

  # Download metadata file
  metadata <- all_metadata()
  write.table(
    metadata,
    paste(base_dir, '/human_general_analysis/data/metadata/recount2_metadata.tsv', sep = ""),
    sep = '\t',
    row.names = FALSE
  )

  # Format of the metadata
  # Based on definitions from NCBI and blog post:
  # https://www.ncbi.nlm.nih.gov/books/NBK56913/#search.what_do_the_different_sra_accessi
  # https://www.ccdatalab.org/blog/2019/3/29/gene-expression-repos-explained
  # Project: A sequencing study (i.e. NCBI sequencing study).
  # Sample: Physical biospecimen on which sequencing was performed, biological source
  #         material (i.e. HeLa cell line). A project contains many samples
  # Experiment: Unique sequencing library for a specific sample.
  #             A sample can have multiple experiments, most have 1 experiment.
  # Run: Sequencing run.  An experiment can contain many runs (i.e. technical replicates)
  # In this case, we want to group runs into projects (for "experiment-level" simulation)

  metadata <- metadata[metadata[,"read_count_as_reported_by_sra"] > 0,]
  project_ids <- unique(metadata$project)

  # Entire recount2 is 8TB, we will only select the subset of studies instead
  selected_project_ids <- sample(project_ids, num_experiments, )

  # Add user selected project id
  selected_project_ids <- append(selected_project_ids, template_project_id)

  # Get data associated with project ids
  # Download the RangedSummarizedExperiment object
  # The RangedSummarizedExperiment object for the counts summarized
  # at the gene level using the Gencode v25 (GRCh38.p7, CHR) annotation as provided by Gencode.
  for (i in 1:length(selected_project_ids)) {
    print(selected_project_ids[i])
    if (!file.exists(file.path(selected_project_ids[i], 'rse_gene.Rdata'))) {
      download_study(selected_project_ids[i])
    }
    load(file.path(selected_project_ids[i], 'rse_gene.Rdata'), verbose = TRUE)

    # Counts are raw read counts from the sequencing run
    # Counts are the number of reads that map that each gene

    # Concatenate counts into one matrix
    if (i == 1) {
      data_counts_all <- assays(scale_counts(rse_gene))$counts
    }
    else {
      data_counts_all <- cbind(data_counts_all, assays(scale_counts(rse_gene))$counts)
    }

    ## Delete it if you don't need it anymore
    unlink(selected_project_ids[i], recursive = TRUE)
  }

  data_counts_all <- t(data_counts_all)

  ## Save counts matrix to file
  write.table(
    data_counts_all,
    paste(local_dir,'/recount2_compendium_data.tsv',sep = ""),
    sep = '\t',
    row.names = TRUE,
    col.names = NA
  )
}

# This function saves recount2 metadata into `metadata_dir`, downloads
# `rse_gene.Rdata` file each recount2 SRA project (experiment) in `download_dir`,
# and saves the corresponding transposed data count file (`t_data_counts.tsv`)
# in `download_dir` as well.
# `overwrite` is a boolean flag to indicate whether the dowloading process will
# overwrite any previously downloaded files. If it is FALSE and the previously
# downloaded files were complete (`validation_filename` exists), no files will
# be downloaded; otherwise previously downloaded files will be overwritten.
# `validation_filename` is the name of an empty file, which will be generated
# in `download_dir` at the end to indicate the completion of the downloading.
#
# Note: Downloading may take 30 minutes to a few hours, depending on the speed
# of your network.
download_recount2_sra <- function (metadata_dir,
                                   download_dir,
                                   overwrite = TRUE,
                                   validation_filename = 'download_completed.txt') {
  # If `overwrite` is FALSE and previous downloading was complete, exit
  if (!overwrite && file.exists(file.path(download_dir, validation_filename))) {
    print("Downloading process skipped because previous one is complete")
    return()
  }

  # Download metadata file
  metadata <- all_metadata()
  write.table(
    metadata,
    paste(metadata_dir, 'recount2_metadata.tsv', sep = "/"),
    sep = '\t',
    row.names = FALSE
  )

  # Format of the metadata:
  # Based on definitions from NCBI and blog post:
  # https://www.ncbi.nlm.nih.gov/books/NBK56913/#search.what_do_the_different_sra_accessi
  # https://www.ccdatalab.org/blog/2019/3/29/gene-expression-repos-explained
  # Project: A sequencing study (i.e. NCBI sequencing study).
  # Sample: Physical biospecimen on which sequencing was performed, biological source
  #         material (i.e. HeLa cell line). A project contains many samples
  # Experiment: Unique sequencing library for a specific sample.
  #             A sample can have multiple experiments, most have 1 experiment.
  # Run: Sequencing run.  An experiment can contain many runs (i.e. technical replicates)
  # In this case, we want to group runs into projects (for "experiment-level" simulation)

  metadata <- metadata[metadata[, "read_count_as_reported_by_sra"] > 0,]
  project_ids <- unique(metadata$project)

  # Save current working directory
  original_wd = getwd()

  # Change working directory to `download_dir`
  setwd(download_dir)

  for (idx in 1:length(project_ids)) {
    pid <- project_ids[idx]
    # As of 09/25/2020, the project ID "SRP053791" (the last entry in `project_ids`)
    # CAN NOT be downloaded due to the following error generated by download_study() function:
    #   Error in download_study(last_p) :
    #     Invalid 'project' argument. There's no such 'project' in the recount_url data.frame.
    # so this project is skipped manually to prevent R session from premature termination.
    if (pid == "SRP053791") {
      next
    }

    download_study(pid)
    load(file.path(pid, 'rse_gene.Rdata'), verbose = TRUE)

    data_counts <- assays(scale_counts(rse_gene))$counts
    t_data_counts <- t(data_counts)
    write.table(
      t_data_counts,
      paste(download_dir, pid, 't_data_counts.tsv', sep = "/"),
      sep = '\t',
      row.names = TRUE,
      col.names = NA
    )
  }

  # Create the validation file to indicate the completion of download
  blank_df <- data.frame()
  write.table(blank_df, file = validation_filename, col.names = FALSE)

  # Set working directory back
  setwd(original_wd)
}
