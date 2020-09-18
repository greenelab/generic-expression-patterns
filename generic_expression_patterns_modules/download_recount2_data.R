## Run this once to setup environment
## Used R 3.6.3

library('recount')

## About
# recount2 contains over 70,000 uniformly processed human RNA-seq samples spanning
# TCGA, SRA and GTEx

set.seed(123)

get_recount2_template_experiment <- function(project_id, local_dir) {
  # This function downloads the gene expression data associated
  # with the input project_id and saves it in the local_dir

  # Reset working directory to `local_dir`
  setwd(local_dir)

  ## Data

  # Get data associated with project ids
  # Download the RSE (RangedSummarizedExperiment) object.
  # The RSE object for the counts summarized at the gene level using the
  # Gencode v25 (GRCh38.p7, CHR) annotation as provided by Gencode.
  if (!file.exists(file.path(project_id, 'rse_gene.Rdata'))) {
    download_study(project_id)
  }
  load(file.path(project_id, 'rse_gene.Rdata'), verbose=TRUE)

  # Counts are the number of reads that map that each gene

  rse_gene_scaled <- scale_counts(rse_gene)
  rse_gene_counts <- assays(rse_gene_scaled)$counts
  data_counts <- t(rse_gene_counts)

  ## Save counts matrix to file
  write.table(
    data_counts,
    paste(local_dir, '/recount2_template_data.tsv', sep = ""),
    sep = '\t',
    row.names = TRUE,
    col.names = NA
  )

}


get_recount2_subset_compendium <- function(template_project_id,
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
    paste(base_dir, '/human_analysis/data/metadata/recount2_metadata.tsv', sep = ""),
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

  # Entire recount2 is 8TB
  # We will only select the subset of studies instead
  selected_project_ids <- sample(project_ids, num_experiments, )

  # Add user selected project id
  selected_project_ids <- append(selected_project_ids, template_project_id)

  # Reset working directory to `local_dir`
  setwd(local_dir)

  ## Data

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

# This function downloads all recount2 `rse_gene.Rdata` files and generates a
# transposed data count file (`t_data_counts.tsv`) based on each `rse_gene.Rdata`.
function download_all_recount2(metadata_dir, download_dir) {
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

  # Reset working directory to `download_dir`
  setwd(download_dir)

  for (idx in 1:length(project_ids)) {
    pid <- project_ids[idx]
    if (!file.exists(file.path(pid, 'rse_gene.Rdata'))) {
      download_study(pid)
    }

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

}
