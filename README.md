# Alzheimer's Disease Transcriptomics Analysis Pipeline

## Project Overview
This computational pipeline performs comprehensive transcriptomic analysis of Alzheimer's disease using RNA-seq data from temporal cortex samples, focusing on differential gene expression, cell type specificity, and pathway enrichment.

## Repository Structure

### Data Files
- `fpkm_table_unnormalized.csv`: Raw gene expression data
- `columns-samples.csv`: Sample metadata
- `DonorInformation.csv`: Donor-specific information
- `row-genes.csv`: Gene annotation information
- `tbi_data_files.csv`: RNA integrity metrics
- `trimmed_means.csv`: Cell type expression data

### Analysis Scripts
- `alzheimers_differential_expression_analysis.py`: Primary differential expression analysis
- `alzheimers_cell_type_and_pathway_analysis.py`: Cell type and pathway enrichment analysis
- `simple-gsea.R`: Gene Set Enrichment Analysis (GSEA)

### Output Files
- `alz_vs_healthy_tcx_results.csv`: Differential expression results
- `top_100_genes.csv`: Top genes by significance
- `ranked_genes_for_gsea.csv`: Ranked gene list for enrichment analysis
- `gsea_results_go.csv`: GO pathway enrichment results
- `gsea_results_reactome.csv`: Reactome pathway enrichment results

### Visualization Outputs
- `volcano_plot.png`: Volcano plot of differential expression
- `pvalue_distribution.png`: P-value distribution
- `logFC_distribution.png`: Log fold change distribution
- `de_genes_cell_type_category_heatmap.png`: Cell type expression heatmap
- `pathway_process_distribution.png`: Pathway distribution (if generated)

## Key Analysis Stages

1. **Differential Expression Analysis**
   - Data preprocessing
   - Filtering of low-expression genes
   - Statistical testing using limma
   - Multiple testing correction

2. **Cell Type Specificity**
   - Mapping genes to cell type expression data
   - Identifying cell type-specific gene expression patterns

3. **Pathway Enrichment**
   - Gene Set Enrichment Analysis (GSEA)
   - Identification of significant biological pathways
   - Exploration of GO and Reactome pathways

## Key Findings Summary
- Identified significant gene expression changes in Alzheimer's disease
- Mapped gene expression to specific neural cell types
- Explored potential molecular mechanisms underlying neurodegeneration

## Required Dependencies
- Python:
  - pandas
  - numpy
  - matplotlib
  - seaborn
  - scipy
  - rpy2
  - patsy
- R:
  - limma
  - fgsea (optional)
- Other:
  - R statistical computing environment

## Installation
1. Clone the repository
2. Install Python dependencies: `pip install -r requirements.txt`
3. Install R dependencies: `R -e "install.packages(c('limma', 'fgsea'))"`

## Usage
```bash
# Run differential expression analysis
python alzheimers_differential_expression_analysis.py

# Run cell type and pathway analysis
python alzheimers_cell_type_and_pathway_analysis.py

# Run GSEA
Rscript simple-gsea.R
```

## Limitations
- Small sample size
- Potential batch effects
- Limited to temporal cortex samples

## Future Directions
- Validate findings in larger cohorts
- Integrate multi-omics data
- Functional validation of key genes
