import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from scipy.spatial import distance

# Load data
fpkm = pd.read_csv("fpkm_table_unnormalized.csv", index_col=0)
samples = pd.read_csv("columns-samples.csv")
donors = pd.read_csv("DonorInformation.csv")

# Merge metadata and ensure rnaseq_profile_id is a string
merged_metadata = pd.merge(samples, donors, on="donor_id")
merged_metadata["rnaseq_profile_id"] = merged_metadata["rnaseq_profile_id"].astype(str)

# Preprocess FPKM
fpkm = fpkm.T
fpkm.index = fpkm.index.astype(str)  # Convert index to string

# Filter genes
min_samples = int(0.2 * len(fpkm))
filtered_fpkm = fpkm.loc[:, (fpkm >= 1).sum(axis=0) >= min_samples]
filtered_fpkm = np.log2(filtered_fpkm + 1)

# Merge with metadata
processed_data = pd.merge(
    merged_metadata.set_index('rnaseq_profile_id'),
    filtered_fpkm,
    left_index=True,
    right_index=True
)

print(f"Original genes: {fpkm.shape[1]}")
print(f"After filtering: {filtered_fpkm.shape[1]}")
print("Final matrix:", processed_data.shape)

from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.packages import importr
import rpy2.robjects as ro

# Activate R-Python conversion
pandas2ri.activate()

# Subset to Temporal Cortex
tcx_data = processed_data[processed_data['structure_acronym'] == 'TCx']

print(f"Samples in TCx: {len(tcx_data)}")
print(tcx_data['dsm_iv_clinical_diagnosis'].value_counts())

# Add RIN Scores
tbi_metrics = pd.read_csv("tbi_data_files.csv")[['rnaseq_profile_id', 'rna_integrity_number']]
tbi_metrics['rnaseq_profile_id'] = tbi_metrics['rnaseq_profile_id'].astype(str)

tcx_data = pd.merge(
    tcx_data.reset_index(),
    tbi_metrics,
    left_on='rnaseq_profile_id',
    right_on='rnaseq_profile_id',
    how='left'
).set_index('rnaseq_profile_id')

# Drop samples with missing covariates
tcx_data = tcx_data.dropna(subset=['dsm_iv_clinical_diagnosis', 'age', 'sex', 'rna_integrity_number'])

# Create binary variable for Alzheimer's vs No Dementia
# This simplifies our model and gives us a clearer contrast
tcx_data = tcx_data[tcx_data['dsm_iv_clinical_diagnosis'].isin(['Alzheimer\'s Disease Type', 'No Dementia'])]
tcx_data['is_alzheimers'] = (tcx_data['dsm_iv_clinical_diagnosis'] == 'Alzheimer\'s Disease Type').astype(int)

print("Sample counts after filtering:")
print(tcx_data['dsm_iv_clinical_diagnosis'].value_counts())

import patsy

# Create design matrix with the binary variable
design_formula = "~ is_alzheimers + age + sex + rna_integrity_number"
design_matrix = patsy.dmatrix(design_formula, data=tcx_data.reset_index(), return_type='dataframe')

# Verify design matrix
print("Design matrix columns:", design_matrix.columns)

# Prepare expression matrix
metadata_cols = tcx_data.columns[:tcx_data.columns.get_loc(list(filtered_fpkm.columns)[0])]
expression_matrix = tcx_data.drop(columns=metadata_cols).T

# Import R packages
limma = importr('limma')
stats = importr('stats')
base = importr('base')

# Convert to R objects with proper type handling
expression_matrix = expression_matrix.apply(pd.to_numeric, errors='coerce')
expression_matrix = expression_matrix.fillna(0)

print("Converting data to R objects...")
with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
    r_expr = ro.conversion.py2rpy(expression_matrix)
    r_design = ro.conversion.py2rpy(design_matrix)

# Voom normalization
print("Performing voom normalization...")
voom_data = limma.voom(r_expr, r_design, plot=False)

# Fit linear model
print("Fitting linear model...")
fit = limma.lmFit(voom_data, r_design)
fit = limma.eBayes(fit)

# Extract results directly for Alzheimer's coefficient (is_alzheimers)
print("Extracting results...")
results = limma.topTable(
    fit, 
    coef="is_alzheimers",  # Use the binary variable coefficient directly
    number=float('inf'), 
    adjust="BH"  # Benjamini-Hochberg correction
)

# Convert back to pandas
with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
    de_results = ro.conversion.rpy2py(results)

# Add gene symbols
gene_annot = pd.read_csv("row-genes.csv")

# Convert both to strings to ensure they can merge properly
de_results.index = de_results.index.astype(str)
gene_annot['gene_id'] = gene_annot['gene_id'].astype(str)

# Merge results with gene annotations
de_results = pd.merge(
    de_results.reset_index().rename(columns={'index': 'gene_id'}),
    gene_annot,
    on='gene_id',
    how='left'
)

# Save results to CSV
de_results.to_csv("alz_vs_healthy_tcx_results.csv")

# Count excluding any non-gene entries (like is_alzheimers)
de_results_genes_only = de_results[de_results['gene_id'] != 'is_alzheimers'] if 'is_alzheimers' in de_results['gene_id'].values else de_results

# Print FDR < 5% count
sig_genes_strict = sum(de_results_genes_only['adj.P.Val'] < 0.05)
print(f"\nFound {sig_genes_strict} significant genes at FDR < 5%")

# Print FDR < 10% count (more relaxed threshold)
sig_genes_relaxed = sum(de_results['adj.P.Val'] < 0.1)
print(f"Found {sig_genes_relaxed} significant genes at FDR < 10%")

# Print nominal p-value < 0.05 count (without multiple testing correction)
nom_sig_genes = sum(de_results['P.Value'] < 0.05)
print(f"Found {nom_sig_genes} genes with nominal p-value < 0.05 (without correction)")

# Plot p-value distribution to check for potential signal
plt.figure(figsize=(10, 6))
plt.hist(de_results['P.Value'], bins=50)
plt.title("Distribution of p-values")
plt.xlabel("p-value")
plt.ylabel("Frequency")
plt.savefig("pvalue_distribution.png")
print("P-value distribution plot saved to 'pvalue_distribution.png'")

# Plot log fold change distribution
plt.figure(figsize=(10, 6))
plt.hist(de_results['logFC'], bins=50)
plt.title("Distribution of log2 fold changes")
plt.xlabel("log2FC")
plt.ylabel("Frequency")
plt.savefig("logFC_distribution.png")
print("LogFC distribution plot saved to 'logFC_distribution.png'")

# Create a volcano plot
plt.figure(figsize=(10, 8))
plt.scatter(
    de_results['logFC'], 
    -np.log10(de_results['P.Value']),
    alpha=0.6, 
    s=4, 
    color='grey'
)

# Highlight significant genes
sig_genes = de_results[de_results['adj.P.Val'] < 0.05]
if len(sig_genes) > 0:
    plt.scatter(
        sig_genes['logFC'],
        -np.log10(sig_genes['P.Value']),
        alpha=1,
        s=25,
        color='red'
    )
    
    # Label significant genes
    for idx, row in sig_genes.iterrows():
        plt.annotate(
            row['gene_symbol'],
            (row['logFC'], -np.log10(row['P.Value'])),
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=8
        )

plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.3)
plt.axvline(x=0, color='black', linestyle='-', alpha=0.3)

# Add title and labels
plt.title("Volcano Plot: Alzheimer's vs No Dementia")
plt.xlabel("log2 Fold Change")
plt.ylabel("-log10 P-value")
plt.savefig("volcano_plot.png")
print("Volcano plot saved to 'volcano_plot.png'")

# Display top genes regardless of significance
print("\nTop 20 genes by p-value:")
top_genes = de_results.sort_values('P.Value').head(20)[['gene_id', 'gene_symbol', 'logFC', 'P.Value', 'adj.P.Val']]
print(top_genes)

# Analyze the significant genes in more detail
if sig_genes_strict > 0:
    print("\nSignificant genes (FDR < 5%):")
    sig_gene_details = de_results[de_results['adj.P.Val'] < 0.05][['gene_id', 'gene_symbol', 'gene_name', 'logFC', 'P.Value', 'adj.P.Val']]
    print(sig_gene_details)

# Check for key Alzheimer's related genes
print("\nChecking for key Alzheimer's-related genes...")
key_genes = ['APP', 'PSEN1', 'PSEN2', 'MAPT', 'APOE', 'TREM2', 'CLU', 'ABCA7', 'BIN1', 'CD33', 'MS4A', 'PICALM', 'SORL1']
ad_related_genes = de_results[de_results['gene_symbol'].isin(key_genes)].sort_values('P.Value')
if not ad_related_genes.empty:
    print(f"Found {len(ad_related_genes)} Alzheimer's-related genes in results:")
    print(ad_related_genes[['gene_symbol', 'logFC', 'P.Value', 'adj.P.Val']])
else:
    print("No known Alzheimer's-related genes found in our results.")

# Check for enrichment in top genes (looking at top 100 by p-value)
print("\nChecking for functional enrichment in top genes...")
top_100_genes = de_results.sort_values('P.Value').head(100)
# Save the top genes to a file for potential external enrichment analysis
top_100_genes[['gene_id', 'gene_symbol', 'logFC', 'P.Value', 'adj.P.Val']].to_csv("top_100_genes.csv")
print("Top 100 genes saved to 'top_100_genes.csv' for enrichment analysis")

# Try an alternative model with no age term (might be collinear with diagnosis)
print("\nTrying alternative model without age term...")
design_formula_noage = "~ is_alzheimers + sex + rna_integrity_number"
design_matrix_noage = patsy.dmatrix(design_formula_noage, data=tcx_data.reset_index(), return_type='dataframe')

with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
    r_design_noage = ro.conversion.py2rpy(design_matrix_noage)

# Fit model 
fit2 = limma.lmFit(voom_data, r_design_noage)
fit2 = limma.eBayes(fit2)

# Extract results
results2 = limma.topTable(fit2, coef="is_alzheimers", number=float('inf'), adjust="BH")

# Convert back to pandas
with ro.conversion.localconverter(ro.default_converter + pandas2ri.converter):
    de_results2 = ro.conversion.rpy2py(results2)

# Process results
de_results2.index = de_results2.index.astype(str)
de_results2 = pd.merge(
    de_results2.reset_index().rename(columns={'index': 'gene_id'}),
    gene_annot,
    on='gene_id',
    how='left'
)

# Save results from the simpler model
de_results2.to_csv("alz_vs_healthy_tcx_noage_results.csv")

# Print statistics
sig_genes_noage = sum(de_results2['adj.P.Val'] < 0.05)
print(f"\nWithout age term: Found {sig_genes_noage} significant genes at FDR < 5%")

# Display top genes from simpler model
print("\nTop 20 genes without age term:")
top_genes_noage = de_results2.sort_values('P.Value').head(20)[['gene_id', 'gene_symbol', 'logFC', 'P.Value', 'adj.P.Val']]
print(top_genes_noage)

# Compare top genes between original and no-age models
common_genes = set(top_genes['gene_id']).intersection(set(top_genes_noage['gene_id']))
print(f"\nNumber of genes in common between top 20 of both models: {len(common_genes)}")

# Try Gene Set Enrichment Analysis (GSEA) approach using fgsea package in R if available
try:
    print("\nAttempting gene set enrichment analysis...")
    # Try to load fgsea package
    try:
        fgsea = importr('fgsea')
        stats_available = True
    except:
        print("fgsea package not available. Skipping GSEA.")
        stats_available = False
    
    if stats_available:
        # Prepare ranked gene list for GSEA (all genes ranked by -log10(p) * sign(logFC))
        ranked_genes = de_results.copy()
        
        # Filter out rows with missing gene symbols or the 'is_alzheimers' row
        ranked_genes = ranked_genes[ranked_genes['gene_symbol'].notna() & 
                                    (ranked_genes['gene_id'] != 'is_alzheimers')]
        
        # Calculate ranking score
        ranked_genes['score'] = -np.log10(ranked_genes['P.Value']) * np.sign(ranked_genes['logFC'])
        
        # Sort by score
        ranked_genes = ranked_genes.sort_values('score', ascending=False)
        
        # Handle potential duplicate gene symbols (use most significant if duplicates exist)
        ranked_genes = ranked_genes.loc[~ranked_genes['gene_symbol'].duplicated()]
        
        # Convert to a named vector for R
        print(f"Creating ranked list of {len(ranked_genes)} genes...")
        
        # Use a dict first, then convert to R vector
        gene_dict = dict(zip(ranked_genes['gene_symbol'].astype(str), ranked_genes['score']))
        
        # Convert to R named vector
        from rpy2.robjects import vectors
        gene_ranks = vectors.FloatVector(list(gene_dict.values()))
        gene_ranks.names = vectors.StrVector(list(gene_dict.keys()))
        
        # Print example of the prepared gene ranks
        print(f"Top 5 genes in ranking: {list(gene_dict.items())[:5]}")
        
        # Simple GSEA with GO terms or KEGG pathways
        print("GSEA analysis requires biological pathway definitions.")
        print("Consider using external tools like DAVID, EnrichR, or WebGestalt with the exported gene lists.")
        
        # Export the ranked gene list for external tools
        ranked_genes[['gene_symbol', 'score']].to_csv("ranked_genes_for_gsea.csv")
        print("Ranked gene list exported to 'ranked_genes_for_gsea.csv' for external GSEA tools")
        
except Exception as e:
    print(f"Error performing gene set analysis: {e}")
    print("Consider using external tools like DAVID, EnrichR, or WebGestalt with the exported gene lists.")