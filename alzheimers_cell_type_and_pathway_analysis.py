import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load your DE results
de_results = pd.read_csv("alz_vs_healthy_tcx_results.csv")

# Load cell type expression data
cell_type_expr = pd.read_csv("trimmed_means.csv")
cell_type_expr = cell_type_expr.set_index('gene_symbol')

# Print some basic stats
print(f"Total genes in cell type data: {len(cell_type_expr)}")
print(f"Total cell types: {len(cell_type_expr.columns)}")

# Check for SLC6A9 specifically - with case insensitive search
slc6a9_candidates = [gene for gene in cell_type_expr.index if 'slc6a9' in gene.lower()]
print(f"Possible SLC6A9 gene matches: {slc6a9_candidates}")

# Try with Ensembl ID or other common gene ID formats if needed
# Sometimes gene symbol mappings vary between datasets

# Get top DE genes (significant or top 50 by p-value)
top_genes = de_results.sort_values('P.Value').head(50)

# Create mapping between your gene IDs and the Allen dataset
# You might need to use a mapping file or try different formats
gene_mapping = {}
for gene_id in top_genes['gene_id']:
    # Extract gene symbol if it's in your results
    if 'gene_symbol' in top_genes.columns:
        symbol = top_genes.loc[top_genes['gene_id'] == gene_id, 'gene_symbol'].values[0]
        if pd.notna(symbol) and symbol in cell_type_expr.index:
            gene_mapping[gene_id] = symbol

print(f"Successfully mapped {len(gene_mapping)} genes to the Allen dataset")

# Analyze cell type specificity for mapped genes
mapped_genes = list(gene_mapping.values())
if mapped_genes:
    # Extract expression data for these genes
    gene_expr_by_celltype = cell_type_expr.loc[mapped_genes]
    
    # Group cell types into major categories
    cell_type_categories = {
        'GABAergic': [col for col in cell_type_expr.columns if any(x in col for x in ['Lamp5', 'Pax6', 'Sncg', 'Vip', 'Sst', 'Pvalb', 'Chandelier'])],
        'Glutamatergic': [col for col in cell_type_expr.columns if any(x in col for x in ['L2.3.IT', 'L4.IT', 'L5.IT', 'L5.ET', 'L6.CT', 'L6b', 'L6.IT', 'L5.6.NP'])],
        'Astrocytes': [col for col in cell_type_expr.columns if 'Astro' in col],
        'Oligodendrocytes': [col for col in cell_type_expr.columns if 'Oligo' in col or 'OPC' in col],
        'Microglia': [col for col in cell_type_expr.columns if 'Micro' in col],
        'Vascular': [col for col in cell_type_expr.columns if any(x in col for x in ['Endo', 'VLMC'])]
    }
    
    # Calculate average expression in each major cell type category
    category_avg_expr = pd.DataFrame(index=mapped_genes)
    for category, cols in cell_type_categories.items():
        if cols:  # Only process non-empty categories
            category_avg_expr[category] = gene_expr_by_celltype[cols].mean(axis=1)
    
    # Create heatmap of average expression by major cell type
    plt.figure(figsize=(10, len(mapped_genes) * 0.3 + 2))
    sns.heatmap(category_avg_expr, cmap="viridis", yticklabels=True, annot=True, fmt=".2f")
    plt.title("Average Expression of DE Genes by Cell Type Category")
    plt.tight_layout()
    plt.savefig("de_genes_cell_type_category_heatmap.png", dpi=300)
    
    # For each gene, determine which cell type has highest expression
    gene_cell_type_enrichment = {}
    for gene in mapped_genes:
        if gene in gene_expr_by_celltype.index:
            expr_across_types = gene_expr_by_celltype.loc[gene]
            max_cell_type = expr_across_types.idxmax()
            max_value = expr_across_types.max()
            gene_cell_type_enrichment[gene] = (max_cell_type, max_value)
            
    # Print cell type enrichment for top genes
    print("\nCell type enrichment for top DE genes:")
    for gene, (cell_type, value) in gene_cell_type_enrichment.items():
        print(f"{gene}: Highest in {cell_type} (expression = {value:.2f})")

# For pathway analysis with GSEA results
try:
    # Load GSEA results - check which file format you have
    gsea_go = pd.read_csv("gsea_results_go.csv")
    
    # Check columns to ensure we're using the right names
    print("\nGSEA GO result columns:", gsea_go.columns.tolist())
    
    # Look at the first few pathway names to see their format
    print("\nExample pathway names:")
    if 'pathway_name' in gsea_go.columns:
        print(gsea_go['pathway_name'].head())
    elif 'description' in gsea_go.columns:
        print(gsea_go['description'].head())
    
    # Fix the pathway mapping using the correct column
    pathway_col = 'description' if 'description' in gsea_go.columns else 'pathway_name'
    
    # Select significant pathways
    sig_pathways = gsea_go[gsea_go['padj'] < 0.1].copy()
    
    # Improved keywords for cellular processes
    cell_type_pathways = {
        'Neuronal': ['neuron', 'synap', 'axon', 'dendrit', 'neurotrans', 'ion channel', 'nerve'],
        'Astrocytes': ['astro', 'glutamate uptake', 'potassium ion trans'],
        'Oligodendrocytes': ['myelin', 'ensheath', 'oligodendro'],
        'Microglia': ['immun', 'inflamm', 'cytokine', 'phago', 'microglia'],
        'RNA Processing': ['rna', 'splicing', 'transcription', 'mrna'],
        'Translation': ['translation', 'ribosom', 'protein synth'],
        'Cell Cycle': ['cell cycle', 'mitosis', 'replication']
    }
    
    # Tag pathways by process type (fix the warning by using .loc)
    for cell_type, keywords in cell_type_pathways.items():
        sig_pathways.loc[:, f'{cell_type}_related'] = sig_pathways[pathway_col].astype(str).str.lower().apply(
            lambda x: any(kw in x.lower() for kw in keywords)
        )
    
    # Count pathways by process
    process_counts = {ct: sig_pathways[f'{ct}_related'].sum() for ct in cell_type_pathways.keys()}
    print("\nNumber of pathways related to each process:")
    for process, count in process_counts.items():
        print(f"{process}: {count}")
    
    # Filter out processes with zero counts
    filtered_counts = {k: v for k, v in process_counts.items() if v > 0}
    
    # Create pie chart if we have matched pathways
    if sum(filtered_counts.values()) > 0:
        plt.figure(figsize=(10, 10))
        plt.pie(filtered_counts.values(), labels=filtered_counts.keys(), autopct='%1.1f%%')
        plt.title("Distribution of Significant Pathways by Process")
        plt.savefig("pathway_process_distribution.png", dpi=300)
        print("Pathway process distribution saved to pathway_process_distribution.png")
    else:
        print("No matching pathways found for the defined categories")
        
except Exception as e:
    print(f"Error in pathway analysis: {e}")