#!/usr/bin/env python3

import os
import argparse
import numpy as np
import pandas as pd
import scipy.io
import anndata as ad
import celltypist


def load_unassigned_folder(base_dir: str) -> ad.AnnData:
    # Expecting base_dir to contain unassigned_logNorm.mtx, genes.txt, cells.txt    
    # These are the default outputs from Seurat's "write10x" function when exporting unassigned cells.
    #Prepare the paths to the required files
    mtx_path = os.path.join(base_dir, "unassigned_logNorm.mtx")
    genes_path = os.path.join(base_dir, "genes.txt")
    cells_path = os.path.join(base_dir, "cells.txt")

    missing = [p for p in [mtx_path, genes_path, cells_path] if not os.path.exists(p)] # Veriify that all required files exist
    if missing:
        raise FileNotFoundError("Missing required files:\n" + "\n".join(missing)) # If any files are missing, raise an error with a message listing the missing files

    X = scipy.io.mmread(mtx_path).tocsr()  # genes x cells # Load the matrix in Matrix Market format and convert it to Compressed Sparse Row format for efficient row slicing
    genes = [g.strip() for g in open(genes_path, encoding="utf-8")] # Read the genes from the genes.txt file, stripping
    cells = [c.strip() for c in open(cells_path, encoding="utf-8")] # Read the cell barcodes from the cells.txt file

    if X.shape != (len(genes), len(cells)): # Make sure the dimensions of the matrix match the number of genes and cells
        raise ValueError(f"Shape mismatch: MTX={X.shape}, genes={len(genes)}, cells={len(cells)}")

    adata = ad.AnnData(X=X.T)  # cells x genes , transpose to have cells as rows and genes as columns, which is the expected format for CellTypist  (Celltypist expects cells as rows and genes as columns, so we transpose the matrix here)
    adata.var_names = genes # Set the variable names (genes) in the AnnData object
    adata.obs_names = cells # Set the observation names (cells) in the AnnData object
    return adata # Return the AnnData object containing the unassigned cells and their gene expression data


def pick_final_labels(pred, use_majority_voting: bool) -> pd.Series:
    """
    Return final 1D Series of labels for each cell.
    """
    # CellTypist's predicted_labels can be either a Series or a DataFrame (if majority voting is used). We want to pick the appropriate column if it's a DataFrame.
    labs = pred.predicted_labels # This can be a Series (if no majority voting) or a DataFrame (if majority voting). We need to pick the right column based on the user's choice.
    if isinstance(labs, pd.DataFrame): # If it's a DataFrame?, we check if the user wants majority voting and if the "majority_voting" column exists. If so, we use that. Otherwise, we fall back to "predicted_labels" or the first column.
        if use_majority_voting and "majority_voting" in labs.columns: # If the user requested majority voting and the column exists, use it.
            return labs["majority_voting"].astype(str)
        if "predicted_labels" in labs.columns: # If the "predicted_labels" column exists, use it as a fallback.
            return labs["predicted_labels"].astype(str)
        return labs.iloc[:, 0].astype(str) # If neither column exists, just take the first column of the DataFrame as a last resort.
    return labs.astype(str)


def compute_confidence(prob: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """
    prob : pred.probability_matrix , where there is the predicted probability for each cell (rows) and each cell type (columns).
    Compute confidence metrics from the probability matrix. Returns (max_prob, margin_top1_top2).
    -- max_prob: the highest predicted probability for each cell.
    -- margin_top1_top2: the difference between the top predicted probability and the second highest,
    """
    
    top1 = prob.max(axis=1) #the max predicted probability for each cell (the confidence of the top prediction)
    top2 = prob.apply(lambda r: np.partition(r.values, -2)[-2], axis=1) # the second highest predicted probability for each cell (the confidence of the second-best prediction)
    margin = top1 - top2 # the margin between the top prediction and the second-best prediction, which can be a measure of confidence
    return top1, margin # Return both the max probability and the margin as confidence metrics for each cell's prediction

def simplify_label(label: str) -> str:
    """
    Map detailed CellTypist labels to broader categories for easier interpretation with my previous annotations. 
    """
    s = str(label) # Convert to string in case it's not already (e.g., if it's a category or something else)

    # cortical excitatory layers
    if s.startswith(("L2", "L3", "L4", "L5", "L6")):
        return "Excitatory Neurons"

    # inhibitory neurons (in this model)
    if s.startswith("InN") or s.startswith("IN"):
        return "Inhibitory Neurons"

    # glia / other
    if s.startswith("Astro"):
        return "Astrocyte"
    if s.startswith("Micro"):
        return "Microglia"
    if s.startswith("Oligo"):
        return "Oligodendrocytes"
    if s.startswith("OPC"):
        return "OPCs"
    if s.startswith("Endo"):
        return "Endothelial"

    # immune broader buckets
    if s.startswith("Macro") or s.startswith("Myeloid") or s.startswith("DC") or s.startswith("B ") or s.startswith("T ") or s.startswith("NK"):
        return "Immune"

    return "Other/Unknown"


def main():
    ap = argparse.ArgumentParser(
        description="Run CellTypist on Seurat-exported unassigned cells (unassigned_logNorm.mtx/genes.txt/cells.txt)."
    )
    ap.add_argument("--base", required=True, help="Folder with unassigned_logNorm.mtx, genes.txt, cells.txt")
    ap.add_argument("--model", required=True, help="Path to CellTypist .pkl model (or model name)")
    ap.add_argument("--majority_voting", action="store_true", help="Use majority voting (recomendation ") # If the model supports majority voting,
    #this will use the majority_voting column for final labels instead of the default predicted_labels. 
    # This can provide more robust predictions by considering multiple models or iterations, but it requires that the model was trained with majority voting enabled.
    # If you want to use majority voting : False  is the default , 
    ap.add_argument("--min_prob", type=float, default=0.0,help="If >0, set labels to 'Unassigned' when max_prob < min_prob") # This allows the user to set a confidence threshold for the predictions. 
    # If the maximum predicted probability for a cell is below this threshold, the cell will be labeled as "Unassigned"
    # instead of being assigned to a specific cell type. This can help improve the reliability of the annotations by filtering out low-confidence predictions.
    ap.add_argument("--out_csv", default=None,help="Output CSV path. Default: <base>/unassigned_celltypist_predictions.csv") 
    args = ap.parse_args()

    adata = load_unassigned_folder(args.base) # Load the unassigned cells from the specified folder into an AnnData object. 
    #This function reads the matrix, genes, and cells files and constructs the AnnData object in the format expected by CellTypist.

    pred = celltypist.annotate(  # RUN CellTypist
        adata, #mi adata 
        model=args.model, #  a path to a .pkl file containing a trained model, 
        majority_voting=args.majority_voting #default False, if True, it will use the majority_voting column 
    )

    prob = pred.probability_matrix  # cells x types matrix of predicted probabilities for each cell and each cell type
    top1, margin = compute_confidence(prob)  # Compute confidence metrics (max probability and margin between top 2 predictions)

    labels = pick_final_labels(pred, use_majority_voting=args.majority_voting) # Get the final predicted labels for each celL
    labels = labels.reindex(prob.index).astype(str)  # Make sure the labels are in the same order as the probability matrix and convert to string

    simple = labels.map(simplify_label) #Aplly the function to simplify the labels 
    # optional filtering by confidence
    if args.min_prob > 0:
        low = top1 < args.min_prob # the probability we assigned in the argumnets, 
        #if the max probability for a cell is below this threshold, we consider it low confidence
        #.loc is used to selecct those cells with low confidence in a pandas series object
        labels.loc[low] = "Unassigned" # For low-confidence cells, we set the label to "Unassigned" instead of the predicted cell type
        simple.loc[low] = "Unassigned" # For low-confidence cells, we set the label to "Unassigned"

    out = pd.DataFrame({
        "cell": prob.index.astype(str),
        "celltypist_label": labels.values,
        "celltype_simple": simple.values,
        "max_prob": top1.values,
        "margin_top1_top2": margin.values
    })

    out_csv = args.out_csv or os.path.join(args.base, "unassigned_celltypist_predictions.csv")
    out.to_csv(out_csv, index=False)

    print("Saved:", out_csv)
    print("Coarse counts:\n", out["celltype_simple"].value_counts().head(20).to_string())


if __name__ == "__main__":
    main()
