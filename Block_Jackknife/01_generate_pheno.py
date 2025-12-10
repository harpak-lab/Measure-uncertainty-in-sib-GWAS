import pandas as pd
import numpy as np
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Generate BJK phenotype files")
    parser.add_argument("--pheno_name", required=True, help="Name of the phenotype (e.g., height)")
    parser.add_argument("--input_file", required=True, help="Path to residual phenotype file")
    parser.add_argument("--out_dir", required=True, help="Directory to save output")
    parser.add_argument("--seed", type=int, default=8675309, help="Random seed")
    args = parser.parse_args()

    np.random.seed(args.seed)
    
    # Load Data
    phenotype_file = pd.read_csv(args.input_file, sep='\t', header=None, names=['FID','IID','pheno'])
    sibpairs = phenotype_file['FID'].unique()

    # Generate 500 permutations
    for iteration in range(1, 501):
        col_name = f'pheno{iteration}'
        phenotype_file[col_name] = phenotype_file['pheno']
        
        # Sample 500 families without replacement
        omission_fams = np.random.choice(sibpairs, 500, replace=False)
        
        # Set omitted families to NA
        phenotype_file.loc[phenotype_file['FID'].isin(omission_fams), col_name] = np.nan

    # Save
    phenotype_file = phenotype_file.drop(columns=['pheno'])
    out_path = os.path.join(args.out_dir, f'{args.pheno_name}.blockjackknife.phenos.txt')
    phenotype_file.to_csv(out_path, sep='\t', index=False, na_rep='NA')
    print(f"Saved to {out_path}")

if __name__ == "__main__":
    main()
