import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import warnings
import os
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
warnings.filterwarnings("ignore")
from tqdm import tqdm
from PIL import Image
import random
import h5py
import cv2

from st_helper import *
Image.MAX_IMAGE_PIXELS = 933120000
# sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 1 #3

def plot_unprocessed(st_dir, img_dir, df_img_name, img_save_dir):
    """
    Here we do not do any processing on the adata/images, 
    we just plot the distribution of gene expression transcripts and spatial overlays on H&E images
    to get a sense of the quality of the data and the gene filter cut-off values (min_counts, max_counts, min_cells)
    """

    print("Plotting plots on unprocessed adata...")
    os.makedirs(img_save_dir, exist_ok=True)
    
    matched_dict = get_matched_dict_from_df(df_img_name, st_dir, img_dir)
    print(f'Number of samples: {len(matched_dict)}')
    sample_names = list(matched_dict.keys())
    st_path_list = [val[0] for val in matched_dict.values()]

    for i, path in enumerate(st_path_list):
        adata = sc.read_visium(
            path, 
            count_file='filtered_feature_bc_matrix.h5', 
            source_image_path='tissue_hires_image.png'
        )
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        print(adata)

        ### 1. Plot Distribution ###
        fig, axs = plt.subplots(1, 4, figsize=(20, 5))
        sns.distplot(adata.obs["total_counts"], kde=False, ax=axs[0])
        # Plot total counts with a threshold
        sns.distplot(adata.obs["total_counts"][adata.obs["total_counts"] < 35000], kde=False, bins=40, ax=axs[1])
        # Plot number of genes by counts for all
        sns.distplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
        # Plot number of genes by counts with a threshold
        sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 6000], kde=False, bins=60, ax=axs[3])
        
        # Adding titles to the leftmost subplots for each sample
        axs[0].set_ylabel(f"Sample {i}\n{sample_names[i]}", rotation=0, labelpad=100, size='large', verticalalignment='center')

        plt.tight_layout()
        plt.savefig(os.path.join(img_save_dir, f"{i}_{sample_names[i]}_dist_plots.png"))
        plt.close()
        
        ### 2. Plot Spatial overlays on H&E Images ######
 
        # Generate spatial plots without showing them
        sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts", "pct_counts_in_top_200_genes"],
                    ncols=3, cmap='plasma', alpha_img=0.5, 
                    title=[f"Sample {i}, {sample_names[i]}, total_counts", "n_genes_by_counts", "pct_counts_in_top_200_genes"], show=False)
        
        # Adjust the layout to make room for the custom titles
        plt.tight_layout()
        # Save the figure
        plt.savefig(os.path.join(img_save_dir, f"{i}_{sample_names[i]}_spatial_plots.png"))
        plt.close()
    print(f"Distribution plots saved in {img_save_dir}")


def load_data(st_dir, img_dir, df_img_name):

    print("Loading data...")
    # Load the data and filter out bad samples
    matched_dict = get_matched_dict_from_df(df_img_name, st_dir, img_dir)

    print(f'Number of samples before filtering: {len(matched_dict)}')
    bad_samples = []
    print(f'Number of bad samples to remove: {len(bad_samples)}')
    # get the matched_dict without the bad_samples
    matched_dict = {k: v for k, v in matched_dict.items() if k not in bad_samples}
    sample_names = list(matched_dict.keys())
    print(f'Number of samples after filtering: {len(matched_dict)}')

    # use a dataframe to store the individaul filtering parameters: sample_name, min_counts, max_counts, min_cells, pct_counts_mt for each sample
    # here we only do basic filtering, more advanced filtering can be done in the future
    df_filter = pd.DataFrame({'sample_name': list(matched_dict.keys()), 
                            'min_counts': [5000]*len(matched_dict), 
                            'max_counts': [35000]*len(matched_dict), 
                            'min_cells': [50]*len(matched_dict), 
                            'pct_counts_mt': [20]*len(matched_dict),
                            'cv_threshold':[210]*len(matched_dict) # cv_threshold is used for computer vision filtering to discard largely empty images
                            })

    adata_list, img_list, hvgs_union = get_data_lists(matched_dict, 
                                                        hvg_list = [], 
                                                        num_hvg = 800, 
                                                        df_filter = df_filter,
                                                        log_norm = False, 
                                                        # min_counts = 5000, max_counts = 35000, 
                                                        # min_cells = 50,
                                                        # pct_counts_mt = 20, 
                                                        )
    print("Loading data complete")
    return adata_list, img_list, hvgs_union, sample_names, matched_dict

###############################################################
#### 1. Plot Distribution of Gene Expression Transcripts ######
# ax1: total_counts
# ax2: n_genes_by_counts
###############################################################
# This can tell the optimal cutoffs for filtering out bad cells (min_counts, max_counts, min_cells)
# and if any abnormal samples need to be removed


def plot_dist_plots(adata_list, sample_names, img_save_dir, start, end):
    print("Plotting distribution plots...")
    os.makedirs(img_save_dir, exist_ok=True)

    # Creating sub-lists of 5 samples each
    for i in range(start, end, 5):
        # Slice the adata_list to get sub-lists of 5 elements each
        sub_list = adata_list[i:i+5]

        # Create a large figure to hold all subplots
        fig, axs = plt.subplots(len(sub_list), 2, figsize=(12, 5*len(sub_list)))  # Adjust the size as needed

        for adata_index, adata in enumerate(sub_list):
            # Plot total counts for all
            sns.distplot(adata.obs["total_counts"], kde=False, ax=axs[adata_index, 0])
            # Plot total counts with a threshold
            #sns.distplot(adata.obs["total_counts"][adata.obs["total_counts"] < 35000], kde=False, bins=40, ax=axs[adata_index, 1])
            # Plot number of genes by counts for all
            sns.distplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[adata_index, 1])
            # Plot number of genes by counts with a threshold
            #sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 6000], kde=False, bins=60, ax=axs[adata_index, 3])
            
            # Adding titles to the leftmost subplots for each sample
            axs[adata_index, 0].set_ylabel(f"Sample { i + adata_index + 1}\n{sample_names[i + adata_index]}", rotation=0, labelpad=100, size='large', verticalalignment='center')

        # Adjust layout for better appearance
        plt.tight_layout()
        # Save the entire figure
        plt.savefig(os.path.join(img_save_dir, f"{i}_{i+len(sub_list)-1}_combined_samples_dist_plots.png"))
        plt.close()  # Close the plot to free memory
    print(f"Distribution plots saved in {img_save_dir}")


###############################################################
#### 2. Plot Spatial overlays on H&E Images ######
# ax1: total_counts
# ax2: n_genes_by_counts
# ax3: pct_counts_in_top_200_genes
###############################################################

def plot_spatial_plots(adata_list, sample_names, img_save_dir, start, end):
    print("Plotting spatial plots...")
    os.makedirs(img_save_dir, exist_ok=True)

    for i, adata in enumerate(tqdm(adata_list[start:end])):
        sample_name = sample_names[start+i]
        #print(sample_name)
        
        # Generate spatial plots without showing them
        sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts", "pct_counts_in_top_200_genes"],
                    ncols=3, cmap='plasma', alpha_img=0.5, title=[f"Sample {start+i}, {sample_name}, total_counts", "n_genes_by_counts", "pct_counts_in_top_200_genes"], show=False)
        
        # Adjust the layout to make room for the custom titles
        plt.tight_layout()
        
        # Save the figure
        plt.savefig(os.path.join(img_save_dir, f"{i}_{sample_name}_spatial_plots.png"))
        plt.close()  # Close the plot to free memory
    print(f"H&E overlay spatial plots saved in {img_save_dir}")




def main():

    st_dir = '/home/mahmoodlab8/workspace/anthony/projects/ST_codebase/data/lymph_23/st_files'
    img_dir = '/home/mahmoodlab8/workspace/anthony/projects/ST_codebase/data/lymph_23/images'
    df_img_name = pd.read_csv('/home/mahmoodlab8/workspace/anthony/projects/ST_codebase/data/lymph_23/df_img_name.csv')
    """
            sample_name	img_name
        0	GSM5924042_frozen_a_1	a_1 Frozen.ndpi
        1	GSM5924044_frozen_a_15	a_15 Frozen.ndpi
        2	GSM5924045_frozen_a_17	a_17 Frozen.ndpi
        3	GSM5924043_frozen_a_3	a_3 Frozen.ndpi
        ...
    """

    plot_unprocessed(st_dir, img_dir, df_img_name, './figures/lymph_23/unprocessed_plots')
    # adata_list, img_list, hvgs_union, sample_names, matched_dict = load_data(st_dir, img_dir, df_img_name)
    # start, end = 0, len(adata_list)
    # # Ensure the directory exists for saving the figure
    # dist_img_save_dir = './figures/lymph_23/combined_dist_plots'
    # spatial_img_save_dir = './figures/lymph_23/spatial_plots'

    # plot_dist_plots(adata_list, sample_names, dist_img_save_dir, start, end)
    #plot_spatial_plots(adata_list, sample_names, spatial_img_save_dir, start, end)


if __name__ == "__main__":
    main()