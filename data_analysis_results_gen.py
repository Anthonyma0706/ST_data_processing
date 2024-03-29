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
import openslide
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

    adata_shape_l = []
    for i, path in enumerate(st_path_list):
        adata = sc.read_visium(
            path, 
            count_file='filtered_feature_bc_matrix.h5', 
            source_image_path='tissue_hires_image.png'
        )
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        adata_shape_l.append(f'{adata.shape}[0] x {adata.shape}[1]')
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
    return adata_shape_l


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


# here we load the raw adata, we do not do any QC filtering now (min_counts, max_counts, min_cells, pct_counts_mt) 
def get_spot_pixel_size_resolution(adata):
    d = list(adata.uns['spatial'].values())[0]
    d_scale_factor = d['scalefactors']
    pix_size = round(d_scale_factor['spot_diameter_fullres'], 2)
    resolution = round(55/ d_scale_factor['spot_diameter_fullres'], 2)
    magnification = '< 5x'
    if resolution < 0.4:
        magnification = '20x'
    elif resolution < 0.8:
        magnification = '10x'

    return pix_size, resolution, magnification

def generate_metrics_df(df_img_name, st_dir, img_dir, load_img = True, cv_thresh = 200, 
                        check_wsi = False,
                        save_path = None):

    matched_dict = get_matched_dict_from_df(df_img_name, st_dir, img_dir)
    print(f'Number of samples: {len(matched_dict)}')
    sample_names = list(matched_dict.keys())
    st_path_list = [val[0] for val in matched_dict.values()]

    df_spatial_l = []
    adata_list = []
    img_list = []
    wsi_list = []
    patches_list = []

    # metrics
    adata_shape_l = []
    pix_size_l = []
    resolution_l = []
    mag_l = []
    img_shape_l = []
    wsi_level0_l = []
    samples_read_l = []

    # pixel coords
    bottom_x_y_l = []
    bottom_x_y_pixel_l = []
    top_x_y_l = []
    top_x_y_pixel_l = []
    max_pixel_l = []
    min_pixel_l = []

    # patch metrics
    bad_patch_percentage_l = []
    good_patch_num_l = []
    bad_patch_num_l = []
    cv_thres_l = []

    for key in tqdm(matched_dict):
        st_path, img_path = matched_dict[key]
        # print(key)
        samples_read_l.append(key)
        adata = sc.read_visium(
            st_path, 
            count_file='filtered_feature_bc_matrix.h5', 
            source_image_path='tissue_hires_image.png'
        )
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        pix_size, resolution, magnification = get_spot_pixel_size_resolution(adata)
        pix_size_l.append(pix_size)
        resolution_l.append(resolution)
        mag_l.append(magnification)
        adata_shape_l.append(f'{adata.shape[0]} x {adata.shape[1]}')

        spatial_coord_path = os.path.join(st_path, "spatial/tissue_positions_list.csv")
        spatial = pd.read_csv(spatial_coord_path, sep="," , na_filter=False, index_col=0, header=None)
        spatial.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
        df_merged = adata.obs.merge(spatial[['pxl_row_in_fullres', 'pxl_col_in_fullres']], 
                                left_index=True, 
                                right_index=True, 
                                how='left')
        # for easy naming in patching functions, while preserving the original names for scanpy
        df_merged['x_array'] = df_merged['array_row']
        df_merged['y_array'] = df_merged['array_col']
        df_merged['x_pixel'] = df_merged['pxl_row_in_fullres']
        df_merged['y_pixel'] = df_merged['pxl_col_in_fullres']

        adata.obs = df_merged
        #print(adata)
        df_spatial_l.append(df_merged)
        adata_list.append(adata)

        # get the spot coordinate at top/bottom of the image and the pixel coordinate, use array_col to locate
        df_merged = df_merged[['array_col', 'array_row', 'pxl_row_in_fullres', 'pxl_col_in_fullres']]
        df_merged = df_merged.sort_values(by=['array_col', 'array_row'])
        # get the first spot when array_col = min
        y, x, x_pixel, y_pixel = df_merged[df_merged['array_col'] == df_merged['array_col'].min()].iloc[0,:4].values
        bottom_x_y_l.append(f'{x} , {y}')
        bottom_x_y_pixel_l.append(f'{x_pixel} , {y_pixel}')

        y, x, x_pixel, y_pixel = df_merged[df_merged['array_col'] == df_merged['array_col'].max()].iloc[0,:4].values
        top_x_y_l.append(f'{x} , {y}')
        top_x_y_pixel_l.append(f'{x_pixel} , {y_pixel}')

        min_pixel_l.append(f'{df_merged["pxl_row_in_fullres"].min()} , {df_merged["pxl_col_in_fullres"].min()}')
        max_pixel_l.append(f'{df_merged["pxl_row_in_fullres"].max()} , {df_merged["pxl_col_in_fullres"].max()}')

        # # convert gene names to captical letters
        # adata.var_names=[name.upper() for name in list(adata.var_names)]
        # adata.var["gene_name"]=adata.var.index.astype("str")

        ### 3. H&E IMAGE
        if load_img:
            img = np.array(Image.open(img_path))
            #print(f'Image shape: {img.shape}')
            img_list.append(img)
            img_shape_l.append(f'{img.shape[0]} x {img.shape[1]} x {img.shape[2]}')
            if check_wsi:
                wsi = openslide.OpenSlide(img_path)
                #print(wsi.level_dimensions)
                wsi_level0_l.append(f'{wsi.level_dimensions[0][0]} x {wsi.level_dimensions[0][1]}')
                wsi_list.append(wsi)
            else:
                wsi_level0_l.append('None')
                wsi_list.append('None')
            
            df = adata.obs[['x_array', 'y_array', 'x_pixel', 'y_pixel']]
            patch_pixel_size = int(pix_size+100) # 20x mag: 142.04 + 100
            try:
                cv_thresh = 200 # 160 turns out to be a good threshold that can filter 'blue spot', 200 is not enough
                img_patches_interest = extract_image_patches(df, img, 
                                                    patch_size = patch_pixel_size,
                                                    )
                bad_patch_indices = identify_bad_patches(img_patches_interest, threshold=cv_thresh) 
                bad_patches_l = list(np.array(img_patches_interest)[bad_patch_indices])
                # Convert bad_patch_indices to a set for faster lookup
                bad_patch_indices_set = set(bad_patch_indices)
                # Use a list comprehension to filter out the bad patches
                good_patches_l = [patch for j, patch in enumerate(img_patches_interest) if j not in bad_patch_indices_set]
                
                #print(f'i = {i}, {sample_names[i]}')
                # print(f'Threshold = {cv_thresh}')
                # print(f'# Good patches: {len(good_patches_l)}, # Bad patches: {len(bad_patches_l)}')
                # print(f'Percentage of bad patches: {len(bad_patch_indices)/len(img_patches_interest) * 100:.1f}%')
                cv_thres_l.append(cv_thresh)
                good_patch_num_l.append(len(good_patches_l))
                bad_patch_num_l.append(len(bad_patches_l))
                bad_patch_percentage_l.append(f'{len(bad_patch_indices)/len(img_patches_interest) * 100:.1f}%')
                # assert len(good_patches_l) + len(bad_patches_l) == len(img_patches_interest)
            except:
                print(f'Error in image')
                img_shape_l.append('None')
                wsi_level0_l.append('None')
                img_list.append('None')
                wsi_list.append('None')
                cv_thres_l.append('None')
                good_patch_num_l.append('None')
                bad_patch_num_l.append('None')
                bad_patch_percentage_l.append('None')
                continue

        else:
            img_shape_l.append('None')
            wsi_level0_l.append('None')
            img_list.append('None')
            wsi_list.append('None')
            cv_thres_l.append('None')
            good_patch_num_l.append('None')
            bad_patch_num_l.append('None')
            bad_patch_percentage_l.append('None')

    # make a dataframe of the metrics
    df_metrics = pd.DataFrame({'sample_name' : samples_read_l,
                                'adata_shape': adata_shape_l,
                                'img_shape (height, width, channels)': img_shape_l,
                                'wsi_level_0 (width, height)': wsi_level0_l,
                                'spot_pixel_size': pix_size_l,
                                'resolution (um/pixel)': resolution_l,
                                'magnification': mag_l,
                                'good_patch_num': good_patch_num_l,
                                'bad_patch_num': bad_patch_num_l,
                                'bad_patch_percentage': bad_patch_percentage_l,
                                'cv_threshold': cv_thres_l,
                                'bottom_x_y': bottom_x_y_l,
                                'bottom_x_y_pixel': bottom_x_y_pixel_l,
                                'top_x_y': top_x_y_l,
                                'top_x_y_pixel': top_x_y_pixel_l,
                                'min_pixel_x_y': min_pixel_l,
                                'max_pixel_x_y': max_pixel_l,
                                })
    
    if save_path is not None:
        df_metrics.to_csv(save_path, index=False)
        print(f'Dataframe saved to {save_path}')

    return df_metrics

def main():
    project = 'crc_14' #'lymph_23'
    if project == 'crc_14':
        st_dir = '/home/mahmoodlab8/workspace/anthony/projects/ST_codebase/data/crc_14/st_files'
        img_dir = '/home/mahmoodlab8/workspace/anthony/projects/ST_codebase/data/crc_14/images'
        df_img_name = pd.read_csv('/home/mahmoodlab8/workspace/anthony/projects/ST_codebase/data/crc_14/df_img_name_crc14.csv')
    elif project == 'lymph_23':
        st_dir = '/home/mahmoodlab8/workspace/anthony/projects/ST_codebase/data/lymph_23/st_files'
        img_dir = '/home/mahmoodlab8/workspace/anthony/projects/ST_codebase/data/lymph_23/images'
        df_img_name = pd.read_csv('/home/mahmoodlab8/workspace/anthony/projects/ST_codebase/data/lymph_23/df_img_name.csv')
    elif project == 'demo':
        st_dir = '/home/mahmoodlab8/workspace/anthony/projects/ST_codebase/data/demo_multi_cancer_data'
        img_dir = '/home/mahmoodlab8/workspace/anthony/projects/ST_codebase/data/demo_multi_cancer_data'
        df_img_name = pd.read_csv('/home/mahmoodlab8/workspace/anthony/projects/ST_codebase/data/demo_multi_cancer_data/df_img_name_demo_0225.csv')
    
    """ df_img_name is a csv that contains the sample names and the image names, os.path.join(img_dir, img_name) will give the full path to the image
            sample_name	img_name
        0	GSM5924042_frozen_a_1	a_1 Frozen.ndpi
        1	GSM5924044_frozen_a_15	a_15 Frozen.ndpi
        2	GSM5924045_frozen_a_17	a_17 Frozen.ndpi
        3	GSM5924043_frozen_a_3	a_3 Frozen.ndpi
        ...
    """

    # adata_shape_l_unprocessed = plot_unprocessed(st_dir, img_dir, df_img_name, f'./figures/{project}/unprocessed_plots')
    save_path_df_metrics = f'./figures/{project}/df_metrics_{project}.csv'
    load_img = True
    if project == 'crc_14':
        load_img = False
    df_metrics = generate_metrics_df(df_img_name, st_dir, img_dir, 
                                     load_img = load_img, 
                                     check_wsi=False,
                                     cv_thresh = 200, 
                                     save_path = save_path_df_metrics)
    
    
    # adata_list, img_list, hvgs_union, sample_names, matched_dict = load_data(st_dir, img_dir, df_img_name)
    # start, end = 0, len(adata_list)
    # # Ensure the directory exists for saving the figure
    # dist_img_save_dir = f'./figures/{project}/combined_dist_plots'
    # spatial_img_save_dir = f'./figures/{project}/spatial_plots'

    # plot_dist_plots(adata_list, sample_names, dist_img_save_dir, start, end)
    #plot_spatial_plots(adata_list, sample_names, spatial_img_save_dir, start, end)


if __name__ == "__main__":
    main()