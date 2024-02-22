# ST_data_processing: Data Analysis of Spatial Transcriptomics 10x Visium Data

### Set up conda environment   
We need to use python 3.8 and `pip install 'scanpy[leiden]'` to install both scanpy, leidenalg (for clustering) and other dependencies.   
`pip install ipykernel --upgrade` is used for running this conda env in jupyternotebook (tested in VS Code).   
```
conda create -n stcode python=3.8
conda activate stcode  

pip install 'scanpy[leiden]'
pip install opencv-python

pip install --upgrade pip   
pip install ipykernel --upgrade 
...
```
### Functionality
1. `st_helper.py` contains all helper functions to load data and process data.   

### Data Folder Format
One must put files into this architecture and accompany a csv file, refered to be read as:   
df_img_name = pd.read_csv('/home/mahmoodlab8/workspace/anthony/projects/ST_codebase/data/lymph_23/df_img_name.csv')   

    sample_name	img_name   
0	GSM5924042_frozen_a_1	a_1 Frozen.ndpi   
1	GSM5924044_frozen_a_15	a_15 Frozen.ndpi   
2	GSM5924045_frozen_a_17	a_17 Frozen.ndpi   
3	GSM5924043_frozen_a_3	a_3 Frozen.ndpi   


Data Folder Path: ~/workspace/anthony/projects/ST_codebase/data/lymph_23:    

    .     
    ├── df_img_name.csv     
    ├── images                                 #contain identifiers that can be mapped with samples in st_files folder   
    │   ├── a_1 Frozen.ndpi   
    │   ├── a_3 Frozen.ndpi   
    │   └── ....ndpi   
    └── st_files   
        ├── GSM5924030_ffpe_c_2   
        │   ├── filtered_feature_bc_matrix.h5   
        │   ├── raw_feature_bc_matrix.h5   
        │   └── spatial   
        │       ├── aligned_fiducials.jpg   
        │       ├── detected_tissue_image.jpg   
        │       ├── scalefactors_json.json   
        │       ├── tissue_hires_image.png   
        │       ├── tissue_lowres_image.png   
        │       └── tissue_positions_list.csv   
        ├── GSM5924031_ffpe_c_3    
        │   [Assuming similar structure to GSM5924030_ffpe_c_2]    
        ...   
        ├── [Other GSM directories following the pattern of GSM5924030_ffpe_c_2]   
        │   ├── filtered_feature_bc_matrix.h5   
        │   ├── raw_feature_bc_matrix.h5    
        │   └── spatial    
        │       ├── aligned_fiducials.jpg   
        │       ├── detected_tissue_image.jpg   
        │       ├── scalefactors_json.json   
        │       ├── tissue_hires_image.png   
        │       ├── tissue_lowres_image.png     
           └── tissue_positions_list.csv    
        └── ...    
    │   
    └── ...    


### Scanpy Data Processing   
In the context of processing spatial transcriptomics data using Scanpy, the filters `total_counts`, `n_genes_by_counts`, and `log1p_n_genes_by_counts`, etc. are used to select and filter cells and genes based on their expression levels and the prevalence across the dataset. Here's a breakdown of what each of these terms means and how the filters are applied:  

**Note: Though the official description refers to cell, but in ST processing each unit is in fact a spot, which can consist of 50~100 cells,
A spot measures ~ 55 micrometers in diameter.**   
  
#### `total_counts`
- **What it is**: `total_counts` refers to the total number of transcript counts (mRNA molecules) detected in each cell (spot). It is a measure of the overall expression level within a cell.
- **How it's used**: Filtering cells based on `total_counts` helps in removing cells with very low or very high total transcript counts. Low counts might indicate dead or damaged cells, or low-quality data, while extremely high counts could be due to doublets (two cells being counted as one) or other artifacts.
    - `sc.pp.filter_cells(adata, min_counts=5000)`: This line filters out cells that have fewer than 5000 total transcript counts, removing cells with low overall expression.
    - `sc.pp.filter_cells(adata, max_counts=35000)`: Conversely, this filters out cells with more than 35000 total transcript counts, targeting cells that may be doublets or otherwise abnormal.

#### `n_genes_by_counts`
- **What it is**: This metric represents the number of genes detected in each cell that have **at least one transcript count**. It provides an idea of the gene diversity within each cell.
- **How it's used**: Filtering by `n_genes_by_counts` can help identify and exclude cells with an unusually low or high number of genes expressed, which can be indicative of poor quality or multiplets, respectively. This filtering criterion isn't explicitly shown in your provided code but would follow a similar syntax to filter cells based on the number of expressed genes.

#### `log1p_n_genes_by_counts`
- **What it is**: This is a logarithmic transformation (log1p) applied to the `n_genes_by_counts`. The `log1p` function applies a natural logarithm (log) transformation after adding 1 to each value, which helps in stabilizing the variance across cells with different expression levels.
- **How it's used**: While not directly used in the filtering steps you've mentioned, this transformation is often applied before further analysis to make the data more amenable to statistical methods that assume normal distribution or similar variance across observations.

#### `pct_counts_in_top_200_genes` and `pct_counts_in_top_500_genes`
- **What they are**: These metrics measure the percentage of total transcript counts that are found in the top 200 or top 500 most expressed genes within a cell, respectively.
- **How they're used**: High values in these metrics might indicate that a few genes are dominating the expression profile of the cell, which could be a sign of technical bias or biological processes like differentiation or stress response. It's also a useful measure to identify cells with an overly skewed expression profile, which could potentially distort downstream analyses.

#### `total_counts_mt`
- **What it is**: This metric represents the total count of mitochondrial transcripts detected in each cell.
- **How it's used**: As mentioned earlier, a high number of mitochondrial transcripts can indicate cell stress or apoptosis, as mitochondrial genes tend to be released into the cytoplasm from damaged mitochondria. Filtering cells based on `total_counts_mt` helps remove potentially compromised cells from the analysis.

#### `log1p_total_counts_mt`
- **What it is**: This is the log-transformed (using log1p, which applies `log(1+x)` to each value) version of `total_counts_mt`.
- **How it's used**: Log transformation is often applied to normalize data and reduce skewness, making the distribution of these counts more uniform and suitable for statistical analysis. It can help in visualizing and identifying outliers in the context of mitochondrial gene expression.

#### `pct_counts_mt`
- **What it is**: This metric calculates the percentage of total transcript counts that are mitochondrial.
- **How it's used**: Similar to `total_counts_mt`, but as a percentage, this measure helps in identifying cells that might be under stress or undergoing apoptosis. Cells with a high percentage of mitochondrial transcripts are often filtered out during QC.

#### `n_counts`
- **What it is**: This is another term for `total_counts`, referring to the total number of transcript counts detected in a cell.
- **How it's used**: It serves as a general measure of the cell's activity or the efficiency of RNA capture and sequencing. Cells with very low `n_counts` might be of low quality, dead, or empty droplets (in droplet-based sequencing), while very high counts might indicate doublets or multiplets.

These metrics are crucial for ensuring the reliability and interpretability of data from scRNA-seq and spatial transcriptomics experiments. They help researchers to filter out low-quality cells, adjust for technical variances, and thus, focus their analyses on healthy, single cells with high-quality gene expression profiles.
#### Other Filters
- `sc.pp.filter_genes(adata, min_cells=10)`: This line filters out genes that are detected in fewer than 10 cells. This step is crucial for removing genes that are rarely expressed across the dataset, which might be noise or specific to certain outliers, thereby focusing the analysis on genes that are more broadly relevant.

In summary, these filtering steps are essential for cleaning spatial transcriptomics data, improving the quality of downstream analyses by removing outliers and focusing on the most informative cells and genes. They help in ensuring that the data set for analysis is of high quality, which is critical for accurate interpretation and discovery in spatial transcriptomics studies.