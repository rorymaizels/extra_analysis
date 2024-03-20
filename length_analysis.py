import sys
import pandas as pd
import numpy as np
import scipy as sp      
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
import pickle

# velocity packages
import scanpy as sc
import scvelo as scv
import anndata as ann
import cellrank as cr

# and
import gget
from tqdm import tqdm
import multiprocessing

import multiprocessing
from tqdm import tqdm
import gget
import numpy as np

def process_gene(gene, sub_data):
    try:
        results = gget.search([gene], 'mus_musculus')
        gene_ids = results[results.gene_name == gene].iloc[0].ensembl_id
        mrna_seq = gget.seq(gene_ids)[1]
        length = len(mrna_seq)

        uns = sub_data[:,gene].layers['unspliced'].sum()
        tot = sub_data[:,gene].layers['spliced'].sum() + sub_data[:,gene].layers['unspliced'].sum()
        up = float(np.nan_to_num(uns/tot))

        return length, up
    except Exception as e:
        print(f"Error processing {gene}: {e}")
        return None, None

    def parallel_gene_processing(gene_list, sub_data):
    with multiprocessing.Pool(processes=int(multiprocessing.cpu_count()-4)) as pool:
        results = list((pool.starmap(process_gene, [(gene, sub_data) for gene in gene_list]), total=len(gene_list)))
    return results

mus_path = "/nemo/lab/briscoej/home/users/maizelr/transcriptomics/mouse_transcriptomics_data/full_data/mouse_full_typed_velocity.loom"
hum_path = "/nemo/lab/briscoej/home/users/maizelr/transcriptomics/human_transcriptomics_analysis/data/human_full_typed_velocity.loom"

mdata = sc.read_loom(mus_path)
hdata = sc.read_loom(hum_path)

# mouse
sub = mdata.copy()
sc.pp.filter_genes(sub, min_cells=100)
gene_list = sub.var_names
sub_data = sub.copy()
results = parallel_gene_processing(gene_list, sub_data)
lengths, unspcts = zip(*results)
np.save('muslen.npy',np.array(lengths))
np.save('musups.npy',np.array(unspcts))

# human
sub = hdata.copy()
sc.pp.filter_genes(sub, min_cells=100)
gene_list = sub.var_names
sub_data = sub.copy()
results = parallel_gene_processing(gene_list, sub_data)
lengths, unspcts = zip(*results)
np.save('humlen.npy',np.array(lengths))
np.save('humups.npy',np.array(unspcts))












