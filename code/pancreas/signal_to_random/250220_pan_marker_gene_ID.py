import numpy as np
import pandas as pd

import requests
def get_ensembl_id(gene_name):
    url = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{gene_name}?content-type=application/json"
    response = requests.get(url)
    if response.ok:
        data = response.json()
        return data[0]["id"] if data else None
    return None

data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/'

genes_mark = pd.read_csv(data_folder+'v4_pancreas/pan_marker_gene_name.csv')
genes_mark = genes_mark['gene_name']
ensembl_ids = {gene: get_ensembl_id(gene) for gene in genes_mark.values.flatten()}.values()
pd.DataFrame({'gene_id': ensembl_ids, 'gene_name': genes_mark.values.flatten()}).to_csv(data_folder+'v4_pancreas/pan_marker_genes.csv')

