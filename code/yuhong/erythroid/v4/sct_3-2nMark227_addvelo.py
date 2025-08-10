## plot
import sctour as sct
import sys
sys.path.append('/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/veloUncertainty')
from v4_functions import *
from sctour_misc import *


def compute_sct_avg_velocity(tnode,timesteps):
    v_shape = tnode.adata.shape
    v = np.zeros(v_shape)
    for t in timesteps:
        v += compute_sctour_velocity(tnode, timestep=t)
    return v/len(timesteps)

def addvelo_sct_ery_nMark(gene_set_name, split_seed):
    import scvelo as scv
    dataset_short = 'ery'
    dataset_long = 'erythroid'
    method_prefix = 'sct'
    method = method_prefix + '_' + gene_set_name
    data_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/out/yuhong/data/v4_'+dataset_long+'/'
    fig_folder = '/home/users/y2564li/kzlinlab/projects/veloUncertainty/git/veloUncertainty/fig/yuhong/v4_'+dataset_long+'/seed'+str(split_seed)+'/'+method+'/'
    total = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total.h5ad')
    split1 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1.h5ad')
    split2 = sc.read_h5ad(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2.h5ad')
    tnode_total = sct.predict.load_model(data_folder+'seed'+str(split_seed)+'/'+method+'/tnode_'+dataset_short+'_'+method+'_total.pth')
    tnode_split1 = sct.predict.load_model(data_folder+'seed'+str(split_seed)+'/'+method+'/tnode_'+dataset_short+'_'+method+'_split1.pth')
    tnode_split2 = sct.predict.load_model(data_folder+'seed'+str(split_seed)+'/'+method+'/tnode_'+dataset_short+'_'+method+'_split2.pth')
    # recover velocity
    timesteps=[i/50 for i in range(1,11)]
    total.layers['velocity'] = compute_sct_avg_velocity(tnode_total, timesteps)
    split1.layers['velocity'] = compute_sct_avg_velocity(tnode_split1, timesteps) 
    split2.layers['velocity'] = compute_sct_avg_velocity(tnode_split2, timesteps)
    print('Velocity computed')
    # compute umap
    get_umap_sct(total)
    get_umap_sct(split1)
    get_umap_sct(split2)
    print('UMAP computed')
    # write data
    total.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_total_outputAdded.h5ad')
    split1.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split1_outputAdded.h5ad')
    split2.write(data_folder+'seed'+str(split_seed)+'/'+method+'/adata_'+dataset_short+'_'+method+'_split2_outputAdded.h5ad')
    print('###################### All done for seed'+str(split_seed)+', '+method)


for i in range(4):
    split_seed = [320, 323, 326, 329][i]
    grid_seed = 227
    gene_set_name = 'nMark' + str(grid_seed)
    addvelo_sct_ery_nMark(gene_set_name, split_seed)
