import numpy as np
from torchdiffeq import odeint
import torch

def compute_sctour_velocity(tnode, timestep=1/100):
    X = tnode.adata.X
    X = np.log1p(X)

    tnode.model.eval()
    X = torch.tensor(X.toarray()).to(tnode.model.device)

    pseudtotime, posterior_mean, _ = tnode.model.encoder(X)
    differences = []

    for i in range(len(pseudtotime)):
        # Define T for the current pseudotime
        T = torch.tensor([pseudtotime[i], pseudtotime[i] + timestep])
        
        # Compute the ODE solutions
        res = odeint(tnode.model.lode_func, 
                     posterior_mean[i,:].to('cpu'), 
                     T.to('cpu'), 
                     method = tnode.model.ode_method).view(-1, tnode.model.n_latent)
        
        # Decode the ODE solutions
        vec1 = tnode.model.decoder(res[0])
        vec2 = tnode.model.decoder(res[1])
        
        # Compute the difference and append to the list
        differences.append((vec2 - vec1).detach().numpy())

    # Convert the list of differences to a numpy array
    differences_matrix = np.array(differences)
    return differences_matrix

