This repository contains all the code used for the analyses in our paper, which involve five RNA velocity methods and three datasets. The `code` folder is organized by dataset and method, and includes scripts for count splitting as well as for computing replicate coherence and signal-to-random coherence scores. The `veloUncertainty` folder contains custom functions used in these analyses. The produced figures can be found in the `fig` folder.

This code is accompanying the paper “Assessing RNA velocity stability across synthetic replicates using count splitting” (https://www.biorxiv.org/content/10.1101/2024.11.23.625009).

# Instructions for reproducing results

1. Install the appropriate RNA velocity method into a virtual environment.
2. Download our GitHub using `git clone https://github.com/linnykos/veloUncertainty/` 
3. After setting up the filepaths to your environment appropriately, you can run the following lines in your script
```
import sys
sys.path.append('FILEPATH-TO-veloUncertainty-FOLDER-INSIDE-GITHUB-REPO')
from countsplit import *
```
Afterwards, you can import the relevant functions for the specific RNA velocity method. (See the list of scripts in `https://github.com/linnykos/veloUncertainty/tree/main/veloUncertainty`. The file suffixes to `v4_functions` denote the helper functions for the corresponding RNA velocity method: sct = scTour, scv = scVelo, velovi = both versions of veloVI, and utv = UniTVelo. 

The analyses details and workflow are documented for the respective datasets under the `code` folder.

# Details for the virtual environments

The following are the dependencies for the scVelo virtual environment:
```
# packages in environment at /home/users/y2564li/miniconda3/envs/scVelo:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main  
_openmp_mutex             5.1                       1_gnu  
absl-py                   2.1.0                    pypi_0    pypi
aiohttp                   3.9.3                    pypi_0    pypi
aiosignal                 1.3.1                    pypi_0    pypi
anndata                   0.10.5.post1             pypi_0    pypi
annoy                     1.17.3                   pypi_0    pypi
antlr-denter              1.3.1                    pypi_0    pypi
antlr4-python3-runtime    4.9.3                    pypi_0    pypi
array-api-compat          1.4.1                    pypi_0    pypi
async-timeout             4.0.3                    pypi_0    pypi
attrs                     23.2.0                   pypi_0    pypi
automata-lib              8.3.0                    pypi_0    pypi
bbknn                     1.6.0                    pypi_0    pypi
ca-certificates           2023.12.12           h06a4308_0  
cached-method             0.1.0                    pypi_0    pypi
cellrank                  2.0.2                    pypi_0    pypi
certifi                   2024.2.2                 pypi_0    pypi
charset-normalizer        3.3.2                    pypi_0    pypi
chex                      0.1.85                   pypi_0    pypi
click                     8.1.7                    pypi_0    pypi
contextlib2               21.6.0                   pypi_0    pypi
contourpy                 1.2.0                    pypi_0    pypi
cycler                    0.12.1                   pypi_0    pypi
cython                    3.0.9                    pypi_0    pypi
docrep                    0.3.2                    pypi_0    pypi
etils                     1.5.2                    pypi_0    pypi
exceptiongroup            1.2.0                    pypi_0    pypi
filelock                  3.13.1                   pypi_0    pypi
flax                      0.8.1                    pypi_0    pypi
fonttools                 4.49.0                   pypi_0    pypi
frozendict                2.4.1                    pypi_0    pypi
frozenlist                1.4.1                    pypi_0    pypi
fsspec                    2024.2.0                 pypi_0    pypi
get-annotations           0.1.2                    pypi_0    pypi
h5py                      3.10.0                   pypi_0    pypi
harmony                   1.2.3587                 pypi_0    pypi
harmonypy                 0.0.9                    pypi_0    pypi
idna                      3.6                      pypi_0    pypi
igraph                    0.10.8                   pypi_0    pypi
importlib-metadata        7.0.1                    pypi_0    pypi
importlib-resources       6.1.1                    pypi_0    pypi
jax                       0.4.24                   pypi_0    pypi
jaxlib                    0.4.24                   pypi_0    pypi
jinja2                    3.0.3                    pypi_0    pypi
joblib                    1.3.2                    pypi_0    pypi
kiwisolver                1.4.5                    pypi_0    pypi
ld_impl_linux-64          2.38                 h1181459_1  
libffi                    3.4.4                h6a678d5_0  
libgcc-ng                 11.2.0               h1234567_1  
libgomp                   11.2.0               h1234567_1  
libstdcxx-ng              11.2.0               h1234567_1  
lightning                 2.1.4                    pypi_0    pypi
lightning-utilities       0.10.1                   pypi_0    pypi
llvmlite                  0.42.0                   pypi_0    pypi
loompy                    3.0.7                    pypi_0    pypi
louvain                   0.8.1                    pypi_0    pypi
markdown-it-py            3.0.0                    pypi_0    pypi
markupsafe                2.1.5                    pypi_0    pypi
matplotlib                3.7.1                    pypi_0    pypi
mdurl                     0.1.2                    pypi_0    pypi
ml-collections            0.1.1                    pypi_0    pypi
ml-dtypes                 0.3.2                    pypi_0    pypi
mpmath                    1.3.0                    pypi_0    pypi
msgpack                   1.0.7                    pypi_0    pypi
mudata                    0.2.3                    pypi_0    pypi
multidict                 6.0.5                    pypi_0    pypi
multipledispatch          1.0.0                    pypi_0    pypi
natsort                   8.4.0                    pypi_0    pypi
ncurses                   6.4                  h6a678d5_0  
nest-asyncio              1.6.0                    pypi_0    pypi
networkx                  3.2.1                    pypi_0    pypi
numba                     0.59.0                   pypi_0    pypi
numpy                     1.26.4                   pypi_0    pypi
numpy-groupies            0.10.2                   pypi_0    pypi
numpyro                   0.13.2                   pypi_0    pypi
nvidia-cublas-cu12        12.1.3.1                 pypi_0    pypi
nvidia-cuda-cupti-cu12    12.1.105                 pypi_0    pypi
nvidia-cuda-nvrtc-cu12    12.1.105                 pypi_0    pypi
nvidia-cuda-runtime-cu12  12.1.105                 pypi_0    pypi
nvidia-cudnn-cu12         8.9.2.26                 pypi_0    pypi
nvidia-cufft-cu12         11.0.2.54                pypi_0    pypi
nvidia-curand-cu12        10.3.2.106               pypi_0    pypi
nvidia-cusolver-cu12      11.4.5.107               pypi_0    pypi
nvidia-cusparse-cu12      12.1.0.106               pypi_0    pypi
nvidia-nccl-cu12          2.19.3                   pypi_0    pypi
nvidia-nvjitlink-cu12     12.3.101                 pypi_0    pypi
nvidia-nvtx-cu12          12.1.105                 pypi_0    pypi
openssl                   3.0.13               h7f8727e_0  
opt-einsum                3.3.0                    pypi_0    pypi
optax                     0.1.9                    pypi_0    pypi
orbax-checkpoint          0.5.3                    pypi_0    pypi
packaging                 23.2                     pypi_0    pypi
pandas                    2.2.0                    pypi_0    pypi
patsy                     0.5.6                    pypi_0    pypi
pillow                    10.2.0                   pypi_0    pypi
pip                       23.3.1           py39h06a4308_0  
progressbar2              4.3.2                    pypi_0    pypi
protobuf                  4.25.3                   pypi_0    pypi
pydot                     2.0.0                    pypi_0    pypi
pygam                     0.9.1                    pypi_0    pypi
pygments                  2.17.2                   pypi_0    pypi
pygpcca                   1.0.4                    pypi_0    pypi
pynndescent               0.5.11                   pypi_0    pypi
pyparsing                 3.1.1                    pypi_0    pypi
pyro-api                  0.1.2                    pypi_0    pypi
pyro-ppl                  1.8.6                    pypi_0    pypi
python                    3.9.18               h955ad1f_0  
python-dateutil           2.8.2                    pypi_0    pypi
python-utils              3.8.2                    pypi_0    pypi
pytorch-lightning         2.2.0.post0              pypi_0    pypi
pytz                      2024.1                   pypi_0    pypi
pyyaml                    6.0.1                    pypi_0    pypi
readline                  8.2                  h5eee18b_0  
requests                  2.31.0                   pypi_0    pypi
rich                      13.7.0                   pypi_0    pypi
scanpy                    1.9.8                    pypi_0    pypi
scikit-learn              1.4.1.post1              pypi_0    pypi
scipy                     1.11.4                   pypi_0    pypi
scvelo                    0.3.1                    pypi_0    pypi
scvi-tools                1.1.0.post2              pypi_0    pypi
seaborn                   0.13.2                   pypi_0    pypi
session-info              1.0.0                    pypi_0    pypi
setuptools                68.2.2           py39h06a4308_0  
six                       1.16.0                   pypi_0    pypi
skeleton-methods          0.0.4                    pypi_0    pypi
sqlite                    3.41.2               h5eee18b_0  
statsmodels               0.14.1                   pypi_0    pypi
stdlib-list               0.10.0                   pypi_0    pypi
sympy                     1.12                     pypi_0    pypi
tensorstore               0.1.53                   pypi_0    pypi
texttable                 1.7.0                    pypi_0    pypi
threadpoolctl             3.3.0                    pypi_0    pypi
tk                        8.6.12               h1ccaba5_0  
toolz                     0.12.1                   pypi_0    pypi
torch                     2.2.0                    pypi_0    pypi
torchmetrics              1.3.1                    pypi_0    pypi
tqdm                      4.66.2                   pypi_0    pypi
triton                    2.2.0                    pypi_0    pypi
typing-extensions         4.9.0                    pypi_0    pypi
tzdata                    2024.1                   pypi_0    pypi
umap-learn                0.5.5                    pypi_0    pypi
urllib3                   2.2.1                    pypi_0    pypi
wheel                     0.41.2           py39h06a4308_0  
wrapt                     1.16.0                   pypi_0    pypi
xz                        5.4.5                h5eee18b_0  
yarl                      1.9.4                    pypi_0    pypi
zipp                      3.17.0                   pypi_0    pypi
zlib                      1.2.13               h5eee18b_0  
```

The following are the dependencies for the UniTVelo virtual environment:
```
# packages in environment at /home/users/y2564li/miniconda3/envs/utvClone:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main  
_openmp_mutex             5.1                       1_gnu  
absl-py                   2.1.0                    pypi_0    pypi
aiohttp                   3.9.3                    pypi_0    pypi
aiosignal                 1.3.1                    pypi_0    pypi
anndata                   0.10.5.post1             pypi_0    pypi
annoy                     1.17.3                   pypi_0    pypi
antlr-denter              1.3.1                    pypi_0    pypi
antlr4-python3-runtime    4.9.3                    pypi_0    pypi
anyio                     4.3.0                    pypi_0    pypi
argon2-cffi               23.1.0                   pypi_0    pypi
argon2-cffi-bindings      21.2.0                   pypi_0    pypi
array-api-compat          1.4.1                    pypi_0    pypi
arrow                     1.3.0                    pypi_0    pypi
asttokens                 2.4.1                    pypi_0    pypi
astunparse                1.6.3                    pypi_0    pypi
async-lru                 2.0.4                    pypi_0    pypi
async-timeout             4.0.3                    pypi_0    pypi
attrs                     23.2.0                   pypi_0    pypi
automata-lib              8.3.0                    pypi_0    pypi
babel                     2.14.0                   pypi_0    pypi
bbknn                     1.6.0                    pypi_0    pypi
beautifulsoup4            4.12.3                   pypi_0    pypi
bleach                    6.1.0                    pypi_0    pypi
ca-certificates           2023.12.12           h06a4308_0  
cached-method             0.1.0                    pypi_0    pypi
cellrank                  2.0.2                    pypi_0    pypi
certifi                   2024.2.2                 pypi_0    pypi
cffi                      1.16.0                   pypi_0    pypi
charset-normalizer        3.3.2                    pypi_0    pypi
chex                      0.1.85                   pypi_0    pypi
click                     8.1.7                    pypi_0    pypi
comm                      0.2.2                    pypi_0    pypi
contextlib2               21.6.0                   pypi_0    pypi
contourpy                 1.2.0                    pypi_0    pypi
cycler                    0.12.1                   pypi_0    pypi
cython                    3.0.9                    pypi_0    pypi
debugpy                   1.8.1                    pypi_0    pypi
decorator                 5.1.1                    pypi_0    pypi
defusedxml                0.7.1                    pypi_0    pypi
docrep                    0.3.2                    pypi_0    pypi
etils                     1.5.2                    pypi_0    pypi
exceptiongroup            1.2.0                    pypi_0    pypi
executing                 2.0.1                    pypi_0    pypi
fastjsonschema            2.19.1                   pypi_0    pypi
filelock                  3.13.1                   pypi_0    pypi
flatbuffers               24.3.25                  pypi_0    pypi
flax                      0.8.1                    pypi_0    pypi
fonttools                 4.49.0                   pypi_0    pypi
fqdn                      1.5.1                    pypi_0    pypi
frozendict                2.4.1                    pypi_0    pypi
frozenlist                1.4.1                    pypi_0    pypi
fsspec                    2024.2.0                 pypi_0    pypi
gast                      0.5.4                    pypi_0    pypi
get-annotations           0.1.2                    pypi_0    pypi
google-pasta              0.2.0                    pypi_0    pypi
grpcio                    1.62.1                   pypi_0    pypi
h11                       0.14.0                   pypi_0    pypi
h5py                      3.10.0                   pypi_0    pypi
harmony                   1.2.3587                 pypi_0    pypi
harmonypy                 0.0.9                    pypi_0    pypi
httpcore                  1.0.5                    pypi_0    pypi
httpx                     0.27.0                   pypi_0    pypi
idna                      3.6                      pypi_0    pypi
igraph                    0.10.8                   pypi_0    pypi
importlib-metadata        7.0.1                    pypi_0    pypi
importlib-resources       6.1.1                    pypi_0    pypi
iprogress                 0.4                      pypi_0    pypi
ipykernel                 6.29.4                   pypi_0    pypi
ipython                   8.18.1                   pypi_0    pypi
ipywidgets                8.1.2                    pypi_0    pypi
isoduration               20.11.0                  pypi_0    pypi
jax                       0.4.24                   pypi_0    pypi
jaxlib                    0.4.24                   pypi_0    pypi
jedi                      0.19.1                   pypi_0    pypi
jinja2                    3.0.3                    pypi_0    pypi
joblib                    1.3.2                    pypi_0    pypi
json5                     0.9.25                   pypi_0    pypi
jsonpointer               2.4                      pypi_0    pypi
jsonschema                4.21.1                   pypi_0    pypi
jsonschema-specifications 2023.12.1                pypi_0    pypi
jupyter                   1.0.0                    pypi_0    pypi
jupyter-client            8.6.1                    pypi_0    pypi
jupyter-console           6.6.3                    pypi_0    pypi
jupyter-core              5.7.2                    pypi_0    pypi
jupyter-events            0.10.0                   pypi_0    pypi
jupyter-lsp               2.2.5                    pypi_0    pypi
jupyter-server            2.14.0                   pypi_0    pypi
jupyter-server-terminals  0.5.3                    pypi_0    pypi
jupyterlab                4.1.6                    pypi_0    pypi
jupyterlab-pygments       0.3.0                    pypi_0    pypi
jupyterlab-server         2.26.0                   pypi_0    pypi
jupyterlab-widgets        3.0.10                   pypi_0    pypi
keras                     3.2.1                    pypi_0    pypi
kiwisolver                1.4.5                    pypi_0    pypi
ld_impl_linux-64          2.38                 h1181459_1  
libclang                  18.1.1                   pypi_0    pypi
libffi                    3.4.4                h6a678d5_0  
libgcc-ng                 11.2.0               h1234567_1  
libgomp                   11.2.0               h1234567_1  
libstdcxx-ng              11.2.0               h1234567_1  
lightning                 2.1.4                    pypi_0    pypi
lightning-utilities       0.10.1                   pypi_0    pypi
llvmlite                  0.42.0                   pypi_0    pypi
loompy                    3.0.7                    pypi_0    pypi
louvain                   0.8.1                    pypi_0    pypi
markdown                  3.6                      pypi_0    pypi
markdown-it-py            3.0.0                    pypi_0    pypi
markupsafe                2.1.5                    pypi_0    pypi
matplotlib                3.7.1                    pypi_0    pypi
matplotlib-inline         0.1.6                    pypi_0    pypi
mdurl                     0.1.2                    pypi_0    pypi
mistune                   3.0.2                    pypi_0    pypi
ml-collections            0.1.1                    pypi_0    pypi
ml-dtypes                 0.3.2                    pypi_0    pypi
mpmath                    1.3.0                    pypi_0    pypi
msgpack                   1.0.7                    pypi_0    pypi
mudata                    0.2.3                    pypi_0    pypi
multidict                 6.0.5                    pypi_0    pypi
multipledispatch          1.0.0                    pypi_0    pypi
namex                     0.0.7                    pypi_0    pypi
natsort                   8.4.0                    pypi_0    pypi
nbclient                  0.10.0                   pypi_0    pypi
nbconvert                 7.16.3                   pypi_0    pypi
nbformat                  5.10.4                   pypi_0    pypi
ncurses                   6.4                  h6a678d5_0  
nest-asyncio              1.6.0                    pypi_0    pypi
networkx                  3.2.1                    pypi_0    pypi
notebook                  7.1.2                    pypi_0    pypi
notebook-shim             0.2.4                    pypi_0    pypi
numba                     0.59.0                   pypi_0    pypi
numpy                     1.26.4                   pypi_0    pypi
numpy-groupies            0.10.2                   pypi_0    pypi
numpyro                   0.13.2                   pypi_0    pypi
nvidia-cublas-cu12        12.1.3.1                 pypi_0    pypi
nvidia-cuda-cupti-cu12    12.1.105                 pypi_0    pypi
nvidia-cuda-nvrtc-cu12    12.1.105                 pypi_0    pypi
nvidia-cuda-runtime-cu12  12.1.105                 pypi_0    pypi
nvidia-cudnn-cu12         8.9.2.26                 pypi_0    pypi
nvidia-cufft-cu12         11.0.2.54                pypi_0    pypi
nvidia-curand-cu12        10.3.2.106               pypi_0    pypi
nvidia-cusolver-cu12      11.4.5.107               pypi_0    pypi
nvidia-cusparse-cu12      12.1.0.106               pypi_0    pypi
nvidia-nccl-cu12          2.19.3                   pypi_0    pypi
nvidia-nvjitlink-cu12     12.3.101                 pypi_0    pypi
nvidia-nvtx-cu12          12.1.105                 pypi_0    pypi
openssl                   3.0.13               h7f8727e_0  
opt-einsum                3.3.0                    pypi_0    pypi
optax                     0.1.9                    pypi_0    pypi
optree                    0.11.0                   pypi_0    pypi
orbax-checkpoint          0.5.3                    pypi_0    pypi
overrides                 7.7.0                    pypi_0    pypi
packaging                 23.2                     pypi_0    pypi
pandas                    2.2.0                    pypi_0    pypi
pandocfilters             1.5.1                    pypi_0    pypi
parso                     0.8.4                    pypi_0    pypi
patsy                     0.5.6                    pypi_0    pypi
pexpect                   4.9.0                    pypi_0    pypi
pillow                    10.2.0                   pypi_0    pypi
pip                       23.3.1           py39h06a4308_0  
platformdirs              4.2.0                    pypi_0    pypi
progressbar2              4.3.2                    pypi_0    pypi
prometheus-client         0.20.0                   pypi_0    pypi
prompt-toolkit            3.0.43                   pypi_0    pypi
protobuf                  4.25.3                   pypi_0    pypi
psutil                    5.9.8                    pypi_0    pypi
ptyprocess                0.7.0                    pypi_0    pypi
pure-eval                 0.2.2                    pypi_0    pypi
pycparser                 2.22                     pypi_0    pypi
pydot                     2.0.0                    pypi_0    pypi
pygam                     0.9.1                    pypi_0    pypi
pygments                  2.17.2                   pypi_0    pypi
pygpcca                   1.0.4                    pypi_0    pypi
pynndescent               0.5.11                   pypi_0    pypi
pyparsing                 3.1.1                    pypi_0    pypi
pyro-api                  0.1.2                    pypi_0    pypi
pyro-ppl                  1.8.6                    pypi_0    pypi
python                    3.9.18               h955ad1f_0  
python-dateutil           2.8.2                    pypi_0    pypi
python-json-logger        2.0.7                    pypi_0    pypi
python-utils              3.8.2                    pypi_0    pypi
pytorch-lightning         2.2.0.post0              pypi_0    pypi
pytz                      2024.1                   pypi_0    pypi
pyyaml                    6.0.1                    pypi_0    pypi
pyzmq                     25.1.2                   pypi_0    pypi
qtconsole                 5.5.1                    pypi_0    pypi
qtpy                      2.4.1                    pypi_0    pypi
readline                  8.2                  h5eee18b_0  
referencing               0.34.0                   pypi_0    pypi
requests                  2.31.0                   pypi_0    pypi
rfc3339-validator         0.1.4                    pypi_0    pypi
rfc3986-validator         0.1.1                    pypi_0    pypi
rich                      13.7.0                   pypi_0    pypi
rpds-py                   0.18.0                   pypi_0    pypi
scanpy                    1.9.8                    pypi_0    pypi
scikit-learn              1.1.3                    pypi_0    pypi
scipy                     1.11.4                   pypi_0    pypi
scvelo                    0.3.1                    pypi_0    pypi
scvi-tools                1.1.0.post2              pypi_0    pypi
seaborn                   0.13.2                   pypi_0    pypi
send2trash                1.8.3                    pypi_0    pypi
session-info              1.0.0                    pypi_0    pypi
setuptools                68.2.2           py39h06a4308_0  
six                       1.16.0                   pypi_0    pypi
skeleton-methods          0.0.4                    pypi_0    pypi
sniffio                   1.3.1                    pypi_0    pypi
soupsieve                 2.5                      pypi_0    pypi
sqlite                    3.41.2               h5eee18b_0  
stack-data                0.6.3                    pypi_0    pypi
statsmodels               0.14.1                   pypi_0    pypi
stdlib-list               0.10.0                   pypi_0    pypi
sympy                     1.12                     pypi_0    pypi
tensorboard               2.16.2                   pypi_0    pypi
tensorboard-data-server   0.7.2                    pypi_0    pypi
tensorflow                2.16.1                   pypi_0    pypi
tensorflow-io-gcs-filesystem 0.36.0                   pypi_0    pypi
tensorstore               0.1.53                   pypi_0    pypi
termcolor                 2.4.0                    pypi_0    pypi
terminado                 0.18.1                   pypi_0    pypi
texttable                 1.7.0                    pypi_0    pypi
tf-keras                  2.16.0                   pypi_0    pypi
threadpoolctl             3.3.0                    pypi_0    pypi
tinycss2                  1.2.1                    pypi_0    pypi
tk                        8.6.12               h1ccaba5_0  
tomli                     2.0.1                    pypi_0    pypi
toolz                     0.12.1                   pypi_0    pypi
torch                     2.2.0                    pypi_0    pypi
torchmetrics              1.3.1                    pypi_0    pypi
tornado                   6.4                      pypi_0    pypi
tqdm                      4.66.2                   pypi_0    pypi
traitlets                 5.14.2                   pypi_0    pypi
triton                    2.2.0                    pypi_0    pypi
types-python-dateutil     2.9.0.20240316           pypi_0    pypi
typing-extensions         4.9.0                    pypi_0    pypi
tzdata                    2024.1                   pypi_0    pypi
umap-learn                0.5.5                    pypi_0    pypi
unitvelo                  0.2.5.2                  pypi_0    pypi
uri-template              1.3.0                    pypi_0    pypi
urllib3                   2.2.1                    pypi_0    pypi
wcwidth                   0.2.13                   pypi_0    pypi
webcolors                 1.13                     pypi_0    pypi
webencodings              0.5.1                    pypi_0    pypi
websocket-client          1.7.0                    pypi_0    pypi
werkzeug                  3.0.2                    pypi_0    pypi
wheel                     0.41.2           py39h06a4308_0  
widgetsnbextension        4.0.10                   pypi_0    pypi
wrapt                     1.16.0                   pypi_0    pypi
xz                        5.4.5                h5eee18b_0  
yarl                      1.9.4                    pypi_0    pypi
zipp                      3.17.0                   pypi_0    pypi
zlib                      1.2.13               h5eee18b_0  
```

The following are the dependencies for the scTour virtual environment:
```
# packages in environment at /home/users/y2564li/miniconda3/envs/sctClone:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main  
_openmp_mutex             5.1                       1_gnu  
absl-py                   2.1.0                    pypi_0    pypi
aiohttp                   3.9.3                    pypi_0    pypi
aiosignal                 1.3.1                    pypi_0    pypi
anndata                   0.10.5.post1             pypi_0    pypi
annoy                     1.17.3                   pypi_0    pypi
antlr-denter              1.3.1                    pypi_0    pypi
antlr4-python3-runtime    4.9.3                    pypi_0    pypi
array-api-compat          1.4.1                    pypi_0    pypi
async-timeout             4.0.3                    pypi_0    pypi
attrs                     23.2.0                   pypi_0    pypi
automata-lib              8.3.0                    pypi_0    pypi
bbknn                     1.6.0                    pypi_0    pypi
ca-certificates           2023.12.12           h06a4308_0  
cached-method             0.1.0                    pypi_0    pypi
cellrank                  2.0.2                    pypi_0    pypi
certifi                   2024.2.2                 pypi_0    pypi
charset-normalizer        3.3.2                    pypi_0    pypi
chex                      0.1.85                   pypi_0    pypi
click                     8.1.7                    pypi_0    pypi
contextlib2               21.6.0                   pypi_0    pypi
contourpy                 1.2.0                    pypi_0    pypi
cycler                    0.12.1                   pypi_0    pypi
cython                    3.0.9                    pypi_0    pypi
docrep                    0.3.2                    pypi_0    pypi
etils                     1.5.2                    pypi_0    pypi
exceptiongroup            1.2.0                    pypi_0    pypi
filelock                  3.13.1                   pypi_0    pypi
flax                      0.8.1                    pypi_0    pypi
fonttools                 4.49.0                   pypi_0    pypi
frozendict                2.4.1                    pypi_0    pypi
frozenlist                1.4.1                    pypi_0    pypi
fsspec                    2024.2.0                 pypi_0    pypi
get-annotations           0.1.2                    pypi_0    pypi
h5py                      3.10.0                   pypi_0    pypi
harmony                   1.2.3587                 pypi_0    pypi
harmonypy                 0.0.9                    pypi_0    pypi
idna                      3.6                      pypi_0    pypi
igraph                    0.10.8                   pypi_0    pypi
importlib-metadata        7.0.1                    pypi_0    pypi
importlib-resources       6.1.1                    pypi_0    pypi
jax                       0.4.24                   pypi_0    pypi
jaxlib                    0.4.24                   pypi_0    pypi
jinja2                    3.0.3                    pypi_0    pypi
joblib                    1.3.2                    pypi_0    pypi
kiwisolver                1.4.5                    pypi_0    pypi
ld_impl_linux-64          2.38                 h1181459_1  
leidenalg                 0.10.2                   pypi_0    pypi
libffi                    3.4.4                h6a678d5_0  
libgcc-ng                 11.2.0               h1234567_1  
libgomp                   11.2.0               h1234567_1  
libstdcxx-ng              11.2.0               h1234567_1  
lightning                 2.1.4                    pypi_0    pypi
lightning-utilities       0.10.1                   pypi_0    pypi
llvmlite                  0.42.0                   pypi_0    pypi
loompy                    3.0.7                    pypi_0    pypi
louvain                   0.8.1                    pypi_0    pypi
markdown-it-py            3.0.0                    pypi_0    pypi
markupsafe                2.1.5                    pypi_0    pypi
matplotlib                3.7.1                    pypi_0    pypi
mdurl                     0.1.2                    pypi_0    pypi
ml-collections            0.1.1                    pypi_0    pypi
ml-dtypes                 0.3.2                    pypi_0    pypi
mpmath                    1.3.0                    pypi_0    pypi
msgpack                   1.0.7                    pypi_0    pypi
mudata                    0.2.3                    pypi_0    pypi
multidict                 6.0.5                    pypi_0    pypi
multipledispatch          1.0.0                    pypi_0    pypi
natsort                   8.4.0                    pypi_0    pypi
ncurses                   6.4                  h6a678d5_0  
nest-asyncio              1.6.0                    pypi_0    pypi
networkx                  3.2.1                    pypi_0    pypi
numba                     0.59.0                   pypi_0    pypi
numpy                     1.26.4                   pypi_0    pypi
numpy-groupies            0.10.2                   pypi_0    pypi
numpyro                   0.13.2                   pypi_0    pypi
nvidia-cublas-cu12        12.1.3.1                 pypi_0    pypi
nvidia-cuda-cupti-cu12    12.1.105                 pypi_0    pypi
nvidia-cuda-nvrtc-cu12    12.1.105                 pypi_0    pypi
nvidia-cuda-runtime-cu12  12.1.105                 pypi_0    pypi
nvidia-cudnn-cu12         8.9.2.26                 pypi_0    pypi
nvidia-cufft-cu12         11.0.2.54                pypi_0    pypi
nvidia-curand-cu12        10.3.2.106               pypi_0    pypi
nvidia-cusolver-cu12      11.4.5.107               pypi_0    pypi
nvidia-cusparse-cu12      12.1.0.106               pypi_0    pypi
nvidia-nccl-cu12          2.19.3                   pypi_0    pypi
nvidia-nvjitlink-cu12     12.3.101                 pypi_0    pypi
nvidia-nvtx-cu12          12.1.105                 pypi_0    pypi
openssl                   3.0.13               h7f8727e_0  
opt-einsum                3.3.0                    pypi_0    pypi
optax                     0.1.9                    pypi_0    pypi
orbax-checkpoint          0.5.3                    pypi_0    pypi
packaging                 23.2                     pypi_0    pypi
pandas                    2.2.0                    pypi_0    pypi
patsy                     0.5.6                    pypi_0    pypi
pillow                    10.2.0                   pypi_0    pypi
pip                       23.3.1           py39h06a4308_0  
progressbar2              4.3.2                    pypi_0    pypi
protobuf                  4.25.3                   pypi_0    pypi
pydot                     2.0.0                    pypi_0    pypi
pygam                     0.9.1                    pypi_0    pypi
pygments                  2.17.2                   pypi_0    pypi
pygpcca                   1.0.4                    pypi_0    pypi
pynndescent               0.5.11                   pypi_0    pypi
pyparsing                 3.1.1                    pypi_0    pypi
pyro-api                  0.1.2                    pypi_0    pypi
pyro-ppl                  1.8.6                    pypi_0    pypi
python                    3.9.18               h955ad1f_0  
python-dateutil           2.8.2                    pypi_0    pypi
python-utils              3.8.2                    pypi_0    pypi
pytorch-lightning         2.2.0.post0              pypi_0    pypi
pytz                      2024.1                   pypi_0    pypi
pyyaml                    6.0.1                    pypi_0    pypi
readline                  8.2                  h5eee18b_0  
requests                  2.31.0                   pypi_0    pypi
rich                      13.7.0                   pypi_0    pypi
scanpy                    1.9.8                    pypi_0    pypi
scikit-learn              1.4.1.post1              pypi_0    pypi
scipy                     1.11.4                   pypi_0    pypi
sctour                    1.0.0                    pypi_0    pypi
scvelo                    0.3.1                    pypi_0    pypi
scvi-tools                1.1.0.post2              pypi_0    pypi
seaborn                   0.13.2                   pypi_0    pypi
session-info              1.0.0                    pypi_0    pypi
setuptools                68.2.2           py39h06a4308_0  
six                       1.16.0                   pypi_0    pypi
skeleton-methods          0.0.4                    pypi_0    pypi
sqlite                    3.41.2               h5eee18b_0  
statsmodels               0.14.1                   pypi_0    pypi
stdlib-list               0.10.0                   pypi_0    pypi
sympy                     1.12                     pypi_0    pypi
tensorstore               0.1.53                   pypi_0    pypi
texttable                 1.7.0                    pypi_0    pypi
threadpoolctl             3.3.0                    pypi_0    pypi
tk                        8.6.12               h1ccaba5_0  
toolz                     0.12.1                   pypi_0    pypi
torch                     2.2.0                    pypi_0    pypi
torchdiffeq               0.2.4                    pypi_0    pypi
torchmetrics              1.3.1                    pypi_0    pypi
tqdm                      4.66.2                   pypi_0    pypi
triton                    2.2.0                    pypi_0    pypi
typing-extensions         4.9.0                    pypi_0    pypi
tzdata                    2024.1                   pypi_0    pypi
umap-learn                0.5.5                    pypi_0    pypi
urllib3                   2.2.1                    pypi_0    pypi
wheel                     0.41.2           py39h06a4308_0  
wrapt                     1.16.0                   pypi_0    pypi
xz                        5.4.5                h5eee18b_0  
yarl                      1.9.4                    pypi_0    pypi
zipp                      3.17.0                   pypi_0    pypi
zlib                      1.2.13               h5eee18b_0
```

The following are the dependencies for the VeloVI virtual environment:
```
# packages in environment at /home/users/y2564li/miniconda3/envs/veloviClone:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                        main  
_openmp_mutex             5.1                       1_gnu  
absl-py                   2.1.0                    pypi_0    pypi
aiohttp                   3.9.3                    pypi_0    pypi
aiosignal                 1.3.1                    pypi_0    pypi
alabaster                 0.7.16                   pypi_0    pypi
anndata                   0.10.5.post1             pypi_0    pypi
annoy                     1.17.3                   pypi_0    pypi
antlr-denter              1.3.1                    pypi_0    pypi
antlr4-python3-runtime    4.9.3                    pypi_0    pypi
array-api-compat          1.4.1                    pypi_0    pypi
async-timeout             4.0.3                    pypi_0    pypi
attrs                     23.2.0                   pypi_0    pypi
automata-lib              8.3.0                    pypi_0    pypi
babel                     2.15.0                   pypi_0    pypi
bbknn                     1.6.0                    pypi_0    pypi
ca-certificates           2023.12.12           h06a4308_0  
cached-method             0.1.0                    pypi_0    pypi
cellrank                  2.0.2                    pypi_0    pypi
certifi                   2024.2.2                 pypi_0    pypi
charset-normalizer        3.3.2                    pypi_0    pypi
chex                      0.1.85                   pypi_0    pypi
click                     8.1.7                    pypi_0    pypi
contextlib2               21.6.0                   pypi_0    pypi
contourpy                 1.2.0                    pypi_0    pypi
cycler                    0.12.1                   pypi_0    pypi
cython                    3.0.9                    pypi_0    pypi
docrep                    0.3.2                    pypi_0    pypi
docutils                  0.21.2                   pypi_0    pypi
etils                     1.5.2                    pypi_0    pypi
exceptiongroup            1.2.0                    pypi_0    pypi
filelock                  3.13.1                   pypi_0    pypi
flax                      0.8.1                    pypi_0    pypi
fonttools                 4.49.0                   pypi_0    pypi
frozendict                2.4.1                    pypi_0    pypi
frozenlist                1.4.1                    pypi_0    pypi
fsspec                    2024.2.0                 pypi_0    pypi
get-annotations           0.1.2                    pypi_0    pypi
h5py                      3.10.0                   pypi_0    pypi
harmony                   1.2.3587                 pypi_0    pypi
harmonypy                 0.0.9                    pypi_0    pypi
idna                      3.6                      pypi_0    pypi
igraph                    0.11.5                   pypi_0    pypi
imagesize                 1.4.1                    pypi_0    pypi
importlib-metadata        7.0.1                    pypi_0    pypi
importlib-resources       6.1.1                    pypi_0    pypi
jax                       0.4.24                   pypi_0    pypi
jaxlib                    0.4.24                   pypi_0    pypi
jinja2                    3.0.3                    pypi_0    pypi
joblib                    1.3.2                    pypi_0    pypi
kiwisolver                1.4.5                    pypi_0    pypi
ld_impl_linux-64          2.38                 h1181459_1  
libffi                    3.4.4                h6a678d5_0  
libgcc-ng                 11.2.0               h1234567_1  
libgomp                   11.2.0               h1234567_1  
libstdcxx-ng              11.2.0               h1234567_1  
lightning                 2.1.4                    pypi_0    pypi
lightning-utilities       0.10.1                   pypi_0    pypi
llvmlite                  0.42.0                   pypi_0    pypi
loompy                    3.0.7                    pypi_0    pypi
louvain                   0.8.1                    pypi_0    pypi
markdown-it-py            3.0.0                    pypi_0    pypi
markupsafe                2.1.5                    pypi_0    pypi
matplotlib                3.7.1                    pypi_0    pypi
mdurl                     0.1.2                    pypi_0    pypi
ml-collections            0.1.1                    pypi_0    pypi
ml-dtypes                 0.3.2                    pypi_0    pypi
mpmath                    1.3.0                    pypi_0    pypi
msgpack                   1.0.7                    pypi_0    pypi
mudata                    0.2.3                    pypi_0    pypi
multidict                 6.0.5                    pypi_0    pypi
multipledispatch          1.0.0                    pypi_0    pypi
natsort                   8.4.0                    pypi_0    pypi
ncurses                   6.4                  h6a678d5_0  
nest-asyncio              1.6.0                    pypi_0    pypi
networkx                  3.2.1                    pypi_0    pypi
numba                     0.59.0                   pypi_0    pypi
numpy                     1.26.4                   pypi_0    pypi
numpy-groupies            0.10.2                   pypi_0    pypi
numpyro                   0.13.2                   pypi_0    pypi
nvidia-cublas-cu12        12.1.3.1                 pypi_0    pypi
nvidia-cuda-cupti-cu12    12.1.105                 pypi_0    pypi
nvidia-cuda-nvrtc-cu12    12.1.105                 pypi_0    pypi
nvidia-cuda-runtime-cu12  12.1.105                 pypi_0    pypi
nvidia-cudnn-cu12         8.9.2.26                 pypi_0    pypi
nvidia-cufft-cu12         11.0.2.54                pypi_0    pypi
nvidia-curand-cu12        10.3.2.106               pypi_0    pypi
nvidia-cusolver-cu12      11.4.5.107               pypi_0    pypi
nvidia-cusparse-cu12      12.1.0.106               pypi_0    pypi
nvidia-nccl-cu12          2.19.3                   pypi_0    pypi
nvidia-nvjitlink-cu12     12.3.101                 pypi_0    pypi
nvidia-nvtx-cu12          12.1.105                 pypi_0    pypi
openssl                   3.0.13               h7f8727e_0  
opt-einsum                3.3.0                    pypi_0    pypi
optax                     0.1.9                    pypi_0    pypi
orbax-checkpoint          0.5.3                    pypi_0    pypi
packaging                 23.2                     pypi_0    pypi
pandas                    2.2.0                    pypi_0    pypi
patsy                     0.5.6                    pypi_0    pypi
pillow                    10.2.0                   pypi_0    pypi
pip                       23.3.1           py39h06a4308_0  
progressbar2              4.3.2                    pypi_0    pypi
protobuf                  4.25.3                   pypi_0    pypi
pydot                     2.0.0                    pypi_0    pypi
pygam                     0.9.1                    pypi_0    pypi
pygments                  2.17.2                   pypi_0    pypi
pygpcca                   1.0.4                    pypi_0    pypi
pynndescent               0.5.11                   pypi_0    pypi
pyparsing                 3.1.1                    pypi_0    pypi
pyro-api                  0.1.2                    pypi_0    pypi
pyro-ppl                  1.8.6                    pypi_0    pypi
python                    3.9.18               h955ad1f_0  
python-dateutil           2.8.2                    pypi_0    pypi
python-igraph             0.11.5                   pypi_0    pypi
python-utils              3.8.2                    pypi_0    pypi
pytorch-lightning         2.2.0.post0              pypi_0    pypi
pytz                      2024.1                   pypi_0    pypi
pyyaml                    6.0.1                    pypi_0    pypi
readline                  8.2                  h5eee18b_0  
requests                  2.31.0                   pypi_0    pypi
rich                      13.7.0                   pypi_0    pypi
scanpy                    1.9.8                    pypi_0    pypi
scikit-learn              1.1.3                    pypi_0    pypi
scipy                     1.11.4                   pypi_0    pypi
scvelo                    0.3.1                    pypi_0    pypi
scvi-tools                1.1.0.post2              pypi_0    pypi
seaborn                   0.13.2                   pypi_0    pypi
session-info              1.0.0                    pypi_0    pypi
setuptools                68.2.2           py39h06a4308_0  
six                       1.16.0                   pypi_0    pypi
skeleton-methods          0.0.4                    pypi_0    pypi
snowballstemmer           2.2.0                    pypi_0    pypi
sphinx                    7.3.7                    pypi_0    pypi
sphinx-autodoc-typehints  2.2.2                    pypi_0    pypi
sphinxcontrib-applehelp   1.0.8                    pypi_0    pypi
sphinxcontrib-devhelp     1.0.6                    pypi_0    pypi
sphinxcontrib-htmlhelp    2.0.5                    pypi_0    pypi
sphinxcontrib-jsmath      1.0.1                    pypi_0    pypi
sphinxcontrib-qthelp      1.0.7                    pypi_0    pypi
sphinxcontrib-serializinghtml 1.1.10                   pypi_0    pypi
sqlite                    3.41.2               h5eee18b_0  
statsmodels               0.14.1                   pypi_0    pypi
stdlib-list               0.10.0                   pypi_0    pypi
sympy                     1.12                     pypi_0    pypi
tensorstore               0.1.53                   pypi_0    pypi
texttable                 1.7.0                    pypi_0    pypi
threadpoolctl             3.3.0                    pypi_0    pypi
tk                        8.6.12               h1ccaba5_0  
tomli                     2.0.1                    pypi_0    pypi
toolz                     0.12.1                   pypi_0    pypi
torch                     2.2.0                    pypi_0    pypi
torchmetrics              1.3.1                    pypi_0    pypi
tqdm                      4.66.2                   pypi_0    pypi
triton                    2.2.0                    pypi_0    pypi
typing-extensions         4.9.0                    pypi_0    pypi
tzdata                    2024.1                   pypi_0    pypi
umap-learn                0.5.5                    pypi_0    pypi
urllib3                   2.2.1                    pypi_0    pypi
velovi                    0.3.1                    pypi_0    pypi
wheel                     0.41.2           py39h06a4308_0  
wrapt                     1.16.0                   pypi_0    pypi
xz                        5.4.5                h5eee18b_0  
yarl                      1.9.4                    pypi_0    pypi
zipp                      3.17.0                   pypi_0    pypi
zlib                      1.2.13               h5eee18b_0
``` 
