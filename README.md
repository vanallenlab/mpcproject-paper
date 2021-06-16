# MPCproject

This repository contains the code and data used for the MPCproject paper.

To protect the privacy of our patient partners, potentially identifying information (e.g. addresses, individual survey answers, germline mutations, exact dates of clinical histories, etc) is provided at a summary level or omitted entirely.

## Creating the figures

For an installation-free way to generate the figures used in this paper, you can use a [Google Colab Notebook](https://colab.research.google.com/drive/1SqX5o3BogoNTp-N7MIOQPCq5rohr3T67?usp=sharing).

To generate the figures on your own computer, the best way is to use a [virtual environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). The easiest way to do this is using `conda` (install [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)).

**Instructions**

1. Clone this repository using `git`:
   ```sh
   git clone https://github.com/vanallenlab/mpcproject-paper.git
   ```
2. Create conda environment
   ```sh
   conda env create --name mpcproject_env --file=environment.yml python=3.8.5
   ```
3. Activate the conda environment:
   ```sh
   source activate mpcproject_env
   ```
4. There are then two options for running the code:

   - Simply run `python generate_figures.py`, which will run all the code and save the figures to the `figures` folder.

   - If you want a more interactive approach, you can run `mpcproject-figures.ipynb` using [Jupyter Notebook](https://jupyter.org/), which is installed by default along with `conda`. To do this, the virtual environment needs to be accessible to Jupyter as a kernel. Run the following lines of code, after which `mpcproject_kernel` will be available as a kernel to select with `Jupyter Notebook` or `Jupyter Lab`.

   ```sh
   (mpcproject_env)$ conda install ipykernel
   (mpcproject_env)$ ipython kernel install --user --name=mpcproject_kernel
   (mpcproject_env)$ conda deactivate
   ```
