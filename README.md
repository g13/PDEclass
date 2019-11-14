# PDE solver
### install python3 (through miniconda, recommended)
miniconda: https://docs.conda.io/en/latest/miniconda.html, this is a slim version of anaconda

In windows:

run Anaconda Prompt and create a new environment(optional)

### install packages: numpy, matplotlib and their depenencies
In the Anaconda Prompt,
go to the newly created environment if you have created one
```bash
activate *new_environment* 
```
and install by:
```bash
conda install *package_name*
```
Jupyter Lab is a recommended environment to present python code, to install:
```bash
conda install -c conda-forge jupyterlab
```
to run Jupyter Lab:
```bash
jupyter lab
```
## to run examples:
### In terminal/cmd
```bash
python PDE2.py
```
### In Jupyter lab, IPython or Jupyter Notebook
```python
run PDE2
```
