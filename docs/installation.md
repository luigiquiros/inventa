<h1>Installation</h1>
 
There are two ways to work with inventa: 

### A) Running inventa with Binder:

Binder allows to run inventa on the cloud with a Binder instance, which is really convenient but you need to save the parameters and results locally as these instances are shutting down after 15 min of inactivity.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/luigiquiros/inventa/main?urlpath=lab/tree/notebook/inventa.ipynb)

(this may take several minutes to build Jupyter assets ... please wait).

#### how to use Binder?

- All the necessary tables described below can be 'drag & drop' in the folder `data/` direclty in the browser, you'll find the folder in the left side.
- The output can be found in the folder `results/` as a TAB separated file.
- If needed, you can modify the jupyter notebook directly in the browser (make sure to save it locally).
- As explained below, if needed, you can modify the `inventa.py` parameters inside `src/inventa.py`, this can be done as well directly in the browser.

### B) Running inventa locally:

### Install the conda environment

First make sure to have [anaconda](https://www.anaconda.com/products/individual) installed.

### Clone the repo locally

First clone the repository using git clone in command line. You may need to install the git package (see [here](https://www.atlassian.com/git/tutorials/install-git)):

```
git clone https://github.com/luigiquiros/inventa.git
```

Create a new conda environment to avoid clashes:

```
conda env create -n inventa -f inventa/environment.yml
```

Then, activate the environment: 

```
conda activate inventa
```

If you need to update the environment run: 
[make sure to run the following to keep the depencies versions]

```
conda env update --file environment.yml
```

If you have an error, try installing `scikit-bio` from `conda-forge` before create the environment:

```
conda install -c conda-forge scikit-bio
```

Also, is possible to have the following error: 'Microsoft Visual C++ 14.0 or greater is required (install Microsoft C++ Build Toold from https://visualstudio.microsoft.com/visual-cpp-build-tools/)'. To solve it please make sure you include the installation of the workload Desktop development with C++


## [Continue to Getting Started](getting-started.md) 

### [Back to home page](index.md)
