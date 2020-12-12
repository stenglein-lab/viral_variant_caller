### This directory contains files necessary to setup a server to run this pipeline

To setup a conda environment containing all of the dependencies necessary to run this pipeline, change to the directory containing these files and run:

```
./setup_environment.sh
```

If you'd like to update this conda environment with additional conda packages, add them to the .yaml file, change to the directory containing these files, and run:

```
conda env update --prefix=$HOME/variant_conda_environment --file variant_conda_environment.yaml
```

