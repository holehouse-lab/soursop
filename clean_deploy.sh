#!/bin/zsh

# source (note this automatically activates the appropriate conda 
# environment!)
source ~/.zshrc

# double check though...
which_pip=$(which pip)

# note will need to update this for linux installation...
if [ "/Users/alex/Python/conda/envs/neuron3/bin/pip" != ${which_pip} ]
then 
    echo "Could not find correct pip..."
    exit 1
else
    echo "Found correct pip..."
fi


echo "Uninstalling old CTraj..."
pip uninstall CTraj

echo "Installing new CTraj..."
pip install .
 
