# Oligomer-Classifier
Oligomer classifier toy model 

This is a simple oligomer classifier. The script ('polymer_generator.py') is used to create random SMILES of three types of oligomers (peptides, plastics and oligosaccharides) and to generate a dataset that includes each oligomer and its corresponding Morgan fingerprint as a 2048 bit vector (check: https://www.rdkit.org/docs/GettingStartedInPython.html#fingerprinting-and-molecular-similarity). 

The notebook ('Oligomer_classifier') contains the model and some additional tests. 

In order to run the script, you have to install the RDKit package (see: https://www.rdkit.org/docs/Install.html#). The neural network is created using Keras and Tensorflow.
