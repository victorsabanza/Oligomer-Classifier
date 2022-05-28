
'''
Polymer Generator script. It generates a .csv file containing the SMILES of random oligomers and their Morgan fingerprint as a 2048 bit vector.
It also includes the class of the oligomer (peptide, plastic or oligosaccharide).
'''

from rdkit.Chem import AllChem
import numpy as np
import pandas as pd
from rdkit import Chem #RDKit Chemistry
from random import choice


#This class is used to generate the polymers
class RandomPolymerGenerator():
    
    '''Create a RandomPolymerGenerator object which creates and stores a 
    list of random polymers from a given list of monomers'''
    
    def __init__(self, monomer_list):
        
        self.monomers = monomer_list
        self.polymers = []
    
    
    def random_polymer(self, pol_len, fix = False):
        
        '''Generate one random polymer of length pol_len from
        the monomers list. It optionally corrects the string adding an
        extra oxigen atom at the beggining if fix = True'''
            
        pol = ''
        
        for i in range(pol_len):
            
            mon = choice(self.monomers)
            
            pol = pol + mon
        
        if fix:
            
            pol = 'O' + pol
            
        return pol
        
          
        
    def generate_polymers(self, length, n, correct = False):
        
        '''Add a n polymers of lengths = length to self.polymers'''
        
        new_polymers = [self.random_polymer(length, fix = correct) for i in range(n)]
        
        self.polymers = self.polymers + new_polymers
    
    
    def filter_duplicates(self):
        
        '''Filter duplicates from the polymers list'''

        self.polymers = set(self.polymers)
        


if __name__ == "__main__":
    
    print('Starting program')

    #Monomer lists
    #Create lists containing monomers (amino acids, plastic monomers and sugars). Concatenation of monomers generates the polymer string 
    #(notice that in the case of amino acids and sugars there is a missing 'O' at the beggining of the string. The purpose of this is concatenating 
    # the strings easily, then it will be corrected).

    aa = ['C(=O)[C@]([H])([H])N', 'C(=O)[C@]([H])(C)N', 'C(=O)[C@]([H])(CCC(O)=O)N', 'C(=O)[C@]([H])(CC(C)C)N', 'C(=O)C(C(C)(C))N', 'C(=O)[C@]([H])(CO)N', 'C(=O)C(CCCNC(N)=N)N', 'C(=O)[C@]([H])(CCC(N)=O)N', 'C(=O)[C@]([H])(C(O)C)N', 'C(=O)[C@]([H])(CC(N)=O)N', 'C(=O)[C@]([H])(CCSC)N', 'C(=O)[C@]([H])(CC1=CNC2=C1C=CC=C2)N', 'C(=O)[C@]([H])(CC(O)=O)N', 'C(=O)[C@]([H])(CC1=CNC=N1)N', 'C(=O)[C@]([H])(CC1=CC=CC=C1)N', 'C(=O)[C@]([H])(CC1=CC=C(O)C=C1)N', 'C(=O)[C@@]([H])(CS)N', 'C(=O)[C@]([H])(C(CC)C)N', 'C(=O)[C@]([H])(C1CCCN1)N', 'C(=O)[C@]([H])(C(C)C)N'] #Code amino acids as SMILES. Notice that we don't include the first -OH in order to make the joining simpler

    monomers = ['CC', 'C(C)C', 'C(Cl)C', 'C(CC)C', 'CC(c1ccccc1)']

    sugars = ['[C@]([C@@]([H])(CO)O[C@@](O)([H])[C@]([H])1O)([H])[C@@]1([H])O', '[C@]([C@@]([H])(CO)O[C@@](O)([H])[C@]([H])1O)([H])[C@@]1([H])N', '[C@]([C@@]([H])(CO)O[C@@](O)([H])[C@]([H])1N)([H])[C@@]1([H])O', '[C@]([C@]([H])(CO)O[C@](O)([H])[C@]([H])1O)([H])[C@@]1([H])O', '[C@]([C@]([H])(CO)O[C@](O)([H])[C@]([H])1O)([H])[C@]1([H])O', '[C@@]([C@]([H])(CO)O[C@](O)([H])[C@@]([H])1O)([H])[C@]1([H])O']


    
    #Generate peptides SMILES list
    
    print("Generating peptides")

    peptides = RandomPolymerGenerator(aa)

    for i in range(3,8):
        
        peptides.generate_polymers(i, 1400, correct=True)

    peptides.filter_duplicates()

    #Transform SMILES into rdkit mol objects
    pep_mols = [Chem.MolFromSmiles(i) for i in peptides.polymers]



    # Generate plastic SMILES. Same process than before.

    print("Generating plastics")

    plastics = RandomPolymerGenerator(monomers)

    for i in range(3,8):
        
        plastics.generate_polymers(i, 2300)

    plastics.filter_duplicates()

    plast_mols = [Chem.MolFromSmiles(i) for i in plastics.polymers]


    # Generate oligosaccharides SMILES
    
    print("Generating oligosaccharides")

    oligosac = RandomPolymerGenerator(sugars)

    for i in range(3,8):
        
        oligosac.generate_polymers(i, 2000, correct=True)

    oligosac.filter_duplicates() 

    oligosach_mols = [Chem.MolFromSmiles(i) for i in oligosac.polymers]


    #Dataset generation

    #Now we generate the dataset containing the Morgan Fingerprint (as a bit vector) for the peptides, oligosaccharides and oligomers 
    # adding their corresponding label (this will be the dependent variable for the algorithm). 

    # Create list with peptide smiles, peptide fingerprints and the corresponding labels and create a dataframe

    print("Generating full dataset")
    
    #convert set into list
    pep_smiles = list(peptides.polymers)
    
    #Create Morgan fingeprint from mol objects
    pep_fp = [np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 4, nBits = 2048)) for mol in pep_mols]

    #Create labels
    labels = ['peptide']*len(pep_fp)
    
    #Create dictionary and then DataFrame 
    pep_data = {'smiles': pep_smiles, 'fp': pep_fp, 'label': labels}

    pep_df = pd.DataFrame(pep_data)


    # Same for plastics

    plastic_smiles = list(plastics.polymers)

    plastic_fp = [np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 4, nBits = 2048)) for mol in plast_mols]

    pl_labels = ['plastic']*len(plastic_fp)

    plastic_data = {'smiles': plastic_smiles, 'fp': plastic_fp, 'label': pl_labels}

    plastic_df = pd.DataFrame(plastic_data)


    # Same for oligosaccharides

    oligosac_smiles = list(oligosac.polymers)

    oligosac_fp = [np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 4, nBits = 2048)) for mol in oligosach_mols]

    olig_labels = ['oligosaccharide']*len(oligosac_fp)

    olig_data = {'smiles': oligosac_smiles, 'fp': oligosac_fp, 'label': olig_labels}

    olig_df = pd.DataFrame(olig_data)



    # Now, combine the three dataframes

    merged_df = pd.concat((pep_df, plastic_df, olig_df), ignore_index=True)


    # Transform fp column into separate columns (one for each digit in the fingerprint) and then drop column with single fingerprint

    fps = pd.DataFrame(merged_df['fp'].tolist())

    merged_df.drop(labels='fp', axis=1)

   
    #Create final dataframe and save it 
    final_df = pd.concat((merged_df.drop(labels='fp', axis=1), fps), axis = 1)
     
    final_df = final_df.sample(frac=1).reset_index(drop=True) 

    print('Dataset generated correctly. Saving as polymer_dataset.csv')

    final_df.to_csv('polymers_dataset.csv')




