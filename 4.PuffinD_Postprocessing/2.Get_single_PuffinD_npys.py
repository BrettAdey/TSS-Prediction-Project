import numpy as np
import pandas as pd
import os
import glob

# Define the base directory and the categories - must match directory structure and number of .npy files inside
# Run once for each category (run HPRT1 variants together)
base_dir = '/puffin_outputs'
categories = {
    'T2T Holdout Reversed': 4315 
    #'T2T Holdout Forward': 4315 
    #'Mononucleotide Global Shuffled Non-repeats': 31161 
    #'Mononucleotide Global Shuffled Repeats': 31161 
    #'Local Mononucleotide Shuffled': 31161 
    #'Global Mononucleotide Shuffled': 31161
    #'Global Dinucleotide Shuffled': 31161
    #'Local Dinucleotide Shuffled': 31161
    #'T2T Reversed Repeats': 31161
    #'T2T Reversed Nonrepeats': 31161
    #'T2T Forward': 31161,
    #'T2T Reverse': 31161,
    #'HPRT1': 1,
    #'HPRT1R': 1,
    #'HPRT1RnoCpG': 1,
}

# Function to load numpy arrays from a given directory and select 'GRO_CAP (forward strand)'
def load_numpy_arrays(category, num_files):
    arrays = []
    filenames = []
    category_dir = os.path.join(base_dir, category)
    
    for npy_file in glob.glob(os.path.join(category_dir, '*.npy'))[:num_files]:
        data = np.load(npy_file)
        selected_assay = data[:, 3, :]  # Selecting 'GRO_CAP (forward strand)' (index 3, ":" for all)
        arrays.append(selected_assay)
        filenames.append(os.path.basename(npy_file))
        
    return arrays, filenames


# Main script
if __name__ == '__main__':
    
    all_data = {} #Create a dictionary to hold all the data
    
    for category, num_files in categories.items():
        # Load the numpy arrays and their filenames
        arrays, filenames = load_numpy_arrays(category, num_files)
        
        #Store data in the dictionary
        all_data[category] = {'arrays': arrays, 'filenames': filenames}

#np.save('HPRT1.npy', all_data)
#np.save('T2T_Forward.npy', all_data)
#np.save('T2T_Reversed.npy', all_data)
#np.save('Reversed_Nonrepeats.npy', all_data)
#np.save('Reversed_Repeats.npy', all_data)
#np.save('Local_Dinuc.npy', all_data)
#np.save('Global_Dinuc.npy', all_data)
#np.save('Global_Mononuc.npy', all_data)
#np.save('Local_Mononuc.npy', all_data)
#np.save('Mononuc_Glob_Repeats.npy', all_data)
#np.save('Mononuc_Glob_Nonrepeats.npy', all_data)
#np.save('Holdout_Forward.npy', all_data)
np.save('Holdout_Reversed.npy', all_data)

