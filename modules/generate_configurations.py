#%%
import numpy as np
import os
from time import time
import serial
import pandas as pd
from itertools import product, permutations


#%%

def generate_unique_combinations(num_vals):
    """
    Generate all possible cominations of current ratios (r1, r2, r3) between the three coils.
    The set of combinations fulfills fulfills following properties:
    - at least one value is 1, all other values are <= 1
    - each triplet of ratios is unique, i.e. the set contains each triplet only once
    - no two triplets are the equivalent up to their sign, i.e. multiplying a triplet by -1 
    does not produce a different element of the set
    - all possible combinations with num_vals values between -1 and 1 are contained in the set

    Arg:
    - num_vals: there are num_vals possible values of r1,r2,r3 that are equally 
    spaced between (incluively) -1 and 1. 
    num_vals must be odd to contain 0.

    Return: 
    -unique_combinations (ndarray): array containing all possible combinations of shape (number of combis, 3)

    Note:
    This may not be the most efficient implementation to generate the set, 
    better do not use it for large values of num_vals.
    """
    # generate all configurations in raw numbers
    generator_combinations = product([1], 
                                    np.linspace(-1,1, num=num_vals), 
                                    np.linspace(-1,1, num=num_vals))
    raw_ratios = np.array(list(generator_combinations))

    combinations_ratios = np.zeros((len(raw_ratios)*6,3))
    for i in range(len(raw_ratios)):
        combinations_ratios[i*6:(i+1)*6] = np.array(list(permutations(raw_ratios[i])))

    # remove duplicates
    unique_combinations = np.unique(combinations_ratios, axis=0)

    # remove combinations that are equal up to the sign
    indices_equivalent_pairs = []
    for i in range(len(unique_combinations)//2):
        for j in range(len(unique_combinations)):
            if np.all(unique_combinations[i] == -1*unique_combinations[j]):
                # print('{},{}: {}, {}'.format(i,j, unique_combinations[i], unique_combinations[j]))
                indices_equivalent_pairs.append(i)
    unique_combinations = np.delete(unique_combinations, indices_equivalent_pairs, axis=0)

    # reverse order for convenience to have 111 at the beginning 
    unique_combinations = np.flip(unique_combinations, axis=0)

    return unique_combinations

#%%
# set the number of values between (incl) -1 and 1 that should be considered 
num_vals = 2

for num_vals in range(2,10):
    # generate the set
    unique_combinations = generate_unique_combinations(num_vals)

    # save the combinations to csv file
    directory = '../config_files'


    df = pd.DataFrame({ 'ratio coil 1': unique_combinations[:,0], 
                        'ratio coil 2': unique_combinations[:,1], 
                        'ratio coil 3': unique_combinations[:,2]})

    output_file_name = 'configs_numvals{}_length{}.csv'.format(num_vals, len(unique_combinations))
    data_filepath = os.path.join(directory, output_file_name)
    df.to_csv(data_filepath, index=False, header=True)

# %%
