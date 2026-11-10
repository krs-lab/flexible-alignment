### Originally developed by Dr. Karl R Schmitz 

import os
import subprocess
import time
import multiprocessing
import math
import io
import os.path
import shutil
from contextlib import redirect_stdout, redirect_stderr
import pandas as pd
from tqdm import tqdm


########################################
### This script will run a pairwise comparison of
###  Lpn effector domains in one folder against those
###  in a second folder. The results will be saved to a 
###  csv file as a half-diagonal matrix.
###
### One important feature of this script is
###  that the output files for each FATCAT run will be
###  deleted after extracting the E-score.
###
### Optionally, one can enter a cutoff E-score
###  below above which the output can be saved in a 
###  separate folder.
###

### Do we have problems reading from explicity paths?
###
### Query domains are read from:
### /home/krs/python/FATCAT/46-lpn_vs_all_lpn/FATCAT_analysis
query_domains_path = "FATCAT_analysis"

### Legionella effector list
### /home/krs/python/FATCAT/lpn_effectors (symlink)
lpn_effectors_path = "lpn_effectors"

### Output csv file is written to:
### /home/krs/python/FATCAT/46-lpn_vs_all_lpn/output
output_csv_path = "output"

### Make a ramdisk beforehand at this point.
### This should limit readwrites to the disk
output_tmp_path = "ramdisk"

### High-scoring files (if flagged) go in:
### 
high_scoring_path = "high_scores"

########################################
print("This program will run a FATCAT comparison of query"
      "pdb files in one target directory to pdb files in"
      "a second directory. The results will be saved"
      "as a matrix in csv format.")


num_cores = int(input("Please enter the number of cores to use:\n"))

#########################################
### Start by getting a list of the query domains
### -  as filenames

query_pdb_filenames = []
count = 0

for filename in os.listdir(query_domains_path):
    if filename.endswith(".pdb"):
        query_pdb_filenames.append(filename)
        count += 1

print("We found",count,"pdb files in",query_domains_path,".")


#########################################
### Start by getting a list of the effector domains
### -  as filenames

effector_pdb_filenames = []
count = 0

for filename in os.listdir(lpn_effectors_path):
    if filename.endswith(".pdb"):
        effector_pdb_filenames.append(filename)
        count += 1

print("We found",count,"pdb files in",lpn_effectors_path,".")


num_of_comparisons = int(len(query_pdb_filenames) * len(effector_pdb_filenames))
print("We will make", num_of_comparisons, "pairwise comparisons.")


##########################################
### Do we want to set a high-scoring cutoff to save log/pdb files?
### If not, just leave flag as "0"
print("Do you want to set a high-scoring cutoff E-value? Copies above the threshold"
"will be saved to the high-scores directory.")
log_score = input("If so, please enter the score as -log10 value:\n")
cutoff_score = 10 ** (-1 * float(log_score))
print("That corresponds to a cutoff of",cutoff_score)

##########################################
### Now we start doing fatcat runs
###  and create the pandas dataframe to hold results.

# Start with an empty DataFrame with named columns, rows
scores_df = pd.DataFrame(index=query_pdb_filenames, columns=effector_pdb_filenames)

####################################
## Define a function to run Fatcat jobs.
## pdb_directory = directory with PDB files
## pdb_1_fn = first PDB filename
## pdb_2_fn = second PDB filename
## 
## 

## Note: this program will capture FATCAT output using StringIO
##  
def run_fatcat(query_pdb_filenames, pdb_1_fn, effector_pdb_filenames, pdb_2_fn, high_score_output_path, ramdisk, cutoff):
    ## strip off ".pdb" to get the prefix
    pdb_1_prefix = pdb_1_fn[:-4]
    pdb_2_prefix = pdb_2_fn[:-4]

    #command = rf"FATCAT -p1 {pdb_directory}/{pdb_1_fn} -p2 {pdb_directory}/{pdb_2_fn} -o {ramdisk}/{pdb_1_prefix}_{pdb_2_prefix} -m -ac -t"
    #command = f"FATCAT -p1 '{pdb_directory}/{pdb_1_fn}' -p2 '{pdb_directory}/{pdb_2_fn}' -m -ac -t"
    command = [
        "FATCAT",
        "-p1", f"{query_pdb_filenames}/{pdb_1_fn}",
        "-p2", f"{effector_pdb_filenames}/{pdb_2_fn}",
        "-o", f"{ramdisk}/{pdb_1_prefix}_{pdb_2_prefix}",
        "-m", "-ac", "-t"
    ]
    cwd = os.getcwd()

    # Execute the command
    subprocess.run(command) 

    # read the output file for the P-value
    filename = f"{ramdisk}/{pdb_1_prefix}_{pdb_2_prefix}.aln"
    #print("Opening the alignment file",filename)
    with open(filename, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if "P-value" in line:
                    p_value = float(line.split()[1])
                    #print("-P-value appears to be",p_value)

    
    
    # If the results met the threshold, copy output to the high score directory
    if (p_value <= cutoff):
        query = pdb_1_prefix + "_" + pdb_2_prefix
        #print("--In IF: Met threshold!!!!!!!!!!!! Looking for",query)
        #if (os.path.isfile(f"{ramdisk}/{query}.aln")):
        #    print("--Which currently exists")
        #else:
        #    print("--But does not exist?")
        
        # Look for all output files that match this basename
        # print("Ramdisk contents:",os.listdir(ramdisk))
        for out_file in os.listdir(ramdisk):
            basename = os.path.splitext(os.path.basename(out_file))[0]
            #print("---Looking for",basename,"in",out_file)
            if query in basename:
                #print("---Found match!!")
                try:
                    #print("Trying to copy then delete...")
                    src_path = os.path.join(ramdisk, out_file)
                    dst_path = os.path.join(high_score_output_path, out_file)
                    shutil.copy2(src_path, dst_path) 
                    os.remove(src_path)
                except Exception as error:
                    print("Error:",error)
                #print("End of if...")

        #print("Outside of for")
        #if (os.path.isfile(f"{ramdisk}/{query}.aln")):
        #    print("....We didn't delete it!")
        #else:
        #    print("....And we managed to delete it!")
        #print("End of score test if")

    # Otherwise, just delete the output files
    else:
        query = pdb_1_prefix + "_" + pdb_2_prefix
        #print("--Didn't meet threshold")
        #if (os.path.isfile(f"{ramdisk}/{query}.aln")):
            #print("--Which currently exists")
        #else:
            #print("--But does not exist?")

        # Look for all output files that match this basename
        for out_file in os.listdir(ramdisk):
            basename = os.path.splitext(os.path.basename(out_file))[0]
            #print("Looking at",basename)
            if query in basename:
                #print("Found match!!")
                src_path = os.path.join(ramdisk, out_file)
                os.remove(src_path)

        #if (os.path.isfile(f"{ramdisk}/{query}.aln")):
            #print("....We didn't delete it!")
        #else:
            #print("....And we managed to delete it!")

    #print("-fun_fatcat is going to return:",p_value)
    return p_value


########################################
### Try and do this multithreaded
### Must first create a queue...
queue = []

### Iterate through the set of comparisons

# First, establish the queue of pdb comparisons
#print("prior to loop")
for row_pdb in query_pdb_filenames:
    #print("loop layer 1, with i",i)
    for col_pdb in effector_pdb_filenames:
        #print("loop layer 2, with j",j)
        #print("Queuing on",i,j,"with",row_pdb,"vs",col_pdb)
        queue.append([row_pdb,col_pdb])


def worker(pdb_1_directory, pdb_1_fn, pdb_2_directory, pdb_2_fn, high_score_output_path, ramdisk, cutoff, 
            return_dict):
    """
    Worker function for the multiprocessing pool.

    Args:
        pdb_1_directory (str): Path to the directory containing PDB files.
        pdb_1_fn (str): Filename of the first PDB file.
        pdb_2_directory (str): Path to directory
        pdb_2_fn (str): Filename of the second PDB file.
        high_score_output_path (str): Path to the permanent output directory.
        ramdisk (str): Path to the ramdisk directory.
        cutoff (float): P-value cutoff for saving output.
        return_dict (dict): Dictionary to store the results.
    """
    score = run_fatcat(pdb_1_directory, pdb_1_fn, pdb_2_directory, pdb_2_fn, 
                       high_score_output_path, ramdisk, cutoff)
    #print(">>>Returned Score is",score)
    return_dict[(pdb_1_fn, pdb_2_fn)] = score

pbar = tqdm(total=len(queue))
manager = multiprocessing.Manager()
results_dict = manager.dict()

#run_fatcat(pdb_directory, pdb_1_fn, pdb_2_fn, output_path, cutoff):
with multiprocessing.Pool(processes=num_cores) as pool:
    for next_pdb_pair in queue:
        pool.apply_async(worker, 
                        args=(query_domains_path, next_pdb_pair[0], lpn_effectors_path, next_pdb_pair[1], 
                                high_scoring_path, output_tmp_path, cutoff_score, 
                                results_dict))
        pbar.update(1)

    pool.close()
    pool.join()

# Populate the scores_df
for (pdb1, pdb2), score in results_dict.items():
    scores_df.loc[pdb1, pdb2] = score

#scores_df.loc[row_pdb,col_pdb] = run_fatcat(effector_domains_path,pair[0],pair[1],high_scoring_path,output_tmp_path,cutoff_score)


scores_df.to_csv(output_csv_path + "/lpn-vs-lpn_scores.csv") 
