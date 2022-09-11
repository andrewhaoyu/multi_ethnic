# Reformat summary statistics and run SBayesR.
# All outputs (including the reformatted summary statistics and a python log file) are in the output_file_dire
# The output result file with posterior effect size is at os.path.join(output_file_dire, output_file_path+".snpRes")

import argparse
import os
import sys
import pandas as pd
import csv

## parse the command line argument
parser = argparse.ArgumentParser(description="Run SBayesR and output the posterior effect sizes.")
parser.add_argument('--race', required=True, type=str, dest="race", 
                    help="the race parameter of summary data")
parser.add_argument('--rho', required=True, type=int, dest="rho", 
                    help="the rho parameter of summary data")
parser.add_argument('--GA', required=True, type=int, dest="GA", 
                    help="the GA parameter of summary data")
parser.add_argument('--size', required=True, type=int, dest="size", 
                    help="the size parameter of summary data")
parser.add_argument('--rep', required=True, type=int, dest="rep", 
                    help="the rep parameter of summary data")
parser.add_argument('--use-public-ldm', default=False, action='store_true', dest="use_public_ldm", 
                    help="If specified, then the program will perform formatting differently")
args = parser.parse_args()

## input parameters
race = args.race
size = args.size
rho = args.rho
GA = args.GA
rep = args.rep

if size == 1:
  N = 15 * 1000
elif size == 2:
  N = 45 * 1000
elif size == 3:
  N = 80 * 1000
else:    
  N = 100 * 1000

## input file path
# summary data (with columns [CHR, POS, SNP_ID, REF, ALT, REF_FRQ PVAL, BETA, SE, N])
# EXAMPLE 1: separate (The following fillers must be included)
is_single_chr_sum = True
sumdata_filepath = "/dcl01/chatterj/data/jin/prs/simulation/{race}/sumdata/megasum-rho{rho}-size{size}-rep{rep}-GA{GA}-chr{chr_num}.txt"
# EXAMPLE 2: combined (absolute filepath with no fillers)
# is_single_chr_sum = False
# sumdata_filepath = "/dcl01/chatterj/data/jin/prs/simulation/{race}/sumdata/megasum-rho{rho}-size{size}-rep{rep}-GA{GA}.txt".format(race=race, rho=rho, size=size, rep=rep, GA=GA) # fake path

## directory of the LD matrix
ld_directory = "/dcs04/nilanjan/data/wlu/PolyPred/SBayesR/ukbEURu_hm3_shrunk_sparse"
ld_mat_file_prefix = "ukbEURu_hm3_{chr_num}_v3_50k.ldm.sparse{file_type}" # (.info / .bin); chr_num = 'all', chr1, chr2
mldm_list_filename = "ukbEURu_hm3_sparse_mldm_list.txt" # will be rewritten with new filepath

## gctb path
gctb_path = "/dcs04/nilanjan/data/wlu/PolyPred/SBayesR/gctb_2.03beta_Linux/gctb"

## output file path
# All outputs (including a python log file) are in the output_file_dire
# The output result file with posterior effect size is at output_file_path+".snpRes"
if args.use_public_ldm:
  print("Use public ldm from {}".format(ld_directory))
output_file_dire = "/dcs04/nilanjan/data/wlu/PolyPred/output/SBayesR_public_ldm/{race}/rho{rho}-size{size}-rep{rep}-GA{GA}".format(race=race,rho=rho,size=size,rep=rep,GA=GA)
else:
  print("Use self ldm from {}".format(ld_directory))
output_file_dire = "/dcs04/nilanjan/data/wlu/PolyPred/output/SBayesR_self_ldm/{race}/rho{rho}-size{size}-rep{rep}-GA{GA}".format(race=race,rho=rho,size=size,rep=rep,GA=GA)

output_file_path = os.path.join(output_file_dire, "SBayesR_result")
python_log = os.path.join(output_file_dire, os.path.basename(__file__[:-3])+".out")


############ Reidrect Output & Sanity Check ###############
# Helper function to run shell command
def run_shell_cmd(shell_cmd):
  print("Running shell command $ " + shell_cmd)
ret = os.popen(shell_cmd)
print(ret.read())

# prep
shell_cmd = "mkdir -p {}".format(output_file_dire)
run_shell_cmd(shell_cmd)

# Redirect stdout
class RedirectOutput:
  
  def __init__(self, filename):
  self.fp = open(filename, 'w')

def write(self, msg):
  self.fp.write(msg)
self.fp.flush()

def flush(self):
  self.fp.flush()

python_out_fp = RedirectOutput(python_log)
sys.stdout = python_out_fp
sys.stderr = python_out_fp

# sanity check
if not os.path.exists(gctb_path):
  msg = "GCTB path" + gctb + " does not exist"
raise Exception(msg)

if is_single_chr_sum:
  for i in range(1, 23):
  curr = sumdata_filepath.format(race=race,rho=rho,size=size,rep=rep,GA=GA,chr_num=i)
if not os.path.exists(curr):
  msg = "Summary stats " + curr + " does not exist"
raise Exception(msg)
else:
  if not os.path.exists(sumdata_filepath):
  msg = "Summary stats " + sumdata_filepath + " does not exist"
raise Exception(msg)

for chr_num in range(1, 23):
  for file_type in [".bin", ".info"]:
  ld_file = ld_mat_file_prefix.format(chr_num="chr"+str(chr_num), file_type=file_type)
ld_path = os.path.join(ld_directory, ld_file)
if not os.path.exists(ld_path):
  msg = "LD matrix file {} does not exist".format(ld_path)
raise Exception(msg)


# ############# Formatting Input Files ################
def read_ld_info(args, path):
  with open(path, 'r') as fp:
  info_lines = fp.readlines()
fp.close()
if not args.use_public_ldm:
  info_lines.pop(0)
return info_lines

print("##################### Formatting Input Files #####################")

all_ld_file = ld_mat_file_prefix.format(chr_num="all", file_type=".info")
all_ld_path = os.path.join(ld_directory, all_ld_file)
if not os.path.exists(all_ld_path):
  print("Create and append to", all_ld_path)
for chr_num in range(1, 23):
  ld_file = ld_mat_file_prefix.format(chr_num="chr"+str(chr_num), file_type=".info")
ld_path = os.path.join(ld_directory, ld_file)
info_lines = read_ld_info(args, ld_path)
with open(all_ld_path, 'a+') as fp:
  for curr_line in info_lines:
  fp.write(curr_line)
if curr_line[-1] != "\n":
  fp.write("\n")
fp.close()
del info_lines

## build the list of LD matrix file paths
mldm_list_path = os.path.join(ld_directory, mldm_list_filename)
with open(mldm_list_path, 'w') as fp:
  for chr_num in range(1, 23):
  ld_file = ld_mat_file_prefix.format(chr_num="chr"+str(chr_num), file_type="")
ld_path = os.path.join(ld_directory, ld_file) + "\n"
fp.write(ld_path)
print("Write the list of LD matrix file paths in {}".format(mldm_list_path))

## format the summary data
# 1. to satisfy SBayesR requirement
# 2. to have consistent SNPID with HapMap3 used for the LD matrix (just in case if this would affect the result)

# sumdata_raw: [CHR, POS, SNP_ID, REF, ALT, REF_FRQ, PVAL, BETA, SE, N, CHR_join, POS_join]
if is_single_chr_sum:
  sumdata_list = []
for i in range(1, 23):
  curr = sumdata_filepath.format(race=race,rho=rho,size=size,rep=rep,GA=GA,chr_num=i)
sumdata_curr = pd.read_table(curr)
sumdata_list.append(sumdata_curr)
sumdata_raw = pd.concat(sumdata_list, axis=0)
else:
  sumdata_raw = pd.read_table(sumdata_filepath)
sumdata_raw["CHR_join"] = sumdata_raw['CHR'].astype(int)
sumdata_raw["POS_join"] = sumdata_raw['POS'].astype(int)
SNP_ID = sumdata_raw["SNP_ID"]
sumdata_raw["SNP_ID"] = [i.strip().split(":")[0].strip() for i in SNP_ID]
del SNP_ID
# sumdata_raw: [CHR, POS, SNP_ID, A2, A1, A1FREQ, PVAL, BETA, SE, N, CHR_join, POS_join]
sumdata_raw = sumdata_raw.rename(columns={"REF": "A2", "ALT": "A1", "REF_FRQ": "A1FREQ"})
sumdata_raw['A1FREQ'] = 1 - sumdata_raw['A1FREQ']
print("Reformat summary data: {} SNPs found in summary data".format(len(sumdata_raw)))

# HapMap info data
# ld_mat_df: ['CHR_join', 'LD_SNP_ID', 'POS_join', "LD_A1", "LD_A2"]
ld_mat_df = pd.read_csv(all_ld_path, sep='\s+', usecols=[0, 1, 3, 4, 5], names=['CHR_join', 'LD_SNP_ID', 'POS_join', "LD_A1", "LD_A2"])
print("Reformat summary data: {} SNPs found in LD matrix".format(len(ld_mat_df)))
ld_mat_df = ld_mat_df.drop(index = ld_mat_df.loc[~ld_mat_df['LD_A1'].isin(['A', 'T', 'G', 'C'])].index)
ld_mat_df = ld_mat_df.drop(index = ld_mat_df.loc[~ld_mat_df['LD_A2'].isin(['A', 'T', 'G', 'C'])].index)
ld_mat_df = ld_mat_df.drop(index = ld_mat_df.loc[ld_mat_df['POS_join'].isna()].index)
ld_mat_df = ld_mat_df.drop(index = ld_mat_df.loc[ld_mat_df['LD_SNP_ID'].isna()].index)
ld_mat_df = ld_mat_df.drop(index = ld_mat_df.loc[ld_mat_df['CHR_join'].isna()].index)
ld_mat_df["CHR_join"] = ld_mat_df['CHR_join'].astype(int)
ld_mat_df["POS_join"] = ld_mat_df['POS_join'].astype(int)
print("Reformat summary data: {} SNPs in LD matrix are valid".format(len(ld_mat_df)))

# perform reformatting
# sumdata_raw: [CHR, POS, SNP_ID, A2, A1, A1FREQ, PVAL, BETA, SE, N, CHR_join, POS_join, LD_SNP_ID, LD_A1, LD_A2]
sumdata_raw = sumdata_raw.merge(ld_mat_df)
del ld_mat_df

print("Reformat summary data: {} SNPs matched with the LD matrix".format(len(sumdata_raw)))
changed_rows = sumdata_raw['SNP_ID']!=sumdata_raw['LD_SNP_ID']
sumdata_raw.loc[changed_rows, 'SNP_ID'] = sumdata_raw.loc[changed_rows, 'LD_SNP_ID']
print("Reformat summary data: {} SNP IDs in the summary data are changed to be consistent with the LD matrix".format(sum(changed_rows)))
del changed_rows

matched_rows = (sumdata_raw['A1']==sumdata_raw['LD_A1']) & (sumdata_raw['A2']==sumdata_raw['LD_A2'])
swapped_rows = (sumdata_raw['A1']==sumdata_raw['LD_A2']) & (sumdata_raw['A2']==sumdata_raw['LD_A1'])
invalid_rows = ~(matched_rows | swapped_rows)
sumdata_raw.loc[swapped_rows, "A1"] = sumdata_raw.loc[swapped_rows, "LD_A1"] 
sumdata_raw.loc[swapped_rows, "A2"] = sumdata_raw.loc[swapped_rows, "LD_A2"] 
sumdata_raw.loc[swapped_rows, "A1FREQ"] = 1 - sumdata_raw.loc[swapped_rows, "A1FREQ"]
sumdata_raw.loc[swapped_rows, "BETA"] = -1 * sumdata_raw.loc[swapped_rows, "BETA"]
sumdata_raw = sumdata_raw.drop(sumdata_raw[invalid_rows].index)
sumdata_raw['Z'] = sumdata_raw['BETA'] / sumdata_raw['SE']
print("Reformat summary data: {} SNPs in summary data were deleted due to unmatched alleles".format(sum(invalid_rows)))
print("Reformat summary data: {} SNPs in summary data were aligned with the LD matrix".format(sum(swapped_rows)))

# save .ma file
sumdata_raw = sumdata_raw[['SNP_ID', 'A1', 'A2', 'A1FREQ', 'BETA', 'SE', 'PVAL', 'N']]
ma_sumstats_file_path = os.path.join(output_file_dire, "sumstats.ma")
sumdata_raw.to_csv(ma_sumstats_file_path, sep=" ", header=True, index=False, doublequote=False, quoting=csv.QUOTE_NONE)

print("Reformat summary data: {} SNPs in summary data remained".format(len(sumdata_raw)))
print("Reformat summary data: stored at {}".format(ma_sumstats_file_path))
print()

del matched_rows, swapped_rows, invalid_rows, sumdata_raw

# ############# Running SBayesR ############
print("##################### Running SBayesR #####################")

final_output_file = output_file_path+"{}.snpRes"
num = 0
gamma_last = 2
while not os.path.exists(final_output_file.format(num)):
  num += 1
gamma_last = gamma_last / 2
if gamma_last < 1e-6:
  break
print("\nNo. {}".format(num))
shell_cmd = """
    {gctb}  --sbayes R \
            --mldm {mldm} \
            --gwas-summary {summary} \
            --out {out} \
            --pi 0.95,0.02,0.02,0.01 \
            --gamma 0,0.01,0.1,{gamma} \
            --seed 0 --chain-length 10000 --burn-in 4000 --out-freq 10 --thin 10 \
             > {log} 2>&1 
    """.format(gctb=gctb_path, 
               mldm=mldm_list_path, 
               summary=ma_sumstats_file_path, 
               out=output_file_path+str(num),
               log=output_file_path+"{}.log".format(num),
               gamma=gamma_last)
run_shell_cmd(shell_cmd)

if not os.path.exists(final_output_file.format(num)):
  print("\nNo. {} (remove the last mixture component)".format(num))
shell_cmd = """
    {gctb}  --sbayes R \
            --mldm {mldm} \
            --gwas-summary {summary} \
            --out {out} \
            --pi 0.95,0.025,0.025 \
            --gamma 0,0.01,0.1 \
            --seed 0 --chain-length 10000 --burn-in 4000 --out-freq 10 --thin 10 \
             > {log} 2>&1 
    """.format(gctb=gctb_path, 
               mldm=mldm_list_path, 
               summary=ma_sumstats_file_path, 
               out=output_file_path+str(num),
               log=output_file_path+"{}.log".format(num))
run_shell_cmd(shell_cmd)

print()
if os.path.exists(final_output_file.format(num)):
  shell_cmd = "cp {orig} {target}".format(orig=final_output_file.format(num), target=final_output_file.format(""))
run_shell_cmd(shell_cmd)
print("Final output", final_output_file.format(""))
else:
  print("Failed")
print("----- END -----")