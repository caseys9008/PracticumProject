# ==================================================================================================================
# This script will take in a FASTA file and a list of residues to be preserved.  A file with the same
# name as the fasta file will be generated that contains fasta sequences of all other of all other possible options
# with the specified conserved amino acids unchanged

# the output files will be shuffled versions of the input fasta with the same residues in the conserved positions

# USAGE: python3 ShufflePeptide.py <<FASTA_file>> <<conserved_residues>> <<output_folder_path>>
# --> conserved_residues should be separated by commas, no spaces (conserved residues are optional)
# --> Ex.) 2,6,7
# ==================================================================================================================
import sys
import os
from random import sample
import math
import itertools
import pprint
# ===================================================================================================================
class FASTA_Handler:

    def parse_fasta(self, infile_fasta,
                    type):  # type = "aa' for amino acid (protein seq) and 'na' for nucleic acid (DNA or RNA)
        valid_na = ('A', 'C', 'G', 'T', 'U', 'R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N', '-')
        valid_aa = (
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'Y', 'Z', 'X', '*', '-')
        description = ""
        sequence = ""
        warnings = ""

        first = True  # keeps track of if we are looking at the first line in the file or not.

        with open(infile_fasta, 'r') as in_file:  # ensures that th file is closed once block is finished...
            for line in in_file:
                if line.startswith(">") and first:  # means it is the first sequence in the file
                    # (only need to set description string)
                    description = line.rstrip()[1:]
                    first = False

                elif first and not line.startswith(">"):  # if the first line is not in the correct format...
                    # need to raise an exception here because parser will not parse correctly if first line is invalid
                    raise Exception("Invalid input: first line must begin with '>'")

                elif line.startswith(
                        ">") and not first:  # if the line is a description but not for the first sequence...
                    # yield the previously recorded description and sequence
                    prev_seq = (description, sequence)
                    yield prev_seq

                    # record the new description line and reset the sequence variable
                    description = line.rstrip()[1:]
                    sequence = ""

                else:
                    # check for valid characters in the sequences
                    for each in line.rstrip():  # rstrip removes the white spaces
                        if type == "na":  # invalid characters will be different depending type of sequence entered
                            if each not in valid_na:
                                warnings += "- Invalid character found: character = '" + each + "'; sequence description = " \
                                            + description + "\n"
                        elif type == "aa":
                            if each not in valid_aa:
                                warnings += "- Invalid character found: character = '" + each + "'; sequence description = " \
                                            + description + "\n"
                    sequence += line.rstrip()  # add this line of sequence to the cumulative sequence string

        # print out all the warnings collected for the file
        last_seq = (description, sequence)  # This makes sure that the last sequence is included as well
        yield last_seq
        print(warnings)

    def format_tuple(self, result):
        x, y = result
        return x + " | " + y

# ==================================================================================================================
myFastaHandler = FASTA_Handler()

# save the command line input to variables
input_file = sys.argv[1]
fasta_desc = ""
input_fasta = ""
conserved_index_list = []
fasta_sequence = ""
output_folder_path = ""

# HANDLE COMMAND LINE ARGUMENTS -----------------------------------------------------------------------------------
if len(sys.argv) == 3:  # could be either only an input and output file or only an input and a list of
    if sys.argv[2].split(',')[0].isnumeric():  # if the second argument was a number
        for each in sys.argv[2].split(','):
            conserved_index_list.append(int(each) - 1)  # turns residue number into residue index
    else:
        output_folder_path = sys.argv[2]

if len(sys.argv) == 4:  # means there will be both conserved residues and an output file to deal with
    for each in sys.argv[2].split(','):
        conserved_index_list.append(int(each) - 1)  # turns residue number into residue index
    output_folder_path = sys.argv[3]

for (desc, seq) in myFastaHandler.parse_fasta(input_file, "aa"):  # there will only be one
    fasta_desc = desc
    input_fasta = seq
input_fasta_list = list(input_fasta)


# MAKE LIST OF ALL RESIDUES TO CONSERVE -----------------------------------------------------------------------------
not_conserved_residues = ''
for i in range(len(input_fasta_list)):
    if i not in conserved_index_list:
        not_conserved_residues += input_fasta_list[i]

# FIND ALL UNIQUE PERMUTATIONS OF THE NON CONSERVED RESIDUES (save to list) -----------------------------------------
all_nc_permutations = set(itertools.permutations(not_conserved_residues))

# LOOP THROUGH EACH PERMUTATION AND ADD BACK IN CONSERVED RESIDUES --------------------------------------------------
all_shuffled_strings = []
for perm in all_nc_permutations: # need to be a small enough number that you can loop through them... (test on RCC for sure)
    perm_string = "".join(perm)
    i = 0
    j = 0
    new_perm_string = ""
    while i < len(input_fasta):
        if i not in conserved_index_list:
            new_perm_string += perm_string[j]
            i += 1
            j += 1
        else: # (when i is in the conserved index list)
            new_perm_string += input_fasta[i]
            i += 1
    all_shuffled_strings.append(new_perm_string)

# SAVE ALL THE SHUFFLED SEQUENCES TO FASTA FILES IN A FOLDER (include reference file for quick lookup) --------------

ref_file_name = output_folder_path + "/" + (sys.argv[1]).split(".")[0] +  "_ALL_OUTPUT.fasta"  # the output folder path will contain the name of the folder
os.makedirs(os.path.dirname(ref_file_name), exist_ok=True)  # this should create the directory correctly
with open(ref_file_name, 'w+') as reference_file:  # have the reference file open the whole time
    for i in range(len(all_shuffled_strings)):
        #gen_file_name = output_folder_path + "/" + (sys.argv[1]).split(".")[0] + "_" + str(i + 1) + ".fasta"
        #with open(gen_file_name, 'w+') as seq_file:
            #seq_file.write(">" + (sys.argv[1]).split(".")[0] + "_" + str(i + 1) + "\n" + all_shuffled_strings[i])
        reference_file.write(">" + (sys.argv[1]).split(".")[0] + "_" + str(i + 1) + "\n" + all_shuffled_strings[i]+ "\n")




# == OTHER STUFF (to maybe use later) ============================================================================
#pprint.pprint(list(itertools.permutations(not_conserved_residues)))

# # calculate the total number of shuffles possible
# letter_counts = []
# for aa in not_conserved_residues:
#     print(aa)
#     count_tuple = (aa, not_conserved_residues.count(aa))
#     print(count_tuple)
#     if count_tuple not in letter_counts:
#         letter_counts.append(count_tuple)
#
# total_poss_wo_repeats = int(math.factorial(len(input_fasta)))
# total_denominator_calc = 1
# for letter_count in letter_counts:
#     aa, count = letter_count
#     total_denominator_calc = int(total_denominator_calc * math.factorial(count))
#
# total_pos_shuffles = total_poss_wo_repeats/total_denominator_calc
# print(int(total_pos_shuffles))  # should I just assume that python can handle a loop of this size?


