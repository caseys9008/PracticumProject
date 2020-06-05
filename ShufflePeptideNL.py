import sys
import os
from random import sample
import math
import itertools
import pprint
import shutil
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
# all_nc_permutations = set(itertools.permutations(not_conserved_residues))
# print("Saved all shuffles to array")


# Create the temp directory
directory = output_folder_path + "/"
print("Temp Directory " + str(directory))
os.makedirs(os.path.dirname(directory), exist_ok=True)  # this should create the directory correctly
print("Created the temp Directory")

# Create the final location directory
final_directory = "/scratch/midway2/caseys/P8_Shuffled/"   # ********************************
os.makedirs(os.path.dirname(final_directory), exist_ok=True)
print("Created the final Directory:  " + final_directory)


max_array_size = 150000000  # < -- the max size of a functioning python array is 536,870,912 elements.
array_below_500mil = []
num_output_files = 1
fastas_per_file = 150000000


# loop through each in the itertools.permutations
for each in itertools.permutations(input_fasta):

    res_conserved = True
    repeat = False

    # check that it has the conserved residues
    for index in conserved_index_list:
        if not each[index] == input_fasta[index]:
            res_conserved = False

    # if it does have the conserved residues and the array does not already contain the value... then add to array
    if res_conserved and not repeat:
        array_below_500mil.append(each)

    # if the array is the max size ...
    if len(array_below_500mil) == max_array_size:
        print("max size array reached... being handled now")
        # convert to a unique value set and loop through and add to files...
        unique_values = list(set(array_below_500mil))
        array_below_500mil = [] # reset the array to be empty (this should help to conserve memory) -- I hope

        # loop through and print the unique values to folders ...

        while len(unique_values) > 0:  # I will be deleting the values as I print them to files
            i = 0
            output_file_name = output_folder_path + "/" + (sys.argv[1]).split(".")[0] + "_OUTPUT_" + str(num_output_files)+ ".fasta"
            with open(output_file_name, 'w+') as output_file:
                for seq in unique_values[:fastas_per_file]:
                    output_file.write(">" + (sys.argv[1]).split(".")[0] + "_" + str((i + 1 + (fastas_per_file * (num_output_files - 1)))) +
                                      "\n" + ''.join(unique_values[i]) + "\n")
                    i += 1
            print("Wrote output file " + str(num_output_files))
            shutil.move(output_file_name, "/scratch/midway2/caseys/P8_Shuffled/"  + (sys.argv[1]).split(".")[0] + "_OUTPUT_" + str(num_output_files)+ ".fasta")
            print("file moved to final location")
            unique_values = unique_values[fastas_per_file:]
            num_output_files += 1


# if there are less than the max array size in the array
    # either to begin with or leftover --> will be handled by this block of code

# get the unique values
print("Done handling all arrays of max size")
unique_values = list(set(array_below_500mil))
array_below_500mil = [] # reset the array to save memory

while len(unique_values) > 0:  # I will be deleting the values as I print them to files
    i = 0
    # start = ((out_file_count - 1) * num_per_file) + 1
    # stop = start + num_per_file - 1
    output_file_name = output_folder_path + "/" + (sys.argv[1]).split(".")[0] + "_OUTPUT_" + str(num_output_files)+ ".fasta"
    with open(output_file_name, 'w+') as output_file:
        for seq in unique_values[:fastas_per_file]:
            output_file.write(">" + (sys.argv[1]).split(".")[0] + "_" +
                                      str((i + 1 + (fastas_per_file * (num_output_files - 1)))) +
                                      "\n" + ''.join(unique_values[i]) + "\n")
            i += 1
    print("Wrote output file " + str(num_output_files))
    shutil.move(output_file_name, "/scratch/midway2/caseys/P8_Shuffled/" + (sys.argv[1]).split(".")[0] + "_OUTPUT_" + str(num_output_files) + ".fasta")
    print("file moved to final location (the end)")
    unique_values = unique_values[fastas_per_file:]
    num_output_files += 1

os.rmdir(directory)
print("Removed the temporary directory")

print("Done Running")




#
# # == OTHER STUFF (to maybe use later) ============================================================================
# #pprint.pprint(list(itertools.permutations(not_conserved_residues)))
#
# # # calculate the total number of shuffles possible
# # letter_counts = []
# # for aa in not_conserved_residues:
# #     print(aa)
# #     count_tuple = (aa, not_conserved_residues.count(aa))
# #     print(count_tuple)
# #     if count_tuple not in letter_counts:
# #         letter_counts.append(count_tuple)
# #
# # total_poss_wo_repeats = int(math.factorial(len(input_fasta)))
# # total_denominator_calc = 1
# # for letter_count in letter_counts:
# #     aa, count = letter_count
# #     total_denominator_calc = int(total_denominator_calc * math.factorial(count))
# #
# # total_pos_shuffles = total_poss_wo_repeats/total_denominator_calc
# # print(int(total_pos_shuffles))  # should I just assume that python can handle a loop of this size?


