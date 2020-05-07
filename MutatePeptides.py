# ==================================================================================================================
# This script will take in a FASTA file and a list of residues to be preserved.  A file with the same
# name as the fasta file will be generated that contains fasta sequences of all other of all other possible options
# with the specified conserved amino acids unchanged

# USAGE: python3 MutatePeptides.py <<FASTA_file>> <<conserved_residues>> <<output_folder_path>>
# --> conserved_residues should be separated by commas, no spaces (conserved residues are optional)
# --> Ex.) 2,6,7
# ==================================================================================================================
import sys
import os
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

#myFastaHandler = fh.FASTA_Handler()  # make sure that FASTA_Handler.py is in the same directory... (change this?)
myFastaHandler = FASTA_Handler()

# save the command line input to variables
input_file = sys.argv[1]
fasta_desc = ""
input_fasta = ""
conserved_index_list = []
fasta_sequence = ""
output_folder_path = ""

if len(sys.argv) == 3:  # could be either only an input and output file or only an input and a list of
    if sys.argv[2].split(',')[0].isnumeric():  # if the second argument was a number
        for each in sys.argv[2].split(','):
            conserved_index_list.append(int(each) - 1)  # turns residue number into residue index
    else:
        output_folder_path = sys.argv[3]

if len(sys.argv) == 4:  # means there will be both conserved residues and an output file to deal with
    for each in sys.argv[2].split(','):
        conserved_index_list.append(int(each) - 1)  # turns residue number into residue index
    output_folder_path = sys.argv[3]

for (desc, seq) in myFastaHandler.parse_fasta(input_file, "aa"):  # there will only be one
    fasta_desc = desc
    input_fasta = seq
input_fasta_list = list(input_fasta)

# save a list of all possible residues (for mutating purposes)
residues = ("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
# should I be including x for selenomethionine?


# alter the residues and save the results.
mutated_sequences = []
for i in range(len(input_fasta_list)):

    # INDEX NOT CONSERVED
    if i not in conserved_index_list:
        if i == 0:
            for each in residues:
                mutated_sequences.append(each)  # need to start with something in the list...
        else:
            new_mutated_list = []
            for residue in residues:
                for seq in mutated_sequences:
                    new_mutated_list.append(seq + residue)
            mutated_sequences = new_mutated_list

    # INDEX CONSERVED
    else:
        if i == 0:
            mutated_sequences.append(input_fasta_list[i])
        else:
            new_mutated_list = []
            for seq in mutated_sequences:
                new_mutated_list.append(seq + input_fasta_list[i])
            mutated_sequences = new_mutated_list

ref_file_name = output_folder_path + "/REFERENCE.txt"  # the output folder path will contain the name of the folder
os.makedirs(os.path.dirname(ref_file_name), exist_ok=True)  # this should create the directory correctly
with open(ref_file_name, 'w+') as reference_file:  # have the reference file open the whole time
    for i in range(len(mutated_sequences)):
        gen_file_name = output_folder_path + "/" + (sys.argv[1]).split(".")[0] + "_" + str(i + 1) + ".fasta"
        with open(gen_file_name, 'w+') as seq_file:
            seq_file.write(">" + (sys.argv[1]).split(".")[0] + "_" + str(i + 1) + "\n" + mutated_sequences[i])
            reference_file.write(">" + (sys.argv[1]).split(".")[0] + "_" + str(i + 1) + "\n" + mutated_sequences[i]+ "\n")





