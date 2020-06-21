import matplotlib.pyplot as plt
import numpy as np
import sys
import os


#takes in the text input files, parses them, and adds them to a graph
def plot(directory):
    print("Plot method called")
    file_names = []
    hydro = []
    SASA = []
    hydroSASA = []
    hbond_acc = []
    hbond_donor = []
    total_energy = []

    # handle the command line input --> will be a directory with a lot of .txt files in it
    for file in os.listdir(directory):
        if file.endswith(".txt"):
            try:
                fileName = directory + "/" + file
                with open(fileName, "r+") as open_text_file:
                    contents = parse_score_string(file, open_text_file.read())
                    file_names.append(file.replace(".txt", ""))
                    hydro.append(float(contents[0]))
                    SASA.append(float(contents[1]))
                    hydroSASA.append(float(contents[2]))
                    hbond_acc.append(float(contents[3]))
                    hbond_donor.append(float(contents[4]))
                    total_energy.append(float(contents[5]))
                    print(contents)
            except FileNotFoundError:
                print("File was not found: " + fileName)


    plt.style.use('ggplot')

    # output hydro plot
    hydro_plot = plt.figure(1)
    x = file_names
    x_pos = [i for i, _ in enumerate(x)]
    plt.bar(x_pos, hydro, color='red')
    plt.xticks(x_pos, x, rotation='vertical')
    rects = plt.patches
    plt.title("Hydro")

    # # output SASA plot
    SASA_plot = plt.figure(2)
    plt.bar(x_pos, SASA, color='orange')
    plt.xticks(x_pos, x, rotation='vertical')
    plt.title("SASA")

    hydoSASA_plot = plt.figure(3)
    plt.bar(x_pos, hydroSASA, color='green')
    plt.xticks(x_pos, x, rotation='vertical')
    plt.title("hydroSASA")

    hbond_acc_plot = plt.figure(4)
    plt.bar(x_pos, hbond_acc, color='blue')
    plt.xticks(x_pos, x, rotation='vertical')
    plt.title("hbond acceptor")

    hbond_donor_plot = plt.figure(5)
    plt.bar(x_pos, hbond_donor, color='purple')
    plt.xticks(x_pos, x, rotation='vertical')
    plt.title("hbond donor")

    total_energy_plot = plt.figure(6)
    plt.bar(x_pos, total_energy, color='black')
    plt.xticks(x_pos, x, rotation='vertical')
    plt.title("Total Energy")




    plt.show()


    input()

    print("ARRAYS")
    print(file_names)
    print(hydro)
    print(SASA)
    print(hydroSASA)
    print(hbond_acc)
    print(hbond_donor)
    print(total_energy)


def parse_score_string(filename, score_string):
    filename = filename.replace(".txt", "")
    score_array = score_string.split(" |")

    #print("File name: " + filename)
    #print(filename)
    array_to_return = []
    #array_to_return.append(str(filename.strip()))
    array_to_return.append(score_array[2].strip())  # hydro
    array_to_return.append(score_array[4].strip())  # SASA
    array_to_return.append(score_array[6].strip())  # hydroSASA
    array_to_return.append(score_array[8].strip())  # hbond-acc
    array_to_return.append(score_array[10].strip())  # hbond-donor
    array_to_return.append(score_array[12].strip())  # total-energy
    return array_to_return






    #print("Score string: " + score_string)





plot("OrigPDB_andStats")






