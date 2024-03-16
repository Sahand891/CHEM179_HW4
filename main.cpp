//
// Created by Sahand Adibnia on 3/3/24.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include <armadillo>
#include "Overlap_Integrals.h"
#include "Fock_Matrix.h"
#include "Diag_Conv.h"

int main() {

    // Empty vector for iteration data
    std::vector<iteration_data> empty_vect;

    // Read in a vector of atoms (in a molecule) from an appropriately formatted text file
    std::string path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao/HW4/sample_input/H2.txt";
    std::vector<Atom> H2_Atoms = read_atoms_from_file(path);

    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao/HW4/sample_input/HF.txt";
    std::vector<Atom> HF_Atoms = read_atoms_from_file(path);

    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao/HW4/sample_input/HO.txt";
    std::vector<Atom> HO_Atoms = read_atoms_from_file(path);


    // Now performing the CNDO/2 energy optimization
    iteration_data H2_initial_it_data = start_CNDO2(H2_Atoms);
    std::vector<iteration_data> H2_conv_it_datas = converge_CNDO2(H2_initial_it_data, empty_vect);

    iteration_data HF_initial_it_data = start_CNDO2(HF_Atoms);
    std::vector<iteration_data> HF_conv_it_datas = converge_CNDO2(HF_initial_it_data, empty_vect);

    iteration_data HO_initial_it_data = start_CNDO2(HO_Atoms);
    std::vector<iteration_data> HO_conv_it_datas = converge_CNDO2(HO_initial_it_data, empty_vect);


    // And now writing them to some nice files!
    write_CNDO2_to_file("/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/my_outputs/H2.txt", H2_initial_it_data, H2_conv_it_datas);
    write_CNDO2_to_file("/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/my_outputs/HF.txt", HF_initial_it_data, HF_conv_it_datas);
    write_CNDO2_to_file("/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/my_outputs/HO.txt", HO_initial_it_data, HO_conv_it_datas);


    return 0;
}
