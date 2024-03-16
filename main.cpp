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

    // Read in a vector of atoms (in a molecule) from an appropriately formatted text file
    std::string path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao_HW4/sample_input/H2.txt";
    std::vector<Atom> H2_Atoms = read_atoms_from_file(path);

    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao_HW4/sample_input/HF.txt";
    std::vector<Atom> HF_Atoms = read_atoms_from_file(path);

    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao_HW4/sample_input/HO.txt";
    std::vector<Atom> HO_Atoms = read_atoms_from_file(path);


    // Now performing the CNDO/2 energy optimization
    iteration_data H2_initial_it_data = start_CNDO2(H2_Atoms);
    std::vector<iteration_data> H2_start_vec = {H2_initial_it_data};
    std::vector<iteration_data> H2_conv_it_datas = converge_CNDO2(H2_initial_it_data, H2_start_vec);

    iteration_data HF_initial_it_data = start_CNDO2(HF_Atoms);
    std::vector<iteration_data> HF_start_vec = {HF_initial_it_data};
    std::vector<iteration_data> HF_conv_it_datas = converge_CNDO2(HF_initial_it_data, HF_start_vec);

    iteration_data HO_initial_it_data = start_CNDO2(HO_Atoms);
    std::vector<iteration_data> HO_start_vec = {HO_initial_it_data};
    std::vector<iteration_data> HO_conv_it_datas = converge_CNDO2(HO_initial_it_data, HO_start_vec);


    // And now writing them to some nice files!
    write_CNDO2_to_file("/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/my_outputs/H2.txt", H2_initial_it_data, H2_conv_it_datas);
    write_CNDO2_to_file("/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/my_outputs/HF.txt", HF_initial_it_data, HF_conv_it_datas);
    write_CNDO2_to_file("/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/my_outputs/HO.txt", HO_initial_it_data, HO_conv_it_datas);



    // For fun - my own inputs!

    // Raw H atom
    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao_HW4/sample_input/H.txt";
    std::vector<Atom> H_Atom = read_atoms_from_file(path);
    iteration_data H_initial_it_data = start_CNDO2(H_Atom);
    std::vector<iteration_data> H_start_vec = {H_initial_it_data};
    std::vector<iteration_data> H_conv_it_datas = converge_CNDO2(H_initial_it_data, H_start_vec, false);

    write_CNDO2_to_file("/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/my_outputs/H.txt", H_initial_it_data, H_conv_it_datas);

    // Raw O atom
    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao_HW4/sample_input/O.txt";
    std::vector<Atom> O_Atom = read_atoms_from_file(path);
    iteration_data O_initial_it_data = start_CNDO2(O_Atom, false, false, 4, 2);
    std::vector<iteration_data> O_start_vec = {O_initial_it_data};
    std::vector<iteration_data> O_conv_it_datas = converge_CNDO2(O_initial_it_data, O_start_vec);

    write_CNDO2_to_file("/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/my_outputs/O.txt", O_initial_it_data, O_conv_it_datas);

    // Raw F atom
    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao_HW4/sample_input/F.txt";
    std::vector<Atom> F_Atom = read_atoms_from_file(path);
    iteration_data F_initial_it_data = start_CNDO2(F_Atom, false, false, 4, 3);
    std::vector<iteration_data> F_start_vec = {HF_initial_it_data};
    std::vector<iteration_data> F_conv_it_datas = converge_CNDO2(F_initial_it_data, F_start_vec);

    write_CNDO2_to_file("/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/my_outputs/F.txt", F_initial_it_data, F_conv_it_datas);


    // Printing out energies of H2, HO, HF relative to two individual H, O, F atoms!!

    double H_energy = compute_total_energy(H_conv_it_datas.back());
    double H2_energy = compute_total_energy(H2_conv_it_datas.back());

    double O_energy = compute_total_energy(O_conv_it_datas.back());
    double F_energy = compute_total_energy(F_conv_it_datas.back());

    double HO_energy = compute_total_energy(HO_conv_it_datas.back());
    double HF_energy = compute_total_energy(HF_conv_it_datas.back());

    std::cout << "Energy of H2 molecule calculated by CNDO/2 theory: " << H2_energy - 2*H_energy << " eV" << std::endl;
    std::cout << "Energy of HF molecule calculated by CNDO/2 theory: " << HF_energy - H_energy - F_energy << " eV" << std::endl;
    std::cout << "Energy of HO molecule calculated by CNDO/2 theory: " << HO_energy - H_energy - O_energy << " eV" << std::endl;


    // Now let's try a nitrogen atom with p=4, q=1!!

    // Raw N atom first (p=4, q=1)
    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao_HW4/sample_input/N.txt";
    std::vector<Atom> N_Atom = read_atoms_from_file(path);
    iteration_data N_initial_it_data = start_CNDO2(N_Atom, false, false, 4, 1);
    std::vector<iteration_data> N_start_vec = {N_initial_it_data};
    std::vector<iteration_data> N_conv_it_datas = converge_CNDO2(N_initial_it_data, N_start_vec);

    write_CNDO2_to_file("/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/my_outputs/N.txt", N_initial_it_data, N_conv_it_datas);

    // Now N2
    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao_HW4/sample_input/N2.txt";
    std::vector<Atom> N2_Atoms = read_atoms_from_file(path);
    iteration_data N2_initial_it_data = start_CNDO2(N2_Atoms);
    std::vector<iteration_data> N2_start_vec = {N2_initial_it_data};
    std::vector<iteration_data> N2_conv_it_datas = converge_CNDO2(N2_initial_it_data, N2_start_vec);

    write_CNDO2_to_file("/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/my_outputs/N2.txt", N2_initial_it_data, N2_conv_it_datas);

    double N_energy = compute_total_energy(N_conv_it_datas.back());
    double N2_energy = compute_total_energy(N2_conv_it_datas.back());

    std::cout << "Energy of N2 molecule calculated by CNDO/2 theory: " << N2_energy - 2*N_energy << " eV" << std::endl;





    return 0;
}
