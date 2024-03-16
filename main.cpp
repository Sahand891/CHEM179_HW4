//
// Created by Sahand Adibnia on 3/3/24.
//


/* HW 4 Notes
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */

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

    // H2 used as an example to demonstrate how each function works

    // Read in a vector of atoms (in a molecule) from an appropriately formatted text file
    std::string path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao/HW4/sample_input/H2.txt";
    std::vector<Atom> H2_Atoms = read_atoms_from_file(path);
//    arma::mat H2_gamma = gamma_matrix(H2_Atoms);
//    //H2_gamma.print();
//
//    arma::mat H2_S = S(H2_Atoms);
//    //H2_S.print();
//
//    arma::mat H2_H_core = h(H2_Atoms);
//    //H2_H_core.print();
//
//    arma::mat H2_C_alpha_1 = find_MO_coefs(H2_H_core);
//    //H2_C_alpha_1.print();
//
//    arma::mat H2_C_beta_1 = find_MO_coefs(H2_H_core);
//    //H2_C_beta_1.print();
//
//    int p = count_alpha_electrons(H2_Atoms);
//    int q = count_beta_electrons(H2_Atoms);
//
//    //std::cout << "p = " << p << std::endl;
//    //std::cout << "q = " << q << std::endl;
//
//    arma::mat H2_C_occ_alpha_1 = C_occ_alpha(H2_C_alpha_1, p);
//    //H2_C_occ_alpha_1.print();
//    arma::mat H2_C_occ_beta_1 = C_occ_alpha(H2_C_beta_1, p);
//    //H2_C_occ_beta_1.print();
//
//    arma::mat H2_P_alpha = P_alpha(H2_C_occ_alpha_1);
//    //H2_P_alpha.print();
//    arma::mat H2_P_beta = P_alpha(H2_C_occ_beta_1);
//    //H2_P_beta.print();
//
//    arma::mat H2_P_AA = P_AA(H2_Atoms, H2_P_alpha, H2_P_beta);
//    //H2_P_AA.print();
//
//    arma::mat H2_F_alpha = F_alpha(H2_Atoms, H2_P_alpha, H2_P_beta);
//    //H2_F_alpha.print();
//
//    arma::mat H2_F_beta = F_alpha(H2_Atoms, H2_P_beta, H2_P_alpha);
//    //H2_F_beta.print();
//
//
//    // Second diagonalization
//    arma::mat H2_C_alpha_2 = find_MO_coefs(H2_F_alpha);
//    //H2_C_alpha_2.print();
//
//    arma::mat H2_C_beta_2 = find_MO_coefs(H2_F_beta);
//    //H2_C_beta_2.print();
//
//    // Compare new vs old density matrices
//    arma::mat H2_C_occ_alpha_2 = C_occ_alpha(H2_C_alpha_2, p);
//    arma::mat H2_P_alpha_new = P_alpha(H2_C_occ_alpha_2);
//    //H2_P_alpha.print();
//    arma::mat H2_C_occ_beta_2 = C_occ_alpha(H2_C_beta_2, q);
//    arma::mat H2_P_beta_new = P_alpha(H2_C_occ_beta_2);
//    //H2_P_beta.print();
//
//    // Confirming that covergence is acheived (after just one iteration lol)!!!!
//    //(H2_P_alpha_new - H2_P_alpha).print();
//    //(H2_P_beta_new - H2_P_beta).print();
//
//
//
//    iteration_data H2_initial_it_data = start_CNDO2(H2_Atoms);
//    iteration_data H2_conv_it_data = converge_CNDO2(H2_initial_it_data);
//
//
//
//
//
//    // Now let's compute the energies



    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao/HW4/sample_input/HF.txt";
    std::vector<Atom> HF_Atoms = read_atoms_from_file(path);
//    arma::mat HF_gamma = gamma_matrix(HF_Atoms);
//    //HF_gamma.print();
//
//    arma::mat HF_S = S(HF_Atoms);
//    //HF_S.print();
//
//    arma::mat HF_H_core = h(HF_Atoms);
//    //HF_H_core.print();
//
//    //iteration_data HF_initial_it_data = start_CNDO2(HF_Atoms);
//    //iteration_data HF_conv_it_data = converge_CNDO2(HF_initial_it_data);






    // So far, all calculations are correct for H2 and HF, but not quite matching Xiao's for HO ... not sure if that's me (HO is the only unrestricted case!) or if Xiao messed up in the HO calculation



    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao/HW4/sample_input/HO.txt";
    std::vector<Atom> HO_Atoms = read_atoms_from_file(path);
    arma::mat HO_gamma = gamma_matrix(HO_Atoms);
    //HO_gamma.print();

    arma::mat HO_S = S(HO_Atoms);
    //HO_S.print();

    arma::mat HO_H_core = h(HO_Atoms);
    //HO_H_core.print();

    //std::cout << count_total_electrons(HO_Atoms) << std::endl;
    //std::cout << count_alpha_electrons(HO_Atoms) << std::endl;

    //iteration_data HO_initial_it_data = start_CNDO2(HO_Atoms);
    //iteration_data HF_conv_it_data = converge_CNDO2(HO_initial_it_data);

//    std::cout << nuc_repl_energy(H2_Atoms) << std::endl;
//    std::cout << nuc_repl_energy(HF_Atoms) << std::endl;
//    std::cout << nuc_repl_energy(HO_Atoms) << std::endl;

    iteration_data H2_initial_it_data = start_CNDO2(H2_Atoms);
    iteration_data H2_conv_it_data = converge_CNDO2(H2_initial_it_data);

    iteration_data HF_initial_it_data = start_CNDO2(HF_Atoms);
    iteration_data HF_conv_it_data = converge_CNDO2(HF_initial_it_data);

    iteration_data HO_initial_it_data = start_CNDO2(HO_Atoms);
    iteration_data HO_conv_it_data = converge_CNDO2(HO_initial_it_data);

    iteration_data A = HF_conv_it_data;

    arma::mat B = h(HF_Atoms);
    std::cout << electron_repl_energy(A.fock_alpha, A.fock_beta, A.P_alpha_new, A.P_beta_new, B) << std::endl;





    return 0;
}
