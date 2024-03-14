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
    arma::mat H2_gamma = gamma_matrix(H2_Atoms);
    //H2_gamma.print();

    arma::mat H2_S = S(H2_Atoms);
    //H2_S.print();

    arma::mat H2_H_core = h(H2_Atoms);
    //H2_H_core.print();

    arma::mat H2_C_alpha_1 = find_MO_coefs(H2_H_core);
    //H2_C_alpha_1.print();

    arma::mat H2_C_beta_1 = find_MO_coefs(H2_H_core);
    //H2_C_beta_1.print();

    int p = count_alpha_electrons(H2_Atoms);
    int q = count_beta_electrons(H2_Atoms);

    std::cout << "p = " << p << std::endl;
    std::cout << "q = " << q << std::endl;

    arma::mat H2_C_occ_alpha_1 = C_occ_alpha(H2_C_alpha_1, p);
    //H2_C_occ_alpha_1.print();
    arma::mat H2_C_occ_beta_1 = C_occ_alpha(H2_C_beta_1, p);
    //H2_C_occ_beta_1.print();

    arma::mat H2_P_alpha = P_alpha(H2_C_occ_alpha_1);
    //H2_P_alpha.print();
    arma::mat H2_P_beta = P_alpha(H2_C_occ_beta_1);
    //H2_P_beta.print();

    arma::mat H2_P_AA = P_AA(H2_Atoms, H2_P_alpha, H2_P_beta);
    //H2_P_AA.print();

    arma::mat H2_F_alpha = F_alpha(H2_Atoms, H2_P_alpha, H2_P_beta);
    //H2_F_alpha.print();

    arma::mat H2_F_beta = F_alpha(H2_Atoms, H2_P_beta, H2_P_alpha);
    //H2_F_beta.print();


    // Second diagonalization
    arma::mat H2_C_alpha_2 = find_MO_coefs(H2_F_alpha);
    //H2_C_alpha_2.print();

    arma::mat H2_C_beta_2 = find_MO_coefs(H2_F_beta);
    //H2_C_beta_2.print();

    // Compare new vs old density matrices
    arma::mat H2_C_occ_alpha_2 = C_occ_alpha(H2_C_alpha_2, p);
    arma::mat H2_P_alpha_new = P_alpha(H2_C_occ_alpha_1);
    //H2_P_alpha.print();
    arma::mat H2_C_occ_beta_2 = C_occ_alpha(H2_C_beta_2, p);
    arma::mat H2_P_beta_new = P_alpha(H2_C_occ_beta_2);
    //H2_P_beta.print();

    // Confirming that covergence is acheived!!!!
    (H2_P_alpha_new - H2_P_alpha).print();
    (H2_P_beta_new - H2_P_beta).print();


    // Now let's compute the energies



    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao/HW4/sample_input/HF.txt";
    std::vector<Atom> HF_Atoms = read_atoms_from_file(path);
    arma::mat HF_gamma = gamma_matrix(HF_Atoms);
    //HF_gamma.print();

    arma::mat HF_S = S(HF_Atoms);
    //HF_S.print();

    arma::mat HF_H_core = h(HF_Atoms);
    //HF_H_core.print();






    // So far, all calculations are correct for H2 and HF, but not quite matching Xiao's for HO ... not sure if that's me (HO is the only unrestricted case!) or if Xiao messed up in the HO calculation



    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao/HW4/sample_input/HO.txt";
    std::vector<Atom> HO_Atoms = read_atoms_from_file(path);
    arma::mat HO_gamma = gamma_matrix(HO_Atoms);
    //HO_gamma.print();

    arma::mat HO_S = S(HO_Atoms);
    //HO_S.print();

    arma::mat HO_H_core = h(HO_Atoms);
    //HO_H_core.print();


    return 0;
}
