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
    H2_H_core.print();

//    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao/HW4/sample_input/HF.txt";
//    std::vector<Atom> HF_Atoms = read_atoms_from_file(path);
//    arma::mat HF_gamma = gamma_matrix(HF_Atoms);
//    HF_gamma.print();
//
//    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW4/Xiao/HW4/sample_input/HO.txt";
//    std::vector<Atom> HO_Atoms = read_atoms_from_file(path);
//    arma::mat HO_gamma = gamma_matrix(HO_Atoms);
//    HO_gamma.print();


    return 0;
}
