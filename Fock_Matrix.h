//
// Created by Sahand Adibnia on 3/10/24.
//


#ifndef HW4_FOCK_MATRIX_H
#define HW4_FOCK_MATRIX_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include <armadillo>
#include "Overlap_Integrals.h"


// Atomic orbital struct
struct AO {
    int A, L;
    contr_Gaussian cGTO;
};

std::vector<AO> atoms_to_AOs(const std::vector<Atom> &Atoms);

// Gives the number of electrons in a molecule
int count_electrons(const std::vector<Atom> &Atoms);

// Compute overlap matrix from a vector of atoms
arma::mat S(const std::vector<Atom> &Atoms);


// Function that takes in MO coefficient matrix, returns density matrix element
double p_uv_alpha(int u, int v, arma::mat c);

double p_uv_tot(int u, int v, arma::mat c_alpha, arma::mat c_beta);




double gamma_AB(Atom A, Atom B);
arma::mat gamma_matrix(std::vector<Atom> atoms);



#endif //HW4_FOCK_MATRIX_H
