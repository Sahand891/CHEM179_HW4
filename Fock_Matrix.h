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
    int A, L; // each AO contains information of which atom is came from!! Best to iterate through AOs!!!
    Atom origin_atom;
    contr_Gaussian cGTO;
};

std::vector<AO> atoms_to_AOs(const std::vector<Atom> &Atoms);

// Gives the number of alpha and beta electrons (p and q) in a molecule
int count_alpha_electrons(const std::vector<Atom> &Atoms);
int count_beta_electrons(const std::vector<Atom> &Atoms);
int count_total_electrons(const std::vector<Atom> &Atoms);

// Return matrix elements and matrix for gamma
double gamma_AB(Atom A, Atom B);
arma::mat gamma_matrix(std::vector<Atom> atoms);

// Compute overlap matrix from a vector of atoms
arma::mat S(const std::vector<Atom> &Atoms);

// Compute a core, 1-electron Hamiltonian; currently only works for two atoms
arma::mat h(const std::vector<Atom> &atoms);


// Function that takes in MO coefficient matrix for alpha, a number of alpha electrons, returns the rectangular occupied coefficient matrix
arma::mat C_occ_alpha(arma::mat c_alpha, int p);
// Calculates density matrix for alpha electrons from occupied coefficient matrix for alpha electrons
arma::mat P_alpha(arma::mat c_occ_alpha);

// Returns vector of total electron density on each atom
arma::vec P_AA(const std::vector<Atom> &atoms, arma::mat P_alpha, arma::mat P_beta);

// Fock matrix for alpha electrons
arma::mat F_alpha(const std::vector<Atom> &atoms, arma::mat P_alpha, arma::mat P_beta);


// Only works for two atoms
double f_uu_alpha(int u, int v, Atom A, Atom B, arma::mat c_alpha, arma::mat c_beta, int num_electrons);



#endif //HW4_FOCK_MATRIX_H
