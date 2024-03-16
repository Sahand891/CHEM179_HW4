//
// Created by Sahand Adibnia on 3/13/24.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include <armadillo>
#include "Overlap_Integrals.h"
#include "Fock_Matrix.h"

#ifndef HW4_DIAGONALIZATION_AND_CONVERGENCE_H
#define HW4_DIAGONALIZATION_AND_CONVERGENCE_H

struct iteration_data {

    std::vector<Atom> atoms;

    arma::mat P_alpha_old;
    arma::mat P_beta_old;

    arma::mat fock_alpha;
    arma::mat fock_beta;

    arma::mat C_alpha_new;
    arma::mat C_beta_new;

    arma::mat P_alpha_new;
    arma::mat P_beta_new;

    int iteration_count;

    // check for convergence just based on alpha density matrix
    bool converged = arma::approx_equal(P_alpha_new,P_alpha_old,"absdiff",1e-6);

    // Other information in case we need it
    arma::vec P_total_new = P_AA(atoms, P_alpha_new, P_beta_new);

    arma::mat gamma = gamma_matrix(atoms);

    int p = count_alpha_electrons(atoms);
    int q = count_beta_electrons(atoms);


};


arma::mat find_MO_coefs(arma::mat F);
arma::vec get_energy_eigs(arma::mat F);


iteration_data start_CNDO2(const std::vector<Atom> &atoms);
iteration_data converge_CNDO2(const iteration_data &it_data);



double nuc_repl_energy(const std::vector<Atom> &atoms);
double electron_repl_energy(const arma::mat &F_alpha, const arma::mat &F_beta, const arma::mat &P_alpha, const arma::mat &P_beta, const arma::mat &H_core);



#endif //HW4_DIAGONALIZATION_AND_CONVERGENCE_H
