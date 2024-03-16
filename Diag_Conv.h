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

    int p;
    int q;

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
    arma::mat overlap = S(atoms);
    arma::mat H_core = h(atoms);



};


arma::mat find_MO_coefs(arma::mat F);
arma::vec get_energy_eigs(arma::mat F);


iteration_data start_CNDO2(const std::vector<Atom> &atoms, bool do_print=false, bool auto_electrons=true, int p_input=0, int q_input=0);
// Stores a vector of iteration data, which can then be put in an output file nicely!!
std::vector<iteration_data> converge_CNDO2(const iteration_data &it_data, std::vector<iteration_data> &it_data_vec, bool do_print=false);


double nuc_repl_energy(const std::vector<Atom> &atoms);
double electron_energy(const arma::mat &F_alpha, const arma::mat &F_beta, const arma::mat &P_alpha, const arma::mat &P_beta, const arma::mat &H_core);
double compute_total_energy(const iteration_data &conv_it_data);


void write_CNDO2_to_file(std::string output_location, iteration_data &starting_it_data, std::vector<iteration_data> &it_data_vec);



#endif //HW4_DIAGONALIZATION_AND_CONVERGENCE_H
