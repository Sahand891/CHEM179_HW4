//
// Created by Sahand Adibnia on 3/10/24.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include <armadillo>
#include "Fock_Matrix.h"

std::vector<AO> atoms_to_AOs(const std::vector<Atom> &Atoms) {

    std::vector<AO> AOs;

    for (auto &atom: Atoms) {
        int atomic_num = atom.A;
        for (auto &cGTO: atom.cGTOs) {
            // just need to take angular momentum of ONE of the primitive Gaussians in a cGTO
            int L = cGTO.primGs[0].L; // will tell me if I have an s or a p orbital
            AO AO_to_append = {atomic_num, L, cGTO};
            AOs.push_back(AO_to_append);
        }
    }

    return AOs;

}

// Diagonal matrix elements of one-electron Hamiltonian matrix; currently only works for 2 atoms
double h_uu(int u, Atom A, std::vector<Atom> atoms) {

    Atom B = atoms[1];

    double term1=0;
    if (A.A == 1) { // if we have an H atom
        term1 = A.IA_s;
    } else { // if have any other atom - C, N, O, F
        if (A.cGTOs[u].primGs[0].L == 0) { // s orbital
            term1 = A.IA_s;
        } else if (A.cGTOs[u].primGs[0].L == 1) { // p orbital
            term1 = A.IA_p;
        }
    }


    double term2 = (A.Z_ - 0.5)* gamma_AB(A,A);

    // This is the part that makes it work for only two-atom molecules
    double term3 = B.Z_* gamma_AB(A,B);

    return -term1-term2-term3;

}

// Off-diagonal matrix elements of 1-electron core hamiltonian
double h_uv(int u, int v, Atom A, std::vector<Atom> atoms) {

    Atom B = atoms[1];

    double S_ele = S_uv(A.cGTOs[u],B.cGTOs[v]);

    return 0.5*(A.beta+B.beta)*S_ele;

}


arma::mat h(std::vector<Atom> atoms) {

    Atom A = atoms[0];
    Atom B = atoms[1];

    // # of AOs = size of Hamiltonian matrix
    int matrix_size = atoms_to_AOs(atoms).size();

    // Initialize a matrix of the appropriate size, just filled with zeroes
    arma::mat h(matrix_size, matrix_size, arma::fill::zeros);


    // Input diagonal components first
    for (size_t i = 0; i < matrix_size; i++) {
        Atom focus_atom = atoms[i];
        h(i,i) = h_uu(i, focus_atom, atoms);
    }

    // Then add off-diagonals
    for (size_t i = 0; i < matrix_size; i++) {
        for (size_t j = 0; j < matrix_size; j++) {
            if (i != j) {
                Atom focus_atom = atoms[i];
                h(i, j) = h_uv(i,j,focus_atom,atoms);
            }
        }
    }

    return h;

}


// Matrix elements of alpha (or beta) density matrix only
double p_uv_alpha(int u, int v, arma::mat c, int num_electrons) {

    // num electrons is p (alpha) or q (beta)
    // assumes c is ordered such that lowest energy orbitals have lower indices

    double sum = 0;
    for (int i = 0; i < num_electrons; i++) {
        double val = c(u, i) * c(v, i);
        sum += val;
    }
    return sum;
}

// Full P_alpha matrix
arma::mat p_uv_alpha_matrix(arma::mat c, int num_electrons) {

    int matrix_size = c.size();

    // Initialize a matrix of the appropriate size, just filled with zeroes
    arma::mat final_mat(matrix_size, matrix_size, arma::fill::zeros);

    // Iterate through indices in matrix, computing appropriate overlap integral!
    for (size_t i = 0; i < matrix_size; i++) {
        for (size_t j = 0; j < matrix_size; j++) {
            final_mat(i, j) = p_uv_alpha(i,j,c,num_electrons);
        }
    }

    return final_mat;

}

// Matrix elements of total density matrix
double p_uv_tot(int u, int v, arma::mat &c_alpha, arma::mat &c_beta, int p, int q) {

    double val1 = p_uv_alpha(u, v, c_alpha, p);
    double val2 = p_uv_alpha(u, v, c_beta, q);

    return val1 + val2;
}

// total electron density on one atom—the sum of diagonal total matrix elements
arma::vec p_AA(const std::vector<Atom> &atoms, arma::mat P_alpha, arma::mat P_beta){

    // First get the total density matrix
    arma::mat P_tot = P_alpha + P_beta;

    // Initialize an arma::vec of the appropriate size
    arma::vec final_vec(atoms.size(), arma::fill::zeros);

    double sum=0;
    // Assumes order of "atoms" vector corresponds to density matrix AO order
    size_t i =0;
    while (i < atoms.size()) {
        if (atoms[i].A == 1) { // if you're dealing with a Hydrogen atom as the atom you're at
            final_vec(i) = P_tot(i,i);
            i += 1;
        } else { // if it's not a hydrogen atom, it's gonna be a 2nd row element atom (C, N, O, F)
            // Now we take the sum of the next FOUR elements
            final_vec(i) = P_tot(i,i) + P_tot(i+1,i+1) + P_tot(i+2,i+2) + P_tot(i+3,i+3);
            i += 4;
        }
    }

    return final_vec;
}






double zero_zero_integral(Atom A, Atom B, int k, int k_prime, int l, int l_prime) {

    double final_val = 0;

    // Only care about valence s orbital, so always take the FIRST cGTO
    double alpha_k_A = A.cGTOs[0].primGs[k].alpha;
    double alpha_k_prime_A = A.cGTOs[0].primGs[k_prime].alpha;
    double alpha_k_B = B.cGTOs[0].primGs[l].alpha;
    double alpha_k_prime_B = B.cGTOs[0].primGs[l_prime].alpha;

    double sigma_A = 1 / (alpha_k_A + alpha_k_prime_A);
    double U_A = pow((M_PI * sigma_A), 1.5);

    double sigma_B = 1 / (alpha_k_B + alpha_k_prime_B);
    double U_B = pow((M_PI * sigma_B), 1.5);

    // Calculating d = (R_A - R_B)**2
    arma::vec d_vec = {(A.X - B.X), (A.Y - B.Y), (A.Z - B.Z)};
    double d = pow(arma::norm(d_vec),2);

    double V_squared = 1 / (sigma_A + sigma_B);
    double T = V_squared * d;


    if (T == 0) {
        final_val = U_A * U_B * sqrt(2 * V_squared) * sqrt(2.0 / M_PI);
    } else {
        final_val = U_A * U_B * sqrt(1.0 / d) * erf(sqrt(T));
    }

    return 27.211*final_val; // convert gamma elements to eV from atomic units

}

double gamma_AB(Atom A, Atom B) {

    // Note: only looks at valence s orbitals of A and B
    // This means we're looking at the FIRST cGTO ALWAYS
    // but at each index i, take the ith primitive Gaussian's contraction coefficient and normalization constant!

    double sum = 0;

    int count = 0;
    for (int k = 0; k < 3; k++) {
        double d_k_A = A.cGTOs[0].contr_coefs[k] * A.cGTOs[0].norm_constants[k];
        for (int k_prime = 0; k_prime < 3; k_prime++) {
            double d_kprime_A = A.cGTOs[0].contr_coefs[k_prime] * A.cGTOs[0].norm_constants[k_prime];
            for (int l = 0; l < 3; l++) {
                double d_l_B = B.cGTOs[0].contr_coefs[l] * B.cGTOs[0].norm_constants[l];
                for (int l_prime = 0; l_prime < 3; l_prime++) {
                    double d_lprime_B = B.cGTOs[0].contr_coefs[l_prime] * B.cGTOs[0].norm_constants[l_prime];

                    double integral = zero_zero_integral(A, B, k, k_prime, l, l_prime);

                    sum += d_k_A * d_kprime_A * d_l_B * d_lprime_B * integral;
                    count += 1;

                }
            }
        }
    }

    //std::cout << count << std::endl;

    return sum;
}

arma::mat gamma_matrix(std::vector<Atom> atoms) {

    int matrix_size = atoms.size();
    arma::mat final_mat(matrix_size, matrix_size, arma::fill::zeros);

    // Iterate through indices in matrix, computing appropriate value of gamma!
    for (size_t i = 0; i < matrix_size; i++) {
        for (size_t j = 0; j < matrix_size; j++) {
            double mat_ij = gamma_AB(atoms[i], atoms[j]);
            final_mat(i, j) = mat_ij;
        }
    }

    return final_mat;

}


arma::mat S(const std::vector<Atom> &Atoms) {

    std::vector<AO> AOs = atoms_to_AOs(Atoms);
    int matrix_size = AOs.size();

    // Initialize a matrix of the appropriate size, just filled with zeroes
    arma::mat S(matrix_size, matrix_size, arma::fill::zeros);

    // Iterate through indices in matrix, computing appropriate overlap integral!
    for (size_t i = 0; i < matrix_size; i++) {
        //std::cout<< i<< std::endl;
        for (size_t j = 0; j < matrix_size; j++) {
            double Sij = S_uv(AOs[i].cGTO, AOs[j].cGTO);
            S(i, j) = Sij;
            //std::cout<< i <<j<< std::endl;
            //S.print();
        }
    }

    return S;
}

int count_total_electrons(const std::vector<Atom> &Atoms) {

    // Here, n is the total number of electrons provided by atomic orbitals of the atoms
    int n = 0;
    for (auto &atom: Atoms) {
        if (atom.A == 1) { // H atom
            n += 1;
        } else { // C, N, O, or F atom
            n += atom.Z_;
        }
    }
    return n;
}

int count_alpha_electrons(const std::vector<Atom> &Atoms) {
    int n = count_total_electrons(Atoms);
    return ceil(n/2);
}

int count_beta_electrons(const std::vector<Atom> &Atoms) {
    int n = count_total_electrons(Atoms);
    return floor(n/2);
}


// Current implementation only works for two-atom molecules
double f_uu_alpha(int u, Atom A, Atom B, arma::mat P_alpha, arma::mat P_beta, int num_electrons) {

    // u is the index of the matrix = the index of the atomic orbital we're looking at

    double term1=0;
    if (A.cGTOs[u].primGs[0].L == 0) { // s orbital
        term1 = A.IA_s;
    } else if (A.cGTOs[u].primGs[0].L == 1) { // p orbital
        term1 = A.IA_p;
    }

    // Compile atoms A and B into a standard vector for formatting compatibility
    std::vector<Atom> atoms = {A,B};

    double term2 = ((p_AA(atoms, P_alpha, P_beta)[u] - A.Z_) - (P_alpha(u,u) - 0.5))*gamma_AB(A,A);

    // This is the part that makes it work for only two-atom molecules
    double term3 = (p_AA(atoms, P_alpha, P_beta)[u] - B.Z_)*gamma_AB(A,B);

    return -term1+term2+term3;

}
