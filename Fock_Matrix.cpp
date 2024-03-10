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

double p_uv_alpha(int u, int v, arma::mat c) {
    double sum = 0;
    for (int i = 0; i < c.size(); i++) {
        double val = c(u, i) * c(v, i);
        sum += val;
    }
    return sum;
}

double p_uv_tot(int u, int v, arma::mat &c_alpha, arma::mat &c_beta) {

    double val1 = p_uv_alpha(u, v, c_alpha);
    double val2 = p_uv_alpha(u, v, c_beta);

    return val1 + val2;
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

int count_electrons(const std::vector<Atom> &Atoms) {
    // Here, n is the total number of electrons
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
