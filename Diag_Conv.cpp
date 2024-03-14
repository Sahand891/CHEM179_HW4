//
// Created by Sahand Adibnia on 3/13/24.
//

#include "Diag_Conv.h"

arma::mat find_MO_coefs(arma::mat F) {

    arma::vec energies;
    arma::mat C;

    arma::eig_sym(energies, C, F); // can use eig_sym because F is guranteed to be symmetric

    return C;
}
