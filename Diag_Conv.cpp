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

// This function starts CNDO2 density matrix optimization with an initial guess of the core Hamiltonian matrix, and also prints out the 0th iteration
iteration_data start_CNDO2(const std::vector<Atom> &atoms) {

    arma::mat gamma = gamma_matrix(atoms);

    arma::mat overlap = S(atoms);

    arma::mat H_core = h(atoms);

    // Compute appropriate matrices to input into iteration_data object

    int p = count_alpha_electrons(atoms);
    int q = count_beta_electrons(atoms);

    // C_alpha = C_beta = C for the zeroth iteration since we're using core Hamiltonian to start

    // 0th iteration: Diagonalizing H_core to get the MO coefficients
    arma::mat C = find_MO_coefs(H_core);

    // Calculating "new" density matrices for the 0th iteration (i.e., starting density matrices)
    arma::mat C_occ_a = C_occ_alpha(C, p);
    arma::mat C_occ_b = C_occ_alpha(C, q);
    arma::mat P_a = P_alpha(C_occ_a);
    arma::mat P_b = P_alpha(C_occ_b);

    // For the first iteration, old matrix just set to be a 1x1 matrix of 0's, and specifically set converged bool to be false (SUPER IMPORTANT!)
    // We won't even use P_alpha_old and P_beta_old in the 0th iteration anyway

    iteration_data it_data = {atoms,
                              arma::mat(1,1,arma::fill::zeros),
                              arma::mat(1,1,arma::fill::zeros),
                              H_core,
                              H_core,
                              C,
                              C,
                              P_a,
                              P_b,
                              1,
                              false};


    // Now print out everything nicely so we can see it!!
    std::cout << "gamma" << std::endl;
    gamma.print();
    std::cout << "Overlap" << std::endl;
    overlap.print();
    std::cout << "p = " << p << " q = " << q << std::endl;
    std::cout << "H_core" << std::endl;
    H_core.print();
    // Printing out zeroth iteration
    std::cout << "Iteration: 0" << std::endl;
    std::cout << "Fa" << std::endl;
    H_core.print();
    std::cout << "Fb" << std::endl;
    H_core.print();
    std::cout << "after solving the eigenvalue equation: " << std::endl;
    std::cout << "Ca" << std::endl;
    C.print();
    std::cout << "Cb" << std::endl;
    C.print();
    std::cout << "p = " << it_data.p << " q = " << it_data.q << std::endl;
    std::cout << "Pa_new" << std::endl;
    P_a.print();
    std::cout << "Pb_new" << std::endl;
    P_b.print();
    std::cout << "P_total" << std::endl;
    P_AA(atoms, P_a, P_b).print();

    return it_data;

}


iteration_data converge_CNDO2(const iteration_data &it_data) {

    // Base case
    if (it_data.converged) {
        return it_data;
    }

    // Extracting constant information for each iteration
    std::vector<Atom> atoms = it_data.atoms;
    int p = it_data.p;
    int q = it_data.q;


    // Old density matrices for this iteration at the new ones from the previous iteration
    arma::mat P_alpha_old = it_data.P_alpha_new;
    arma::mat P_beta_old = it_data.P_beta_new;


    // Make the new fock matrix based on the density matrices from the previous iteration
    arma::mat Fock_alpha = F_alpha(atoms, P_alpha_old, P_beta_old);
    arma::mat Fock_beta = F_alpha(atoms, P_beta_old, P_alpha_old);


    // Diagonalize the fock matrices and extra MO coefficients
    arma::mat C_alpha_new = find_MO_coefs(Fock_alpha);
    arma::mat C_beta_new = find_MO_coefs(Fock_beta);


    // Calculate the new density matrices based on the new MO coefficient matrices
    // Compare new vs old density matrices
    arma::mat C_occ_alpha_new = C_occ_alpha(C_alpha_new, p);
    arma::mat P_alpha_new = P_alpha(C_occ_alpha_new);

    arma::mat C_occ_beta_new = C_occ_alpha(C_beta_new, q);
    arma::mat P_beta_new = P_alpha(C_occ_beta_new);


    // Compile all the appropriate info into a new iteration_data object
    iteration_data new_it_data = {atoms,
                                  P_alpha_old,
                                  P_beta_old,
                                  Fock_alpha,
                                  Fock_beta,
                                  C_alpha_new,
                                  C_beta_new,
                                  P_alpha_new,
                                  P_beta_new,
                                  it_data.iteration_count+1};

    std::cout << "Iteration: " << it_data.iteration_count << std::endl;
    std::cout << "Fa" << std::endl;
    Fock_alpha.print();
    std::cout << "Fb" << std::endl;
    Fock_beta.print();
    std::cout << "after solving the eigenvalue equation: " << std::endl;
    std::cout << "Ca" << std::endl;
    C_alpha_new.print();
    std::cout << "Cb" << std::endl;
    C_beta_new.print();
    std::cout << "p = " << p << " q = " << q<< std::endl;
    std::cout << "Pa_new" << std::endl;
    P_alpha_new.print();
    std::cout << "Pb_new" << std::endl;
    P_beta_new.print();
    std::cout << "P_total" << std::endl;
    new_it_data.P_total_new.print();

    return converge_CNDO2(new_it_data);

}

// Calculate energy eigenvalues of a converged Fock matrix
arma::vec get_energy_eigs(arma::mat F) {
    arma::vec energies;
    arma::mat C;
    arma::eig_sym(energies, C, F); // can use eig_sym because F is guaranteed to be symmetric

    return energies;
}


// Calculate distance between two atoms
double atom_distance(const Atom &A, const Atom &B) {
    return sqrt(pow((A.X - B.X),2) + pow((A.Y - B.Y),2) + pow((A.Z - B.Z),2));
}

// Calculate nuclear repulsion energy for set of atoms
double nuc_repl_energy(const std::vector<Atom> &atoms) {
    double sum=0;
    for (auto& A : atoms) {
        for (auto& B : atoms) {
            if (A.X == B.X and A.Y == B.Y and A.Z == B.Z) { // check if they're the same atom
                sum += 0;
            } else {
                double R_AB = atom_distance(A,B);
                sum += A.Z_*B.Z_ / R_AB;
            }
        }
    }
    return sum * 27.211 / 2; // multiply by eV/a.u. ratio, divide by 2 to eliminate double counting
}


// Calculate electron energy from converged density and Fock matrices
double electron_repl_energy(const arma::mat &F_alpha, const arma::mat &F_beta, const arma::mat &P_alpha, const arma::mat &P_beta, const arma::mat &H_core) {

    int matrix_length = sqrt(F_alpha.size());

    // Iterating through indices of the matrices
    double alpha_term=0, beta_term=0;
    for (int u=0; u < matrix_length; u++) {
        for (int v=0; v < matrix_length; v++) {
            alpha_term += P_alpha(u,v)*(H_core(u,v)+F_alpha(u,v));
            beta_term += P_beta(u,v)*(H_core(u,v)+F_beta(u,v));
        }
    }

    return (alpha_term+beta_term) / 2;
}


// Calculate total energy of a converged system (sum of nuclear repulsion and electron energies)
double total_energy(double nuc_repl_e, double elec_repl_e) {
    return nuc_repl_e+elec_repl_e;
}