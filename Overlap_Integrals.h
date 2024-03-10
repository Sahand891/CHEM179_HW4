//
// Created by Sahand Adibnia on 3/10/24.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include <armadillo>

#ifndef HW4_OVERLAP_INTEGRALS_H
#define HW4_OVERLAP_INTEGRALS_H

struct prim_Gaussian {

    /* Every unnormalized primitive Gaussian (3D) has:
     * 1. A set of coordinates: X, Y, Z
     * 2. An exponent (alpha)
     * 3. Dimensional angular momentum (exponents for X, Y, Z coords)
     * 4. Total angular momentum (sum of dimensional ang momentums)
     */

    double X, Y, Z, alpha;
    int l=0, m=0, n=0;
    int L = l+m+n;

};

double normalizer(prim_Gaussian primG);

struct contr_Gaussian {
    std::vector<prim_Gaussian> primGs;
    std::vector<double> contr_coefs;
    std::vector<double> norm_constants;
};

struct STO3G {
    std::vector<double> exps;
    std::vector<double> con_coef_s;
    std::vector<double> con_coef_p;
};


struct Atom {
    // General atom data type

    double X,Y,Z;
    std::vector<contr_Gaussian> cGTOs;
    int A; // atomic number

    // Semi-emiprical parameters
    double IA_s;
    double IA_p;
    double beta;

    // valence atomic number
    int Z_ = std::max(1, A-2);
};


Atom construct_Atom(double X, double Y, double Z, std::string atom_symbol);

std::vector<Atom> read_atoms_from_file(const std::string& path);


// Overlap integral between two primitive Gaussians, based on HW2 implementation
double S_ab(prim_Gaussian primG1, prim_Gaussian primG2);

// Overlap integral between two contracted Gaussians, what I will input to overlap matrix S
double S_uv(contr_Gaussian contrG1, contr_Gaussian contrG2);



#endif //HW4_OVERLAP_INTEGRALS_H
