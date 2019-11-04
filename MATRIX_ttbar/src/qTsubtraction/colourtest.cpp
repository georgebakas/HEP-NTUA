
#include <cstddef>
#include <iostream>
#include <cassert>
#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include "ColorFull/Core/Trace_basis.h"
#include "ColorFull/Core/Orthogonal_basis.h"
#include "ColorFull/Core/Tree_level_gluon_basis.h"


extern "C" {
  void ol_setparameter_int(const char* param, int val);
  int ol_register_process(const char* process, int amptype);
  int ol_n_external(int id);
  void ol_tree_colbasis_dim(int id, int* ncolb, int* colelemsz, int* nheltot);
  void ol_tree_colbasis(int id, int* basis, int* needed);
  void ol_tree_colourflow(int id, int* flowbasis);
  void ol_phase_space_point(int id, double energy, double* pp);
  void ol_evaluate_tree(int id, double* pp, double* m2);
  void ol_evaluate_tree_colvect(int id, double* pp, double* amp, int* nhel);
  void ol_evaluate_tree_colvect2(int id, double* pp, double* m2arr );
  void ol_start();
  void ol_finish();
}

int factorial(int n) {
  return (n == 1 || n == 0) ? 1 : factorial(n-1)*n;
}

template<typename T> class Tensor3;

// A Dynamic rank 2 array which is stored in a (contiguous) rank 1 array.
// Can be used to exchange data with Fortran code (automatic cast to T*):
// type, dimension(y,x) ^= Matrix<type>(x,y)
// The order of the indices is reversed between C++ and Fortran
// as with static C++ arrays.
// Note that there is no copy constructor / assignment operator or whatsoever.
template<typename T>
class Matrix {
  friend class Tensor3<T>;
public:
  Matrix(std::size_t xdim, std::size_t ydim):
    allocated(true), xdim_(xdim), ydim_(ydim), arr(new T[xdim*ydim]) {}
  ~Matrix() {
    if (allocated) {
      delete[] arr;
    }
  }
  std::size_t xdim() const {
    return xdim_;
  }
  std::size_t ydim() const {
    return ydim_;
  }
  T* operator[](std::size_t x) const {
    assert(x < xdim_);
    return &arr[x*ydim_];
  }
  operator T*() const {
    return &arr[0];
  }
private:
  Matrix(): allocated(false) {}
  bool allocated;
  std::size_t xdim_, ydim_;
  T* arr;
  void point_to(T* ar, std::size_t xdim, std::size_t ydim) {
    xdim_ = xdim;
    ydim_ = ydim;
    arr = ar;
  }
};


// A Dynamic rank 3 array which is stored in a (contiguous) rank 1 array.
// Can be used to exchange data with Fortran code (automatic cast to T*):
// type, dimension(z,x,y) ^= Tensor3<type>(x,y,z)
// The order of the indices is reversed between C++ and Fortran
// as with static C++ arrays.
// Note that there is no copy constructor / assignment operator or whatsoever.
template<typename T>
class Tensor3
{
public:
  Tensor3(std::size_t xdim, std::size_t ydim, std::size_t zdim)
    : xdim_(xdim), ydim_(ydim), zdim_(zdim),
      arr(new T[xdim*ydim*zdim]), matrices(new Matrix<T>[xdim]) {
    for (std::size_t n = 0; n < xdim; ++n) {
      matrices[n].point_to(&arr[n*ydim*zdim], ydim, zdim);
    }
  }
  ~Tensor3() {
    delete[] arr;
    delete[] matrices;
  }
  std::size_t xdim() const {
    return xdim_;
  }
  std::size_t ydim() const {
    return ydim_;
  }
  std::size_t zdim() const {
    return zdim_;
  }
  Matrix<T> operator[](std::size_t x) const {
    assert(x < xdim_);
    return matrices[x];
  }
  operator T*() const {
    return &arr[0];
  }
private:
  std::size_t xdim_, ydim_, zdim_;
  T* arr;
  Matrix<T>* matrices;
};


// Calculate the product of the spin and colour average factors
// of the (two) initial state particles and the symmetry factor
// of the outgoing particles.
int avg_factor(std::vector<int> particles) {
  int avgfac = 1;
  // initial state gluon: colour average factor 8
  if (particles[0] == 21) avgfac *= 8;
  if (particles[1] == 21) avgfac *= 8;
  // initial state quark: colour average factor 3
  if (std::abs(particles[0]) >= 1 and std::abs(particles[0]) <= 6) avgfac *= 3;
  if (std::abs(particles[1]) >= 1 and std::abs(particles[1]) <= 6) avgfac *= 3;
  // spin average
  // fermion, photon, gluon: 2
  if (std::abs(particles[0]) <= 22) avgfac *= 2;
  if (std::abs(particles[1]) <= 22) avgfac *= 2;
  // Z,W+,W-: 3
  if (particles[0] == 23 or std::abs(particles[0]) == 24) avgfac *= 3;
  if (particles[1] == 23 or std::abs(particles[1]) == 24) avgfac *= 3;
  // symmetry factor of outgoing particles
  std::map<int,int> outmap;
  for (std::vector<int>::size_type p = 2; p < particles.size(); ++p) {
    outmap[particles[p]] += 1;
  }
  for (std::map<int,int>::iterator it = outmap.begin();
       it != outmap.end(); ++it) {
    avgfac *= factorial(it->second);
  }
  return avgfac;
}


// Import a trace basis as returned by OpenLoops as Matrix() into ColorFull,
// creating a Col_basis object. The scalar product matrix and the leading nc
// scalar product matrix are calculated and stored inside the Col_basis
// object. Return the Col_basis object.
ColorFull::Col_basis convert_colbasis(const Matrix<int>& olbasis,
                                      const std::vector<int>& particles) {
  // import the trace basis with ColorFull
  ColorFull::Col_basis colbasis;
  int ncolb, colelemsz, col;
  std::vector<int> qline;
  ncolb = olbasis.xdim();
  colelemsz = olbasis.ydim();
  for (int r = 0; r < ncolb; ++r) {
    ColorFull::Col_str basisel;
    for (int c = 0; c < colelemsz; ++c) {
      col = olbasis[r][c];
      if ((col == 0 or (col != 0 and c == colelemsz-1)) and !qline.empty()) {
        // end of qline
        ColorFull::Quark_line ql;
        if (particles[qline.back()-1] == 21) {
          // is a trace
          // Tr(a,b,c) -> (a,b,c)
          for (std::vector<int>::iterator it = qline.begin();
               it != qline.end(); ++it) {
            ql.append(*it);
          }
          ql.open = false;
        }
        else {
          // is a chain
          // T(a,b,c)_ij -> {i,a,b,c,j}
          ql.append(qline.end()[-2]);
          for (std::vector<int>::iterator it = qline.begin();
               it != qline.end()-2; ++it) {
            ql.append(*it);
          }
          ql.append(qline.back());
          ql.open = true;
        }
        basisel = basisel * ql;
        qline.clear();
      }
      else {
        qline.push_back(col);
      }
    }
    colbasis.append(basisel);
  }
  // calculate the colour matrix
  colbasis.scalar_product_matrix();
  colbasis.leading_scalar_product_matrix();
  return colbasis;
}


int colourtest()
{
  int id, ncolb, colelemsz, nheltot, nex;
  std::vector<int> particles;

  ol_setparameter_int("no_splash", 1);
  ol_setparameter_int("verbose", 1);
  ol_setparameter_int("order_ew", 0);


  particles.push_back(21);
  particles.push_back(21);
  particles.push_back(6);
  particles.push_back(-6);
  // Build a string "a b -> c d ..." (using a stream) and register the process.
  std::stringstream procstrstream;
  procstrstream << particles[0] << " " << particles[1] << " ->";
  for (std::vector<int>::iterator it = particles.begin() + 2;
        it != particles.end(); ++it) {
    procstrstream << " " << *it;
  }
  id = ol_register_process(procstrstream.str().c_str(), 1);

  ol_start();

  // get the trace basis from OpenLoops
  nex = ol_n_external(id);
  ol_tree_colbasis_dim(id, &ncolb, &colelemsz, &nheltot);
  Matrix<int> basis = Matrix<int>(ncolb,colelemsz);
  Matrix<int> needed = Matrix<int>(ncolb,ncolb);
  ol_tree_colbasis(id, basis, needed);

  // import the colour basis with ColorFull and calculate the colour matrix
  ColorFull::Col_basis colbasis = convert_colbasis(basis, particles);

  // calculate the amplitude for a random phase space point
  int nhel;
  double energy = 500.;
  Matrix<double> pp = Matrix<double>(nex,4);
  Matrix<double> amp = Matrix<double>(nheltot,2*ncolb);
  double m2;
  ol_phase_space_point(id, energy, pp);
  ol_evaluate_tree_colvect(id, pp, amp, &nhel);
  // The array amp[h] is the amplitude for helicity configuration h
  // with the real and imaginary parts both as double
  // {re1,im1,re2,im2,...} for each colour basis element.
  // Only the first nhel helicity configurations are non-zero.
  // For comparison: calculate the squared matrix element with OpenLoops
  ol_evaluate_tree(id, pp, &m2);
  // Calculate the matrix element from the amplitude using the colour matrix:
  // me = sum_h amp*.colmatrix.amp (h numbers the helicity configurations)
  int avgfac;
  double me, diagme;
  std::vector<double> diagmes;
  me = 0;
  avgfac = avg_factor(particles);
  for (int h = 0; h < nhel; ++h) {
    for (int r = 0; r < ncolb; ++r) {
      for (int c = 0; c < ncolb; ++c) {
        me += colbasis.d_spm[r][c]
            * (amp[h][2*r]*amp[h][2*c] + amp[h][2*r+1]*amp[h][2*c+1]);
      }
    }
  }
  // Apply average and symmetry factors.
  me /= avgfac;
  // Full colour diagonal elements of the colour interference matrix
  // for comparison with leading colour.
  for (int r = 0; r < ncolb; ++r) {
    diagme = 0;
    for (int h = 0; h < nhel; ++h) {
      diagme += (amp[h][2*r]*amp[h][2*r] + amp[h][2*r+1]*amp[h][2*r+1])
              * colbasis.d_spm[r][r];
    }
    diagmes.push_back(diagme/avgfac);
  }

  // Calculate the helicity summed diagonal elements of the colour
  // interference matrix in leading colour with OpenLoops. Colour and
  // average factors are not included. The colour factor is the same for all
  // elements, namely (TF*Nc)^ng*(Nc)^nqq) where ng=(number of external
  // gluons) and nqq=(number of quark lines).
  std::vector<double> m2arr;
  m2arr.reserve(ncolb);
  ol_evaluate_tree_colvect2(id, pp, &m2arr[0]);

  // colour basis in colour flow notation:
  // The two integers in flowbasis[ncolb][nex] give the incoming and outgoing
  // fundamental colour through external leg 'nex' in basis element
  // number 'ncolb'.
  Tensor3<int> flowbasis(ncolb,nex,2);
  ol_tree_colourflow(id, flowbasis);

  ol_finish();

  // screen output
  std::cout.precision(15);
  std::cout << std::endl << "Colour basis and diagonal elements of the colour"
            " interference matrix:" << std::endl;
  std::cout << std::endl << "         full colour |      leading colour"
            " |       ratio lc/full |" << std::setw(4*nex) << "colour flow basis"
            << " | trace colour basis" << std::endl;
  for (int n = 0; n < ncolb; ++n) {
    std::cout << std::setw(21) << diagmes[n] << std::setw(22)
              << colbasis.leading_d_spm[n][n]*m2arr[n]/avgfac << std::setw(22)
              << colbasis.leading_d_spm[n][n]*m2arr[n]/(avgfac*diagmes[n])
              << "  ";
    for (int c = 0; c < nex; ++c) {
      std::cout << "(" << flowbasis[n][c][0]
                << "," << flowbasis[n][c][1] << ")";
    }
    std::cout << "  " << colbasis.at(n) << std::endl;
  }
  std::cout << std::endl
    << "OpenLoops amplitude squared with ColorFull: " << me << std::endl
    << "Squared matrix element from OpenLoops:      " << m2 << std::endl
    << "ratio:                                      " << me/m2 << std::endl;

  return 0;
}
