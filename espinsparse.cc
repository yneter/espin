#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>
#include <cmath>
#include <vector>

#ifndef ESPINSPARSE
#define ESPINSPARSE

#include "espingen.cc"

struct SpinSparse { 
    typedef SparseMatrix<complexg> SpinMatrix;
    typedef SparseMatrix<double> SpinMatrixReal;
private: 
    int matrix_size;

    SpinMatrixReal rsx;
    SpinMatrixReal isy;
    SpinMatrixReal rsz;
    SpinMatrixReal rsp;
    SpinMatrixReal rsm;
    SpinMatrixReal rid;

    SpinMatrix sx;
    SpinMatrix sy;
    SpinMatrix sz;
    SpinMatrix sp;
    SpinMatrix sm;
    SpinMatrix id;
    SpinMatrix Hfull;

    void init_matrix(void) {     
       id.setIdentity();
       rid.setIdentity();

       sx.setZero();
       sy.setZero();
       sz.setZero();
       sp.setZero();
       sm.setZero();

       rsx.setZero();
       isy.setZero();
       rsz.setZero();
       rsp.setZero();
       rsm.setZero();


       double S = ((double)matrix_size - 1.)/2.;
       std::vector< Triplet<double> > coef;      
       for (int i = 0; i < matrix_size; i++) { 
	  double mz = S - (double)i;
	  coef.push_back( Triplet<double>( i, i, mz ) );
       }
       rsz.setFromTriplets(coef.begin(), coef.end());              
       sz.setFromTriplets(coef.begin(), coef.end());              


       coef.clear();
       for (int i = 0; i < matrix_size-1; i++) { 
	  double mz = S - (double)i;
	  double nz = S - (double)(i+1);
	  coef.push_back( Triplet<double>( i, i+1, sqrt(S*(S+1) - mz * nz) ) );
       }
       sp.setFromTriplets(coef.begin(), coef.end());              
       rsp.setFromTriplets(coef.begin(), coef.end());              


       coef.clear();
       for (int i = 0; i < matrix_size-1; i++) { 
	  double mz = S - (double)i;
	  double nz = S - (double)(i+1);
	  coef.push_back( Triplet<double>( i+1, i, sqrt(S*(S+1) - mz * nz) ) );
       }
       sm.setFromTriplets(coef.begin(), coef.end());              
       rsm.setFromTriplets(coef.begin(), coef.end());              
       
       rsx = (rsp + rsm)/2.;
       isy = (rsp - rsm)/2.;

       sx = (sp + sm)/2.;
       sy = (sp - sm)/(2.*iii);
  
    }


    void resize_matrix(int N) { 
       matrix_size = N;
       sx.resize(matrix_size, matrix_size);
       sy.resize(matrix_size, matrix_size);
       sz.resize(matrix_size, matrix_size);
       id.resize(matrix_size, matrix_size);
       sp.resize(matrix_size, matrix_size);
       sm.resize(matrix_size, matrix_size);
       rsx.resize(matrix_size, matrix_size);
       isy.resize(matrix_size, matrix_size);
       rsz.resize(matrix_size, matrix_size);
       rid.resize(matrix_size, matrix_size);
       rsp.resize(matrix_size, matrix_size);
       rsm.resize(matrix_size, matrix_size);
       Hfull.resize(matrix_size, matrix_size);
    }

public :       
    double D; // D 
    double E; // E
    Rotation rot; // rotation from laboratory frame to molecular frame 
    Vector3d g3; // g factors in molecular frame 
    Vector3d B;  // magnetic field in laboratory frame 


    SpinSparse(int N): matrix_size(N), 
		       sx(matrix_size, matrix_size),
		       sy(matrix_size, matrix_size),
		       sz(matrix_size, matrix_size),
		       id(matrix_size, matrix_size),
		       sp(matrix_size, matrix_size),
		       sm(matrix_size, matrix_size),
		       rsx(matrix_size, matrix_size),
		       isy(matrix_size, matrix_size),
		       rsz(matrix_size, matrix_size),
		       rid(matrix_size, matrix_size),
		       rsp(matrix_size, matrix_size),
		       rsm(matrix_size, matrix_size),
		       Hfull(matrix_size, matrix_size),
		       D(0),
		       E(0),
		       rot ( Matrix3d() = Matrix3d::Identity() ),
		       g3 ( Vector3d() = Vector3d::Constant(1.0) ),
		       B (  Vector3d() = Vector3d::Constant(0.0) )

    { 

       init_matrix();
    }

    SpinSparse(void) : matrix_size(0), 
		       D(0),
		       E(0),
		       rot ( Matrix3d() = Matrix3d::Identity() ),
		       g3 ( Vector3d() = Vector3d::Constant(1.0) ),
		       B (  Vector3d() = Vector3d::Constant(0.0) )
    {
    }

    const SpinMatrixReal& rSx(void) const { return rsx; }
    const SpinMatrixReal& iSy(void) const { return isy; }
    const SpinMatrixReal& rSz(void) const { return rsz; }
    const SpinMatrixReal& rSp(void) const { return rsp; }
    const SpinMatrixReal& rSm(void) const { return rsm; }
    const SpinMatrixReal& rId(void) const { return rid; }

    const SpinMatrix& Sx(void) const { return sx; }
    const SpinMatrix& Sy(void) const { return sy; }
    const SpinMatrix& Sz(void) const { return sz; }
    const SpinMatrix& Sp(void) const { return sp; }
    const SpinMatrix& Sm(void) const { return sm; }
    const SpinMatrix& Id(void) const { return id; }

    int size(void) const { return matrix_size; }

    void resize(int N) { resize_matrix(N); init_matrix(); }

    const SpinMatrix& update_hamiltonian(void) { 
       Matrix3d r_matrix = rot.matrix();
       SpinMatrix rotSx = r_matrix(0, 0) * Sx() - (iii * r_matrix(0, 1)) * iSy() + r_matrix(0, 2) * Sz();
       SpinMatrix rotSy = r_matrix(1, 0) * Sx() - (iii * r_matrix(1, 1)) * iSy() + r_matrix(1, 2) * Sz();
       SpinMatrix rotSz = r_matrix(2, 0) * Sx() - (iii * r_matrix(2, 1)) * iSy() + r_matrix(2, 2) * Sz();
       Vector3d rBvec = r_matrix * B;
       Hfull = D * (rotSz * rotSz - 2.0*Id()/3.0) + E * (rotSx * rotSx -  rotSy * rotSy) 
	   + g3[0] * rotSx * rBvec[0] + g3[1] * rotSy * rBvec[1] + g3[2] * rotSz * rBvec[2];
       return Hfull;
    }

    const SpinMatrix& hamiltonian(void) const { 
       return Hfull;
    }

    SpinSparse& operator=(const SpinSparse &S) { 
       this->D = S.D;
       this->E = S.E;
       this->rot = S.rot;
       this->g3 = S.g3;
       this->B = S.B;
       update_hamiltonian();
       return *this;
    }  

};

struct SpinHalfSparse : public SpinSparse { 
   enum { matrix_size = 2 };
   SpinHalfSparse() : SpinSparse(matrix_size) {  }
};

struct SpinTupleSparse { 
    typedef SparseMatrix< complexg > SpinMatrix;
    typedef SparseMatrix<double> SpinMatrixReal;

private:
    int matrix_size;
    std::vector< SpinSparse > Svec;
    int compute_size(void) { 
       int s = 1;
       for (int j = 0; j < Svec.size(); j++) { 
	  s *= Svec[j].size();
       }
       return s;
    }
    SpinMatrixReal rid;
    SpinMatrixReal rzero;
    SpinMatrix id;
    SpinMatrix zero;
   
    void resize_matrix(void) { 
       rid.resize(matrix_size, matrix_size);
       rid.setIdentity();
       rzero.resize(matrix_size, matrix_size);
       rzero.setZero();
       id.resize(matrix_size, matrix_size);
       id.setIdentity();
       zero.resize(matrix_size, matrix_size);
       zero.setZero();
    }
public :
    int size(void) { return matrix_size; }

    int nspins(void) const { return Svec.size(); }

    std::vector<int> spin_size_list(void) const { 
       std::vector<int> ss( nspins() );
       for (int j = 0; j < nspins(); j++) { 
	  ss[j] = Svec[j].size();
       }
       return ss;
    }

    void spin_size_list(const std::vector<int> &spins) { 
       Svec.clear();
       for (int i = 0; i < spins.size(); i++) { 
	  Svec.push_back( SpinSparse(spins[i]) );
       }
       matrix_size = compute_size();
       resize_matrix();
    }


    SpinTupleSparse(void) { 

    }

    SpinTupleSparse(const std::vector<int> &spins) { 
       spin_size_list(spins);
    } 

    int find_left_matrix_size(int i) { 
       int s = 1;
       for (int j = 0; j < i; j++) { 
	  s *= Svec[j].size();
       }
       return s;
    };

    int find_item_matrix_size(int i) { 
       return Svec[i].size();
    }

    template <typename Matrix> Matrix make_matrix_Hi(int i, const Matrix &Hi) { 
       int size_left = find_left_matrix_size(i);
       int size_i = find_item_matrix_size(i);
       int size_right = matrix_size / (size_left * size_i);

       SpinMatrixReal id_left(size_left, size_left);
       id_left.setIdentity();
       SpinMatrixReal id_right(size_right, size_right);
       id_right.setIdentity();

       return kroneckerProduct(id_left, kroneckerProduct(Hi, id_right)).eval();
    } 


    SpinMatrix make_matrix_HiHj(int I, SpinMatrix Hi, int J, SpinMatrix Hj ) { 
       int size_I_left = find_left_matrix_size(I);
       int size_I = find_item_matrix_size(I);
       int size_J_left = find_left_matrix_size(J);
       int size_J = find_item_matrix_size(J);
       int size_J_right = matrix_size / (size_J_left * size_J);
       int size_center = size_J_left / (size_I_left * size_I);

       SpinMatrixReal id_left(size_I_left, size_I_left);
       id_left.setIdentity();
       SpinMatrixReal id_right(size_J_right, size_J_right);
       id_right.setIdentity();
       SpinMatrixReal id_center(size_center, size_center);
       id_center.setIdentity();

       return kroneckerProduct(id_left, 
			       kroneckerProduct(Hi, 
						kroneckerProduct(id_center, 
								 kroneckerProduct(Hj, id_right)
								 ))).eval();
    }

    SpinSparse &S(int i) { return Svec[i]; }
    const SpinSparse &Sconst(int i) const { return Svec[i]; }

    SpinMatrixReal rSx(int i) { return make_matrix_Hi<SpinMatrixReal>(i, Svec[i].rSx()); }
    SpinMatrixReal iSy(int i) { return make_matrix_Hi<SpinMatrixReal>(i, Svec[i].iSy()); }
    SpinMatrixReal rSz(int i) { return make_matrix_Hi<SpinMatrixReal>(i, Svec[i].rSz()); }
    SpinMatrixReal rSp(int i) { return make_matrix_Hi<SpinMatrixReal>(i, Svec[i].rSp()); }
    SpinMatrixReal rSm(int i) { return make_matrix_Hi<SpinMatrixReal>(i, Svec[i].rSm()); }
    SpinMatrix Sx(int i) { return make_matrix_Hi<SpinMatrix>(i, Svec[i].Sx()); }
    SpinMatrix Sy(int i) { return make_matrix_Hi<SpinMatrix>(i, Svec[i].Sy()); }
    SpinMatrix Sz(int i) { return make_matrix_Hi<SpinMatrix>(i, Svec[i].Sz()); }
    SpinMatrix Sp(int i) { return make_matrix_Hi<SpinMatrix>(i, Svec[i].Sp()); }
    SpinMatrix Sm(int i) { return make_matrix_Hi<SpinMatrix>(i, Svec[i].Sm()); }

    const SpinMatrixReal& rId(void) { return rid; } 
    const SpinMatrixReal& rZero(void) { return rzero; } 
    const SpinMatrix& Id(void) { return id; } 
    const SpinMatrix& Zero(void) { return zero; } 
    
    SpinMatrixReal rSx(void) { SpinMatrixReal srSx = rZero(); for (int s = 0; s < nspins(); s++) srSx += rSx(s); return srSx; }
    SpinMatrixReal iSy(void) { SpinMatrixReal siSy = rZero(); for (int s = 0; s < nspins(); s++) siSy += iSy(s); return siSy; }
    SpinMatrixReal rSz(void) { SpinMatrixReal srSz = rZero(); for (int s = 0; s < nspins(); s++) srSz += rSz(s); return srSz; }
    SpinMatrix Sx(void) { SpinMatrix sSx = Zero(); for (int s = 0; s < nspins(); s++) sSx += Sx(s); return sSx; }
    SpinMatrix Sy(void) { SpinMatrix sSy = Zero(); for (int s = 0; s < nspins(); s++) sSy += Sy(s); return sSy; }
    SpinMatrix Sz(void) { SpinMatrix sSz = Zero(); for (int s = 0; s < nspins(); s++) sSz += Sz(s); return sSz; }



    SpinTupleSparse& operator=(const SpinTupleSparse &S) { 
       this->spin_size_list( S.spin_size_list() );
       for (int s = 0; s < nspins(); s++) { 
	  Svec[s] = S.Sconst(s);
       }
       return *this;
    }  
};


#endif

