#ifndef ESPINCAVITY
#define ESPINCAVITY

#include "espinsparse.cc"
#include <iostream>
#include <tuple>

class CavitySparse { 
public:
    typedef double ecavity_float;
    typedef SparseMatrix<ecavity_float> CavityMatrix;
private :
    CavityMatrix a_matrix;
    CavityMatrix ap_matrix;
    CavityMatrix N_matrix;
    CavityMatrix Id_matrix;
    CavityMatrix Hfull;

public :
    double omega_c;
    double Fx;
  
    CavitySparse(int N) : 
      a_matrix(N, N),
      ap_matrix(N, N),
      N_matrix(N, N),
      Id_matrix(N, N),
      Hfull(N, N)
    { 
       Fx = 0;
       std::vector< Triplet<ecavity_float> > coef;      
       for (int i = 0; i < N-1; i++) { 
	 coef.push_back( Triplet<ecavity_float>( i, i+1, static_cast<ecavity_float>(sqrt(i+1)) ) );
       }
       a_matrix.setFromTriplets(coef.begin(), coef.end());

       coef.clear();
       for (int i = 0; i < N-1; i++) { 
	 coef.push_back( Triplet<ecavity_float>( i+1, i, static_cast<ecavity_float>(sqrt(i+1)) ) );
       }
       ap_matrix.setFromTriplets(coef.begin(), coef.end());


       coef.clear();
       for (int i = 0; i < N; i++) { 
	  coef.push_back( Triplet<ecavity_float>( i, i, static_cast<ecavity_float>(i) ) );
       }
       N_matrix.setFromTriplets(coef.begin(), coef.end());
       
       Id_matrix.setIdentity();
    }

    int size(void) const { 
       return Id_matrix.cols();
    }

    const CavityMatrix& a(void) const { 
       return a_matrix;
    }

    const CavityMatrix& ap(void) const { 
       return ap_matrix;
    }

    const CavityMatrix& N(void) const { 
       return N_matrix;
    }

    const CavityMatrix& Id(void) const {
       return Id_matrix;
    }

    VectorXcd coherent(complex<double> alpha) { 
       VectorXcd c(size());
       c.setZero();
       if (abs(alpha) < 1e-15 ) {
	 c(0) = 1;
       } else { 
	 double an0 = norm(alpha)/2.0;
	 for (int i = 0; i < size(); i++) {
	   c(i) = exp( -an0 + log(alpha) * (double)i - 0.5 * std::lgamma(i+1) );
	 }
       }
       return c;
    }

    VectorXcd vac_n(int n) { 
	VectorXcd v( size() );
	v.setZero();
	v(n) = 1;
	return v;
    }
  
    VectorXcd coherent(double alphax, double alphay) { 
       std::complex<double> alpha(alphax, alphay);
       return coherent(alpha);
    }
    MatrixXcd coherent_rho(complex<double> alpha) { 
       VectorXcd c = coherent(alpha);
       return c * c.adjoint();
    }

    void update_hamiltonian(void) { 
       Hfull = omega_c * N_matrix + Fx * (a_matrix + ap_matrix);
    }

    const SparseMatrix<double>& hamiltonian(void) const { 
       return Hfull;
    }

};



//
// too much code overlap with SpinTupleSparse

struct CavityTupleSparse : public QuantTupleSparse< CavitySparse > { 
    typedef CavitySparse::CavityMatrix  CavityMatrix;
private:
    CavityMatrix Hfull; 
    CavityMatrix Id_matrix;
public: 

    CavityTupleSparse(const std::vector<int> &vec) : QuantTupleSparse< CavitySparse >(vec) { 
       Hfull.resize( size(), size() );
       Id_matrix.resize( size(), size() );
       Id_matrix.setIdentity();
    } 

    CavityMatrix a(int i) const { return make_matrix_Hi<CavityMatrix>(i, quant_sys_const(i).a()); }
    CavityMatrix ap(int i) const { return make_matrix_Hi<CavityMatrix>(i, quant_sys_const(i).ap()); }
    CavityMatrix N(int i) const { return make_matrix_Hi<CavityMatrix>(i, quant_sys_const(i).N()); }
    const CavityMatrix &Id(void) const { return Id_matrix; }
    
    CavitySparse &cavity(int i) { return quant_sys(i); }
    const CavitySparse &cavity_const(int i) const { return quant_sys_const(i); }

    void load_uncoupled_hamiltonian(void) { 
       Hfull.setZero();
       for (int s = 0; s < nsys(); s++) {
	  cavity(s).update_hamiltonian();
	  Hfull += make_matrix_Hi<CavityMatrix>(s, cavity(s).hamiltonian() );
       }       
    }

    void update_hamiltonian(void) { 
       load_uncoupled_hamiltonian();
    }

    const CavityMatrix& hamiltonian(void) const { 
       return Hfull;
    }
};


struct SpinTupleInCavity  { 
    std::vector<double> Bz;
    std::vector<double> Omega;
    double omegaC;
    double F;

    typedef SparseMatrix<double> SIC_Matrix;
    SpinTupleSparse S;
    CavitySparse os;
    bool rwa;
private:
    int matrix_size;
    SIC_Matrix Hfull;
public :
    int size(void) const { return matrix_size; }
    int nspins(void) const { return S.nspins(); } 

private:
    void resize_nspins(void) { 
       Bz.resize( nspins() );
       Omega.resize( nspins() );
    }
public: 

    SpinTupleInCavity(int N, std::vector<int> spin_size_list) : S(spin_size_list),
			       os(N)
    { 
       matrix_size = S.size() * os.size();
       Hfull.setZero();
       resize_nspins();
       rwa = true;
    }


    int enc_nl_s(int nl, int s) { 
       return nl * S.size() + s;
    }


    int dec_nl(int n) {
       return n / S.size();
    }

    int dec_s(int n) {
       return n % S.size();
    }


    const SIC_Matrix& update_hamiltonian(void) { 
       Hfull = omegaC * kroneckerProduct(os.N(), S.rId()) + F * ( kroneckerProduct(os.a(), S.rId()) + kroneckerProduct(os.ap(), S.rId()) );
       for (int i = 0; i < S.nspins(); i++) { 
	  Hfull += ( Bz[i] * kroneckerProduct(os.Id(), S.rSz(i)) + Omega[i] *  ( kroneckerProduct(os.a(), S.rSp(i)) + kroneckerProduct(os.ap(), S.rSm(i))) );
	  if (!rwa) { 
	     Hfull += Omega[i] *  ( kroneckerProduct(os.ap(), S.rSp(i)) + kroneckerProduct(os.a(), S.rSm(i))) ;
	  }
       }
       return Hfull;
    }

    const SIC_Matrix& hamiltonian(void) const { 
       return Hfull;
    }

    void operator()(const VectorXcd &x , VectorXcd &dxdt, double t)
    {
        dxdt = -iii * hamiltonian() * x;
    }


    complexg trace_AB(const SparseMatrix<double>& A, const MatrixXcd &B) { return (A * B).eval().trace(); }
    complexg trace_cAB(const SparseMatrix<complexg>& A, const MatrixXcd &B) { return (A * B).eval().trace(); }
    double re_trace_AB(const SparseMatrix<double>& A, const MatrixXcd &B) { return real( trace_AB(A, B) ); }
    complexg trA(const MatrixXcd &rho) { return trace_AB( kroneckerProduct(os.a(), S.rId()), rho); }
    double trN(const MatrixXcd &rho) { return re_trace_AB( kroneckerProduct(os.N(), S.rId()), rho); }
    double trS2(const MatrixXcd &rho) { return re_trace_AB( kroneckerProduct(os.Id(), S.rSx() * S.rSx() - S.iSy() * S.iSy() + S.rSz() * S.rSz()), rho); }
    double trSx(const MatrixXcd &rho) { return re_trace_AB( kroneckerProduct(os.Id(), S.rSx()), rho); }
    double trSy(const MatrixXcd &rho) { return imag ( trace_AB( kroneckerProduct(os.Id(), S.iSy()), rho) ); }
    double trSz(const MatrixXcd &rho) { return re_trace_AB( kroneckerProduct(os.Id(), S.rSz()), rho); }
    double trSx(int i, const MatrixXcd &rho) { return re_trace_AB( kroneckerProduct(os.Id(), S.rSx(i)), rho); }
    double trSy(int i, const MatrixXcd &rho) { return imag ( trace_AB( kroneckerProduct(os.Id(), S.iSy(i)), rho) ); }
    double trSz(int i, const MatrixXcd &rho) { return re_trace_AB( kroneckerProduct(os.Id(), S.rSz(i)), rho); }
    double trBell_spin_cavity(int i, const MatrixXcd &rho) { 
       SparseMatrix<complexg> MQ = 2. * kroneckerProduct(os.Id(), S.Sz(i));
       SparseMatrix<complexg> MR = 2. * kroneckerProduct(os.Id(), S.Sx(i));
       double E = sqrt( trN(rho) + 0.5 );
       SparseMatrix<complexg> MS = kroneckerProduct(os.a() + os.ap(), S.Id())  / (2. * E); 
       SparseMatrix<complexg> MT = kroneckerProduct(iii * (os.a() - os.ap()), S.Id())  / (2. * E);
       return real( trace_cAB(MQ * MS + MR * MS + MR * MT - MQ * MT, rho) );
    }


    double trBell_spin(int i, int j, 
		       const Vector3d &s1, const Vector3d &s2, 
		       const Vector3d &s3, const Vector3d &s4, const MatrixXcd &rho) { 
          Vector3d v1, v2, v3, v4;
	  v1 = s1;
	  v2 = s2;
	  v3 = s3;
	  v4 = s4;
	  v1.normalize();
	  v2.normalize();
	  v3.normalize();
	  v4.normalize();
	  SparseMatrix<complexg> MQ = 2. * kroneckerProduct(os.Id(), v1(0) * S.Sx(i) + v1(1) * S.Sy(i) + v1(2) * S.Sz(i));
	  SparseMatrix<complexg> MR = 2. * kroneckerProduct(os.Id(), v2(0) * S.Sx(i) + v2(1) * S.Sy(i) + v2(2) * S.Sz(i));
	  SparseMatrix<complexg> MS = 2. * kroneckerProduct(os.Id(), v3(0) * S.Sx(j) + v3(1) * S.Sy(j) + v3(2) * S.Sz(j));
	  SparseMatrix<complexg> MT = 2. * kroneckerProduct(os.Id(), v4(0) * S.Sx(j) + v4(1) * S.Sy(j) + v4(2) * S.Sz(j));
	  return fabs( real( trace_cAB(MQ * MS + MR * MS + MR * MT - MQ * MT, rho) ) );
    }

    double trBell_spin_random(int i, int j, int Nsample, const MatrixXcd &rho) { 
       Rotation rot1, rot2;
       Rotation rot1max, rot2max;
       Matrix3d r1, r2;
       double bmax = 0;
       double b; 

       for (int k=0; k < Nsample; k++) {
	  rot1.random();
	  rot2.random();
	  r1 = rot1.matrix();
	  r2 = rot2.matrix();

	  b = trBell_spin(i, j, r1.row(0), r1.row(1), r2.row(0), r2.row(1), rho);
	  if (b > bmax) { 
	     bmax = b;
	     rot1max = rot1;
	     rot2max = rot2;
	  }

       }

       r1 = rot1max.matrix();
       r2 = rot2max.matrix();
       std::cout << "# bmax " << bmax << std::endl;
       std::cout << "# vector1a " << r1(0,0) << "   " << r1(0,1) << "     " << r1(0,2) << std::endl;
       std::cout << "# vector1b " << r1(1,0) << "   " << r1(1,1) << "     " << r1(1,2) << std::endl;
       std::cout << "# vector2a " << r2(0,0) << "   " << r2(0,1) << "     " << r2(0,2) << std::endl;
       std::cout << "# vector2b " << r2(1,0) << "   " << r2(1,1) << "     " << r2(1,2) << std::endl;
       return bmax;
    }



    void print_info(void) { 
       std::cout << "# Bz ";
       for (int s = 0; s < nspins(); s++) std::cout << Bz[s] << "   ";
       std::cout << std::endl;
       std::cout << "# Omega ";
       for (int s = 0; s < nspins(); s++) std::cout << Omega[s] << "   ";
       std::cout << std::endl;
       std::cout << "# Spin size ";
       for (int s = 0; s < nspins(); s++) std::cout << S.S(s).size() << "   ";
       std::cout << std::endl;
       std::cout << "# omegaC " << omegaC << std::endl;
       std::cout << "# F " << F << std::endl;
       std::cout << "# os size " << os.size() << std::endl;
       std::cout << "# spin size " << S.size() << std::endl;
       std::cout << "# rwa " << rwa << std::endl;
    }
};


struct SpinCavityTuple  { 
    SpinTupleSparse S;
    CavityTupleSparse os;

    int size(void) { return S.size() * os.size(); }
private :
    SparseMatrix<complexg> Hfull;
public : 
  typedef SparseMatrix< std::complex<double> > SpinCavityTupleMatrix;

    SpinCavityTuple(std::vector<int> spins, std::vector<int> cavities) : S(spins), os(cavities) {        
       Hfull.resize( size(), size() );
    }
    SpinCavityTupleMatrix Sx(int i) const { return kroneckerProduct(S.Sx(i), os.Id()); }
    SpinCavityTupleMatrix Sy(int i) const { return kroneckerProduct(S.Sy(i), os.Id()); }
    SpinCavityTupleMatrix Sz(int i) const { return kroneckerProduct(S.Sz(i), os.Id()); }
    SpinCavityTupleMatrix Sx2(int i) const { return kroneckerProduct(S.Sx2(i), os.Id()); }
    SpinCavityTupleMatrix Sy2(int i) const { return kroneckerProduct(S.Sy2(i), os.Id()); }
    SpinCavityTupleMatrix Sz2(int i) const { return kroneckerProduct(S.Sz2(i), os.Id()); }
    SpinCavityTupleMatrix N(int i) const { return kroneckerProduct(S.Id(), os.N(i)); }

    double trSx(int i, const VectorXcd &psi) { return real( std::complex<double>( psi.adjoint() * Sx(i) * psi) ); }
    double trSy(int i, const VectorXcd &psi) { return real( std::complex<double>( psi.adjoint() * Sy(i) * psi) ); }
    double trSz(int i, const VectorXcd &psi) { return real( std::complex<double>( psi.adjoint() * Sz(i) * psi) ); }
    double trSx2(int i, const VectorXcd &psi) { return real( std::complex<double>( psi.adjoint() * Sx2(i) * psi) ); }
    double trSy2(int i, const VectorXcd &psi) { return real( std::complex<double>( psi.adjoint() * Sy2(i) * psi) ); }
    double trSz2(int i, const VectorXcd &psi) { return real( std::complex<double>( psi.adjoint() * Sz2(i) * psi) ); }
    double trN(int i, const VectorXcd &psi) { return real( std::complex<double>( psi.adjoint() * N(i) * psi) ); }

  
    void update_hamiltonian(void) {
       S.update_hamiltonian();
       os.update_hamiltonian();
       Hfull = kroneckerProduct( S.hamiltonian(), os.Id() ).eval() + kroneckerProduct( S.Id(), os.hamiltonian() ).eval();
    }

    SparseMatrix<complexg> &hamiltonian(void) { 
       return Hfull; 
    }

};

struct SpinTupleInCavityMaster;

#ifdef ESPIN_NOT_WORKING

namespace Eigen {
namespace internal {
  template <> struct traits< SpinTupleInCavityMaster > : public Eigen::internal::traits<Eigen::SparseMatrix<complexg> >
  {};

}
}

#endif

struct SpinTupleInCavityMaster;

struct SpinTupleInCavityMasterDt { 
    const SpinTupleInCavityMaster &spin_in_cavity;
    void operator()(const MatrixXcd &x , MatrixXcd &dxdt, double t);
  
    SpinTupleInCavityMasterDt(const SpinTupleInCavityMaster &spin_in_cavity_ref) : 
       spin_in_cavity(spin_in_cavity_ref) { 
    }
};

#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

#include "espinsylvester.cc"

struct SpinTupleInCavityMaster : public SpinTupleInCavity
{ 
    double gamma;
    double nav;
    double gammaS;
    VectorXd gammaSvec;
    double gammaSup;
    double temp;
private:
    enum { Lrate, Lop, Ladjoint, Ladjoint_Lop };
    enum { gamma_os, gamma_os_up, gamma_cavity }; 
    int number_of_lindblad_operators;

    typedef std::tuple< double, SparseMatrix< double >, SparseMatrix< double >, SparseMatrix< double > > LindbladOperator;
    std::vector < LindbladOperator > lindblad_operators;    

    int master_size;
    MatrixXcd Jacobi_for_Liouville;
    SparseMatrix<complexg> Master_full;

    MatrixXcd rho_matr;
public:
    SpinTupleInCavityMaster(int N, std::vector<int> spin_size_list) : SpinTupleInCavity(N, spin_size_list),
			  rho_matr(SpinTupleInCavity::size(), SpinTupleInCavity::size()),
			  master_size(SpinTupleInCavity::size() * SpinTupleInCavity::size()),
			  number_of_lindblad_operators( gamma_cavity + 2 * S.nspins() ),
			  lindblad_operators( number_of_lindblad_operators ),
			  Master_full(master_size, master_size)
    { 
      //       Jacobi_for_Liouville = MatrixXcd::Constant(SpinTupleInCavity::size() , SpinTupleInCavity::size(), 1.0);
       gammaS = 0;
       gammaSup = 0;
       gamma = 0;
       nav = 0;

       std::get<Lop>( lindblad_operators[gamma_os] ) = kroneckerProduct(os.a(), S.rId());
       std::get<Ladjoint>( lindblad_operators[gamma_os] ) = kroneckerProduct(os.ap(), S.rId());
       std::get<Ladjoint_Lop>( lindblad_operators[gamma_os] ) = kroneckerProduct(os.N(), S.rId());       

       std::get<Lop>( lindblad_operators[gamma_os_up] ) = kroneckerProduct(os.ap(), S.rId());
       std::get<Ladjoint>( lindblad_operators[gamma_os_up] ) = kroneckerProduct(os.a(), S.rId());
       std::get<Ladjoint_Lop>( lindblad_operators[gamma_os_up] ) = kroneckerProduct(os.a() * os.ap(), S.rId());       

       for (int s = 0; s < S.nspins(); s++) { 
	  int gamma_spin = gamma_cavity + 2 * s;
	  int gamma_spin_up = gamma_spin + 1;

	  std::get<Lop>( lindblad_operators[gamma_spin] ) = kroneckerProduct(os.Id(), S.rSm(s));
	  std::get<Ladjoint>( lindblad_operators[gamma_spin] ) = kroneckerProduct(os.Id(), S.rSp(s));
	  std::get<Ladjoint_Lop>( lindblad_operators[gamma_spin] ) = kroneckerProduct(os.Id(), S.rSp(s) * S.rSm(s) );       
	  
	  std::get<Lop>( lindblad_operators[gamma_spin_up] ) = kroneckerProduct(os.Id(), S.rSp(s));
	  std::get<Ladjoint>( lindblad_operators[gamma_spin_up] ) = kroneckerProduct(os.Id(), S.rSm(s));
	  std::get<Ladjoint_Lop>( lindblad_operators[gamma_spin_up] ) = kroneckerProduct(os.Id(), S.rSm(s) * S.rSp(s) );       
       }
    }


    typedef complexg Scalar;
    typedef double RealScalar;
    typedef int StorageIndex;

    enum {
       ColsAtCompileTime = Eigen::Dynamic,
       MaxColsAtCompileTime = Eigen::Dynamic,
       IsRowMajor = false
    };
  
    Index rows() const { return master_size; }
    Index cols() const { return master_size; } 

    template<typename Rhs>
    Eigen::Product<SpinTupleInCavityMaster,Rhs,Eigen::AliasFreeProduct> 
    operator*(const Eigen::MatrixBase<Rhs>& x) const {
       return Eigen::Product<SpinTupleInCavityMaster,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
    }



    int enc_master(int n, int m) { 
       return n * SpinTupleInCavity::size() + m;
    }

    template <typename Type> MatrixXcd map_to_mat(const Type &vec) const { 
       return Map< const MatrixXcd >(vec.data(), SpinTupleInCavity::size(), SpinTupleInCavity::size());
    }

    template <typename Type> VectorXcd map_to_vec(const Type &mat) const { 
       return Map< const VectorXcd > (mat.data(), master_size);
    }

    const SparseMatrix<complexg>& master(void) { 
       return Master_full;
    }

    double gammaSpin(int i) { 
       if (gammaSvec.size() == S.nspins()) {
	 return gammaSvec(i);
       } else {
	 return gammaS;
       }
    }
  
    void update_jacobi_preconditioner(void) { 
      /***
       int matrix_size = SpinTupleInCavity::size(); 
       MatrixXcd rho_tmp( matrix_size, matrix_size );
       for (int i = 0; i < matrix_size ; i++) { 
	  rho_tmp.setZero();
	  rho_tmp(i, i) = 1.0;
	  rho_tmp = liouvillian( rho_tmp );
	  if (abs(rho_tmp(i, i)) > 1e-2) { 
	    Jacobi_for_Liouville(i, i) = 1.0/ rho_tmp(i, i);
	  }
       }  
      ***/
    }

    void update_lindblad(void) {  
       std::get<Lrate>( lindblad_operators[gamma_os] ) = gamma * ( nav + 1. );
       std::get<Lrate>( lindblad_operators[gamma_os_up] ) = gamma * nav;
       for (int s = 0; s < S.nspins(); s++) { 
	  int gamma_spin = gamma_cavity + 2 * s;
	  int gamma_spin_up = gamma_spin + 1;
	  if (gammaSvec.size() == S.nspins()) {
	     std::get<Lrate>( lindblad_operators[gamma_spin] ) = gammaSvec(s);
	  } else {
	     std::get<Lrate>( lindblad_operators[gamma_spin] ) = gammaS;
	  }
	  std::get<Lrate>( lindblad_operators[gamma_spin_up] ) = gammaSup;
       }
    }

    const SparseMatrix< double >& update_hamiltonian(void) { 
       update_lindblad();
       return SpinTupleInCavity::update_hamiltonian();
    }
    
    const SparseMatrix<complexg>& update_master(void) { 
       update_hamiltonian();
       Master_full.setZero();

       std::vector< Triplet<complexg> > coef;      

       for (int k=0; k<hamiltonian().outerSize(); k++) { 
 	  for (SparseMatrix<double>::InnerIterator it(hamiltonian(),k); it; ++it) {
	     int n = it.row();   // row index
	     int m = it.col();   // col index (here it is equal to k)
	     for (int c = 0; c < SpinTupleInCavity::size(); c++) { 

	       // d rho(n, c)/dt = -i H(n, m) rho(m, c) 

	       coef.push_back( Triplet<complexg>( enc_master(n, c), enc_master(m, c) , - iii * it.value() ) );

	       // d rho(c, m)/dt = i rho(c, n) H(n, m) 
	       coef.push_back( Triplet<complexg>( enc_master(c, m), enc_master(c, n) ,  iii * it.value() ) );

	     }
	  }
       }


       for (int n1 = 0; n1 < SpinTupleInCavity::size(); n1++) { 
	  for (int m1 = 0; m1 < SpinTupleInCavity::size(); m1++) { 
	     int nl = this->dec_nl(n1);
	     int ml = this->dec_nl(m1);
	     coef.push_back( Triplet<complexg>( enc_master(n1, m1), enc_master(n1, m1) , -0.5 * gamma * (double) ( nl + ml ) ) );
	    
	     if (nl+1 < this->os.size() && ml+1 < this->os.size() ) { 
	        int sn = this->dec_s(n1);
		int sm = this->dec_s(m1);
		int n2 = this->enc_nl_s(nl+1, sn);
		int m2 = this->enc_nl_s(ml+1, sm);
		coef.push_back( Triplet<complexg>( enc_master(n1, m1), enc_master(n2, m2) , gamma * sqrt( (nl+1) * (ml+1) ) ) );
	     }
	  }
       }

       Master_full.setFromTriplets(coef.begin(), coef.end());
       return Master_full;
    }

    template <class MatrixIn, class MatrixOut> MatrixOut liouvillian_lindbladt_gen(const MatrixIn &rho) const { 
       MatrixOut Lrho(SpinTupleInCavity::size(), SpinTupleInCavity::size());
       Lrho.setZero();
       for (int n = 0; n < number_of_lindblad_operators; n++) { 
	 Lrho += std::get<Lrate>(lindblad_operators[n]) * ( 
		  std::get<Lop>(lindblad_operators[n]) * rho * std::get<Ladjoint>(lindblad_operators[n]) 
		  - 0.5 * std::get<Ladjoint_Lop>(lindblad_operators[n]) * rho 
		  - 0.5 * rho * std::get<Ladjoint_Lop>(lindblad_operators[n]) );
       }
       return Lrho;
    }

  // -iii H rho - 1/2 L^+ L rho = -iii * (H - iii L^+ L / 2) rho 
    template <class MatrixIn, class MatrixOut> MatrixOut liouvillian_gen(const MatrixIn &rho) const { 
       return -iii * hamiltonian() * rho + iii * rho * hamiltonian() + liouvillian_lindbladt_gen<MatrixIn,MatrixOut>(rho);
    }

    MatrixXcd liouvillian(const MatrixXcd &rho) const { return liouvillian_gen<MatrixXcd, MatrixXcd>(rho); }
    MatrixXcd liouvillian_lindbladt(const MatrixXcd &rho) const { return liouvillian_lindbladt_gen<MatrixXcd, MatrixXcd>(rho); }
    MatrixXcd liouvillian_hamiltonian(const MatrixXcd &rho) const { return -iii * hamiltonian() * rho + iii * rho * hamiltonian(); }

    SparseMatrix<complexg> liouvillian_sparse(const SparseMatrix< complexg > &rho) const { 
       return liouvillian_gen< SparseMatrix<complexg>, SparseMatrix<complexg> >(rho);
    }

    const SparseMatrix<complexg>&  update_master_from_liouville(void) { 
       update_hamiltonian();
       
       Master_full.setZero();
       std::vector< Triplet<complexg> > coef;      


       for (int k=0; k<hamiltonian().outerSize(); k++) { 
 	  for (SparseMatrix<double>::InnerIterator it(hamiltonian(),k); it; ++it) {
	     int n = it.row();   // row index
	     int m = it.col();   // col index (here it is equal to k)
	     for (int c = 0; c < SpinTupleInCavity::size(); c++) { 

	       // d rho(n, c)/dt = -i H(n, m) rho(m, c) 

	       coef.push_back( Triplet<complexg>( enc_master(n, c), enc_master(m, c) , - iii * it.value() ) );

	       // d rho(c, m)/dt = i rho(c, n) H(n, m) 
	       coef.push_back( Triplet<complexg>( enc_master(c, m), enc_master(c, n) ,  iii * it.value() ) );

	     }
	  }
       }

       for (int l = 0; l < number_of_lindblad_operators; l++) { 
	  double gamma = std::get<Lrate>(lindblad_operators[l]);
	  if (gamma == 0) continue;

	  SparseMatrix< double > &ladjlop = std::get<Ladjoint_Lop>(lindblad_operators[l]);

	  for (int k=0; k< ladjlop.outerSize(); k++) { 
	     for (SparseMatrix<double>::InnerIterator it(ladjlop, k); it; ++it) {

	        int n = it.row();   // row index
		int m = it.col();   // col index 
		// (L^+ L)_{nm} rho_{mp}
		// rho_{pn} (L^+ L)_{nm} 
		for (int p = 0; p < SpinTupleInCavity::size(); p++) {
		   coef.push_back( Triplet<complexg>( enc_master(n, p), enc_master(m, p) , -0.5 * gamma * it.value() ) );
		   coef.push_back( Triplet<complexg>( enc_master(p, m), enc_master(p, n) , -0.5 * gamma * it.value() ) );
		}
	     }
	  }

	  SparseMatrix< double > &lop = std::get<Lop>(lindblad_operators[l]);
	  SparseMatrix< double > &ladj = std::get<Ladjoint>(lindblad_operators[l]);

	  for (int k1=0; k1< lop.outerSize(); k1++) { 
	     for (SparseMatrix<double>::InnerIterator it1(lop, k1); it1; ++it1) {
	        for (int k2=0; k2< ladj.outerSize(); k2++) { 
		   for (SparseMatrix<double>::InnerIterator it2(ladj, k2); it2; ++it2) {
		      int n = it1.row();   // row index
		      int m = it1.col();   // col index
		      int p = it2.row();   // row index
		      int q = it2.col();   // col index
		      // L_{nm} rho_{mp} Ladj_{pq} 
		      coef.push_back( Triplet<complexg>( enc_master(n, q), enc_master(m, p) , gamma * it1.value() * it2.value() ) );
		   }
		}
	     }
	  }

       }

       Master_full.setFromTriplets(coef.begin(), coef.end());
       return Master_full;

    }

    // updates the master equation looping explicitly through density matrices 
    // it is  constructed explicitely using the liouvillian_sparse function but it is slow
    const SparseMatrix<complexg>&  update_master_from_liouville_explicit(void) { 
       update_hamiltonian();
       Master_full.setZero();
       std::vector< Triplet<complexg> > coef;      

       int matrix_size = SpinTupleInCavity::size();
       for (int n = 0; n < matrix_size; n++) { 
	  for (int m = 0; m < matrix_size; m++) { 
	     SparseMatrix< complexg > rhoit(matrix_size, matrix_size);
	     rhoit.insert(n, m) = 1.0;
	     SparseMatrix< complexg > Lrho(matrix_size, matrix_size);
	     Lrho = liouvillian_sparse( rhoit );

	     for (int k=0; k<Lrho.outerSize(); k++) { 
	       for (SparseMatrix< complexg >::InnerIterator it(Lrho,k); it; ++it) {
		  int nf = it.row();   // row index
		  int mf = it.col();   // col index (here it is equal to k)
		  coef.push_back( Triplet<complexg>( enc_master(nf, mf), enc_master(n, m) , it.value() ) );
	       }
	     }  

	  }
       }

       Master_full.setFromTriplets(coef.begin(), coef.end());
       return Master_full;
    }

    MatrixXcd liouvillian_jacobi(const MatrixXcd &rho) const { 
       return Jacobi_for_Liouville.cwiseProduct( liouvillian(rho) );
    }

    MatrixXcd inverse_jacobi(const MatrixXcd &rho) const { 
       return (rho.array() / Jacobi_for_Liouville.array()).matrix();
    }

    MatrixXcd liouvillian(void) const { 
       return liouvillian(rho_matr);
    }

    MatrixXcd& rho(void) { return rho_matr; }
    complexg rho(int i, int j) { return rho_matr(i, j); } 
    double trId(void) { return real( rho_matr.trace() ); }
    complexg trA(void) { return SpinTupleInCavity::trA(rho_matr); }
    double trN(void) { return SpinTupleInCavity::trN(rho_matr); }
    double trS2(void) { return SpinTupleInCavity::trS2(rho_matr); }
    double trSz(void) { return SpinTupleInCavity::trSz(rho_matr); }
    double trSx(void) { return SpinTupleInCavity::trSx(rho_matr); }
    double trSy(void) { return SpinTupleInCavity::trSy(rho_matr); }
    double trSz(int s) { return SpinTupleInCavity::trSz(s, rho_matr); }
    double trSx(int s) { return SpinTupleInCavity::trSx(s, rho_matr); }
    double trSy(int s) { return SpinTupleInCavity::trSy(s, rho_matr); }
    double trBell_spin_cavity(int s) { return SpinTupleInCavity::trBell_spin_cavity(s, rho_matr); }
    double trBell_spin(int i, int j, const Vector3d &v1, const Vector3d &v2, const Vector3d &v3, const Vector3d &v4) {
       return SpinTupleInCavity::trBell_spin(i, j, v1, v2, v3, v4, rho_matr);
    }
    double trBell_spin_random(int i, int j, int N) { return SpinTupleInCavity::trBell_spin_random(i, j, N, rho_matr); }

    void print_info(void) { 
       SpinTupleInCavity::print_info();
       std::cout << "# gamma " << gamma << std::endl;
       std::cout << "# nav " << nav << std::endl;
	  if (gammaSvec.size() == S.nspins()) {
	     for (int n = 0; n < S.nspins(); n++) { 
	       std::cout << "# gammaS spin" << n << "  " << gammaSvec(n) << std::endl;
	     }
	  } else {
	     std::cout << "# gammaS " << gammaS << std::endl;
	  }
       
       std::cout << "# gammaSup " << gammaSup << std::endl;
    }
    
    double find_rho_eq(int niter=10) { 
       update_hamiltonian();
       SparseMatrix<complexg, ColMajor> Hdiff = Master_full;
       SparseLU< SparseMatrix<complexg, ColMajor>, COLAMDOrdering<int> > solver;
       solver.analyzePattern(Hdiff); 
       solver.factorize(Hdiff); 

       VectorXcd xxx(master_size); 
       VectorXcd yyy(master_size);

       int NS = S.size();
       yyy.setZero();
       for (int s = 0; s < NS; s++) { 
	  yyy( enc_master( enc_nl_s(0, s), enc_nl_s(0, s) ) ) = 1./(double) NS;
       }

       for (int k=0; k<niter; k++) { 
	  xxx = solver.solve(yyy);
	  double n = xxx.norm();
	  xxx /= n;
	  yyy = xxx;
       } 

       yyy = Master_full * xxx;
       double err = yyy.norm();

       int imax = 0;
       for (int i = 1; i < SpinTupleInCavity::size(); i++) { 
 	  if ( norm( xxx( enc_master(i, i) ) ) >  norm( xxx( enc_master(imax, imax) ) ) ) { 
	     imax = i;
	  }
       }       
       
       for (int i = 0; i < SpinTupleInCavity::size(); i++) { 
	  for (int j = 0; j < SpinTupleInCavity::size(); j++) { 
	     rho_matr(i, j) = xxx( enc_master(i, j) ) * exp( -iii * arg( xxx(enc_master(imax, imax) )) );
	  } 
       }

       // for s != 0 all density matricies are traceless
       double sum = real( rho_matr.trace() );
       rho_matr /= sum;
       return err;
    }

    double find_rho_rate(int itermax, bool checkerr) {       
       MatrixXd Hfull(SpinTupleInCavity::size(), SpinTupleInCavity::size());
       Hfull = hamiltonian();
       
       MatrixXd evec;
       VectorXd eval;
       SelfAdjointEigenSolver<MatrixXd> eigensolver( Hfull);
       if (eigensolver.info() != Success) abort();
       eval = eigensolver.eigenvalues();
       evec = eigensolver.eigenvectors();

       int matrix_size = SpinTupleInCavity::size();
       MatrixXd K(matrix_size, matrix_size);
       K.setZero();      

       VectorXcd Ln;
       //
       // K(m, n) = evec.col(m).adjoin() * liouvillian_lindbladt( evec.col(n) * evec.col(n).adjoint() ) * evec.col(m)
       // using the definition of lindbladt operators is more efficient 
       //
       for (int l = 0; l < number_of_lindblad_operators; l++) { 
	  double gamma = std::get<Lrate>(lindblad_operators[l]);
	  if (gamma == 0) continue;

	  SparseMatrix< double > &lop = std::get<Lop>(lindblad_operators[l]);
	  for (int n = 0; n < matrix_size; n++) {
	     Ln = lop * evec.col(n);
	     for (int m = 0; m < matrix_size; m++) {
	        complexg Lmn = evec.col(m).adjoint() * Ln;
		//	        complexg Lnm = Ln.adjoint() * evec.col(m);
		K(m, n) += gamma * norm(Lmn);
	     }
	     K(n, n) -= gamma * Ln.squaredNorm();
	  }
       }

       //       VectorXd v1 = VectorXd::Constant(matrix_size, 1.0);
       //       std::cerr << "# K error " << (K.transpose() * v1).norm() << std::endl;
       MatrixXd ker = K.fullPivLu().kernel();
       
       if (ker.cols() > 1) { 
	  std::cerr << "# kernel > 1 " << ker.cols() << std::endl;
       }

       VectorXd Kvec = ker.col(0) / ker.col(0).sum();

       MatrixXcd rho0( matrix_size, matrix_size );
       // 
       // rho0 solves 
       // L_hamiltonian rho0 = 0 
       // P_{diag} L_lindblad rho0 = 0 
       // Tr [ rho0 ] = 1 
       //
       rho0.setZero();
       for (int i = 0; i < matrix_size; i++) { 
	  rho0 += Kvec(i) * evec.col(i) * evec.col(i).adjoint();
       }

       MatrixXcd Lrho( matrix_size, matrix_size );

       Lrho = evec.adjoint() * liouvillian_lindbladt(rho0) * evec;
       std::cerr << "# error 1 diag L_lindblad rho0 " << Lrho.diagonal().norm() << std::endl;

       rho_matr = rho0;

       // equations to be solved iteratively are 
       //
       // L_lindblad rho_n + L_hamiltonian rho_{n+1} = 0
       // P_{diag} L_lindblad rho_{n+1} = 0 
       // Tr [ rho_{n+1} ] = 0 
       //

       //       MatrixXcd rho_n( matrix_size, matrix_size );
       MatrixXcd rho_next( matrix_size, matrix_size );
       
       rho_next = rho0;

       double err;
       for (int t = 0; t < itermax; t++) {
	  double err_in = liouvillian().norm();
	  std::cerr << "# liouvillian norm step " << t << " : " << err_in << std::endl;
	  Lrho = evec.adjoint() * liouvillian_lindbladt(rho_next) * evec;

	  //	  std::cerr << "# diagonal error " << Lrho.diagonal().norm() << std::endl;
	  for (int i = 0; i < matrix_size; i++) { 
	     for (int j = 0; j < matrix_size; j++) { 
	        if (i != j) { 
		   rho_next(i, j) = -iii * Lrho(i, j) / (eval(i) - eval(j));
		} else { 
		   rho_next(i, j) = 0;
		}
	     }
	  }

	  rho_next = evec * rho_next * evec.adjoint(); 
	  // rho_n solves 
	  // L_lindblad rho_n + L_hamiltonian rho_{n+1} = 0
	  //	  std::cerr << "# error 1 L_lindblad rho_n + L_hamiltonian rho_{n+1} " << (liouvillian_hamiltonian(rho_next) + liouvillian_lindbladt(rho_n)).norm() << std::endl;
	  
	  Lrho = evec.adjoint() * liouvillian_lindbladt(rho_next) * evec;
	  //	  std::cerr << "# error 1 diag L_lindblad rho_{n+1} " << Lrho.diagonal().norm() << std::endl;
	  VectorXd rhodiag;
	  rhodiag = -K.bdcSvd(ComputeThinU | ComputeThinV).solve( Lrho.diagonal().real() );

	  // at this stage rho_next has trace 0 
	  for (int i = 0; i < matrix_size; i++) { 
	     rho_next += rhodiag(i) * evec.col(i) * evec.col(i).adjoint();
	  }       
	  // we add rhodiag which gives a contribution of order \epsilon^n
	  
	  rho_next -= rho_next.trace() * rho0;
	  // where does this come from ??? 
	  // justification seems to be that if rho_next.trace() is of order \epsilon^n
	  // lindbladt( rho_next.trace() * rho0 ) is of order \espilon^{n+1}
	  // so we can use it as a counter term to fix the trace to 1 ...
	  	  
	  //	  std::cerr << "# error 2 L_lindblad rho_n + L_hamiltonian rho_{n+1} " << (liouvillian_hamiltonian(rho_next) + liouvillian_lindbladt(rho_n)).norm() << std::endl;

	  Lrho = evec.adjoint() * liouvillian_lindbladt(rho_next) * evec;

	  //	  std::cerr << "# error 2 diag L_lindblad rho_{n+1} " << Lrho.diagonal().norm() << std::endl;
	  
	  double err_out = liouvillian(rho_matr + rho_next).norm();
	  
	  if (err_out < err_in || !checkerr) { 
	    rho_matr += rho_next;
	    //	    rho_n = rho_next;
	    err = err_out;
	  } else { 
	    err = err_in;
	    break;
	  }
       }
       return err;
    }



    double find_rho_rate2(int itermax, bool checkerr) {       
       MatrixXd Hfull(SpinTupleInCavity::size(), SpinTupleInCavity::size());
       Hfull = hamiltonian();
       
       MatrixXd evec;
       VectorXd eval;
       SelfAdjointEigenSolver<MatrixXd> eigensolver( Hfull);
       if (eigensolver.info() != Success) abort();
       eval = eigensolver.eigenvalues();
       evec = eigensolver.eigenvectors();

       int matrix_size = SpinTupleInCavity::size();
       MatrixXd K(matrix_size, matrix_size);
       MatrixXcd Kc(matrix_size, matrix_size);
       K.setZero();      
       Kc.setZero();
       
       VectorXcd Ln;
       VectorXcd Lm;
       //
       // Kc(m, n) = evec.col(m).adjoint() * liouvillian_lindbladt( evec.col(m) * evec.col(n).adjoint() ) * evec.col(n)
       // computation using the definition of lindbladt operators
       // use K as a temporary for the calculation of Kc
       for (int l = 0; l < number_of_lindblad_operators; l++) { 
	  double gamma = std::get<Lrate>(lindblad_operators[l]);
	  if (gamma == 0) continue;
	  SparseMatrix< double > &lop = std::get<Lop>(lindblad_operators[l]);
	  K = lop * evec;
	  for (int n = 0; n < matrix_size; n++) {
	     Ln = K.col(n);
	     complexg Lnn = evec.col(n).adjoint() * Ln;
	     complexg L2nn = evec.col(n).adjoint() * lop.adjoint() * Ln;
	     for (int m = 0; m < matrix_size; m++) {
	        Lm = K.col(m);
	        complexg Lmm = evec.col(m).adjoint() * Lm;
		complexg L2mm = evec.col(m).adjoint() * lop.adjoint() * Lm;

		Kc(m, n) += gamma * Lmm * conj(Lnn);
		Kc(m, n) -= gamma * 0.5 * L2nn;
		Kc(m, n) -= gamma * 0.5 * L2mm;
	     }
	  }
       }


       //
       // K(m, n) = evec.col(m).adjoin() * liouvillian_lindbladt( evec.col(n) * evec.col(n).adjoint() ) * evec.col(m)
       // using the definition of lindbladt operators is more efficient 
       //
       K.setZero();      
       for (int l = 0; l < number_of_lindblad_operators; l++) { 
	  double gamma = std::get<Lrate>(lindblad_operators[l]);
	  if (gamma == 0) continue;

	  SparseMatrix< double > &lop = std::get<Lop>(lindblad_operators[l]);
	  for (int n = 0; n < matrix_size; n++) {
	     Ln = lop * evec.col(n);
	     for (int m = 0; m < matrix_size; m++) {
	        complexg Lmn = evec.col(m).adjoint() * Ln;
		//	        complexg Lnm = Ln.adjoint() * evec.col(m);
		K(m, n) += gamma * norm(Lmn);
	     }
	     K(n, n) -= gamma * Ln.squaredNorm();
	  }
       }

       
       //       VectorXd v1 = VectorXd::Constant(matrix_size, 1.0);
       //       std::cerr << "# K error " << (K.transpose() * v1).norm() << std::endl;
       MatrixXd ker = K.fullPivLu().kernel();
       
       if (ker.cols() > 1) { 
	  std::cerr << "# kernel > 1 " << ker.cols() << std::endl;
       }

       VectorXd Kvec = ker.col(0) / ker.col(0).sum();

       MatrixXcd rho0( matrix_size, matrix_size );
       // 
       // rho0 solves 
       // L_hamiltonian rho0 = 0 
       // P_{diag} L_lindblad rho0 = 0 
       // Tr [ rho0 ] = 1 
       //
       rho0.setZero();
       for (int i = 0; i < matrix_size; i++) { 
	  rho0 += Kvec(i) * evec.col(i) * evec.col(i).adjoint();
       }

       MatrixXcd Lrho( matrix_size, matrix_size );

       Lrho = evec.adjoint() * liouvillian_lindbladt(rho0) * evec;
       std::cerr << "# error 1 diag L_lindblad rho0 " << Lrho.diagonal().norm() << std::endl;

       rho_matr = rho0;

       // equations to be solved iteratively are 
       //
       // L_lindblad rho_n + L_hamiltonian rho_{n+1} = 0
       // P_{diag} L_lindblad rho_{n+1} = 0 
       // Tr [ rho_{n+1} ] = 0 
       //

       //       MatrixXcd rho_n( matrix_size, matrix_size );
       MatrixXcd rho_next( matrix_size, matrix_size );
       MatrixXcd rho_next_evecbasis( matrix_size, matrix_size );
       
       rho_next = rho0;

       double err;
       for (int t = 0; t < itermax; t++) {
	  double err_in = liouvillian().norm();
	  std::cerr << "# liouvillian norm step " << t << " : " << err_in << std::endl;
	  Lrho = evec.adjoint() * liouvillian_lindbladt(rho_next) * evec;
	  rho_next_evecbasis = evec.adjoint() * rho_next * evec;
	  
	  //	  std::cerr << "# diagonal error " << Lrho.diagonal().norm() << std::endl;
	  for (int i = 0; i < matrix_size; i++) { 
	     for (int j = 0; j < matrix_size; j++) { 
	        if (i != j) { 
		   rho_next(i, j) = -(Lrho(i, j) - Kc(i, j) * rho_next_evecbasis(i, j) )/ ( -iii * (eval(i) - eval(j)) + Kc(i, j) );
		} else { 
		   rho_next(i, j) = 0;
		}
	     }
	  }

	  rho_next = evec * rho_next * evec.adjoint(); 
	  // rho_n solves 
	  // L_lindblad rho_n + L_hamiltonian rho_{n+1} = 0
	  //	  std::cerr << "# error 1 L_lindblad rho_n + L_hamiltonian rho_{n+1} " << (liouvillian_hamiltonian(rho_next) + liouvillian_lindbladt(rho_n)).norm() << std::endl;
	  
	  Lrho = evec.adjoint() * liouvillian_lindbladt(rho_next) * evec;
	  //	  std::cerr << "# error 1 diag L_lindblad rho_{n+1} " << Lrho.diagonal().norm() << std::endl;
	  VectorXd rhodiag;
	  rhodiag = -K.bdcSvd(ComputeThinU | ComputeThinV).solve( Lrho.diagonal().real() );
	  
	  for (int i = 0; i < matrix_size; i++) { 
	     rho_next += rhodiag(i) * evec.col(i) * evec.col(i).adjoint();
	  }       

	  rho_next -= rho_next.trace() * rho0;

	  //	  std::cerr << "# error 2 L_lindblad rho_n + L_hamiltonian rho_{n+1} " << (liouvillian_hamiltonian(rho_next) + liouvillian_lindbladt(rho_n)).norm() << std::endl;

	  Lrho = evec.adjoint() * liouvillian_lindbladt(rho_next) * evec;

	  //	  std::cerr << "# error 2 diag L_lindblad rho_{n+1} " << Lrho.diagonal().norm() << std::endl;
	  
	  double err_out = liouvillian(rho_matr + rho_next).norm();
	  
	  if (err_out < err_in || !checkerr) { 
	    rho_matr += rho_next;
	    //	    rho_n = rho_next;
	    err = err_out;
	  } else { 
	    err = err_in;
	    break;
	  }
       }
       return err;
    }


    double find_rho_rate3(int itermax, bool checkerr) {       
       MatrixXcd Hfull(SpinTupleInCavity::size(), SpinTupleInCavity::size());
       Hfull = hamiltonian();
       
       MatrixXcd evec;
       VectorXd eval;
       SelfAdjointEigenSolver<MatrixXcd> eigensolver( Hfull);
       if (eigensolver.info() != Success) abort();
       eval = eigensolver.eigenvalues();
       evec = eigensolver.eigenvectors();

       int matrix_size = SpinTupleInCavity::size();
       MatrixXd K(matrix_size, matrix_size);
       MatrixXcd Kc(matrix_size, matrix_size);
       MatrixXcd Ktmp(matrix_size, matrix_size);
       MatrixXcd rho0( matrix_size, matrix_size );
       
       VectorXcd Ln;
       VectorXcd Lm;

       double err;
       for (int t = 0; t < itermax; t++) {
	  K.setZero();      
	  Kc.setZero();

	  //
	  // Kc(m, n) = evec.col(m).adjoint() * liouvillian( evec.col(m) * evec.col(n).adjoint() ) * evec.col(n)
	  // computation using the definition of lindbladt operators
	  // use K as a temporary for the calculation of Kc
	  for (int l = 0; l < number_of_lindblad_operators; l++) { 
	     double gamma = std::get<Lrate>(lindblad_operators[l]);
	     if (gamma == 0) continue;
	     SparseMatrix< double > &lop = std::get<Lop>(lindblad_operators[l]);
	     Ktmp = lop * evec;
	     for (int n = 0; n < matrix_size; n++) {
	        Ln = Ktmp.col(n);
		complexg Lnn = evec.col(n).adjoint() * Ln;
		complexg L2nn = evec.col(n).adjoint() * lop.adjoint() * Ln;
		for (int m = 0; m < matrix_size; m++) {
		   if (m == n) continue;

		   Lm = Ktmp.col(m);
		   complexg Lmm = evec.col(m).adjoint() * Lm;
		   complexg L2mm = evec.col(m).adjoint() * lop.adjoint() * Lm;
		   
		   Kc(m, n) += gamma * Lmm * conj(Lnn);
		   Kc(m, n) -= gamma * 0.5 * L2nn;
		   Kc(m, n) -= gamma * 0.5 * L2mm;
		}
	     }
	  }
	  
	  //
	  // K(m, n) = evec.col(m).adjoin() * liouvillian_lindbladt( evec.col(n) * evec.col(n).adjoint() ) * evec.col(m)
	  // using the definition of lindbladt operators is more efficient 
	  //
	  K.setZero();      
	  for (int l = 0; l < number_of_lindblad_operators; l++) { 
	     double gamma = std::get<Lrate>(lindblad_operators[l]);
	     if (gamma == 0) continue;

	     SparseMatrix< double > &lop = std::get<Lop>(lindblad_operators[l]);
	     for (int n = 0; n < matrix_size; n++) {
	       Ln = lop * evec.col(n);
	       for (int m = 0; m < matrix_size; m++) {
		 complexg Lmn = evec.col(m).adjoint() * Ln;
		 //	        complexg Lnm = Ln.adjoint() * evec.col(m);
		 K(m, n) += gamma * norm(Lmn);
	       }
	       K(n, n) -= gamma * Ln.squaredNorm();
	     }
	  }
	  
	  std::cerr << "# here " << std::endl;
	  MatrixXd ker = K.fullPivLu().kernel();
	  
	  if (ker.cols() > 1) { 
	    std::cerr << "# kernel > 1 " << ker.cols() << std::endl;
	  }
	  
	  VectorXd Kvec = ker.col(0) / ker.col(0).sum();

	  // 
	  // rho0 solves 
	  // L_hamiltonian rho0 = 0 
	  // P_{diag} L_lindblad rho0 = 0 
	  // Tr [ rho0 ] = 1 
	  //
	  rho0.setZero();
	  for (int i = 0; i < matrix_size; i++) { 
	    rho0 += Kvec(i) * evec.col(i) * evec.col(i).adjoint();
	  }

	  MatrixXcd Lrho( matrix_size, matrix_size );

	  Lrho = evec.adjoint() * liouvillian(rho0) * evec;
	  std::cerr << "# error 1 diag L_lindblad rho0 " << Lrho.diagonal().norm() << std::endl;

	  rho_matr = rho0;

	  MatrixXcd rho_next( matrix_size, matrix_size );
	  double err_in = liouvillian().norm();
	  std::cerr << "# liouvillian norm step " << t << " : " << err_in << std::endl;
	  
	  //	  std::cerr << "# diagonal error " << Lrho.diagonal().norm() << std::endl;
	  for (int i = 0; i < matrix_size; i++) { 
	     for (int j = 0; j < matrix_size; j++) { 
	        if (i != j) { 
		   rho_next(i, j) = -Lrho(i, j)/ ( -iii * (eval(i) - eval(j)) + Kc(i, j) );
		} else { 
		   rho_next(i, j) = 0;
		}
	     }
	  }

	  rho_next = evec * rho_next * evec.adjoint(); 
	  // rho_n solves 
	  // L_lindblad rho_n + L_hamiltonian rho_{n+1} = 0
	  //	  std::cerr << "# error 1 L_lindblad rho_n + L_hamiltonian rho_{n+1} " << (liouvillian_hamiltonian(rho_next) + liouvillian_lindbladt(rho_n)).norm() << std::endl;
	  
	  Lrho = evec.adjoint() * liouvillian(rho_next) * evec;
	  //	  std::cerr << "# error 1 diag L_lindblad rho_{n+1} " << Lrho.diagonal().norm() << std::endl;
	  VectorXd rhodiag;
	  rhodiag = -K.bdcSvd(ComputeThinU | ComputeThinV).solve( Lrho.diagonal().real() );
	  
	  for (int i = 0; i < matrix_size; i++) { 
	     rho_next += rhodiag(i) * evec.col(i) * evec.col(i).adjoint();
	  }       
	  rho_matr = rho0 + rho_next - rho_next.trace() * rho0;
	  double err_out = liouvillian(rho_matr).norm();
	  std::cerr << "# liouvillian norm step " << t << " : " << err_out << std::endl;
	  
	  if (err_out < err_in || !checkerr) { 
	    rho_matr += rho_next;
	    //	    rho_n = rho_next;
	    err = err_out;
	  } else { 
	    err = err_in;
	    break;
	  }

	  SelfAdjointEigenSolver<MatrixXcd> esolver( rho_matr );
	  if (esolver.info() != Success) abort();
	  evec = esolver.eigenvectors();
	  for (int n = 0; n < matrix_size; n++) {
	     complexg Hnn = evec.col(n).adjoint() * Hfull * evec.col(n);
	     eval(n) = real( Hnn );
	  }
	   
       }
       return err;
    }

  

  
    double find_rho_approx(int itermax, bool checkerr) { return find_rho_rate(itermax, checkerr); } 
    double find_rho_approx(void) { return find_rho_approx(10, true); } 
  
    typedef runge_kutta_dopri5<MatrixXcd, double, MatrixXcd, double,vector_space_algebra> Master_stepper;

    double find_rho_relaxation(void) { 
       int matrix_size = SpinTupleInCavity::size();
       Master_stepper master_stepper;
       SpinTupleInCavityMasterDt sict(*this);
       MatrixXcd rhot = rho_matr;
       MatrixXcd rho0( matrix_size, matrix_size );
       double DT = 1.0;


       double err;
       double scale_up = 1.3;
       double scale_dn = 2.0;
       for (int n = 0; n < 1000; n++) { 
	  rho0 = rhot;
	  double err_in = liouvillian(rho0).norm();
	  double dt = 0.1 * DT;
	  integrate_adaptive(master_stepper, sict, rhot, 0., DT, dt);
	  double err_out = liouvillian(rhot).norm();
	  std::cerr << "# find_rho_relaxation " << err_in << "     " << err_out << std::endl;
	  if (err_out < err_in) { 
	     DT *= scale_up;
	     err = err_out;
	  } else { 
	     rhot = rho0;
	     DT /= scale_dn;
	     err = err_in;
	  }
	  if (err < 1e-4) { 
	     break;
	  }
       }
       return err;
    }    

    double find_rho_approx_nondiag(int itermax, double errmax) {       
       double k = 0;
       int matrix_size = SpinTupleInCavity::size();
       MatrixXcd iHeff(SpinTupleInCavity::size(), SpinTupleInCavity::size());
       iHeff = iii * hamiltonian();
       for (int n = 0; n < number_of_lindblad_operators; n++) { 
	  iHeff += 0.5 * std::get<Lrate>(lindblad_operators[n]) * std::get<Ladjoint_Lop>(lindblad_operators[n]);
       }
       iHeff += (k/2.) * MatrixXd::Identity(matrix_size, matrix_size);

       ComplexSchur<MatrixXcd> schur_iHadj(matrix_size);
       schur_iHadj.compute(iHeff.adjoint());
       
       //
       // -iii Heff rho + iii rho Heff^+ - k rho + ( \sum_n L_n \rho L_n^+ + k \rho) = 0 
       // [ iii Heff . - iii . Heff^+ + k ] rho = ( \sum_n L_n . L_n^+ + k ) \rho 
       // 
       // \rho_{n+1} = [ iii Heff . - iii . Heff^+ + k ]^{-1}  ( \sum_n L_n . L_n^+ + k ) \rho_n
       // 


       MatrixXcd rho_next(matrix_size, matrix_size);
       MatrixXcd rho_n(matrix_size, matrix_size);

       rho_n = rho_matr;

       //
       // compute ( \sum_n L_n . L_n^+ + k ) \rho_n
       // store result to rho_next 
       //
       double err;
       for (int count = 0; count < itermax; count++) { 
	  rho_next.setZero();
	  for (int l = 0; l < number_of_lindblad_operators; l++) { 
	    double gamma = std::get<Lrate>(lindblad_operators[l]);
	    if (gamma == 0) continue;
	    rho_next += gamma * std::get<Lop>(lindblad_operators[l])*rho_n*std::get<Ladjoint>(lindblad_operators[l]);
	  }
	  rho_next += k * rho_n;
	  
	  // solve 
	  // [ iii Heff . - iii . Heff^+ + k ] rho_n = rho_next 
	  sylvester_adjoint( iHeff, schur_iHadj, rho_next, rho_n );

	  rho_n /= rho_n.trace();

	  err = liouvillian(rho_n).norm();
	  if (!(count % 20)) { 
	     std::cerr << "# Sylvester error " << err << " step " << count << std::endl;
	  }
	  if (err < errmax) { break; }
       }
       rho_matr = rho_n;
       std::cerr << "# Sylvester final error " << err << std::endl;
       return err;
    }
};

void SpinTupleInCavityMasterDt::operator() (const MatrixXcd &x , MatrixXcd &dxdt, double t) { 
   dxdt = spin_in_cavity.liouvillian(x);
}


int main_test(int argc, char **argv) { 
    double errmax = 1e-6;
    const static int NL = 70;
    std::vector<int> spin_size = { 2, 2 };
    SpinTupleInCavityMaster sic(NL, spin_size);
    std::vector<double> Delta = { -0.2, 0.2, 0.25 };
    bool first = true;
    MatrixXcd rho_prev;
    for (int s = 0; s < sic.nspins(); s++) { 
       std::cout << "# Delta" << s << " " << Delta[s] << std::endl;
    }
    //    for (double domega = -0.1; domega <= 0.1; domega += 0.001) {
    for (double domega = -0.2; domega <= 0.2; domega += 0.0017) {
       for (int s = 0; s < sic.nspins(); s++) { 
	  sic.Bz[s] = Delta[s] - domega;
	  sic.Omega[s]  = 0.1;
       }
       sic.omegaC = -domega;
       sic.gamma = 0.02;
       sic.nav = 0;
       sic.gammaS = 0.02;
       double n = 20.0;
       sic.F = sic.gamma * sqrt(n) / 2.0;
       sic.gammaSup = 0.0;

       // sic.update_master();
       sic.update_master_from_liouville();
       double err, err_prev;
       //       err = sic.find_rho_eq();

       if (!first) { 
          rho_prev = sic.rho();
          err_prev = sic.liouvillian().norm();
       }

       err = sic.find_rho_approx();
       if (!first && err_prev < err) { 
          sic.rho() = rho_prev;
          std::cerr << "# reset to rho_prev with error " << err_prev << " instead of " << err << std::endl;
       }
       if (err > errmax) { 
	  err = sic.find_rho_approx_nondiag(1000, errmax);
       }

       if (first) { sic.print_info(); first = false; }

       complexg alpha = sic.trA();
       std::cout << "# " << domega << "  ";
       std::cout << sic.trN() << "  " << sic.trSx() << "  " << sic.trSy() << "   " << sic.trSz() << "  " << sic.trS2() << "  ";
       std::cout << real(alpha) << "   " << imag(alpha) << "   ";
       std::cout << err << "  ";
       for (int s = 0; s < sic.nspins(); s++) { 
	  std::cout << sic.trSx(s) << "  " << sic.trSy(s) << "   " << sic.trSz(s) << "  ";
       }
       std::cout << " % " << std::endl;

       for (int n = 0; n < sic.os.size(); n++) { 
	  std::cout << domega << "   "  << n << "   "; 
	  for (int m = 0; m < sic.S.size(); m++) { 
	     int nm = sic.enc_nl_s(n, m);
	     std::cout << norm( sic.rho(nm, nm) ) << "   "; 
	  }
	  std::cout << std::endl;
       }
       std::cout << std::endl;
    }
    return 0;
}

#endif
