#ifndef ESPINCAVITY
#define ESPINCAVITY

#include "espinsparse.cc"
#include <iostream>
#include <tuple>

class CavitySparse { 
    SparseMatrix<double> a_matrix;
    SparseMatrix<double> ap_matrix;
    SparseMatrix<double> N_matrix;
    SparseMatrix<double> Id_matrix;
public :
    CavitySparse(int N) : 
      a_matrix(N, N),
      ap_matrix(N, N),
      N_matrix(N, N),
      Id_matrix(N, N)
    { 
      std::vector< Triplet<double> > coef;      
       for (int i = 0; i < N-1; i++) { 
	  coef.push_back( Triplet<double>( i, i+1, sqrt(i+1) ) );
       }
       a_matrix.setFromTriplets(coef.begin(), coef.end());

       coef.clear();
       for (int i = 0; i < N-1; i++) { 
	  coef.push_back( Triplet<double>( i+1, i, sqrt(i+1) ) );
       }
       ap_matrix.setFromTriplets(coef.begin(), coef.end());


       coef.clear();
       for (int i = 0; i < N; i++) { 
	  coef.push_back( Triplet<double>( i, i, (double)i ) );
       }
       N_matrix.setFromTriplets(coef.begin(), coef.end());
       
       Id_matrix.setIdentity();
    }

    int size(void) { 
       return Id_matrix.cols();
    }

    const SparseMatrix<double>& a(void) const { 
       return a_matrix;
    }

    const SparseMatrix<double>& ap(void) const { 
       return ap_matrix;
    }

    const SparseMatrix<double>& N(void) const { 
       return N_matrix;
    }

    const SparseMatrix<double>& Id(void) const {
       return Id_matrix;
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
    double re_trace_AB(const SparseMatrix<double>& A, const MatrixXcd &B) { return real( trace_AB(A, B) ); }
    double trN(const MatrixXcd &rho) { return re_trace_AB( kroneckerProduct(os.N(), S.rId()), rho); }
    double trS2(const MatrixXcd &rho) { return re_trace_AB( kroneckerProduct(os.Id(), S.rSx() * S.rSx() - S.iSy() * S.iSy() + S.rSz() * S.rSz()), rho); }
    double trSx(const MatrixXcd &rho) { return re_trace_AB( kroneckerProduct(os.Id(), S.rSx()), rho); }
    double trSy(const MatrixXcd &rho) { return imag ( trace_AB( kroneckerProduct(os.Id(), S.iSy()), rho) ); }
    double trSz(const MatrixXcd &rho) { return re_trace_AB( kroneckerProduct(os.Id(), S.rSz()), rho); }
    double trSx(int i, const MatrixXcd &rho) { return re_trace_AB( kroneckerProduct(os.Id(), S.rSx(i)), rho); }
    double trSy(int i, const MatrixXcd &rho) { return imag ( trace_AB( kroneckerProduct(os.Id(), S.iSy(i)), rho) ); }
    double trSz(int i, const MatrixXcd &rho) { return re_trace_AB( kroneckerProduct(os.Id(), S.rSz(i)), rho); }



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

struct SpinTupleInCavityMasterTimeDerivative { 
    const SpinTupleInCavityMaster &spin_in_cavity;
    void operator()(const MatrixXcd &x , MatrixXcd &dxdt, double t);
  
    SpinTupleInCavityMasterTimeDerivative(const SpinTupleInCavityMaster &spin_in_cavity_ref) : 
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
    double gammaSup;
    double temp;
private:
    enum { Lrate, Lop, Ladjoint, Ladjoint_Lop };
    enum { gamma_os, gamma_os_up, gamma_cavity }; 
    int number_of_lindblad_operators;

    typedef std::tuple< double, SparseMatrix< double >, SparseMatrix< double >, SparseMatrix< double > > LindbladOperator;
    std::vector < LindbladOperator > lindblad_operators;    

    int master_size;
    MatrixXcd rho_matr;
    MatrixXcd Jacobi_for_Liouville;
    SparseMatrix<complexg> Master_full;
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
	  std::get<Lrate>( lindblad_operators[gamma_spin] ) = gammaS;
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
    double trN(void) { return SpinTupleInCavity::trN(rho_matr); }
    double trS2(void) { return SpinTupleInCavity::trS2(rho_matr); }
    double trSz(void) { return SpinTupleInCavity::trSz(rho_matr); }
    double trSx(void) { return SpinTupleInCavity::trSx(rho_matr); }
    double trSy(void) { return SpinTupleInCavity::trSy(rho_matr); }
    double trSz(int s) { return SpinTupleInCavity::trSz(s, rho_matr); }
    double trSx(int s) { return SpinTupleInCavity::trSx(s, rho_matr); }
    double trSy(int s) { return SpinTupleInCavity::trSy(s, rho_matr); }

    void print_info(void) { 
       SpinTupleInCavity::print_info();
       std::cout << "# gamma " << gamma << std::endl;
       std::cout << "# nav " << nav << std::endl;
       std::cout << "# gammaS " << gammaS << std::endl;
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

    double find_rho_approx(void) {       
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
       for (int t = 0; t < 10; t++) {
	  double err_in = liouvillian().norm();
	  std::cerr << "# liouvillian norm step " << t << " : " << err_in << std::endl;
	  Lrho = evec.adjoint() * liouvillian_lindbladt(rho_next) * evec;

	  //	  std::cerr << "# diagonal error " << Lrho.diagonal().norm() << std::endl;
	  for (int i = 0; i < matrix_size; i++) { 
	     for (int j = 0; j < matrix_size; j++) { 
	        if (i != j) { 
		   rho_next(i, j) = -iii * Lrho(i, j) / ( eval(i) - eval(j) );
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
	  
	  if (err_out < err_in) { 
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

    typedef runge_kutta_dopri5<MatrixXcd, double, MatrixXcd, double,vector_space_algebra> Master_stepper;

    double find_rho_relaxation(void) { 
       int matrix_size = SpinTupleInCavity::size();
       Master_stepper master_stepper;
       SpinTupleInCavityMasterTimeDerivative sict(*this);
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

void SpinTupleInCavityMasterTimeDerivative::operator() (const MatrixXcd &x , MatrixXcd &dxdt, double t) { 
   dxdt = spin_in_cavity.liouvillian(x);
}


main() { 
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
    for (double domega = -0.2; domega <= 0.2; domega += 0.0033) {
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
       if (err > 1e-6) { 
	  err = sic.find_rho_approx_nondiag(1000, errmax);
       }

       if (first) { sic.print_info(); first = false; }

       std::cout << "# " << domega << "  ";
       std::cout << sic.trN() << "  " << sic.trSx() << "  " << sic.trSy() << "   " << sic.trSz() << "  " << sic.trS2() << "  ";
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
}

#endif
