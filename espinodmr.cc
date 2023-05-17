#include "espindense.cc"
#include <Eigen/Sparse>


#ifndef ESPINODMR
#define ESPINODMR


/*
 * ESR_Signal
 *
 * Output : Computes ESR signal
 *
 * Input : spins, a reference on SpinSystem object
 * SpinSystem should define 
 * spins.matrix_size
 * spins.evec
 * spins.eval
 *
 * spins.Sx
 * spins.Sy
 * spins.Sz
 *
 * no longer needed : spins.Bac_field_basis_matrix()
 */ 
template <class SpinSystem> class ESR_Signal { 
public :
    typedef Matrix<double, SpinSystem::matrix_size, 1> SpinVectorReal;
    typedef Matrix<complexg, SpinSystem::matrix_size, 1> SpinVector;
    typedef Matrix<complexg, SpinSystem::matrix_size, SpinSystem::matrix_size> SpinMatrix;
protected:
    SpinMatrix Vx;
    SpinSystem &spins;
public:
    SpinVectorReal rho0;
    double gamma;
    double gamma_diag;
    Vector3d Bac;  // magnetic field in laboratory frame 


    ESR_Signal(SpinSystem &spin_system) : spins(spin_system){
       Bac << 1.0, 0.0, 0.0;
    }

    void update_from_spin_hamiltonian(void) { 
       Vx = spins.evec.adjoint() * (Bac(0) * spins.Sx() + Bac(1) * spins.Sy() + Bac(2) * spins.Sz() )* spins.evec; 
      //	Vx = spins.Bac_field_basis_matrix(); 
    }

    double omega_nm(int n, int m) { 
       return spins.eval(n) - spins.eval(m);
    }

    void load_rho0_thermal(double Temp) { 
       for (int i = 0; i < spins.matrix_size ; i++) { 
	  rho0(i) = exp(- spins.eval(i) / Temp);
       }
       double t = rho0.sum();
       rho0 /= t;
    }

    void load_rho0_from_state_projections(const SpinVector &state) { 
       for (int i = 0; i < spins.matrix_size ; i++) { 
	  complexg ci = state.adjoint() * spins.evec.col(i);
	  rho0(i) = norm(ci);
       }
    }

    void load_rho0_from_projector(const SpinMatrix &proj) { 
       for (int i = 0; i < spins.matrix_size ; i++) { 
	  complexg ci = spins.evec.col(i).adjoint() * proj * spins.evec.col(i);
	  rho0(i) = real(ci);
       }
    }

    void load_rho0(const std::vector<double> &values) { 
       for (int i = 0; i < spins.matrix_size ; i++) { 
	  rho0(i) = values[i];
       }
       double t = rho0.sum();
       rho0 /= t;
    }    

    void load_rho0(const double values[]) { 
       double t = 0.0;
       for (int i = 0; i < spins.matrix_size ; i++) { 
	  rho0(i) = values[i];
	  t += values[i];
       }
       rho0 /= t;
    }    


    complexg chi1(double omega) { 
       complexg c1 = 0.0;
       for (int m = 0; m < spins.matrix_size ; m++) { 
	  for (int n = 0; n < spins.matrix_size ; n++) { 
	     // 
	     // the contribution to chi1 vanishes for n == m, whether gamma is the same for diagonal and non diagonal elements is not relvant here 
	     // 
	     c1 -= (rho0(m) - rho0(n)) * norm(Vx(n, m)) / ( omega_nm(n, m) - omega - iii * gamma );
	  }
       }
       return c1;
    }    

    complexg sq(complexg x) { return x*x; } 
  
    complexg dchi1_domega(double omega) { 
       complexg c1 = 0.0;
       for (int m = 0; m < spins.matrix_size ; m++) { 
	  for (int n = 0; n < spins.matrix_size ; n++) { 
	     // 
	     // the contribution to chi1 vanishes for n == m, whether gamma is the same for diagonal and non diagonal elements is not relvant here 
	     // 
	     c1 -= (rho0(m) - rho0(n)) * norm(Vx(n, m)) / sq( omega_nm(n, m) - omega - iii * gamma );
	  }
       }
       return c1;
    }    
  
};


/*
 * ODMR_Signal
 *
 * Output : Computes ODMR signal
 *
 * Input : spins, a reference on SpinSystem object
 * SpinSystem should define 
 * for base class ESR_Signal
 * spins.matrix_size
 * spins.evec
 * spins.eval
 * and for ODMR signal 
 * spins.singlet_projector()
 */ 
template <class SpinSystem> class ODMR_Signal : public ESR_Signal<SpinSystem> { 
    typedef Matrix<double, SpinSystem::matrix_size, 1> SpinVectorReal;
    typedef Matrix<complexg, SpinSystem::matrix_size, SpinSystem::matrix_size> SpinMatrix;
    SpinMatrix rho2;
    SpinMatrix Sproj_eig_basis;
public :
    ODMR_Signal(SpinSystem &spin_system): ESR_Signal<SpinSystem>(spin_system) { 

    }

    void load_rho0_from_singlet(void) { 
       for (int i = 0; i < this->spins.matrix_size ; i++) { 
	  this->rho0(i) = real(Sproj_eig_basis(i, i));
       }
       double t = this->rho0.sum();
       this->rho0 /= t;
    }    

    void update_from_spin_hamiltonian(void) { 
        Sproj_eig_basis = this->spins.evec.adjoint() * this->spins.singlet_projector() * this->spins.evec;
	this->ESR_Signal<SpinSystem>::update_from_spin_hamiltonian();
    }

    void update_from_spin_hamiltonian_local_basis(void) { 
        Sproj_eig_basis = this->spins.singlet_projector();
	this->ESR_Signal<SpinSystem>::update_from_spin_hamiltonian();
    }


    // explicit calculation of rho2 - close to analytical formula but slow
    void find_rho2_explicit(double omega) { 
       for (int m = 0; m < this->spins.matrix_size ; m++) { 
	  for (int n = 0; n < this->spins.matrix_size ; n++) { 
	     complexg rrr = 0.0;
	     for (int nu = 0; nu < this->spins.matrix_size ; nu++) { 
	        for (int p = -1; p <= 1; p += 2) { 
		  // Vtmp(nu, m) = (rho0(m) - rho0(nu)) * V(nu, m) / ( omega_nm(nu, m) - omega * (double) p - iii * gamma )
		   rrr += this->Vx(n, nu) * (this->rho0(m) - this->rho0(nu)) * this->Vx(nu, m) / ( this->omega_nm(nu, m) - omega * (double) p - iii * this->gamma );
		  //  nu->n and m->nu : Vtmp(n, nu)  
		   rrr -= ((this->rho0(nu) - this->rho0(n)) * this->Vx(n, nu) / ( this->omega_nm(n, nu) - omega * (double) p - iii * this->gamma )) * this->Vx(nu, m);
		}
	     }
	     // relaxation may be different for diaganonal and non diagonal terms
	     double gamma_nm = (n == m) ? this->gamma_diag : this->gamma;
	     rho2(n, m) = rrr / ( this->omega_nm(n, m) - iii * gamma_nm );
	  }
       }

    }


    // optimized calculation of rho2
    void find_rho2(double omega) { 
       SpinMatrix Vtmp = SpinMatrix::Zero();
       for (int m = 0; m < this->spins.matrix_size ; m++) { 
	  for (int nu = 0; nu < this->spins.matrix_size ; nu++) { 
	     for (int p = -1; p <= 1; p += 2) { 
	        Vtmp(nu, m) += (this->rho0(m) - this->rho0(nu)) * this->Vx(nu, m) / (this->omega_nm(nu, m) - omega * (double) p - iii * this->gamma);
	     }
	  }
       }      
       rho2 = this->Vx * Vtmp - Vtmp * this->Vx;
       for (int m = 0; m < this->spins.matrix_size ; m++) { 
	  for (int n = 0; n < this->spins.matrix_size ; n++) { 
	     // relaxation may be different for diaganonal and non diagonal terms
	     double gamma_nm = (n == m) ? this->gamma_diag : this->gamma;
	     rho2(n, m) /= ( this->omega_nm(n, m) - iii * gamma_nm );
	  }
       }
    }


    double odmr(double omega) { 
       double odmr_amp = 0.0;
       find_rho2(omega);
       
       for (int m = 0; m < this->spins.matrix_size ; m++) { 
	  for (int n = 0; n < this->spins.matrix_size ; n++) { 
	     odmr_amp += real( rho2(m , n) * Sproj_eig_basis(n, m) );
	  }
       }

       return odmr_amp;
    }
};


/*
 * Merrifield
 *
 * Computes steady state spin density matrix from the master equation 
 * using eigen's matrix free iterative solvers.
 * Matrix free solver code borrowed from  
 * https://eigen.tuxfamily.org/dox/group__MatrixfreeSolverExample.html
 *
 * Input : spins, a reference on SpinSystem object
 * SpinSystem should define 
 * spins.matrix_size
 * spins.singlet_projector()
 * spins.hamiltonian()
 *
 */ 
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

template<typename SpinSystem> class Merrifield;

namespace Eigen {
namespace internal {
  // make Merrifield look like Matrix<complexg, SpinSystem::matrix_size^2, SpinSystem::matrix_size^2> 
  template<class SpinSystem>
  //  struct traits< Merrifield<SpinSystem> > : public Eigen::internal::traits< Matrix<complexg, SpinSystem::matrix_size*SpinSystem::matrix_size, SpinSystem::matrix_size*SpinSystem::matrix_size> >
  struct traits< Merrifield<SpinSystem> > : public Eigen::internal::traits<Eigen::SparseMatrix<complexg> >
  {};

}
}

template <class SpinSystem> class Merrifield :public Eigen::EigenBase< Merrifield<SpinSystem> >  { 
    SpinSystem &spins;
    bool rho_initialized;
public:
    double gammaS;
    double gamma;
    double alphaS;
    double alpha;

    // Required typedefs, constants, and method:
    typedef complexg Scalar;
    typedef double RealScalar;
    typedef int StorageIndex;

    enum { master_size = SpinSystem::matrix_size * SpinSystem::matrix_size };
    enum { matrix_size = SpinSystem::matrix_size };
    typedef Matrix<complexg, SpinSystem::matrix_size, SpinSystem::matrix_size> SpinMatrix;
    typedef Matrix<complexg, master_size, 1> SpinMatrixVecForm;
    SpinMatrix rho;
    SpinMatrix Ps;


    enum {
       ColsAtCompileTime = Eigen::Dynamic,
       MaxColsAtCompileTime = Eigen::Dynamic,
       IsRowMajor = false
    };
  
    Index rows() const { return SpinSystem::matrix_size * SpinSystem::matrix_size; }
    Index cols() const { return SpinSystem::matrix_size * SpinSystem::matrix_size; } 

    typedef Matrix<complexg, master_size, master_size> LiouvilleMatrix;
    int enc_master(int n, int m) { 
       return m * SpinSystem::matrix_size + n;
    }
    LiouvilleMatrix liouville_full_matrix;    
  //    SparseMatrix<complexg>  liouville_sparse_matrix;    

    template<typename Rhs>
    Eigen::Product<Merrifield<SpinSystem>,Rhs,Eigen::AliasFreeProduct> 
    operator*(const Eigen::MatrixBase<Rhs>& x) const {
       return Eigen::Product<Merrifield<SpinSystem>,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
    }


    Merrifield(SpinSystem &spin_system) : spins(spin_system) {
       Ps = spins.singlet_projector();
       rho_initialized = false;       
       //       liouville_sparse_matrix.resize(master_size, master_size);
    }
private : 
    double trace_rho_Ps(const SpinMatrix &mrho) { 
	Matrix<complexg, 1, 1> sum;
	for (int i = 0; i < SpinSystem::matrix_size; i++) {
 	   sum += Ps.row(i) * mrho.col(i);
	}
	return real(sum(0));
    }



public  : 

    SpinMatrix Liouvillian(const SpinMatrix &rho) const { 
      return -iii * ( spins.hamiltonian() * rho - rho * spins.hamiltonian() )
	- gamma * rho
	- 0.5 * gammaS * (Ps * rho + rho * Ps);
    }


    SpinMatrix map_to_mat(SpinMatrixVecForm vec) const { 
        return Map< SpinMatrix >(vec.data());
    }

    SpinMatrixVecForm map_to_vec(SpinMatrix mat) const { 
        return Map< SpinMatrixVecForm > (mat.data());
    }

    SpinMatrixVecForm Ps_to_vec(void) { 
        return map_to_vec(Ps);
    }

    void find_rho(bool use_previous_as_guess = true) { 
      //        Eigen::BiCGSTAB< Merrifield<SpinSystem>, Eigen::DiagonalPreconditioner<complexg> > bicg;
        Eigen::BiCGSTAB< Merrifield<SpinSystem> , Eigen::IdentityPreconditioner > bicg;
      //        Eigen::ConjugateGradient< Merrifield<SpinSystem> , Lower|Upper, Eigen::IdentityPreconditioner > bicg;
	bicg.compute(*this);
	SpinMatrixVecForm x;
	SpinMatrixVecForm y = -Ps_to_vec();
	if (!rho_initialized || !use_previous_as_guess) { 
	   x = bicg.solve(y);    
	} else { 
	   x = bicg.solveWithGuess(y, map_to_vec(rho));
	}
	//	std::cout << "BiCGSTAB: #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << std::endl;
	std::cerr << "#iterations:     " << bicg.iterations() << std::endl;
	std::cerr << "estimated error: " << bicg.error()      << std::endl;

	rho = map_to_mat(x);
	rho_initialized = true;
    }


    void find_rho(SpinMatrix guess) { 
      //        Eigen::BiCGSTAB< Merrifield<SpinSystem>, Eigen::DiagonalPreconditioner<complexg> > bicg;
        Eigen::BiCGSTAB< Merrifield<SpinSystem> , Eigen::IdentityPreconditioner > bicg;
	//       Eigen::ConjugateGradient< Merrifield<SpinSystem> , Lower|Upper, Eigen::IdentityPreconditioner > bicg;
	SpinMatrixVecForm x;
	SpinMatrixVecForm y = -Ps_to_vec();
	bicg.compute(*this);
	x = bicg.solveWithGuess(y, map_to_vec(guess));
	//	x = bicg.solve(y);    
	//	std::cout << "BiCGSTAB: #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << std::endl;
	std::cerr << "#iterations:     " << bicg.iterations() << std::endl;
	std::cerr << "#estimated error: " << bicg.error()      << std::endl;
	rho = map_to_mat(x);
	rho_initialized = true;
    }

    const LiouvilleMatrix &update_liouville_matrix(void) { 
        spins.update_hamiltonian();
	liouville_full_matrix.setZero();
	//	liouville_sparse_matrix.setZero();
	//	std::vector< Triplet<complexg> > coef;      
       
	const SpinMatrix &H = spins.hamiltonian();
	for (int n=0; n< matrix_size; n++) { 
	   for (int m=0; m< matrix_size; m++) { 
	      for (int c=0; c< matrix_size; c++) { 
	       // d rho(n, c)/dt = -i H(n, m) rho(m, c) 
		liouville_full_matrix( enc_master(n, c), enc_master(m, c)) -= iii * H(n, m);
	       // d rho(c, m)/dt = i rho(c, n) H(n, m) 
		liouville_full_matrix( enc_master(c, m), enc_master(c, n)) += iii * H(n, m);
		/**
		if ( H(n, m) != complex<double>(0, 0) ) { 
		   coef.push_back( Triplet<complexg>( enc_master(n, c), enc_master(m, c) , - iii * H(n, m) ) );
		   coef.push_back( Triplet<complexg>( enc_master(c, m), enc_master(c, n) , + iii * H(n, m) ) );
		}
		**/
	       // d rho(n, c)/dt = -gammaS PS(n, m) rho(m, c) 
		liouville_full_matrix( enc_master(n, c), enc_master(m, c)) -= 0.5 * gammaS * Ps(n, m);
	       // d rho(c, m)/dt = -gammaS rho(c, n) PS(n, m) 
		liouville_full_matrix( enc_master(c, m), enc_master(c, n)) -= 0.5 * gammaS * Ps(n, m);
		/**
		if (Ps(n, m) != complex<double>(0, 0)) { 
		   coef.push_back( Triplet<complexg>( enc_master(n, c), enc_master(m, c) , - gammaS * Ps(n, m) ) );
		   coef.push_back( Triplet<complexg>( enc_master(c, m), enc_master(c, n) , - gammaS * Ps(n, m) ) );
		}
		**/

	      }	      
	      //  drho(n, m)/dt = -gamma  rho(n, m)
	      liouville_full_matrix( enc_master(n, m), enc_master(n, m)) -= gamma;
	      //	      coef.push_back( Triplet<complexg>( enc_master(n, m), enc_master(n, m) , - gamma ) );
	   }
	}      

	//	liouville_sparse_matrix.setFromTriplets(coef.begin(), coef.end());       
	return liouville_full_matrix;
    }

    void find_rho_exact(bool sparse = false) { 
       SpinMatrixVecForm y;
       y.setZero();
       for (int n = 0; n < matrix_size; n++) y( enc_master(n, n) ) = -alpha;
       y -= alphaS * Ps_to_vec();

       SpinMatrixVecForm x;

       update_liouville_matrix();

       x = liouville_full_matrix.lu().solve(y);

       rho = map_to_mat(x);       
       rho_initialized = true;
    }

    double odmr(double omega, Vector3d Bac) { 
       SpinMatrix Vac = Bac(0) * spins.Sx() + Bac(1) * spins.Sy() + Bac(2) * spins.Sz();
       // solve \partial_t \rho = L \rho - i [ Vac, \rho ] 
       SpinMatrix rho2; 
       SpinMatrixVecForm z1;
       z1 = map_to_vec( -iii * Vac * rho + iii * rho * Vac );
       LiouvilleMatrix Lp, Lm;
       Lp = liouville_full_matrix - iii * omega * LiouvilleMatrix::Identity();
       Lm = liouville_full_matrix + iii * omega * LiouvilleMatrix::Identity();
       SpinMatrixVecForm xp, xm;
       xp = Lp.lu().solve(z1);
       xm = Lm.lu().solve(z1);
       SpinMatrixVecForm zp2, zm2;
       zp2 = map_to_vec( -iii * Vac * map_to_mat(xp) + iii * map_to_mat(xp) * Vac );
       zm2 = map_to_vec( -iii * Vac * map_to_mat(xm) + iii * map_to_mat(xm) * Vac );
       xp = liouville_full_matrix.lu().solve(zp2);
       xm = liouville_full_matrix.lu().solve(zm2);
       rho2 = map_to_mat( xp ) + map_to_mat( xm );
       return trace_rho_Ps(rho2);
    }

    double rho_error(void) { 
       return (Liouvillian(rho) + alphaS * Ps + alpha * SpinMatrix::Identity()).norm();
    }

    double PL(void) { 
       return trace_rho_Ps(rho);
    }

    


};

/*
 * Merrifield from rate equations 
 * requires 
 * 
 * Input : spins, a reference on SpinSystem object
 * SpinSystem should define 
 * spins.matrix_size
 * spins.singlet_projector 
 */
template <class SpinSystem> class MerrifieldRate { 
    SpinSystem &spins;
public:
    typedef typename SpinSystem::SpinMatrix SpinMatrix;
    typedef Matrix<complexg, SpinSystem::matrix_size*SpinSystem::matrix_size, 1> SpinMatrixVecForm;
    SpinMatrix rho;

    MerrifieldRate(SpinSystem &spin_system) : spins(spin_system) {
    }

    double gamma;
    double gammaS;
    
    double singlet_content(int i) {
      //       Matrix<complexg, 1, 1> iProj = evec.block(0, i, matrix_size, 1).adjoint() * singlet_projector() * evec.block(0, i, matrix_size, 1);
        complexg iProj = spins.evec.col(i).adjoint() * spins.singlet_projector() * spins.evec.col(i);
	return real(iProj);
    }

    void find_rho(void) { 
        rho = SpinMatrix::Zero();
	for (int i = 0; i < SpinSystem::matrix_size; i++) {
	   double alpha_i = singlet_content(i);
	   rho(i, i) = alpha_i / (gamma + 2.0 * gammaS * alpha_i);
	}	
	rho = spins.evec * rho * spins.evec.adjoint();
    }

    double PL(int npow = 2) { 
        double sum = 0.0;
	for (int i = 0; i < SpinSystem::matrix_size; i++) {
	   double alpha_n = singlet_content(i);
	   sum += pow(alpha_n, npow) / (gamma + 2.0 * gammaS * alpha_n);
	}
	return sum;
    }
};

namespace Eigen {
namespace internal {
  template<typename Rhs, class SpinSystem>
  struct generic_product_impl<Merrifield<SpinSystem>, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<Merrifield<SpinSystem>,Rhs,generic_product_impl<Merrifield<SpinSystem>,Rhs> >
  {
    typedef typename Product<Merrifield<SpinSystem>,Rhs>::Scalar Scalar;
    template<typename Dest>
    static void scaleAndAddTo(Dest& dst, const Merrifield<SpinSystem>& lhs, const Rhs& rhs, const Scalar& alpha)
    {
      // This method should implement "dst += alpha * lhs * rhs" inplace,
      // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
      assert(alpha==Scalar(1) && "scaling is not implemented");
      typename Merrifield<SpinSystem>::SpinMatrix rho = lhs.map_to_mat(rhs);     
      typename Merrifield<SpinSystem>::SpinMatrix L = lhs.Liouvillian(rho);
      typename Merrifield<SpinSystem>::SpinMatrixVecForm lhs_x_rhs = lhs.map_to_vec(L);
      dst += lhs_x_rhs;
    }
  };
}
}

#endif
