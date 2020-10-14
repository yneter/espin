#ifndef ESPINTUPLE 
#define ESPINTUPLE

#include "espindense.cc"

typedef std::unique_ptr<GenericSpinBase> SpinBasePtr;
// typedef std::shared_ptr<GenericSpinBase> SpinBasePtr;
namespace SpinTupleAux { 
  template <typename T> constexpr int find_matrix_size(void) { return T::matrix_size; }

  template<typename T, typename... Tp> 
  constexpr inline typename std::enable_if< sizeof...(Tp) >= 1, int>::type find_matrix_size(void) { 
    return T::matrix_size * find_matrix_size<Tp...>();
  }

  template <typename T> void fill_S(std::vector< SpinBasePtr > &S) { S.push_back(SpinBasePtr(new T)); }

  template<typename T, typename... Tp> 
  constexpr inline typename std::enable_if< sizeof...(Tp) >= 1, void>::type fill_S(std::vector< SpinBasePtr > &S) { 
    S.push_back(SpinBasePtr(new T));
    fill_S<Tp...>(S);
  }

  template<int I, int J, typename T, typename... Tp> 
  constexpr inline typename std::enable_if< I == J , int>::type find_item_matrix_size(void) { 
    return T::matrix_size;
  }

  template<int I, int J, typename T, typename... Tp> 
  constexpr inline typename std::enable_if< J < I, int>::type find_item_matrix_size(void) { 
    return find_item_matrix_size<I, J+1, Tp...>();
  }

  template <int I, typename... Tp>  constexpr inline int find_item_matrix_size(void) { 
    return find_item_matrix_size<I,0,Tp...>();
  }

  template<int I, int J, typename T, typename... Tp> 
  constexpr inline typename std::enable_if< I == J , int>::type find_left_matrix_size(void) { 
    return 1;
  }

  template<int I, int J, typename T, typename... Tp> 
  constexpr inline typename std::enable_if< J < I, int>::type find_left_matrix_size(void) { 
    return T::matrix_size * find_left_matrix_size<I, J+1, Tp...>();
  }

  template <int I, typename... Tp>  constexpr inline int find_left_matrix_size(void) { 
    return find_left_matrix_size<I,0,Tp...>();
  }

}






template <typename... Tp> struct SpinTuple { 
    enum { matrix_size = SpinTupleAux::find_matrix_size<Tp...>() };
    enum { spin_number = sizeof...(Tp) };
    typedef Matrix<complexg, matrix_size, matrix_size> SpinMatrix;
    typedef Matrix<complexg, matrix_size, 1> SpinVector;
    typedef Matrix<double, matrix_size, 1> SpinVectorReal;
    typedef Matrix<double, 3, spin_number> SpinPositionMatrix;
    typedef Matrix<double, spin_number, spin_number> SpinExchangeMatrix;

private:
    SpinMatrix Hfull;
    std::vector< SpinBasePtr > Svec;
    bool dipole_dipole_enabled;
    double Gamma_dip;
    bool exchange_enabled;
    SpinPositionMatrix spin_positions;
    SpinExchangeMatrix spin_exchange_matrix;

public:

    SpinTuple() 
    {
        SpinTupleAux::fill_S<Tp...>(Svec);
	dipole_dipole_enabled = false;
	exchange_enabled = false;
    }

private :
    template <int I> SpinMatrix make_matrix_Hi(TransferSpinMatrix Hi) { 
       constexpr static const int size_left = SpinTupleAux::find_left_matrix_size<I, Tp...>();
       constexpr static const int size_i = SpinTupleAux::find_item_matrix_size<I, Tp...>();
       constexpr static const int size_right = matrix_size / (size_left * size_i);
       typedef Matrix<complexg, size_left, size_left> MatrixLeft;
       typedef Matrix<complexg, size_i, size_i> MatrixItem;
       typedef Matrix<complexg, size_right, size_right> MatrixRight;
       return kroneckerProduct(MatrixLeft::Identity(), kroneckerProduct(Map< MatrixItem > (Hi), MatrixRight::Identity())).eval();
    } 


    template <int I, int J> inline typename std::enable_if< I < J, SpinMatrix>::type make_matrix_HiHj(TransferSpinMatrix Hi, TransferSpinMatrix Hj) { 
       constexpr static const int size_I_left = SpinTupleAux::find_left_matrix_size<I, Tp...>();
       constexpr static const int size_I = SpinTupleAux::find_item_matrix_size<I, Tp...>();
       constexpr static const int size_J_left = SpinTupleAux::find_left_matrix_size<J, Tp...>();
       constexpr static const int size_J = SpinTupleAux::find_item_matrix_size<J, Tp...>();
       constexpr static const int size_J_right = matrix_size / (size_J_left * size_J);
       constexpr static const int size_center = size_J_left / (size_I_left * size_I);
       
       typedef Matrix<complexg, size_I_left, size_I_left> MatrixLeft;
       typedef Matrix<complexg, size_I, size_I> MatrixI;
       typedef Matrix<complexg, size_center, size_center> MatrixCenter;
       typedef Matrix<complexg, size_J, size_J> MatrixJ;
       typedef Matrix<complexg, size_J_right, size_J_right> MatrixRight;
       
       return kroneckerProduct(MatrixLeft::Identity(), 
			       kroneckerProduct(Map< MatrixI > (Hi), 
						kroneckerProduct(MatrixCenter::Identity(), 
								 kroneckerProduct(Map< MatrixJ > (Hj), MatrixRight::Identity())
								 ))).eval();
    } 

    template <int I, int J> inline typename std::enable_if< J < I, SpinMatrix>::type make_matrix_HiHj(TransferSpinMatrix Hi, TransferSpinMatrix Hj) { 
        return make_matrix_HiHj<J, I>(Hj, Hi);
    } 


    template <int I> inline typename std::enable_if< I < sizeof...(Tp), void>::type uncoupled_hamiltonian(void) { 
        Hfull += make_matrix_Hi<I>( Svec[I]->hamiltonian_gen() );
	uncoupled_hamiltonian<I+1>();
    }

    template <int I> inline typename std::enable_if< I == sizeof...(Tp), void>::type uncoupled_hamiltonian(void) { }


public : 
    void load_uncoupled_hamiltonian(void) { 
       Hfull = SpinMatrix::Zero();
       uncoupled_hamiltonian<0>(); 
    }


    template <int I, int J> void add_exchange(double Jij) { 
       Hfull += Jij * make_matrix_HiHj<I, J> ( Svec[I]->Sx_gen(), Svec[J]->Sx_gen() );
       Hfull += Jij * make_matrix_HiHj<I, J> ( Svec[I]->Sy_gen(), Svec[J]->Sy_gen() );
       Hfull += Jij * make_matrix_HiHj<I, J> ( Svec[I]->Sz_gen(), Svec[J]->Sz_gen() );
    }

    void add_matrix(const SpinMatrix &M) { 
       Hfull += M;
    }

    // normalizes uvec to 1 
    template <int I, int J> void add_dipole_dipole(double Jij, Vector3d uvec) { 
       constexpr static const int size_I = SpinTupleAux::find_item_matrix_size<I, Tp...>();
       constexpr static const int size_J = SpinTupleAux::find_item_matrix_size<J, Tp...>();
       typedef Matrix<complexg, size_I, size_I> MatrixI;
       typedef Matrix<complexg, size_J, size_J> MatrixJ;
       double unorm = uvec.norm();
       uvec /= unorm;
       MatrixI uSI = uvec(0) * Map< MatrixI > ( Svec[I]->Sx_gen() ) 
	           + uvec(1) * Map< MatrixI > ( Svec[I]->Sy_gen() )  
	           + uvec(2) * Map< MatrixI > ( Svec[I]->Sz_gen() );       

       MatrixJ uSJ = uvec(0) * Map< MatrixJ > ( Svec[J]->Sx_gen() ) 
	           + uvec(1) * Map< MatrixJ > ( Svec[J]->Sy_gen() )  
	           + uvec(2) * Map< MatrixJ > ( Svec[J]->Sz_gen() );       
       add_exchange<I, J>(Jij);
       Hfull -= 3.0 * Jij * make_matrix_HiHj<I, J> ( uSI.data(), uSJ.data() );
    }

    GenericSpinBase &S(int i) { return *Svec[i]; }

    SpinMatrix hamiltonian(void) const { 
       return Hfull;
    }


    SpinVectorReal eval;   // eigenvalues
    SpinMatrix evec; // eigenvectors
    void diag(void) { 
       SelfAdjointEigenSolver<SpinMatrix> eigensolver(Hfull);
       if (eigensolver.info() != Success) abort();
       eval = eigensolver.eigenvalues();
       evec = eigensolver.eigenvectors();
    }

private : 
    template <int I, int J> inline typename std::enable_if< I < J, void>::type add_dipole_dipole_I_less_J(double Jij, Vector3d uvec) {
       //       std::cerr << "# dipole loop " << I << "     " << J << "      " << Jij << std::endl;
       add_dipole_dipole<I,J>(Jij, uvec);
    }

    template <int I, int J> inline typename std::enable_if< I >= J, void>::type add_dipole_dipole_I_less_J(double Jij, Vector3d uvec) {
    }

    template <int K> inline typename std::enable_if< K != 0, void>::type add_dipole_dipole_loop(double Gamma, const SpinPositionMatrix &spin_position_matrix) { 
       constexpr static const int I = K / spin_number;
       constexpr static const int J = K % spin_number;      

       Vector3d r12 = spin_position_matrix.col(I) - spin_position_matrix.col(J);
       Vector3d uvec = r12.normalized();
       double d12 = r12.norm();
       double Jij = Gamma / (d12 * d12 * d12);

       add_dipole_dipole_I_less_J<I,J>(Jij, uvec);
       add_dipole_dipole_loop<K-1>(Gamma, spin_position_matrix);
    }

    template <int K> inline typename std::enable_if< K == 0, void>::type add_dipole_dipole_loop(double Gamma, const SpinPositionMatrix &spin_position_matrix) { }


    template <int I, int J> inline typename std::enable_if< I < J, void>::type add_exchange_I_less_J(double Jex) {
       add_exchange<I, J>(Jex);
    }

    template <int I, int J> inline typename std::enable_if< I >= J, void>::type add_exchange_I_less_J(double Jex) {
    }
    
    template <int K> inline typename std::enable_if< K != 0, void>::type add_exchange_loop(void) { 
       constexpr static const int I = K / spin_number;
       constexpr static const int J = K % spin_number;      
       add_exchange_I_less_J<I,J> ( spin_exchange_matrix(I, J) );
       add_exchange_loop<K-1>();
    }

    template <int K> inline typename std::enable_if< K == 0, void>::type add_exchange_loop(void) { 
    }

public : 

    void dipole_dipole_interaction(double Gamma, const SpinPositionMatrix &spin_positions_matrix) { 
       dipole_dipole_enabled = true;
       Gamma_dip = Gamma;
       spin_positions = spin_positions_matrix;
    }

    void disable_dipole_dipole_interaction(void) { 
       dipole_dipole_enabled = false;
    }

    void exchange_interaction(const SpinExchangeMatrix &exchange_interaction_matrix) { 
       exchange_enabled = true;
       spin_exchange_matrix = exchange_interaction_matrix;
    }

    SpinExchangeMatrix exchange_interaction(void) { 
       return spin_exchange_matrix;
    }

    void disable_exchange_interaction(void) { 
       exchange_enabled = false;
    }

private : 
    template <int I> inline typename std::enable_if< I < sizeof...(Tp), SpinMatrix>::type Sx_loop(void) { 
        return make_matrix_Hi<I>( Svec[I]->Sx_gen() ) + Sx_loop<I+1>();
    }

    template <int I> inline typename std::enable_if< I == sizeof...(Tp), SpinMatrix>::type Sx_loop(void) { 
        return SpinMatrix::Zero();
    }

    template <int I> inline typename std::enable_if< I < sizeof...(Tp), SpinMatrix>::type Sy_loop(void) { 
        return make_matrix_Hi<I>( Svec[I]->Sy_gen() ) + Sy_loop<I+1>();
    }

    template <int I> inline typename std::enable_if< I == sizeof...(Tp), SpinMatrix>::type Sy_loop(void) { 
        return SpinMatrix::Zero();
    }

    template <int I> inline typename std::enable_if< I < sizeof...(Tp), SpinMatrix>::type Sz_loop(void) { 
        return make_matrix_Hi<I>( Svec[I]->Sz_gen() ) + Sz_loop<I+1>();
    }

    template <int I> inline typename std::enable_if< I == sizeof...(Tp), SpinMatrix>::type Sz_loop(void) { 
        return SpinMatrix::Zero();
    }

public :
    int size(void) const { return matrix_size; }
    int nspins(void) const { return spin_number; } 

    SpinMatrix Sx(void) { return Sx_loop<0>(); }
    SpinMatrix Sy(void) { return Sy_loop<0>(); }  
    SpinMatrix Sz(void) { return Sz_loop<0>(); }
    SpinMatrix S2(void) { return Sx()*Sx() + Sy()*Sy() + Sz()*Sz(); }
    SpinMatrix Jproj(double J) { 
       SpinVectorReal Jeval;   // eigenvalues
       SpinMatrix Jevec; // eigenvectors
       SelfAdjointEigenSolver<SpinMatrix> eigensolver( S2() );
       double J2 = J*(J+1.);
       if (eigensolver.info() != Success) abort();
       Jeval = eigensolver.eigenvalues();
       Jevec = eigensolver.eigenvectors();
       //       std::cerr << "# S2 diag " << (S2() - Jevec * Jeval.asDiagonal() * Jevec.adjoint()).norm() << std::endl;
       //       std::cerr << Jeval << std::endl;
       for (int n = 0; n < size(); n++) { 
	  if (fabs(Jeval(n) - J2) < 1e-5) { 
	     Jeval(n) = 1.0;
	  } else { 
	     Jeval(n) = 0.0;
	  }
       }
       return Jevec * Jeval.asDiagonal() * Jevec.adjoint();
    }

    complexg Sx(int n, int m) { 
       return evec.col(n).adjoint() * Sx() * evec.col(m);
    }

    complexg Sy(int n, int m) { 
       return evec.col(n).adjoint() * Sy() * evec.col(m);
    }

    complexg Sz(int n, int m) { 
       return evec.col(n).adjoint() * Sz() * evec.col(m);
    }

    complexg matrix_elem(const SpinMatrix &M, int n, int m) { 
       return evec.col(n).adjoint() * M * evec.col(m);
    }

    void update_hamiltonian(void) { 
       load_uncoupled_hamiltonian();
       if (dipole_dipole_enabled) { 
	  add_dipole_dipole_loop<spin_number * spin_number-1>(Gamma_dip, spin_positions);
       }
       if (exchange_enabled) { 
	  add_exchange_loop<spin_number * spin_number-1>();
       }
    }

};


typedef SpinTuple< TripletSpin, TripletSpin, SpinHalf > HFE_SpinTuple;

struct HFE : public HFE_SpinTuple { 
   typedef HFE_SpinTuple::SpinMatrix  SpinMatrix;

   double J;
   double dJ;
   double t;
   double Jdip;
   Vector3d r12;

   HFE() { 
      J = dJ = t = Jdip = 0.0;
      r12 << 0, 0, 1;
   }

   SpinMatrix update_hamiltonian(void) { 
       S(2).B << t, 0, 0;
       load_uncoupled_hamiltonian();
       add_exchange<0,1>(J);
       add_dipole_dipole<0,1>(Jdip, r12);
       add_matrix( dJ * kroneckerProduct( TripletPair::exchange_matrix(), SpinHalf::sz ).eval() );
       return hamiltonian();
   }

    static const SpinMatrix singlet_projector(void) { 
       return kroneckerProduct( TripletPair::singlet_projector(), SpinHalf::id ).eval();
    }

    double singlet_content(int i) {
      //       Matrix<complexg, 1, 1> iProj = evec.block(0, i, matrix_size, 1).adjoint() * singlet_projector() * evec.block(0, i, matrix_size, 1);
       Matrix<complexg, 1, 1> iProj = evec.col(i).adjoint() * singlet_projector() * evec.col(i);
       return real(iProj(0,0));
    }


    double PLa(const SpinMatrix &rho) {
       static const SpinMatrix Pa = kroneckerProduct( TripletPair::singlet_projector(), (SpinHalf::SpinMatrix() << 1.0, 0.0, 0.0, 0.0).finished() ).eval();
       Matrix<complexg, 1, 1> trace;
       trace << 0.0;
       for (int i = 0; i < rho.rows(); i++) { 
	 trace += Pa.row(i) * rho.col(i);
       }
       return real(trace(0,0));
    }

    double PLb(const SpinMatrix &rho) {
       static const SpinMatrix Pb = kroneckerProduct( TripletPair::singlet_projector(), (SpinHalf::SpinMatrix() << 0.0, 0.0, 0.0, 1.0).finished() ).eval();
       Matrix<complexg, 1, 1> trace;
       trace << 0.0;
       for (int i = 0; i < rho.rows(); i++) { 
	 trace += Pb.row(i) * rho.col(i);
       }
       return real(trace(0,0));
    }

};


#endif 
