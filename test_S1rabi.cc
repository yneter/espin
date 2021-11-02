#include "espincavity.cc"
#include "espindense.cc"
#include "espintensor.cc"

struct Spin1CavityRabi : public SpinCavityTuple {
   double lambda;
   double sigma;
   double omega;
   double Bac;
   // only a single triplet spin 
   // Ncavity - number of states for oscillator 
   Spin1CavityRabi(int Ncavity) : SpinCavityTuple( std::vector<int>( { 3 } ) ,
						   std::vector<int>( { Ncavity } ) )
   {
      lambda = 0;
      sigma = 0;
   }

   void update_hamiltonian(void) { 
      SpinCavityTuple::update_hamiltonian();
      SpinCavityTuple::hamiltonian() += lambda * kroneckerProduct( S.Sx(0), os.a(0) + os.ap(0) ).eval(); // Sx (a + a^+) 
      SpinCavityTuple::hamiltonian() += sigma * kroneckerProduct( S.Sz2(0), os.a(0) + os.ap(0) ).eval(); // Sz^2 (b + b^+)
   }

   SparseMatrix<complexg> hamiltonian(double time) { 
      return hamiltonian() + Bac * Sx(0) * cos(omega * time);
   }

   SparseMatrix<complexg> hamiltonian(void) { return SpinCavityTuple::hamiltonian(); }

   void operator()(const VectorXcd &x , VectorXcd &dxdt, double t);
};

void Spin1CavityRabi::operator() (const VectorXcd &x , VectorXcd &dxdt, double t) { 
   dxdt = -iii * this->hamiltonian(t) * x;
}


int main(int argc, char **argv)
{
   Spin1CavityRabi s1a(50);
   double Ds = 1;  // D parameter for spin 
   s1a.S.S(0).D = Ds; 
   s1a.S.S(0).B << 0., 0.5, 0.; // Bx, By, Bz for spin 

   SingleSpin1 triplet;
   triplet.S = s1a.S.S(0);
   triplet.update_hamiltonian();
   triplet.diag();
   std::cout << "# triplet eval " << triplet.eval(0) << "  " << triplet.eval(1) << std::endl;
   // in zero magnetic field - spin eigenvalues are split from the ground state by D +/- E 
   s1a.omega = triplet.eval(1) - triplet.eval(0);
   s1a.Bac = 0.02; 
   s1a.update_hamiltonian();
   SparseMatrix< std::complex<double> > H0 = s1a.hamiltonian();

   std::complex< double > sx01 = triplet.evec.col(1).adjoint() * triplet.S.Sx() * triplet.evec.col(0);
   double rabi_frequency = s1a.Bac * abs(sx01);
   std::cout << "# sx01 " << rabi_frequency << std::endl; // Rabi frequency 

   s1a.os.cavity(0).omega_c = atof(argv[1]) * rabi_frequency;    // frequency of first cavity set equal to Rabi frequency    
   s1a.lambda = 0.0;   
   s1a.sigma = 0.1 * s1a.os.cavity(0).omega_c;      // coupling strength 
   s1a.update_hamiltonian();
   
   typedef runge_kutta_dopri5<VectorXcd, double, VectorXcd, double,vector_space_algebra> H_stepper;
   H_stepper hstep;

   VectorXcd psi0(s1a.size());
   psi0 = espin_tensor( VectorXcd () = triplet.evec.col(0) , s1a.os.cavity(0).vac_n(0) );
   VectorXcd psi1(s1a.size());
   psi1 = espin_tensor( VectorXcd () = triplet.evec.col(1) , s1a.os.cavity(0).vac_n(0) );

   
   double DT = 0.1;
   double t = 0;
   while (t <= 5000) { 
      int ndt = 10;
      double dt = DT / (double) ndt;
      integrate_n_steps(hstep, s1a, psi0, t, dt, ndt);
      integrate_n_steps(hstep, s1a, psi1, t, dt, ndt);
      std::cout << t
		<< "   " << real( std::complex<double> ( psi0.adjoint() * H0 * psi0 ) )
		<< "   " << real( std::complex<double> ( psi1.adjoint() * H0 * psi1 ) )
		<< "   " << real( std::complex<double> ( psi1.adjoint() * psi0 )  )
		<< "   " << s1a.trSz(0, psi0) << "    " << s1a.trSz2(0, psi0) << "    "
		<< "   " << s1a.trN(0, psi0) << "   " << psi0.norm()-1. 
		<< "   " << s1a.trN(0, psi1) << "   " << psi1.norm()-1. << std::endl;
      t += DT;
   }
}
