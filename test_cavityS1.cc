#include "espincavity.cc"

struct Spin1w2Cavities : public SpinCavityTuple {
   double lambda0;
   double sigma1;
   // only a single triplet spin 
   // Ncavity1 - number of states for oscillator 1 
   // Ncavity2 - number of states for oscillator 2
   Spin1w2Cavities(int Ncavity1, int Ncavity2) : SpinCavityTuple( std::vector<int>( { 3 } ) ,
								  std::vector<int>( { Ncavity1, Ncavity2 } ) )
   {
   }

   void update_hamiltonian(void) { 
      SpinCavityTuple::update_hamiltonian();
      SpinCavityTuple::hamiltonian() += lambda0 * kroneckerProduct( S.Sx(0), os.a(0) + os.ap(0) ).eval(); // Sx (a + a^+) 
      //      SpinCavityTuple::hamiltonian() += sigma1 * kroneckerProduct( S.Sz2(0), os.a(1) + os.ap(1) ).eval(); // Sz^2 (b + b^+)
      SpinCavityTuple::hamiltonian() += sigma1 * kroneckerProduct( S.Sx2(0) - S.Sy2(0), os.a(1) + os.ap(1) ).eval(); // Sz^2 (b + b^+)
      //      SpinCavityTuple::hamiltonian() += sigma1 * kroneckerProduct( S.Sy(0), os.a(1) + os.ap(1) ).eval(); // Sz^2 (b + b^+)            
   }
};

double ipr(const VectorXd &vec) {
   double s = 0;
   for (int m = 0;m < vec.size(); m++) { 
      s += pow( vec(m), 4.0 );
   }
   return s;
}

int main(int argc, char **argv)
{
  //   Spin1w2Cavities s1ab(20, 50);
   Spin1w2Cavities s1ab(10, 30);  
   double Ds = 1;  // D parameter for spin 
   double Es = 0.1; // E parameter for spin 
   s1ab.S.S(0).D = Ds; 
   s1ab.S.S(0).E = Es;  
   s1ab.S.S(0).B << 0, 0, 0; // Bx, By, Bz for spin 

   std::cout << "# D " << s1ab.S.S(0).D << std::endl;
   
   // in zero magnetic field - spin eigenvalues are split from the ground state by D +/- E 
   s1ab.os.cavity(0).omega_c = Ds+Es;    // frequency of first cavity 
   s1ab.os.cavity(1).omega_c = 2. * Es;  // frequency of second cavity 

   for (double lambda = 1e-3; lambda <= s1ab.os.cavity(1).omega_c; lambda += 1e-3) { 

   s1ab.lambda0 = lambda;    // coupling between spin and first cavity Sx.(a + a^+) type 
   s1ab.sigma1 =  lambda;     // coupling between spin and second cavity Sz^2.(a + a^+) type 

   s1ab.update_hamiltonian();


   //   MatrixXcd H = s1ab.hamiltonian(); // converts sparse matrix hamiltonian into dense matrix for dense diagonalization 
   //   SelfAdjointEigenSolver<MatrixXcd> eigensolver(H);
   MatrixXd H = s1ab.hamiltonian().real();
   SelfAdjointEigenSolver<MatrixXd> eigensolver(H);   
   VectorXd eval;
   if (eigensolver.info() != Success) abort();
   eval = eigensolver.eigenvalues();
   MatrixXd evec;
   evec = eigensolver.eigenvectors();
   std::cout << "# evec trace " << evec.trace() << std::endl;
   for (int n = 0; n < s1ab.size()-1; n++) {
     std::cout << lambda << "  " << n << "   " << eval(n) << "    " << eval(n+1)-eval(n) << "   " << ipr(evec.col(n)) << std::endl;
   }
   std::cout << std::endl;
   }
}
