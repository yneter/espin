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
      SpinCavityTuple::hamiltonian() += sigma1 * kroneckerProduct( S.Sz2(0), os.a(1) + os.ap(1) ).eval(); // Sz^2 (b + b^+)
   }
};



int main(int argc, char **argv)
{
   Spin1w2Cavities s1ab(3, 3);
   double Ds = 1;  // D parameter for spin 
   double Es = 0.1; // E parameter for spin 
   s1ab.S.S(0).D = Ds; 
   s1ab.S.S(0).E = Es;  
   s1ab.S.S(0).B << 0, 0, 0; // Bx, By, Bz for spin 

   // in zero magnetic field - spin eigenvalues are split from the ground state by D +/- E 
   s1ab.os.cavity(0).omega_c = Ds+Es;    // frequency of first cavity 
   s1ab.os.cavity(1).omega_c = 2. * Es;  // frequency of second cavity 
   s1ab.lambda0 = 0.03;    // coupling between spin and first cavity Sx.(a + a^+) type 
   s1ab.sigma1 = 0.03;     // coupling between spin and second cavity Sz^2.(a + a^+) type 

   s1ab.update_hamiltonian();

   MatrixXcd H = s1ab.hamiltonian(); // converts sparse matrix hamiltonian into dense matrix for dense diagonalization 
   VectorXd eval;
   SelfAdjointEigenSolver<MatrixXcd> eigensolver(H);
   if (eigensolver.info() != Success) abort();
   eval = eigensolver.eigenvalues();
   for (int n = 0; n < s1ab.size(); n++) {
     std::cout << n << "   " << eval(n) << std::endl;
   }
}
