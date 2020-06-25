#include <iostream>
#include <vector>
#include "espindense.cc"

using namespace Eigen;
using namespace std;

int main() 
{
   SpinPair<Spin3, Spin2p5> sp;
   double s = 0.0;
   for (double D = 0.3; D <= 1.5; D += 1e-5) { 
     sp.S1.D = D;
     sp.S2.D = 1.2;
     sp.S1.B << 0.2, 0.4, 0.9;
     sp.S2.B << 0.2, -0.9, 0.9;
     sp.J = 0.5;
     sp.update_hamiltonian();
     sp.diag();
     s += sp.eval(0);
   }
   std::cerr << s << std::endl;
}

