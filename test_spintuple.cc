#include <iostream>
#include <vector>
#include "espintuple.cc"

using namespace Eigen;
using namespace std;


int main() 
{
   SpinTuple<Spin3, Spin2p5> sp2;
   typename SpinTuple<Spin1, Spin2p5>::SpinExchangeMatrix J;
   double s2 = 0.0;
   for (double D = 0.3; D <= 1.5; D += 1e-5) { 
     sp2.S(0).D = D;
     sp2.S(1).D = 1.2;
     sp2.S(0).B << 0.2, 0.4, 0.9;
     sp2.S(1).B << 0.2, -0.9, 0.9;
     J.setZero();
     J(0,1) = 0.5;
     sp2.exchange_interaction(J);
     sp2.update_hamiltonian();
     //     std::cout << sp2.hamiltonian() << std::endl;
     sp2.diag();
     s2 += sp2.eval(0);
   }
   std::cerr << s2 << std::endl;
}
