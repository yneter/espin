#include "espinsparse.cc"
#include "espincavity.cc"



int main(int argc, char **argv) { 
    double errmax = 1e-6;
    const static int NL = 60;
    std::vector<int> spin_size = { 2 };
    SpinTupleInCavityMaster sic(NL, spin_size);
    std::vector<double> Delta = { 0.03 };
    bool first = true;
    MatrixXcd rho_prev;
    for (int s = 0; s < sic.nspins(); s++) { 
       std::cout << "# Delta" << s << " " << Delta[s] << std::endl;
    }
    sic.gamma = 0.02;
    sic.nav = 0;
    sic.gammaS = 0.0;
    double n = 50.0;
    //    sic.F = sic.gamma * sqrt(n) / 2.0;
    sic.F = 0.05;    
    std::cout << "# F " << sic.F << std::endl;
    std::cout << "# gamma " << sic.gamma << std::endl;
    double domega = 0.005;
    for (double domega = -0.3; domega <= 0.3; domega += 0.00053) {
    //    for (double domega = -0.05; domega <= 0.05; domega += 0.0033) {
       for (int s = 0; s < sic.nspins(); s++) { 
	  sic.Bz[s] = Delta[s] + domega;
	  sic.Omega[s]  = 0.05;
       }
       sic.omegaC = domega;
       sic.gammaSup = 0.0;
       SparseMatrix<complexg> master_matrix;
       master_matrix = sic.update_master_from_liouville();

       /** eigenvalues 
       MatrixXcd master_matrix_dense = master_matrix;
       VectorXcd eval;
       ComplexEigenSolver<MatrixXcd> eigensolver(master_matrix_dense);
       if (eigensolver.info() != Success) abort();
       eval = eigensolver.eigenvalues();
       for (int i = 0; i < eval.size(); i++) {
	  if (real(eval(i)) > -3.0 * sic.gamma) { 
	     std::cout << domega << "  " << i << "   " << real(eval(i)) << "  " << imag(eval(i)) << std::endl;
	  }
       }
       **/

       double err, err_prev;

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
	  err = sic.find_rho_approx_nondiag(3000, errmax);
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
