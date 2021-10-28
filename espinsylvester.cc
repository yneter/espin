#ifndef ESPIN_SYLVESTER
#define ESPIN_SYLVESTER

using namespace Eigen;
using namespace std;

//
// solves A X + X B = C
// 
void sylvester(const MatrixXcd &A, const MatrixXcd &B, const MatrixXcd &C, MatrixXcd &X) { 
    int matrix_size = A.rows();
    ComplexSchur<MatrixXcd> schurAadj(matrix_size);
    schurAadj.compute(A.adjoint());
    ComplexSchur<MatrixXcd> schurB(matrix_size);
    schurB.compute(B );

    std::cerr << "# finished Schur " << std::endl;

    // (U Taa U^+)^+ X + X V Tb V^+ = C 
    // Taa^+ U^+ X V + U^+ X V Tb = U^+ C V 

    MatrixXcd Cuv = schurAadj.matrixU().adjoint() * C * schurB.matrixU();
    // Taa^+ Y + Y Tb = Cuv 

    for (int n = 0; n < matrix_size; n++) { 

       for (int m = n; m < matrix_size; m++) { 
	  // conj(Ta_{k n}) Y_{k m} + Y_{n k} Tb_{k, m}      = C_{n m}
	  complex<double> Snm = 0.0;
	  for (int k = 0; k < n; k++) { 
	     Snm += conj( schurAadj.matrixT()(k, n) ) * X(k, m);
	  }	  
	  for (int k = 0; k < m; k++) { 
	     Snm += X(n, k) * schurB.matrixT()(k, m);
	  }

	  X(n, m) = ( Cuv(n, m) - Snm ) / ( conj( schurAadj.matrixT()(n, n) ) + schurB.matrixT()(m, m) );

	  if (m != n) { 
	     complex<double> Smn = 0.0;
	     for (int k = 0; k < m; k++) { 
	        Smn += conj( schurAadj.matrixT()(k, m) ) * X(k, n);
	     }
	     for (int k = 0; k < n; k++) { 
	        Smn += X(m, k) * schurB.matrixT()(k, n);
	     }
	     X(m, n) = ( Cuv(m, n) - Smn ) / ( conj( schurAadj.matrixT()(m, m) ) + schurB.matrixT()(n, n) );
	  }
       }
    }

    X = schurAadj.matrixU() * X * schurB.matrixU().adjoint();

    //    std::cout << "# error " << ( A * X + X * B - C ).norm() << std::endl;

}

//
// special case B = adjoint(A)
// 

void sylvester_adjoint(const MatrixXcd &A, const ComplexSchur<MatrixXcd> &schurAadj, const MatrixXcd &C, MatrixXcd &X) { 
    int matrix_size = A.rows();

    MatrixXcd Cuv = schurAadj.matrixU().adjoint() * C * schurAadj.matrixU();

    for (int n = 0; n < matrix_size; n++) { 

       for (int m = n; m < matrix_size; m++) { 
	  complex<double> Snm = 0.0;
	  for (int k = 0; k < n; k++) { 
	     Snm += conj( schurAadj.matrixT()(k, n) ) * X(k, m);
	  }	  
	  for (int k = 0; k < m; k++) { 
	     Snm += X(n, k) * schurAadj.matrixT()(k, m);
	  }

	  X(n, m) = ( Cuv(n, m) - Snm ) / ( conj( schurAadj.matrixT()(n, n) ) + schurAadj.matrixT()(m, m) );

	  if (m != n) { 
	     complex<double> Smn = 0.0;
	     for (int k = 0; k < m; k++) { 
	        Smn += conj( schurAadj.matrixT()(k, m) ) * X(k, n);
	     }
	     for (int k = 0; k < n; k++) { 
	        Smn += X(m, k) * schurAadj.matrixT()(k, n);
	     }
	     X(m, n) = ( Cuv(m, n) - Smn ) / ( conj( schurAadj.matrixT()(m, m) ) + schurAadj.matrixT()(n, n) );
	  }
       }
    }

    X = schurAadj.matrixU() * X * schurAadj.matrixU().adjoint();

    //    std::cout << "# error " << ( A * X + X * A.adjoint() - C ).norm() << std::endl;
}

#endif

/***
main()
{
  int N = 10;
  MatrixXcd A = MatrixXcd::Random(N, N);
  MatrixXcd B = MatrixXcd::Random(N, N);
  MatrixXcd C = MatrixXcd::Random(N, N);
  MatrixXcd X = MatrixXcd::Random(N, N);
  ComplexSchur<MatrixXcd> schurAadj(N);
  schurAadj.compute(A.adjoint());

  sylvester(A,B,C,X);

  sylvester_adjoint(A,schurAadj,C,X);
}
***/
