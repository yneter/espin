#include "espindense.cc"
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

struct SingleSpin1Dt : public SingleSpin1 {
    typedef SingleSpin1::SpinMatrix SpinMatrix;
    double Bac;
    double omega;
  
    SpinMatrix hamiltonian(double time) { 
       return SingleSpin1::hamiltonian() + Bac * S.Sx() * cos(omega * time);
    } 

    void operator()(const Vector3cd &x , Vector3cd &dxdt, double t);    
};

void SingleSpin1Dt::operator() (const Vector3cd &x , Vector3cd &dxdt, double t) { 
    dxdt = -iii * this->hamiltonian(t) * x;
}


int main(int argc, char **argv)
{
    SingleSpin1Dt s1;
    s1.S.D = 1;
    s1.S.B << 0.0, 0.5, 0;
    s1.update_hamiltonian();
    s1.diag();
    s1.Bac = 0.05;
    s1.omega = s1.eval(1) - s1.eval(0);

    typedef runge_kutta_dopri5<Vector3cd, double, Vector3cd, double,vector_space_algebra> H_stepper;
    H_stepper hstep;
    SingleSpin1Dt::SpinVector psi0 = s1.evec.col(0);

    std::complex< double > sx01 = s1.evec.col(1).adjoint() * s1.S.Sx() * psi0;
    std::cout << "# sx01 " << s1.Bac * abs(sx01) << std::endl;
    
    double DT = 0.1;
    double t = 0;
    while (t <= 1000) { 
      int ndt = 10;
      double dt = DT / (double) ndt;
      integrate_n_steps(hstep, s1, psi0, t, dt, ndt);
      Vector3cd w = s1.evec.adjoint() * psi0;
      std::cout << t << "   " << norm( w(0) ) << "    " << norm( w(1) ) << "    " << norm( w(2) ) << "    " << psi0.norm() << std::endl;
      t += DT;
    }
}
