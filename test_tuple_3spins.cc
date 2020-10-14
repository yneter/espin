#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/Geometry>
#include <vector>
#include <boost/scoped_ptr.hpp>

#include "espin/espintuple.cc"

class tipsADT {
    Vector3d C;
    Vector3d B;
    Vector3d A;

    Vector3d M1X;
    Vector3d M1Y;
    Vector3d M1Z;

    Vector3d M12;
public :
    Vector3d a(void) { return A; } 
    Vector3d b(void) { return B; } 
    Vector3d c(void) { return C; } 


    Vector3d mx(void) { return M1X; } 
    Vector3d my(void) { return M1Y; } 
    Vector3d mz(void) { return M1Z; } 
    Vector3d r12(void) { return M12; }

    tipsADT() { 
       Vector3d C55;
       C55 << -2.567, -3.232, 7.571;
       Vector3d C149;
       C149 << -3.305, -6.421, 23.386;
       Vector3d C431;
       C431 << 5.01, -3.232, 7.571;
       Vector3d C243;
       C243 << -3.817, 4.852, 7.571;
       
       C = C149-C55; 
       A = C431-C55;
       B = C243-C55;

       Vector3d M1C96;
       M1C96 << 7.163, 2.283, 6.631;
       Vector3d M1C143;
       M1C143 << 6.001, 2.612, 9.184;
       Vector3d M1C97;
       M1C97 << 7.247, 3.536, 7.270;
       
       M1X = (M1C143 - M1C96).normalized();
       M1Z = M1X.cross(M1C97-M1C96).normalized();
       M1Y = M1Z.cross(M1X);

       Vector3d M1C49;
       M1C49 << -1.575, 2.612, 9.184;
       Vector3d M2C143;
       M2C143 << 6.001, 2.612, 9.184;
       M12 = M2C143 - M1C49;
    };

    Matrix3d rot (void) const { 
       
       Matrix3d R; 
       R.col(0) = M1X;
       R.col(1) = M1Y;
       R.col(2) = M1Z;
       return R.transpose();
    }

};


#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

struct TripletPairTime : public TripletPair { 
    double rate;
    double dt;
    double J0;
    bool pure_singlet; 

    typedef Eigen::Matrix<complexg, TripletPair::matrix_size, 1> spin_state;
    spin_state ut;


    double Jt(double t) {
       if (J0 > 0) 
	 return J0 - rate * t;
       else 
	 return J0 + rate * t;
    }

    TripletPair::SpinMatrix hamiltonian_at_time(double time) { 
       this->J = Jt(time);
       return update_hamiltonian();
    }

    void operator()(const spin_state &x , spin_state &dxdt , double t )
    {
       dxdt = -iii * hamiltonian_at_time(t) * x;
    }

    void run_to_crossing_jumps(int N) { 
       this->hamiltonian_at_time(0.0);
       this->diag();
       int i_singlet;
       for (int i = 0; i < TripletPair::matrix_size; i++) {
	 if (this->singlet_content(i) > 0.9) { 
	   i_singlet = i;
	    std::cerr << "# singlet content " << this->singlet_content(i) << std::endl;
	 }
       }

       ut = this->evec.col(i_singlet);

       double time = fabs(J0/rate);
       double Dt = time / (double) N;
       double t = 0.0;
       for (int n = 0; n < N; n++)  { 
	  this->hamiltonian_at_time(t + Dt/2.0);
	  this->diag();
	  ut = this->evec.adjoint() * ut;
	  for (int i = 0; i < ut.size(); i++) { 
	     ut(i) *= exp(-iii * this->eval(i) * Dt);
	  }
	  ut = this->evec * ut;
	  t += Dt;
       }       

       this->hamiltonian_at_time(time);
       std::cerr << "# norm ut-1: " <<  ut.norm() - 1.0 << std::endl;
       std::cerr << "# ";
       for (int i = 0; i < ut.size(); i++) std::cerr << norm(ut(i)) << " ";
       std::cerr << std::endl;
    }
  
    void run_to_crossing(void) { 
       this->hamiltonian_at_time(0.0);
       this->diag();
       int i_singlet;
       for (int i = 0; i < TripletPair::matrix_size; i++) {
	 if (this->singlet_content(i) > 0.9) { 
	   i_singlet = i;
	 }
       }

       typedef runge_kutta_dopri5<spin_state, double, spin_state, double,vector_space_algebra> rk5_stepper;
       rk5_stepper rk5;

       ut = this->evec.col(i_singlet);
       double time = fabs(J0/rate);
       double t0 = 0.0;
       
       int N = 50;
       double Dt = time / (double) N;
       for (int n = 0; n < N; n++)  { 
	 integrate_adaptive(rk5, *this, ut,  t0, t0 + Dt, dt);
	 t0 += Dt;
         std:cerr << "# norm ut-1: " <<  ut.norm() - 1.0 << " time " << t0 << std::endl;
       }
       this->hamiltonian_at_time(time);
    }

    void load_singlet(void) { 
       ut = singlet();
    }

    const SpinMatrix singlet_projector(void) { 
       return ut * ut.adjoint();
    }

};


int main(int argc, char **argv)
{
    tipsADT tips;

    Vector3d uz_triplet = tips.mz();
    Vector3d uz;
    uz << 0.0, 0.0, 1.0;
    Vector3d uy = uz.cross(uz_triplet).normalized();
    Vector3d ux_perp = uy.cross(uz_triplet);

    Vector3d r1 = atof(argv[1]) * tips.a() + atof(argv[2]) * tips.b() + atof(argv[3]) * tips.c();
    Vector3d r2 = atof(argv[4]) * tips.a() + atof(argv[5]) * tips.b() + atof(argv[6]) * tips.c();

    typedef SpinTuple< TripletSpin, TripletSpin, TripletSpin > TripletTriad;

    //    TripletPairTime triplet_pair;
    TripletTriad triplet3;
    Vector3d Bvec;
    double D = atof(argv[7]);
    double E = atof(argv[8]);
    for (int s = 0; s < triplet3.nspins(); s++) { 
       triplet3.S(s).D = D;
       triplet3.S(s).E = E;
       triplet3.S(s).rot = tips.rot();
    }
    double GammaD = 51.92;    
    TripletTriad::SpinPositionMatrix rspins;
    rspins.col(0).setZero();
    rspins.col(1) = 0.1 * r1;
    rspins.col(2) = 0.1 * r2;
    triplet3.dipole_dipole_interaction(GammaD, rspins);


    double ax = atof(argv[9]);
    double ay = atof(argv[10]);
    double az = atof(argv[11]);

    TripletTriad::SpinMatrix proj3 = triplet3.Jproj( 3. );
    TripletTriad::SpinMatrix proj1 = triplet3.Jproj( 1. );

    for (double Btot = 0.0; Btot <= 4000.0; Btot += 2.0) { 
      Vector3d uB = (ax * ux_perp +  ay * uy + az * uz_triplet );
      uB.normalize();
      Vector3d B = Btot * uB;
      for (int s = 0; s < triplet3.nspins(); s++) { 
	 triplet3.S(s).B = B;
      }
      triplet3.update_hamiltonian();
      triplet3.diag(); 

      std::vector<int> selected_states;
      for (int i = 0; i < triplet3.size(); i++) { 
	 if ( real( triplet3.matrix_elem(proj1, i, i) ) < 0.1 ) { 
	//	 if ( real( triplet3.matrix_elem(proj3, i, i) ) > 0.5 ) { 
	   selected_states.push_back(i);
	 }

	 for (int i = 0; i < selected_states.size(); i++) { 
	    for (int j = 0; j < i; j++) { 
	       int a = selected_states[i];
	       int b = selected_states[j];
	       double Amp = norm(triplet3.Sx(a,b)) + norm(triplet3.Sy(a,b));	     
	       double epsilon = triplet3.eval(a) - triplet3.eval(b);
	       double dd = 60.0;
	       double D0 = 1500.0;
	       string tag; 
	       if ( fabs( Btot + epsilon - D0 ) < dd || fabs( Btot - epsilon + D0 ) < dd ) { 
		  tag = " SM "; 
	       } else if ( fabs( Btot - epsilon - D0 ) < dd ) { 
		  tag = " SP ";
	       } else { 
		  tag = " S0 "; 
	       }
	       double Ta = real( triplet3.matrix_elem(proj1, a, a) );
	       double Tb = real( triplet3.matrix_elem(proj1, b, b) );
	       double Qa = real( triplet3.matrix_elem(proj3, a, a) );
	       double Qb = real( triplet3.matrix_elem(proj3, b, b) );
	       std::cout << "   " << epsilon << "   " << 0.357 * Btot <<  "   " << Ta << "   " <<  Tb << "   " << Qa << "    " << Qb << "    " << Amp << tag << std::endl; 
	    }
	 }
	 /**
	 for (int i = 0; i < triplet3.size(); i++) { 
	    for (int j = 0; j < triplet3.size(); j++) { 
	       std::cout << "   " << triplet3.eval(i) - triplet3.eval(j) << "   " << 0.357 * Btot << "    " << triplet3.eval(i) << " Z " << std::endl;
	    }
	 }
	 **/
       }
    }
}


 

