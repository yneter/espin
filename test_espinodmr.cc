#include "espinodmr.cc"
#include "espintuple.cc"

int main_hfe() 
{
    HFE hfe_spins;
    hfe_spins.S(0).D = hfe_spins.S(1).D = 0.1;
    hfe_spins.S(0).rot = Rotation::X(0).eval();
    hfe_spins.S(1).rot = Rotation::X(M_PI/2.0).eval();
    hfe_spins.S(0).g3 << 1.0, 1.1, 1.2;
    hfe_spins.S(1).g3 << 1.3, 1.0, 1.0;
    hfe_spins.J = 5.0/3.0;
    hfe_spins.dJ = 1.0/3.0;
    hfe_spins.t = 0.3;
    
    Merrifield<HFE> merrifield(hfe_spins);
    merrifield.gammaS = -0.001;
    merrifield.gamma = 0.003;

    MerrifieldRate<HFE> mr(hfe_spins);
    mr.gammaS = merrifield.gammaS;
    mr.gamma = merrifield.gamma;

    for (double B = 0.0; B < 20; B += 0.003) { 
       hfe_spins.S(0).B  << 0, 0, B;
       hfe_spins.S(1).B  << 0, 0, B;
       hfe_spins.update_hamiltonian();
       hfe_spins.diag(); // needed for PL_from_rate()
       merrifield.find_rho();
       mr.find_rho();

       double PLa = hfe_spins.PLa( merrifield.rho );
       double PLb = hfe_spins.PLb( merrifield.rho );
       double PLar = hfe_spins.PLa( mr.rho );
       double PLbr = hfe_spins.PLb( mr.rho );

       cout << B << "     " 
	    << merrifield.PL() << "    " << mr.PL() << "     " 
	    << PLa << "    " << PLb << "     " 
	    << PLar << "    " << PLbr << "     " 
	    << merrifield.rho_error() << endl;
    }
    return 0;
}



//
// demonstration code for Merrifield class with both Liouville and rate equation versions
//
int main_merrifield()
{
    TripletPair triplet_pair;
    triplet_pair.S1.D = triplet_pair.S2.D = 0.1;
    triplet_pair.S1.E = triplet_pair.S2.E = 0.0;    
    triplet_pair.S1.B  << 0, 0, 1.0;
    triplet_pair.S2.B  << 0, 0, 1.0;
    triplet_pair.S1.g3 << 1.0, 1.0, 1.01; // g factors in molecular frame 
    triplet_pair.S1.g3 << 1.02, 1.0, 1.00; // g factors in molecular frame 
    //    triplet_pair.S1.rot = Rotation::X(0).eval();
    //    triplet_pair.S2.rot = Rotation::X(M_PI/2.0).eval();
    triplet_pair.S1.rot.random();
    triplet_pair.S2.rot.random();

    triplet_pair.r12 = random_unit_vector();
    //    triplet_pair.J = 5.0/3.0;
    triplet_pair.J = 0.0;
    triplet_pair.Jdip = 0.0;
    triplet_pair.update_hamiltonian();

    double t = 3.0*5.0/3.0;
    TripletPair::SpinMatrix tex = TripletPair::SpinMatrix::Zero();
    for (int i = 0; i < 3; i++) { 
       tex(i*3+i, i*3+i) = t;
    }
    for (int i = 0; i < 3; i++) { 
       for (int j = 0; j < i; j++) { 
	  tex(i*3+j, j*3+i) = 0.5*t;
	  tex(j*3+i, i*3+j) = 0.5*t;
       }
    } 
    triplet_pair.diag(); // needed for PL_from_rate()

    Merrifield<TripletPair> merrifield(triplet_pair);
    merrifield.gammaS = 0.001;
    merrifield.gamma = 0.003;


    MerrifieldRate<TripletPair> mr(triplet_pair);
    mr.gammaS = merrifield.gammaS;
    mr.gamma = merrifield.gamma;
    
    for (double B = 0.0; B < 20; B += 0.003) { 
       triplet_pair.S1.B  << 0, 0, B;
       triplet_pair.S2.B  << 0, 0, B;
       triplet_pair.update_hamiltonian();
       triplet_pair.add_matrix(tex);
       triplet_pair.diag(); // needed for MerrifieldRate
       merrifield.find_rho();
       double PL = merrifield.PL();
       cout << B << "     " << PL << "    " << mr.PL() << "     " << merrifield.rho_error() << endl;
    }
    return 0;
}



int main_diag()
{
    TripletPair triplet_pair;
    triplet_pair.J = -1.0;
    triplet_pair.update_hamiltonian();
    triplet_pair.diag(); // needed for PL_from_rate()

    for (int i = 0; i < TripletPair::matrix_size; i++) { 
      cout << i << "   " << triplet_pair.eval[i] << "    " << triplet_pair.singlet_content(i) << "   " << triplet_pair.triplet_content(i) << "   " << triplet_pair.quintet_content(i) << endl;
    }
    cout << "#" << endl;

    TripletPair t2;
    double t = 2.0;
    TripletPair::SpinMatrix tex = TripletPair::SpinMatrix::Zero();
    for (int i = 0; i < 3; i++) { 
       tex(i*3+i, i*3+i) = t;
    }

    double s = 1.0;
    for (int i = 0; i < 3; i++) { 
       for (int j = 0; j < i; j++) { 
	  tex(i*3+j, j*3+i) = s;
	  tex(j*3+i, i*3+j) = s;
       }
    }

    t2.add_matrix(tex);
    t2.diag(); // needed for PL_from_rate()
    for (int i = 0; i < TripletPair::matrix_size; i++) { 
      cout << i << "   " << t2.eval[i] << "    " << t2.singlet_content(i) << "   " << t2.triplet_content(i) << "   " << t2.quintet_content(i) << endl;
    }

    return 1.0;



}


//
// demonstration code for MR/ODMR colormap as function of B and frequency 
//
int main_odmr()
{
    TripletPair triplet_pair;
    Vector3d Bvec;

    triplet_pair.S1.D = triplet_pair.S2.D = 1.0;
    triplet_pair.S1.E = triplet_pair.S2.E = 0.15;
    triplet_pair.J = 0.0;
    triplet_pair.Jdip = 0.03;
    double Bz = 5.0;
    srand(1);

    double quintet_max = 0.0;

    Rotation quintet_t1_rot;
    Rotation quintet_t2_rot;
    Vector3d quintet_rdip;

    int Nsamples = 5000;
    std::vector<double> slist(Nsamples);

    double theta = 1.1;
    for (int count = 0; count < Nsamples; count++) { 
       Rotation triplet1_rot, triplet2_rot;
       triplet_pair.S1.rot = (Rotation::Y(0) * Rotation::X(theta)).eval();
       triplet_pair.S2.rot.random();
       triplet_pair.r12 = random_unit_vector();
       triplet_pair.S1.B << 0, 0, Bz;
       triplet_pair.S2.B << 0, 0, Bz;
       triplet_pair.update_hamiltonian();
       triplet_pair.diag(); 
       double si = 0.0;
       for (int i = 0; i < triplet_pair.matrix_size ; i++) { 
	  double qi = triplet_pair.quintet_content(i);
	  si += pow(qi, 4.0);
       }
       if (si > quintet_max) { 
	  quintet_max = si;
	  quintet_t1_rot = triplet_pair.S1.rot;
	  quintet_t2_rot = triplet_pair.S2.rot;
	  quintet_rdip = triplet_pair.r12;
       }
       slist[count] = si;
    }

    cout << "# quintet_max " << quintet_max << endl;
    cout << "# triplet Euler angles :" << endl;
    cout << "# ";
    for (int i = 0; i < 3; i++) cout << quintet_t1_rot.euler_angles()[i] << "   ";
    cout << endl;
    cout << "# ";
    for (int i = 0; i < 3; i++) cout << quintet_t2_rot.euler_angles()[i] << "   ";
    cout << endl;
    cout << "# Rdip :" << endl;
    cout << "# ";
    for (int i = 0; i < 3; i++) cout << quintet_rdip[i] << "   ";
    cout << endl;

    Bvec << 0, 0, Bz;
    triplet_pair.load_field_basis_Hamiltonian(quintet_t1_rot, quintet_t2_rot, quintet_rdip, Bvec );
    triplet_pair.diag();
    cout << "# qunitet/triplet projections at B = " << Bz << endl;
    for (int i = 0; i < triplet_pair.matrix_size ; i++) 
      cout << "# " << triplet_pair.quintet_content(i) << "    " << triplet_pair.triplet_content(i) << "    " << triplet_pair.singlet_content(i) << "    " << triplet_pair.sz_elem(i) << endl;


    double cos1z = quintet_t1_rot.matrix()(2,2);
    double cos2z = quintet_t2_rot.matrix()(2,2);
    cout << "# angles to field " << endl;
    cout << "# " << theta << "   " << "    " << quintet_max << "   "  << cos1z << "    " << cos2z << "    " << endl;

    ODMR_Signal<TripletPair> odmr_from_triplets(triplet_pair);    

    const int N_averages = 10000;
    double omega_span = 10.0;     
    const int n_omega_samples = 1000;

    double B = 5.0;
    //    for (double B = 0; B < 3.0; B += 0.01) { 

       vector<complexg> chi_B(n_omega_samples, 0.0);
       vector<double> odmr_B(n_omega_samples, 0.0);
    
       for (int sample = 0; sample < N_averages; sample++) {        
	  Rotation rot1;
	  rot1.random();
	  Rotation r1_sample = (rot1.matrix() * quintet_t1_rot.matrix()).eval();
	  Rotation r2_sample = (rot1.matrix() * quintet_t2_rot.matrix()).eval();
	  Vector3d rdip_sample = rot1.matrix() * random_unit_vector(); // rot1.matrix() * quintet_rdip;
	  Bvec << 0, 0, B;
	  triplet_pair.load_field_basis_Hamiltonian(r1_sample, r2_sample, rdip_sample, Bvec);
	  triplet_pair.diag();

	  odmr_from_triplets.update_from_spin_hamiltonian();
	  odmr_from_triplets.load_rho0_thermal(10.0);
	  //	  odmr_from_triplets.load_rho0_from_singlet();
	  odmr_from_triplets.gamma = 3e-3;
	  odmr_from_triplets.gamma_diag = 3e-3;

	  for (int omega_index = 0; omega_index < n_omega_samples; omega_index++) { 
	    double omega = omega_span * (double)omega_index/(double)(n_omega_samples-1);
	    chi_B[omega_index] += odmr_from_triplets.chi1(omega)  / (double) N_averages;
	    odmr_B[omega_index] += odmr_from_triplets.odmr(omega) / (double) N_averages;
	  }
       }


       for (int omega_index = 0; omega_index < n_omega_samples; omega_index++) { 
	 double omega = omega_span * (double)omega_index/(double)(n_omega_samples-1);
	 cout  << B << "    " << omega << "   " << real(chi_B[omega_index]) << "    " << imag(chi_B[omega_index]) << "     " << odmr_B[omega_index] << endl;
       }

       cout << endl;
       return 0;
}

int main() {
    return main_hfe();
}


