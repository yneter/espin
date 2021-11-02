#ifndef ESPINTENSOR
#define ESPINTENSOR

#include <iostream>
#include <string>
#include <initializer_list>
#include <vector>
#include <algorithm>
#include <math.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
using namespace Eigen;

class eSpinSparseTensor {
   std::vector<int> dim_list;
   std::vector<int> dim_prod;
  
   int enc_run(int dimindex, int t) 
   {
      if (dimindex == 0 && t >= 0 && t < dim_list[0]) {
	return t; 
      } else {
	std::cerr << "# indexing error in enc_run " << dimindex << std::endl;	
	return 0;
      }
   }

   template<typename... Args> int enc_run(int dimindex, int t, Args... args) 
   {      
      if (t >= 0 && t < dim_list[dimindex]) { 
	 return t * dim_prod[dimindex-1]  + enc_run(dimindex-1, args...);
      } else {
	 std::cerr << "# indexing error in enc_run dimindex" << dimindex << " n=" << t << " but " << "dim="<< dim_list[dimindex] << std::endl;
	 return 0;
      }
   }


   void setvec_run(std::vector<double> &v, int dimindex, double x) 
   {
     //      v[v.size()-1] = x;
      v[0] = x;      
   }

   template<typename... Args> void setvec_run(std::vector<double> &v, int dimindex, double x, Args... args) 
   {      
     //      v[v.size()-1-dimindex] = x;
      v[dimindex] = x;
      setvec_run(v, dimindex-1, args...);
   }


    double x_from_index(int n, int i) {
       return xmin[n] + (xmax[n] - xmin[n]) * (double) i / (double) ( dim_list[n] - 1 );
    }

    int index_from_x_nocut(int n, double x) {
       return ceil( (x - xmin[n])/(xmax[n] - xmin[n]) * (double)( dim_list[n] - 1 ) );
    }

    int index_from_x(int n, double x) {
      if (x >= xmax[n]) { return dim_list[n] - 1; }
      else if (x <= xmin[n]) { return 0; }
      return index_from_x_nocut(n, x);
    }


public :
   std::vector<double> xmax;
   std::vector<double> xmin;

   template<typename... Args> int enc(Args... args) 
   {      
      return enc_run(dim_prod.size() - 1, args...);
   }

   int enc(std::vector<int> vi, bool reverse = true)
   {
      int N = dim_prod.size();
      int s;
      if (reverse) { 
	 s = vi[N-1];
	 for (int n = 0; n < vi.size()-1; n++) {
	    s += vi[n] * dim_prod[N-2-n];
	 }
      } else {
	 s = vi[0];
	 for (int n = 1; n < vi.size(); n++) {
	    s += vi[n] * dim_prod[n-1];
	 }
      }
      return s;
   }
  
   template<typename... Args> void setmax(Args... args) 
   {      
      xmax.resize( dim_list.size() );
      setvec_run(xmax, dim_prod.size() - 1, args...);
   }

   template<typename... Args> void setmin(Args... args) 
   {      
      xmin.resize( dim_list.size() );
      setvec_run(xmin, dim_prod.size() - 1, args...);
   }

   template<typename... Args> int encx(Args... args) 
   {      
      std::vector<double> v( dim_list.size() );
      std::vector<int> out( dim_list.size() );


      setvec_run(v, v.size() - 1, args...);
      for (int n = 0; n < v.size(); n++) {
	 out[n] = index_from_x(n, v[n]);
      }
      return enc(out, false);
   }

   int encx(std::vector<double> x) { 
      std::vector<int> out( dim_list.size() );
      for (int n = 0; n < out.size(); n++) {
	 out[n] = index_from_x(n, x[x.size()-1-n]);
      }
      return enc(out, false);      
   }
  
   void dim(std::vector<int> d) {
      dim_list = d;
      std::reverse(dim_list.begin(), dim_list.end());
      dim_prod.resize( dim_list.size() );
      dim_prod[0] = dim_list[0];
      for (int n = 1; n < dim_list.size(); n++) { 
	 dim_prod[n] = dim_list[n] * dim_prod[n-1];
      }
   }

   // output of dec is not reversed ...
   std::vector<int> dec(int n) { 
      std::vector<int> vd;
      int N = dim_prod.size();
      vd.resize(N);
      for (int s=0; s < N-1; s++) {
	// t * dim_prod[dimindex-1]  + enc_run(dimindex-1, args...);
	 vd[s] = n / dim_prod[N-2-s];
	 n -= vd[s] * dim_prod[N-2-s];
      }
      vd[N-1] = n;
      return vd;
   }

   std::vector<double> decx(int n) { 
      std::vector<int> vd = dec(n);
      std::vector<double> vx( vd.size() );
      for (int n=0; n<vd.size(); n++) {
	 vx[n] = x_from_index(vd.size()-1-n, vd[n]);
      }
      return vx;
   }
};


VectorXcd espintensor(std::vector< std::reference_wrapper< VectorXcd > > &vlist)  {   
  VectorXcd tmp;
  if (vlist.size() == 1) tmp = vlist[0];
  if (vlist.size() >= 2) {
    tmp = kroneckerProduct( static_cast<VectorXcd &>(vlist[0]), static_cast<VectorXcd &>(vlist[1]) ).eval();
  }
  for (int n = 2; n < vlist.size(); n++) {
    tmp = kroneckerProduct( tmp, static_cast<VectorXcd &>(vlist[n]) ).eval();
  }
  return tmp;
}

VectorXcd espintensor(std::vector< std::reference_wrapper< VectorXcd > > &vlist, VectorXcd &v) { 
   vlist.push_back(v);
   return espintensor(vlist);
}
  
template<typename... Args> VectorXcd espintensor(std::vector< std::reference_wrapper< VectorXcd > > &vlist, VectorXcd &v, Args... args) { 
   vlist.push_back(v);
   return espintensor(vlist, args...);
}
  
template<typename... Args> VectorXcd espintensor(VectorXcd &v, Args... args) { 
   std::vector< std::reference_wrapper< VectorXcd > > vlist;
   vlist.push_back(v);
   return espintensor(vlist, args...);
}


int main_test_espintensor()
{
  eSpinSparseTensor a;
  std::vector<int> adim = { 11, 8, 23 };
  a.dim(adim);
  

  a.setmax(10, 7, 22);
  a.setmin(0, 0, 0);

  
  std::cout << a.encx( 3, 1, 22 ) << "  " <<  a.encx( { 3., 1., 22. } ) << "   " << a.enc(3, 1, 22) << "   " << a.enc( {3, 1, 22} )<< std::endl;

  for (double x : a.decx( a.encx( { 9., 7., 20.9 } ) ) ) {
    std::cout << x << "  ";
  }
  std::cout << std::endl;
  return 0;
} 


int main_test_tensor()
{
   VectorXcd a(2);
   a << 1, 2;
   VectorXcd b(2);
   b << 1, 3;
   std::cout << espintensor(a, b, a) << std::endl;
   return 0;
}
#endif
