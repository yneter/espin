#include "espindense.cc"



int main_check_spin()
{
  Spin0p5 S0p5;
  std::cout << (S0p5.sx * S0p5.sy - S0p5.sy * S0p5.sx - iii * S0p5.sz).norm() << std::endl;
  Spin1 S1;
  std::cout << (S1.sx * S1.sy - S1.sy * S1.sx - iii * S1.sz).norm() << std::endl;
  Spin1p5 S1p5;
  std::cout << (S1p5.sx * S1p5.sy - S1p5.sy * S1p5.sx - iii * S1p5.sz).norm() << std::endl;

  Spin2 S2;
  std::cout << (S2.sx * S2.sy - S2.sy * S2.sx - iii * S2.sz).norm() << std::endl;
  Spin2p5 S2p5;
  std::cout << (S2p5.sx * S2p5.sy - S2p5.sy * S2p5.sx - iii * S2p5.sz).norm() << std::endl;

  Spin3 S3;
  std::cout << (S3.sx * S3.sy - S3.sy * S3.sx - iii * S3.sz).norm() << std::endl;
  Spin3p5 S3p5;
  std::cout << (S3p5.sx * S3p5.sy - S3p5.sy * S3p5.sx - iii * S3p5.sz).norm() << std::endl;
  Spin4 S4;
  std::cout << (S4.sx * S4.sy - S4.sy * S4.sx - iii * S4.sz).norm() << std::endl;
  Spin4p5 S4p5;
  std::cout << (S4p5.sx * S4p5.sy - S4p5.sy * S4p5.sx - iii * S4p5.sz).norm() << std::endl;
  Spin5 S5;
  std::cout << (S5.sx * S5.sy - S5.sy * S5.sx - iii * S5.sz).norm() << std::endl;
  return 0;
}



int main() {
    return main_check_spin();
}

// #endif 


