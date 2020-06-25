#include "espindense.cc"



int main_check_spin()
{
  Spin1 S1;
  std::cout << (S1.sx * S1.sy - S1.sy * S1.sx + iii * S1.sz).norm() << std::endl;
  Spin1p5 S1p5;
  std::cout << (S1p5.sx * S1p5.sy - S1p5.sy * S1p5.sx + iii * S1p5.sz).norm() << std::endl;

  Spin2 S2;
  std::cout << (S2.sx * S2.sy - S2.sy * S2.sx + iii * S2.sz).norm() << std::endl;
  Spin2p5 S2p5;
  std::cout << (S2p5.sx * S2p5.sy - S2p5.sy * S2p5.sx + iii * S2p5.sz).norm() << std::endl;

  Spin3 S3;
  std::cout << (S3.sx * S3.sy - S3.sy * S3.sx + iii * S3.sz).norm() << std::endl;
  Spin3p5 S3p5;
  std::cout << (S3p5.sx * S3p5.sy - S3p5.sy * S3p5.sx + iii * S3p5.sz).norm() << std::endl;

  Spin4 S4;
  std::cout << (S4.sx * S4.sy - S4.sy * S4.sx + iii * S4.sz).norm() << std::endl;
}



int main() {
    return main_check_spin();
}

// #endif 


