

#include <iostream>
#include "MD_GB_hpp.hpp"


double density = 0.28; //Before it was 0.33, now it is 0.23
double temperature = 1;

int main(int argc, const char * argv[]) {

  int is_init_file = 0;
  std::cout << "Hello, World!\n";

  GB* a = new GB(density, temperature);
  a->outputit();
  //  a->loaddata();
  //  exit(0);  

  if (is_init_file==0) {
    a->initialisation();
    a->run(0);
  }else a->run(1);
  a->savedata();

   // a->g_corr();
   // a->g_orient();

  return 0;
}
