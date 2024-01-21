

#include <iostream>
#include "MD_GB.hpp"

double density = 0.28;
double temperature = 1;

int main(int argc, const char * argv[]) {

  int is_init_file = 0;
  std::cout << "Statistical mechanics is only for cool people!\n";

  GB* a = new GB(density, temperature);
  a->outputit();
  //  a->loaddata();
  //  exit(0);  

  if (is_init_file==0) {
    a->initialisation();
    a->run(0);
  }else a->run(1);
  a->savedata();

 
  //  a->g_corr();

  return 0;
}
