#include <gsl/gsl_errno.h>
#include "gsl_wrapper.h"
#include "gsl_wrapper_tests.h"
#include "recon_tests.h"
#include "vdi2048.h"

#include <iostream>


int main (void)
{
  dvrlib::gsl_enable_exceptions();

  try{
    //gsl_wrapper_test_suite();
    //recon_test_suite();
	  dvrlib::example_VDI2048();
  }
  catch(const dvrlib::gsl_exception& e){
    std::cout << "Caught GSL exception" << std::endl;
    std::cout << "  reason: " << e.reason << std::endl;
    std::cout << "  file:   " << e.file << std::endl;
    std::cout << "  line:   " << e.line << std::endl;
    std::cout << "  gsl_errno: " << e.gsl_errno << std::endl;
  }
    
  
  return 0;
}
