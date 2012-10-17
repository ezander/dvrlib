#include <vector>
#include <string>
#include <cmath>
#include "gsl_wrapper.h"
#include "recon.h"
#include "utils.h"


using std::string;

#ifdef FOO
void example_VDI2048_1(){
  recon_system system;
  
  system.add_var("m_FDKeI",  46.241, 0.800);
  system.add_var("m_FDKeII", 45.668, 0.790);
  system.add_var("m_SpI",    44.575, 0.535);
  system.add_var("m_SpII",   44.319, 0.532);
  system.add_var("m_V",       0.525, 0.105);
  system.add_var("m_HK",     69.978, 0.854);
  system.add_var("m_A7",     10.364, 0.168);
  system.add_var("m_A6",      3.744, 0.058);
  system.add_var("m_A5",      4.391, 0.058);
  system.add_var("m_HDNK",   18.498, 0.205);
  system.add_var("m_D",       2.092, 0.272);

  double rho_FDKeI_FDKeII = 0.2;
  double rho_SpI_SpII = 0.4;

  matrix S_x = system.getCovarianceMatrix();
  system.addCovariance(S_x, system.find_var("m_FDKeI"), 
		       system.find_var("m_FDKeII"), rho_FDKeI_FDKeII);
  system.addCovariance(S_x, system.find_var("m_SpI"), 
  		       system.find_var("m_SpII"), rho_SpI_SpII);

  PRINT_TITLE("Matrix S_X");
  PRINT(S_x);

  //m_FD1 = m_FDKeI + m_FDKeII - 0.2 * m_V
  //m_FD2 = m_SpI + m_SpII - 0.6 * m_V
  //m_FD3 = m_HK + m_A7 + m_A6 + m_A5 + 0.4 * m_V
  //m_HDAnz = m_A7 + m_A6 + m_A5
  //M_FD1 - M_FD2 = 0
  //M_FD2 - M_FD3 = 0
  //M_HDAnz - M_HDNK = 0
  double Fc[][11]={
    {1, 1, -1, -1, 0.4,  0,  0,  0,  0,  0, 0},
    {0, 0,  1,  1,  -1, -1, -1, -1, -1,  0, 0},
    {0, 0,  0,  0,   0,  0,  1,  1,  1, -1, 0},
  };
  matrix F(3, 11, Fc);

  // 
  vector x=system.getValues();

  PRINT_TITLE("Matrix F + (+extra cols)");
  PRINT(F);
  
  // A15
  PRINT_TITLE("Matrix F*S_X (A15)");
  PRINT(F*S_x);

  PRINT_TITLE("Matrix F*X (A17)");
  PRINT(F*x);
  
  PRINT_SEP;

  vector v(x.size());
  lin_recon(F*x, S_x, F, v);

  PRINT_TITLE("Vector v");
  PRINT(v);
  PRINT(F*(x+v));
   
}
#endif
