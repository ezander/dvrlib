
#include <iostream>

#include <string.h>
#include <vector>
#include <stdio.h>

#include "gsl_wrapper.h"
#include "recon.h"
#include "utils.h"
#include "recon_system.h"
#include "vdi2048.h"


void example_VDI2048_1(){
  dvrlib::recon_system system;

  system.add_var("m_FDKeI",  46.241, 2.500);
  system.add_var("m_FDKeII", 45.668, 2.500);
  system.add_var("m_SpI",    44.575, 0.535);
  system.add_var("m_SpII",   44.319, 0.532);
  system.add_var("m_V",       0.525, 0.105);
  system.add_var("m_HK",     69.978, 0.854);
  system.add_var("m_A7",     10.364, 0.168);
  system.add_var("m_A6",      3.744, 0.058);
  system.add_var("m_A5",      4.391, 0.058);
  system.add_var("m_HDNK",   18.498, 0.205);
  system.add_var("m_D",       2.092, 0.272);

  system.add_covariance_coeff("m_FDKeI", "m_FDKeII", 0.2);
  system.add_covariance_coeff("m_SpI", "m_SpII", 0.4);

  // get and print the covariance matrix
  dvrlib::matrix S_x = system.get_covariance_matrix();
  PRINT_TITLE("Matrix S_X (compare A9, values from A6 and A8)");
  PRINT(S_x);

  // define the constraints
  //   m_FD1 = m_FDKeI + m_FDKeII - 0.2 * m_V
  //   m_FD2 = m_SpI + m_SpII - 0.6 * m_V
  //   m_FD3 = m_HK + m_A7 + m_A6 + m_A5 + 0.4 * m_V
  //   m_HDAnz = m_A7 + m_A6 + m_A5
  // plus
  //   M_FD1 - M_FD2 = 0
  //   M_FD2 - M_FD3 = 0
  //   M_HDAnz - M_HDNK = 0
  double Fc[][11]={
    {1, 1, -1, -1, 0.4,  0,  0,  0,  0,  0, 0},
    {0, 0,  1,  1,  -1, -1, -1, -1, -1,  0, 0},
    {0, 0,  0,  0,   0,  0,  1,  1,  1, -1, 0},
  };
  dvrlib::matrix F(3, 11, Fc);

  dvrlib::vector x=system.get_values();
  dvrlib::vector v(x.size());

  dvrlib::lin_recon(F*x, S_x, F, v);

  PRINT_TITLE("Vector v (compare A20)");
  PRINT(v);

  // test that the constraints are fulfilled
  PRINT_TITLE("Constraint fulfillment F*(v+x) (should be close to zero)");
  PRINT(F*(x+v));

  // compute and print covariance S_v
  dvrlib::matrix S_v(S_x.size1(), S_x.size2());
  dvrlib::lin_cov_update(S_x, F, S_v);

  PRINT_TITLE("Matrix S_v (compare A21)");
  PRINT(S_v);

  // compute and print covariance update S_X_new
  dvrlib::matrix S_xnew(S_x.size1(), S_x.size2());
  S_xnew = S_x - S_v;

  PRINT_TITLE("Matrix S_xnew (compare A33)");
  PRINT(S_xnew);


}

void example_VDI2048_2(){
    dvrlib::recon_system system;
  
  //system.add_var("m_FDKeI",  46.241, 0.800);
  //system.add_var("m_FDKeII", 45.668, 0.790);
  system.add_var("m_FDKeI",  46.241, 2.500);
  system.add_var("m_FDKeII", 45.668, 2.500);
  system.add_var("m_SpI",    44.575, 0.535);
  system.add_var("m_SpII",   44.319, 0.532);
  system.add_var("m_V",       0.525, 0.105);
  system.add_var("m_HK",     69.978, 0.854);
  system.add_var("m_A7",     10.364, 0.168);
  system.add_var("m_A6",      3.744, 0.058);
  system.add_var("m_A5",      4.391, 0.058);
  system.add_var("m_HDNK",   18.498, 0.205);
  system.add_var("m_D",       2.092, 0.272);

  system.add_var("m_FD1",   0, -1);
  system.add_var("m_FD2",   0, -1);
  system.add_var("m_FD3",   0, -1);
  system.add_var("m_HDAnz", 0, -1);

  system.add_covariance_coeff("m_FDKeI", "m_FDKeII", 0.2);
  system.add_covariance_coeff("m_SpI", "m_SpII", 0.4);


  // get and print the covariance matrix
  dvrlib::matrix S_x = system.get_covariance_matrix();
  PRINT_TITLE("Matrix S_X (compare A9, values from A6 and A8)");
  PRINT(S_x);

  // define the constraints
  //   m_FD1 = m_FDKeI + m_FDKeII - 0.2 * m_V
  //   m_FD2 = m_SpI + m_SpII - 0.6 * m_V
  //   m_FD3 = m_HK + m_A7 + m_A6 + m_A5 + 0.4 * m_V
  //   m_HDAnz = m_A7 + m_A6 + m_A5
  // plus
  //   M_FD1 - M_FD2 = 0
  //   M_FD2 - M_FD3 = 0
  //   M_HDAnz - M_HDNK = 0
  double Fc[][15]={
    {1, 1, 0, 0, -0.2, 0, 0, 0, 0,  0, 0, -1,  0,  0,  0},
    {0, 0, 1, 1, -0.6, 0, 0, 0, 0,  0, 0,  0, -1,  0,  0},
    {0, 0, 0, 0,  0.4, 1, 1, 1, 1,  0, 0,  0,  0, -1,  0},
    {0, 0, 0, 0,    0, 0, 1, 1, 1,  0, 0,  0,  0,  0, -1},
    {0, 0, 0, 0,    0, 0, 0, 0, 0,  0, 0,  1, -1,  0,  0},
    {0, 0, 0, 0,    0, 0, 0, 0, 0,  0, 0,  0,  1, -1,  0},
    {0, 0, 0, 0,    0, 0, 0, 0, 0, -1, 0,  0,  0,  0,  1},
  };
  dvrlib::matrix F(7, 15, Fc);
  PRINT_TITLE("Constraints (compare A1-A4 and A12)");
  system.print_constraints(F);

  PRINT_TITLE("Constraints as Matrix F");
  PRINT(F);
  
  // compute and print reconciliation vector v
  dvrlib::vector x=system.get_values();
  dvrlib::vector v(x.size());

  dvrlib::lin_recon(F*x, S_x, F, v);

  PRINT_TITLE("Vector v (compare A20)");
  PRINT(v);

  // test that the constraints are fulfilled
  PRINT_TITLE("Constraint fulfillment F*(v+x) (should be close to zero)");
  PRINT(F*(x+v));

  // compute and print covariance S_v
  dvrlib::matrix S_v(S_x.size1(), S_x.size2());
  dvrlib::lin_cov_update(S_x, F, S_v);

  PRINT_TITLE("Matrix S_v (compare A21)");
  PRINT(S_v);
  
}


void example_allgemein(
	const char* name[],
	const double wert[],
	const double konfidenz[],
	const char* cov1[],
	const char* cov2[],
	const double wert_cov[],
	gsl_matrix *zzz
	                   )
    {
    dvrlib::recon_system system;

  std::cout <<  "example_allgemein " << std::endl;
  std::cout <<  name[0] << std::endl;

/*
 for (int i = 0; i< name.size(); i++)
      {
	  system.add_var(name,  wert, konfidenz);
	  //system.add_var("m_FDKeI",  46.241, 0.800);
      }
*/
 /*
  for (int i = 0; i< cov1.size(); i++)
      {
	  system.add_covariance_coeff(cov1,  cov2, wert_cov);
	  //system.add_covariance_coeff("m_FDKeI", "m_FDKeII", 0.2);
      }

  // get and print the covariance matrix
  matrix S_x = system.get_covariance_matrix();
  PRINT_TITLE("Matrix S_X (compare A9, values from A6 and A8)");
  PRINT(S_x);*/

  // define the constraints
  //   m_FD1 = m_FDKeI + m_FDKeII - 0.2 * m_V
  //   m_FD2 = m_SpI + m_SpII - 0.6 * m_V
  //   m_FD3 = m_HK + m_A7 + m_A6 + m_A5 + 0.4 * m_V
  //   m_HDAnz = m_A7 + m_A6 + m_A5
  // plus
  //   M_FD1 - M_FD2 = 0
  //   M_FD2 - M_FD3 = 0
  //   M_HDAnz - M_HDNK = 0
  /*double Fc[][15]={
    {1, 1, 0, 0, -0.2, 0, 0, 0, 0,  0, 0, -1,  0,  0,  0},
    {0, 0, 1, 1, -0.6, 0, 0, 0, 0,  0, 0,  0, -1,  0,  0},
    {0, 0, 0, 0,  0.4, 1, 1, 1, 1,  0, 0,  0,  0, -1,  0},
    {0, 0, 0, 0,    0, 0, 1, 1, 1,  0, 0,  0,  0,  0, -1},
    {0, 0, 0, 0,    0, 0, 0, 0, 0,  0, 0,  1, -1,  0,  0},
    {0, 0, 0, 0,    0, 0, 0, 0, 0,  0, 0,  0,  1, -1,  0},
    {0, 0, 0, 0,    0, 0, 0, 0, 0, -1, 0,  0,  0,  0,  1},
  };*/
  /*matrix F(zeilen_F_Matrix, spalten_F_Matrix, Fc);

  PRINT_TITLE("Constraints (compare A1-A4 and A12)");
  system.print_constraints(F);

  PRINT_TITLE("Constraints as Matrix F");
  PRINT(F);

  // compute and print reconciliation vector v
  vector x=system.get_values();
  vector v(x.size());

  lin_recon(F*x, S_x, F, v);

  PRINT_TITLE("Vector v (compare A20)");
  PRINT(v);

  // test that the constraints are fulfilled
  PRINT_TITLE("Constraint fulfillment F*(v+x) (should be close to zero)");
  PRINT(F*(x+v));

  // compute and print covariance S_v
  matrix S_v(S_x.size1(), S_x.size2());
  lin_cov_update(S_x, F, S_v);

  PRINT_TITLE("Matrix S_v (compare A21)");
  PRINT(S_v);
  */

}


void example_VDI2048() {
    //example_VDI2048_2();
    example_VDI2048_1();

    /*
    gsl_matrix *zzz;
    zzz = gsl_matrix_calloc(7,15);
    gsl_matrix_set (zzz, 0,0, 1.0 );
    gsl_matrix_set (zzz, 0,1, 1.0 );
    gsl_matrix_set (zzz, 0,4, -0.2 );
    gsl_matrix_set (zzz, 0,11, -1.0 );

    gsl_matrix_set (zzz, 1,2, 1.0 );
    gsl_matrix_set (zzz, 1,3, 1.0 );
    gsl_matrix_set (zzz, 1,4, -0.6 );
    gsl_matrix_set (zzz, 1,12, -1.0 );

    gsl_matrix_set (zzz, 2,4, 0.4 );
    gsl_matrix_set (zzz, 2,5, 1.0 );
    gsl_matrix_set (zzz, 2,6, 1.0 );
    gsl_matrix_set (zzz, 2,7, 1.0 );
    gsl_matrix_set (zzz, 2,8, 1.0 );
    gsl_matrix_set (zzz, 2,13, -1.0 );

    gsl_matrix_set (zzz, 3,6, 1.0 );
    gsl_matrix_set (zzz, 3,7, 1.0 );
    gsl_matrix_set (zzz, 3,8, 1.0 );
    gsl_matrix_set (zzz, 3,14, -1.0 );

    gsl_matrix_set (zzz, 4,12, 1.0 );
    gsl_matrix_set (zzz, 4,13, -1.0 );

    gsl_matrix_set (zzz, 5,13, 1.0 );
    gsl_matrix_set (zzz, 5,14, -1.0 );

    gsl_matrix_set (zzz, 6,9, -1.0 );
    gsl_matrix_set (zzz, 6,14, 1.0 );

    const char* name[] = {"m_FDKeI", "m_FDKeII",
	    "m_SpI", "m_SpII", "m_V", "m_HK",
	    "m_A7", "m_A6", "m_A5", "m_HDNK",
	    "m_D"};

    const double wert[] = 	{46.241,45.668,
	    44.575, 44.319, 0.525, 69.978,
	    10.364, 3.744, 4.391, 18.498,
	    2.092};

    const double konfidenz[] = {2.500, 2.500,
	    0.535, 0.532, 0.105, 0.854,
	    0.168, 0.058, 0.058 , 0.205,
	    0.272  };

    const char* cov1[] = {"m_FDKeI", "m_SpI"};
    const char* cov2[] = {"m_FDKeII", "m_SpII"};
    const double wert_cov[] = {0.2 , 0.4};

    example_allgemein(
	    name,
	    wert,
	    konfidenz,
	    cov1 ,
	    cov2 ,
	    wert_cov,
	    zzz);
    */
}  


void schnittstelle::example_VDI2048(){
    //example_VDI2048_2();
    example_VDI2048_1();
}



