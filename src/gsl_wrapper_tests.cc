#include <cassert>
#include <cmath>
#include <iostream>

#include "gsl_wrapper.h"

namespace dvrlib{

void test_vector() {
    vector v(3, 4.0);
    assert(v.get(0)==4.0);
    assert(v.get(1)==4.0);
    assert(v.get(2)==4.0);

    for (int i = 0; i < 3; i++) {
	    v.set(i, 1.23 + i);
    }
    assert(v.get(0)==1.23);
    assert(v.get(1)==2.23);
    assert(v.get(2)==3.23);

    vector w(3, 5.0);
    assert(w.get(0)==5);

    w+=v;
    assert(w.get(0)==6.23);
    assert(w.get(2)==8.23);
    assert(v.get(0)==1.23);

    v = w * 3;
    assert(v.get(0)==18.69);
    assert(w.get(0)==6.23);

    double src[] = {3, 4, 10, 7};
    vector x(4, src);
    assert(x.get(0)==3);
    assert(x.get(3)==7);

    x-=x;
    assert(x.get(0)==0);
    assert(x.get(1)==0);

    vector z(4, src);
    assert(z.get(0)==3);
    assert(z.get(3)==7);

    x = -z;
    assert(x.get(0)==-3);
    assert(x.get(3)==-7);

    assert(x[0] ==-3);
    assert(x[3] ==-7);

}


void test_vector_view() {
    vector v(10, 4.0);
    vector_view vv = v.subvector(3, 2);
    vv.set(0,2);
    vv.set(1,5);
    assert(v.get(0)==4);
    assert(v.get(3)==2);
    assert(v.get(4)==5);
    vector x(vv);
    x.set(0, 7);
    assert(x.get(0)==7);
    assert(x.get(1)==5);
    assert(vv.get(0)==2);
    assert(v.get(3)==2);
    x.set(1,8);
    vv=x;
    assert(v.get(3)==7);
    assert(v.get(4)==8);
    x=vv;
    vv.set(0,-10);
    assert(x.get(0)==7);
    assert(vv.get(0)==-10);

    vector v3(3);
    v3.set(0,3);
    v3.set(1,-4);
    v3.set(2,12);
    assert(v3.norm1()==19);
    assert(v3.norm2()==13);
}


void test_matrix() {
    matrix m(5, 6, true);
    assert(m.get(0,0)==1);
    assert(m.get(0,1)==0);
    assert(m.get(4,4)==1);
    assert(m.get(4,5)==0);

    m.set(3,4,7);
    m.set(1,3,9);
    matrix m2 = m.transpose();
    assert(m.get(3,4)==7);
    assert(m.get(1,3)==9);
    assert(m2.get(4,3)==7);
    assert(m2.get(3,1)==9);

    matrix A(2,2);
    A.set(0,0,4);
    A.set(0,1,6);
    A.set(1,0,10);
    A.set(1,1,3);
    vector b(2);
    b.set(0,38);
    b.set(1,59);
    vector c=A.linsolve(b);
    assert(c.get(0)==5);
    assert(c.get(1)==3);
    assert(b.get(1)==59);
    assert(A.get(0,0)==4);
    assert(A.get(1,1)==3);

    matrix B(A.inverse());
    double det=-48;
    assert(fabs(B.get(0,0)- +3.0/det)<1e-9);
    assert(fabs(B.get(0,1)- -6.0/det)<1e-9);
    assert(fabs(B.get(1,0)- -10.0/det)<1e-9);
    assert(fabs(B.get(1,1)- +4.0/det)<1e-9);


    matrix A2(2,3);
    A2.set(0,0,4);
    A2.set(0,1,6);
    A2.set(0,2,7);
    A2.set(1,0,10);
    A2.set(1,1,3);
    A2.set(1,2,5);
    vector x(3);
    x.set(0, 1);
    x.set(1, 2);
    x.set(2, 3);
    vector y=A2*x;
    assert(y.get(0)==37);
    assert(y.get(1)==31);

    double src[][3] = {{3, 4, 10}, {7, 5, 11}};
    matrix A3(2, 3, src);
    assert(A3.get(0,0)==3);
    assert(A3.get(1,0)==7);
    assert(A3.get(0,1)==4);
    assert(A3.get(1,2)==11);

    matrix A4(2, 3, src);
    A4 = A3 + A3;
    assert(A4.get(0,0)==6);
    assert(A4.get(1,0)==14);
    assert(A4.get(0,1)==8);
    assert(A4.get(1,2)==22);

    A4 += A3;
    assert(A4.get(0,0)==9);
    assert(A4.get(1,0)==21);
    assert(A4.get(0,1)==12);
    assert(A4.get(1,2)==33);

    A4 -= A3;
    assert(A4.get(0,0)==6);
    assert(A4.get(1,0)==14);
    assert(A4.get(0,1)==8);
    assert(A4.get(1,2)==22);

    A4 = A4 - A3;
    assert(A3.get(0,0)==3);
    assert(A3.get(1,0)==7);
    assert(A3.get(0,1)==4);
    assert(A3.get(1,2)==11);

    A4 = -A3;
    assert(A3.get(0,0)==3);
    assert(A3.get(1,0)==7);
    assert(A3.get(0,1)==4);
    assert(A3.get(1,2)==11);

    assert(A4.get(0,0)==-3);
    assert(A4.get(1,0)==-7);
    assert(A4.get(0,1)==-4);
    assert(A4.get(1,2)==-11);

    A4 *= -2.0;
    assert(A4.get(0,0)==6);
    assert(A4.get(1,0)==14);
    assert(A4.get(0,1)==8);
    assert(A4.get(1,2)==22);

    A4 = A4 * -2.0;
    assert(A4.get(0,0)==-12);
    assert(A4.get(1,0)==-28);
    assert(A4.get(0,1)==-16);
    assert(A4.get(1,2)==-44);

    vector d(A4[1]);
    assert(d.get(0)==-28);
    assert(d.get(1)==-20);
    assert(d.get(2)==-44);
}


void test_matrix_view(){
    matrix m(10, 10, true);
    matrix_view mv = m.submatrix(3,3,2,2);
    assert(mv.get(0,0)==1);
    assert(mv.get(0,1)==0);

    mv.set(0,1,2);
    mv.set(1,0,5);
    assert(mv.get(0,1)==2);
    assert(mv.get(1,0)==5);
    assert(m.get(3,4)==2);
    assert(m.get(4,3)==5);

    matrix m3(2, 2, true);
    m3.set(0,0,5);
    m3.set(1,0,6);
    m3.set(0,1,7);
    m3.set(1,1,8);
    m.submatrix(2,4,2,2) = m3;
    assert(m.get(2,4)==5);
    assert(m.get(3,4)==6);
    assert(m.get(2,5)==7);
    assert(m.get(3,5)==8);
}


void gsl_wrapper_test_suite() {
    test_vector();
    test_vector_view();
    test_matrix();
    test_matrix_view();
}

}
