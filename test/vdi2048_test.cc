#include <catch2/catch_all.hpp>
#include "gsl_wrapper.h"
#include "recon.h"
#include "recon_system.h"

using namespace dvrlib;

static void setup_system(recon_system& sys) {
    sys.add_var("m_FDKeI",  46.241, 0.800);
    sys.add_var("m_FDKeII", 45.668, 0.790);
    sys.add_var("m_SpI",    44.575, 0.535);
    sys.add_var("m_SpII",   44.319, 0.532);
    sys.add_var("m_V",       0.525, 0.105);
    sys.add_var("m_HK",     69.978, 0.854);
    sys.add_var("m_A7",     10.364, 0.168);
    sys.add_var("m_A6",      3.744, 0.058);
    sys.add_var("m_A5",      4.391, 0.058);
    sys.add_var("m_HDNK",   18.498, 0.205);
    sys.add_var("m_D",       2.092, 0.272);
    sys.add_var("m_FD1",     0,    -1);
    sys.add_var("m_FD2",     0,    -1);
    sys.add_var("m_FD3",     0,    -1);
    sys.add_var("m_HDAnz",   0,    -1);
    sys.add_covariance_coeff("m_FDKeI",  "m_FDKeII", 0.2);
    sys.add_covariance_coeff("m_SpI",    "m_SpII",   0.4);
}

static matrix make_F() {
    double Fc[][15] = {
        {1, 1, 0, 0, -0.2, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0},
        {0, 0, 1, 1, -0.6, 0, 0, 0, 0, 0, 0,  0,-1, 0, 0},
        {0, 0, 0, 0,  0.4, 1, 1, 1, 1, 0, 0,  0, 0,-1, 0},
        {0, 0, 0, 0,  0,   0, 1, 1, 1, 0, 0,  0, 0, 0,-1},
        {0, 0, 0, 0,  0,   0, 0, 0, 0, 0, 0,  1,-1, 0, 0},
        {0, 0, 0, 0,  0,   0, 0, 0, 0, 0, 0,  0, 1,-1, 0},
        {0, 0, 0, 0,  0,   0, 0, 0, 0,-1, 0,  0, 0, 0, 1},
    };
    return matrix(7, 15, Fc);
}

TEST_CASE("VDI2048 regression") {
    recon_system sys;
    setup_system(sys);
    matrix S_x = sys.get_covariance_matrix();
    matrix F   = make_F();
    vector x   = sys.get_values();

    vector v(x.size());
    lin_recon(F * x, S_x, F, v);

    matrix S_v(S_x.size1(), S_x.size2());
    lin_cov_update(S_x, F, S_v);

    matrix S_x_new = S_x - S_v;
    vector x_new   = x + v;
    vector confint_new(S_x.size1());
    extract_confidence(S_x_new, confint_new);

    const double tol = 1e-6;

    SECTION("corrections v") {
        double expected[] = {
            -1.27651079357085, -1.25002785826582,
             0.350904047191697, 0.348096449940544,
             0.00134787242233519, 0.575060829315,
             0.0117829833264503, 0.00140440603423252,
             0.00140440603423252, 0.0155917953949165,
             0, 89.2771917736789, 89.2771917736789,
             89.2771917736788, 18.5135917953949
        };
        for (int i = 0; i < v.size(); i++)
            CHECK(v.get(i) == Catch::Approx(expected[i]).epsilon(tol));
    }

    SECTION("constraints satisfied after correction") {
        vector residual = F * (x + v);
        for (int i = 0; i < residual.size(); i++)
            CHECK(residual.get(i) == Catch::Approx(0.0).margin(1e-10));
    }

    SECTION("updated values x_new") {
        double expected[] = {
            44.9644892064291, 44.4179721417342,
            44.9259040471917, 44.6670964499405,
             0.526347872422335, 70.553060829315,
            10.3757829833265,  3.74540440603423,
             4.39240440603423, 18.5135917953949,
             2.092, 89.2771917736789, 89.2771917736789,
            89.2771917736788, 18.5135917953949
        };
        for (int i = 0; i < x_new.size(); i++)
            CHECK(x_new.get(i) == Catch::Approx(expected[i]).epsilon(tol));
    }

    SECTION("updated confidence intervals") {
        double expected[] = {
            0.575626570589129, 0.572817587480311,
            0.404463177421466, 0.402919370199657,
            0.104623410451298, 0.559872848891984,
            0.133003389109488, 0.056695251568643,
            0.056695251568643, 0.137102511396334,
            0.272
        };
        for (int i = 0; i < confint_new.size(); i++)
            CHECK(confint_new.get(i) == Catch::Approx(expected[i]).epsilon(tol));
    }
}
