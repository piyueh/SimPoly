/**
 * \file evaluate.cpp
 * \brief Unit tests for polynomial evaluation.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-02-01
 */


# include <valarray>
# include <complex>
# include <random>

# include <gtest/gtest.h>

# include "operations.h"
# include "exceptions.h"


extern std::default_random_engine generator;


# ifndef NDEBUG
TEST(PolynomialEvaluation, WrongLengthTests)
{
    try { simpoly::op::evaluate(nullptr, 0, 0.0); }
    catch (simpoly::exceptions::ZeroCoeffsLength e) {};
    
    try { simpoly::op::evaluate(nullptr, -1, 0.0); }
    catch (simpoly::exceptions::NegativeCoeffsLength e) {};
    
    try { simpoly::op::evaluate(std::valarray<double>(0), 0.0); }
    catch (simpoly::exceptions::ZeroCoeffsLength e) {};
    
    try { simpoly::op::evaluate(nullptr, 0, std::complex<double>(0.0)); }
    catch (simpoly::exceptions::ZeroCoeffsLength e) {};
    
    try { simpoly::op::evaluate(std::valarray<std::complex<double>>(0), 0.0); }
    catch (simpoly::exceptions::ZeroCoeffsLength e) {};
    
    try { simpoly::op::evaluate_root(1.0, nullptr, -1, 0.0); }
    catch (simpoly::exceptions::NegativeDegree e) {};
}
# endif

TEST(PolynomialEvaluation, RealPolyPointer1)
{
    double *c, *x, *expect;
    c = new double[16]{
        0.5147664494717302, 0.6377586019590598, 0.8693376634434886,
        0.3616343298342888, 0.3671257345690335, 0.9879013608892221, 
        0.057927626742414, 0.4893598360407232, 0.4563590270223324, 
        0.846309210807406, 0.3207035239642677, 0.7058682864690925, 
        0.4044647334180196, 0.8381610445132728, 0.2320406628170849, 
        0.6094284242837951};
    
    x = new double[11]{
        -0.4228109374105191, -0.8916690947701997, 0.3004964434901052,
        0.5289741804671604, 0.391125993645951, -0.9705540830999013, 
        0.140811530330192, -0.081151319191916, -0.6987384555851917, 
        -0.7946392224478713, 1.};
    
    expect = new double[11]{
        0.3708304280699384, -0.5004742070702876, 0.8003385938623538,
        1.232717689712284, 0.9378580252923071, -1.5694372775647332, 
        0.6230171267614112, 0.4685557372060891, 0.2439320727821255, 
        0.040352420677454, 8.6991465162452322};
    
    for(unsigned i=0; i<11; ++i)
    {
        double result = simpoly::op::evaluate(c, 16, x[i]);
        ASSERT_NEAR(expect[i], result, 1e-12);
    }
    
    delete[] c; delete[] x; delete[] expect;
}

TEST(PolynomialEvaluation, RealPolyPointer2)
{
    double *c, *x, *expect;
    c = new double[10]{
        0.2592584530394499, -0.7514681609411962, 0.5091058938166564,
        -0.1583353677036581, -0.482985886510811, -0.0168285696265613,
        0.1508398332186933, -0.6081712532961356, -0.5204924007984593,
        -0.8852066676326502};
    
    x = new double[11]{
        -0.7534486553162925, -0.0347858160484482, -0.736302652337975,
        -0.2999595980829783, -0.7307329134372322, -0.7458658132191887,
        -0.8773789778035546, -0.4204977672318795, 0.9370281988772768,
        0.6570947501205981, 0.1976671115913535};
    
    expect = new double[11]{
        1.157259553740849, 0.2860208900847923, 1.1201937097154451,
        0.5311058905544308, 1.1087026242468756, 1.1405386945874234,
        1.5421396366024629, 0.6642644256506487, -1.5985969097982302,
        -0.2100670419848677, 0.1286447515096008};
    
    for(unsigned i=0; i<11; ++i)
    {
        double result = simpoly::op::evaluate(c, 10, x[i]);
        ASSERT_NEAR(expect[i], result, 1e-12);
    }
    
    delete[] c; delete[] x; delete[] expect;
}

TEST(PolynomialEvaluation, RealPolyValarray1)
{
    std::valarray<double> c, x, expect; 
    c = std::valarray<double>({
        0.5147664494717302, 0.6377586019590598, 0.8693376634434886,
        0.3616343298342888, 0.3671257345690335, 0.9879013608892221, 
        0.057927626742414, 0.4893598360407232, 0.4563590270223324, 
        0.846309210807406, 0.3207035239642677, 0.7058682864690925, 
        0.4044647334180196, 0.8381610445132728, 0.2320406628170849, 
        0.6094284242837951});
    
    x = std::valarray<double>({
        -0.4228109374105191, -0.8916690947701997, 0.3004964434901052,
        0.5289741804671604, 0.391125993645951, -0.9705540830999013, 
        0.140811530330192, -0.081151319191916, -0.6987384555851917, 
        -0.7946392224478713, 1.});
    
    expect = std::valarray<double>({
        0.3708304280699384, -0.5004742070702876, 0.8003385938623538,
        1.232717689712284, 0.9378580252923071, -1.5694372775647332, 
        0.6230171267614112, 0.4685557372060891, 0.2439320727821255, 
        0.040352420677454, 8.6991465162452322});
    
    for(unsigned i=0; i<11; ++i)
    {
        double result = simpoly::op::evaluate(c, x[i]);
        ASSERT_NEAR(expect[i], result, 1e-12);
    }
}

TEST(PolynomialEvaluation, RealPolyValarray2)
{
    std::valarray<double> c, x, expect; 
    c = std::valarray<double>({
        0.2592584530394499, -0.7514681609411962, 0.5091058938166564,
        -0.1583353677036581, -0.482985886510811, -0.0168285696265613,
        0.1508398332186933, -0.6081712532961356, -0.5204924007984593,
        -0.8852066676326502});
    
    x = std::valarray<double>({
        -0.7534486553162925, -0.0347858160484482, -0.736302652337975,
        -0.2999595980829783, -0.7307329134372322, -0.7458658132191887,
        -0.8773789778035546, -0.4204977672318795, 0.9370281988772768,
        0.6570947501205981, 0.1976671115913535});
    
    expect = std::valarray<double>({
        1.157259553740849, 0.2860208900847923, 1.1201937097154451,
        0.5311058905544308, 1.1087026242468756, 1.1405386945874234,
        1.5421396366024629, 0.6642644256506487, -1.5985969097982302,
        -0.2100670419848677, 0.1286447515096008});
    
    for(unsigned i=0; i<11; ++i)
    {
        double result = simpoly::op::evaluate(c, x[i]);
        ASSERT_NEAR(expect[i], result, 1e-12);
    }
}

TEST(PolynomialEvaluation, ComplexPolyPointer)
{
    std::complex<double> *c, *x, *expect;
    c = new std::complex<double>[8]{
        std::complex<double>(0.759414027655581, 0.5692286126871828),
        std::complex<double>(-0.3909314128950405, -0.5208053215270787),
        std::complex<double>(-0.4159999364898488, -0.929797800647657),
        std::complex<double>(-0.05289305288679613, 0.7623008176129775),
        std::complex<double>(0.3810204528687382, 0.5443453788782495),
        std::complex<double>(-0.182901423575166, 0.5832026946660676),
        std::complex<double>(0.5132115703426237, -0.2746174894763567),
        std::complex<double>(-0.05772567034291853, 0.2078453112062857)};

    x = new std::complex<double>[5]{
        std::complex<double>(0.5022381178753841, 0.6970911204735433),
        std::complex<double>(-0.2660269666516346, 0.3408308392987081),
        std::complex<double>(0.07225847489108794, -0.4890141431874784),
        std::complex<double>(-0.6984620632833527, 0.3385783677682219),
        std::complex<double>(0.6724492933394415, 0.6822524332030151)};
    
    expect = new std::complex<double>[5]{
        std::complex<double>(1.910143035375531, -0.9284375830984708),
        std::complex<double>(0.8535065511576807, 0.7411792224928896),
        std::complex<double>(0.4293431036176729, 0.9751281975153519),
        std::complex<double>(0.1349251253981901, 0.4761603912848932),
        std::complex<double>(1.294817005520558, -1.875701442672196)};
    
    for(unsigned i=0; i<5; ++i)
    {
        std::complex<double> result = simpoly::op::evaluate(c, 8, x[i]);
        ASSERT_NEAR(expect[i].real(), result.real(), 1e-12);
        ASSERT_NEAR(expect[i].imag(), result.imag(), 1e-12);
    }
    
    delete[] c; delete[] x; delete[] expect;
}

TEST(PolynomialEvaluation, ComplexPolyValarray)
{
    std::valarray<std::complex<double>> 
        
        c({std::complex<double>(0.759414027655581, 0.5692286126871828),
            std::complex<double>(-0.3909314128950405, -0.5208053215270787),
            std::complex<double>(-0.4159999364898488, -0.929797800647657),
            std::complex<double>(-0.05289305288679613, 0.7623008176129775),
            std::complex<double>(0.3810204528687382, 0.5443453788782495),
            std::complex<double>(-0.182901423575166, 0.5832026946660676),
            std::complex<double>(0.5132115703426237, -0.2746174894763567),
            std::complex<double>(-0.05772567034291853, 0.2078453112062857)}), 
        
        x({std::complex<double>(0.5022381178753841, 0.6970911204735433),
            std::complex<double>(-0.2660269666516346, 0.3408308392987081),
            std::complex<double>(0.07225847489108794, -0.4890141431874784),
            std::complex<double>(-0.6984620632833527, 0.3385783677682219),
            std::complex<double>(0.6724492933394415, 0.6822524332030151)}), 
        
        expect({std::complex<double>(1.910143035375531, -0.9284375830984708),
            std::complex<double>(0.8535065511576807, 0.7411792224928896),
            std::complex<double>(0.4293431036176729, 0.9751281975153519),
            std::complex<double>(0.1349251253981901, 0.4761603912848932),
            std::complex<double>(1.294817005520558, -1.875701442672196)});

    for(unsigned i=0; i<5; ++i)
    {
        std::complex<double> result = simpoly::op::evaluate(c, x[i]);
        ASSERT_NEAR(expect[i].real(), result.real(), 1e-12);
        ASSERT_NEAR(expect[i].imag(), result.imag(), 1e-12);
    }
}

TEST(PolynomialEvaluation, RealRootPolyPointer)
{
    double l, *r, *x, *expect;
    
    l = 5.0;
    
    r = new double[13]{
        0.7634365506829128, 0.8181683656769851, 0.6127459448202521,
        -0.8553740751871957, -0.118458881248837, 0.2596055040782295,
        0.4167095094451332, -0.2111594728172452, 0.8450265251779554,
        -0.0581022083352503, 0.8087840971799483, -0.9677953376811022,
        0.6071235649521618};
    
    x = new double[11]{-1., -0.8, -0.6, -0.4, -0.2,  0., 0.2, 0.4, 0.6, 0.8, 1.}; 
    
    expect = new double[11]{
        -7.5489115691444941e-01, -2.3735165519244153e-01, -2.3998048722027623e-01,
        -2.7540763589693359e-02, 6.2463983606482018e-05, 1.0334842643644593e-04,
        6.1825068446831833e-05, -3.5617200624308186e-06, 4.5214620084522330e-08,
        -2.2918937481954420e-08, 2.1914364795650341e-03};
    
    for(unsigned i=0; i<11; ++i)
    {
        double result = simpoly::op::evaluate_root(l, r, 13, x[i]);
        ASSERT_NEAR(expect[i], result, 1e-12);
    }
    
    delete[] r; delete[] x; delete[] expect;
}

TEST(PolynomialEvaluation, RealRootPolyValarray)
{
    double l = 5.0;
    
    std::valarray<double>
        
        r({0.7634365506829128, 0.8181683656769851, 0.6127459448202521,
            -0.8553740751871957, -0.118458881248837, 0.2596055040782295,
            0.4167095094451332, -0.2111594728172452, 0.8450265251779554,
            -0.0581022083352503, 0.8087840971799483, -0.9677953376811022,
            0.6071235649521618}),
        
        x({-1., -0.8, -0.6, -0.4, -0.2,  0., 0.2, 0.4, 0.6, 0.8, 1.}),
        
        expect({
            -7.5489115691444941e-01, -2.3735165519244153e-01, -2.3998048722027623e-01,
            -2.7540763589693359e-02, 6.2463983606482018e-05, 1.0334842643644593e-04,
            6.1825068446831833e-05, -3.5617200624308186e-06, 4.5214620084522330e-08,
            -2.2918937481954420e-08, 2.1914364795650341e-03});
    
    for(unsigned i=0; i<11; ++i)
    {
        double result = simpoly::op::evaluate_root(l, r, x[i]);
        ASSERT_NEAR(expect[i], result, 1e-12);
    }
}
