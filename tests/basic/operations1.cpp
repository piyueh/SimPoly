/**
 * \file operations1.cpp
 * \brief Unit tests for polynomial basic operations.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-01-29
 */


# include <random>

# include <gtest/gtest.h>

# include "basic.h"


extern std::default_random_engine generator;

using namespace simpoly;


TEST(PolynomialCoeffs, GetCoeffsFromRoots1)
{
    basic::DArry roots({
            0.8538721859301963146,  0.2921220802626576241,
            -0.9125881264026460826,  0.2467728385548548786,
            -0.4448851063326948463, -0.6172124325075281082,
            -0.6301262014641872966});
    
    basic::DArry expect({
            -0.0097194008368678246,  0.0203713227801400679,
            0.1774520322384245952, -0.1187869060899048868,
            -1.1062246986215185096, -0.6075137057194068824,
            1.2120447619593475164,  1.});

    basic::DArry result = basic::to_coefficients(1.0, roots);

    ASSERT_EQ(expect.size(), result.size());
    
    for(unsigned i=0; i<expect.size(); ++i)
        ASSERT_NEAR(expect[i], result[i], 1e-10);
}


static auto f = [](basic::Cmplx i, basic::Cmplx j)->bool {
    if (std::abs(std::abs(i)-std::abs(j)) <= 1e-11)
    {
        if (std::abs(i.real()-j.real()) <= 1e-11)
            if (i.imag() <= j.imag()) return true;
        
        if ((i.real()-j.real()) < -1e-11) return true;
    }
    if ((std::abs(i)-std::abs(j)) < -1e-11) return true;
    return false;
};
TEST(PolynomialCoeffs, GetCoeffsFromRoots2)
{
    basic::CArry roots({
            basic::Cmplx(0.0523249771240981154, 0.7007289936213532133),
            basic::Cmplx(-0.1511969052827171200, 0.4108412398766685403),
            basic::Cmplx(0.7769313311349028695, 0.7301249050268194818),
            basic::Cmplx(0.0046666287791925498, 0.9696847924026805732),
            basic::Cmplx(-0.5789703445402629711, -0.6424983236922838792),
            basic::Cmplx(-0.8725637786525710649, -0.2471912341041830885),
            basic::Cmplx(-0.8482223096380121508, 0.1376940842539751131),
            basic::Cmplx(0.0332852932796565693, -0.8528866044880742958),
            basic::Cmplx(0.8522720871683582278, 0.1953603750595562794),
            basic::Cmplx(-0.8680725141564455960, 0.423548363231143199) ,
            basic::Cmplx(0.8824640585760064049, -0.7545152381037976141),
            basic::Cmplx(0.8001160654174328535, 0.5327952053290956336),
            basic::Cmplx(0.7524249058508991617, -0.4850694648111775198),
            basic::Cmplx(0.4670358901871698443, 0.7246115099814800509),
            basic::Cmplx(0.6425671238360994852, -0.8227184096042670092),
            basic::Cmplx(-0.8746249058224726536, 0.4526170928635313562),
            basic::Cmplx(-0.8405665222993630170, 0.4624841604488800773)});
    
    basic::CArry expect({
            basic::Cmplx(0.1284473674676091248, 0.027106308015679835) ,
            basic::Cmplx(-0.0055811935483588704, 0.5895676656919293102),
            basic::Cmplx(-1.0697675991580708477, 0.0739888122261974646),
            basic::Cmplx(0.0350214941923967515, -1.2581110962554780119),
            basic::Cmplx(1.1311480794477093426, 0.5697505628564978331),
            basic::Cmplx(-0.9578963956878190356, 0.8762473636503139707),
            basic::Cmplx(-0.7142421431728066539, -0.9120235859686047775),
            basic::Cmplx(0.6191826581229090687, 0.1287334984037127317),
            basic::Cmplx(-1.0386382947640331320, 0.6191766562897030868),
            basic::Cmplx(-1.4622204324047241020, -2.0911509295099053318),
            basic::Cmplx(2.5643980136130921643, -1.3562694162887596327),
            basic::Cmplx(0.6477469434997339404, 3.3575574217272889932),
            basic::Cmplx(-4.6940592271922572110, -0.9211194673159723845),
            basic::Cmplx(1.7891645824502104833, -4.8182053599055656434),
            basic::Cmplx(3.5022002736720923011, 2.2497280480793735435),
            basic::Cmplx(-2.7137599509433454337, 1.2366520250938943803),
            basic::Cmplx(-0.2298710809619715079, -1.9356114472913998892),
            basic::Cmplx(1.0000000000000000000, 0.)});

    basic::CArry result = basic::to_coefficients(basic::Cmplx(1.0), roots);

    ASSERT_EQ(expect.size(), result.size());
    
    for(unsigned i=0; i<expect.size(); ++i)
    {
        ASSERT_NEAR(expect[i].real(), result[i].real(), 1e-10);
        ASSERT_NEAR(expect[i].imag(), result[i].imag(), 1e-10);
    }
}


TEST(PolynomialAdd, EqualLength)
{
    std::uniform_int_distribution<int> dist1(1, 20);
    std::uniform_real_distribution<double> dist2(-5.0, 5.0);
    
    int len = dist1(generator);
    
    std::valarray<double> p1(len), p2(len), result(len);
    
    for(unsigned i=0; i<len; ++i)
    {
        p1[i] = dist2(generator);
        p2[i] = dist2(generator);
    }
    
    result = simpoly::basic::add(p1, p2);
    
    for(unsigned i=0; i<len; ++i)
        ASSERT_NEAR(p1[i]+p2[i], result[i], 1e-12);
}

TEST(PolynomialAdd, FirstOneLonger)
{
    std::uniform_int_distribution<int> dist1(1, 20);
    
    int len = dist1(generator);
    std::valarray<double> p1(len), result(len);
    
    std::uniform_int_distribution<int> dist2(0, len-1);
    std::valarray<double> p2(dist2(generator));
    
    std::uniform_real_distribution<double> rdouble(-5.0, 5.0);
    for(unsigned i=0; i<len; ++i) p1[i] = rdouble(generator);
    for(unsigned i=0; i<p2.size(); ++i) p2[i] = rdouble(generator);
    
    result = simpoly::basic::add(p1, p2); 
    for(unsigned i=0; i<p2.size(); ++i) ASSERT_NEAR(p1[i]+p2[i], result[i], 1e-12);
    for(unsigned i=p2.size(); i<len; ++i) ASSERT_NEAR(p1[i], result[i], 1e-12);
}

TEST(PolynomialAdd, SecondOneLonger)
{
    std::uniform_int_distribution<int> dist1(1, 20);
    
    int len = dist1(generator);
    std::valarray<double> p1(len), result(len);
    
    std::uniform_int_distribution<int> dist2(0, len-1);
    std::valarray<double> p2(dist2(generator));
    
    std::uniform_real_distribution<double> rdouble(-5.0, 5.0);
    for(unsigned i=0; i<len; ++i) p1[i] = rdouble(generator);
    for(unsigned i=0; i<p2.size(); ++i) p2[i] = rdouble(generator);
    
    result = simpoly::basic::add(p2, p1); 
    for(unsigned i=0; i<p2.size(); ++i) ASSERT_NEAR(p1[i]+p2[i], result[i], 1e-12);
    for(unsigned i=p2.size(); i<len; ++i) ASSERT_NEAR(p1[i], result[i], 1e-12);
}


TEST(PolynomialGCD, GCD1)
{
    basic::DArry p1({
        6.48000000e+02,   5.40000000e+02,  -8.10000000e+02,
        -1.11900000e+03,  -1.48000000e+02,   4.43000000e+02,
        3.26000000e+02,   1.03000000e+02,   1.60000000e+01,
        1.00000000e+00});
    
    basic::DArry q1_expect({-108., -108., 45., 104., 54., 12., 1.});
    basic::DArry q2_expect({18., 21., 8., 1.});
    basic::DArry q3_expect({3., 1.});
    basic::DArry q4_expect({1.});
    
    basic::DArry q1 = basic::GCD(p1, basic::derivative(p1));
    basic::DArry q2 = basic::GCD(q1, basic::derivative(q1));
    basic::DArry q3 = basic::GCD(q2, basic::derivative(q2));
    basic::DArry q4 = basic::GCD(q3, basic::derivative(q3));
    
    for(unsigned i=0; i<q1_expect.size(); ++i)
        ASSERT_NEAR(q1_expect[i], q1[i], 1e-8);
    
    for(unsigned i=0; i<q2_expect.size(); ++i)
        ASSERT_NEAR(q2_expect[i], q2[i], 1e-8);
    
    for(unsigned i=0; i<q3_expect.size(); ++i)
        ASSERT_NEAR(q3_expect[i], q3[i], 1e-8);
    
    for(unsigned i=0; i<q4_expect.size(); ++i)
        ASSERT_NEAR(q4_expect[i], q4[i], 1e-8);
    
}

TEST(PolynomialGCD, GCD2)
{
    basic::DArry p1({
            -7.5202173502138876682e-04,  -8.6835988738648978158e-03,
            5.5348919281959091387e-02,   1.4315233176654140745e-01,
            -4.3769681353357492437e-01,  -1.2069048603868517411e+00,
            -2.8638873139566145554e-02,   1.6438835517368903805e+00,
            1.0000000000000000000e+00});
    
    basic::DArry q1_expect({-0.2431327906295455976, 1.0});
    basic::DArry q2_expect({1.});
    
    basic::DArry q1 = basic::GCD(p1, basic::derivative(p1));
    basic::DArry q2 = basic::GCD(q1, basic::derivative(q1));
    
    for(unsigned i=0; i<q1_expect.size(); ++i)
        ASSERT_NEAR(q1_expect[i], q1[i], 1e-8);
    
    for(unsigned i=0; i<q2_expect.size(); ++i)
        ASSERT_NEAR(q2_expect[i], q2[i], 1e-8);
    
}

TEST(PolynomialGCD, GCD3)
{
    basic::DArry p1({
            0.0000000000000000000e+00,  -2.9491200000000012543e-02,
            8.8473600000000013344e-02,   9.6112640000000015839e-01,
            -3.1193088000000002147e+00,  -5.5872000000000010544e+00,
            2.5158400000000003871e+01,  -2.2080000000000006288e+00,
            -6.3264000000000010004e+01,   5.2000000000000014211e+01,
            3.5999999999999992895e+01,  -6.0000000000000000000e+01,
            2.0000000000000000000e+01});
    
    basic::DArry q1_expect({1.0, -2.0, 1.0});
    basic::DArry q2_expect({-1.0, 1.0});
    basic::DArry q3_expect({1.});
    
    basic::DArry q1 = basic::GCD(p1, basic::derivative(p1));
    basic::DArry q2 = basic::GCD(q1, basic::derivative(q1));
    basic::DArry q3 = basic::GCD(q2, basic::derivative(q2));
    
    for(unsigned i=0; i<q1_expect.size(); ++i)
        ASSERT_NEAR(q1_expect[i], q1[i], 1e-8);
    
    for(unsigned i=0; i<q2_expect.size(); ++i)
        ASSERT_NEAR(q2_expect[i], q2[i], 1e-8);
    
    for(unsigned i=0; i<q3_expect.size(); ++i)
        ASSERT_NEAR(q3_expect[i], q3[i], 1e-8);
    
}
