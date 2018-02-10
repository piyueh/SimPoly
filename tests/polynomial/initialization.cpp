/**
 * \file tests/polynomial/initialization.cpp
 * \brief Unit tests for initialization of Polynomial class.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-02-08
 */

# include <algorithm>
# include <cmath>
# include <gtest/gtest.h>

# include "exceptions.h"
# include "polynomial.h"

using namespace simpoly;

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

void run(basic::DArry &c, basic::DArry &rr, basic::CArry &cr, basic::CArry &allr,
        const poly::Polynomial &p, poly::PolyType type, bool use_r, unsigned n,
        unsigned nr, unsigned nc)
{
    basic::DArry pc = p.coef();
    basic::DArry prr = p.real_roots();
    basic::CArry pcr = p.cmplx_roots();
    basic::CArry pallr = p.roots();

    std::sort(rr.begin(), rr.end());
    std::sort(prr.begin(), prr.end());
    std::sort(cr.begin(), cr.end(), f);
    std::sort(pcr.begin(), pcr.end(), f);
    std::sort(allr.begin(), allr.end(), f);
    std::sort(pallr.begin(), pallr.end(), f);

    ASSERT_EQ(type, p.type());
    ASSERT_EQ(use_r, p.use_roots());
    ASSERT_EQ(n, p.degree());
    ASSERT_EQ(nr, p.n_real_roots());
    ASSERT_EQ(nc, p.n_cmplx_roots());

    ASSERT_EQ(c.size(), pc.size());
    for(unsigned i=0; i<c.size(); ++i)
        ASSERT_NEAR(c[i], pc[i], 1e-10);

    ASSERT_EQ(rr.size(), prr.size());
    for(unsigned i=0; i<rr.size(); ++i)
        ASSERT_NEAR(rr[i], prr[i], 1e-10);

    ASSERT_EQ(cr.size(), pcr.size());
    for(unsigned i=0; i<cr.size(); ++i)
    {
        ASSERT_NEAR(cr[i].real(), pcr[i].real(), 1e-10);
        ASSERT_NEAR(cr[i].imag(), pcr[i].imag(), 1e-10);
    }

    ASSERT_EQ(allr.size(), pallr.size());
    for(unsigned i=0; i<allr.size(); ++i)
    {
        ASSERT_NEAR(allr[i].real(), pallr[i].real(), 1e-10);
        ASSERT_NEAR(allr[i].imag(), pallr[i].imag(), 1e-10);
    }
}

TEST(PolynomialConstruct, Coefficient1)
{
    basic::DArry c({
            -0.1195107918695455,  0.8447947406345584, -0.6030868103450979,
            -0.9766533913297424, -0.2677767116634977,  0.3804660521221173,
            -0.4716581713021364, -0.8085712194756074,  0.1556975750852925,
            -0.6073461984864883,  0.8943598261002668, -0.0478803277253317,
            -0.5541643941753234, -0.008585713188056 ,  0.9873624852167924,
            0.1332436172903213 });

    basic::DArry rr({
            -7.3453301351775968, -0.9283764426295956, 0.1669466647980708,
            0.5160978326951132, 1.0922685323108878});

    basic::CArry cr({
            basic::Cmplx(-1.0053226251035317, -0.6685050805086484),
            basic::Cmplx(-1.0053226251035317, 0.6685050805086484),
            basic::Cmplx(-0.6566634869174144, -0.6613997230665936),
            basic::Cmplx(-0.6566634869174144, 0.6613997230665936),
            basic::Cmplx(-0.1313113999816178, -1.013006128623221),
            basic::Cmplx(-0.1313113999816178, 1.013006128623221),
            basic::Cmplx(0.5190399155424927, -0.8510559355096275),
            basic::Cmplx(0.5190399155424927, 0.8510559355096275),
            basic::Cmplx(0.8183519974697149, -0.62842193876576),
            basic::Cmplx(0.8183519974697149, 0.62842193876576)});

    basic::CArry allroots({
            basic::Cmplx(-7.3453301351775968, 0.),
            basic::Cmplx(-1.0053226251035317, -0.6685050805086484),
            basic::Cmplx(-1.0053226251035317, 0.6685050805086484),
            basic::Cmplx(-0.9283764426295956, 0.),
            basic::Cmplx(-0.6566634869174144, -0.6613997230665936),
            basic::Cmplx(-0.6566634869174144, 0.6613997230665936),
            basic::Cmplx(-0.1313113999816178, -1.013006128623221),
            basic::Cmplx(-0.1313113999816178, 1.013006128623221),
            basic::Cmplx(0.1669466647980708, 0.),
            basic::Cmplx(0.5160978326951132, 0.),
            basic::Cmplx(0.5190399155424927, -0.8510559355096275),
            basic::Cmplx(0.5190399155424927, 0.8510559355096275),
            basic::Cmplx(0.8183519974697149, -0.62842193876576),
            basic::Cmplx(0.8183519974697149, 0.62842193876576),
            basic::Cmplx(1.0922685323108878, 0.) });

    poly::Polynomial p(c);

    run(c, rr, cr, allroots, p, poly::General, false, 15, 5, 10);
}

TEST(PolynomialConstruct, Coefficient2)
{
    basic::DArry c({
            0.0254213669105730726, -0.2534500155209851724, 0.8062095562787726655,
            -0.6357752008483941175, -1.5610318058296166477, 3.2409361412883592735,
            -0.7130473867752948891, -3.1340615462415519765, 3.3223860099444975091,
            0.9010891620340866215, -8.1217453795715854881, 12.2404621022169397548,
            -8.2173667074445813086, 2.1000000000000000000});

    basic::DArry rr({
            0.7141171158689039178,  0.2116130677773324642, -0.6427357434615013076,
            0.99051581546728662  , 0.6675917803331954392,  0.92450641235501263  ,
            0.5897280649985707779});

    basic::CArry cr({
            basic::Cmplx(0.0470049721208805060, -0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, -0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, 0.3828854883388757013),
            basic::Cmplx(0.0470049721208805060, 0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, 0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, -0.3828854883388757013)});

    basic::CArry allr({
            basic::Cmplx(0.7141171158689039178, 0.0),
            basic::Cmplx(0.2116130677773324642, 0.0),
            basic::Cmplx(-0.6427357434615013076, 0.0),
            basic::Cmplx(0.99051581546728662, 0.0),
            basic::Cmplx(0.6675917803331954392, 0.0),
            basic::Cmplx(0.92450641235501263, 0.0),
            basic::Cmplx(0.5897280649985707779, 0.0),
            basic::Cmplx(0.0470049721208805060, -0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, -0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, 0.3828854883388757013),
            basic::Cmplx(0.0470049721208805060, 0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, 0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, -0.3828854883388757013)});

    poly::Polynomial p(c);

    run(c, rr, cr, allr, p, poly::General, false, 13, 7, 6);
}

TEST(PolynomialConstruct, RealRoots)
{
    basic::DArry c({
        -0.0030656367973936837, -0.0546778459167431147, -0.1104907299430930162,
        0.9134673653237410651, 0.232604797603839647 , -2.2719417419792904411,
        -0.096235546433642627 ,  1.5});

    basic::DArry rr({
        -0.793002440830254951 , -0.1450108736699442513, 0.89637937289469094,
        -0.943023789745914609 , 0.3940559055598014115, -0.0732864366680110457,
        0.7280452934153942568});

    basic::CArry cr({});
    basic::CArry allr({
        -0.793002440830254951 , -0.1450108736699442513, 0.89637937289469094,
        -0.943023789745914609 , 0.3940559055598014115, -0.0732864366680110457,
        0.7280452934153942568});

    poly::Polynomial p(1.5, rr);

    run(c, rr, cr, allr, p, poly::General, true, 7, 7, 0);
}

TEST(PolynomialConstruct, CmplxRoots)
{
    basic::DArry c({
            0.0294157315584267429, 0.1828545566550552615, 0.7752484826061304402,
            2.1749217360894204099, 4.7522887221290321946, 7.8363735137657037910,
            10.2309579575511087057, 9.9964047794244574874, 7.3390840252023901158,
            3.4616726139023166198, 1.0000000000000000000});

    basic::DArry rr({});

    basic::CArry cr({
            basic::Cmplx(-0.2210988873678600974, -0.4981159220444357771),
            basic::Cmplx(-0.6320197898502710743, -0.9979011634849095369),
            basic::Cmplx(-0.0868288229751668084, -0.5910951393758168138),
            basic::Cmplx(-0.2514232122885062815, 0.5350583485005697693),
            basic::Cmplx(-0.5394655944693540484, 0.5272551705386470111),
            basic::Cmplx(-0.2210988873678600974, 0.4981159220444357771),
            basic::Cmplx(-0.6320197898502710743, 0.9979011634849095369),
            basic::Cmplx(-0.0868288229751668084, 0.5910951393758168138),
            basic::Cmplx(-0.2514232122885062815, -0.5350583485005697693),
            basic::Cmplx(-0.5394655944693540484, -0.5272551705386470111)});

    basic::CArry allr({
            basic::Cmplx(-0.2210988873678600974, -0.4981159220444357771),
            basic::Cmplx(-0.6320197898502710743, -0.9979011634849095369),
            basic::Cmplx(-0.0868288229751668084, -0.5910951393758168138),
            basic::Cmplx(-0.2514232122885062815, 0.5350583485005697693),
            basic::Cmplx(-0.5394655944693540484, 0.5272551705386470111),
            basic::Cmplx(-0.2210988873678600974, 0.4981159220444357771),
            basic::Cmplx(-0.6320197898502710743, 0.9979011634849095369),
            basic::Cmplx(-0.0868288229751668084, 0.5910951393758168138),
            basic::Cmplx(-0.2514232122885062815, -0.5350583485005697693),
            basic::Cmplx(-0.5394655944693540484, -0.5272551705386470111)});

    poly::Polynomial p(1.0, cr);

    run(c, rr, cr, allr, p, poly::General, false, 10, 0, 10);
}

TEST(PolynomialConstruct, MixedRoots)
{
    basic::DArry c({
            4.8761816514376562538e-04, -3.0680255180199938331e-02, -1.9399430784745028133e-01,
            2.2327406840241159625e-01, 8.5211236507700105491e-01, -7.1470499725845870742e-01,
            -1.8098649104020723133e+00, 7.2797782400452593077e-01, 1.9871178294829077160e+00,
            -8.4472365176380781460e-01, -8.8108740335678636946e-01, 1.0000000000000000000e+00});

    basic::DArry rr({
            0.0145742365199481583, -0.1634630878288063638, -0.6482642852285591673,
            0.8651125370780239354, -0.8332453037107427907});

    basic::CArry cr({
            basic::Cmplx(0.5716851193463614589, 0.2604540748186188193),
            basic::Cmplx(0.8863017933854371311, 0.952203188776464815),
            basic::Cmplx(-0.6348002594683372912, 0.5028578564320600464),
            basic::Cmplx(0.5716851193463614589, -0.2604540748186188193),
            basic::Cmplx(0.8863017933854371311, -0.952203188776464815),
            basic::Cmplx(-0.6348002594683372912, -0.5028578564320600464)});

    basic::CArry allr({
            basic::Cmplx(0.0145742365199481583, 0.0),
            basic::Cmplx(-0.1634630878288063638, 0.0),
            basic::Cmplx(-0.6482642852285591673, 0.0),
            basic::Cmplx(0.8651125370780239354, 0.0),
            basic::Cmplx(-0.8332453037107427907, 0.0),
            basic::Cmplx(0.5716851193463614589, 0.2604540748186188193),
            basic::Cmplx(0.8863017933854371311, 0.952203188776464815),
            basic::Cmplx(-0.6348002594683372912, 0.5028578564320600464),
            basic::Cmplx(0.5716851193463614589, -0.2604540748186188193),
            basic::Cmplx(0.8863017933854371311, -0.952203188776464815),
            basic::Cmplx(-0.6348002594683372912, -0.5028578564320600464)});

    poly::Polynomial p(1.0, rr, cr);

    run(c, rr, cr, allr, p, poly::General, false, 11, 5, 6);
}

TEST(PolynomialConstruct, CoefficientRealRoots)
{
    basic::DArry c({
        8.8116801897574424103e-08,  -1.5755432076005173523e-05,
        3.2755888134778607694e-04,   1.0390721688220840024e-02,
        3.0402292514555038366e-02,  -1.0868964479964349223e-01,
        -3.5559953729810578338e-01,   4.0262150963981052643e-01,
        1.3375567521249724923e+00,  -6.4769378204836203228e-01,
        -1.9747141187711405586e+00,   3.7393558660794212756e-01,
        1.0000000000000000000e+00});

    basic::DArry rr({
        -0.4511966608803148482,  0.0067448609489018096,
        0.0207296190229941946, -0.9665704195214050998,
        0.8617202469880191895, -0.0723920445878822871,
        -0.8553373479037251759, -0.2703413922920632206,
        -0.4764932045263101656,  0.5856198013651621181,
        0.7879411424312969814,  0.4556398123473843764});

    basic::CArry cr({});

    basic::CArry allr({
        -0.4511966608803148482,  0.0067448609489018096,
        0.0207296190229941946, -0.9665704195214050998,
        0.8617202469880191895, -0.0723920445878822871,
        -0.8553373479037251759, -0.2703413922920632206,
        -0.4764932045263101656,  0.5856198013651621181,
        0.7879411424312969814,  0.4556398123473843764});

    poly::Polynomial p(c, rr);

    run(c, rr, cr, allr, p, poly::General, true, 12, 12, 0);
}

TEST(PolynomialConstruct, CoefficientCmplxRoots)
{
    basic::DArry c({
            0.8518872382911975016,  1.1637936035152987024, 2.5451161666633632663,
            2.5567555527462921816, 3.3173821473742890475,  2.1100889511742426663,
            2.4532504244953017292,  1.5634925815615026146, 1.0000000000000000000});

    basic::DArry rr({});

    basic::CArry cr({
            basic::Cmplx(-0.8016639666671099285, 0.8357146561432540199),
            basic::Cmplx(0.5924357957191641599, 0.8630938655638200618),
            basic::Cmplx(-0.6149591360259734074, -0.7153509672829487798),
            basic::Cmplx(0.0424410161931678687, 0.8059405380525002105),
            basic::Cmplx(-0.8016639666671099285, -0.8357146561432540199),
            basic::Cmplx(0.5924357957191641599, -0.8630938655638200618),
            basic::Cmplx(-0.6149591360259734074, 0.7153509672829487798),
            basic::Cmplx(0.0424410161931678687, -0.8059405380525002105)});

    basic::CArry allr({
            basic::Cmplx(-0.8016639666671099285, 0.8357146561432540199),
            basic::Cmplx(0.5924357957191641599, 0.8630938655638200618),
            basic::Cmplx(-0.6149591360259734074, -0.7153509672829487798),
            basic::Cmplx(0.0424410161931678687, 0.8059405380525002105),
            basic::Cmplx(-0.8016639666671099285, -0.8357146561432540199),
            basic::Cmplx(0.5924357957191641599, -0.8630938655638200618),
            basic::Cmplx(-0.6149591360259734074, 0.7153509672829487798),
            basic::Cmplx(0.0424410161931678687, -0.8059405380525002105)});

    poly::Polynomial p(c, cr);

    run(c, rr, cr, allr, p, poly::General, false, 8, 0, 8);
}

TEST(PolynomialConstruct, CoefficientMixedRoots)
{
    basic::DArry c({
            0.0254213669105730726, -0.2534500155209851724, 0.8062095562787726655,
            -0.6357752008483941175, -1.5610318058296166477, 3.2409361412883592735,
            -0.7130473867752948891, -3.1340615462415519765, 3.3223860099444975091,
            0.9010891620340866215, -8.1217453795715854881, 12.2404621022169397548,
            -8.2173667074445813086, 2.1000000000000000000});

    basic::DArry rr({
            0.7141171158689039178,  0.2116130677773324642, -0.6427357434615013076,
            0.99051581546728662  , 0.6675917803331954392,  0.92450641235501263  ,
            0.5897280649985707779});

    basic::CArry cr({
            basic::Cmplx(0.0470049721208805060, -0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, -0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, 0.3828854883388757013),
            basic::Cmplx(0.0470049721208805060, 0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, 0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, -0.3828854883388757013)});

    basic::CArry allr({
            basic::Cmplx(0.7141171158689039178, 0.0),
            basic::Cmplx(0.2116130677773324642, 0.0),
            basic::Cmplx(-0.6427357434615013076, 0.0),
            basic::Cmplx(0.99051581546728662, 0.0),
            basic::Cmplx(0.6675917803331954392, 0.0),
            basic::Cmplx(0.92450641235501263, 0.0),
            basic::Cmplx(0.5897280649985707779, 0.0),
            basic::Cmplx(0.0470049721208805060, -0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, -0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, 0.3828854883388757013),
            basic::Cmplx(0.0470049721208805060, 0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, 0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, -0.3828854883388757013)});

    poly::Polynomial p(c, rr, cr);

    run(c, rr, cr, allr, p, poly::General, false, 13, 7, 6);
}

TEST(PolynomialConstruct, MoveConstruction)
{
    basic::DArry c({
            0.0254213669105730726, -0.2534500155209851724, 0.8062095562787726655,
            -0.6357752008483941175, -1.5610318058296166477, 3.2409361412883592735,
            -0.7130473867752948891, -3.1340615462415519765, 3.3223860099444975091,
            0.9010891620340866215, -8.1217453795715854881, 12.2404621022169397548,
            -8.2173667074445813086, 2.1000000000000000000});

    basic::DArry rr({
            0.7141171158689039178,  0.2116130677773324642, -0.6427357434615013076,
            0.99051581546728662  , 0.6675917803331954392,  0.92450641235501263  ,
            0.5897280649985707779});

    basic::CArry cr({
            basic::Cmplx(0.0470049721208805060, -0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, -0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, 0.3828854883388757013),
            basic::Cmplx(0.0470049721208805060, 0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, 0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, -0.3828854883388757013)});

    basic::CArry allr({
            basic::Cmplx(0.7141171158689039178, 0.0),
            basic::Cmplx(0.2116130677773324642, 0.0),
            basic::Cmplx(-0.6427357434615013076, 0.0),
            basic::Cmplx(0.99051581546728662, 0.0),
            basic::Cmplx(0.6675917803331954392, 0.0),
            basic::Cmplx(0.92450641235501263, 0.0),
            basic::Cmplx(0.5897280649985707779, 0.0),
            basic::Cmplx(0.0470049721208805060, -0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, -0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, 0.3828854883388757013),
            basic::Cmplx(0.0470049721208805060, 0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, 0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, -0.3828854883388757013)});

    poly::Polynomial p(std::move(poly::Polynomial(c)));

    run(c, rr, cr, allr, p, poly::General, false, 13, 7, 6);
}

TEST(PolynomialConstruct, CopyConstruction)
{
    basic::DArry c({
            0.0254213669105730726, -0.2534500155209851724, 0.8062095562787726655,
            -0.6357752008483941175, -1.5610318058296166477, 3.2409361412883592735,
            -0.7130473867752948891, -3.1340615462415519765, 3.3223860099444975091,
            0.9010891620340866215, -8.1217453795715854881, 12.2404621022169397548,
            -8.2173667074445813086, 2.1000000000000000000});

    basic::DArry rr({
            0.7141171158689039178,  0.2116130677773324642, -0.6427357434615013076,
            0.99051581546728662  , 0.6675917803331954392,  0.92450641235501263  ,
            0.5897280649985707779});

    basic::CArry cr({
            basic::Cmplx(0.0470049721208805060, -0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, -0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, 0.3828854883388757013),
            basic::Cmplx(0.0470049721208805060, 0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, 0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, -0.3828854883388757013)});

    basic::CArry allr({
            basic::Cmplx(0.7141171158689039178, 0.0),
            basic::Cmplx(0.2116130677773324642, 0.0),
            basic::Cmplx(-0.6427357434615013076, 0.0),
            basic::Cmplx(0.99051581546728662, 0.0),
            basic::Cmplx(0.6675917803331954392, 0.0),
            basic::Cmplx(0.92450641235501263, 0.0),
            basic::Cmplx(0.5897280649985707779, 0.0),
            basic::Cmplx(0.0470049721208805060, -0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, -0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, 0.3828854883388757013),
            basic::Cmplx(0.0470049721208805060, 0.9345941739461847142),
            basic::Cmplx(0.8211549515211351569, 0.1913435896808370895),
            basic::Cmplx(-0.6393122975865159763, -0.3828854883388757013)});

    poly::Polynomial p(c);
    poly::Polynomial pp(p);

    run(c, rr, cr, allr, pp, poly::General, false, 13, 7, 6);
}
