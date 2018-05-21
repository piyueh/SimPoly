/**
 * \file find_roots.cpp
 * \brief Unit tests for roots finding functions.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-02-01
 */


# include <random>
# include <algorithm>
# include <cmath>

# include <gtest/gtest.h>

# include "basic.h"
# include "exceptions.h"

using namespace simpoly;

extern std::default_random_engine generator;
static std::uniform_real_distribution<double> drand(-1.0, 1.0);

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

static auto f2 = [](basic::Cmplx i, basic::Cmplx j)->bool{return i.real()<j.real();};

# ifndef NDEBUG
TEST(PolynomialRoots, CheckExceptions)
{
    ASSERT_THROW(basic::aberth(basic::CArry(0), basic::CArry(0)),
            simpoly::exceptions::ZeroCoeffsLength);

    ASSERT_THROW(basic::aberth(basic::DArry(0), basic::CArry(0)),
            simpoly::exceptions::ZeroCoeffsLength);

    ASSERT_THROW(basic::aberth(basic::CArry(0), basic::DArry(0)),
            simpoly::exceptions::ZeroCoeffsLength);

    ASSERT_THROW(basic::aberth(basic::DArry(0), basic::DArry(0)),
            simpoly::exceptions::ZeroCoeffsLength);

    ASSERT_THROW(basic::aberth(basic::CArry(0)),
            simpoly::exceptions::ZeroCoeffsLength);

    ASSERT_THROW(basic::aberth(basic::DArry(0)),
            simpoly::exceptions::ZeroCoeffsLength);

    ASSERT_THROW(basic::yan_and_chieng_2006(basic::CArry(0)),
            simpoly::exceptions::ZeroCoeffsLength);

    ASSERT_THROW(basic::yan_and_chieng_2006(basic::DArry(0)),
            simpoly::exceptions::ZeroCoeffsLength);
}
# endif

TEST(PolynomialRoots, ComplexRootsAberth1)
{

    basic::CArry c({basic::Cmplx(-0.2333161111467611, -0.07203119056462484),
            basic::Cmplx(0.4691992245557093, 0.1069593270345827),
            basic::Cmplx(-0.4626667261153909, 0.01132740583617468),
            basic::Cmplx(-0.2089349052104357, 1.204569153622907),
            basic::Cmplx(1.937034853087647, -0.8696243603297597),
            basic::Cmplx(0.1479889416963552, -1.895113546067258),
            basic::Cmplx(-2.340724865590584, -0.01971983376224859),
            basic::Cmplx(-0.5151163499273337, 1.464560317090324),
            basic::Cmplx(1.0, 0.0)});

    basic::CArry expect({basic::Cmplx(0.7005658123938592, -0.8463540437130455),
            basic::Cmplx(0.9446989441874694, -0.257812273394918),
            basic::Cmplx(0.05705818190287304, 0.4239984931397236),
            basic::Cmplx(-0.8895720736786421, -0.4724245379748282),
            basic::Cmplx(-0.9000656197820271, 0.5501970894798538),
            basic::Cmplx(-0.8014849228035563, -0.5048110577579954),
            basic::Cmplx(0.6918133537758795, -0.1193762463708401),
            basic::Cmplx(0.7121026739314782, -0.2379777404982741)});

    basic::CArry guess(expect.size());
    for(auto &it: guess) it = basic::Cmplx(drand(generator), drand(generator));

    basic::CArry result = basic::aberth(c, guess);

    // sort `expect` and `result` so we can compare them term by term
    std::sort(std::begin(expect), std::end(expect), f);
    std::sort(std::begin(result), std::end(result), f);

    for(unsigned i=0; i<expect.size(); ++i)
    {
        ASSERT_NEAR(expect[i].real(), result[i].real(), 1e-10);
        ASSERT_NEAR(expect[i].imag(), result[i].imag(), 1e-10);
    }
}

TEST(PolynomialRoots, ComplexRootsAberth2)
{

    basic::CArry c({
            basic::Cmplx(-0.005728250342466052, -0.002271398296873492),
            basic::Cmplx(-0.01308912745807032, -0.01213438492830167),
            basic::Cmplx(-0.02955995754847203, 0.03399218074497463),
            basic::Cmplx(-0.2140782539763119, 0.1236268405829473),
            basic::Cmplx(-0.3066913953594199, 0.01961305798663569),
            basic::Cmplx(0.100694771640935, 0.5612758300940467),
            basic::Cmplx(-0.6193124230980729, 1.405897414495),
            basic::Cmplx(-0.2255050490430183, 0.4788858014628719),
            basic::Cmplx(1.784790975021791, 1.09577302430897),
            basic::Cmplx(2.041509181070611, 1.483192570436406),
            basic::Cmplx(3.610434184456693, -0.326446012171619),
            basic::Cmplx(5.304415179297847, -0.370402492941317),
            basic::Cmplx(5.745823527966784, -2.0847187371134),
            basic::Cmplx(7.311086993663613, -1.626690508281606),
            basic::Cmplx(4.782488494653972, -2.104480116649998),
            basic::Cmplx(5.855325915119935, -3.399342743824328),
            basic::Cmplx(4.255103489250693, -2.165302523845158),
            basic::Cmplx(3.91570950673342, -1.740942782152447),
            basic::Cmplx(0.9524783367812539, -0.6506310156921662),
            basic::Cmplx(1.141297932726919, -1.94837189217718),
            basic::Cmplx(1.0, 0.0) });

    basic::CArry expect({
            basic::Cmplx(0.7975362840622007, -0.8462970398548786),
            basic::Cmplx(-0.3807754021788712, -0.06812598591443542),
            basic::Cmplx(-0.1304125017095847, 0.7766448933011634),
            basic::Cmplx(-0.001963077096688526, 0.8622475739085034),
            basic::Cmplx(0.7175359397080345, 0.8037881512142599),
            basic::Cmplx(0.2107341754008711, -0.2366132065261823),
            basic::Cmplx(0.6555109430788004, -0.1170830430835776),
            basic::Cmplx(-0.3270742988659125, -0.9078065531974187),
            basic::Cmplx(-0.2645322491092501, 0.9236502191723146),
            basic::Cmplx(-0.7199783731914202, -0.1537600692786125),
            basic::Cmplx(-0.1084740382683762, 0.6145022185388467),
            basic::Cmplx(-0.354895207020788, 0.1991856712898046),
            basic::Cmplx(-0.7567110900552072, -0.7691744505329461),
            basic::Cmplx(0.317398215607551, -0.4452176445162543),
            basic::Cmplx(0.2834616418719842, -0.8552686943685404),
            basic::Cmplx(-0.08733263110724621, 0.9277137258818431),
            basic::Cmplx(0.6090163019740282, -0.8409782066445819),
            basic::Cmplx(0.1121993093826672, 0.941367896408785),
            basic::Cmplx(-0.7887725695474679, 0.3323586067639508),
            basic::Cmplx(-0.9237693056622434, 0.8072378296151359) });

    basic::DArry guess(expect.size());
    for(auto &it: guess) it = drand(generator);

    basic::CArry result = basic::aberth(c, guess);

    // sort `expect` and `result` so we can compare them term by term
    std::sort(std::begin(expect), std::end(expect), f);
    std::sort(std::begin(result), std::end(result), f);

    for(unsigned i=0; i<expect.size(); ++i)
    {
        ASSERT_NEAR(expect[i].real(), result[i].real(), 1e-10);
        ASSERT_NEAR(expect[i].imag(), result[i].imag(), 1e-10);
    }
}

TEST(PolynomialRoots, ComplexRootsAberth3)
{

    basic::DArry c({
            -0.1195107918695455,  0.8447947406345584, -0.6030868103450979,
            -0.9766533913297424, -0.2677767116634977,  0.3804660521221173,
            -0.4716581713021364, -0.8085712194756074,  0.1556975750852925,
            -0.6073461984864883,  0.8943598261002668, -0.0478803277253317,
            -0.5541643941753234, -0.008585713188056 ,  0.9873624852167924,
            0.1332436172903213 });

    basic::CArry expect({
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

    basic::CArry guess(expect.size());
    for(auto &it: guess) it = basic::Cmplx(drand(generator), drand(generator));

    basic::CArry result = basic::aberth(c, guess);

    // sort `expect` and `result` so we can compare them term by term
    std::sort(std::begin(expect), std::end(expect), f);
    std::sort(std::begin(result), std::end(result), f);

    for(unsigned i=0; i<expect.size(); ++i)
    {
        ASSERT_NEAR(expect[i].real(), result[i].real(), 1e-10);
        ASSERT_NEAR(expect[i].imag(), result[i].imag(), 1e-10);
    }
}

TEST(PolynomialRoots, ComplexRootsAberth4)
{

    basic::DArry c({
            0.5808674560998166, -0.6495173545341641,  0.2415758641182864,
            -0.5399634531874522, -0.1987466764202912,  0.524098092459949 ,
            0.9492147992555258, -0.8465445984978208, -0.2634181317553206 });

    basic::CArry expect({
            basic::Cmplx(-3.9897450303048148, 0.),
            basic::Cmplx(-0.8870172577239701, -0.566285872801541),
            basic::Cmplx(-0.8870172577239701, 0.566285872801541),
            basic::Cmplx(-0.0269354166777918, -0.7935763762618989),
            basic::Cmplx(-0.0269354166777918, 0.7935763762618989),
            basic::Cmplx(0.8661900415687693, 0.),
            basic::Cmplx(0.8688844424396298, -0.3985703106289106),
            basic::Cmplx(0.8688844424396298, 0.3985703106289106)});

    basic::DArry guess(expect.size());
    for(auto &it: guess) it = drand(generator);

    basic::CArry result = basic::aberth(c, guess);

    // sort `expect` and `result` so we can compare them term by term
    std::sort(std::begin(expect), std::end(expect), f);
    std::sort(std::begin(result), std::end(result), f);

    for(unsigned i=0; i<expect.size(); ++i)
    {
        ASSERT_NEAR(expect[i].real(), result[i].real(), 1e-10);
        ASSERT_NEAR(expect[i].imag(), result[i].imag(), 1e-10);
    }
}

TEST(PolynomialRoots, ComplexRootsAberth5)
{

    basic::CArry c({
            basic::Cmplx(-0.006895471860191621, -0.001561652498943983),
            basic::Cmplx(-0.004665280949253687, 0.01159701097464432),
            basic::Cmplx(0.04893628621350306, 0.002564810956982647),
            basic::Cmplx(0.1021068565648835, -0.0778142575160364),
            basic::Cmplx(-0.07586962266273767, -0.3021235965460747),
            basic::Cmplx(-0.3012987800316255, -0.3799700159537112),
            basic::Cmplx(-0.5084001253398745, 0.1911652072168439),
            basic::Cmplx(-0.7708335428363232, 0.6687961657562093),
            basic::Cmplx(0.505299934562367, 0.1622859246626371),
            basic::Cmplx(1.0, 0.0) });


    basic::CArry expect({
            basic::Cmplx(0.1771003961465119, -0.4981138554161064),
            basic::Cmplx(0.4633804950170575, -0.4035761261924611),
            basic::Cmplx(-0.1712007201504315, 0.6307368750421651),
            basic::Cmplx(-0.5929071284868146, -0.1075285838151725),
            basic::Cmplx(-0.9567684509846892, 0.1545082143067493),
            basic::Cmplx(-0.4307735278681735, 0.3268702315765142),
            basic::Cmplx(0.2943806862289293, 0.06505809552997199),
            basic::Cmplx(0.9608825921769193, -0.06627893962635278),
            basic::Cmplx(-0.2493942766416761, -0.263961836067945) });


    basic::CArry result = basic::aberth(c);

    // sort `expect` and `result` so we can compare them term by term
    std::sort(std::begin(expect), std::end(expect), f);
    std::sort(std::begin(result), std::end(result), f);

    for(unsigned i=0; i<expect.size(); ++i)
    {
        ASSERT_NEAR(expect[i].real(), result[i].real(), 1e-10);
        ASSERT_NEAR(expect[i].imag(), result[i].imag(), 1e-10);
    }
}

TEST(PolynomialRoots, ComplexRootsAberth6)
{

    basic::DArry c({
            -0.0728414666575747, -0.1931113936126789, -0.0998162346855869,
            -0.6720636118447425, -0.6376528554841387,  0.897012721205759 ,
            0.3974499550174013, -0.174367579389703 , -0.5953274416841128,
            0.4356308848726764, -0.0437916918310677,  0.923789663992717 ,
            -0.5211615410858714,  0.1200906586123869,  0.3267892130542112,
            0.4086253995278857,  0.508166061783101 });


    basic::CArry expect({
            basic::Cmplx(-1.4408869505091435, 0.),
            basic::Cmplx(-0.7299902298252351, -0.5361193861781802),
            basic::Cmplx(-0.7299902298252351, 0.5361193861781802),
            basic::Cmplx(-0.6102568526266997, 0.),
            basic::Cmplx(-0.4774192936120940, -1.1167377487906553),
            basic::Cmplx(-0.4774192936120940, 1.1167377487906553),
            basic::Cmplx(-0.3613144393856812, 0.),
            basic::Cmplx(-0.1326344194445007, -1.168802820704177),
            basic::Cmplx(-0.1326344194445007, 1.168802820704177),
            basic::Cmplx(0.1384183504785491, -0.4838537114619603),
            basic::Cmplx(0.1384183504785491, 0.4838537114619603),
            basic::Cmplx(0.6470536898794871, -0.8326642583192213),
            basic::Cmplx(0.6470536898794871, 0.8326642583192213),
            basic::Cmplx(0.8910630266770523, -0.478482517291651),
            basic::Cmplx(0.8910630266770523, 0.478482517291651),
            basic::Cmplx(0.9353581471985004, 0.) });


    basic::CArry result = basic::aberth(c);

    // sort `expect` and `result` so we can compare them term by term
    std::sort(std::begin(expect), std::end(expect), f);
    std::sort(std::begin(result), std::end(result), f);

    for(unsigned i=0; i<expect.size(); ++i)
    {
        ASSERT_NEAR(expect[i].real(), result[i].real(), 1e-10);
        ASSERT_NEAR(expect[i].imag(), result[i].imag(), 1e-10);
    }
}

TEST(PolynomialRoots, RealRootsAberth1)
{
    basic::DArry c({
            0.0005010250303020029,  0.006234914567457598,
            -0.008217029116396798, -0.1730879244322208,  0.2766216841964072,
            1.248746139278724, -3.203690116155496, -0.1925619244634082,
            6.985683557019738, -7.280199120562615,  2.34});

    basic::DArry expect({
            0.7602860830452414, -0.1591090490225655,  0.9977072170854013,
            0.357143338844115,  0.7414751420913581, -0.0927725821643739,
            -0.5554400402569037,  0.8916535577665332,  0.4763459267936514,
            -0.3060933888138182});

    basic::DArry guess(expect.size());
    for(auto &it: guess) it = drand(generator);

    basic::CArry result = basic::aberth(c, guess);

    // sort `expect` and `result` so we can compare them term by term
    std::sort(std::begin(expect), std::end(expect));
    std::sort(std::begin(result), std::end(result), f2);

    for(unsigned i=0; i<expect.size(); ++i)
    {
        ASSERT_NEAR(expect[i], result[i].real(), 1e-10);
        ASSERT_NEAR(0.0, result[i].imag(), 1e-10);
    }
}

TEST(PolynomialRoots, RealRootsAberth2)
{
    basic::DArry c({
            -1.196351529195439e-07,  4.919973730289604e-06,
            -1.986062723842243e-05, -0.0002469080505527796,
            0.001284978870504016,  0.00348963837303726, -0.02404374066899382,
            -0.01347466727223127,  0.200125699333026, -0.06577075073777044,
            -0.810921913206468,  0.6347898119672782,  1.538343111992133,
            -1.605177209592742, -1.062861119607442,  1.278 });

    basic::DArry expect({
            -0.5036077278270599,  0.5253156753043096, -0.1863051985638844,
            0.4066469386958105, -0.327348593000119,  0.472961157808012,
            -0.5504948721890652,  0.1513160141154197,  0.2625784242407205,
            0.5203189352440989,  0.4459543391846599, -0.3634891145277974,
            -0.7818520862771727,  0.7310561595323164,  0.02860966626244443 });

    basic::DArry guess(expect.size());
    for(auto &it: guess) it = drand(generator);

    basic::CArry result = basic::aberth(c, guess, 1e-10);

    // sort `expect` and `result` so we can compare them term by term
    std::sort(std::begin(expect), std::end(expect));
    std::sort(std::begin(result), std::end(result), f2);

    for(unsigned i=0; i<expect.size(); ++i)
    {
        ASSERT_NEAR(expect[i], result[i].real(), 1e-10);
        ASSERT_NEAR(0.0, result[i].imag(), 1e-10);
    }
}

TEST(PolynomialRoots, ComplexRootsYC1)
{

    basic::CArry c({basic::Cmplx(-0.2333161111467611, -0.07203119056462484),
            basic::Cmplx(0.4691992245557093, 0.1069593270345827),
            basic::Cmplx(-0.4626667261153909, 0.01132740583617468),
            basic::Cmplx(-0.2089349052104357, 1.204569153622907),
            basic::Cmplx(1.937034853087647, -0.8696243603297597),
            basic::Cmplx(0.1479889416963552, -1.895113546067258),
            basic::Cmplx(-2.340724865590584, -0.01971983376224859),
            basic::Cmplx(-0.5151163499273337, 1.464560317090324),
            basic::Cmplx(1.0, 0.0)});

    basic::CArry expect({basic::Cmplx(0.7005658123938592, -0.8463540437130455),
            basic::Cmplx(0.9446989441874694, -0.257812273394918),
            basic::Cmplx(0.05705818190287304, 0.4239984931397236),
            basic::Cmplx(-0.8895720736786421, -0.4724245379748282),
            basic::Cmplx(-0.9000656197820271, 0.5501970894798538),
            basic::Cmplx(-0.8014849228035563, -0.5048110577579954),
            basic::Cmplx(0.6918133537758795, -0.1193762463708401),
            basic::Cmplx(0.7121026739314782, -0.2379777404982741)});

    basic::CArry result = basic::yan_and_chieng_2006(c);

    // sort `expect` and `result` so we can compare them term by term
    std::sort(std::begin(expect), std::end(expect), f);
    std::sort(std::begin(result), std::end(result), f);

    for(unsigned i=0; i<expect.size(); ++i)
    {
        ASSERT_NEAR(expect[i].real(), result[i].real(), 1e-10);
        ASSERT_NEAR(expect[i].imag(), result[i].imag(), 1e-10);
    }
}

TEST(PolynomialRoots, ComplexRootsYC2)
{

    basic::DArry c({
            -0.1195107918695455,  0.8447947406345584, -0.6030868103450979,
            -0.9766533913297424, -0.2677767116634977,  0.3804660521221173,
            -0.4716581713021364, -0.8085712194756074,  0.1556975750852925,
            -0.6073461984864883,  0.8943598261002668, -0.0478803277253317,
            -0.5541643941753234, -0.008585713188056 ,  0.9873624852167924,
            0.1332436172903213 });

    basic::CArry expect({
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

    basic::CArry result = basic::yan_and_chieng_2006(c);

    // sort `expect` and `result` so we can compare them term by term
    std::sort(std::begin(expect), std::end(expect), f);
    std::sort(std::begin(result), std::end(result), f);

    for(unsigned i=0; i<expect.size(); ++i)
    {
        ASSERT_NEAR(expect[i].real(), result[i].real(), 1e-10);
        ASSERT_NEAR(expect[i].imag(), result[i].imag(), 1e-10);
    }
}

TEST(PolynomialRoots, RealRootsYC)
{
    basic::DArry c({
            0.0005010250303020029,  0.006234914567457598,
            -0.008217029116396798, -0.1730879244322208,  0.2766216841964072,
            1.248746139278724, -3.203690116155496, -0.1925619244634082,
            6.985683557019738, -7.280199120562615,  2.34});

    basic::DArry expect({
            0.7602860830452414, -0.1591090490225655,  0.9977072170854013,
            0.357143338844115,  0.7414751420913581, -0.0927725821643739,
            -0.5554400402569037,  0.8916535577665332,  0.4763459267936514,
            -0.3060933888138182});

    basic::CArry result = basic::yan_and_chieng_2006(c);

    // sort `expect` and `result` so we can compare them term by term
    std::sort(std::begin(expect), std::end(expect));
    std::sort(std::begin(result), std::end(result), f2);

    for(unsigned i=0; i<expect.size(); ++i)
    {
        ASSERT_NEAR(expect[i], result[i].real(), 1e-10);
        ASSERT_NEAR(0.0, result[i].imag(), 1e-10);
    }
}

TEST(PolynomialRoots, MultiRealRootsYC1)
{
    basic::DArry c({
            -7.5202173502138876682e-04,  -8.6835988738648978158e-03,
            5.5348919281959091387e-02,   1.4315233176654140745e-01,
            -4.3769681353357492437e-01,  -1.2069048603868517411e+00,
            -2.8638873139566145554e-02,   1.6438835517368903805e+00,
            1.0000000000000000000e+00});


    basic::DArry expect({
            -0.7592589202499554490, -0.7515504281961997712,
            -0.7121750519637163324, -0.6197568280813230412,
            -0.0649616158126642695, 0.2431327906295455976,
            0.2431327906295455976, 0.7775537113078748730});


    basic::CArry result = basic::yan_and_chieng_2006(c);

    // sort `expect` and `result` so we can compare them term by term
    std::sort(std::begin(expect), std::end(expect));
    std::sort(std::begin(result), std::end(result), f2);

    for(unsigned i=0; i<expect.size(); ++i)
    {
        ASSERT_NEAR(expect[i], result[i].real(), 1e-8);
        ASSERT_NEAR(0.0, result[i].imag(), 1e-10);
    }
}

TEST(PolynomialRoots, MultiRealRootsYC2)
{
    basic::DArry c({
            0.0, -2.94912e-02, 8.84736e-02, 9.611264e-01, -3.1193088e+00,
            -5.5872e+00, 2.51584e+01, -2.208e+00, -6.3264e+01, 5.2e+01, 3.6e+01,
            -6.0e+01, 2.0e+01});


    basic::DArry expect({
            1., 1., 1., 0.8, 0.6, 0.4, 0.2, 0., -0.2, -0.4, -0.6, -0.8 });


    basic::CArry result = basic::yan_and_chieng_2006(c);

    // sort `expect` and `result` so we can compare them term by term
    std::sort(std::begin(expect), std::end(expect));
    std::sort(std::begin(result), std::end(result), f2);

    for(unsigned i=0; i<expect.size(); ++i)
    {
        ASSERT_NEAR(expect[i], result[i].real(), 1e-8);
        ASSERT_NEAR(0.0, result[i].imag(), 1e-10);
    }
}

// TODO: add tests for roots of low-degree polynomials
