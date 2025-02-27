#include <Interpolation.hpp>
#include "include/Assertion.hpp"
#include "include/Vec.hpp"
#include "include/rel_err.hpp"

#include <algorithm>
#include <iostream>

int main() {
    using namespace std;
    using namespace intp;
    // These tests are done by comparing interpolation results with that from
    // MMA.

    Assertion assertion;
    constexpr double tol = 1e-14;

    // 1D interpolation test

    array<double, 13> f{{0.905515, 0.894638, -0.433134, 0.43131, -0.131052,
                         0.262974, 0.423888, -0.562671, -0.915567, -0.261017,
                         -0.47915, -0.00939326, -0.445962}};

    InterpolationFunction1D<0> interp1_segment(util::get_range(f));

    InterpolationFunction1D<1> interp1_linear(util::get_range(f));

    InterpolationFunction1D<3> interp1{std::make_pair(0, .5 * (f.size() - 1)),
                                       util::get_range(f)};

    auto coords_1d_half = {1.968791374707961,  0.23397112295047862,
                           4.183162505803409,  5.744438451300649,
                           0.322817339924379,  5.431730686042901,
                           1.2745684388796121, 0.21339197545741406,
                           1.5321444095845136, 0.64165470553346};

    // some random points
    auto coords_1d = {3.937582749415922,  0.46794224590095723,
                      8.366325011606818,  11.488876902601298,
                      0.645634679848758,  10.863461372085801,
                      2.5491368777592243, 0.4267839509148281,
                      3.064288819169027,  1.28330941106692};

    std::vector<double> vals_1d_segment;
    // linear interpolation results
    std::vector<double> vals_1d_linear;
    for (auto&& coord : coords_1d) {
        const std::size_t idx = static_cast<std::size_t>(std::floor(coord));
        const double val =
            f[idx] + (f[idx + 1] - f[idx]) * (coord - static_cast<double>(idx));
        vals_1d_linear.push_back(val);
        vals_1d_segment.push_back(
            f[static_cast<std::size_t>(std::round(coord))]);
    }

    // Values pre-computed by MMA
    auto vals_1d = {-0.10342680140577024, 1.4404772594187227,
                    -0.6740383362867217,  0.059044924309338276,
                    1.3546807267401677,   -0.08084215107000194,
                    0.01725408002459637,  1.4414547842595904,
                    0.444175804071652,    0.39522879177749654};

    std::cout << "\n1D Interpolation Test:\n";

    double d = rel_err(interp1_segment, util::get_range(coords_1d),
                       util::get_range(vals_1d_segment));
    assertion(d < tol);
    std::cout << "\n1D test (segment) "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    d = rel_err(interp1_linear, util::get_range(coords_1d),
                util::get_range(vals_1d_linear));
    assertion(d < tol);
    std::cout << "\n1D test (linear) "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    d = rel_err(interp1, util::get_range(coords_1d_half),
                util::get_range(vals_1d));
    assertion(d < tol);
    std::cout << "\n1D test (cubic) "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    assertion(std::abs(interp1(-.5) - (-6.3167718907512755)) < tol,
              "Out of left boundary extrapolation did not work as expected.");
    if (assertion.last_status() == 0) {
        std::cout << "Out of left boundary extrapolation work as expected.\n";
    }
    assertion(std::abs(interp1(6.5) - (-4.508470210464194)) < tol,
              "Out of right boundary extrapolation did not work as expected.");
    if (assertion.last_status() == 0) {
        std::cout << "Out of right boundary extrapolation work as expected.\n";
    }

    // 2D interpolation test

    std::cout << "\n2D Interpolation Test:\n";

    // random 5x5 mesh grid
    constexpr array<array<double, 5>, 5> f2{
        {{0.1149376591802711, 0.366986491879822, 0.8706293477610783,
          0.35545268902617266, 0.5618622985788111},
         {0.20103081669915168, -0.35136997700640915, 0.3459158000280991,
          -0.6411854570589961, -0.924275445068762},
         {0.0609591391461084, -0.9863975752412686, 0.6961713590119816,
          0.8400181479240318, 0.3758703506622121},
         {-0.849674573235994, 0.7037626673100483, 0.5167482535099701,
          -0.961602882417603, 0.9016008505150341},
         {-0.7763245335851283, 0.729592963642911, -0.5861424925445622,
          -0.3508132480272552, 0.7075736670162986}}};
    Mesh<double, 2> f2d{f2.size(), f2[0].size()};
    for (size_t i = 0; i < f2d.dim_size(0); ++i) {
        for (size_t j = 0; j < f2d.dim_size(1); ++j) { f2d(i, j) = f2[i][j]; }
    }

    InterpolationFunction<double, 2, 3> interp2{
        f2d, make_pair(0., static_cast<double>(f2d.dim_size(0)) - 1.),
        make_pair(0., static_cast<double>(f2d.dim_size(1)) - 1.)};

    // some random points
    constexpr array<array<double, 2>, 10> coords_2d{
        {{2.7956403386992656, 0.16702594760696154},
         {0.14620918612600242, 1.80574833501798},
         {0.912690726300375, 1.617752882497217},
         {3.4743760585635837, 3.235474073888084},
         {3.2819951814307826, 1.2227470066832273},
         {2.6455842429367955, 2.5079261616463135},
         {2.4398083528195187, 2.6353068861468163},
         {0.3906420637429102, 0.4383471604184015},
         {3.702137078579508, 3.366239031714958},
         {1.0831817822035914, 1.7656688800854026}}};

    // and their spline value, pre-computed by MMA
    constexpr array<double, 10> vals_2d{
        {-0.5407817111102959, 0.724190993527158, 0.25189765520792373,
         -1.4791573194664718, 1.1465605774177101, 0.3098484329528253,
         0.5464869245836507, 0.0474466547034479, -1.0511757647770876,
         0.2755974099701995}};

    assertion(interp2.uniform(0) && interp2.uniform(1),
              "Uniform properties check \033[1;91mfailed\033[0m.");

    d = rel_err(interp2, util::get_range(coords_2d), util::get_range(vals_2d));
    assertion(d < tol);
    std::cout << "\n2D test "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    try {
        interp2.at(-1, 1);
        assertion(false, "Out of boundary check \033[1;91mfailed\033[0m.\n");
    } catch (const std::domain_error&) {
        std::cout << "Out of boundary check succeed.\n";
    }

    // 3D interpolation test

    std::cout << "\n3D Interpolation Test:\n";

    // random 5x6x7 mesh grid
    constexpr array<array<array<double, 7>, 6>, 5> f3{
        {{{{0.7227454556577735, 0.12338574902709754, 0.3377337587737288,
            0.8091623938771417, 0.9293380971320708, -0.5741047532957699,
            -0.2200926742338356},
           {-0.10763006732480607, 0.0013429236738118355, 0.4513036540324049,
            -0.026003741605878705, -0.3695306079209315, 0.8136434953103922,
            -0.5116150484018025},
           {-0.8677823131366349, 0.8365100339234104, 0.7501947382164902,
            0.29910348046032365, 0.6232336678925097, 0.5087192462972618,
            -0.893065394522043},
           {-0.4851529104440644, 0.03957209074264867, 0.2539529968314116,
            -0.34512413139965936, -0.15618833950707955, -0.41941677707757874,
            0.4458983469783453},
           {-0.7660093051011847, 0.5482763224961875, -0.9987913478497763,
            -0.9842994610123474, -0.7801902911749243, 0.5613574074106569,
            -0.49371188067385186},
           {0.5583049111784795, -0.9511960118489133, -0.7610577450906089,
            0.48188343529691346, 0.23547454702291803, 0.20415133442747058,
            0.9181145422295938}}},
         {{{0.14911700752757273, -0.05837746001531219, 0.5368187139585459,
            -0.8940613010669436, -0.5018302442340068, 0.27020723614181774,
            0.12359775508403592},
           {-0.3683916194955037, -0.5525267297116154, -0.5872030595156685,
            -0.8061235674442466, -0.23709048741408623, 0.42772609271483164,
            0.4827241146251704},
           {-0.8997166066727429, 0.09739710605428176, 0.5320698398177925,
            0.05469688654269378, 0.4503706682941915, -0.7821018114726392,
            -0.3379863042117526},
           {0.5105732505948533, 0.6195932582442767, 0.6841027629496503,
            -0.024765852986073256, 0.10834864531079891, -0.49612833775897114,
            0.6109637757093997},
           {-0.28292428408539205, 0.30657791183926, -0.4741374456208911,
            0.714641208996071, 0.8309281186389432, 0.20447906224501944,
            0.1498769698488771},
           {0.38814726304764857, -0.43209235857228956, -0.8375165379497882,
            -0.9320039920129055, -0.5820061765266624, -0.6009461842921451,
            0.8425964963548309}}},
         {{{0.6848105040087451, -0.09424447387569135, 0.630497576004454,
            0.319121226100175, 0.15547601325199478, 0.37766869124805513,
            -0.418251240744814},
           {-0.5841177409538716, -0.07943995214433475, 0.6405419622330082,
            -0.18822915511374738, 0.9682166115359014, -0.24310913205955398,
            0.4207346330747752},
           {0.45689131291795304, 0.5592009847191926, -0.5794736694118892,
            0.2777598519942086, 0.06893779977151926, -0.10558235108004288,
            0.24127228408049373},
           {0.8653133697361124, 0.8808125543307121, 0.013003929742553488,
            -0.25819110436923287, 0.7150667606571139, -0.8474607538945635,
            -0.21182942625066659},
           {0.7953892915393812, -0.5146297490600227, 0.9797805357099336,
            -0.19418087953367502, 0.30615862986502584, 0.9621224147788383,
            0.591540018663443},
           {-0.0326501323196986, 0.4340496651711856, -0.6274862226468154,
            0.43246185178053453, -0.6752962418803996, -0.11639979303524317,
            0.04992688161509218}}},
         {{{0.7246870093148794, 0.41031155832285204, 0.6075964984936815,
            -0.9309929613876653, -0.5305318949868476, -0.684400029780897,
            -0.03189572751676284},
           {0.3890456863699665, 0.055814217385724785, -0.028984551458518748,
            0.1946456227803317, -0.28468283193227384, 0.07443098789465097,
            -0.3397710281207207},
           {-0.3252622666425089, -0.7764904739232343, 0.659017163729533,
            -0.6314623347234112, -0.4459102255849552, -0.8305995999617672,
            -0.6736542039955138},
           {-0.09946179372459873, 0.48571832213388744, 0.06431964245524391,
            0.9248711036318769, 0.27818775843144383, 0.06436195186180793,
            0.4804631389346339},
           {0.8394854989160261, 0.46911286378594497, -0.3880347646613327,
            -0.8793296857106343, 0.7535849141300499, -0.14621751679049622,
            -0.24084757493862208},
           {0.263291050906274, 0.2426439223615553, 0.024235715636066857,
            -0.3441033743486446, -0.6157381917061411, 0.8654487680234961,
            -0.015712235651818673}}},
         {{{0.25941803193881485, -0.5643528433065192, -0.6939218246816452,
            -0.6573164675211882, 0.3044833933508735, -0.15696192470423354,
            -0.7292088799733678},
           {0.7157059796941971, 0.5010086085806371, 0.42635964783799274,
            -0.6918122089292549, 0.5343027744642965, -0.3177068701933763,
            0.7728881262187897},
           {0.1919928613372428, 0.9481231231381191, -0.8495829859085204,
            -0.5016908143642169, -0.25281765568651915, -0.8041546515214235,
            -0.9379455299920623},
           {0.9756710429966389, 0.002338101461675013, 0.5512942665704528,
            0.11255265169286277, -0.4446511594248159, 0.923624877032104,
            -0.35888035429907994},
           {-0.4993622134330429, -0.41713411728302896, 0.5241211722633352,
            -0.8565133189111758, 0.009666595018221535, 0.0024308669289925255,
            0.21168620320785259},
           {0.7819003932609805, 0.9688407967537689, -0.438602023010886,
            0.6460148816650522, -0.463700457054959, 0.7497559992824492,
            -0.6604977204100679}}}}};

    Mesh<double, 3> f3d{f3.size(), f3[0].size(), f3[0][0].size()};
    for (size_t i = 0; i < f3d.dim_size(0); ++i) {
        for (size_t j = 0; j < f3d.dim_size(1); ++j) {
            for (size_t k = 0; k < f3d.dim_size(2); ++k) {
                f3d(i, j, k) = f3[i][j][k];
            }
        }
    }

    InterpolationFunction<double, 3, 3> interp3(
        f3d, make_pair(0., static_cast<double>(f3d.dim_size(0)) - 1.),
        make_pair(0., static_cast<double>(f3d.dim_size(1)) - 1.),
        make_pair(0., static_cast<double>(f3d.dim_size(2)) - 1.));

    // some random points
    constexpr array<array<double, 3>, 10> coords_3d{
        {{3.702442458842895, 0.3823502932775078, 2.0413783528852125},
         {0.5158748832148454, 3.998837357745451, 3.063577604910418},
         {0.9488026738796771, 1.0559217201840418, 2.5167646532589334},
         {3.787943594036949, 1.0660346792124562, 0.39476961351333895},
         {2.5437164557222864, 2.1812414553861963, 3.194304608069536},
         {0.8943435648540952, 3.274226260967546, 0.3988216219352534},
         {2.607886328617485, 1.718155000146739, 3.5736295561232927},
         {3.0725177261490373, 1.054538982002727, 0.4944889286172529},
         {1.7731799360168479, 1.237539640884, 1.0101529018305797},
         {0.8655632442457968, 1.3379555874162659, 1.7103327524709568}}};

    // and their spline value, pre-computed by MMA
    constexpr array<double, 10> vals_3d{
        {-0.20015704400375117, 0.5183267778903129, -0.7197026899005371,
         0.5183443130595233, -0.07245353025037997, 0.5534353986280456,
         -0.2674229002109916, 0.1673843822053797, -0.021928200124974297,
         -0.260677062462001}};

    d = rel_err(interp3, util::get_range(coords_3d), util::get_range(vals_3d));
    assertion(d < tol);
    std::cout << "\n3D test "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    assertion(!interp3.periodicity(0) && !interp3.periodicity(1) &&
              !interp3.periodicity(2));

    // 1D interpolation test with periodic boundary

    std::cout << "\n1D Interpolation with Periodic Boundary Test:\n";

#ifdef INTP_PERIODIC_NO_DUMMY_POINT
    const auto f_periodic_range = std::make_pair(f.begin(), std::prev(f.end()));
    InterpolationFunction1D<0> interp1_periodic_segment(f_periodic_range, true);

    InterpolationFunction1D<1> interp1_periodic_linear(f_periodic_range, true);

    InterpolationFunction1D<4> interp1_periodic(f_periodic_range, true);

#else
    // Since the InterpolationFunction ignores the last term along periodic
    // dimension, thus we can use the origin non-periodic data to interpolate a
    // periodic spline function

    InterpolationFunction1D<1> interp1_periodic_segment(util::get_range(f),
                                                        true);
    InterpolationFunction1D<1> interp1_periodic_linear(util::get_range(f),
                                                       true);

    InterpolationFunction1D<4> interp1_periodic(util::get_range(f), true);
#endif

    // linear interpolation results
    std::vector<double> vals_1d_periodic_segment;
    std::vector<double> vals_1d_periodic_linear;
    for (auto&& coord : coords_1d) {
        const std::size_t idx = static_cast<std::size_t>(std::floor(coord));
        const double val = f[idx] + (f[(idx + 1) % (f.size() - 1)] - f[idx]) *
                                        (coord - static_cast<double>(idx));
        vals_1d_periodic_linear.push_back(val);
        vals_1d_periodic_segment.push_back(
            f[static_cast<std::size_t>(std::round(coord))]);
    }

    auto vals_1d_periodic = {
        -0.09762254647017743, 1.168800757853312,  -0.6682906902062101,
        0.44128062499073417,  1.1642814012848224, -0.12481864202139212,
        -0.00445415461594469, 1.159234754035029,  0.4508845136779133,
        0.45967108162968584,
    };

    d = rel_err(interp1_periodic_segment, util::get_range(coords_1d),
                util::get_range(vals_1d_periodic_segment));
    assertion(d < tol);
    std::cout << "\n1D test (segment) with periodic boundary "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    d = rel_err(interp1_periodic_linear, util::get_range(coords_1d),
                util::get_range(vals_1d_periodic_linear));
    assertion(d < tol);
    std::cout << "\n1D test (linear) with periodic boundary "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    d = rel_err(interp1_periodic, util::get_range(coords_1d),
                util::get_range(vals_1d_periodic));
    assertion(d < tol);
    std::cout << "\n1D test (quartic) with periodic boundary "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    assertion(std::abs(interp1_periodic(*coords_1d.begin() - (f.size() - 1)) -
                       *vals_1d_periodic.begin()) < tol,
              "Out of left periodic boundary did not work as expected.\n");
    assertion(std::abs(interp1_periodic(*coords_1d.begin() + (f.size() - 1)) -
                       *vals_1d_periodic.begin()) < tol,
              "Out of Right periodic boundary did not work as expected.\n");

    // 2D interpolation test with one dimension being periodic boundary
    std::cout << "\n2D Interpolation with one Periodic Boundary Test:\n";

#ifdef INTP_PERIODIC_NO_DUMMY_POINT
    Mesh<double, 2> f2d_periodic_y{f2.size(), f2[0].size() - 1};
    for (size_t i = 0; i < f2d_periodic_y.dim_size(0); ++i) {
        for (size_t j = 0; j < f2d_periodic_y.dim_size(1); ++j) {
            f2d_periodic_y(i, j) = f2[i][j];
        }
    }
    InterpolationFunction<double, 2, 3> interp2_periodic(
        {false, true}, f2d_periodic_y,
        make_pair(0., static_cast<double>(f2d.dim_size(0)) - 1.),
        make_pair(0., static_cast<double>(f2d.dim_size(1)) - 1.));

#else
    InterpolationFunction<double, 2, 3> interp2_periodic(
        {false, true}, f2d,
        make_pair(0., static_cast<double>(f2d.dim_size(0)) - 1.),
        make_pair(0., static_cast<double>(f2d.dim_size(1)) - 1.));
#endif

    auto vals_2d_periodic = {-0.5456439415470818, 0.7261218483070795,
                             0.21577722210958022, -1.6499933881987376,
                             1.224021619908732,   0.34969937176489274,
                             0.5845798304532657,  0.2875923130734858,
                             -1.4740569960870218, 0.258215214830246};

    d = rel_err(interp2_periodic, util::get_range(coords_2d),
                util::get_range(vals_2d_periodic));
    assertion(d < tol);
    std::cout << "\n2D test with periodic boundary "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

#ifdef _DEBUG
    if (assertion.last_status() != 0) {
        interp2_periodic.spline().debug_output();
    }
#endif

    // 1D interpolation derivative test

    std::cout << "\n1D Interpolation Derivative Test:\n";

    auto vals_1d_derivative_1 = {-1.027459647923634,  -0.23460057495062184,
                                 1.626359518442344,   -0.5964266359199579,
                                 -1.6344806285675368, 1.1365216604820132,
                                 2.451880834458231,   0.14293981479176837,
                                 0.16561298233542285, -3.58205140200432};

    d = rel_err(
        [&](double x) { return interp1.derivative_at(std::make_pair(x, 1)); },
        util::get_range(coords_1d_half), util::get_range(vals_1d_derivative_1));
    assertion(d < tol);
    std::cout << "\n1D test of derivative "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    auto vals_1d_derivative_periodic = {
        -0.5956189309592923, 0.18879244402223305, 0.823542093101076,
        0.9423837958144585,  -0.2574829498865512, 0.8141621781881445,
        1.297682043678497,   0.2748905754634491,  0.20278909403301104,
        -1.744727637863168};

    d = rel_err(
        [&](double x) {
            return interp1_periodic.derivative_at(std::make_pair(x, 1));
        },
        util::get_range(coords_1d),
        util::get_range(vals_1d_derivative_periodic));
    assertion(d < tol);
    std::cout << "\n1D test of derivative in periodic case "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    // 2D interpolation derivative test

    std::cout << "\n2D Interpolation Derivative Test:\n";

    auto vals_2d_derivative_x2_y1 = {-0.7543201698562789, 7.386529487473519,
                                     1.4374808443975349,  0.015013232428628703,
                                     0.4101648846318855,  0.3989275680682247,
                                     -1.2969356922668822, -3.0725969893922946,
                                     -3.4258871562563455, 0.43843555890798447};

    d = rel_err(
        [&](std::array<double, 2> coord) {
            return interp2.derivative_at(coord, 2, 1);
        },
        util::get_range(coords_2d), util::get_range(vals_2d_derivative_x2_y1));
    assertion(d < tol);
    std::cout << "\n2D test of derivative (2,1) "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    // 3D interpolation derivative test

    std::cout << "\n3D Interpolation Derivative Test:\n";

    auto vals_3d_derivative_x1_y0_z3 = {
        26.856939941150102, -2.3883090024921603, 8.857003271951815,
        -5.988828759251878, -9.656188308304278,  -3.8790002779026658,
        -8.950315245767326, 2.859291392187189,   -1.0605065454924238,
        0.9708580137688204};

    d = rel_err(
        [&](std::array<double, 3> coord) {
            return interp3.derivative_at(coord, {1, 0, 3});
        },
        util::get_range(coords_3d),
        util::get_range(vals_3d_derivative_x1_y0_z3));
    assertion(d < tol);
    std::cout << "\n3D test of derivative (1,0,3) "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    // 1D non-uniform interpolation test

    std::cout << "\n1D nonuniform Interpolation Test:\n";

    auto input_coords_1d = {0.,
                            0.9248468590818175,
                            2.270222487146226,
                            3.498804659217474,
                            3.7530418178453955,
                            4.791302870662092,
                            6.3336649913937375,
                            7.041926963118405,
                            7.612472537805882,
                            9.034787492568283,
                            10.078257485188475,
                            11.440973163521294,
                            12.};

    // InterpolationFunction1D<0> interp1_nonuniform_segment(
    //     util::get_range(input_coords_1d), util::get_range(f));

    InterpolationFunction1D<1> interp1_nonuniform_linear(
        util::get_range(input_coords_1d), util::get_range(f));

    InterpolationFunction1D<3> interp1_nonuniform(
        util::get_range(input_coords_1d), util::get_range(f));

    // linear interpolation results
    std::vector<double> vals_1d_nonuniform_linear;
    for (auto&& coord : coords_1d) {
        const auto iter = std::lower_bound(input_coords_1d.begin(),
                                           input_coords_1d.end(), coord);
        const std::size_t idx =
            static_cast<std::size_t>(iter - input_coords_1d.begin());
        const double val = f[idx - 1] + (f[idx] - f[idx - 1]) *
                                            (coord - *std::prev(iter)) /
                                            (*iter - *std::prev(iter));
        vals_1d_nonuniform_linear.push_back(val);
    }

    auto vals_1d_nonuniform = {-0.4562057772431492, 1.3471094229928755,
                               -0.6079379355534298, -0.016699918339397407,
                               1.2453576785633913,  -0.17792750766792387,
                               0.05227508784539216, 1.3530673230583836,
                               0.877847815343884,   0.27427682690484234};

    d = rel_err(interp1_nonuniform_linear, util::get_range(coords_1d),
                util::get_range(vals_1d_nonuniform_linear));
    assertion(d < tol);
    std::cout << "\n1D nonuniform test (linear) "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    d = rel_err(interp1_nonuniform, util::get_range(coords_1d),
                util::get_range(vals_1d_nonuniform));
    assertion(d < tol);
    std::cout << "\n1D nonuniform test (cubic) "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    std::cout << "\n1D nonuniform Interpolation Test with Periodic Boundary:\n";

    // 1D non-uniform periodic interpolation test

#ifdef INTP_PERIODIC_NO_DUMMY_POINT
    InterpolationFunction1D<1> interp1_nonuniform_periodic_linear(
        util::get_range(input_coords_1d), f_periodic_range, true);

    InterpolationFunction1D<4> interp1_nonuniform_periodic(
        util::get_range(input_coords_1d), f_periodic_range, true);
#else
    InterpolationFunction1D<1> interp1_nonuniform_periodic_linear(
        util::get_range(input_coords_1d), util::get_range(f), true);

    InterpolationFunction1D<4> interp1_nonuniform_periodic(
        util::get_range(input_coords_1d), util::get_range(f), true);
#endif

    // linear interpolation results
    std::vector<double> vals_1d_nonuniform_periodic_linear;
    for (auto&& coord : coords_1d) {
        const auto iter = std::lower_bound(input_coords_1d.begin(),
                                           input_coords_1d.end(), coord);
        const std::size_t idx =
            static_cast<std::size_t>(iter - input_coords_1d.begin());
        const double val = f[idx - 1] + (f[idx % (f.size() - 1)] - f[idx - 1]) *
                                            (coord - *std::prev(iter)) /
                                            (*iter - *std::prev(iter));
        vals_1d_nonuniform_periodic_linear.push_back(val);
    }

    auto vals_1d_nonuniform_periodic = {
        -0.4464803199487985,    1.2298382629311808, -0.6268441539418139,
        0.06800541396058792,    1.1730724308997087, -0.6075709757451075,
        0.00021566856863144274, 1.2284708203600656, 0.8052326379455318,
        0.3011733005197456};

    d = rel_err(interp1_nonuniform_periodic_linear, util::get_range(coords_1d),
                util::get_range(vals_1d_nonuniform_periodic_linear));
    assertion(d < tol);
    std::cout << "\n1D nonuniform test (linear) with periodic boundary "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    d = rel_err(interp1_nonuniform_periodic, util::get_range(coords_1d),
                util::get_range(vals_1d_nonuniform_periodic));
    assertion(d < tol);
    std::cout << "\n1D nonuniform test (quartic) with periodic boundary "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    // 2D non-uniform with one dimension being periodic boundary

    std::cout << "\n2D nonuniform Interpolation Test with Periodic Boundary:\n";

    auto nonuniform_coord_for_2d = {0., 1.329905262345947, 2.200286645200226,
                                    3.1202827237815516, 4.};

#ifdef INTP_PERIODIC_NO_DUMMY_POINT
    Mesh<double, 2> f2d_periodic_x{f2.size() - 1, f2[0].size()};
    for (size_t i = 0; i < f2d_periodic_x.dim_size(0); ++i) {
        for (size_t j = 0; j < f2d_periodic_x.dim_size(1); ++j) {
            f2d_periodic_x(i, j) = f2[i][j];
        }
    }
    InterpolationFunction<double, 2, 3> interp2_X_periodic_Y_nonuniform(
        {true, false}, f2d_periodic_x,
        std::make_pair(0., static_cast<double>(f2d.dim_size(0)) - 1.),
        util::get_range(nonuniform_coord_for_2d));
#else
    InterpolationFunction<double, 2, 3> interp2_X_periodic_Y_nonuniform(
        {true, false}, f2d,
        std::make_pair(0., static_cast<double>(f2d.dim_size(0)) - 1.),
        util::get_range(nonuniform_coord_for_2d));
#endif

    auto vals_2d_X_periodic_Y_nonuniform = {
        -0.7160424002258807,  0.6702846215275077,  0.04933393138376427,
        -0.48309957376784946, 0.7622619905083738,  0.28742506313682975,
        0.4762197211423307,   -0.2924876444803407, -0.05355290342562366,
        0.08938131323971249};

    d = rel_err(interp2_X_periodic_Y_nonuniform, util::get_range(coords_2d),
                util::get_range(vals_2d_X_periodic_Y_nonuniform));
    assertion(d < tol);
    std::cout << "\n2D test with x-periodic and y-nonuniform "
              << (assertion.last_status() == 0 ? "succeed"
                                               : "\033[1;91mfailed\033[0m")
              << '\n';
    std::cout << "Relative Error = " << d << '\n';

    // interpolation on vector

    std::vector<Vec<2, float>> pts;
    constexpr std::size_t pts_n = 31;
#ifdef INTP_PERIODIC_NO_DUMMY_POINT
    for (std::size_t i = 0; i < pts_n; ++i) {
#else
    for (std::size_t i = 0; i <= pts_n; ++i) {
#endif
        double theta = 2. * M_PI * static_cast<double>(i) / pts_n;
        pts.emplace_back(std::cos(theta), std::sin(theta));
    }
    InterpolationFunction1D<3, Vec<2, float>, float> circle(
        util::get_range(pts), true);
    {
        constexpr std::size_t sample_n = 1024;
        float err = 0;
        for (std::size_t i = 0; i < sample_n; ++i) {
#ifndef M_PIf
#define M_PIf 3.14159265358979323846f
#endif
            float theta = 2.f * M_PIf * static_cast<float>(i) / sample_n;
            err += std::abs(circle(theta).mag() - 1.f);
        }
        err /= sample_n;
        assertion(err < 1.e-5f);
        std::cout << "\nInterpolation of 2d points on a circle "
                  << (assertion.last_status() == 0 ? "succeed"
                                                   : "\033[1;91mfailed\033[0m")
                  << '\n';
        std::cout << "Relative Error = " << err << '\n';
    }

    return assertion.status();
}
