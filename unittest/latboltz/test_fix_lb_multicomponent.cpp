/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 Contributing authors: Abhin Suresh (University of Delaware)
------------------------------------------------------------------------- */

// unit tests for fix style multicomponent fluid dynamics

#include "test_main.h"
#include "yaml_writer.h"

#include "gtest/gtest.h"

#include "input.h"
#include "fix.h"
#include "fix_lb_multicomponent.h"
#include "lammps.h"
#include "modify.h"
#include "latboltz_const.h"

#include <array>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <string>
#include <tuple>

using namespace LAMMPS_NS;

namespace LAMMPS_NS {

// Keep access to lattice internals local to this test executable.
class FixLbMulticomponentTestAccess {
public:
    using Populations = std::array<double, 19>;
    struct Site { int x, y, z; };
    struct Moments { double rho, phi, psi, pressure;
                     std::array<double, 3> velocity; };

    // return a interior site of subdomain, avoiding halo sites
    static Site interior_site(const FixLbMulticomponent &fix){
        return {fix.halo_extent[0], fix.halo_extent[1], fix.halo_extent[2]};
    }

    // Populate fix object f_lb, g_lb, k_lb
    static void set_populations(FixLbMulticomponent &fix, Site site,
                                const Populations &f, const Populations &g,
                                const Populations &k)
    {
        for (int i = 0; i < 19; ++i) {
            fix.f_lb[site.x][site.y][site.z][i] = f[i];
            fix.g_lb[site.x][site.y][site.z][i] = g[i];
            fix.k_lb[site.x][site.y][site.z][i] = k[i];
        }
    }

    // Call the unit calc_moments of the fix object
    static Moments calculate_moments(FixLbMulticomponent &fix, Site site)
    {
        fix.calc_moments(site.x, site.y, site.z);
        const auto *u = fix.u_lb[site.x][site.y][site.z];
        return {fix.density_lb[site.x][site.y][site.z],
                fix.phi_lb[site.x][site.y][site.z],
                fix.psi_lb[site.x][site.y][site.z],
                fix.pressure_lb[site.x][site.y][site.z],
                {u[0], u[1], u[2]}};
    }

    // Scalar fields passed to calc_gradient_laplacian by calc_equilibrium
    enum class Field { Rho, Phi, Psi };
    struct Derivatives { std::array<double, 3> gradient;
                         double laplacian; };
    // polynomial q(X, Y, Z) in lattice coordinates relative to the tested site
    using ScalarFunction = double (*)(double, double, double);

    // Fill the field on the 3x3x3 block around the site, which covers the D3Q19 stencil
    static void set_scalar_field(FixLbMulticomponent &fix, Field field, Site site,
                                 ScalarFunction q)
    {
        double ***values = scalar_field(fix, field);
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz)
                    values[site.x + dx][site.y + dy][site.z + dz] = q(dx, dy, dz);
    }

    // Overwrite gradient and laplacian outputs on the 3x3x3 block around the site
    static void fill_derivatives(FixLbMulticomponent &fix, Field field, Site site, double value)
    {
        double ****gradient = gradient_field(fix, field);
        double ***laplacian = laplacian_field(fix, field);
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    const int x = site.x + dx, y = site.y + dy, z = site.z + dz;
                    for (int dir = 0; dir < 3; ++dir) gradient[x][y][z][dir] = value;
                    laplacian[x][y][z] = value;
                }
    }

    // Return copies of the gradient and laplacian stored at the site
    static Derivatives read_derivatives(FixLbMulticomponent &fix, Field field, Site site)
    {
        const double *g = gradient_field(fix, field)[site.x][site.y][site.z];
        return {{g[0], g[1], g[2]}, laplacian_field(fix, field)[site.x][site.y][site.z]};
    }

    // Call the unit calc_gradient_laplacian with the field/output triple used in calc_equilibrium
    static Derivatives calculate_derivatives(FixLbMulticomponent &fix, Field field, Site site)
    {
        fix.calc_gradient_laplacian(site.x, site.y, site.z, scalar_field(fix, field),
                                    gradient_field(fix, field), laplacian_field(fix, field));
        return read_derivatives(fix, field, site);
    }

private:
    static double ***scalar_field(FixLbMulticomponent &fix, Field field)
    {
        switch (field) {
            case Field::Rho: return fix.density_lb;
            case Field::Phi: return fix.phi_lb;
            default: return fix.psi_lb;
        }
    }

    static double ****gradient_field(FixLbMulticomponent &fix, Field field)
    {
        switch (field) {
            case Field::Rho: return fix.density_gradient;
            case Field::Phi: return fix.phi_gradient;
            default: return fix.psi_gradient;
        }
    }

    static double ***laplacian_field(FixLbMulticomponent &fix, Field field)
    {
        switch (field) {
            case Field::Rho: return fix.laplace_rho;
            case Field::Phi: return fix.laplace_phi;
            default: return fix.laplace_psi;
        }
    }
};

} // namespace LAMMPS_NS

class FixLbMulticomponentTest : public ::testing::Test {
protected:
    std::unique_ptr<LAMMPS> lmp;
    FixLbMulticomponent *fix = nullptr;

    void SetUp() override{
        const char *args[] = {"FixLbMulticomponentTest", "-log", "none", "-screen", "none",
                              "-nocite"};
        lmp.reset(new LAMMPS(sizeof(args) / sizeof(args[0]),
                            const_cast<char **>(args), MPI_COMM_WORLD));
        lmp->input->one("boundary p p p");
        lmp->input->one("region fluid block 0 8 0 8 0 8");
        lmp->input->one("create_box 0 fluid");
        lmp->input->one("timestep 1.0");
        // zero kappas reduce the pressure to rho * c_s^2, an independent reference value
        lmp->input->one("fix mcmp all lb/multicomponent 1 0.166667 1.0 D3Q19 dx 1 "
                        "kappa1 0 kappa2 0 kappa3 0 init mixture");
        auto fixes = lmp->modify->get_fix_by_style("lb/multicomponent");
        ASSERT_EQ(fixes.size(), 1);
        fix = dynamic_cast<FixLbMulticomponent *>(fixes[0]);
        ASSERT_NE(fix, nullptr);
    }

    void TearDown() override { 
        lmp.reset(); 
    }
};

// TEST_F(FixtureClassName, TestName)
TEST_F(FixLbMulticomponentTest, ReconstructsDensityCompositionAndVelocity)
{
    using Access = FixLbMulticomponentTestAccess;
    Access::Populations f{}, g{}, k{};
    
    std::array<int, 19> hits{};
    struct Entry { int ex, ey, ez; double f, g, k; };
    // Hand-chosen populations keyed by lattice vector, independent of the e19 ordering.
    // Diagonals are nonzero and opposite pairs are unequal so every direction contributes.
    const Entry table[] = {
        { 0, 0, 0, 1.30,  0.10, 0.20},
        { 1, 0, 0, 0.40,  0.05, 0.00}, {-1, 0, 0, 0.10,  0.00, 0.03},
        { 0, 1, 0, 0.20,  0.00, 0.00}, { 0,-1, 0, 0.50, -0.04, 0.00},
        { 0, 0, 1, 0.60,  0.00, 0.05}, { 0, 0,-1, 0.20,  0.00, 0.00},
        { 1, 1, 0, 0.10,  0.03, 0.00}, {-1,-1, 0, 0.05,  0.00, 0.00},
        { 1,-1, 0, 0.07,  0.00, 0.00}, {-1, 1, 0, 0.03,  0.00, 0.04},
        { 1, 0, 1, 0.08,  0.00, 0.00}, {-1, 0,-1, 0.02,  0.14, 0.00},
        { 1, 0,-1, 0.04,  0.00, 0.06}, {-1, 0, 1, 0.06,  0.00, 0.00},
        { 0, 1, 1, 0.09,  0.00, 0.02}, { 0,-1,-1, 0.01,  0.00, 0.00},
        { 0, 1,-1, 0.05,  0.00, 0.00}, { 0,-1, 1, 0.10,  0.02, 0.00},
    };

    for (const auto &t : table)
        for (int i = 0; i < 19; ++i)
            if (e19[i][0] == t.ex && e19[i][1] == t.ey && e19[i][2] == t.ez) {
                f[i] = t.f;
                g[i] = t.g;
                k[i] = t.k;
                ++hits[i];
            }
    // every D3Q19 direction must be assigned exactly once
    for (int i = 0; i < 19; ++i) ASSERT_EQ(hits[i], 1) << "direction " << i;

    // Get interior lattice site of the subdomain
    const auto site = Access::interior_site(*fix);
    // Set fix object f_lb, g_lb, k_lb values
    Access::set_populations(*fix, site, f, g, k);
    // Call the unit fix.calc_moments and returns density, phi, psi, pressure, velocity
    const auto calculated = Access::calculate_moments(*fix, site);

    constexpr double tolerance = 1e-13;
    // rho = 1.3 (rest) + 2.0 (axial) + 0.7 (diagonal)
    EXPECT_NEAR(calculated.rho, 4.0, tolerance);
    EXPECT_NEAR(calculated.phi, 0.30, tolerance);
    EXPECT_NEAR(calculated.psi, 0.40, tolerance);
    // j = (0.3 + 0.13, -0.3 + 0.04, 0.4 + 0.21), u = j / rho
    EXPECT_NEAR(calculated.velocity[0],  0.1075, tolerance);
    EXPECT_NEAR(calculated.velocity[1], -0.065,  tolerance);
    EXPECT_NEAR(calculated.velocity[2],  0.1525, tolerance);
    // all kappa = 0 in the fixture, so p = rho * c_s^2
    EXPECT_NEAR(calculated.pressure, 4.0 / 3.0, tolerance);
}

// Tests of calc_gradient_laplacian, kept in their own suites for the CTest filter
class FixLbDerivativeTest : public FixLbMulticomponentTest {};

// Quadratic polynomial with its analytic gradient and laplacian at the origin.
// The D3Q19 stencil is exact for polynomials up to second order.
struct PolynomialCase {
    const char *name;
    FixLbMulticomponentTestAccess::ScalarFunction q;
    std::array<double, 3> gradient;
    double laplacian;
};

void PrintTo(const PolynomialCase &c, std::ostream *os) { *os << c.name; }

using Field = FixLbMulticomponentTestAccess::Field;

std::string field_name(Field field)
{
    switch (field) {
        case Field::Rho: return "rho";
        case Field::Phi: return "phi";
        default: return "psi";
    }
}

const PolynomialCase polynomial_cases[] = {
    {"Constant", [](double, double, double) { return 5.0; }, {0, 0, 0}, 0},
    {"X", [](double X, double, double) { return X; }, {1, 0, 0}, 0},
    {"Y", [](double, double Y, double) { return Y; }, {0, 1, 0}, 0},
    {"Z", [](double, double, double Z) { return Z; }, {0, 0, 1}, 0},
    {"Linear", [](double X, double Y, double Z) { return 5 + 2 * X - 3 * Y + 4 * Z; },
     {2, -3, 4}, 0},
    {"XX", [](double X, double, double) { return X * X; }, {0, 0, 0}, 2},
    {"YY", [](double, double Y, double) { return Y * Y; }, {0, 0, 0}, 2},
    {"ZZ", [](double, double, double Z) { return Z * Z; }, {0, 0, 0}, 2},
    {"UnequalCurvature", [](double X, double Y, double Z) { return X * X + 2 * Y * Y + 3 * Z * Z; },
     {0, 0, 0}, 12},
    {"MixedXY", [](double X, double Y, double) { return (X + 1) * (Y + 2); }, {2, 1, 0}, 0},
    {"MixedYZ", [](double, double Y, double Z) { return (Y + 1) * (Z + 2); }, {0, 2, 1}, 0},
    {"MixedZX", [](double X, double, double Z) { return (Z + 1) * (X + 2); }, {1, 0, 2}, 0},
    {"Combined",
     [](double X, double Y, double Z) {
         return 5 + 2 * X - 3 * Y + 4 * Z + X * X + 2 * Y * Y + 3 * Z * Z;
     },
     {2, -3, 4}, 12},
};

class FixLbDerivativePolynomialTest
    : public FixLbDerivativeTest,
      public ::testing::WithParamInterface<std::tuple<Field, PolynomialCase>> {};

TEST_P(FixLbDerivativePolynomialTest, MatchesAnalyticDerivatives)
{
    using Access = FixLbMulticomponentTestAccess;
    const Field field = std::get<0>(GetParam());
    const PolynomialCase &c = std::get<1>(GetParam());
    const auto site = Access::interior_site(*fix);

    // Sentinel outputs: the site must be overwritten (not accumulated into),
    // and its neighbours must be left untouched.
    constexpr double sentinel = 1.0e3;
    Access::set_scalar_field(*fix, field, site, c.q);
    Access::fill_derivatives(*fix, field, site, sentinel);
    const auto calculated = Access::calculate_derivatives(*fix, field, site);

    constexpr double tolerance = 1e-13;
    EXPECT_NEAR(calculated.gradient[0], c.gradient[0], tolerance);
    EXPECT_NEAR(calculated.gradient[1], c.gradient[1], tolerance);
    EXPECT_NEAR(calculated.gradient[2], c.gradient[2], tolerance);
    EXPECT_NEAR(calculated.laplacian, c.laplacian, tolerance);

    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz) {
                if (dx == 0 && dy == 0 && dz == 0) continue;
                const auto neighbour = Access::read_derivatives(
                    *fix, field, {site.x + dx, site.y + dy, site.z + dz});
                EXPECT_EQ(neighbour.gradient, (std::array<double, 3>{sentinel, sentinel, sentinel}))
                    << "neighbour (" << dx << "," << dy << "," << dz << ")";
                EXPECT_EQ(neighbour.laplacian, sentinel)
                    << "neighbour (" << dx << "," << dy << "," << dz << ")";
            }
}

INSTANTIATE_TEST_SUITE_P(
    Polynomials, FixLbDerivativePolynomialTest,
    ::testing::Combine(::testing::Values(Field::Rho, Field::Phi, Field::Psi),
                       ::testing::ValuesIn(polynomial_cases)),
    [](const ::testing::TestParamInfo<FixLbDerivativePolynomialTest::ParamType> &info) {
        return field_name(std::get<0>(info.param)) + "_" + std::get<1>(info.param).name;
    });

// Negative control: for q = X^3 the continuum gradient at the origin is zero, but the
// second-order stencil gives 3 * sum_i w_i e_ix^4 = 3 * (2/18 + 8/36) = 1. This pins the
// expected truncation error and shows the analytic cases are not trivially satisfied.
TEST_F(FixLbDerivativeTest, CubicShowsSecondOrderTruncation)
{
    using Access = FixLbMulticomponentTestAccess;
    const auto site = Access::interior_site(*fix);
    Access::set_scalar_field(*fix, Field::Phi, site,
                             [](double X, double, double) { return X * X * X; });
    const auto calculated = Access::calculate_derivatives(*fix, Field::Phi, site);

    constexpr double tolerance = 1e-13;
    EXPECT_NEAR(calculated.gradient[0], 1.0, tolerance);
    EXPECT_NEAR(calculated.gradient[1], 0.0, tolerance);
    EXPECT_NEAR(calculated.gradient[2], 0.0, tolerance);
    EXPECT_NEAR(calculated.laplacian, 0.0, tolerance);
}

LAMMPS *init_lammps()
{
    const char *args[] = {"FixLBMulticomponent", "-log", "none", 
                          "-echo", "screen", "-nocite"};
    
    char **argv = (char **)args;
    int argc = sizeof(args) / sizeof(char *);

    LAMMPS *lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);

    // utility lambda to improve readability
    auto command = [&](const std::string &line){
        lmp->input->one(line); 
    };
    
    command("region fluid block -16 16 -16 16 -3 3");
    command("create_box 0 fluid");
    command("timestep 1.0");
    command("fix mcmp all lb/multicomponent "
            "1 0.166667 1.0 D3Q19 dx 1 "
            "C1 0.333333 C2 0.333333 C3 0.333334 "
            "kappa1 0.01 kappa2 0.02 kappa3 0.05 "
            "init mixture");

    return lmp;
}

void run_lammps(LAMMPS *lmp)
{
    lmp->input->one("run 100");
}

void generate_yaml_file(const char *outfile, const TestConfig &config)
{
    LAMMPS *lmp = init_lammps();

    YamlWriter writer(outfile);

    write_yaml_header(&writer, &test_config, lmp->version);

    delete lmp;
} 
// TEST(TestSuitName, TestName)
TEST(FixLBMulticomponent, plain)
{
    LAMMPS *lmp = init_lammps();
    
    // Check for LAMMPS initailization with the fix, else skip the test
    if (!lmp) {                                                       
      std::cerr << "One or more prerequisite styles are not available "
                   "in this LAMMPS configuration:\n";               
      for (auto &prerequisite : test_config.prerequisites)         
        std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
                                                                    
      GTEST_SKIP();                                                 
    }                                                                 
                                                                   
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // get the base-class Fix * matching the style
    auto fixes = lmp->modify->get_fix_by_style("lb/multicomponent");    
        
    ASSERT_EQ(fixes.size(), 1);
    ASSERT_NE(fixes[0], nullptr);

    // cast to FixLbMulticomponent * is it exists
    auto *fix = dynamic_cast<FixLbMulticomponent *>(fixes[0]);
    
    // verify the cast to FixLbMulticomponent succeeded
    ASSERT_NE(fix, nullptr);
    
    // obtain initial total momentum
    double jx0, jy0, jz0, m0;
    jx0 = jy0 = jz0 = m0 = 0.0;
    fix->get_total_momentum(jx0, jy0, jz0);
    fix->get_total_mass(m0);
    if (rank == 0){
      std::cout << "Initial P: " << jx0 << " " << jy0 << " " << jz0 << std::endl;
      std::cout << "Initial M: " << m0 << std::endl;
      EXPECT_NEAR(jx0, 0, 1e-6);
      EXPECT_NEAR(jy0, 0, 1e-6);
      EXPECT_NEAR(jz0, 0, 1e-6);
    }
    
    // runs lammps time steps invoking fix_lb_multicomponent
    run_lammps(lmp); 

    // obtain final total momentum
    double jx1, jy1, jz1, m1;
    jx1 = jy1 = jz1 = m1 = 0.0;
    fix->get_total_momentum(jx1, jy1, jz1);
    fix->get_total_mass(m1);
    
    // verify conservation total momentum
    if (rank == 0){
      std::cout << "Final P: " << jx1 << " " << jy1 << " " << jz1 << std::endl;
      std::cout << "Final M: " << m1 << std::endl;
      EXPECT_NEAR(jx1, 0, 1e-6);
      EXPECT_NEAR(jy1, 0, 1e-6);
      EXPECT_NEAR(jz1, 0, 1e-6);
      EXPECT_NEAR(m0, m1, 1e-6);
    }
    delete lmp;
}
