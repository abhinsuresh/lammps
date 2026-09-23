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

using namespace LAMMPS_NS;

namespace LAMMPS_NS {

// Keep access to lattice internals local to this test executable.
class FixLbMulticomponentTestAccess {
public:
    using Populations = std::array<double, 19>;
    struct Site { int x, y, z; };
    struct Moments {
        double rho, phi, psi;
        std::array<double, 3> velocity;
    };

    static Site owned_site(const FixLbMulticomponent &fix)
    {
        return {fix.subNbx / 2, fix.subNby / 2, fix.subNbz / 2};
    }

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

    static Moments reconstruct(FixLbMulticomponent &fix, Site site)
    {
        fix.calc_moments(site.x, site.y, site.z);
        const auto *u = fix.u_lb[site.x][site.y][site.z];
        return {fix.density_lb[site.x][site.y][site.z],
                fix.phi_lb[site.x][site.y][site.z],
                fix.psi_lb[site.x][site.y][site.z], {u[0], u[1], u[2]}};
    }
};

} // namespace LAMMPS_NS

class LBMomentReconstruction : public ::testing::Test {
protected:
    std::unique_ptr<LAMMPS> lmp;
    FixLbMulticomponent *fix = nullptr;

    void SetUp() override
    {
        const char *args[] = {"LBMomentReconstruction", "-log", "none", "-screen", "none",
                              "-nocite"};
        lmp.reset(new LAMMPS(sizeof(args) / sizeof(args[0]),
                            const_cast<char **>(args), MPI_COMM_WORLD));
        lmp->input->one("boundary p p p");
        lmp->input->one("region fluid block 0 8 0 8 0 8");
        lmp->input->one("create_box 0 fluid");
        lmp->input->one("timestep 1.0");
        lmp->input->one("fix mcmp all lb/multicomponent 1 0.166667 1.0 D3Q19 dx 1 init mixture");
        auto fixes = lmp->modify->get_fix_by_style("lb/multicomponent");
        ASSERT_EQ(fixes.size(), 1);
        fix = dynamic_cast<FixLbMulticomponent *>(fixes[0]);
        ASSERT_NE(fix, nullptr);
    }

    void TearDown() override { lmp.reset(); }
};

TEST_F(LBMomentReconstruction, ReconstructsDensityCompositionAndVelocity)
{
    using Access = FixLbMulticomponentTestAccess;
    Access::Populations f{}, g{}, k{};
    // Resolve direction indices from D3Q19, while keeping expected moments explicit.
    for (int i = 0; i < 19; ++i) {
        const int x = e19[i][0], y = e19[i][1], z = e19[i][2];
        if (x == 0 && y == 0 && z == 0) {
            f[i] = 2.0;
            g[i] = 0.3;
            k[i] = 0.4;
        } else if (y == 0 && z == 0) {
            f[i] = x == 1 ? 0.4 : 0.1;
        } else if (x == 0 && z == 0) {
            f[i] = y == 1 ? 0.2 : 0.5;
        } else if (x == 0 && y == 0) {
            f[i] = z == 1 ? 0.6 : 0.2;
        }
    }
    const auto site = Access::owned_site(*fix);
    Access::set_populations(*fix, site, f, g, k);
    const auto actual = Access::reconstruct(*fix, site);

    constexpr double tolerance = 1e-13;
    EXPECT_NEAR(actual.rho, 4.0, tolerance);
    EXPECT_NEAR(actual.phi, 0.3, tolerance);
    EXPECT_NEAR(actual.psi, 0.4, tolerance);
    EXPECT_NEAR(actual.velocity[0], 0.075, tolerance);
    EXPECT_NEAR(actual.velocity[1], -0.075, tolerance);
    EXPECT_NEAR(actual.velocity[2], 0.100, tolerance);
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

/*
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    MPI_Finalize();
    return result;
}
*/
