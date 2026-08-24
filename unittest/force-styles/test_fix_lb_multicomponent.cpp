#include "gtest/gtest.h"

#include "input.h"
#include "lammps.h"

#include <mpi.h>

#include <string>

using namespace LAMMPS_NS;

TEST(FixLBMulticomponent, plain)
{
    const char *args[] = {"FixLBMulticomponent", "-log", "none", 
                          "-echo", "screen", "-nocite"};

    char **argv = (char **)args;
    int argc = sizeof(args) / sizeof(char *);

    LAMMPS *lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);

    auto command = [&](const std::string &line)
    {
        lmp->input->one(line);
    };
    
    command("region fluid block -16 16 -16 16 -3 3");
    command("create_box 0 fluid");

    command("timestep 1.0");

    command(
        "fix mcmp all lb/multicomponent "
        "1 0.166667 1.0 D3Q19 dx 1 "
        "C1 0.333333 C2 0.333333 C3 0.333334 "
        "kappa1 0.01 kappa2 0.02 kappa3 0.05 "
        "init mixture"
    );

    command("run 10");

    delete lmp;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    MPI_Finalize();
    return result;
}
