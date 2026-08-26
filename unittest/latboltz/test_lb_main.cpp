#include "test_lb_main.h"

#include "test_config_reader.h"
#include "yaml_writer.h"

#include "gtest/gtest.h"

#include <iostream>
#include <mpi.h>
#include <string>

TestConfig test_config;

bool read_yaml_file(const char *infile, TestConfig &config)
{
    auto reader = TestConfigReader(config);

    if (reader.parse_file(infile))
        return false;

    config.basename = reader.get_basename();
    return true;
}


void write_yaml_header(YamlWriter *writer,
                       const TestConfig *config,
                       const char *version)
{
    writer->emit("lammps_version", version);
    writer->emit("epsilon", config->epsilon);

    std::string block;

    for (auto &command : config->pre_commands)
        block += command + "\n";

    writer->emit_block("pre_commands", block);

    block.clear();

    for (auto &command : config->post_commands)
        block += command + "\n";

    writer->emit_block("post_commands", block);
}


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    ::testing::InitGoogleTest(&argc, argv);

    if (argc == 1){
        int result = RUN_ALL_TESTS();
        MPI_Finalize();
        return result;
    }

    if (std::string(argv[1]) == "--generate") {
        if (argc < 3) {
            std::cerr << "Missing YAML output filename\n";
            MPI_Finalize();
            return 1;
        }

        generate_yaml_file(argv[2], test_config);

        MPI_Finalize();
        return 0;
    }

    if (!read_yaml_file(argv[1], test_config)) {
        std::cerr << "Error reading YAML file: "
                  << argv[1] << "\n";

        MPI_Finalize();
        return 2;
    }

    int result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}
