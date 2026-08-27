#include "test_lb_main.h"

#include "test_config_reader.h"
#include "yaml_writer.h"
#include "gtest/gtest.h"
#include "utils.h"

#include <iostream>
#include <mpi.h>
#include <string>
#include <ctime>

using LAMMPS_NS::utils::trim;

// test configuration settings read from yaml file
TestConfig test_config;

// common read_yaml_file function
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
    // lammps version
    writer->emit("lammps_version", version);
    
    // date_generated
    std::time_t now   = time(nullptr);
    std::string block = trim(ctime(&now));
    writer->emit("date_generated", block);
    writer->emit("epsilon", config->epsilon);

    /*
    block.clear();
    std::string block;
    
    // pres_commands
    for (auto &command : config->pre_commands)
        block += command + "\n";

    writer->emit_block("pre_commands", block);

    block.clear();
    
    // post_commands
    for (auto &command : config->post_commands)
        block += command + "\n";
26-2004080931
    writer->emit_block("post_commands", block);
    */
    writer->emit("input_file", config->input_file);
}


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    ::testing::InitGoogleTest(&argc, argv);

    // if no args provided
    if (argc == 1){
        int result = RUN_ALL_TESTS();
        MPI_Finalize();
        return result;
    }

    // if generate yaml is the 2nd arg, need 3 args
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

    // if yaml file provided but can't read
    if (!read_yaml_file(argv[1], test_config)) {
        std::cerr << "Error reading YAML file: "
                  << argv[1] << "\n";

        MPI_Finalize();
        return 2;
    }
    
    // if one arg provided, which can be read by yaml
    int result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}
