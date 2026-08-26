#ifndef LMP_TEST_LB_MAIN_H
#define LMP_TEST_LB_MAIN_H

#include "test_config.h"

class YamlWriter;

extern TestConfig test_config;
// extern std::string INPUT_FOLDER;

bool read_yaml_file(const char *, TestConfig &);

void write_yaml_header(YamlWriter *, const TestConfig *, const char *);

void generate_yaml_file(const char *, const TestConfig &);

#endif
