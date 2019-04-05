#pragma once
#include <inttypes.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int * currentTime();
void make_directory(const char* name);
void print_path();
