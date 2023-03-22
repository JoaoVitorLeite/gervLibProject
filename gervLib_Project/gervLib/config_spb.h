#ifndef CONFIG_SPB_H
#define CONFIG_SPB_H

#include <cstddef>
#include <string>
#include <filesystem>

typedef unsigned long long ull;
typedef long long ll;

static double MAXDIST = 0.0;
static double EPS = 0.0;
static double MAXINT = 0.0;
static ull p = 0;
static ull PIVOT_NUM = 0; //numero de pivos
static ull GRID_L = ((1ull << p) - 1);
static ull PREC = 1e5;
static const size_t PAGE_SIZE = 256;
static const size_t BTREE_PAGE_SIZE = 4096;
static const size_t MIN_BTREE_LEAF_NUM = 55+1;
static const size_t MIN_INNER_SIZE = 2048;
static long long IOwrite = 0;
static long long IOread = 0;
//static std::string baseFilePath = std::filesystem::current_path().string();
static std::string baseFilePath = "../SPB/spb_files";


#endif // CONFIG_SPB_H
