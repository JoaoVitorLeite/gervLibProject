#ifndef UTIL_H
#define UTIL_H

#include <cstdlib>
#include <stdexcept>
#include <unordered_set>
#include <iostream>
#include <vector>
#include <dirent.h>
#include <algorithm>
#include <regex>
#include <queue>
#include <vector>

size_t* uniqueRandomNumber(size_t start, size_t end, size_t nNumber, size_t seed = 0);



size_t* randomNumber(size_t start, size_t end, size_t nNumber, size_t seed = 0);



size_t* shuffleIndex(size_t start, size_t end, size_t seed = 0);



std::string selectFolderName(std::string directory, std::string startwith);



#endif // UTIL_H
