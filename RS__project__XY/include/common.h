//
// Created by 11379 on 2020-10-30.
//

#ifndef FRAME_COMMON_H
#define FRAME_COMMON_H

#include "classAndVarDefine.h"
#include <stdlib.h>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <string>
#include <map>
#include <string.h>
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <numeric>
#include <random>
#define PrecisionValue 1e-6
#define InfInt 9999999999
#define VALUE 1000000
#define XY_MIN(i, j)   (((i) > (j)) ? (j) : (i))
#define XY_MAX(i, j)   (((i) < (j)) ? (j) : (i))
#define XY_SWAP(x, y, type) {type tmp = (x); (x) = (y); (y) = (tmp);}
#endif //FRAME_COMMON_H
