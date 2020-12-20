//
// Created by Zhuoqing Fang on 7/8/20.
//

#include <numeric>
#include <algorithm>
#include <functional>
#include "linalg.h"

// vector math
float dot_product(std::vector<float> & v1, std::vector<float> & v2) {
    assert(v1.size() == v2.size());
    float res = std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0F);
    return res;
}


std::vector<float> add_vectors(std::vector<float> & v1, std::vector<float> & v2) {
    assert(v1.size() == v2.size());
    std::vector<float> res(v1.size(), 0.0F);
    // res[0] = v1[0] + v2[0]
    // ...
    std::transform(v1.begin(), v1.end(), v2.begin(), res.begin(),
                   [](int x, int y){return x + y;});
    return res;
}

std::vector<float> sub_vectors(std::vector<float> & v1, float c) {
    std::vector<float> v2(v1.size(),0.0);
    std::transform(v1.begin(), v1.end(), v2.begin(),
                   [&c](float x){return x - c;}); // [&c] or [c], capture variable for closure
    return v2;
}
