#pragma once

#include <math.h>

namespace cvl{

template<class T> T cubick_ding(T c2, T c1, T c0, T & root) {
    
    T a = c1 - c2 * c2 / 3.0;
    T b = (2.0 * c2 * c2 * c2 - 9.0 * c2 * c1) / 27.0 + c0;
    T c = b * b / 4.0 + a * a * a / 27.0;
    if (c != 0) {
        if (c > 0) {
            c = std::sqrt(c);
            b *= -0.5;
            root = std::cbrt(b + c) + std::cbrt(b - c) - c2 / 3.0;
            return true;
        } else {
            c = 3.0 * b / (2.0 * a) * std::sqrt(-3.0 / a);
            root = 2.0 * std::sqrt(-a / 3.0) * std::cos(std::acos(c) / 3.0) - c2 / 3.0;
        }
    } else {
        root = -c2 / 3.0 + (a != 0 ? (3.0 * b / a) : 0);
    }
    
    return false;
}


}
