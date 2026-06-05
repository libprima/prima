#ifndef PRIMA_CPP_TEST_COMMON_HPP
#define PRIMA_CPP_TEST_COMMON_HPP

#include <iostream>
#include <string>
#include <utility>

namespace prima::test {

struct TestRunner {
    int failures = 0;

    void check(bool cond, const std::string& msg) {
        if (!cond) {
            std::cerr << "FAIL: " << msg << std::endl;
            ++failures;
        } else {
            std::cout << "PASS: " << msg << std::endl;
        }
    }

    void check_close(double actual, double expected, double tol, const std::string& msg) {
        if (std::abs(actual - expected) < tol) {
            std::cout << "PASS: " << msg << std::endl;
        } else {
            std::cerr << "FAIL: " << msg << " (got " << actual << ", expected " << expected << ", tol " << tol << ")" << std::endl;
            ++failures;
        }
    }

    int done() {
        if (failures > 0) {
            std::cerr << "\n" << failures << " test(s) FAILED!" << std::endl;
            return 1;
        }
        std::cout << "\nAll tests PASSED." << std::endl;
        return 0;
    }
};

} // namespace prima::test

#endif // PRIMA_CPP_TEST_COMMON_HPP
