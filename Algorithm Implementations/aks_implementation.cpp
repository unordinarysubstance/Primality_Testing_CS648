#include <iostream>
#include <vector>
#include <cmath>
#include <gmp.h>
#include <gmpxx.h>
#include<numeric>
#include <cstring>

bool perfectPower(const mpz_class& n) {
    if (n < 2) return false;

    unsigned long max_b = mpz_sizeinbase(n.get_mpz_t(), 2); // max b = log2(n)

    for (unsigned long b = 2; b <= max_b; ++b) {
        mpz_class root;
        bool is_exact = mpz_root(root.get_mpz_t(), n.get_mpz_t(), b);

        // Check if root^b == n (for edge cases where mpz_root isn't exact)
        mpz_class power;
        mpz_pow_ui(power.get_mpz_t(), root.get_mpz_t(), b);
        if (power == n) {
            return true;
        }

        mpz_class root_plus_1 = root + 1;
        // Also check (root + 1)^b to catch rounding issues
        mpz_pow_ui(power.get_mpz_t(), root_plus_1.get_mpz_t(), b);
        if (power == n) {
            return true;
        }
    }

    return false;
}

mpz_class fastMod(mpz_class base, mpz_class power, const mpz_class& mod) {
    mpz_class result = 1;
    base = base % mod;

    while (power > 0) {
        if (power % 2 == 1) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        power /= 2;
    }

    return result;
}

unsigned long findR(const mpz_class& n) {
    double logn = mpz_sizeinbase(n.get_mpz_t(), 2);
    double maxK = pow(logn, 2);
    double maxR = pow(logn, 5);  // unused in this code, but in original Python
    bool nexR = true;
    unsigned long r = 1;

    while (nexR) {
        r++;
        nexR = false;

        for (unsigned long k = 1; k <= static_cast<unsigned long>(maxK); ++k) {
            mpz_class result = fastMod(n, k, r);
            if (result == 0 || result == 1) {
                nexR = true;
                break;
            }
        }
    }

    return r;
}

std::vector<mpz_class> multi(const std::vector<mpz_class>& a, const std::vector<mpz_class>& b, const mpz_class& n, size_t r) {
    size_t len = a.size() + b.size() - 1;
    std::vector<mpz_class> x(r, 0);  // result always has length r due to mod x^r - 1

    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            size_t idx = (i + j) % r;
            x[idx] = (x[idx] + a[i] * b[j]) % n;
        }
    }

    return x;
}

// --- Fast modular exponentiation for polynomials
std::vector<mpz_class> fastPoly(std::vector<mpz_class> base,mpz_class power,size_t r) {
    std::vector<mpz_class> x(r, 0);
    x[0] = 1;

    mpz_class n = power;
    mpz_class a = base[0];  // constant term of the input polynomial

    while (power > 0) {
        if (power % 2 == 1) {
            x = multi(x, base, n, r);
        }
        base = multi(base, base, n, r);
        power /= 2;
    }

    // x[0] -= a
    x[0] = (x[0] - a) % n;
    if (x[0] < 0) x[0] += n;

    // x[n % r] -= 1
    size_t idx = mpz_class(n % r).get_ui();
    x[idx] = (x[idx] - 1) % n;
    if (x[idx] < 0) x[idx] += n;

    return x;
}

mpz_class gcd(mpz_class a, mpz_class b) {
    mpz_class g;
    mpz_gcd(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return g;
}

int eulerPhi(int r) {
    int count = 0;
    for (int i = 1; i <= r; ++i) {
        if (gcd(i, r) == 1)
            ++count;
    }
    return count;
}

std::string aks(mpz_class n) {
    // Step 1: Check if n is a perfect power
    if (perfectPower(n)) {
        return "Composite";
    }

    // Step 2: Find the smallest r such that order_n(r) > log2(n)^2
    unsigned long r = findR(n);

    // Step 3: Check GCD(a, n) for 2 ≤ a ≤ min(r, n)
    for (mpz_class a = 2; a < std::min((mpz_class)r, n); ++a) {
        if (gcd(a, n) > 1) {
            return "Composite";
        }
    }

    // Step 4: If n ≤ r, then n is prime
    if (n <= r) {
        return "Prime";
    }

    // Step 5: Polynomial congruence check
    mpz_class phi = eulerPhi(r);
    double logn = mpz_sizeinbase(n.get_mpz_t(), 2);
    unsigned long limit = static_cast<unsigned long>(std::floor(std::sqrt(phi.get_d()) * logn));

    for (mpz_class a = 1; a <= limit; ++a) {
        std::vector<mpz_class> poly = {a, 1}; // x + a
        std::vector<mpz_class> result = fastPoly(poly, n, r);

        bool any_nonzero = false;
        for (const auto& coeff : result) {
            if (coeff != 0) {
                any_nonzero = true;
                break;
            }
        }

        if (any_nonzero) {
            return "Composite";
        }
    }

    // Step 6: Passed all tests
    return "Prime";
}

int main() {
    std::string input;
    std::cout << "Enter a number: ";
    std::cin >> input;
    mpz_class n(input);  // Change this to test other numbers
    std::cout << aks(n) << std::endl;
    return 0;
}
