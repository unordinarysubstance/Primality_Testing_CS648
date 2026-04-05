#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <vector>
#include <cmath>
#include <gmp.h>
#include <gmpxx.h>
#include <numeric>
#include <cstring>

using namespace std;
using namespace std::chrono;

// ==========================================
// 1. MILLER-RABIN IMPLEMENTATION (Fast, Probabilistic)
// ==========================================

uint64_t power(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t res = 1;
    base = base % mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (__int128)res * base % mod;
        exp = exp >> 1;
        base = (__int128)base * base % mod;
    }
    return res;
}

bool millerRabinTest(uint64_t d, uint64_t n) {
    uint64_t a = 2 + rand() % (n - 4);
    uint64_t x = power(a, d, n);
    if (x == 1 || x == n - 1) return true;
    while (d != n - 1) {
        x = (__int128)x * x % n;
        d *= 2;
        if (x == 1) return false;
        if (x == n - 1) return true;
    }
    return false;
}

bool isPrimeMillerRabin(uint64_t n, int k) {
    if (n <= 1 || n == 4) return false;
    if (n <= 3) return true;
    uint64_t d = n - 1;
    while (d % 2 == 0) d /= 2;
    for (int i = 0; i < k; i++) {
        if (!millerRabinTest(d, n)) return false;
    }
    return true;
}

// ==========================================
// 2. AKS IMPLEMENTATION (Slow, Deterministic)
// ==========================================

bool perfectPower(const mpz_class& n) {
    if (n < 2) return false;
    unsigned long max_b = mpz_sizeinbase(n.get_mpz_t(), 2); 
    for (unsigned long b = 2; b <= max_b; ++b) {
        mpz_class root;
        bool is_exact = mpz_root(root.get_mpz_t(), n.get_mpz_t(), b);
        mpz_class power;
        mpz_pow_ui(power.get_mpz_t(), root.get_mpz_t(), b);
        if (power == n) return true;
        mpz_class root_plus_1 = root + 1;
        mpz_pow_ui(power.get_mpz_t(), root_plus_1.get_mpz_t(), b);
        if (power == n) return true;
    }
    return false;
}

mpz_class fastMod(mpz_class base, mpz_class power, const mpz_class& mod) {
    mpz_class result = 1;
    base = base % mod;
    while (power > 0) {
        if (power % 2 == 1) result = (result * base) % mod;
        base = (base * base) % mod;
        power /= 2;
    }
    return result;
}

unsigned long findR(const mpz_class& n) {
    double logn = mpz_sizeinbase(n.get_mpz_t(), 2);
    double maxK = pow(logn, 2);
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
    std::vector<mpz_class> x(r, 0); 
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            size_t idx = (i + j) % r;
            x[idx] = (x[idx] + a[i] * b[j]) % n;
        }
    }
    return x;
}

std::vector<mpz_class> fastPoly(std::vector<mpz_class> base, mpz_class power, size_t r) {
    std::vector<mpz_class> x(r, 0);
    x[0] = 1;
    mpz_class n = power;
    mpz_class a = base[0]; 
    while (power > 0) {
        if (power % 2 == 1) x = multi(x, base, n, r);
        base = multi(base, base, n, r);
        power /= 2;
    }
    x[0] = (x[0] - a) % n;
    if (x[0] < 0) x[0] += n;
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
        if (gcd(i, r) == 1) ++count;
    }
    return count;
}

// Wrapped to return bool for easier benchmarking
bool isPrimeAKS(mpz_class n) {
    if (perfectPower(n)) return false;
    unsigned long r = findR(n);
    for (mpz_class a = 2; a < std::min((mpz_class)r, n); ++a) {
        if (gcd(a, n) > 1) return false;
    }
    if (n <= r) return true;
    mpz_class phi = eulerPhi(r);
    double logn = mpz_sizeinbase(n.get_mpz_t(), 2);
    unsigned long limit = static_cast<unsigned long>(std::floor(std::sqrt(phi.get_d()) * logn));

    for (mpz_class a = 1; a <= limit; ++a) {
        std::vector<mpz_class> poly = {a, 1}; 
        std::vector<mpz_class> result = fastPoly(poly, n, r);
        bool any_nonzero = false;
        for (const auto& coeff : result) {
            if (coeff != 0) {
                any_nonzero = true;
                break;
            }
        }
        if (any_nonzero) return false;
    }
    return true;
}


// ==========================================
// 3. BENCHMARKING ENGINE
// ==========================================

uint64_t generateRandomBitLength(int bits) {
    uint64_t min_val = 1ULL << (bits - 1);
    uint64_t max_val = (1ULL << bits) - 1;
    std::random_device rd;
    std::mt19937_64 eng(rd());
    std::uniform_int_distribution<uint64_t> distr(min_val, max_val);
    return distr(eng);
}

int main() {
    ofstream file("benchmark_data.csv");
    file << "bit_length,iterations_k,miller_rabin_time_us,aks_time_us\n";

    cout << "Starting Benchmark..." << endl;

    // WARNING: AKS is O(b^7.5). Keep the max bits low initially (e.g., up to 20 or 24 bits) 
    // to test if it completes in a reasonable time. You can increase it later.
    for (int b = 4; b <= 20; b += 2) { 
        uint64_t test_num = generateRandomBitLength(b);
        mpz_class test_num_gmp(std::to_string(test_num)); // Convert via string to avoid GMP ambiguity
        
        cout << "Testing " << b << "-bit number: " << test_num << "..." << flush;

        for (int k = 1; k <= 40; k += 5) {
            
            // Measure Miller-Rabin
            auto start_mr = high_resolution_clock::now();
            isPrimeMillerRabin(test_num, k);
            auto end_mr = high_resolution_clock::now();
            auto duration_mr = duration_cast<microseconds>(end_mr - start_mr).count();

            // Measure AKS (We only really need to run this once per number since it doesn't use 'k', 
            // but running it inside the loop simulates the direct comparison for the CSV)
            auto start_aks = high_resolution_clock::now();
            isPrimeAKS(test_num_gmp);
            auto end_aks = high_resolution_clock::now();
            auto duration_aks = duration_cast<microseconds>(end_aks - start_aks).count();

            file << b << "," << k << "," << duration_mr << "," << duration_aks << "\n";
        }
        cout << " Done." << endl;
    }

    file.close();
    cout << "Benchmarking complete. Data saved to benchmark_data.csv" << endl;
    return 0;
}