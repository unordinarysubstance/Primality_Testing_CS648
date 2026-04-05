#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <gmp.h>
#include <gmpxx.h>
#include <numeric>
#include <cstring>
#include <ctime>

// ============================================================================
// AKS IMPLEMENTATION (From your aks_implementation.cpp)
// ============================================================================

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

bool is_prime_aks(mpz_class n) {
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
        for (const auto& coeff : result) {
            if (coeff != 0) return false;
        }
    }
    return true;
}


// ============================================================================
// MILLER-RABIN IMPLEMENTATION (From your miller_rabin.cpp)
// ============================================================================

void mod_exp(mpz_t result, const mpz_t base, const mpz_t exp, const mpz_t mod) {
    mpz_powm(result, base, exp, mod);
}

bool miller_test(const mpz_t n, const mpz_t d, const mpz_t a) {
    mpz_t x;
    mpz_init(x);
    mod_exp(x, a, d, n);

    mpz_t n_minus_1;
    mpz_init(n_minus_1);
    mpz_sub_ui(n_minus_1, n, 1);

    if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, n_minus_1) == 0) {
        mpz_clear(x);
        mpz_clear(n_minus_1);
        return true;
    }

    mpz_t temp_d;
    mpz_init_set(temp_d, d);

    while (mpz_cmp(temp_d, n) < 0) {
        mpz_mul_ui(temp_d, temp_d, 2);
        mod_exp(x, a, temp_d, n);

        if (mpz_cmp_ui(x, 1) == 0) {
            mpz_clear(x);
            mpz_clear(temp_d);
            mpz_clear(n_minus_1);
            return false;
        }

        if (mpz_cmp(x, n_minus_1) == 0) {
            mpz_clear(x);
            mpz_clear(temp_d);
            mpz_clear(n_minus_1);
            return true;
        }
    }

    mpz_clear(x);
    mpz_clear(temp_d);
    mpz_clear(n_minus_1);
    return false;
}

bool is_prime_miller_rabin(mpz_class n_class, int k = 20) {
    mpz_t n;
    mpz_init_set(n, n_class.get_mpz_t());

    if (mpz_cmp_ui(n, 2) == 0 || mpz_cmp_ui(n, 3) == 0) { mpz_clear(n); return true; }
    if (mpz_cmp_ui(n, 1) <= 0 || mpz_even_p(n)) { mpz_clear(n); return false; }

    mpz_t d;
    mpz_init_set(d, n);
    mpz_sub_ui(d, d, 1);

    while (mpz_even_p(d)) mpz_divexact_ui(d, d, 2);

    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, std::time(nullptr));

    mpz_t a;
    mpz_init(a);

    for (int i = 0; i < k; ++i) {
        mpz_sub_ui(a, n, 3);
        mpz_urandomm(a, state, a);
        mpz_add_ui(a, a, 2);

        if (!miller_test(n, d, a)) {
            mpz_clear(d);
            mpz_clear(a);
            gmp_randclear(state);
            mpz_clear(n);
            return false;
        }
    }

    mpz_clear(d);
    mpz_clear(a);
    gmp_randclear(state);
    mpz_clear(n);
    return true;
}


// ============================================================================
// RSA PRIME GENERATION BENCHMARK
// ============================================================================

// Generates a random n-bit integer, ensuring the highest and lowest bits are 1
mpz_class generate_random_odd_candidate(int bits, gmp_randstate_t state) {
    mpz_class num;
    mpz_urandomb(num.get_mpz_t(), state, bits);
    mpz_setbit(num.get_mpz_t(), bits - 1); // Ensure it's exactly 'bits' long
    mpz_setbit(num.get_mpz_t(), 0);        // Ensure it's odd
    return num;
}

void run_benchmark(int bit_length, bool test_aks = false) {
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, std::time(nullptr));

    std::cout << "=================================================\n";
    std::cout << " Generating " << bit_length << "-bit RSA Prime\n";
    std::cout << "=================================================\n";

    // --- MILLER-RABIN BENCHMARK ---
    int mr_candidates = 0;
    mpz_class p_mr;
    auto start_mr = std::chrono::high_resolution_clock::now();
    
    while (true) {
        mr_candidates++;
        p_mr = generate_random_odd_candidate(bit_length, state);
        if (is_prime_miller_rabin(p_mr, 20)) {
            break;
        }
    }
    auto end_mr = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_mr = end_mr - start_mr;

    std::cout << "[Miller-Rabin]\n";
    std::cout << "  Found Prime: " << p_mr.get_str(10) << "\n";
    std::cout << "  Candidates Tested: " << mr_candidates << "\n";
    std::cout << "  Time Taken: " << elapsed_mr.count() << " seconds\n\n";

    // --- AKS BENCHMARK ---
    if (test_aks) {
        int aks_candidates = 0;
        mpz_class p_aks;
        auto start_aks = std::chrono::high_resolution_clock::now();
        
        while (true) {
            aks_candidates++;
            p_aks = generate_random_odd_candidate(bit_length, state);
            if (is_prime_aks(p_aks)) {
                break;
            }
        }
        auto end_aks = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_aks = end_aks - start_aks;

        std::cout << "[AKS]\n";
        std::cout << "  Found Prime: " << p_aks.get_str(10) << "\n";
        std::cout << "  Candidates Tested: " << aks_candidates << "\n";
        std::cout << "  Time Taken: " << elapsed_aks.count() << " seconds\n";
        std::cout << "  Speedup (MR / AKS): " << elapsed_aks.count() / elapsed_mr.count() << "x\n";
    } else {
        std::cout << "[AKS]\n  Skipped (Enable via test_aks boolean). Warning: Only enable for <= 16 bits!\n";
    }

    gmp_randclear(state);
}

int main() {
    std::cout << "RSA Prime Generation Benchmark\n\n";

    // We start with extremely small sizes because your AKS algorithm 
    // calculates polynomial expansions without Fast Fourier Transforms (FFT).
    
    run_benchmark(16, true);  // 16-bit: AKS might take a few seconds
    
    // WARNING: Changing the below flag to true will likely lock up your CPU for a long time.
    run_benchmark(32, false); // 32-bit: AKS will be extremely slow
    
    run_benchmark(512, false); // 512-bit (Real RSA-1024 prime): MR takes milliseconds, AKS will not finish.

    return 0;
}