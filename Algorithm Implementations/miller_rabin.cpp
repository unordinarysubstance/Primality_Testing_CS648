#include <iostream>
#include <cmath>
#include <chrono>
#include <gmp.h>
#include <ctime>

// Perform (base^exp) % mod using GMP
void mod_exp(mpz_t result, const mpz_t base, const mpz_t exp, const mpz_t mod) {
    mpz_powm(result, base, exp, mod);
}

// Returns true if n passes one Miller-Rabin test with base a
bool miller_test(const mpz_t n, const mpz_t d, const mpz_t a) {
    mpz_t x;
    mpz_init(x);
    mod_exp(x, a, d, n);

    mpz_t n_minus_1;
    mpz_init(n_minus_1);
    mpz_sub_ui(n_minus_1, n, 1);

    if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, n_minus_1) == 0) {
        mpz_clear(x);
        // mpz_clear(n_minus_1);
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
            return false;
        }

        if (mpz_cmp(x, n_minus_1) == 0) {
            mpz_clear(x);
            mpz_clear(temp_d);
            return true;
        }
    }

    mpz_clear(x);
    mpz_clear(temp_d);
    mpz_clear(n_minus_1);
    return false;
}

// Custom Miller-Rabin with log-log number of rounds
bool is_prime_miller_rabin(const mpz_t n, int k = -1) {
    size_t num_digits = mpz_sizeinbase(n, 10);
    if (k == -1) {
        double loglog = std::log2(static_cast<double>(num_digits));
        k = std::max(5, 2 * static_cast<int>(std::ceil(loglog)));
    }

    if (mpz_cmp_ui(n, 2) == 0 || mpz_cmp_ui(n, 3) == 0)
        return true;
    if (mpz_cmp_ui(n, 1) <= 0 || mpz_even_p(n))
        return false;

    mpz_t d;
    mpz_init_set(d, n);
    mpz_sub_ui(d, d, 1);  // d = n - 1

    while (mpz_even_p(d))
        mpz_divexact_ui(d, d, 2);

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
            return false;
        }
    }

    mpz_clear(d);
    mpz_clear(a);
    gmp_randclear(state);
    return true;
}

int main() {
    mpz_t num;
    
    std::string str;
    std::cout << "Enter a number to test for primality: ";
    std::getline(std::cin, str);
    // str = str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
    mpz_init_set_str(num, str.c_str(), 10);

    // Custom Miller-Rabin Benchmark
    auto start_custom = std::chrono::high_resolution_clock::now();
    bool result_custom = is_prime_miller_rabin(num);
    auto end_custom = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_custom = end_custom - start_custom;

    std::cout << "[Custom Miller-Rabin] ";
    std::cout << (result_custom ? "Probably Prime!" : "Composite.") << "\n";
    std::cout << "Time: " << elapsed_custom.count() << " seconds\n\n";

    size_t num_digits = mpz_sizeinbase(num, 10);
    double loglog = std::log2(num_digits);
    int k = std::max(5, 2*static_cast<int>(std::ceil(loglog)));

    // GMP Built-in Benchmark
    auto start_gmp = std::chrono::high_resolution_clock::now();
    int result_gmp = mpz_probab_prime_p(num, k);
    auto end_gmp = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_gmp = end_gmp - start_gmp;

    std::cout << "[GMP Built-in] ";
    if (result_gmp == 2)
        std::cout << "Definitely Prime!";
    else if (result_gmp == 1)
        std::cout << "Probably Prime!";
    else
        std::cout << "Composite.";
    std::cout << "\nTime: " << elapsed_gmp.count() << " seconds\n";

    mpz_clear(num);
    return 0;
}