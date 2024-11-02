#ifndef LAT_FACT
#define LAT_FACT

#pragma GCC target("avx2")

#include <iostream>
#include <tuple>
#include <vector>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ.h>

typedef long long ll;

std::vector<double> zero_vector_RR(const long n);
std::vector<long> zero_vector_ZZ(const long n);
std::vector<std::vector<long>> vector_all(const long n);
void GenSemiPrime(NTL::ZZ& n);
bool is_zero(const std::vector<long> v);
std::vector<std::vector<double>> identity_mat(const long n);
std::vector<std::vector<long>> identity_mat_long(const long n);
ll dot(const std::vector<ll> x, const std::vector<ll> y);
double dot(const std::vector<ll> x, const std::vector<double> y);
double dot(const std::vector<double> x, const std::vector<double> y);
std::vector<std::vector<long>> trans(const std::vector<std::vector<long>> A);
std::vector<std::vector<ll>> trans(const std::vector<std::vector<ll>> A);
std::vector<double> mul_vec_mat(const std::vector<ll> x, const std::vector<std::vector<double>> A);
std::vector<std::vector<double>> mul_mat(const std::vector<std::vector<ll>> A, const std::vector<std::vector<ll>> B);
std::vector<std::vector<double>> mul_mat(const std::vector<std::vector<ll>> A, const std::vector<std::vector<double>> B);
bool have_1(const std::vector<long> v, const long j);
void aug_gauss_mod2(std::vector<std::vector<long>>& A, std::vector<std::vector<long>>& B);
std::vector<std::vector<long>> kernel_basis(const std::vector<std::vector<long>> A);
std::vector<std::vector<long>> kernel_mod2(const std::vector<std::vector<long>> A);
std::vector<std::vector<double>> InverseMatrix(std::vector<std::vector<double>> a);
std::vector<long> factor_basis(const long n);
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> Gram_Schmidt(const std::vector<std::vector<ll>> b);
std::tuple<std::vector<double>, std::vector<std::vector<double>>> Gram_Schmidt_squared(const std::vector<std::vector<ll>> b);
double babai_error(const std::vector<std::vector<ll>> b, const std::vector<ll> w);
std::vector<double> coef2lat(const std::vector<ll> v, const std::vector<std::vector<ll>> b);
std::vector<double> lat2coef(const std::vector<double> v, const std::vector<std::vector<ll>> b);
std::vector<ll> ENUM_CVP(const std::vector<std::vector<double>> mu, const std::vector<double> B, double R, const std::vector<double> a);
std::vector<std::vector<ll>> ENUM_CVP_all(const std::vector<std::vector<double>> mu, const std::vector<double> B, double R, const std::vector<double> a, const std::vector<ll> t, const std::vector<std::vector<ll>> b);
std::tuple<std::vector<std::vector<ll>>, NTL::mat_ZZ> imp_prime_mat(const double c, const std::vector<long> p, const long n);
std::vector<ll> target(const NTL::ZZ N, const double c, const long n);
bool sr_test(const NTL::ZZ u, const NTL::ZZ v, const NTL::ZZ N, const std::vector<long> p, const long n);
extern "C" ll *lattice_factorization(const int bit_flag, const char* NN, const int info_flag);

#endif // !LAT_FACT