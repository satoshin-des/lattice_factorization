#include <iostream>
#include <vector>
#include <random>
#include <iomanip>
#include <tuple>
#include <algorithm>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

#include "lat_fact.hpp"

/**
 * @brief Creates a zero-vector.
 * 
 * @param n Length of the vector
 * @return std::vector<double> zero-vector whose length is n
 */
std::vector<double> zero_vector_RR(const long n){const std::vector<double> z(n); return z;}
std::vector<long> zero_vector_ZZ(const long n){const std::vector<long> z(n); return z;}


/**
 * @brief Generates all combinations of n 0 and 1
 * 
 * @param n length
 * @return std::vector<std::vector<long>> All conbinations
 */
std::vector<std::vector<long>> vector_all(const long n){
    std::vector<long> tmp(n), v(n);
    std::vector<std::vector<long>> v_all;
    for(int i = 0; i < n; ++i){
        tmp[i] = 1; v = tmp;
        do v_all.emplace_back(v); while(std::next_permutation(v.begin(), v.end()));
    }
    return v_all;
}


/**
 * @brief Generates a random semi-prime number
 * 
 * @param n 
 */
void GenSemiPrime(NTL::ZZ& n){
    const int m = NTL::to_int(n) / 2, k = NTL::to_int(n);
    for(;;){
        n = NTL::RandomPrime_ZZ(m) * NTL::RandomPrime_ZZ(NTL::to_int(k) - m + 1);
        if(NTL::NumBits(n) == k) break;
    }
}


/**
 * @brief Tests v is zero-vector or not
 * 
 * @param v Vector
 * @return true if v is zero-vector 
 * @return false if v is non-zero-vector
 */
bool is_zero(const std::vector<long> v){
    const int n = v.size();
    for(int i = 0; i < n; ++i) if(v[i] != 0) return false;
    return true;
}


/**
 * @brief Generates an identity matrix
 * 
 * @param n size
 * @return std::vector<std::vector<double>> identity matrix
 */
std::vector<std::vector<double>> identity_mat(const long n){
    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    for(int i = 0; i < n; ++i) A[i][i] = 1.0;
    return A;
}
std::vector<std::vector<long>> identity_mat_long(const long n){
    std::vector<std::vector<long>> A(n, std::vector<long>(n));
    for(int i = 0; i < n; ++i) A[i][i] = 1;
    return A;
}


/**
 * @brief Computes inner product of vectors x and y
 * 
 * @param x a vector
 * @param y a vector
 * @return double inner product
 */
ll dot(const std::vector<ll> x, const std::vector<ll> y){
    ll z = 0;
    const int n = x.size();
    for(int i = 0; i < n; ++i) z += x[i] * y[i];
    return z;
}
double dot(const std::vector<ll> x, const std::vector<double> y){
    double z = 0.0;
    const int n = x.size();
    for(int i = 0; i < n; ++i) z += x[i] * y[i];
    return z;
}
double dot(const std::vector<double> x, const std::vector<double> y){
    double z = 0.0;
    const int n = x.size();
    for(int i = 0; i < n; ++i) z += x[i] * y[i];
    return z;
}


/**
 * @brief Computes transpose matrix of A
 * 
 * @param A matrix
 * @return std::vector<std::vector<long>> transposed matrix
 */
std::vector<std::vector<long>> trans(const std::vector<std::vector<long>> A){
    const int n = A[0].size(), m = A.size(); int i, j;
    std::vector<std::vector<long>> B(n, std::vector<long>(m));
    for(i = 0; i < n; ++i) for(j = 0; j < m; ++j) B[i][j] = A[j][i];
    return B;
}
std::vector<std::vector<ll>> trans(const std::vector<std::vector<ll>> A){
    const int n = A[0].size(), m = A.size(); int i, j;
    std::vector<std::vector<ll>> B(n, std::vector<ll>(m));
    for(i = 0; i < n; ++i) for(j = 0; j < m; ++j) B[i][j] = A[j][i];
    return B;
}


/**
 * @brief Computes product of matrices and vectors
 * 
 * @param x vector
 * @param A matrix
 * @return std::vector<double> x * A
 */
std::vector<double> mul_vec_mat(const std::vector<ll> x, const std::vector<std::vector<double>> A){
    const int m = A[0].size(), n = x.size(); int i, j;
    std::vector<double> v(m);
    for(i = 0; i < m; ++i) for(j = 0; j < n; ++j) v[i] += x[j] * A[j][i];
    return v;
}

/*Computes product of two matriceses*/
std::vector<std::vector<double>> mul_mat(const std::vector<std::vector<ll>> A, const std::vector<std::vector<ll>> B){
    const int n = A.size(), m = B[0].size(), l = B.size(); int i, j, k;
    std::vector<std::vector<double>> C(n, std::vector<double>(m));
    for(i = 0; i < n; ++i){
        for(j = 0; j < m; ++j){
            for(k = 0; k < l; ++k) C[i][j] += A[i][k] * B[k][j];
        }
    }
    return C;
}
std::vector<std::vector<double>> mul_mat(const std::vector<std::vector<ll>> A, const std::vector<std::vector<double>> B){
    const int n = A.size(), m = B[0].size(), l = B.size(); int i, j, k;
    std::vector<std::vector<double>> C(n, std::vector<double>(m));
    for(i = 0; i < n; ++i){
        for(j = 0; j < m; ++j){
            for(k = 0; k < l; ++k) C[i][j] += A[i][k] * B[k][j];
        }
    }
    return C;
}

/*Preparation of Gaussian elimination over Z/2Z*/
bool have_1(const std::vector<long> v, const long j){
    for(int i = 0; i < j; ++i) if(v[i] == 1) return true;
    return false;
}

/*Gaussian elimination of augumented matrices over Z/2Z*/
void aug_gauss_mod2(std::vector<std::vector<long>>& A, std::vector<std::vector<long>>& B){
    const int n = A.size(), m = A[0].size(); int i, j, k, s;
    bool flag;
    for(i = 0; i < n; ++i){
        for(j = 0; j < m; ++j){
            A[i][j] = std::abs(A[i][j] % 2); B[i][j] = std::abs(B[i][j] % 2);
        }
    }
    for(j = 0; j < m; ++j){
        flag = false;
        for(i = 0; i < n; ++i){
            if(A[i][j] == 1 && (! flag) && (! have_1(A[i], j))){
                flag = true; s = i;
            }
        }
        for(i = 0; i < n; ++i){
            if(A[i][j] == 1 && flag && i != s){
                for(k = 0; k < m; ++k){
                    A[i][k] ^= A[s][k]; B[i][k] ^= B[s][k];
                }
            }
        }
    }
}

/*Computes basises of kaernel of matrices A mod 2*/
std::vector<std::vector<long>> kernel_basis(const std::vector<std::vector<long>> A){
    const int n = A[0].size();
    std::vector<std::vector<long>> B = trans(A), E = identity_mat_long(n), ker;
    aug_gauss_mod2(B, E);
    
    for(int i = 0; i < n; ++i) if(is_zero(B[i])) ker.emplace_back(E[i]);
    return ker;
}

/*Computes all elements in kernel of A*/
std::vector<std::vector<long>> kernel_mod2(const std::vector<std::vector<long>> A){
    std::vector<std::vector<long>> v_all, ker;
    const std::vector<std::vector<long>> ker_basis = kernel_basis(A);
    
    /*Computes all basis vectors of a kernel of A*/
    const int d = ker_basis.size(), m = ker_basis[0].size(); int i, j, k;
    std::vector<long> tmp(m);
    v_all = vector_all(d);
    const int s = v_all.size();

    /*Computes all linear combinations of basis vector of kernel of A modulo 2*/
    for(i = 0; i < s; ++i){
        for(j = 0; j < d; ++j){
            for(k = 0; k < m; ++k) tmp[k] ^= v_all[i][j] * ker_basis[j][k];
        }
        ker.emplace_back(tmp);
        tmp = zero_vector_ZZ(m);
    }
    return ker;
}


/**
 * @brief Computes inverse matrix
 * 
 * @param a Matrix
 * @return std::vector<std::vector<double>> inverse matrix
 */
std::vector<std::vector<double>> InverseMatrix(std::vector<std::vector<double>> a){
    double buf;
    const int n = a.size();
    int i,j,k;
    std::vector<std::vector<double>> inv_a(n, std::vector<double>(n));

    for(i = 0; i < n; ++i) inv_a[i][i] = 1.0;

    for(i = 0; i < n; ++i){
        buf = 1.0 / a[i][i];
        for(j = 0; j < n; ++j){
            a[i][j] *= buf;
            inv_a[i][j] *= buf;
        }
        for(j = 0; j < n; ++j){
            if(i != j){
                buf = a[j][i];
                for(k = 0; k < n; ++k){
                    a[j][k] -= a[i][k] * buf;
                    inv_a[j][k] -= inv_a[i][k] * buf;
                }
            }
        }
    }
    return inv_a;
}

/**
 * @brief Generates factor basis
 * 
 * @param n 
 * @return std::vector<long> factor basis
 */
std::vector<long> factor_basis(const long n){
    NTL::PrimeSeq s;
    long p;
    std::vector<long> p_list = {-1};

    while(p_list.size() <= n) p_list.emplace_back(s.next());
    return p_list;
}


/**
 * @brief Gram_Schmidt's orthogonalization algorithm
 * 
 * @param b 
 * @return std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> 
 */
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> Gram_Schmidt(const std::vector<std::vector<ll>> b){
    const int n = b.size(), m = b[0].size(); int i, j, k;
    std::vector<std::vector<double>> GSOb(n, std::vector<double>(m)), mu = identity_mat(n);
    for(i = 0; i < n; ++i){
        for(j = 0; j < m; ++j) GSOb[i][j] = b[i][j];
        for(j = 0; j < i; ++j){
            mu[i][j] = dot(b[i], GSOb[j]) / dot(GSOb[j], GSOb[j]);
            for(k = 0; k < m; ++k) GSOb[i][k] -= mu[i][j] * GSOb[j][k];
        }
    }
    return std::forward_as_tuple(GSOb, mu);
}
std::tuple<std::vector<double>, std::vector<std::vector<double>>> Gram_Schmidt_squared(const std::vector<std::vector<ll>> b){
    const int n = b.size(), m = b[0].size(); int i, j, k;
    std::vector<double> B(n);
    std::vector<std::vector<double>> GSOb(n, std::vector<double>(m)), mu = identity_mat(n);
    for(i = 0; i < n; ++i){
        for(j = 0; j < m; ++j) GSOb[i][j] = b[i][j];
        for(j = 0; j < i; ++j){
            mu[i][j] = dot(b[i], GSOb[j]) / dot(GSOb[j], GSOb[j]);
            for(k = 0; k < m; ++k) GSOb[i][k] -= mu[i][j] * GSOb[j][k];
        }
        B[i] = dot(GSOb[i], GSOb[i]);
    }
    return std::forward_as_tuple(B, mu);
}


/**
 * @brief Babai's nearest plane algorithm.
 * 
 * @param b lattice basis matrix
 * @param w vector
 * @return std::vector<ll> a vector that is approx. solution of CVP for target w.
 */
double babai_error(const std::vector<std::vector<ll>> b, const std::vector<ll> w){
    std::vector<std::vector<double>> GSOb, mu;
    const int n = b.size(), m = w.size(); int i, j;
    std::vector<ll> t = w;
    double c;
    std::tie(GSOb, mu) = Gram_Schmidt(b);
    for(i = n - 1; i >= 0; --i){
        c = round(dot(t, GSOb[i]) / dot(GSOb[i], GSOb[i]));
        for(j = 0; j < m; ++j)t[j] -= c * GSOb[i][j];
    }
    return dot(t, t);
}


/**
 * @brief Conversion of coefficient and lattice.
 * 
 * @param v 
 * @param b 
 * @return std::vector<double> 
 */
std::vector<double> coef2lat(const std::vector<ll> v, const std::vector<std::vector<ll>> b){
    const int n = b.size(), m = b[0].size(); int i, j;
    std::vector<double> x(m);
    for(i = 0; i < n; ++i){
        for(j = 0; j < m; ++j) x[j] += v[i] * b[i][j];
    }
    return x;
}
std::vector<double> lat2coef(const std::vector<double> v, const std::vector<std::vector<ll>> b){
    const int n = b.size();
    std::vector<double> x(n);
    for(int i = 0; i < n; ++i) x[i] = v[i] / b[i][i];
    return x;
}


/**
 * @brief Find a close vector
 * 
 * @param mu GSO-coefficient
 * @param B squared norm of GSO-vector
 * @param R bound of enumerate
 * @param a coefficient of target vector
 * @return std::vector<ll> 
 */
std::vector<ll> ENUM_CVP(const std::vector<std::vector<double>> mu, const std::vector<double> B, double R, const std::vector<double> a){
    const int n = B.size(); int i, k;
    double tmp;
    std::vector<std::vector<double>> sigma(n + 1, std::vector<double>(n)), CVP_list;
    std::vector<long> r(n), w(n);
    std::vector<double> c(n), rho(n + 1), u;
    std::vector<ll> v(n), x(n);
    v[0] = 1;
    for(i = 0; i < n; ++i) r[i] = i;
    for(k = n - 1; k >= 0; --k){
        for(i = n - 1; i >= k + 1; --i) sigma[i][k] = sigma[i + 1][k] + (a[i] - v[i]) * mu[i][k];
        c[k] = a[k] + sigma[k + 1][k];
        v[k] = round(c[k]);
        w[k] = 1;
        tmp = c[k] - v[k]; tmp *= tmp;
        rho[k] = rho[k + 1] + tmp * B[k];
    }
    k = 0;
    for(k = 0; k <= n;){
        tmp = c[k] - v[k]; tmp *= tmp;
        rho[k] = rho[k + 1] + tmp * B[k];
        if(rho[k] <= R){
            if(k == 0) return v;
            --k;
            r[k] = (r[k] > r[k + 1] ? r[k]: r[k + 1]);
            for(i = r[k]; i > k; --i) sigma[i][k] = sigma[i + 1][k] + (a[i] - v[i]) * mu[i][k];
            c[k] = a[k] + sigma[k + 1][k];
            v[k] = round(c[k]);
            w[k] = 1;
        }else{
            ++k;
            if(k == n){x.clear(); return x;}
            r[k - 1] = k;
            if(v[k] > c[k]) v[k] -= w[k];else v[k] += w[k];
            ++w[k];
        }
    }
}


/**
 * @brief Enumerates the closest vector
 * 
 * @param mu GSO-coefficient
 * @param B squared norm of GSO-vector
 * @param R bound of enumerate
 * @param a coefficient of target vector
 * @param t target vector
 * @param b lattice
 * @return std::vector<std::vector<ll>> the closest vector for target a
 */
std::vector<std::vector<ll>> ENUM_CVP_all(const std::vector<std::vector<double>> mu, const std::vector<double> B, double R, const std::vector<double> a, const std::vector<ll> t, const std::vector<std::vector<ll>> b){
    const int n = B.size(), m = t.size(); int i;
    std::vector<std::vector<ll>> CVP_list;
    std::vector<ll> ENUM_CVP_v(n), pre_ENUM_CVP_v(n), x(n);
    for(;;){
        pre_ENUM_CVP_v = ENUM_CVP_v;
        ENUM_CVP_v = ENUM_CVP(mu, B, R, a);
        if(ENUM_CVP_v.empty()) return CVP_list;
        if(std::find(CVP_list.begin(), CVP_list.end(), ENUM_CVP_v) == CVP_list.end()) CVP_list.emplace_back(ENUM_CVP_v);

        R *= 0.9;
    }
}


/**
 * @brief Generates a prime matrix
 * 
 * @param c precision parameter
 * @param p primes
 * @param n rank of lattices
 * @return std::tuple<std::vector<std::vector<ll>>, NTL::mat_ZZ> a prime matrix
 */
std::tuple<std::vector<std::vector<ll>>, NTL::mat_ZZ> imp_prime_mat(const double c, const std::vector<long> p, const long n){
    std::vector<std::vector<ll>> A(n, std::vector<ll>(n + 1));
    NTL::mat_ZZ B; B.SetDims(n, n + 1);
    std::vector<ll> diag(n);//entry
    int i;
    for(i = 0; i < n; ++i) diag[i] = ceil((i + 1.0) * 0.6);
    std::random_device seed_gen;
    std::mt19937_64 engine(seed_gen());
    std::shuffle(diag.begin(), diag.end(), engine); //shuffle randomly

    if(c == 5.0){
        for(i = 0; i < n; ++i){
            A[i][i] = diag[i];
            A[i][n] = round(100000 * log(p[i + 1]));
            B[i][i] = diag[i];
            B[i][n] = A[i][n];
        }
    }else{
        for(i = 0; i < n; ++i){
            A[i][i] = diag[i];
            A[i][n] = round(pow(10, c) * log(p[i + 1]));
            B[i][i] = diag[i];
            B[i][n] = A[i][n];
        }
    }

    return std::forward_as_tuple(A, B);
}


/**
 * @brief Generates a target vector
 * 
 * @param N 
 * @param c 
 * @param n 
 * @return std::vector<ll> 
 */
std::vector<ll> target(const NTL::ZZ N, const double c, const long n){
    std::vector<ll> t(n + 1);
    if(c == 5.0) t[n] = round(100000 * NTL::to_double(NTL::log(N)));
    else t[n] = round(pow(10, c) * NTL::to_double(NTL::log(N)));
    return t;
}


/**
 * @brief Test if (u, v) is sr-pair or not.
 * 
 * @param u smooth number
 * @param v smooth number
 * @param N a composite number
 * @param p primes
 * @param n rank of lattice
 * @return true if (u, v) is sr-pair
 * @return false if (u, v) is nor sr-pair
 */
bool sr_test(const NTL::ZZ u, const NTL::ZZ v, const NTL::ZZ N, const std::vector<long> p, const long n){
    NTL::ZZ a = u - N * v, r;
    for(int i = 0; i < n; ++i){
        r = p[i + 1];
        while(a % r == 0 && a != 0) a /= r;
    }
    if(a == 1 || a == -1) return true;
    return false;
}

extern "C" ll *lattice_factorization(const int bit_flag, const double NN, const int info_flag){
     NTL::RR::SetPrecision(60);

    bool new_pair = false;
    long loop_times = 0, num = 0, e1, e2, s, n;
    int i, j, k;
    float c = 5.0;
    ll* rtn = (ll *)malloc(2 * sizeof(ll));
    NTL::ZZ q, a, u, v, U, S = NTL::to_ZZ(0), N = NTL::to_ZZ(NN);
    NTL::ZZ_p X, Y;
    std::vector<double> B, w, e;
    std::vector<std::vector<double>> mu;
    std::vector<std::vector<ll>> L, close_vecs;
    NTL::mat_ZZ LL;

    if(bit_flag) GenSemiPrime(N);
    //if(argc >= 4){n = atoi(argv[3]); if(argc >= 5) c = atof(argv[4]); }else
    n = 2.2 * NTL::NumBits(N) / log(NTL::NumBits(N)) - 8;
    const long K = 2 * n * n, J = K * 0.66, m = n + 1;

    X.init(N); Y.init(N);

    std::vector<long> vec(K), ee(K);
    const std::vector<long> p = factor_basis(K);
    std::vector<std::vector<long>> kernel, prekernel, A(K, std::vector<long>(K)), already_ee = {};
    const std::vector<ll> t = target(N, c, n);
    std::vector<std::vector<ll>> BB(n, std::vector<ll>(m));

    for(;;){
        ++loop_times;
        w = zero_vector_RR(n);

        /*Composition of CVP*/
        std::tie(L, LL) = imp_prime_mat(c, p, n);

        /*Reduces lattice basis*/
        NTL::LLL(u, LL);
        for(i = 0; i < n; ++i) for(j = 0; j < m; ++j) BB[i][j] = NTL::to_double(LL[i][j]);

        std::tie(B, mu) =    Gram_Schmidt_squared(BB);
        
        /*Composition of coefficient vectors of target vectors*/
        w = mul_vec_mat(t, mul_mat(trans(BB), InverseMatrix(mul_mat(BB, trans(BB)))));

        /*Computes approximate solutions of CVP*/
        close_vecs = ENUM_CVP_all(mu, B, babai_error(BB, t), w, t, BB);
        S += close_vecs.size();
        for(i = 0; i < close_vecs.size(); ++i){
            e = lat2coef(coef2lat(close_vecs[i], BB), L);

            /*Composition of a pair of integers (u, v) correspond to close_vec[i]*/
            u = v = 1;
            for(j = 0; j < n; ++j){
                if(e[j] > 0) u *= NTL::to_ZZ(NTL::power(NTL::to_ZZ(p[j + 1]), e[j]));
                else if(e[j] < 0) v *= NTL::to_ZZ(NTL::power(NTL::to_ZZ(p[j + 1]), -e[j]));
            }

            if(u > 0 && v > 0 && sr_test(u, v, N, p, K)){
                vec = zero_vector_ZZ(K);
                for(j = 0; j < K; ++j){
                    e1 = e2 = 0;
                    U = u;
                    q = p[j + 1];
                    
                    if(j < n) while(U % q == 0 && U != 0){++e1; U /= q;}
                    
                    a = u - v * N;
                    while(a % q == 0 && a != 0){++e2; a /= q;}
                    vec[j] = e2 - e1;
                }

                if(std::find(A.begin(), A.end(), vec) == A.end()){
                    if(num < K) A[num] = vec;else A.emplace_back(vec);
                    ++num; new_pair = true;
                    if(info_flag) printf("%ld pairs were found. (%ld pairs are needed.)\n",num, K);
                }
            }
        }

        if(num >= J){
            if(new_pair){
                ee = zero_vector_ZZ(K);
                kernel = kernel_mod2(trans(A)); s = kernel.size();
                for(i = 0; i < prekernel.size(); ++i)prekernel[i].resize(kernel[0].size());

                if(kernel != prekernel){
                    for(i = 0; i < s; ++i){
                        if(std::find(prekernel.begin(), prekernel.end(), kernel[i]) == prekernel.end()){
                            for(j = 0; j < K; ++j) if(kernel[i][j] == 1) for(k = 0; k < K; ++k) ee[k] += A[j][k];

                            if(std::find(already_ee.begin(), already_ee.end(), ee) == already_ee.end()){
                                already_ee.emplace_back(ee);
                                X = Y = 1;
                                for(j = 0; j < K; ++j){
                                    if(ee[j] > 0) X *= NTL::conv<NTL::ZZ_p>(NTL::PowerMod(NTL::to_ZZ(p[j + 1]), NTL::to_ZZ(ee[j] / 2), N));
                                    else if(ee[j] < 0) Y *= NTL::conv<NTL::ZZ_p>(NTL::PowerMod(NTL::to_ZZ(p[j + 1]), NTL::to_ZZ(-ee[j] / 2), N));
                                }

                                if(X != 1 || Y != 1){
                                    if(info_flag) std::cout << "X = " << X << ", Y = " << Y << std::endl;
                                    if(X != Y) q = NTL::GCD(NTL::rep(X - Y), N);else q = NTL::GCD(NTL::rep(X + Y), N);

                                    if(q != 1 && q != N){
                                        if(info_flag){
                                            printf("==============================\nN = "); std::cout << N;
                                            printf(", bit size = %ld\nc = %.1f, beta = 2.0\nN = %lu * %lu\nloop times = %ld\nnumber of sr-pairs = %ld\n==============================", NTL::NumBits(N), c, NTL::to_long(q), NTL::to_long(N / q), loop_times, num);
                                            std::cout << "\n#Vector = " << S << std::endl;
                                        }
                                        rtn[0] = NTL::to_ulong(q);
                                        rtn[1] = NTL::to_ulong(N / q);
                                        return rtn;
                                    }
                                }
                            }
                            ee = zero_vector_ZZ(K);
                        }
                    }
                    prekernel = kernel;
                }
                new_pair = false;
            }
        }
    }
}
