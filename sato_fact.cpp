#include<iostream>
#include<vector>
#include<cmath>
#include<random>
#include<iomanip>
#include<tuple>
#include<algorithm>
#include<cstdlib>
#include<string>
#include<cfloat>
#include<NTL/ZZ.h>
#include<NTL/ZZ_p.h>
#include<NTL/RR.h>
#include<NTL/mat_ZZ.h>
#include<NTL/matrix.h>
#include<NTL/vector.h>
#include<NTL/LLL.h>
#pragma GCC target("avx2")
#define EPSILON 0.0000000000000000000000000000001

bool vector_not_equal(const std::vector<double> x, const std::vector<double> y){
  const int n = x.size();
  for(int i = 0; i < n; ++i) if(fabs(x.at(i) - y.at(i)) > DBL_MIN) return true;
  return false;
}

bool not_included(const std::vector<double> x, const std::vector<std::vector<double>> vecs){
  const int n = vecs.size();
  if(! n) return true;
  for(int i = 0; i < n; ++i) if(vector_not_equal(x, vecs.at(i))) return true;
  return false;
}

std::vector<double> zero_vector_RR(const long n){const std::vector<double> z(n); return z;}
std::vector<long> zero_vector_ZZ(const long n){const std::vector<long> z(n); return z;}

/*Generates all combinations of n 0 and 1*/
std::vector<std::vector<long>> vector_all(const long n){
  std::vector<long> tmp(n), v(n);
  std::vector<std::vector<long>> v_all;
  for(int i = 0; i < n; ++i){
    tmp.at(i) = 1; v = tmp;
    do v_all.push_back(v); while(std::next_permutation(v.begin(), v.end()));
  }
  return v_all;
}

void GenSemiPrime(NTL::ZZ& n){
  const int m = NTL::to_int(n) * 0.5, k = NTL::to_int(n);
  while(true){
    n = NTL::RandomPrime_ZZ(m) * NTL::RandomPrime_ZZ(NTL::to_int(k) - m + 1);
    if(NTL::NumBits(n) == k) break;
  }
}

/*Tests v is zero-vector or not*/
bool is_zero(const std::vector<long> v){
  const int n = v.size();
  for(int i = 0; i < n; ++i) if(v.at(i) != 0) return false;
  return true;
}

/*Generates an identity matrices*/
std::vector<std::vector<double>> identity_mat(const long n){
  std::vector<std::vector<double>> A(n, std::vector<double>(n));
  for(int i = 0; i < n; ++i) A.at(i).at(i) = 1.0;
  return A;
}

/*Generates an identity matrices*/
std::vector<std::vector<long>> identity_mat_long(const long n){
  std::vector<std::vector<long>> A(n, std::vector<long>(n));
  for(int i = 0; i < n; ++i) A.at(i).at(i) = 1;
  return A;
}

#if 0
/*Print a matrices A*/
void print_mat(const std::vector<std::vector<auto>> A){
  puts("[");
  for(long i = 0; i < A.size(); ++i){
    printf("[");
    for(long j = 0; j < A.at(0).size(); ++j) std::cout << A.at(i).at(j) << " ";
    puts("]");
  }
  puts("]");
}

/*Prints a vector v*/
void print_vec(const std::vector<auto> v){
  printf("[");
  for(long i = 0; i < v.size(); ++i) std::cout << v.at(i) << " ";
  puts("]");
}
#endif

/*Computes inner product of vectors x and y*/
double dot(const std::vector<double> x, const std::vector<double> y){
  double z = 0.0;
  const int n = x.size();
  for(int i = 0; i < n; ++i) z += x.at(i) * y.at(i);
  return z;
}

/*Computes transpose matrices of A*/
std::vector<std::vector<long>> trans(const std::vector<std::vector<long>> A){
  const int n = A.at(0).size(), m = A.size(); int i, j;
  std::vector<std::vector<long>> B(n, std::vector<long>(m));
  for(i = 0; i < n; ++i) for(j = 0; j < m; ++j) B.at(i).at(j) = A.at(j).at(i);
  return B;
}

/*Computes transpose matrices of A*/
std::vector<std::vector<double>> trans_RR(const std::vector<std::vector<double>> A){
  const int n = A.at(0).size(), m = A.size(); int i, j;
  std::vector<std::vector<double>> B(n, std::vector<double>(m));
  for(i = 0; i < n; ++i) for(j = 0; j < m; ++j) B.at(i).at(j) = A.at(j).at(i);
  return B;
}

/*Computes product of matrices and vectors*/
std::vector<double> mul_vec_mat(const std::vector<double> x, const std::vector<std::vector<double>> A){
  const int m = A.at(0).size(), n = x.size(); int i, j;
  std::vector<double> v(m);
  for(i = 0; i < m; ++i) for(j = 0; j < n; ++j) v.at(i) += x.at(j) * A.at(j).at(i);
  return v;
}

/*Computes product of two matriceses*/
std::vector<std::vector<double>> mul_mat(const std::vector<std::vector<double>> A, const std::vector<std::vector<double>> B){
  const int n = A.size(), m = B.at(0).size(), l = B.size(); int i, j, k;
  std::vector<std::vector<double>> C(n, std::vector<double>(m));
  for(i = 0; i < n; ++i){
    for(j = 0; j < m; ++j){
      for(k = 0; k < l; ++k) C.at(i).at(j) += A.at(i).at(k) * B.at(k).at(j);
    }
  }
  return C;
}

/*Preparation of Gaussian elimination over Z/2Z*/
bool have_1(const std::vector<long> v, const long j){
  for(int i = 0; i < j; ++i) if(v.at(i) == 1) return true;
  return false;
}

/*Gaussian elimination of augumented matrices over Z/2Z*/
void aug_gauss_mod2(std::vector<std::vector<long>>& A, std::vector<std::vector<long>>& B){
  const int n = A.size(), m = A.at(0).size(); int i, j, k, s;
  bool flag;
  for(i = 0; i < n; ++i){
    for(j = 0; j < m; ++j){
      A.at(i).at(j) = fabs(A.at(i).at(j) % 2); B.at(i).at(j) = fabs(B.at(i).at(j) % 2);
    }
  }
  for(j = 0; j < m; ++j){
    flag = false;
    for(i = 0; i < n; ++i){
      if(A.at(i).at(j) == 1 && (! flag) && (! have_1(A.at(i), j))){
        flag = true; s = i;
      }
    }
    for(i = 0; i < n; ++i){
      if(A.at(i).at(j) == 1 && flag && i != s){
        for(k = 0; k < m; ++k){
          A.at(i).at(k) ^= A.at(s).at(k); B.at(i).at(k) ^= B.at(s).at(k);
        }
      }
    }
  }
}

/*Computes basises of kaernel of matrices A mod 2*/
std::vector<std::vector<long>> kernel_basis(const std::vector<std::vector<long>> A){
  const int n = A.at(0).size();
  std::vector<std::vector<long>> B = trans(A), E = identity_mat_long(n), ker;
  aug_gauss_mod2(B, E);
  
  for(int i = 0; i < n; ++i) if(is_zero(B.at(i))) ker.push_back(E.at(i));
  return ker;
}

/*Computes all elements in kernel of A*/
std::vector<std::vector<long>> kernel_mod2(const std::vector<std::vector<long>> A){
  std::vector<std::vector<long>> v_all, ker;
  const std::vector<std::vector<long>> ker_basis = kernel_basis(A);
  
  /*Computes all basis vectors of a kernel of A*/
  const int d = ker_basis.size(), m = ker_basis.at(0).size(); int i, j, k;
  std::vector<long> tmp(m);
  v_all = vector_all(d);
  const int s = v_all.size();

  /*Computes all linear combinations of basis vector of kernel of A modulo 2*/
  for(i = 0; i < s; ++i){
    for(j = 0; j < d; ++j){
      for(k = 0; k < m; ++k) tmp.at(k) ^= v_all.at(i).at(j) * ker_basis.at(j).at(k);
    }
    ker.push_back(tmp);
    tmp = zero_vector_ZZ(m);
  }
  return ker;
}

/*Computes inverse matrices of M*/
std::vector<std::vector<double>> mat_inv(std::vector<std::vector<double>> M){
  std::vector<std::vector<double>> X;
  const int n = M.size(), m = M.at(0).size(); int i, j;
  NTL::mat_ZZ A, A_inv;
  NTL::ZZ d;
  A.SetDims(n, m);
  for(i = 0; i < n; ++i) for(j = 0; j < m; ++j) A[i][j] = NTL::to_ZZ(M.at(i).at(j));
  NTL::inv(d, A_inv, A);
  const double e = 1.0 / NTL::to_double(d);
  for(i = 0; i < n; ++i) for(j = 0; j < m; ++j) M.at(i).at(j) = NTL::to_double(A_inv[i][j]) * e;
  return M;
}

std::vector<long> factor_basis(const long n){
  NTL::PrimeSeq s;
  long p;
  std::vector<long> p_list = {-1};

  while(p_list.size() <= n) p_list.push_back(s.next());
  return p_list;
}

/*Gram_Schmidt's orthogonalization algorithm*/
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> Gram_Schmidt(const std::vector<std::vector<double>> b){
  const int n = b.size(), m = b.at(0).size(); int i, j, k;
  std::vector<std::vector<double>> GSOb(n, std::vector<double>(m)), mu = identity_mat(n);
  for(i = 0; i < n; ++i){
    for(j = 0; j < m; ++j) GSOb.at(i).at(j) = b.at(i).at(j);
    for(j = 0; j < i; ++j){
      mu.at(i).at(j) = dot(b.at(i), GSOb.at(j)) / dot(GSOb.at(j), GSOb.at(j));
      for(k = 0; k < m; ++k) GSOb.at(i).at(k) -= mu.at(i).at(j) * GSOb.at(j).at(k);
    }
  }
  return std::forward_as_tuple(GSOb, mu);
}

/*Gram_Schmidt's orthogonalization algorithm*/
std::tuple<std::vector<double>, std::vector<std::vector<double>>> Gram_Schmidt_squared(const std::vector<std::vector<double>> b){
  const int n = b.size(), m = b.at(0).size(); int i, j, k;
  std::vector<double> B(n);
  std::vector<std::vector<double>> GSOb(n, std::vector<double>(m)), mu = identity_mat(n);
  for(i = 0; i < n; ++i){
    for(j = 0; j < m; ++j) GSOb.at(i).at(j) = b.at(i).at(j);
    for(j = 0; j < i; ++j){
      mu.at(i).at(j) = dot(b.at(i), GSOb.at(j)) / dot(GSOb.at(j), GSOb.at(j));
      for(k = 0; k < m; ++k) GSOb.at(i).at(k) -= mu.at(i).at(j) * GSOb.at(j).at(k);
    }
    B.at(i) = dot(GSOb.at(i), GSOb.at(i));
  }
  return std::forward_as_tuple(B, mu);
}

std::vector<std::vector<double>> LLL_reduce(const std::vector<std::vector<double>> b, const float delta = 0.99){
  const int n = b.size(), m = b.at(0).size(); int i, j;
  NTL::ZZ det;
  std::vector<std::vector<double>> B(n, std::vector<double>(m));
  NTL::Mat<NTL::ZZ> c; c.SetDims(n, m);
  for(i = 0; i < n; ++i) for(j = 0; j < m; ++j) c[i][j] = NTL::to_ZZ(b.at(i).at(j));
  NTL::LLL(det, c, delta);
  for(i = 0; i < n; ++i) for(j = 0; j < m; ++j) B.at(i).at(j) = NTL::to_double(c[i][j]);
  return B;
}

std::vector<double> Babai(const std::vector<std::vector<double>> b, const std::vector<double> w){
  std::vector<std::vector<double>> GSOb, mu;
  const int n = b.size(), m = w.size(); int i, j;
  std::vector<double> t = w, v(m);
  double c;
  std::tie(GSOb, mu) = Gram_Schmidt(b);
  for(i = n - 1; i >= 0; --i){
    c = round(dot(t, GSOb.at(i)) / dot(GSOb.at(i), GSOb.at(i)));
    for(j = 0; j < m; ++j)t.at(j) -= c * GSOb.at(i).at(j);
  }
  for(i = 0; i < m; ++i)v.at(i) = w.at(i) - t.at(i);
  return v;
}

double babai_error(const std::vector<std::vector<double>> b, const std::vector<double> w){
  std::vector<double> t = Babai(b, w), error_vec(w.size());
  const int n = error_vec.size();
  for(int i = 0; i < n; ++i) error_vec.at(i) = t.at(i) - w.at(i);
  return sqrt(dot(error_vec, error_vec));
}

/*Converts coefficient vector to lattice vector*/
std::vector<double> coef2lat(const std::vector<auto> v, const std::vector<std::vector<auto>> b){
  const int n = b.size(), m = b.at(0).size(); int i, j;
  std::vector<double> x(m);
  for(i = 0; i < n; ++i){
    for(j = 0; j < m; ++j) x.at(j) += v.at(i) * b.at(i).at(j);
  }
  return x;
}

std::vector<double> coef2ldlat(const std::vector<auto> v, const std::vector<std::vector<auto>> b){
  const int n = b.size(), m = b.at(0).size(); int i, j;
  std::vector<double> x(m);
  for(i = 0; i < n; ++i){
    for(j = 0; j < m; ++j) x.at(j) += v.at(i) * b.at(i).at(j);
  }
  return x;
}
std::vector<double> lat2coef(const std::vector<auto> v, const std::vector<std::vector<auto>> b){
  const int n = b.size();
  std::vector<double> x(n);
  for(int i = 0; i < n; ++i) x.at(i) = v.at(i) / b.at(i).at(i);
  return x;
}

/*Enumerates closest vectors*/
std::vector<double> ENUM_CVP(const std::vector<std::vector<double>> mu, const std::vector<double> B, const double R, const std::vector<double> a){
  const int n = B.size(); int i, k;
  std::vector<std::vector<double>> sigma(n + 2, std::vector<double>(n + 1)), CVP_list;
  std::vector<long> r(n + 1), w(n + 1);
  std::vector<double> c(n + 1), rho(n + 2), v(n + 1), x(n), u;
  double R2 = R * R;
  v.at(1) = 1;
  for(i = 0; i <= n; ++i) r.at(i) = i;
  for(k = n; k >= 1; --k){
    for(i = n; i >= k + 1; --i) sigma.at(i).at(k) = sigma.at(i + 1).at(k) + (a.at(i - 1) - v.at(i)) * mu.at(i - 1).at(k - 1);
    c.at(k) = a.at(k - 1) + sigma.at(k + 1).at(k);
    v.at(k) = round(c.at(k));
    w.at(k) = 1;
    rho.at(k) = rho.at(k + 1) + (c.at(k) - v.at(k)) * (c.at(k) - v.at(k)) * B.at(k - 1);
  }
  k = 1;
  while(true){
    rho.at(k) = rho.at(k + 1) + (c.at(k) - v.at(k)) * (c.at(k) - v.at(k)) * B.at(k - 1);
    if(rho.at(k) < R2 || fabs(rho.at(k) - R2) < DBL_MIN){
      if(k == 1){
        for(i = 1; i <= n; ++i) x.at(i - 1) = v.at(i);
        return x;
      }
      --k;
      r.at(k - 1) = fmax(r.at(k - 1), r.at(k));
      for(i = r.at(k); i > k; --i) sigma.at(i).at(k) = sigma.at(i + 1).at(k) + (a.at(i - 1) - v.at(i)) * mu.at(i - 1).at(k - 1);
      c.at(k) = a.at(k - 1) + sigma.at(k + 1).at(k);
      v.at(k) = round(c.at(k));
      w.at(k) = 1;
    }else{
      ++k;
      if(k == n + 1){x.clear(); return x;}
      r.at(k - 1) = k;
      if(v.at(k) > c.at(k)) v.at(k) -= w.at(k);else v.at(k) += w.at(k);
      ++w.at(k);
    }
  }
}

/*Generates coefficient vectors of close vectors*/
std::vector<std::vector<double>> ENUM_CVP_all(const std::vector<std::vector<double>> mu, const std::vector<double> B, double R, const std::vector<double> a, const std::vector<double> t, const std::vector<std::vector<double>> b){
  const int n = B.size(), m = t.size(); int i;
  std::vector<std::vector<double>> CVP_list;
  std::vector<double> ENUM_CVP_v(n), pre_ENUM_CVP_v(n), x(n), c, d(m);
  while(true){
    for(i = 0; i < n; ++i) pre_ENUM_CVP_v.at(i) = ENUM_CVP_v.at(i);
    ENUM_CVP_v = ENUM_CVP(mu, B, R, a);
    if(ENUM_CVP_v.empty()) return CVP_list;
    if(not_included(ENUM_CVP_v, CVP_list))CVP_list.push_back(ENUM_CVP_v);

    c = coef2lat(ENUM_CVP_v, b);
    for(i = 0; i < m; ++i) d.at(i) = c.at(i) - t.at(i);
    R = fmin(sqrt(dot(d, d)), R * 0.9);
  }
}


/******************************************
main sub-routine starts
*******************************************/
std::vector<std::vector<double>> imp_prime_mat(const double c, const std::vector<long> p, const long n){
  std::vector<std::vector<double>> A(n, std::vector<double>(n + 1));
  std::vector<double> diag(n);//entry
  int i;
  for(i = 0; i < n; ++i) diag.at(i) = ceil((i + 1.0) * 0.5);

  std::random_device seed_gen;
  std::mt19937_64 engine(seed_gen());
  std::shuffle(diag.begin(), diag.end(), engine); //shuffle randomly

  if(c == 5.0){
    for(i = 0; i < n; ++i){
      A.at(i).at(i) = diag.at(i);
      A.at(i).at(n) = round(100000 * log(p.at(i + 1)));
    }
  }else{
    for(i = 0; i < n; ++i){
      A.at(i).at(i) = diag.at(i);
      A.at(i).at(n) = round(pow(10, c) * log(p.at(i + 1)));
    }
  }

  return A;
}

/*Generates a target vector*/
std::vector<double> target(const NTL::ZZ N, const double c, const long n){
  std::vector<double> t(n + 1);
  if(c == 5.0) t.at(n) = round(100000 * NTL::to_double(NTL::log(N)));
  else t.at(n) = round(pow(10, c) * NTL::to_double(NTL::log(N)));
  return t;
}

bool sr_test(const NTL::ZZ u, const NTL::ZZ v, const NTL::ZZ N, const std::vector<long> p, const long n){
  NTL::ZZ a = u - N * v, r;
  for(int i = 0; i < n; ++i){
    r = NTL::to_ZZ((long)p.at(i + 1));
    while(a % r == 0 && a != 0) a /= r;
  }
  if(a == 1 || a == -1) return true;
  return false;
}

/******************************************
main sub-routine ends
*******************************************/


/*factorization*/
int main(int argc, char **argv){
  bool new_pair = false; const bool bit_input = atoi(argv[1]);
  long loop_times = 0, num = 0, e1, e2, s, n;
  int i, j, k;
  float c = 5.0;
  NTL::ZZ q, a, u, v, U, N = NTL::to_ZZ(argv[2]), S = NTL::to_ZZ(0), XX, YY; XX = YY = 1;
  NTL::ZZ_p X, Y;
  std::vector<double> B, w, e;
  std::vector<std::vector<double>> BB, mu, close_vecs, L;

  if(bit_input) GenSemiPrime(N);
  if(argc >= 4){n = atoi(argv[3]); if(argc >= 5) c = atof(argv[4]); }else n = 2.2 * NTL::NumBits(N) / log(NTL::NumBits(N)) - 8;
  const long K = 2 * n * n, J = K * 0.66;

  X.init(N); Y.init(N);

  std::vector<long> vec(K), ee(K);
  const std::vector<long> p = factor_basis(K);
  std::vector<std::vector<long>> kernel, prekernel, A(K, std::vector<long>(K));
  const std::vector<double> t = target(N, c, n);

  while(true){
    ++loop_times;
    w = zero_vector_RR(n);

    /*Composition of CVP*/
    L = imp_prime_mat(c, p, n);

    /*Reduces lattice basis*/
    BB = LLL_reduce(L);

    std::tie(B, mu) =  Gram_Schmidt_squared(BB);
    
    /*Composition of coefficient vectors of target vectors*/
    w = mul_vec_mat(t, mul_mat(trans_RR(BB), mat_inv(mul_mat(BB, trans_RR(BB)))));

    /*Computes approximate solutions of CVP*/
    close_vecs = ENUM_CVP_all(mu, B, babai_error(BB, t), w, t, BB);
    S += close_vecs.size();
    for(i = 0; i < close_vecs.size(); ++i){
      e = lat2coef(coef2ldlat(close_vecs.at(i), BB), L);

      /*Composition of a pair of integers (u, v) correspond to close_vec[i]*/
      u = v = 1;
      for(j = 0; j < n; ++j){
        if(e.at(j) > 0) u *= NTL::to_ZZ(NTL::power(NTL::to_ZZ(p.at(j + 1)), e.at(j)));
        else if(e.at(j) < 0) v *= NTL::to_ZZ(NTL::power(NTL::to_ZZ(p.at(j + 1)), -e.at(j)));
      }

      if(u > 0 && v > 0 && sr_test(u, v, N, p, K)){
        vec = zero_vector_ZZ(K);
        for(j = 0; j < K; ++j){
          e1 = e2 = 0;
          U = u;
          q = NTL::to_ZZ(p.at(j + 1));
          
          if(j < n) while(U % q == 0 && U != 0){++e1; U /= q;}
          
          a = u - v * N;
          while(a % q == 0 && a != 0){++e2; a /= q;}
          vec.at(j) = e2 - e1;
        }

        if(std::find(A.begin(), A.end(), vec) == A.end()){
          if(num < K) A.at(num) = vec;else A.push_back(vec);
          ++num; new_pair = true;
          printf("%ld pairs were found. (%ld pairs are needed.)\n",num, K);
        }
      }
    }

    if(num >= J){
      if(new_pair){
        ee = zero_vector_ZZ(K);
        kernel = kernel_mod2(trans(A)); s = kernel.size();
        for(i = 0; i < prekernel.size(); ++i)prekernel.at(i).resize(kernel.at(0).size());

        if(kernel != prekernel){
          for(i = 0; i < s; ++i){
            if(std::find(prekernel.begin(), prekernel.end(), kernel.at(i)) == prekernel.end()){
              for(j = 0; j < K; ++j) if(kernel.at(i).at(j) == 1) for(k = 0; k < K; ++k) ee.at(k) += A.at(j).at(k);

              X = Y = 1;
              for(j = 0; j < K; ++j){
                if(ee.at(j) > 0) X *= NTL::conv<NTL::ZZ_p>(NTL::PowerMod(NTL::to_ZZ(p.at(j + 1)), NTL::to_ZZ(ee.at(j) / 2), N));
                else if(ee.at(j) < 0) Y *= NTL::conv<NTL::ZZ_p>(NTL::PowerMod(NTL::to_ZZ(p.at(j + 1)), NTL::to_ZZ(-ee.at(j) / 2), N));
              }

              if((X != 1 || Y != 1) && (XX != NTL::rep(X) && YY != NTL::rep(Y))){
                std::cout << "X = " << X << ", Y = " << Y << std::endl;
                if(X != Y) q = NTL::GCD(NTL::rep(X - Y), N);else q = NTL::GCD(NTL::rep(X + Y), N);

                if(q != 1 && q != N){
                  printf("==============================\nN = "); std::cout << N;
                  printf(", bit size = %ld\nc = %.1f, beta = 2.0\nN = %lu * %lu\nloop times = %ld\nnumber of sr-pairs = %ld\n==============================", NTL::NumBits(N), c, NTL::to_long(q), NTL::to_long(N / q), loop_times, num);
                  std::cout << "\n#Vector = " << S << std::endl;
                  return 0;
                }
                XX = NTL::rep(X); YY = NTL::rep(Y);
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
