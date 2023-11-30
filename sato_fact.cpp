#include<iostream>
#include <vector>
#include <cmath>
#include<random>
#include <iomanip>
#include<tuple>
#include<algorithm>
#include<cstdlib>
#include<string>
#include<NTL/ZZ.h>
#include<NTL/mat_ZZ.h>
#include<NTL/matrix.h>
#include<NTL/vector.h>
#include<NTL/LLL.h>
#define EPSILON 1e-6
#define BKZ_REDUCE 1

using namespace std;
using namespace NTL;

bool vector_is_equal(const vector<double> x, const vector<double> y){
  for(long i = 0; i < x.size(); ++i){
    if(fabs(x.at(i) - y.at(i)) > EPSILON)return false;
  }
  return true;
}

bool is_included(const vector<double> x, const vector<vector<double>> vecs){
  if(vecs.size() == 0)return false;
  for(long i = 0; i < vecs.size(); ++i){
    if(not vector_is_equal(x, vecs.at(i)))return false;
  }
  return true;
}

vector<double> zero_vector_RR(const int n){
  vector<double> z(n);
  return z;
}

vector<long long> zero_vector_ZZ(const int n){
  vector<long long> z(n);
  return z;
}

/*Generates all combinations of n 0 and 1*/
vector<vector<long long>> vector_all(const int n){
  vector<long long> tmp(n), v(n);
  vector<vector<long long>> v_all;
  for(long i = 0; i < n; ++i){
    tmp.at(i) = 1; v = tmp;
    do v_all.push_back(v); while(next_permutation(v.begin(), v.end()));
  }
  return v_all;
}

/*Tests v is zero-vector or not*/
bool is_zero(const vector<long long> v){
  for(long i = 0; i < v.size(); ++i) if(v.at(i) != 0) return false;
  return true;
}

/*Generates an identity matrices*/
vector<vector<double>> identity_mat(const int n){
  vector<vector<double>> A(n, vector<double>(n));
  for(long i = 0; i < n; ++i) A.at(i).at(i) = 1.0;
  return A;
}

/*Generates an identity matrices*/
vector<vector<long long>> identity_mat_int(const int n){
  vector<vector<long long>> A(n, vector<long long>(n));
  for(long i = 0; i < n; ++i) A.at(i).at(i) = 1;
  return A;
}

/*Prints a matrices A*/
void print_mat(const vector<vector<auto>> A){
  cout << "[\n";
  for(long i = 0; i < A.size(); ++i){
    cout << "[";
    for(long j = 0; j < A.at(0).size(); ++j) cout << A.at(i).at(j) << " ";
    cout << "]" << endl;
  }
  cout << "]\n";
}

/*Prints a vector v*/
void print_vec(const vector<auto> v){
  cout << "[";
  for(long i = 0; i < v.size(); ++i) cout << v.at(i) << " ";
  cout << "]\n";
}

/*Computes inner product of vectors x and y*/
float dot(const vector<auto> x, const vector<auto> y)
{
  float z = 0.0;
  for(long i = 0; i < x.size(); ++i) z += x.at(i) * y.at(i);
  return z;
}

/*Computes transpose matrices of A*/
vector<vector<long long>> trans(const vector<vector<auto>> A)
{
  int n = A.at(0).size(), m = A.size();
  vector<vector<long long>> B(n, vector<long long>(m));
  for(long i = 0; i < n; ++i){
    for(long j = 0; j < m; ++j) B.at(i).at(j) = A.at(j).at(i);
  }
  return B;
}

/*Computes transpose matrices of A*/
vector<vector<double>> trans_RR(const vector<vector<double>> A)
{
  int n = A.at(0).size(), m = A.size();
  vector<vector<double>> B(n, vector<double>(m));
  for(long i = 0; i < n; ++i){
    for(long j = 0; j < m; ++j) B.at(i).at(j) = A.at(j).at(i);
  }
  return B;
}

/*Computes product of matrices and vectors*/
vector<double> mul_vec_mat(const vector<auto> x, const vector<vector<double>> A){
  vector<double> v(A.at(0).size());
  for(long i = 0; i < A.at(0).size(); ++i){
    for(long j = 0; j < x.size(); ++j)v.at(i) += x.at(j) * A.at(j).at(i);
  }
  return v;
}

/*Computes product of two matriceses*/
vector<vector<double>> mul_mat(const vector<vector<double>> A, const vector<vector<double>> B)
{
  long n = A.size(), m = B.at(0).size();
  vector<vector<double>> C(n, vector<double>(m));
  for(long i = 0; i < n; ++i){
    for(long j = 0; j < m; ++j){
      for(long k = 0; k < B.size(); ++k) C.at(i).at(j) += A.at(i).at(k) * B.at(k).at(j);
    }
  }
  return C;
}

/*Preparation of Gaussian elimination over Z/2Z*/
bool have_1(const vector<long long> v, const int j){
  for(long i = 0; i < j; ++i) if(v.at(i) == 1) return true;
  return false;
}

/*Gaussian elimination of augumented matrices over Z/2Z*/
tuple<vector<vector<long long>>, vector<vector<long long>>> aug_gauss_mod2(vector<vector<long long>> A, vector<vector<long long>> B){
  long s;
  bool flag;
  for(long i = 0; i < A.size(); ++i){
    for(long j = 0; j < A.at(0).size(); ++j){
      A.at(i).at(j) = (A.at(i).at(j) % 2 + 2) % 2;
      B.at(i).at(j) = (B.at(i).at(j) % 2 + 2) % 2;
    }
  }
  for(long j = 0; j < A.at(0).size(); ++j){
    flag = false;
    for(long i = 0; i < A.size(); ++i){
      if(A.at(i).at(j) == 1 && (not flag) && (not have_1(A.at(i), j))){
        flag = true;
        s = i;
      }
    }
    for(long i = 0; i < A.size(); ++i){
      if(A.at(i).at(j) == 1 && flag && i != s){
        for(long k = 0; k < A.at(0).size(); ++k){
          A.at(i).at(k) ^= A.at(s).at(k);
          B.at(i).at(k) ^= B.at(s).at(k);
        }
      }
    }
  }
  return forward_as_tuple(A, B);
}

/*Computes basises of kaernel of matrices A mod 2*/
vector<vector<long long>> kernel_basis(const vector<vector<long long>> A){
  vector<vector<long long>> A_T, I, B, E, ker;
  A_T = trans(A);
  I = identity_mat_int(A_T.size());
  tie(B, E) = aug_gauss_mod2(A_T, I);
  for(long i = 0; i < B.size(); ++i){
    if(is_zero(B.at(i)))ker.push_back(E.at(i));
  }
  return ker;
}

/*Computes all elements in kernel of A*/
vector<vector<long long>> kernel_mod2(const vector<vector<long long>> A){
  vector<vector<long long>> ker_basis, v_all, ker;
  long long d;
  
  /*Computes all basis vectors of a kernel of A*/
  ker_basis = kernel_basis(A);
  vector<long long> tmp(ker_basis.at(0).size());
  d = ker_basis.size(); v_all = vector_all(d);
  
  /*Computes all linear combinations of basis vector of kernel of A modulo 2*/
  for(long i = 0; i < v_all.size(); ++i){
    for(long j = 0; j < d; ++j){
      for(long k = 0; k < ker_basis.at(0).size(); ++k) tmp.at(k) ^= v_all.at(i).at(j) * ker_basis.at(j).at(k);
    }
    ker.push_back(tmp);
    tmp = zero_vector_ZZ(tmp.size());
  }
  return ker;
}

/*Computes pseudo-inverse matrices of M*/
vector<vector<double>> mat_inv(vector<vector<double>> M){
  vector<vector<double>> X;
  int n = M.size(), m = M.at(0).size();
  mat_ZZ A, A_inv;
  ZZ d;
  A.SetDims(n, m);
  for(long i = 0; i < n; ++i) for(long j = 0; j < m; ++j) A[i][j] = conv<ZZ>((double)M.at(i).at(j));
  inv(d, A_inv, A);
  for(long i = 0; i < n; ++i) for(long j = 0; j < m; ++j) M.at(i).at(j) = conv<double>(A_inv[i][j]) / conv<double>(d);
  return M;
}

/*Eratosthenes' sieve algorithm*/
vector<long> esieve(const int M){
  vector<long> v(M); v.at(0) = 1;
  long long s = 1, p = s + 1, t, k;
  while(p * p <= M){
    t = s;
    while(1){
      t += p;
      if(t >= M) break;
      v.at(t) = 1;
    }
    for(k = 1; k < M; ++k) if(v.at(s + k) == 0) break;
    s += k; p = s + 1;
  }
  vector<long> p_vec {};
  vector<long>::iterator it;
  it = p_vec.begin();
  for(long i = M - 1; i >= 0; --i){
    if(v.at(i) == 0) it = p_vec.insert(it, i + 1);
  }
  return p_vec;
}

/*Generates factor basis*/
vector<long> factor_basis(const int n){
  int M = round(n * log((float)n));
  while(true){
    vector<long> p_list = esieve(M);
    if(p_list.size() >= n){
      p_list.insert(p_list.begin(), -1);
      return p_list;
    }
    ++M;
  }
}

/*Gram_Schmidt's orthogonalization algorithm*/
tuple<vector<vector<double>>, vector<vector<double>>> Gram_Schmidt(const vector<vector<double>> b)
{
  long n = b.size(), m = b.at(0).size();
  vector<vector<double>> GSOb(n, vector<double>(m));
  vector<vector<double>> mu = identity_mat(n);
  for(long i = 0; i < n; ++i){
    for(long j = 0; j < m; ++j) GSOb.at(i).at(j) = b.at(i).at(j);
    for(long j = 0; j < i; ++j){
      mu.at(i).at(j) = dot(b.at(i), GSOb.at(j)) / dot(GSOb.at(j), GSOb.at(j));
      for(long k = 0; k < m; ++k) GSOb.at(i).at(k) -= mu.at(i).at(j) * GSOb.at(j).at(k);
    }
  }
  return forward_as_tuple(GSOb, mu);
}

/*Generates a square norm vector of GSOb*/
vector<double> make_square_norm(const vector<vector<double>> GSOb){
  vector<double> B(GSOb.size());
  for(long i = 0; i < GSOb.size(); ++i) B.at(i) = dot(GSOb.at(i), GSOb.at(i));
  return B;
}

vector<vector<double>> LLL_reduce(const vector<vector<double>> b, const float delta){
  long n = b.size(), m = b.at(0).size();
  vector<vector<double>> B(n, vector<double>(m));
  Mat<ZZ> c; c.SetDims(n, m);
  for(long i = 0; i < n; ++i) for(long j = 0; j < m; ++j) c[i][j] = conv<ZZ>((double)b.at(i).at(j));
  LLL_FP(c, delta);
  for(long i = 0; i < n; ++i) for(long j = 0; j < m; ++j) B.at(i).at(j) = conv<double>(c[i][j]);
  return B;
}

vector<vector<double>> BKZ_reduce(const vector<vector<double>> b, const int beta, const float delta){
  long n = b.size(), m = b.at(0).size();
  Mat<ZZ> c; c.SetDims(n, m);
  vector<vector<double>> B(n, vector<double>(m));
  for(long i = 0; i < n; ++i) for(long j = 0; j < m; ++j) c[i][j] = conv<ZZ>((double)b.at(i).at(j));
  BKZ_FP(c, delta, beta);
  for(long i = 0; i < n; ++i) for(long j = 0; j < m; ++j) B.at(i).at(j) = conv<double>(c[i][j]);
  return B;
}

vector<double> Babai(const vector<vector<double>> b, const vector<double> w){
  vector<vector<double>> GSOb, mu;
  vector<double> t, v(w.size());
  float c;
  tie(GSOb, mu) = Gram_Schmidt(b);
  t = w;
  for(long i = b.size() - 1; i >= 0; --i){
    c = round(dot(t, GSOb.at(i)) / dot(GSOb.at(i), GSOb.at(i)));
    for(long j = 0; j < t.size(); ++j)t.at(j) -= c * GSOb.at(i).at(j);
  }
  for(long i = 0; i < v.size(); ++i)v.at(i) = w.at(i) - t.at(i);
  return v;
}

double babai_error(const vector<vector<double>> b, const vector<double> w){
  vector<double> t;
  vector<double> error_vec(w.size());
  t = Babai(b, w);
  for(long i = 0; i < error_vec.size(); ++i) error_vec.at(i) = fabs(t.at(i) - w.at(i));
  return sqrt(dot(error_vec, error_vec));
}

/*Converts coefficient vector to lattice vector*/
vector<double> coef2lat(const vector<auto> v, const vector<vector<auto>> b){
  vector<double> x(b.at(0).size());
  for(long i = 0; i < b.size(); ++i){
    for(long j = 0; j < b.at(0).size(); ++j) x.at(j) += v.at(i) * b.at(i).at(j);
  }
  return x;
}

vector<double> coef2ldlat(const vector<auto> v, const vector<vector<auto>> b){
  vector<double> x(b.at(0).size());
  for(long i = 0; i < b.size(); ++i){
    for(long j = 0; j < b.at(0).size(); ++j) x.at(j) += v.at(i) * b.at(i).at(j);
  }
  return x;
}

vector<double> lat2coef(const vector<auto> v, const vector<vector<auto>> b){
  long n = b.size();
  vector<double> x(n);
  for(long i = 0; i < n; ++i) x.at(i) = v.at(i) / b.at(i).at(i);
  return x;
}

/*Enumerates closest vectors*/
vector<double> ENUM_CVP(const vector<vector<double>> mu, const vector<double> B, const float R, const vector<double> a){
  long n = B.size();
  vector<vector<double>> sigma(n + 2, vector<double>(n + 1)), CVP_list;
  vector<long long> r(n + 1), w(n + 1);
  vector<double> c(n + 1), rho(n + 2), v(n + 1), x(n), u;
  long k;
  v.at(1) = 1;
  for(long i = 0; i <= n; ++i) r.at(i) = i;
  for(k = n; k >= 1; --k){
    for(long i = n; i >= k + 1; --i) sigma.at(i).at(k) = sigma.at(i + 1).at(k) + (a.at(i - 1) - v.at(i)) * mu.at(i - 1).at(k - 1);
    c.at(k) = a.at(k - 1) + sigma.at(k + 1).at(k);
    v.at(k) = round(c.at(k));
    w.at(k) = 1;
    rho.at(k) = rho.at(k + 1) + (c.at(k) - v.at(k)) * (c.at(k) - v.at(k)) * B.at(k - 1);
  }
  k = 1;
  while(true){
    rho.at(k) = rho.at(k + 1) + (c.at(k) - v.at(k)) * (c.at(k) - v.at(k)) * B.at(k - 1);
    if(rho.at(k) < R * R || fabs(rho.at(k) - R * R) < EPSILON){
      if(k == 1){
        for(long i = 1; i <= n; ++i) x.at(i - 1) = v.at(i);
        return x;
      }
      --k;
      r.at(k - 1) = max(r.at(k - 1), r.at(k));
      for(long i = r.at(k); i >= k + 1; --i) sigma.at(i).at(k) = sigma.at(i + 1).at(k) + (a.at(i - 1) - v.at(i)) * mu.at(i - 1).at(k - 1);
      c.at(k) = a.at(k - 1) + sigma.at(k + 1).at(k);
      v.at(k) = round(c.at(k));
      w.at(k) = 1;
    }else{
      ++k;
      if(k == n + 1){
        x.clear();
        return x;
      }
      r.at(k - 1) = k;
      if(v.at(k) > c.at(k)) v.at(k) -= w.at(k);else v.at(k) += w.at(k);
      ++w.at(k);
    }
  }
}

/*Generates coefficient vectors of close vectors*/
vector<vector<double>> ENUM_CVP_all(const vector<vector<double>> mu, const vector<double> B, float R, const vector<double> a, const vector<double> t, const vector<vector<double>> b){
  long n = B.size();
  vector<vector<double>> CVP_list;
  vector<double> ENUM_CVP_v(n), pre_ENUM_CVP_v(n), x(n), c, d(t.size());
  while(true){
    for(long i = 0; i < n; ++i) pre_ENUM_CVP_v.at(i) = ENUM_CVP_v.at(i);
    ENUM_CVP_v = ENUM_CVP(mu, B, R, a);
    if(ENUM_CVP_v.empty()) return CVP_list;
    if(not is_included(ENUM_CVP_v, CVP_list))CVP_list.push_back(ENUM_CVP_v);

    /*distance between target vector and close vector*/
    c = coef2lat(ENUM_CVP_v, b);
    for(long i = 0; i < t.size(); ++i) d.at(i) = c.at(i) - t.at(i);
    R = fmin(sqrt(dot(d, d)), R * 0.9);
  }
}


/******************************************
主要なsub-routine始め
*******************************************/
vector<vector<double>> imp_prime_mat(const float c, const vector<long> p, const int n){
  vector<vector<double>> A(n, vector<double>(n + 1));
  vector<double> diag(n);//対角成分
  mt19937_64 get_rand_mt;

  for(long i = 0; i < n; ++i) diag.at(i) = ceil((i + 1) / 2.0);

  random_device seed_gen;
  mt19937 engine(seed_gen());
  shuffle(diag.begin(), diag.end(), engine); //乱数的に並び替え

  for(long i = 0; i < n; ++i)A.at(i).at(i) = diag.at(i);
  
  for(long i = 0; i < n; ++i)A.at(i).at(n) = round(pow(10, c) * log((float)p.at(i + 1)));

  return A;
}

/*Generates a target vector*/
vector<double> target(const ZZ N, const float c, const int n){
  vector<double> t(n + 1);
  t.at(n) = round(pow(10, c) * log(conv<ZZ>(N)));
  return t;
}

bool sr_test(const ZZ u, const ZZ v, const ZZ N, const vector<long> p, const int n){
  ZZ a = u - N * v, r;
  for(long i = 0; i < n; ++i){
    r = conv<ZZ>(p.at(i + 1));
    while(a % r == 0 && a != 0) a /= r;
  }
  if(a == 1 || a == -1) return true;
  return false;
}

/******************************************
主要なsub-routine終わり
*******************************************/



/******************************************
素因数分解
*******************************************/
int main(int argc, char **argv){
  int loop_times = 0, num = 0, n, K, e1, e2, s;
  float c = 5.0;
  ZZ X, Y, q, N, a, u, v, U, r;
  vector<double> B, w;
  vector<vector<double>> BB, GSOb, mu, close_vecs, L;
  vector<double> t, e;
  N = conv<ZZ>(argv[1]);
  n = ceil(log2(conv<double>(N))) / (1.5 * log(ceil(log2(conv<double>(N))) / 1.5));
  if(argc == 3) n = atoi(argv[2]);
  if(argc == 4) c = atof(argv[3]);

  K = round(2.0 * pow(n, 2.0));
  vector<long long> vec(K), ee(K);
  vector<vector<long long>> kernel, A, prekernel;
  vector<long> p = factor_basis(K);

  cout << "Dimension = " << n << "\nNumber of primes = " << K << "\nN = " << N << "\nbit size = " << ceil(log2(conv<double>(N))) << endl;

  while(true){
    ++loop_times;

    w = zero_vector_RR(n);

    /*Composition of CVP*/
    L = imp_prime_mat(c, p, n); t = target(N, c, n);

    /*Reduces lattice basis*/
    #if BKZ_REDUCE
    if(n >= 30)BB = BKZ_reduce(L, 20, 0.99);else BB = LLL_reduce(L, 0.99);
    #endif
    #if not BKZ_REDUCE
    BB = LLL_reduce(L, 0.99);
    #endif
    tie(GSOb, mu) = Gram_Schmidt(BB); B = make_square_norm(GSOb);
    
    /*Composition of coefficient vectors of target vectors*/
    w = mul_vec_mat(t, mul_mat(trans_RR(BB), mat_inv(mul_mat(BB, trans_RR(BB)))));

    /*Computes approximate solutions of CVP*/
    close_vecs = ENUM_CVP_all(mu, B, babai_error(BB, t), w, t, BB);

    for(long i = 0; i < close_vecs.size(); ++i){
      e = lat2coef(coef2ldlat(close_vecs.at(i), BB), L);

      /*Composition of a pair of integers (u, v) correspond to close_vec[i]*/
      u = v = 1;
     
     for(long j = 0; j < n; ++j){
        if(e.at(j) > 0)u *= conv<ZZ>(pow(p.at(j + 1), e.at(j)));
        else if(e.at(j) < 0)v *= conv<ZZ>(pow(p.at(j + 1), -e.at(j)));
      }

      if(u > 0 && v > 0 && sr_test(u, v, N, p, K)){//if (u, v) is sr-pair, then
        vec = zero_vector_ZZ(K);
        for(long j = 0; j < K; ++j){
          e1 = e2 = 0;
          U = u;
          r = conv<ZZ>(p.at(j + 1));
          
          /*Computes exponents of p_{j} that is prime factor of u*/
          if(j < n)while(U % r == 0 && U != 0){++e1; U /= r;}
          
          /*Computes exponents of p_{j} that is prime factor of a=u-vN*/
          a = u - v * N;
          while(a % r == 0 && a != 0){++e2; a /= r;}
          vec.at(j) = e2 - e1;//exponential vector of each prime factor
        }
        
        /*Appends exponential vector of each prime factor to A*/
        if(find(A.begin(), A.end(), vec) == A.end()){
          ++num;
          cout << num << " sr-pairs were found. (" << K << "pairs are needed)" << endl;
          A.push_back(vec);//Append to A
        }
      }
    }

    /*prime factorization*/
    if(A.size() >= K){
      ee = zero_vector_ZZ(K);
      kernel = kernel_mod2(trans(A));
      for(long i = 0; i < prekernel.size(); ++i)prekernel.at(i).resize(kernel.at(0).size());

      if(kernel != prekernel){
      /*Verifies any t = T that satisfies\sum_{j=1}^{ee.size()}t_j(e_{i.j}-e'_{i.j})=0(mod\ 2)*/
        for(long i = 0; i < kernel.size(); ++i){
          if(find(prekernel.begin(), prekernel.end(), kernel.at(i)) == prekernel.end()){
            for(long j = 0; j < kernel.at(0).size(); ++j){
              if(kernel.at(i).at(j) == 1) for(long k = 0; k < K; ++k) ee.at(k) += A.at(j).at(k);
            }
            X = Y = 1;
            for(long j = 0; j < K; ++j){
              if(ee.at(j) > 0) X = MulMod(X, PowerMod(conv<ZZ>(p.at(j + 1)), conv<ZZ>((long)ee.at(j) / 2), N), N);
              if(ee.at(j) < 0) Y = MulMod(Y, PowerMod(conv<ZZ>(p.at(j + 1)), conv<ZZ>(-(long)ee.at(j) / 2), N), N);
            }
            cout << "X = " << X << ", Y = " << Y << endl;
            if(X != Y) q = GCD(X - Y, N);else q = GCD(X + Y, N);

            if(1 < q && q < N){
              cout << "c = " << c << ", beta = " << 2.0 << endl << N << " = " << q << " * " << N / q << endl;
              cout << "loop times = " << loop_times << endl;
              return 0;
            }
            ee = zero_vector_ZZ(K);
          }
        }
        prekernel = kernel;
      }
    }
  }
}
