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

using namespace std;

bool vector_not_equal(const vector<double> x, const vector<double> y){
  int n = x.size();
  for(int i = 0; i < n; ++i){
    if(fabs(x.at(i) - y.at(i)) > EPSILON)return true;
  }
  return false;
}

bool not_included(const vector<double> x, const vector<vector<double>> vecs){
  int n = vecs.size();
  if(n == 0)return true;
  for(int i = 0; i < n; ++i){
    if(vector_not_equal(x, vecs.at(i)))return true;
  }
  return false;
}

vector<double> zero_vector_RR(const int n){
  vector<double> z(n);
  return z;
}

vector<int> zero_vector_ZZ(const int n){
  vector<int> z(n);
  return z;
}

/*Generates all combinations of n 0 and 1*/
vector<vector<int>> vector_all(const int n){
  vector<int> tmp(n), v(n);
  vector<vector<int>> v_all;
  for(int i = 0; i < n; ++i){
    tmp.at(i) = 1; v = tmp;
    do v_all.push_back(v); while(next_permutation(v.begin(), v.end()));
  }
  return v_all;
}

/*Tests v is zero-vector or not*/
bool is_zero(const vector<int> v){
  int n = v.size();
  for(int i = 0; i < n; ++i) if(v.at(i) != 0) return false;
  return true;
}

/*Generates an identity matrices*/
vector<vector<double>> identity_mat(const int n){
  vector<vector<double>> A(n, vector<double>(n));
  for(int i = 0; i < n; ++i) A.at(i).at(i) = 1.0;
  return A;
}

/*Generates an identity matrices*/
vector<vector<int>> identity_mat_int(const int n){
  vector<vector<int>> A(n, vector<int>(n));
  for(int i = 0; i < n; ++i) A.at(i).at(i) = 1;
  return A;
}

#if 0
/*Prints a matrices A*/
void print_mat(const vector<vector<auto>> A){
  cout << "[\n";
  for(int i = 0; i < A.size(); ++i){
    cout << "[";
    for(int j = 0; j < A.at(0).size(); ++j) cout << A.at(i).at(j) << " ";
    cout << "]" << endl;
  }
  cout << "]\n";
}

/*Prints a vector v*/
void print_vec(const vector<auto> v){
  cout << "[";
  for(int i = 0; i < v.size(); ++i) cout << v.at(i) << " ";
  cout << "]\n";
}
#endif

/*Computes inner product of vectors x and y*/
double dot(const vector<double> x, const vector<double> y)
{
  double z = 0.0;
  int n = x.size();
  for(int i = 0; i < n; ++i) z += x.at(i) * y.at(i);
  return z;
}

/*Computes transpose matrices of A*/
vector<vector<int>> trans(const vector<vector<int>> A)
{
  int n = A.at(0).size(), m = A.size();
  vector<vector<int>> B(n, vector<int>(m));
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < m; ++j) B.at(i).at(j) = A.at(j).at(i);
  }
  return B;
}

/*Computes transpose matrices of A*/
vector<vector<double>> trans_RR(const vector<vector<double>> A)
{
  int n = A.at(0).size(), m = A.size();
  vector<vector<double>> B(n, vector<double>(m));
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < m; ++j) B.at(i).at(j) = A.at(j).at(i);
  }
  return B;
}

/*Computes product of matrices and vectors*/
vector<double> mul_vec_mat(const vector<double> x, const vector<vector<double>> A){
  int m = A.at(0).size(), n = x.size();
  vector<double> v(m);
  for(int i = 0; i < m; ++i){
    for(int j = 0; j < n; ++j)v.at(i) += x.at(j) * A.at(j).at(i);
  }
  return v;
}

/*Computes product of two matriceses*/
vector<vector<double>> mul_mat(const vector<vector<double>> A, const vector<vector<double>> B)
{
  int n = A.size(), m = B.at(0).size(), l = B.size();
  vector<vector<double>> C(n, vector<double>(m));
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < m; ++j){
      for(int k = 0; k < l; ++k) C.at(i).at(j) += A.at(i).at(k) * B.at(k).at(j);
    }
  }
  return C;
}

/*Preparation of Gaussian elimination over Z/2Z*/
bool have_1(const vector<int> v, const int j){
  for(int i = 0; i < j; ++i) if(v.at(i) == 1) return true;
  return false;
}

/*Gaussian elimination of augumented matrices over Z/2Z*/
void aug_gauss_mod2(vector<vector<int>>& A, vector<vector<int>>& B){
  int s, n = A.size(), m = A.at(0).size();
  bool flag;
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < m; ++j){
      A.at(i).at(j) = fabs(A.at(i).at(j) % 2);
      B.at(i).at(j) = fabs(B.at(i).at(j) % 2);
    }
  }
  for(int j = 0; j < m; ++j){
    flag = false;
    for(int i = 0; i < n; ++i){
      if(A.at(i).at(j) == 1 && (! flag) && (! have_1(A.at(i), j))){
        flag = true;
        s = i;
      }
    }
    for(int i = 0; i < n; ++i){
      if(A.at(i).at(j) == 1 && flag && i != s){
        for(int k = 0; k < m; ++k){
          A.at(i).at(k) ^= A.at(s).at(k);
          B.at(i).at(k) ^= B.at(s).at(k);
        }
      }
    }
  }
}

/*Computes basises of kaernel of matrices A mod 2*/
vector<vector<int>> kernel_basis(const vector<vector<int>> A){
  int n = A.at(0).size();
  vector<vector<int>> B = trans(A), E = identity_mat_int(n), ker;
  aug_gauss_mod2(B, E);
  
  for(int i = 0; i < n; ++i){
    if(is_zero(B.at(i))) ker.push_back(E.at(i));
  }
  return ker;
}

/*Computes all elements in kernel of A*/
vector<vector<int>> kernel_mod2(const vector<vector<int>> A){
  vector<vector<int>> ker_basis = kernel_basis(A), v_all, ker;
  
  /*Computes all basis vectors of a kernel of A*/
  int d = ker_basis.size(), m = ker_basis.at(0).size();
  vector<int> tmp(m);
  v_all = vector_all(d);
  
  /*Computes all linear combinations of basis vector of kernel of A modulo 2*/
  for(int i = 0; i < v_all.size(); ++i){
    for(int j = 0; j < d; ++j){
      for(int k = 0; k < m; ++k) tmp.at(k) ^= v_all.at(i).at(j) * ker_basis.at(j).at(k);
    }
    ker.push_back(tmp);
    tmp = zero_vector_ZZ(tmp.size());
  }
  return ker;
}

/*Computes inverse matrices of M*/
vector<vector<double>> mat_inv(vector<vector<double>> M){
  vector<vector<double>> X;
  int n = M.size(), m = M.at(0).size();
  NTL::mat_ZZ A, A_inv;
  NTL::ZZ d;
  A.SetDims(n, m);
  for(int i = 0; i < n; ++i) for(int j = 0; j < m; ++j) A[i][j] = NTL::conv<NTL::ZZ>(M.at(i).at(j));
  NTL::inv(d, A_inv, A);
  double e = 1.0 / conv<double>(d);
  for(int i = 0; i < n; ++i) for(int j = 0; j < m; ++j) M.at(i).at(j) = NTL::conv<double>(A_inv[i][j]) * e;
  return M;
}

/*Eratosthenes' sieve algorithm*/
vector<int> esieve(const int M){
  vector<int> v(M); v.at(0) = 1;
  int s = 1, p = s + 1, t, k;
  while(p * p <= M){
    t = s;
    while(true){
      t += p;
      if(t >= M) break;
      v.at(t) = 1;
    }
    for(k = 1; k < M; ++k) if(v.at(s + k) == 0) break;
    s += k; p = s + 1;
  }
  vector<int> p_vec {};
  vector<int>::iterator it;
  it = p_vec.begin();
  for(int i = M - 1; i >= 0; --i){
    if(v.at(i) == 0) it = p_vec.insert(it, i + 1);
  }
  return p_vec;
}

/*Generates factor basis*/
vector<int> factor_basis(const int n){
  int M = round(n * log(n));
  vector<int> p_list;
  while(true){
    p_list = esieve(M);
    if(p_list.size() >= n){
      p_list.insert(p_list.begin(), -1);
      return p_list;
    }
    ++M;
  }
}

/*Gram_Schmidt's orthogonalization algorithm*/
tuple<vector<vector<double>>, vector<vector<double>>> Gram_Schmidt(const vector<vector<double>> b){
  int n = b.size(), m = b.at(0).size();
  vector<vector<double>> GSOb(n, vector<double>(m)), mu = identity_mat(n);
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < m; ++j) GSOb.at(i).at(j) = b.at(i).at(j);
    for(int j = 0; j < i; ++j){
      mu.at(i).at(j) = dot(b.at(i), GSOb.at(j)) / dot(GSOb.at(j), GSOb.at(j));
      for(int k = 0; k < m; ++k) GSOb.at(i).at(k) -= mu.at(i).at(j) * GSOb.at(j).at(k);
    }
  }
  return forward_as_tuple(GSOb, mu);
}

/*Gram_Schmidt's orthogonalization algorithm*/
tuple<vector<double>, vector<vector<double>>> Gram_Schmidt_squared(const vector<vector<double>> b){
  int n = b.size(), m = b.at(0).size();
  vector<double> B(n);
  vector<vector<double>> GSOb(n, vector<double>(m)), mu = identity_mat(n);
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < m; ++j) GSOb.at(i).at(j) = b.at(i).at(j);
    for(int j = 0; j < i; ++j){
      mu.at(i).at(j) = dot(b.at(i), GSOb.at(j)) / dot(GSOb.at(j), GSOb.at(j));
      for(int k = 0; k < m; ++k) GSOb.at(i).at(k) -= mu.at(i).at(j) * GSOb.at(j).at(k);
    }
    B.at(i) = dot(GSOb.at(i), GSOb.at(i));
  }
  return forward_as_tuple(B, mu);
}

vector<vector<double>> LLL_reduce(const vector<vector<double>> b, const float delta){
  int n = b.size(), m = b.at(0).size();
  vector<vector<double>> B(n, vector<double>(m));
  NTL::Mat<NTL::ZZ> c; c.SetDims(n, m);
  for(int i = 0; i < n; ++i) for(int j = 0; j < m; ++j) c[i][j] = NTL::conv<NTL::ZZ>(b.at(i).at(j));
  NTL::LLL_FP(c, delta);
  for(int i = 0; i < n; ++i) for(int j = 0; j < m; ++j) B.at(i).at(j) = NTL::conv<double>(c[i][j]);
  return B;
}

vector<vector<double>> BKZ_reduce(const vector<vector<double>> b, const int beta, const float delta){
  int n = b.size(), m = b.at(0).size();
  NTL::Mat<NTL::ZZ> c; c.SetDims(n, m);
  vector<vector<double>> B(n, vector<double>(m));
  for(int i = 0; i < n; ++i) for(int j = 0; j < m; ++j) c[i][j] = NTL::conv<NTL::ZZ>(b.at(i).at(j));
  NTL::BKZ_FP(c, delta, beta);
  for(int i = 0; i < n; ++i) for(int j = 0; j < m; ++j) B.at(i).at(j) = NTL::conv<double>(c[i][j]);
  return B;
}

vector<double> Babai(const vector<vector<double>> b, const vector<double> w){
  vector<vector<double>> GSOb, mu;
  int n = b.size(), m = w.size();
  vector<double> t = w, v(m);
  double c;
  tie(GSOb, mu) = Gram_Schmidt(b);
  for(int i = n - 1; i >= 0; --i){
    c = round(dot(t, GSOb.at(i)) / dot(GSOb.at(i), GSOb.at(i)));
    for(int j = 0; j < m; ++j)t.at(j) -= c * GSOb.at(i).at(j);
  }
  for(int i = 0; i < m; ++i)v.at(i) = w.at(i) - t.at(i);
  return v;
}

double babai_error(const vector<vector<double>> b, const vector<double> w){
  vector<double> t = Babai(b, w), error_vec(w.size());
  int n = error_vec.size();
  for(int i = 0; i < n; ++i) error_vec.at(i) = t.at(i) - w.at(i);
  return sqrt(dot(error_vec, error_vec));
}

/*Converts coefficient vector to lattice vector*/
vector<double> coef2lat(const vector<auto> v, const vector<vector<auto>> b){
  int n = b.size(), m = b.at(0).size();
  vector<double> x(m);
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < m; ++j) x.at(j) += v.at(i) * b.at(i).at(j);
  }
  return x;
}

vector<double> coef2ldlat(const vector<auto> v, const vector<vector<auto>> b){
  int n = b.size(), m = b.at(0).size();
  vector<double> x(m);
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < m; ++j) x.at(j) += v.at(i) * b.at(i).at(j);
  }
  return x;
}

vector<double> lat2coef(const vector<auto> v, const vector<vector<auto>> b){
  int n = b.size();
  vector<double> x(n);
  for(int i = 0; i < n; ++i) x.at(i) = v.at(i) / b.at(i).at(i);
  return x;
}

/*Enumerates closest vectors*/
vector<double> ENUM_CVP(const vector<vector<double>> mu, const vector<double> B, const double R, const vector<double> a){
  int n = B.size(), k;
  vector<vector<double>> sigma(n + 2, vector<double>(n + 1)), CVP_list;
  vector<int> r(n + 1), w(n + 1);
  vector<double> c(n + 1), rho(n + 2), v(n + 1), x(n), u;
  v.at(1) = 1;
  for(int i = 0; i <= n; ++i) r.at(i) = i;
  for(k = n; k >= 1; --k){
    for(int i = n; i >= k + 1; --i) sigma.at(i).at(k) = sigma.at(i + 1).at(k) + (a.at(i - 1) - v.at(i)) * mu.at(i - 1).at(k - 1);
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
        for(int i = 1; i <= n; ++i) x.at(i - 1) = v.at(i);
        return x;
      }
      --k;
      r.at(k - 1) = max(r.at(k - 1), r.at(k));
      for(int i = r.at(k); i >= k + 1; --i) sigma.at(i).at(k) = sigma.at(i + 1).at(k) + (a.at(i - 1) - v.at(i)) * mu.at(i - 1).at(k - 1);
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
vector<vector<double>> ENUM_CVP_all(const vector<vector<double>> mu, const vector<double> B, double R, const vector<double> a, const vector<double> t, const vector<vector<double>> b){
  int n = B.size(), m = t.size();
  vector<vector<double>> CVP_list;
  vector<double> ENUM_CVP_v(n), pre_ENUM_CVP_v(n), x(n), c, d(m);
  while(true){
    for(int i = 0; i < n; ++i) pre_ENUM_CVP_v.at(i) = ENUM_CVP_v.at(i);
    ENUM_CVP_v = ENUM_CVP(mu, B, R, a);
    if(ENUM_CVP_v.empty()) return CVP_list;
    if(not_included(ENUM_CVP_v, CVP_list))CVP_list.push_back(ENUM_CVP_v);

    /*distance between target vector and close vector*/
    c = coef2lat(ENUM_CVP_v, b);
    for(int i = 0; i < m; ++i) d.at(i) = c.at(i) - t.at(i);
    R = fmin(sqrt(dot(d, d)), R * 0.9);
  }
}


/******************************************
主要なsub-routine始め
*******************************************/
vector<vector<double>> imp_prime_mat(const double c, const vector<int> p, const int n){
  vector<vector<double>> A(n, vector<double>(n + 1));
  vector<double> diag(n);//対角成分

  for(int i = 0; i < n; ++i) diag.at(i) = ceil((i + 1) * 0.5);

  random_device seed_gen;
  mt19937_64 engine(seed_gen());
  shuffle(diag.begin(), diag.end(), engine); //乱数的に並び替え

  if(c == 5.0){
    for(int i = 0; i < n; ++i){
      A.at(i).at(i) = diag.at(i);
      A.at(i).at(n) = round(100000 * log(p.at(i + 1)));
    }
  }else{
    for(int i = 0; i < n; ++i){
      A.at(i).at(i) = diag.at(i);
      A.at(i).at(n) = round(pow(10, c) * log(p.at(i + 1)));
    }
  }

  return A;
}

/*Generates a target vector*/
vector<double> target(const NTL::ZZ N, const double c, const int n){
  vector<double> t(n + 1);
  if(c == 5.0){
    t.at(n) = round(100000 * log(N));
  }else{
    t.at(n) = round(pow(10, c) * log(N));
  }
  return t;
}

bool sr_test(const NTL::ZZ u, const NTL::ZZ v, const NTL::ZZ N, const vector<int> p, const int n){
  NTL::ZZ a = u - N * v, r;
  for(int i = 0; i < n; ++i){
    r = NTL::conv<NTL::ZZ>(p.at(i + 1));
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
  double c = 5.0;
  NTL::ZZ X, Y, q, N, a, u, v, U;
  vector<double> B, w, t, e;
  vector<vector<double>> BB, mu, close_vecs, L;
  N = NTL::conv<NTL::ZZ>(argv[1]);
  n = ceil(log2(NTL::conv<double>(N))) / (1.5 * log(ceil(log2(NTL::conv<double>(N))) * 0.6666667));
  if(argc >= 3) n = atoi(argv[2]);
  if(argc >= 4) c = atof(argv[3]);

  K = 2 * n * n;
  vector<int> vec(K), ee(K), p = factor_basis(K);
  vector<vector<int>> kernel, A, prekernel;

  cout << "Dimension = " << n << "\nNumber of primes = " << K << "\nN = " << N << "\nbit size = " << ceil(log2(conv<double>(N))) << endl;

  while(true){
    ++loop_times;
    w = zero_vector_RR(n);

    /*Composition of CVP*/
    L = imp_prime_mat(c, p, n); t = target(N, c, n);

    /*Reduces lattice basis*/
    if(n > 30) BB = BKZ_reduce(L, 20, 0.99);else BB = LLL_reduce(L, 0.99);

    tie(B, mu) =  Gram_Schmidt_squared(BB);
    
    /*Composition of coefficient vectors of target vectors*/
    w = mul_vec_mat(t, mul_mat(trans_RR(BB), mat_inv(mul_mat(BB, trans_RR(BB)))));

    /*Computes approximate solutions of CVP*/
    close_vecs = ENUM_CVP_all(mu, B, babai_error(BB, t), w, t, BB);

    for(int i = 0; i < close_vecs.size(); ++i){
      e = lat2coef(coef2ldlat(close_vecs.at(i), BB), L);

      /*Composition of a pair of integers (u, v) correspond to close_vec[i]*/
      u = v = 1;
      for(int j = 0; j < n; ++j){
        if(e.at(j) > 0) u *= NTL::conv<NTL::ZZ>(pow(p.at(j + 1), e.at(j)));
        else if(e.at(j) < 0) v *= NTL::conv<NTL::ZZ>(pow(p.at(j + 1), -e.at(j)));
      }

      if(u > 0 && v > 0 && sr_test(u, v, N, p, K)){//if (u, v) is sr-pair, then
        vec = zero_vector_ZZ(K);
        for(int j = 0; j < K; ++j){
          e1 = e2 = 0;
          U = u;
          q = NTL::conv<NTL::ZZ>(p.at(j + 1));
          
          /*Computes exponents of p_{j} that is prime factor of u*/
          if(j < n)while(U % q == 0 && U != 0){++e1; U /= q;}
          
          /*Computes exponents of p_{j} that is prime factor of a=u-vN*/
          a = u - v * N;
          while(a % q == 0 && a != 0){++e2; a /= q;}
          vec.at(j) = e2 - e1;//exponential vector of each prime factor
        }
        
        /*Appends exponential vector of each prime factor to A*/
        if(find(A.begin(), A.end(), vec) == A.end()){
          ++num;
          printf("%d spair were found. (%d pairs are needed.)\n",num, K);
          A.push_back(vec);//Append to A
        }
      }
    }

    /*prime factorization*/
    if(A.size() >= K){
      ee = zero_vector_ZZ(K);
      kernel = kernel_mod2(trans(A)); s = kernel.size();
      for(int i = 0; i < prekernel.size(); ++i)prekernel.at(i).resize(kernel.at(0).size());

      if(kernel != prekernel){
      /*Verifies any t = T that satisfies\sum_{j=1}^{ee.size()}t_j(e_{i.j}-e'_{i.j})=0(mod\ 2)*/
        for(int i = 0; i < s; ++i){
          if(find(prekernel.begin(), prekernel.end(), kernel.at(i)) == prekernel.end()){
            for(int j = 0; j < K; ++j){
              if(kernel.at(i).at(j) == 1) for(int k = 0; k < K; ++k) ee.at(k) += A.at(j).at(k);
            }
            X = Y = 1;
            for(int j = 0; j < K; ++j){
              if(ee.at(j) > 0) X = NTL::MulMod(X, NTL::PowerMod(NTL::conv<NTL::ZZ>(p.at(j + 1)), NTL::conv<NTL::ZZ>((int)ee.at(j) / 2), N), N);
              else if(ee.at(j) < 0) Y = NTL::MulMod(Y, NTL::PowerMod(NTL::conv<NTL::ZZ>(p.at(j + 1)), NTL::conv<NTL::ZZ>(-(int)ee.at(j) / 2), N), N);
            }
            cout << "X = " << X << ", Y = " << Y << endl;
            if(X != Y) q = NTL::GCD(X - Y, N);else q = NTL::GCD(X + Y, N);

            if(1 < q && q < N){
              cout << "c = " << c << ", beta = 2\nN = " << q << " * " << N / q << "\nloop times = " << loop_times << "\nnumber of sr-pairs = " << num << endl;
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
