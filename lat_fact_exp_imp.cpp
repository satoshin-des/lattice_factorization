#include<iostream>
#include <vector>
#include <cmath>
#include<random>
#include <iomanip>
#include<time.h>
#include<tuple>
#include<algorithm>
#include<ctime>
#include<cstdlib>
#include<string>
#include<NTL/ZZ.h>
#include<NTL/mat_ZZ.h>
#include<NTL/matrix.h>
#include<NTL/vector.h>
#include<NTL/LLL.h>
#define EPSILON 1e-6
#define original 0

/*
~$ g++ -O3 -o lat_fact_imp.exe lat_fact_imp.cpp -lntl
~$ ./lat_fact_imp.exe [N] [c] [alpha] [beta] [TIME] [SHIFT]
*/

using namespace std;
using namespace NTL;

/******************************************
準備始め
*******************************************/

/*Generates all combinations of n 0 and 1*/
vector<vector<long long>> vector_all(long long n){
  vector<long long> tmp(n), v(n);
  vector<vector<long long>> v_all;
  for(long i = 0; i < n; i++){
    tmp.at(i) = 1;
    v = tmp;
    do{
      v_all.push_back(v);
    }while(next_permutation(v.begin(), v.end()));
  }
  return v_all;
}

/*Computes m*n mod N*/
long long mulm(long long m, long long n, long long N){
  unsigned long long s = 0, qm, qn, rm, rn, tmp1, tmp2, tmp3, tmp4 = 0, tmp5 = 0, tmp6, qq, rq, sum;
  qm = m / 10000; rm = m % 10000;
  qn = n / 10000; rn = n % 10000; 
  qq = qm /1000; rq = qm % 1000; //qm = 1000 * qq + rq

  s += (rn * rm) % N;

  /*Computes 100000000 * qn * qm mod N*/
  tmp1 = qn;
  for(int i = 0; i < 8; i++)tmp1 = (tmp1 * 10) % N;//tmp1 = 100000000 * qn mod N
  for(long long i = 0; i < qq; i++)tmp5 = (tmp5 + tmp1) % N;//tmp5 = tmp1 * qq mod N
  for(int i = 0; i < 3; i++)tmp5 = (10 * tmp5) % N;//tmp5 = 1000 * tmp1 * qq mod N
  sum = 0; for(long long i = 0; i < rq; i++)sum = (sum + tmp1) % N;//sum = tmp1 * rq mod N
  tmp4 = (tmp5 + sum) % N;//tmp4 = tmp1 * (1000 * qq + rq) mod N = 100000000 * qn * qm mod N
  s = (s + tmp4) % N;

  /*Computes 10000 * qn * rm, 10000 * qm * rn mod N*/
  tmp2 = qm; tmp3 = qn;
  for(int i = 0; i < 4; i++){tmp2 = (tmp2 * 10) % N; tmp3 = (tmp3 * 10) % N;}
  sum = 0; for(long long i = 0; i < rn; i++){sum = (sum + tmp2) % N;} tmp2 = sum;
  sum = 0; for(long long i = 0; i < rm; i++){sum = (sum + tmp3) % N;} tmp3 = sum;
  s = (s + tmp2) % N; s = (s + tmp3) % N;

  return s;
}

/*Computes GCD of a and b*/
long long gcd(long long a, long long b) {
  if (a % b == 0) {return b;} else {return gcd(b, a % b);}
}

/*Computes number of times that d can divid N*/
long long div_count(long long N, long long d){
  long long count = 0;
  if(d == -1){
    if(N < 0){return 1;}else{return 0;}
  }
  while(N % d == 0){
    count++;
    N /= d;
  }
  return count;
}

/*Computes A^M mod N*/
long long powm(long A, long M, long N){
  unsigned long long S = 1;
  unsigned long long T = A;
  while(M > 0){
    if(M % 2 == 1){S = mulm(S, T, N);}
    T = mulm(T, T, N);
    M = ceil(M / 2);
  }
  return S;
}

/*Tests v is zero-vector or not*/
bool is_zero(vector<long long> v){
  for(long i = 0; i < v.size(); i++){
    if(v.at(i) != 0) return false;
  }
  return true;
}

/*Generates an identity matrices*/
vector<vector<long double>> identity_mat(long n){
  vector<vector<long double>> A(n, vector<long double>(n));
  for(long i = 0; i < n; i++){
    A.at(i).at(i) = 1.0;
  }
  return A;
}

/*Generates an identity matrices*/
vector<vector<long long>> identity_mat_int(long n){
  vector<vector<long long>> A(n, vector<long long>(n));
  for(long i = 0; i < n; i++){
    A.at(i).at(i) = 1;
  }
  return A;
}

/*Prints a matrices A*/
void print_mat(vector<vector<long long>> A){
  cout << "[\n";
  for(long i = 0; i < A.size(); i++){
    cout << "[";
    for(long j = 0; j < A.at(0).size(); j++) cout << A.at(i).at(j) << " ";
    cout << "]" << endl;
  }
  cout << "]\n";
}

/*Prints a vector v*/
void print_vec(vector<long long> v){
  cout << "[";
  for(long i = 0; i < v.size(); i++) cout << v.at(i) << " ";
  cout << "]\n";
}

/*Computes inner product of vectors x and y*/
long double dot(vector<long double> x, vector<long double> y)
{
  long double z = 0.0;
  for(long i = 0; i < x.size(); i++){
    z += x.at(i) * y.at(i);
  }
  return z;
}

/*Computes transpose matrices of A*/
vector<vector<long long>> trans(vector<vector<long long>> A)
{
  vector<vector<long long>> B(A.at(0).size(), vector<long long>(A.size()));
  for(long i = 0; i < A.at(0).size(); i++){
    for(long j = 0; j < A.size(); j++) B.at(i).at(j) = A.at(j).at(i);
  }
  return B;
}

/*Computes transpose matrices of A*/
vector<vector<long double>> trans_double(vector<vector<long double>> A)
{
  vector<vector<long double>> B(A.at(0).size(), vector<long double>(A.size()));
  for(long i = 0; i < A.at(0).size(); i++){
    for(long j = 0; j < A.size(); j++) B.at(i).at(j) = A.at(j).at(i);
  }
  return B;
}

/*Computes determinant of matrices A*/
long double det(vector<vector<long double>> A){
  long n = A.size();
  long double det = 1.0, buf;
  for(long i = 0; i < n; i++){
    for(long j = 0; j < n; j++){
      if(i < j){
        buf = A.at(j).at(i) / A.at(i).at(i);
        for(long k = 0; k < n; k++) A.at(j).at(k) -= A.at(i).at(k) * buf;
      }
    }
  }
  for(long i = 0; i < n; i++) det *= A.at(i).at(i);
  return det;
}

/*Computes product of matrices and vectors*/
vector<long double> mul_mat_vec(vector<vector<long double>> A, vector<long double> x){
  vector<long double> v(A.size());
  for(long i = 0; i < A.size(); i++){
    for(long j = 0; j < x.size(); j++)v.at(i) += A.at(i).at(j) * x.at(j);
  }
  return v;
}

/*Computes product of two matriceses*/
vector<vector<long double>> mul_mat(vector<vector<long double>> A, vector<vector<long double>> B)
{
  long n = A.size(), m = B.at(0).size();
  vector<vector<long double>> C(n, vector<long double>(m));
  for(long i = 0; i < n; i++){
    for(long j = 0; j < m; j++){
      for(long k = 0; k < B.size(); k++) C.at(i).at(j) += A.at(i).at(k) * B.at(k).at(j);
    }
  }
  return C;
}

/*Preparation of Gaussian elimination over Z/2Z*/
bool have_1(vector<long long> v, long j){
  for(long i = 0; i < j; i++){
    if(v.at(i) == 1) return true;
  }
  return false;
}

/*Gaussian elimination over Z/2Z*/
vector<vector<long long>> gauss_mod2(vector<vector<long long>> A){
  long s;
  bool flag;
  for(long i = 0; i < A.size(); i++){
    for(long j = 0; j < A.at(0).size(); j++) A.at(i).at(j) = (A.at(i).at(j) % 2 + 2) % 2;
  }
  for(long j = 0; j < A.at(0).size(); j++){
    flag = false;
    for(long i = 0; i < A.size(); i++){
      if(A.at(i).at(j) == 1 && (not flag) && (not have_1(A.at(i), j))){
        flag = true;
        s = i;
      }
    }
    for(long i = 0; i < A.size(); i++){
      if(A.at(i).at(j) == 1 && flag && i != s){
        for(long k = 0; k < A.at(0).size(); k++) A.at(i).at(k) ^= A.at(s).at(k);
      }
    }
  }
  return A;
}

/*Gaussian elimination of augumented matrices over Z/2Z*/
tuple<vector<vector<long long>>, vector<vector<long long>>> aug_gauss_mod2(vector<vector<long long>> A, vector<vector<long long>> B){
  long s;
  bool flag;
  for(long i = 0; i < A.size(); i++){
    for(long j = 0; j < A.at(0).size(); j++){
      A.at(i).at(j) = (A.at(i).at(j) % 2 + 2) % 2;
      B.at(i).at(j) = (B.at(i).at(j) % 2 + 2) % 2;
    }
  }
  for(long j = 0; j < A.at(0).size(); j++){
    flag = false;
    for(long i = 0; i < A.size(); i++){
      if(A.at(i).at(j) == 1 && (not flag) && (not have_1(A.at(i), j))){
        flag = true;
        s = i;
      }
    }
    for(long i = 0; i < A.size(); i++){
      if(A.at(i).at(j) == 1 && flag && i != s){
        for(long k = 0; k < A.at(0).size(); k++){
          A.at(i).at(k) ^= A.at(s).at(k);
          B.at(i).at(k) ^= B.at(s).at(k);
        }
      }
    }
  }
  return forward_as_tuple(A, B);
}

/*Computes basises of kaernel of matrices A mod 2*/
vector<vector<long long>> kernel_basis(vector<vector<long long>> A){
  vector<vector<long long>> A_T, I, B, E, ker;
  A_T = trans(A);
  I = identity_mat_int(A_T.size());
  tie(B, E) = aug_gauss_mod2(A_T, I);
  for(long i = 0; i < B.size(); i++){
    if(is_zero(B.at(i)))ker.push_back(E.at(i));
  }
  return ker;
}

/*Computes all elements in kernel of A*/
vector<vector<long long>> kernel_mod2(vector<vector<long long>> A){
  vector<vector<long long>> ker_basis, v_all, ker;
  long long d;
  
  /*Computes all basis vectors of a kernel of A*/
  ker_basis = kernel_basis(A);
  vector<long long> tmp(ker_basis.at(0).size());
  d = ker_basis.size(); v_all = vector_all(d);
  
  /*Computes all linear combinations of basis vector of kernel of A modulo 2*/
  for(long i = 0; i < v_all.size(); i++){
    for(long j = 0; j < d; j++){
      for(long k = 0; k < ker_basis.at(0).size(); k++){
        tmp.at(k) ^= v_all.at(i).at(j) * ker_basis.at(j).at(k);
      }
    }
    ker.push_back(tmp);
    fill(tmp.begin(), tmp.end(), 0);
  }
  return ker;
}

/*Computes pseudo-inverse matrices of M*/
vector<vector<long double>> p_inv(vector<vector<long double>> M){
  vector<vector<long double>> X;
  X = mul_mat(trans_double(M), M);
  mat_ZZ A, A_inv;
  ZZ d;
  A.SetDims(X.size(), X.at(0).size());
  for(long i = 0; i < X.size(); i++){
    for(long j = 0; j < X.at(0).size(); j++)A[i][j] = to_ZZ((double)X.at(i).at(j));
  }
  inv(d, A_inv, A);
  for(long i = 0; i < X.size(); i++){
    for(long j = 0; j < X.at(0).size(); j++)X.at(i).at(j) = to_double(A_inv[i][j]) / to_double(d);
  }
  return mul_mat(X, trans_double(M));
}

/*Solve system of equation Mx = 0 for a vector x*/
vector<long double> solve_eq(vector<vector<long double>> M, vector<long double> x){
  vector<vector<long double>> A;
  A = p_inv(M);
  return mul_mat_vec(A, x);
}

/*Eratosthenes' sieve algorithm*/
vector<long long> esieve(long long M){
  vector<long long> v(M); v.at(0) = 1;
  long long s = 1, p = s + 1, t, k;
  while(p * p <= M){
    t = s;
    while(1){
      t += p;
      if(t >= M) break;
      v.at(t) = 1;
    }
    for(k = 1; k < M; k++){
      if(v.at(s + k) == 0) break;
    }
    s += k;
    p = s + 1;
  }
  vector<long long> p_vec {};
  vector<long long>::iterator it;
  it = p_vec.begin();
  for(long i = M - 1; i >= 0; i--){
    if(v.at(i) == 0) it = p_vec.insert(it, i + 1);
  }
  return p_vec;
}

/*Generates factor basis*/
vector<long long> factor_basis(long long n){
  long M = round((long double)n * log((long double)n));
  while(1){
    vector<long long> p_list = esieve(M);
    if(p_list.size() >= n){
      p_list.insert(p_list.begin(), -1);
      return p_list;
    }
    M++;
  }
}

/*Gram_Schmidt's orthogonalization algorithm*/
tuple<vector<vector<long double>>, vector<vector<long double>>> Gram_Schmidt(vector<vector<long double>> b)
{
  long n = b.size(), m = b.at(0).size();
  vector<vector<long double>> GSOb(n, vector<long double>(m));
  vector<vector<long double>> mu = identity_mat(n);
  for(long i = 0; i < n; i++){
    for(long j = 0; j < m; j++) GSOb.at(i).at(j) = b.at(i).at(j);
    for(long j = 0; j < i; j++){
      mu.at(i).at(j) = dot(b.at(i), GSOb.at(j)) / dot(GSOb.at(j), GSOb.at(j));
      for(long k = 0; k < m; k++) GSOb.at(i).at(k) -= mu.at(i).at(j) * GSOb.at(j).at(k);
    }
  }
  return forward_as_tuple(GSOb, mu);
}

/*Generates a square norm vector of GSOb*/
vector<long double> make_square_norm(vector<vector<long double>> GSOb){
  vector<long double> B(GSOb.size());
  for(long i = 0; i < GSOb.size(); i++) B.at(i) = dot(GSOb.at(i), GSOb.at(i));
  return B;
}

/*Computes a volume of a lattice b*/
long double vol(vector<vector<long double>> b){
  long double v = 1.0;
  vector<vector<long double>> GSOb, mu;
  tie(GSOb, mu) = Gram_Schmidt(b);
  for(long i = 0; i < b.size(); i++) v *= sqrt(dot(GSOb.at(i), GSOb.at(i)));
  return v;
}

/*Computes partial size reduced basis*/
tuple<vector<vector<long double>>, vector<vector<long double>>> size_reduce(vector<vector<long double>> b, vector<vector<long double>> mu, long i, long j){
  long double q;
  long m = b.at(0).size();
  if((fabs(mu.at(i).at(j)) > 0.5) && (fabs(fabs(mu.at(i).at(j)) - 0.5) >= EPSILON)){
    q = round(mu.at(i).at(j));
    for(long k = 0; k < m; k++) b.at(i).at(k) -= q * b.at(j).at(k);
    for(long k = 0; k <= j; k++){
      mu.at(i).at(k) -= q * mu.at(j).at(k);
      mu.at(j).at(j) = 1.0;
    }
  }
  return forward_as_tuple(b, mu);
}

/*GSO update*/
tuple<vector<vector<long double>>, vector<long double>> GSOupdate(vector<vector<long double>> mu, vector<long double> B, long k){
  long double nu = mu.at(k).at(k - 1), t;
  long double BB = B.at(k) + nu * nu * B.at(k - 1);
  long n = mu.size();
  mu.at(k).at(k - 1) = nu * B.at(k - 1) / BB;
  B.at(k) *= B.at(k - 1) / BB;
  B.at(k - 1) = BB;
  for(long j = 0; j <= k - 2; j++){
    t = mu.at(k - 1).at(j);
    mu.at(k - 1).at(j) = mu.at(k).at(j);
    mu.at(k).at(j) = t;
  }
  for(long i = k + 1; i < n; i++){
    t = mu.at(i).at(k);
    mu.at(i).at(k) = mu.at(i).at(k - 1) - nu * t;
    mu.at(i).at(k - 1) = t + mu.at(k).at(k - 1) * mu.at(i).at(k);
  }
  return forward_as_tuple(mu, B);
}

/*Computes LLL redued basis*/
vector<vector<long double>> LLL(vector<vector<long double>> b, long double d = 0.75){
  long n = b.size();
  vector<vector<long double>> GSOb, mu; tie(GSOb, mu) = Gram_Schmidt(b);
  vector<long double> B = make_square_norm(GSOb);
  long k = 1;
  while(k < n){
    for(long j = k - 1; j >= 0; j--) tie(b, mu) = size_reduce(b, mu, k, j);
    if((B.at(k) > (d - mu.at(k).at(k - 1) * mu.at(k).at(k - 1)) * B.at(k - 1)) || (fabs(B.at(k) - (d - mu.at(k).at(k - 1) * mu.at(k).at(k - 1)) * B.at(k - 1)) < EPSILON)){
      k++;
    }else{
      vector<long double> v = b.at(k - 1);
      b.at(k - 1) = b.at(k);
      b.at(k) = v;
      tie(mu, B) = GSOupdate(mu, B, k);
      k = max(k - 1, 1);
    }
  }
  return b;
}

/*LLL reduce (using NTL)*/
vector<vector<long double>> LLL_reduce(vector<vector<long double>> b, long double delta){
  long n = b.size(), m = b.at(0).size();
  long v = vol(b) * vol(b);
  ZZ det2 = to_ZZ(v);
  Mat<ZZ> c; c.SetDims(n, m);
  for(long i = 0; i < n; i++){
    for(long j = 0; j < m; j++) c[i][j] = to_ZZ((double)b.at(i).at(j));
  }
  LLL_FP(c, delta);
  for(long i = 0; i < n; i++){
    for(long j = 0; j < m; j++) b.at(i).at(j) = to_double(c[i][j]);
  }
  return b;
}

/*BKZ reduce (using NTL)*/
vector<vector<long double>> BKZ_reduce(vector<vector<long double>> b, long beta, long double delta){
  long n = b.size(), m = b.at(0).size();
  long v = vol(b) * vol(b);
  ZZ det2 = to_ZZ(v);
  Mat<ZZ> c; c.SetDims(n, m);
  for(long i = 0; i < n; i++){
    for(long j = 0; j < m; j++) c[i][j] = to_ZZ((double)b.at(i).at(j));
  }
  BKZ_FP(c, delta, beta);
  for(long i = 0; i < n; i++){
    for(long j = 0; j < m; j++) b.at(i).at(j) = to_double(c[i][j]);
  }
  return b;
}

/*Converts coefficient vector to lattice vector*/
vector<long double> coef2lat(vector<long double> v, vector<vector<long double>> b){
  vector<long double> x(b.at(0).size());
  for(long i = 0; i < b.size(); i++){
    for(long j = 0; j < b.at(0).size(); j++) x.at(j) += v.at(i) * b.at(i).at(j);
  }
  return x;
}

/*Converts lattice vector to lattice vector*/
vector<long double> lat2coef(vector<long double> v, vector<vector<long double>> b){
  long n = b.size();
  vector<long double> x(n);
  for(long i = 0; i < n; i++){
    x.at(i) = v.at(i) / b.at(i).at(i);
  }
  return x;
}

/*Enumerates closest vectors*/
vector<long double> ENUM_CVP(vector<vector<long double>> mu, vector<long double> B, double R, vector<long double> a){
  long n = B.size();
  vector<vector<long double>> sigma(n + 2, vector<long double>(n + 1));
  vector<long long> r(n + 1), w(n + 1);
  vector<long double> c(n + 1), rho(n + 2), v(n + 1), x(n), u;
  long k;
  v[1] = 1;
  for(long i = 0; i <= n; i++) r.at(i) = i;
  for(k = n; k >= 1; k--){
    for(long i = n; i >= k + 1; i--)sigma.at(i).at(k) = sigma.at(i + 1).at(k) + (a.at(i - 1) - v.at(i)) * mu.at(i - 1).at(k - 1);
    c.at(k) = a.at(k - 1) + sigma.at(k + 1).at(k);
    v.at(k) = round(c.at(k));
    w.at(k) = 1;
    rho.at(k) = rho.at(k + 1) + (c.at(k) - v.at(k)) * (c.at(k) - v.at(k)) * B.at(k - 1);
  }
  k = 1;
  while(1){
    rho.at(k) = rho.at(k + 1) + (c.at(k) - v.at(k)) * (c.at(k) - v.at(k)) * B.at(k - 1);
    if(rho.at(k) < R * R || fabs(rho.at(k) - R * R) < EPSILON){
      if(k == 1){
        for(long i = 1; i <= n; i++)x.at(i - 1) = v.at(i);
        return x;
      }
      k--;
      r.at(k - 1) = max(r.at(k - 1), r.at(k));
      for(long i = r.at(k); i >= k + 1; i--){
        sigma.at(i).at(k) = sigma.at(i + 1).at(k) + (a.at(i - 1) - v.at(i)) * mu.at(i - 1).at(k - 1);
      }
      c.at(k) = a.at(k - 1) + sigma.at(k + 1).at(k);
      v.at(k) = round(c.at(k));
      w.at(k) = 1;
    }else{
      k++;
      if(k == n + 1){
        x.clear();
        return x;
      }
      r.at(k - 1) = k;
      if(v.at(k) > c.at(k)){
        v.at(k) -= w.at(k);
      }else{
        v.at(k) += w.at(k);
      }
      w.at(k)++;
    }
  }
}

/*Generates coefficient vectors of close vectors*/
vector<vector<long double>> ENUM_CVP_all(vector<vector<long double>> mu, vector<long double> B, double R, vector<long double> a, double eps){
  long n = B.size();
  vector<vector<long double>> CVP_list;
  vector<long double> ENUM_CVP_v(n), pre_ENUM_CVP_v(n), x(n), v;
  while(1){
    for(long i = 0; i < n; i++) pre_ENUM_CVP_v.at(i) = ENUM_CVP_v.at(i);
    ENUM_CVP_v = ENUM_CVP(mu, B, R, a);
    if(ENUM_CVP_v.empty()) return CVP_list;
    if(find(CVP_list.begin(), CVP_list.end(), ENUM_CVP_v) == CVP_list.end())CVP_list.push_back(ENUM_CVP_v);
    R *= eps;
  }
}

/******************************************
準備終わり
*******************************************/

/******************************************
主要なsub-routine始め
*******************************************/

/*整数vector 𝑨の指数行列を作成する*/
vector<vector<long long>> exp_mat(vector<long long> A, vector<long long> p, long long n){
  vector<vector<long long>> E(n + 1, vector<long long>(A.size()));
  for(long i = 0; i < E.size(); i++){
    for(long j = 0; j < E.at(0).size(); j++){
      E.at(i).at(j) = div_count(A.at(j), p.at(i));
      if(i == n)E.at(i).at(j) = div_count(A.at(j), -1);
    }
  }
  return E;
}


/*Generates prime basis matrices*/
vector<vector<long double>> prime_mat(long long N, long double c, vector<long long> p, long n){
  vector<vector<long double>> A(2 * n, vector<long double>(n + 1));
  vector<long double> diag(n);//対角成分
  mt19937_64 get_rand_mt;

  for(long i = 0; i < n; i++) diag.at(i) = ceil((i + 1) / 2.0);

  random_device seed_gen;
  mt19937 engine(seed_gen());
  shuffle(diag.begin(), diag.end(), engine); //乱数的に並び替え
  for(long i = 0; i < n; i++){A.at(i).at(i) = diag.at(i); A.at(i + n).at(i) = diag.at(i);}

  for(long i = 0; i < n; i++){
    A.at(i).at(n) = round(pow(10, c) * log((long double)p.at(i + 1)));
    A.at(i + n).at(n) = round(pow(10, c) * log((long double)p.at(n - i - 1)));
  }

  return A;
}

vector<vector<vector<long double>>> ext_mat(vector<vector<long double>> B){
  long m = B.at(0).size(), n;
  n = m - 1;
  vector<vector<vector<long double>>> C(n, vector<vector<long double>>(n, vector<long double>(n + 1)));
  
  for(long k = 0; k < n; k++){
    for(long i = 0; i < n; i++){
      for(long j = 0; j < n + 1; j++){
        if(j == n){
          C.at(k).at(i).at(j) = B.at(i + k).at(j);
        }else{
          C.at(k).at(i).at(((j - k) + n + 1) % (n + 1)) = B.at(i + k).at(j);
        }
      }
    }
  }
  return C;
}

/*Generates a target vector*/
vector<long double> target(long long N, long double c, long n){
  vector<long double> t(n + 1);

  #if original == 0
    t.at(n) = round(pow(10, c) * log(N));
  #endif
  #if original == 1
    t.at(n) = pow(N, c) * log(N);
  #endif
  #if original == 2
    t.at(n) = pow(10, c) * log(N);
  #endif

  return t;
}

/*Tests whether pairs of integer (u, v) are sr-pair or not*/
bool sr_test(long long u, long long v, long long N, vector<long long> p, long long n){
  long long a = fabs(u - N * v);
  for(long i = 0; i < n; i++){
    while(a % p.at(i + 1) == 0){
      a /= p.at(i + 1);
    }
  }
  if(a == 1) return true;
  return false;
}

/******************************************
主要なsub-routine終わり
*******************************************/


int main(int argc, char **argv){
  bool loop_end, flag = false;
  long long a, e1, e2, num = 0, r, K, kernel_size = 0, SHIFT;
  unsigned long long u, v, U, N, q, qq, n, X, Y;
  long double c, alpha, TIME, cc, beta;
  vector<vector<vector<long double>>> S, SS, b, L, GSOb(SHIFT, vector<vector<long double>>(0, vector<long double>(0))), mu(SHIFT, vector<vector<long double>>(0, vector<long double>(0))), close_vecs(SHIFT, vector<vector<long double>>(0, vector<long double>(0)));
  vector<vector<long double>> L_T, b_T, B(SHIFT, vector<long double>(0));
  vector<long double> t, tt, w, x, y, tmp, e;
  if(argc == 7){
    N = (unsigned long long)atof(argv[1]);
    c = atof(argv[2]); cc = c;
    alpha = atof(argv[3]);
    beta = atof(argv[4]);
    TIME = atof(argv[5]);
    SHIFT = atof(argv[6]);
  }else{
    cout << "Errors: Comandline argument is incorrect." << endl;
    return 0;
  }
  n = round(pow(log2(N), alpha) / (alpha * log2(log2(N))));
  K = round(beta * pow(n, beta));
  cout << "Dimension = " << n << endl;
  cout << "Number of primes = " << K << endl;
  vector<long long> vec(K), ee(K), T, p;
  vector<vector<long long>> kernel, A, A_T, prekernel;
  w.resize(n);
  p = factor_basis(K);
  clock_t start = clock();
  time_t start_alg; time(&start_alg);
  time_t start_alg_c; time(&start_alg_c);
  while(1){
    /*After TIME [secs] have elapsed since the discovery of the sr-pairs, updates dimension*/
    if(time(NULL) - start_alg >= TIME || flag){
      flag = false;
      n++;
      K = round(beta * pow(n, beta));
      p = factor_basis(K);
      vec.resize(K);
      ee.resize(K);
      w.resize(n);
      for(long i = 0; i < A.size(); i++)A.at(i).resize(K);
      c = cc;
      time(&start_alg); time(&start_alg_c);
    }
    if(29 >= TIME && TIME >= 8){
      if(time(NULL) - start_alg_c >= TIME / 4){c += 0.03; time(&start_alg_c);}
    }else if(8 > TIME){
      if(time(NULL) - start_alg_c >= TIME / 2){c += 0.03; time(&start_alg_c);}
    }else if(TIME >= 30){
      if(time(NULL) - start_alg_c >= 10){c += 0.03; time(&start_alg_c);}
    }

    /*Composition of CVP*/
    S = ext_mat(prime_mat(N, c, p, n));
    t = target(N, c, n);//target vector

    /*For every lattice, find approximately closest vectors to a target vector t.*/
    for(long g = 0; g < SHIFT; g++){
      fill(w.begin(), w.end(), 0);//initialize a vector w

      /*Reduces lattice basis*/
      b.at(g) = L.at(g) = S.at(g);
      if(n >= 40){b.at(g) = BKZ_reduce(b.at(g), 20, 0.99);}else{b.at(g) = LLL(b.at(g), 0.99);}
      tie(GSOb.at(g), mu.at(g)) = Gram_Schmidt(b.at(g));
      B.at(g) = make_square_norm(GSOb.at(g));

      /*Composition of coefficient vectors of target vectors*/
      /*as target vectors of original lattice L*/
      loop_end = false;
      for(long i = 0; i <  L.size(); i++){
        if(L.at(g).at(i).at(i) == 1){
          w.at(i) = t.back() / L.at(g).at(i).back();
          break;
        }
      }

      /*as target vectors of basis reduced lattice b*/
      L_T = trans_double(L.at(g));
      tt = mul_mat_vec(L_T, w);
      b_T = trans_double(b.at(g));
      w = solve_eq(b_T, tt);

      /*Computes approximate solutions of CVP*/
      close_vecs.at(g) = ENUM_CVP_all(mu.at(g), B.at(g), 100, w, 0.99);
    }

    /*For every lattice matrices, find candidate pairs of integers for integer factorization.*/
    for(long g = 0; g < SHIFT; g++){
      for(long i = 0; i < close_vecs.at(g).size(); i++){
        x = close_vecs.at(g).at(i);
        tmp = coef2lat(x, b.at(g));
        e = lat2coef(tmp, L.at(g));

        /*Composition of a pair of integers (u, v) correspond to close_vec[i]*/
        u = v = 1;
        for(long j = 0; j < n; j++){
          if(e.at(j) > 0)u *= pow(p.at(j + 1), (long long)e.at(j));
          if(e.at(j) < 0)v *= pow(p.at(j + 1), -(long long)e.at(j));
        }

        if(u != 0 && v != 0 && sr_test(u, v, N, p, K)){//if (u, v) is sr-pair, then
          fill(vec.begin(), vec.end(), 0);
          for(long j = 0; j < K; j++){
            e1 = 0; e2 = 0;
            U = u;
          
            /*Computes exponents of p_{j} that is prime factor of u*/
            if(j < n)while(U % p.at(j + 1) == 0){e1++; U /= p.at(j + 1);}
          
            /*Computes exponents of p_{j} that is prime factor of a=u-vN*/
            a = u - v * N;
            while(a % p.at(j + 1) == 0 && a != 0){e2++; a /= p.at(j + 1);}
            vec.at(j) = e2 - e1;//exponential vector of each prime factor
          }
        
          /*Appends exponential vector of each prime factor to A*/
          if(find(A.begin(), A.end(), vec) == A.end()){
            num++;
            cout << num << " sr-pairs were found. (" << K << "pairs are needed)\n";
            A.push_back(vec);//Append to A
            time(&start_alg); time(&start_alg_c);//時間の数えなおし
          }
        }
      }

      /*prime factorization*/
      if(A.size() != 0 && A.size() >= A.at(0).size()){
        fill(ee.begin(), ee.end(), 0);
        A_T = trans(A);
        kernel = kernel_mod2(A_T);
        for(long i = 0; i < prekernel.size(); i++)prekernel.at(i).resize(kernel.at(0).size());

        if(kernel != prekernel){
        /*Verifies any t = T that satisfies\sum_{j=1}^{ee.size()}t_j(e_{i.j}-e'_{i.j})=0(mod\ 2)*/
          for(long i = 0; i < kernel.size(); i++){
            if(find(prekernel.begin(), prekernel.end(), kernel.at(i)) == prekernel.end()){
              cout << i + 1 << " / " << kernel.size() << endl;
              T = kernel.at(i);
              for(long j = 0; j < T.size(); j++){
                if(T.at(j) == 1){
                  for(long k = 0; k < K; k++)ee.at(k) += A.at(j).at(k);
                }
              }
              X = Y = 1;
              for(long j = 0; j < K; j++){
                if(ee.at(j) > 0){X = mulm(X, powm(p.at(j + 1), ee.at(j) / 2, N), N); X %= N;}
                if(ee.at(j) < 0){Y = mulm(Y, powm(p.at(j + 1), -ee.at(j) / 2, N), N); Y %= N;}
              }
              cout << "X = " << X << ", Y = " << Y << endl;
              if(X != Y){q = gcd(X - Y, N); qq = gcd(X + Y, N);}

              if(1 < q && q < N){
                clock_t end = clock();
                cout << "c = " << cc << " -> " << c << ", alpha = " << alpha << ", beta = " << beta << endl;
                cout << N << " = " << q << " * " << N / q << endl;
                cout << "Run time = " << (long double)(end - start) / CLOCKS_PER_SEC << "[secs]" << endl;
                return 0;
              }
        
              if(1 < qq && qq < N){
                clock_t end = clock();
                cout << "c = " <<  cc << " -> " << c << ", alpha = " << alpha << ", beta = " << beta << endl;
                cout << N << " = " << qq << " * " << N / qq << endl;
                cout << "Run time = " << (long double)(end - start) / CLOCKS_PER_SEC << "[secs]" << endl;
                return 0;
              }
              fill(ee.begin(), ee.end(), 0);
            }
          }
          prekernel = kernel;
          flag = true;
        }
      }
    }
  }
}

/*
A list of results
----11 bits---------------------
c = 4 -> 4, alpha = 1.5, beta = 1.1
1763 = 43 * 41
Run time = 0.251809[secs]

----12 bits---------------------
c = 4 -> 4, alpha = 1.5, beta = 1.1
3599 = 61 * 59
Run time = 0.017995[secs]

----13 bits---------------------
c = 4 -> 4, alpha = 1.5, beta = 1.1
7811 = 73 * 107
Run time = 0.08957[secs]

----14 bits---------------------
c = 4 -> 4, alpha = 1.5, beta = 1.1
9047 = 83 * 109
Run time = 0.061336[secs]

----15 bits---------------------
c = 5.4 -> 5.4, alpha = 1.8, beta = 1.1
28459 = 149 * 191
Run time = 0.467635[secs]

method changed result
c = 4 -> 4, alpha = 1.5, beta = 1.1
28459 = 191 * 149
Run time = 0.259203[secs]

----16 bits---------------------
c = 4 -> 4, alpha = 1.5, beta = 1.1
49163 = 233 * 211
Run time = 0.218718[secs]

----17 bits---------------------
c = 4 -> 4, alpha = 1.5, beta = 1.1
108539 = 311 * 349
Run time = 0.220647[secs]

----18 bits---------------------
206711 = 421 * 491
Run time = 0.659499[secs]

----19 bits---------------------
c = 4 -> 4, alpha = 1.5, beta = 1.1
342029 = 599 * 571
Run time = 1.3286[secs]

----20 bits---------------------
c = 5.4 -> 5.4, alpha = 1.7, beta = 1.1
695683 = 757 * 919
Run time = 1.80727[secs]

c = 4 -> 4, alpha = 1.5, beta = 1.1
695683 = 919 * 757
Run time = 0.871962[secs]

----21 bits---------------------
c = 5.4 -> 5.4, alpha = 1.7, beta = 1.1
1970551 = 359 * 5489
Run time = 2.58173[secs]

c = 4 -> 4, alpha = 1.5, beta = 1.1
1970551 = 3949 * 499
Run time = 1.3365[secs]

----22 bits---------------------
c = 5.4 -> 5.4, alpha = 1.7, beta = 1.1
2381663 = 1427 * 1669
Run time = 1.98352[secs]

c = 4 -> 4, alpha = 1.5, beta = 1.1
2381663 = 1669 * 1427
Run time = 1.95376[secs]

----23 bits---------------------
c = 5.4 -> 5.4, alpha = 1.7, beta = 1.1
8350847 = 2617 * 3191
Run time = 4.68485[secs]

c = 4 -> 4, alpha = 1.5, beta = 1.1
8350847 = 2617 * 3191
Run time = 3.00512[secs]

----24 bits---------------------
c = 5.4 -> 5.4, alpha = 1.7, beta = 1.1
10586183 = 3931 * 2693
Run time = 5.98235[secs]

c = 4 -> 4, alpha = 1.5, beta = 1.1
10586183 = 3931 * 2693
Run time = 5.76747[secs]

----25 bits---------------------
c = 5.4 -> 5.4, alpha = 1.7, beta = 1.1
20903063 = 4561 * 4583
Run time = 8.71293[secs]

c = 4 -> 4.03, alpha = 1.5, beta = 1.1
20903063 = 4583 * 4561
Run time = 7.1407[secs]

----26 bits---------------------
c = 5.4 -> 5.4, alpha = 1.7, beta = 1.1
45021971 = 1957477 * 23
Run time = 13.1473[secs]

c = 4 -> 4, alpha = 1.5, beta = 1.1
45021971 = 1957477 * 23
Run time = 4.48952[secs]

----27 bits---------------------
c = 5.4 -> 5.43, alpha = 1.7, beta = 1.1
90476623 = 9551 * 9473
Run time = 16.1835[secs]

c = 4 -> 4.03, alpha = 1.5, beta = 1.1
90476623 = 9551 * 9473
Run time = 17.5321[secs]

----28 bits---------------------
c = 5.4 -> 5.55, alpha = 1.7, beta = 1.1
153713251 = 12953 * 11867
Run time = 31.1185[secs]

c = 4 -> 4.03, alpha = 1.5, beta = 1.1
153713251 = 11867 * 12953
Run time = 30.2523[secs]

----29 bits---------------------
c = 5.4 -> 5.46, alpha = 1.7, beta = 1.1
269419387 = 16417 * 16411
Run time = 35.068[secs]

c = 4 -> 4, alpha = 1.5, beta = 1.1
269419387 = 16417 * 16411
Run time = 34.935[secs]

----30 bits---------------------
c = 5.4 -> 5.43, alpha = 1.7, beta = 1.1
586016227 = 4614301 * 127
Run time = 44.536[secs]

c = 4 -> 4.03, alpha = 1.5, beta = 1.4
586016227 = 35687 * 16421
Run time = 32.4403[secs]

----31 bits---------------------
c = 5.4 -> 5.4, alpha = 1.7, beta = 1.1
2047104727 = 38713 * 52879
Run time = 10.2485[secs]

c = 4 -> 4.09, alpha = 1.5, beta = 1.4
2047104727 = 38713 * 52879
Run time = 71.2341[secs]

----32 bits---------------------
c = 5.4 -> 5.4, alpha = 1.7, beta = 1.1
2529975709 = 42073 * 60133
Run time = 9.95274[secs]

c = 4 -> 4.21, alpha = 1.5, beta = 1.5
2047104727 = 38713 * 52879
Run time = 67.4204[secs]

----33 bits---------------------
c = 5.4 -> 5.4, alpha = 1.7, beta = 1.1
5726447347 = 76001 * 75347
Run time = 17.8248[secs]

c = 4 -> 4.15, alpha = 1.5, beta = 1.45
5726447347 = 76001 * 75347
Run time = 112.346[secs]

----34 bits---------------------
c = 5.4 -> 5.4, alpha = 1.7, beta = 1.1
12604829417 = 96337 * 130841
Run time = 29.7365[secs]

c = 5.4 -> 5.52, alpha = 1.5, beta = 1.1
12604829417 = 96337 * 130841
Run time = 97.6522[secs]

----35 bits---------------------
c = 5.4 -> 5.4, alpha = 1.7, beta = 1.1
32465228347 = 199603 * 162649
Run time = 38.2366[secs]

c = 3.7 -> 3.82, alpha = 1.64, beta = 1.1
32465228347 = 31 * 1047265430
Run time = 52.1529[secs]

----36 bits---------------------
c = 5.4 -> 5.4, alpha = 1.7, beta = 1.1
40949546243 = 244313 * 167611
Run time = 50.6907[secs]

c = 3.7 -> 3.94, alpha = 1.63, beta = 1.1
40949546243 = 31 * 1320953104
Run time = 70.2147[secs]

----37 bits---------------------
c = 3.7 -> 3.7, alpha = 1.615, beta = 1.1
114391859609 = 343073 * 333433
Run time = 127.231[secs]

----38 bits---------------------
c = 3.7 -> 3.76, alpha = 1.61, beta = 1.1
222974630273 = 492421 * 452813
Run time = 199.449[secs]

----39 bits---------------------
c = 3.97 -> 4.87, alpha = 1.6, beta = 1.1
516842741119 = 932177 * 554447
Run time = 167.662[secs]

----40 bits---------------------
c = 4 -> 4.18, alpha = 1.58, beta = 1.1
769319877019 = 884293 * 869983
Run time = 233.661[secs]

----41 bits---------------------
c = 3.6 -> 4.32, alpha = 1.57, beta = 1.1
1436884384709 = 1707067 * 841727
Run time = 430.634[secs]

----42 bits---------------------
c = 3.7 -> 3.73, alpha = 1.57, beta = 1.05
2267839867523 = 1370933 * 1654231
Run time = 887.263[secs]

----50 bits---------------------
708343001672333
*/