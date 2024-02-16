from fpylll import *
from time import perf_counter
from datetime import datetime
from random import seed, shuffle
import numpy as np

#=============================
# Generates an n-bits semiprime
#=============================
def gen_semiprime(n):
	while 1:
		p = next_prime(randint(2 ^ ((n - 1) // 2), 2 ^ ((n + 1) // 2)))
		q = next_prime(randint(2 ^ ((n - 3) // 2), 2 ^ ((n + 1) // 2)))
		N = p * q
		if N.nbits() == n and p != q: return N

#===================
# Input an integer N
#===================
N = ZZ(input("N = "))
if N < 100: N = gen_semiprime(N)
R = Integers(N)

#==========================
# Set a rank of the lattice
#==========================
bit_size = N.nbits()
if bit_size >= 40:
	n = round(2.2 * bit_size / log(bit_size) - 8); m = n + 1
else:
	n = 15; m = 16
J = n * n
K = J + J
J = (K + K) // 3

#=================================
# Generate a list of prime numbers
#=================================
Prime = Primes()
P = np.array([-1] + prime_range(Prime.unrank(K)))
Q = P; P = vector(ZZ, P)

#================================
# Construction of a target vector
#================================
t = vector(ZZ, m); t[n] = ceil(100000 * log(N))

B_right = np.array([np.ceil(100000 * np.log(P[1 : n + 1]))]).T
S = np.ceil(np.arange(1, m) * 0.6)
#S = np.ones(n); S[0 : n // 8] = 4; S[n // 8 : n // 4] = 3; S[n // 4 : n // 2] = 2
#S = np.ceil(np.sqrt(10 * np.log(np.arange(1, m))))
List = []; num = 0; loop_times = 0; already_list = []; f: bool = False; r = 0

start = perf_counter()
#=============================
# Main factorization algorithm
#=============================
for _ in range(factorial(n) // 2):
	loop_times += 1
	# Construction of lattice basis
	seed(datetime.now().timestamp()); shuffle(S)
	C = Matrix(ZZ, np.hstack([np.diag(S), B_right]))
	B = IntegerMatrix.from_matrix(C)

	# Lattice basis reduction
	LLL.reduction(B)
	M = GSO.Mat(B); _ = M.update_gso()
	
	# Babai's algorithm 
	w = M.babai(t)
	# Solving approximate CVP
	radius = (vector(ZZ, B.multiply_left(w)) - t).norm()
	enum = Enumeration(M, strategy = EvaluatorStrategy.BEST_N_SOLUTIONS)
	solutions = enum.enumerate(0, n, radius * radius + 1, 0, M.from_canonical(t))
	for _, b in solutions:
		b = IntegerMatrix.from_iterable(1, B.nrows, map(lambda x: int(round(x)), b)); w = b * B
		e = C.solve_left(vector(w[0])); e = np.array(e)
		
		# Construction of (u, v)-pairs
		PositiveIndex = np.where(e > 0)
		NegativeIndex = np.where(e < 0)
		primes_for_u = Q[PositiveIndex[0] + 1]
		primes_for_v = Q[NegativeIndex[0] + 1]
		u_factor = primes_for_u ^ e[PositiveIndex[0]]
		v_factor = primes_for_v ^ (-e[NegativeIndex[0]])
		u = prod(vector(ZZ, u_factor))
		v = prod(vector(ZZ, v_factor))

		#T = u - v * N
		LT = np.array(factor(u - v * N))
		Lu = np.array(factor(u))
		L = LT.T; M = Lu.T
		if len(L) > 0:
			if set(L[0]) <= set(P):# Smoothness check
				vec = zero_vector(ZZ, K)
				for p in np.r_[L[0], M[0]]:
					e1 = e2 = 0
					j = np.where(P == p); j = j[0][0]
					if p in M[0]:
						index = np.where(Lu == p)
						e1 = Lu[index[0][0]][1]
					if p in L[0]:
						index = np.where(LT == p)
						e2 = LT[index[0][0]][1]
					vec[j - 1] = e2 - e1
				if vec not in List:
					num += 1
					print(num, "sr-pairs are found. u =", u,"v =", v)
					List.append(vec)
					AA = Matrix(GF(2), List)
					r = AA.rank()

	if len(List) > r:
		V = AA.kernel();
		for tt in V:
			if not tt.is_zero():
				OneIndex = np.where(tt); OneIndex = OneIndex[0]
				Array = np.array(List); Array = Array[OneIndex]
				ee = vector(ZZ, Array.sum(axis = 0))
				if ee not in already_list:
					already_list.append(ee)

					ee = np.array(ee) // 2
					PositiveIndex = np.where(ee > 0)
					NegativeIndex = np.where(ee < 0)
					#primes_for_X = Q[PositiveIndex[0] + 1]
					#primes_for_Y = Q[NegativeIndex[0] + 1]
					#X_factor = primes_for_X ^ ee[PositiveIndex[0]]
					#Y_factor = primes_for_Y ^ (-ee[NegativeIndex[0]])
					#X = prod(vector(R, X_factor))
					#Y = prod(vector(R, Y_factor))

					X = Y = R(1)
					for j in PositiveIndex[0]:
						X *= R(P[j + 1]) ^ ee[j]
					for j in NegativeIndex[0]:
						Y *= R(P[j + 1]) ^ (-ee[j])

					print("X =", X, "Y =", Y)
					if X != Y:
						p = gcd(ZZ(X - Y), N)
					else:
						p = gcd(ZZ(X + Y), N)
					if 1 < p < N:
						print("==============================\nN =", N, ", bit size = ", bit_size, "\nc = 7, beta = 2.0\nN = ", p, "*", N / p, "\nloop times = ", loop_times, "\nnumber of sr-pairs =", num, "\n==============================")
						f = True; break
				if f:
					break
		if f:
			break

end = perf_counter()
if f:
	print(end - start, "[secs]")
else:
	print("Factorization failed.")
