from fpylll import *
from time import perf_counter
from datetime import datetime
import random
import numpy as np

#=============================
#Generates an n-bits semiprime
#=============================
def gen_semiprime(n):
	while 1:
		p = next_prime(randint(2 ^ ((n - 1) // 2), 2 ^ ((n + 1) // 2)))
		q = next_prime(randint(2 ^ ((n - 3) // 2), 2 ^ ((n + 1) // 2)))
		N = p * q
		if(N.nbits() == n):
			return N

def DivisionCount(N, p):
	if mod(N, p) != 0:
		return 0
	if mod(N, p) == 0:
		return 1 + DivisionCount(N // p, p)

N = ZZ(input("N = "))
if N < 100:
	N = gen_semiprime(N)
R = Integers(N)

start = perf_counter()

f = False
bit_size = N.nbits()
if bit_size >= 40:
	n = round(2.2 * bit_size / log(bit_size) - 8); m = n + 1
else:
	n = 15; m = 16
J = n * n
K = J + J
J = K * 0.666

#A list of prime numbers.
Prime = Primes()
P = np.array([-1] + prime_range(Prime.unrank(K)))
Q = P; P = vector(ZZ, P)
B_right = np.array([np.ceil(100000 * np.log(P[1 : n + 1]))]).T
S = np.ceil(np.arange(1, m) * 0.6)

#Construction of a target vector
t = vector(ZZ, m); t[n] = ceil(100000 * log(N))

List = []; r = 0; num = 0; loop_times = 0; already_list = []
while True:
	loop_times += 1
	# Construction of lattice basis
	random.seed(datetime.now().timestamp()); random.shuffle(S)
	C = Matrix(ZZ, np.hstack([np.diag(S), B_right]))
	B = IntegerMatrix.from_matrix(copy(C))

	# Lattice basis reduction
	LLL.reduction(B)
	M = GSO.Mat(B); _ = M.update_gso()
	
	# Babai's algorithm 
	w = M.babai(t)

	# Solving approximate CVP
	radius = (vector(ZZ, B.multiply_left(w)) - t).norm() ^ 2
	solutions = []
	enum = Enumeration(M, strategy = EvaluatorStrategy.BEST_N_SOLUTIONS, nr_solutions = 20)
	solutions = enum.enumerate(0, n, radius, 0, M.from_canonical(t))
	
	for a, b in solutions:
		b = IntegerMatrix.from_iterable(1, B.nrows, map(lambda x: int(round(x)), b)); w = b * B
		e = C.solve_left(vector(w[0]))
		e = np.array(e)
		
		# Construction of (u, v)-pairs
		PositiveIndex = np.where(e > 0)
		NegativeIndex = np.where(e < 0)
		primes_for_u = Q[PositiveIndex[0] + 1]
		primes_for_v = Q[NegativeIndex[0] + 1]
		u_factor = primes_for_u ^ e[PositiveIndex[0]]
		v_factor = primes_for_v ^ (-e[NegativeIndex[0]])
		u = prod(vector(ZZ, u_factor))
		v = prod(vector(ZZ, v_factor))
		
		T = u - v * N
		L = np.array(factor(T)).T
		if len(set(L[0]) - set(P)) == 0:# Smoothness check
			vec = vector(ZZ, K)
			for j in range(K):
				p = P[j + 1]
				e1 = DivisionCount(u, p)
				e2 = DivisionCount(T, p)
				vec[j] = e2 - e1
			if vec not in List:
				num += 1
				print(num, "sr-pairs are found. (K =", K, ")")
				List.append(vec)
				#AA = Matrix(GF(2), List)
				#r = AA.rank()

	if num >= J:
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
						primes_for_X = Q[PositiveIndex[0] + 1]
						primes_for_Y = Q[NegativeIndex[0] + 1]
						X_factor = (primes_for_X ^ ee[PositiveIndex[0]]) % N
						Y_factor = primes_for_Y ^ (-ee[NegativeIndex[0]]) % N
						X = prod(vector(R, X_factor))
						Y = prod(vector(R, Y_factor))
						#X = XY_product(P, ee, K - 1, R, positive = True)
						#Y = XY_product(P, ee, K - 1, R, positive = False)
						print("X =", X, ", Y =", Y)
						if X != Y:
							p = gcd(ZZ(X - Y), N)
						else:
							p = gcd(ZZ(X + Y), N)
						if 1 < p < N:
							print("==============================\nN = ", N, ", bit size = ", N.nbits(), "\nc = 5, beta = 2.0\nN = ", p, " * ", N / p, "\nloop times = ", loop_times, "\nnumber of sr-pairs = ", num, "\n==============================")
							f = True; break
					if f:
						break
			if f:
				break

end = perf_counter()
print(end - start, "[secs]")
