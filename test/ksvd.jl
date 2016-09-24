# basic one
Y = [
     0  3;
     1  0;
     1 -3;
    -2  3;
     0  0;
]

D, X = ksvd(Y, 20, tolerance_nonzero = -1)
@test norm(Y-D*X) < 1e-6

# max_iter must be > 0
@test_throws ArgumentError ksvd(Y, 2, max_iter = 0)
@test_throws ArgumentError ksvd(Y, 2, max_iter = -1)

# n_atoms must be larger than the signal dimensions (same as the dimensions of
# atoms) since K-SVD is an algorithm for designing overcomplete dictionaries
@test_throws ArgumentError ksvd(rand(4, 3), 2)

# But should work well when size(Y, 1) == n_atoms
Y = [
    -1 1 2;
     1 0 1
]
D, X = ksvd(Y, 2, max_iter_mp = 800)
@test norm(Y-D*X) < 0.001  # relax the constraint since the dictionary is small

# Return only if X is sparse enough
Y = [
	0  2  3 -1  1;
	1 -3  1  3  0
]

# size(X) == (K, N)
K = 5
N = size(Y, 2)

# more than 20% of elements must be zeros
tolerance_nonzero = 0.8 * K * N
D, X = ksvd(Y, K, max_iter = Int(1e10), tolerance_nonzero = tolerance_nonzero)
@test sum(X .!= 0) <= tolerance_nonzero
