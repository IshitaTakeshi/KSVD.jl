# basic one
Y = [
     0  3;
     1  0;
     1 -3;
    -2  3;
     0  0;
]

D, X = ksvd(Y, 8, max_iter_mp = 600, sparsity_allowance = 1.0)
@test norm(Y-D*X) < 1e-4

# sparsity_allowance must be > 0
@test_throws ArgumentError ksvd(Y, 20, sparsity_allowance = -0.1)
@test_throws ArgumentError ksvd(Y, 20, sparsity_allowance = 1.1)

# should work normally
ksvd(Y, 20, max_iter = 1)
ksvd(Y, 20, sparsity_allowance = 0.0)
ksvd(Y, 20, sparsity_allowance = 1.0)

# But should work well when size(Y, 1) == n_atoms
Y = [
    -1 1 2;
     1 0 1
]
D, X = ksvd(Y, 2, max_iter_mp = 4000)
@test norm(Y-D*X) < 0.001  # relax the constraint since the dictionary is small

# Return only if X is sparse enough
Y = [
	0  2  3 -1  1;
	1 -3  1  3  0
]

# More than 20% of elements in X must be zeros
sparsity_allowance = 0.2
D, X = ksvd(Y, 5, max_iter = Int(1e10), sparsity_allowance = sparsity_allowance)
@test sum(iszero, X) / length(X) > sparsity_allowance
