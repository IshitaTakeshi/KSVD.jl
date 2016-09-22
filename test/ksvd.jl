# basic one
Y = [
     0  3;
     1  0;
     1 -3;
    -2  3;
     0  0;
]
D, X = ksvd(Y, 20)
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
