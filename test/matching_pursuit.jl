D = [
	-5  -9  -9   7   0;
	 7   0   4   0   9;
	 8   0   5   4  -6
]


# when data is representad as a vector
Y = [
    -1;
     0;
     5
]

X = matching_pursuit(Y, D)
@test norm(Y-D*X) < 1e-6


# when data is representad as a matrix
Y = [
    -1  2;
     0  3;
     5 -5
]

X = matching_pursuit(Y, D)
@test norm(Y-D*X) < 1e-6


# when max_iter < 0
@test_throws ArgumentError matching_pursuit(Y, D; max_iter = 0)
@test_throws ArgumentError matching_pursuit(Y, D; max_iter = -1)


# when tolerance < 0
@test_throws ArgumentError matching_pursuit(Y, D, tolerance = 0.)
@test_throws ArgumentError matching_pursuit(Y, D; tolerance = -1.0)


# when dimensions of data and atoms do not match

Y = [
     2;
    -5
]

@test_throws ArgumentError matching_pursuit(Y, D)
