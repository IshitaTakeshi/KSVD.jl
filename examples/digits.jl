# Use Python3
ENV["PYTHON"] = strip(readstring(`which python3`))

using PyCall
using KSVD


# Load digits dataset from scikt-learn
@pyimport sklearn.datasets as datasets
digits = datasets.load_digits()

# Each column of Y is a flattened character image
Y = digits["data"]'

# Derive D and X such that DX ≈ Y
D, X = ksvd(
    Y,
    # the number of atoms in D
    256,
    # max iterations of K-SVD
    max_iter = 200,
    # max iterations of matching pursuit called in K-SVD
    max_iter_mp = 40,
    # stop iteration if more than 96% of elements in X are zeros
    sparsity_allowance = 0.96
)

# Find D and X such that Y ≈ DX
println("||Y - D * X|| = $(norm(Y - D * X))")

println("The ratio of zero elemnts in the matrix X: ",
        sum(X .== 0) / length(X))

writedlm("digits256.dlm", D)
