using PyCall
using KSVD

# Use Python3
ENV["PYTHON"] = strip(readall(`which python3`))

# Load digits dataset from scikt-learn
@pyimport sklearn.datasets as datasets
digits = datasets.load_digits()

# Each column of Y is a flattened character image
Y = digits["data"]'

D, X = ksvd(Y, 256, max_iter_mp = 20, tolerance_zeros = 1.0)

# Find D and X such that Y ≈ DX
assert(norm(Y - D * X))

#
println("sparsity: $(sum(X .== 0) / length(X))")

writedlm("digits_D.dlm", D)
