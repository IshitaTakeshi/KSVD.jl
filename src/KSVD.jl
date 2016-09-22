module KSVD

# This is an implementation of the K-SVD algorithm.
# The original paper:
# K-SVD: An Algorithm for Designing Overcomplete Dictionaries
# for Sparse Representation
# http://www.cs.technion.ac.il/~freddy/papers/120.pdf

# Variable names are based on the original paper.
# If you try to read the code, I recommend you to see Figure 2 first.
#

export ksvd, matching_pursuit

include("matching_pursuit.jl")

default_max_iter = 200

srand(1234)  # for stability of tests


function error_matrix(Y, D, X, k)
    Eₖ = Y
    for j in 1:size(D, 2)
        if j != k
            Eₖ -= D[:, j] * X[j, :]
        end
    end
    return Eₖ
end


function init_dictionary(n::Int, K::Int)
    # D must be a rank-n matrix
    assert(n <= K)
    D = rand(n, K)
    while rank(D) != n
        D = rand(n, K)
    end

    for k in 1:K
        D[:, k] /= norm(D[:, k])
    end
    return D
end


function ksvd(Y::AbstractArray, D::AbstractArray, X::AbstractArray)
    N = size(Y, 2)
    for k in 1:size(X, 1)
        xₖ = X[k, :]

        # ignore if the k-th row is zeros
        if all(xₖ .== 0)
            continue
        end

        # wₖ is the indices where xₖ is non-zero,
        # which is equivalent to [i for i in N if xₖ[i] != 0]
        _, wₖ, _ = findnz(xₖ)

        # Eₖ * Ωₖ implies a selection of error columns that
        # correspond to examples that use the atom D[:, k]
        Eₖ = error_matrix(Y, D, X, k)
        Ωₖ = sparse(wₖ, 1:length(wₖ), ones(length(wₖ)), N, length(wₖ))
        # Note that S is a vector that contains diagonal elements of Δ
        # such that Eₖ * Ωₖ == U * Δ * V .
        # Non-zero entries of X are set to
        # the first column of V multiplied by Δ(1, 1)
        U, S, V = svd(Eₖ * Ωₖ)
        D[:, k] = U[:, 1]
        X[k, wₖ] = V[:, 1] * S[1]
    end
    return D, X
end


# TODO add tolerance of D as an argument
function ksvd(Y::Matrix, n_atoms::Int; max_iter::Int = default_max_iter)
    """
    K-SVD derives the most efficient dictionary D
    """

    K = n_atoms
    n, N = size(Y)

    if K < n
        throw(ArgumentError("size(Y, 1) must be >= K"))
    end
    if max_iter <= 0
        throw(ArgumentError("`max_iter` must be > 0"))
    end

    # D is a dictionary matrix that contains atoms for columns.
    D = init_dictionary(n, K)  # size(D) == (n, K)

    X = spzeros(K, N)  # just for making X global in this function
    for i in 1:max_iter
        X = matching_pursuit(Y, D, max_iter = 200)
        D, X = ksvd(Y, D, X)
    end
    return D, X
end

end # module
