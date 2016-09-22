import Base: sparsevec
using DataStructures

# The implementation is referencing the wikipedia page
# https://en.wikipedia.org/wiki/Matching_pursuit#The_algorithm

default_max_iter = 20
default_tolerance = 1e-16


function sparsevec(d::DefaultDict, m::Int)
    sparsevec(collect(keys(d)), collect(values(d)), m)
end


function matching_pursuit_(data::Vector, dictionary::Matrix,
                           max_iter::Int, tolerance::Float64)
    n_atoms = size(dictionary, 2)

    residual = copy(data)

    xdict = DefaultDict(Int, Float64, 0)
    for i in 1:max_iter
        if norm(residual) < tolerance
            return sparsevec(xdict, n_atoms)
        end

        # find an atom with maximum inner product
        products = dictionary' * residual
        _, maxindex = findmax(abs(products))
        maxval = products[maxindex]
        atom = dictionary[:, maxindex]

        # c is the length of the projection of data onto atom
        a = maxval / sum(atom.^2)  # equivalent to maxval / norm(atom)^2
        residual -= atom * a

        xdict[maxindex] += a
    end
    return sparsevec(xdict, n_atoms)
end


function matching_pursuit(data::Vector, dictionary::Matrix;
                          max_iter::Int = default_max_iter,
                          tolerance::Float64 = default_tolerance)
    """
    Finds x such that Dx approximates y where y is `data` and D is `dictionary`.
    """

    if tolerance <= 0
        throw(ArgumentError("`tolerance` must be > 0"))
    end

    if max_iter <= 0
        throw(ArgumentError("`max_iter` must be > 0"))
    end

    if size(data, 1) != size(dictionary, 1)
        throw(ArgumentError(
            "Dimensions must match: `size(data, 1)` and `size(dictionary, 1)`."
        ))
    end

    if rank(dictionary) != size(data, 1)
        throw(ArgumentError(
            "The rank of the dictionary is too small. " *
            "`rank(dictionary)` must be equal to `size(data, 1)`."
        ))
    end

    matching_pursuit_(data, dictionary, max_iter, tolerance)
end


function matching_pursuit(data::Matrix, dictionary::Matrix;
                          max_iter::Int = default_max_iter,
                          tolerance::Float64 = default_tolerance)
    K = size(dictionary, 2)
    N = size(data, 2)

    X = spzeros(K, N)

    for i in 1:N
        X[:, i] = matching_pursuit(
            vec(data[:, i]),
            dictionary,
            max_iter = max_iter,
            tolerance = tolerance
        )
    end
    return X
end
