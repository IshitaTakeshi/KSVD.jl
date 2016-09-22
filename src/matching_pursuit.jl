import Base: sparsevec
using DataStructures


default_max_iter = 20
tolerance = default_tolerance = 1e-16


function sparsevec(d::DefaultDict, m::Int)
    sparsevec(collect(keys(d)), collect(values(d)), m)
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
            tolerance = default_tolerance
        )
    end
    return X
end


function matching_pursuit(data::Vector, dictionary::Matrix;
                          max_iter::Int = default_max_iter,
                          tolerance::Float64 = default_tolerance)
    """
    Finds x such that y = Dx where y is `data` and D is `dictionary`.
    """

    assert(tolerance > 0)
    assert(max_iter > 0)
    # TODO check whether matrix dimensions match

    n, K = size(dictionary)

    # TODO test
    for i in 1:K
        assert(any(dictionary[:, i] .!= 0))
    end

    residual = copy(data)
    dictionary = copy(dictionary)

    xdict = DefaultDict(Int, Float64, 0)
    for i in 1:max_iter
        if norm(residual) < tolerance
            return sparsevec(xdict, size(dictionary, 2))
        end

        # find an atom with maximum inner product
        products = dictionary' * residual
        _, maxindex = findmax(abs(products))
        maxval = products[maxindex]
        atom = dictionary[:, maxindex]

        # c is the length of the projection of data onto atom
        c = maxval / sum(atom.^2)  # equivalent to maxval / norm(atom)^2
        residual -= atom * c

        xdict[maxindex] += c
    end
    return sparsevec(xdict, size(dictionary, 2))
end
