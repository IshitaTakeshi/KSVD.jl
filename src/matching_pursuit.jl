using DataStructures

# The implementation is referencing the wikipedia page
# https://en.wikipedia.org/wiki/Matching_pursuit#The_algorithm

const default_max_iter = 20
const default_tolerance = 1e-6


function SparseArrays.sparsevec(d::DefaultDict, m::Int)
    SparseArrays.sparsevec(collect(keys(d)), collect(values(d)), m)
end


function matching_pursuit_(data::AbstractVector, dictionary::AbstractMatrix,
                           max_iter::Int, tolerance::Float64)
    n_atoms = size(dictionary, 2)

    residual = copy(data)

    xdict = DefaultDict{Int, Float64}(0.)
    for i in 1:max_iter
        if norm(residual) < tolerance
            return sparsevec(xdict, n_atoms)
        end

        # find an atom with maximum inner product
        products = dictionary' * residual
        _, maxindex = findmax(abs.(products))
        maxval = products[maxindex]
        atom = dictionary[:, maxindex]

        # c is the length of the projection of data onto atom
        a = maxval / sum(abs2, atom)  # equivalent to maxval / norm(atom)^2
        residual -= atom * a

        xdict[maxindex] += a
    end
    return sparsevec(xdict, n_atoms)
end


"""
    matching_pursuit(data::Vector, dictionary::AbstractMatrix;
                     max_iter::Int = $default_max_iter,
                     tolerance::Float64 = $default_tolerance)

Find ``x`` such that ``Dx = y`` or ``Dx ≈ y`` where y is `data` and D is `dictionary`.
```
# Arguments
* `max_iter`: Hard limit of iterations
* `tolerance`: Exit when the norm of the residual < tolerance
```
"""
function matching_pursuit(data::AbstractVector, dictionary::AbstractMatrix;
                          max_iter::Int = default_max_iter,
                          tolerance = default_tolerance)

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

    matching_pursuit_(data, dictionary, max_iter, tolerance)
end


"""
    matching_pursuit(data::AbstractMatrix, dictionary::AbstractMatrix;
                     max_iter::Int = $default_max_iter,
                     tolerance::Float64 = $default_tolerance)

Find ``X`` such that ``DX = Y`` or ``DX ≈ Y`` where Y is `data` and D is `dictionary`.
```
# Arguments
* `max_iter`: Hard limit of iterations
* `tolerance`: Exit when the norm of the residual < tolerance
```
"""
function matching_pursuit(data::AbstractMatrix, dictionary::AbstractMatrix;
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
