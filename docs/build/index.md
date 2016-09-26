
<a id='KSVD-Documentation-1'></a>

# KSVD Documentation

<a id='KSVD.ksvd-Tuple{Array{T,2},Int64}' href='#KSVD.ksvd-Tuple{Array{T,2},Int64}'>#</a>
**`KSVD.ksvd`** &mdash; *Method*.



```
ksvd(Y::Matrix, n_atoms::Int;
     tolerance_zeros::Float64 = 0.9,
     max_iter::Int = 200,
     max_iter_mp::Int = 200)
```

Run K-SVD that designs an efficient dictionary for sparse representations.

<a id='KSVD.matching_pursuit-Tuple{Array{T,1},Array{T,2}}' href='#KSVD.matching_pursuit-Tuple{Array{T,1},Array{T,2}}'>#</a>
**`KSVD.matching_pursuit`** &mdash; *Method*.



```
matching_pursuit(data::Vector, dictionary::Matrix;
                 max_iter::Int = 20,
                 tolerance::Float64 = 1.0e-6)
```

Find `x` such that `Dx = y` or `Dx ≈ y` where y is `data` and D is `dictionary`.

<a id='KSVD.matching_pursuit-Tuple{Array{T,2},Array{T,2}}' href='#KSVD.matching_pursuit-Tuple{Array{T,2},Array{T,2}}'>#</a>
**`KSVD.matching_pursuit`** &mdash; *Method*.



```
matching_pursuit(data::Matrix, dictionary::Matrix;
                 max_iter::Int = 20,
                 tolerance::Float64 = 1.0e-6)
```

Find `X` such that `DX = Y` or `DX ≈ Y` where Y is `data` and D is `dictionary`.

