# KSVD Documentation

```@docs
ksvd(Y::Matrix, n_atoms::Int;
     tolerance_zeros::Float64 = default_tolerance_zeros,
     max_iter::Int = default_max_iter,
     max_iter_mp::Int = default_max_iter_mp)
matching_pursuit(data::Vector, dictionary::Matrix;
                 max_iter::Int = default_max_iter,
                 tolerance::Float64 = default_tolerance)
matching_pursuit(data::Matrix, dictionary::Matrix;
                 max_iter::Int = default_max_iter,
                 tolerance::Float64 = default_tolerance)
```
