# KSVD.jl

[K-SVD](https://en.wikipedia.org/wiki/K-SVD) is an algorithm for creating overcomplete dictionaries for sparse representations.  

This package implements:

* K-SVD as described in the original paper: [K-SVD: An Algorithm for Designing Overcomplete Dictionaries for Sparse Representation](http://www.cs.technion.ac.il/~freddy/papers/120.pdf)
* [Matching Pursuit](https://en.wikipedia.org/wiki/Matching_pursuit) for representing signals using a given dictionary.

# Installation
Launch Julia and type

```julia
Pkg.add("KSVD")
```

# Usage

Assume that each column of Y represents a feature vector (or an input signal from some system).  
D is a dictionary. Each column of D represents an atom.  
K-SVD derives D and X such that DX ≈ Y from only Y.  

```julia
D, X = ksvd(
    Y,
    256,  # the number of atoms in D
    max_iter = 200,  # max iterations of K-SVD
    max_iter_mp = 40,  # max iterations of matching pursuit called in K-SVD
    sparsity_allowance = 0.96  # stop iteration when more than 96% of elements in X become zeros
)
```

[Matching Pursuit](https://en.wikipedia.org/wiki/Matching_pursuit) derives X from D and Y such that DX = Y in constraint that X be as sparse as possible.

```julia
X_sparse = matching_pursuit(Y, D, max_iter = 200)
```

# Example
Samples of [the digits dataset in scikit-learn ](http://scikit-learn.org/stable/auto_examples/datasets/plot_digits_last_image.html) and the obtained dictionary of 256 atoms.

![images](examples/digit_images.png)
![digits256](examples/digits256.png)

See [examples](examples) for more details.

# Provided functions

Only a few functions are provided: `ksvd` and `matching_pursuit`.
See [the documentation](docs/build/index.md).
