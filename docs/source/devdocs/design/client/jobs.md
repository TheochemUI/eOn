# Job Design

The key insight is that every job has roughly the same size, and also that they
make construct subsequent jobs which have varying arguments (especially in terms
of the number of matter objects). This means that the Jobs are less [LSP
compatible](https://www.wikiwand.com/en/Liskov_substitution_principle). To this
end, a modern C++ approach avoiding passing `void*` data and also keep value
semantics, it is best to utilize `std::variant` and `std::visit`. Note that:

- The size of `std::variant` depends on the largest member of the variant
- The corresponding visitor can invoke objects not constrained to an inheritance structure
 + Which is convinient for duck typing
 
The type of job to be constructed is taken from a configuration object, so it is
to be known at runtime only.

Another design is to leverage CRTP similar to the `Potential` class, and simply
forward all arguments.
