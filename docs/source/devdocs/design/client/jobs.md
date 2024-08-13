# Job Design

The key insight is that every job has roughly the same size, and that they may
construct subsequent jobs with varying arguments (especially in terms of the
number of `Matter` objects). This variability in arguments means that the jobs
are less compatible with the [Liskov Substitution Principle
(LSP)](https://www.wikiwand.com/en/Liskov_substitution_principle). Specifically,
the differences in the types and numbers of arguments passed to each job make it
difficult to substitute one job type for another in a way that respects LSP.

To address these challenges, a modern C++ approach that avoids passing `void*`
data and maintains value semantics is to utilize `std::variant` and
`std::visit`. This approach offers several benefits:

- **Type Safety**: `std::variant` ensures that only the types defined in the
  variant can be used, providing strong type safety compared to traditional
  polymorphism.
- **Efficiency**: The size of `std::variant` depends on the largest member of
  the variant, and by using `std::visit`, we avoid the overhead of virtual
  function calls and dynamic memory allocation.
- **Flexibility**: The corresponding visitor can invoke objects not constrained
  to an inheritance structure, which is convenient for scenarios requiring duck
  typing or when extending behavior without modifying base classes.

The type of job to be constructed is determined at runtime, based on a
configuration object. This dynamic selection aligns well with the use of
`std::variant`, allowing for compile-time type safety while still supporting
runtime decisions.

### Alternative Approaches

Another design option is to leverage the Curiously Recurring Template Pattern
(CRTP), similar to the `Potential` class, and forward all arguments directly to
the specific job type. However, the main point is that virtual inheritance is a
poor model for the reality of `Job` instances. Each job exposes a similar
interface but requires differing arguments, making traditional polymorphism less
suitable.

### Considerations

For cases where more control over the job execution is needed, especially if
different job types require distinct handling based on runtime conditions or
configuration, `std::get_if` can be used instead of `std::visit`. This allows
for greater flexibility, enabling custom logic before or after the job’s
`runImpl` method is called.

While `std::variant` and `std::visit` provide strong guarantees and efficiency,
they can increase compile-time complexity, particularly as the number of job
types grows. It’s essential to balance these trade-offs based on the specific
requirements of the project.
