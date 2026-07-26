// automatically generated 
// vesin version: 0.6.0

#include <cassert>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <numeric>
#include <tuple>
#include <type_traits>

#ifndef VESIN_CPU_CELL_LIST_HPP
#define VESIN_CPU_CELL_LIST_HPP

#include <cstddef>
#include <vector>

#include "vesin.h"

#ifndef VESIN_TYPES_HPP
#define VESIN_TYPES_HPP

#include <array>
#include <cassert>
#include <string>

#ifndef VESIN_MATH_HPP
#define VESIN_MATH_HPP

#include <array>
#include <cmath>
#include <stdexcept>

namespace vesin {
struct Vector;

Vector operator*(Vector vector, double scalar);

struct Vector: public std::array<double, 3> {
    double dot(Vector other) const {
        return (*this)[0] * other[0] + (*this)[1] * other[1] + (*this)[2] * other[2];
    }

    double norm() const {
        return std::sqrt(this->dot(*this));
    }

    Vector normalize() const {
        return *this * (1.0 / this->norm());
    }

    Vector cross(Vector other) const {
        return Vector{
            (*this)[1] * other[2] - (*this)[2] * other[1],
            (*this)[2] * other[0] - (*this)[0] * other[2],
            (*this)[0] * other[1] - (*this)[1] * other[0],
        };
    }
};

inline Vector operator+(Vector u, Vector v) {
    return Vector{
        u[0] + v[0],
        u[1] + v[1],
        u[2] + v[2],
    };
}

inline Vector operator-(Vector u, Vector v) {
    return Vector{
        u[0] - v[0],
        u[1] - v[1],
        u[2] - v[2],
    };
}

inline Vector operator*(double scalar, Vector vector) {
    return Vector{
        scalar * vector[0],
        scalar * vector[1],
        scalar * vector[2],
    };
}

inline Vector operator*(Vector vector, double scalar) {
    return Vector{
        scalar * vector[0],
        scalar * vector[1],
        scalar * vector[2],
    };
}

struct Matrix: public std::array<std::array<double, 3>, 3> {
    Matrix():
        std::array<std::array<double, 3>, 3>{{{0.0}}} {}

    Matrix(std::array<std::array<double, 3>, 3> values):
        std::array<std::array<double, 3>, 3>(values) {}

    double determinant() const {
        // clang-format off
        return (*this)[0][0] * ((*this)[1][1] * (*this)[2][2] - (*this)[2][1] * (*this)[1][2])
             - (*this)[0][1] * ((*this)[1][0] * (*this)[2][2] - (*this)[1][2] * (*this)[2][0])
             + (*this)[0][2] * ((*this)[1][0] * (*this)[2][1] - (*this)[1][1] * (*this)[2][0]);
        // clang-format on
    }

    Matrix inverse() const {
        auto det = this->determinant();

        if (std::abs(det) < 1e-30) {
            throw std::runtime_error("this matrix is not invertible");
        }

        auto inverse = Matrix();
        inverse[0][0] = ((*this)[1][1] * (*this)[2][2] - (*this)[2][1] * (*this)[1][2]) / det;
        inverse[0][1] = ((*this)[0][2] * (*this)[2][1] - (*this)[0][1] * (*this)[2][2]) / det;
        inverse[0][2] = ((*this)[0][1] * (*this)[1][2] - (*this)[0][2] * (*this)[1][1]) / det;
        inverse[1][0] = ((*this)[1][2] * (*this)[2][0] - (*this)[1][0] * (*this)[2][2]) / det;
        inverse[1][1] = ((*this)[0][0] * (*this)[2][2] - (*this)[0][2] * (*this)[2][0]) / det;
        inverse[1][2] = ((*this)[1][0] * (*this)[0][2] - (*this)[0][0] * (*this)[1][2]) / det;
        inverse[2][0] = ((*this)[1][0] * (*this)[2][1] - (*this)[2][0] * (*this)[1][1]) / det;
        inverse[2][1] = ((*this)[2][0] * (*this)[0][1] - (*this)[0][0] * (*this)[2][1]) / det;
        inverse[2][2] = ((*this)[0][0] * (*this)[1][1] - (*this)[1][0] * (*this)[0][1]) / det;
        return inverse;
    }
};

inline Vector operator*(Matrix matrix, Vector vector) {
    return Vector{
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    };
}

inline Vector operator*(Vector vector, Matrix matrix) {
    return Vector{
        vector[0] * matrix[0][0] + vector[1] * matrix[1][0] + vector[2] * matrix[2][0],
        vector[0] * matrix[0][1] + vector[1] * matrix[1][1] + vector[2] * matrix[2][1],
        vector[0] * matrix[0][2] + vector[1] * matrix[1][2] + vector[2] * matrix[2][2],
    };
}

} // namespace vesin

#endif

namespace vesin {

class BoundingBox {
public:
    BoundingBox(const BoundingBox&) = delete;
    BoundingBox& operator=(const BoundingBox&) = delete;

    BoundingBox(BoundingBox&&) = default;
    BoundingBox& operator=(BoundingBox&&) = default;

    BoundingBox(Matrix matrix, const bool periodic[3]):
        matrix_(matrix),
        periodic_({periodic[0], periodic[1], periodic[2]}),
        max_positions_({-1e300, -1e300, -1e300}),
        min_positions_({1e300, 1e300, 1e300}) {

        // find number of periodic directions and their indices
        int n_periodic = 0;
        int periodic_idx_1 = -1;
        int periodic_idx_2 = -1;
        for (int i = 0; i < 3; ++i) {
            if (periodic_[i]) {
                n_periodic += 1;
                if (periodic_idx_1 == -1) {
                    periodic_idx_1 = i;
                } else if (periodic_idx_2 == -1) {
                    periodic_idx_2 = i;
                }
            }
        }

        // adjust the box matrix to have a simple orthogonal dimension along
        // non-periodic directions
        if (n_periodic == 0) {
            matrix_ = Matrix({{
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, 1},
            }});
        } else if (n_periodic == 1) {
            assert(periodic_idx_1 != -1);
            // Make the two non-periodic directions orthogonal to the periodic one
            auto a = Vector{matrix_[periodic_idx_1]};
            auto b = Vector{0, 1, 0};
            if (std::abs(a.normalize().dot(b)) > 0.9) {
                b = Vector{0, 0, 1};
            }
            auto c = a.cross(b).normalize();
            b = c.cross(a).normalize();

            // Assign back to the matrix picking the "non-periodic" indices without ifs
            matrix_[(periodic_idx_1 + 1) % 3] = b;
            matrix_[(periodic_idx_1 + 2) % 3] = c;
        } else if (n_periodic == 2) {
            assert(periodic_idx_1 != -1 && periodic_idx_2 != -1);
            // Make the one non-periodic direction orthogonal to the two periodic ones
            auto a = Vector{matrix_[periodic_idx_1]};
            auto b = Vector{matrix_[periodic_idx_2]};
            auto c = a.cross(b).normalize();

            // Assign back to the matrix picking the "non-periodic" index without ifs
            matrix_[(3 - periodic_idx_1 - periodic_idx_2)] = c;
        }

        // precompute the inverse matrix
        auto det = matrix_.determinant();
        if (std::abs(det) < 1e-30) {
            throw std::runtime_error("the box matrix is not invertible");
        }

        this->inverse_ = matrix_.inverse();

        // precompute distances between faces of the bounding box
        auto a = Vector{matrix_[0]};
        auto b = Vector{matrix_[1]};
        auto c = Vector{matrix_[2]};

        // Plans normal vectors
        auto na = b.cross(c).normalize();
        auto nb = c.cross(a).normalize();
        auto nc = a.cross(b).normalize();

        distances_between_faces_ = Vector{
            periodic_[0] ? std::abs(na.dot(a)) : max_positions_[0] - min_positions_[0],
            periodic_[1] ? std::abs(nb.dot(b)) : max_positions_[1] - min_positions_[1],
            periodic_[2] ? std::abs(nc.dot(c)) : max_positions_[2] - min_positions_[2],
        };
    }

    const Matrix& matrix() const {
        return this->matrix_;
    }

    bool periodic(size_t spatial) const {
        return this->periodic_[spatial];
    }

    /// Convert a vector from cartesian coordinates to fractional coordinates
    ///
    /// For non-periodic dimensions, the fractional coordinates are not wrapped
    /// inside [0, 1], but are normalized by the corresponding box length.
    Vector cartesian_to_fractional(Vector cartesian) const {
        auto fractional = cartesian * inverse_;
        if (!periodic_[0]) {
            fractional[0] = (cartesian[0] - min_positions_[0]) / distances_between_faces_[0];
        }

        if (!periodic_[1]) {
            fractional[1] = (cartesian[1] - min_positions_[1]) / distances_between_faces_[1];
        }

        if (!periodic_[2]) {
            fractional[2] = (cartesian[2] - min_positions_[2]) / distances_between_faces_[2];
        }

        return fractional;
    }

    /// Convert a vector from fractional coordinates to cartesian coordinates
    Vector fractional_to_cartesian(Vector fractional) const {
        auto cartesian = fractional * matrix_;

        if (!periodic_[0]) {
            cartesian[0] *= distances_between_faces_[0];
            cartesian[0] += min_positions_[0];
        }

        if (!periodic_[1]) {
            cartesian[1] *= distances_between_faces_[1];
            cartesian[1] += min_positions_[1];
        }

        if (!periodic_[2]) {
            cartesian[2] *= distances_between_faces_[2];
            cartesian[2] += min_positions_[2];
        }

        return cartesian;
    }

    /// Get the three distances between faces of the bounding box
    Vector distances_between_faces() const {
        return distances_between_faces_;
    }

    void make_bounding_for(const double (*points)[3], size_t n_points) {
        make_bounding_for(reinterpret_cast<const Vector*>(points), n_points);
    }

    void make_bounding_for(const Vector* points, size_t n_points) {
        // find the min and max coordinates along each axis
        for (size_t i = 0; i < n_points; i++) {
            for (size_t spatial = 0; spatial < 3; spatial++) {
                if (!std::isfinite(points[i][spatial])) {
                    throw std::runtime_error(
                        "point " + std::to_string(i) + " has non-finite coordinate " +
                        "along axis " + std::to_string(spatial) + ": " +
                        std::to_string(points[i][spatial])
                    );
                }

                if (points[i][spatial] < min_positions_[spatial]) {
                    min_positions_[spatial] = points[i][spatial];
                }
                if (points[i][spatial] > max_positions_[spatial]) {
                    max_positions_[spatial] = points[i][spatial];
                }
            }
        }

        for (int dim = 0; dim < 3; dim++) {
            // if all atoms have the same coordinate in this dimension, pretend
            // that the bounding box is at least 1 unit wide to avoid numerical issues
            if (max_positions_[dim] - min_positions_[dim] < 1e-6) {
                max_positions_[dim] = min_positions_[dim] + 1;
            }

            if (!periodic_[dim]) {
                // add a 1% margin to make sure all points are strictly inside the
                // bounding box
                distances_between_faces_[dim] = (max_positions_[dim] - min_positions_[dim]) * 1.01;

                // make sure the distance is not too small (to prevent searching
                // too many cells down the line). This can happen if all point
                // are in the same plane in this direction
                if (distances_between_faces_[dim] < 1.0) {
                    distances_between_faces_[dim] = 1.0;
                }
            }
        }
    }

private:
    Matrix matrix_;
    std::array<bool, 3> periodic_;

    Matrix inverse_;
    Vector min_positions_;
    Vector max_positions_;
    Vector distances_between_faces_;
};

/// A cell shift represents the displacement along cell axis between the actual
/// position of an atom and a periodic image of this atom.
///
/// The cell shift can be used to reconstruct the vector between two points,
/// wrapped inside the unit cell.
struct CellShift: public std::array<int32_t, 3> {
    CellShift(int32_t a, int32_t b, int32_t c):
        std::array<int32_t, 3>({a, b, c}) {}

    CellShift(std::array<int32_t, 3> shifts):
        std::array<int32_t, 3>(shifts) {}

    CellShift(int32_t shifts[3]):
        std::array<int32_t, 3>({shifts[0], shifts[1], shifts[2]}) {}

    CellShift():
        std::array<int32_t, 3>({0, 0, 0}) {}

    /// Compute the shift vector in cartesian coordinates, using the given cell
    /// matrix (stored in row major order).
    Vector cartesian(const BoundingBox& box) const {
        assert(box.periodic(0) || (*this)[0] == 0);
        assert(box.periodic(1) || (*this)[1] == 0);
        assert(box.periodic(2) || (*this)[2] == 0);

        auto vector = Vector{
            static_cast<double>((*this)[0]),
            static_cast<double>((*this)[1]),
            static_cast<double>((*this)[2]),
        };
        return vector * box.matrix();
    }
};

inline CellShift operator+(CellShift a, CellShift b) {
    return CellShift{
        a[0] + b[0],
        a[1] + b[1],
        a[2] + b[2],
    };
}

inline CellShift operator-(CellShift a, CellShift b) {
    return CellShift{
        a[0] - b[0],
        a[1] - b[1],
        a[2] - b[2],
    };
}

} // namespace vesin

#endif

namespace vesin {
namespace cpu {

struct VerletList;

/// Extra CPU allocation metadata stored in `VesinNeighborList::opaque`.
struct ExtraDataCpu {
    /// Initialize empty CPU-side metadata.
    ExtraDataCpu() = default;
    /// Release optional Verlet cache state.
    ~ExtraDataCpu();

    /// Disallow copy construction; this object owns CPU-side cache metadata.
    ExtraDataCpu(const ExtraDataCpu&) = delete;
    /// Disallow copy assignment; this object owns CPU-side cache metadata.
    ExtraDataCpu& operator=(const ExtraDataCpu&) = delete;
    /// Disallow move construction; the C API stores this object by pointer.
    ExtraDataCpu(ExtraDataCpu&&) = delete;
    /// Disallow move assignment; the C API stores this object by pointer.
    ExtraDataCpu& operator=(ExtraDataCpu&&) = delete;

    /// Persisted GrowableNeighborList output capacity.
    size_t capacity = 0;
    /// Optional cached Verlet state for `skin > 0` calculations.
    VerletList* verlet_state = nullptr;
};

class ThreadPool;

void free_neighbors(VesinNeighborList& neighbors);

void cell_list_neighbors(
    const Vector* points,
    size_t n_points,
    BoundingBox box,
    VesinOptions options,
    VesinNeighborList& neighbors,
    size_t& capacity
);

// eOn local patch (upstream candidate): CPU implementation of
// VesinBruteForce. O(n^2) minimum-image pair search — one image per pair.
// Faster than the cell list when the cutoff is comparable to the box (the
// cell grid degenerates to a couple of cells and re-enumerates every pair
// per periodic shift) and for small systems. Requires every periodic box
// width to be at least the cutoff.
//
// When `visitor` is non-null it is invoked in the same scan for every pair
// with distance2 <= visit_cutoff^2 (see vesin_neighbors_visit).
void brute_force_neighbors(
    const Vector* points,
    size_t n_points,
    BoundingBox box,
    VesinOptions options,
    VesinNeighborList& neighbors,
    size_t& capacity,
    VesinPairVisitor visitor = nullptr,
    void* visitor_data = nullptr,
    double visit_cutoff = 0.0
);

void neighbors(
    const Vector* points,
    size_t n_points,
    BoundingBox box,
    VesinOptions options,
    VesinNeighborList& neighbors
);

/// The cell list is used to sort atoms inside bins/cells.
///
/// The list of potential pairs is then constructed by looking through all
/// neighboring cells (the number of cells to search depends on the cutoff and
/// the size of the cells) for each atom to create pair candidates.
class CellList {
public:
    /// Create a new `CellList` for the given bounding box and cutoff,
    /// determining all required parameters.
    CellList(BoundingBox box, double cutoff);

    /// Add a single point to the cell list at the given `position`. The point
    /// is uniquely identified by its `index`.
    void add_point(size_t index, Vector position);

    /// Iterate over all possible pairs, calling the given callback every time
    template <typename Function>
    void foreach_pair(ThreadPool& thread_pool, size_t n_threads, Function callback);

    /// Number of cells in this list.
    size_t n_cells() const {
        return cells_shape_[0] * cells_shape_[1] * cells_shape_[2];
    }

private:
    /// How many cells do we need to look at when searching neighbors to include
    /// all neighbors below cutoff
    std::array<int32_t, 3> n_search_;

    /// the cells themselves are a list of points & corresponding
    /// shift to place the point inside the cell
    struct Point {
        size_t index;
        CellShift shift;
    };
    struct Cell: public std::vector<Point> {};

    // raw data for the cells
    std::vector<Cell> cells_;
    // shape of the cell array
    std::array<size_t, 3> cells_shape_;

    BoundingBox box_;

    Cell& get_cell(std::array<int32_t, 3> index);
    const Cell& get_cell(std::array<int32_t, 3> index) const;
};

/// Wrapper around `VesinNeighborList` that behaves like a std::vector,
/// automatically growing memory allocations.
class GrowableNeighborList {
public:
    VesinNeighborList& neighbors;
    size_t capacity;
    VesinOptions options;

    size_t length() const {
        return neighbors.length;
    }

    void increment_length() {
        neighbors.length += 1;
    }

    void set_pair(size_t index, size_t first, size_t second);
    void set_shift(size_t index, vesin::CellShift shift);
    void set_distance(size_t index, double distance);
    void set_vector(size_t index, vesin::Vector vector);

    // reset length to 0, and allocate/deallocate members of
    // `neighbors` according to `options`
    void reset();

    // allocate more memory & update capacity, `new_size = 0` means doubling the
    // current capacity, otherwise grow to at least `new_size` capacity. This
    // function does not update the `length` of `neighbors`.
    void grow(size_t new_size = 0);

    // sort the pairs currently in the neighbor list
    void sort();
};

} // namespace cpu
} // namespace vesin

#endif
#ifndef VESIN_CPU_THREADPOOL_HPP
#define VESIN_CPU_THREADPOOL_HPP

#include <algorithm>
#include <exception>
#include <vector>

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

namespace vesin {
namespace cpu {

/// Small reusable thread pool for CPU parallel regions.
///
/// The pool always creates one worker slot per available CPU core (according to
/// `std::thread::hardware_concurrency()`), and reuses worker threads across
/// multiple `run` calls. Thread 0 is always the caller of `run`, while workers
/// handle thread IDs 1..max_threads()-1.
///
/// A `run` call is identified by a monotonically increasing `generation_`:
/// workers remember the last generation they have executed and only wake when a
/// new one appears.
class ThreadPool {
public:
    ThreadPool():
        n_threads_(std::max<size_t>(1, static_cast<size_t>(std::thread::hardware_concurrency()))) {
        // eOn local patch: workers spawn lazily in run() on the first
        // parallel dispatch. Serial consumers (options.n_threads == 1) and
        // short-lived processes must not pay hardware_concurrency() thread
        // creations just for constructing the pool.
        workers_.reserve(n_threads_ - 1);
    }

    ~ThreadPool() {
        {
            auto lock = std::lock_guard<std::mutex>(mutex_);
            stopping_ = true;
            generation_ += 1;
        }
        start_cv_.notify_all();

        for (auto& worker_thread : workers_) {
            worker_thread.join();
        }
    }

    /// Maximum number of threads available in this pool.
    ///
    /// This is equal to the number of hardware threads detected at
    /// construction time, clamped to at least 1.
    size_t max_threads() const {
        return n_threads_;
    }

    /// Execute `n_tasks` invocations of `task(task_i, thread_id)`.
    ///
    /// `active_threads` selects how many threads from the pool should
    /// participate in this run. It is clamped to `[1, min(max_threads(), n_tasks)]`.
    ///
    /// Tasks are split deterministically into contiguous chunks based on
    /// `thread_id` among active threads; this keeps ordering stable across
    /// runs.
    ///
    /// Exceptions are captured and rethrown on the caller thread after workers
    /// finish/abort.
    template <typename Function>
    void run(size_t active_threads, size_t n_tasks, Function task) {
        if (n_tasks == 0) {
            return;
        }

        auto run_lock = std::unique_lock<std::mutex>(run_mutex_);

        active_threads = std::max<size_t>(1, active_threads);
        active_threads = std::min(active_threads, n_threads_);
        active_threads = std::min(active_threads, n_tasks);

        if (active_threads <= 1 || n_tasks <= 1) {
            for (size_t task_i = 0; task_i < n_tasks; task_i++) {
                task(task_i, 0);
            }
            return;
        }

        // Lazy worker spawn (see constructor). Workers start with
        // seen_generation = 0, so any worker that begins after the
        // generation bump below still picks up this run's tasks; run()
        // blocks until every active worker has executed its chunk.
        if (workers_.empty() && n_threads_ > 1) {
            for (size_t thread_id = 1; thread_id < n_threads_; thread_id++) {
                workers_.emplace_back([this, thread_id]() {
                    this->worker(thread_id);
                });
            }
        }

        {
            auto lock = std::lock_guard<std::mutex>(mutex_);
            n_tasks_ = n_tasks;
            task_data_ = &task;
            task_function_ = [](void* data, size_t task_i, size_t thread_id) {
                auto* function = static_cast<Function*>(data);
                (*function)(task_i, thread_id);
            };
            first_exception_ = nullptr;
            has_exception_.store(false);
            active_threads_ = active_threads;
            running_workers_ = active_threads_ - 1;
            generation_ += 1;
        }
        start_cv_.notify_all();

        /// The caller also participates in the work, so we execute its chunk
        /// directly here while workers are running their chunks in the background.
        this->execute_assigned_tasks(0, active_threads);

        auto lock = std::unique_lock<std::mutex>(mutex_);
        done_cv_.wait(lock, [this]() {
            return running_workers_ == 0;
        });

        if (first_exception_ != nullptr) {
            std::rethrow_exception(first_exception_);
        }
    }

private:
    /// Main loop for background workers.
    ///
    /// Workers block on `start_cv_` while idle, wake for each new generation,
    /// execute their chunk, then decrement `running_workers_`.
    void worker(size_t thread_id) {
        size_t seen_generation = 0;

        auto lock = std::unique_lock<std::mutex>(mutex_);
        while (true) {
            start_cv_.wait(lock, [this, seen_generation]() {
                return stopping_ || generation_ != seen_generation;
            });

            if (stopping_) {
                return;
            }

            seen_generation = generation_;
            auto is_active = thread_id < active_threads_;
            auto active_threads = active_threads_;
            lock.unlock();
            if (is_active) {
                this->execute_assigned_tasks(thread_id, active_threads);
            }
            lock.lock();

            if (is_active) {
                running_workers_ -= 1;
                if (running_workers_ == 0) {
                    done_cv_.notify_one();
                }
            }
        }
    }

    /// Execute this thread's deterministic chunk for the current generation.
    void execute_assigned_tasks(size_t thread_id, size_t active_threads) {
        auto begin = (thread_id * n_tasks_) / active_threads;
        auto end = ((thread_id + 1) * n_tasks_) / active_threads;

        for (size_t task_i = begin; task_i < end; task_i++) {
            if (has_exception_.load()) {
                return;
            }

            try {
                task_function_(task_data_, task_i, thread_id);
            } catch (...) {
                auto lock = std::lock_guard<std::mutex>(mutex_);
                if (first_exception_ == nullptr) {
                    first_exception_ = std::current_exception();
                    has_exception_.store(true);
                }
                return;
            }
        }
    }

    size_t n_threads_;
    std::vector<std::thread> workers_;

    /// Serialize `run` calls on this pool.
    ///
    /// The current implementation stores one shared run context; allowing
    /// concurrent `run` calls would race on that state.
    std::mutex run_mutex_;

    /// Protects shared run state (`stopping_`, `generation_`, `running_workers_`,
    /// `task_`, `n_tasks_`, and `first_exception_`).
    std::mutex mutex_;
    /// Signaled when a new run generation starts or shutdown begins.
    std::condition_variable start_cv_;
    /// Signaled when all background workers completed the current run.
    std::condition_variable done_cv_;

    /// Are we trying to stop the pool? If true, workers should exit as soon as
    /// possible
    bool stopping_ = false;
    /// Incremented for each new `run` and for shutdown.
    ///
    /// Together with each worker's local `seen_generation`, this separates real
    /// work starts from spurious wakeups and prevents re-running the same batch.
    size_t generation_ = 0;
    /// Number of background workers (thread IDs >= 1) still active in the
    /// current generation.
    size_t running_workers_ = 0;
    /// Number of threads participating in the current generation.
    size_t active_threads_ = 1;

    using task_function_t = void (*)(void*, size_t, size_t);

    size_t n_tasks_ = 0;
    void* task_data_ = nullptr;
    task_function_t task_function_ = nullptr;

    /// Store the first exception thrown by any worker, if any, to rethrow on
    /// the caller thread after workers finish/abort.
    std::exception_ptr first_exception_;
    std::atomic<bool> has_exception_ = false;
};

} // namespace cpu
} // namespace vesin

#endif // VESIN_CPU_THREADPOOL_HPP
#ifndef VESIN_VERLET_HPP
#define VESIN_VERLET_HPP

#include <array>
#include <cstddef>
#include <vector>


namespace vesin {
namespace cpu {

/// On-CPU Verlet neighbor list.
///
/// This class manages an over-complete candidate neighbor list, computed for
/// `cutoff + skin`, and the associated reference state for validating
/// neighbor-list reuse.
struct VerletList {
    /// Initialize an empty Verlet list with no cached candidates.
    VerletList() = default;
    ~VerletList();

    VerletList(const VerletList&) = delete;
    VerletList& operator=(const VerletList&) = delete;
    VerletList(VerletList&&) = delete;
    VerletList& operator=(VerletList&&) = delete;

    /// Copy options relevant to cache reuse from the user options.
    ///
    /// If the cache-driving options changed (cutoff, skin, full-list flag),
    /// clear existing cached candidates.
    void set_options(VesinOptions options);

    /// Return `true` if the cache should be rebuilt for the current `points`.
    ///
    /// The caller passes the current points and box; a rebuild is required
    /// after any structural change that could invalidate cached candidates (box
    /// topology/periodicity changes, points count changes, or large
    /// displacement from the reference points).
    bool needs_rebuild(
        const Vector* points,
        size_t n_points,
        const BoundingBox& box
    ) const;

    /// Build the over-complete candidate list at `cutoff + skin`.
    ///
    /// The rebuild operation stores candidates in full `VesinNeighborList` form
    /// and captures the state used to validate future `needs_rebuild` checks.
    void rebuild(
        const Vector* points,
        size_t n_points,
        const BoundingBox& box
    );

    /// Filter cached candidates at the exact user cutoff for the current
    /// `points`.
    void filter(
        const Vector* points,
        const BoundingBox& box,
        VesinOptions options,
        VesinNeighborList& neighbors,
        size_t& output_capacity
    ) const;

    /// Number of pairs currently stored in the cached candidate list.
    size_t candidate_count() const {
        return candidates_.length;
    }

private:
    /// Release candidate buffers and reset cache metadata.
    void clear_candidates();

    /// Reference points at the time the candidates were built.
    std::vector<Vector> ref_points_;
    /// Box matrix used for candidate generation and displacement validation.
    Matrix ref_matrix_;
    /// Periodicity flags for the cached box.
    std::array<bool, 3> ref_periodic_ = {false, false, false};

    /// Over-complete candidate list generated at `cutoff + skin`.
    ///
    /// The list is kept in normal neighbor-list representation so rebuild and
    /// recompute paths can share storage and filtering logic.
    VesinNeighborList candidates_;

    /// Options used to build the current cache.
    VesinOptions options_ = {};
    /// Rebuild threshold used to invalidate a cache (`(skin/2)^2`).
    double half_skin_sq_ = 0.0;

    /// Whether a usable cache is currently available.
    bool has_cache_ = false;
};

} // namespace cpu
} // namespace vesin

#endif

using namespace vesin::cpu;

static void free_extra_data(VesinNeighborList& neighbors) {
    delete static_cast<ExtraDataCpu*>(neighbors.opaque);
    neighbors.opaque = nullptr;
}

static ExtraDataCpu& extra_data(VesinNeighborList& neighbors) {
    if (neighbors.opaque == nullptr) {
        neighbors.opaque = new ExtraDataCpu();
    }

    return *static_cast<ExtraDataCpu*>(neighbors.opaque);
}

static void validate_algorithm(VesinOptions options) {
    if (options.algorithm != VesinAutoAlgorithm &&
        options.algorithm != VesinCellList &&
        options.algorithm != VesinBruteForce) {
        throw std::runtime_error("unknown algorithm on CPU");
    }
}

ExtraDataCpu::~ExtraDataCpu() {
    delete this->verlet_state;
}

struct ThreadLocalNeighborLists {
    ThreadLocalNeighborLists(size_t n_threads):
        raw(n_threads - 1) {
        for (auto& list : raw) {
            list.device = {VesinCPU, 0};
        }
    }

    std::vector<VesinNeighborList> raw;

    ~ThreadLocalNeighborLists() {
        for (auto& list : raw) {
            vesin::cpu::free_neighbors(list);
        }
    }
};

static ThreadPool& global_thread_pool() {
    static auto pool = ThreadPool();
    return pool;
}

void vesin::cpu::brute_force_neighbors(
    const Vector* points,
    size_t n_points,
    BoundingBox box,
    VesinOptions options,
    VesinNeighborList& raw_neighbors,
    size_t& capacity,
    VesinPairVisitor visitor,
    void* visitor_data,
    double visit_cutoff
) {
    auto visit_cutoff2 = visit_cutoff * visit_cutoff;
    // Exactly one pair per (i, j) — the nearest periodic image. Callers own
    // minimum-image validity: with a periodic width between one and two
    // cutoffs, second images inside the cutoff are NOT reported (use the
    // cell list when all images matter).
    auto widths = box.distances_between_faces();
    for (size_t k = 0; k < 3; k++) {
        if (box.periodic(k) && widths[k] < options.cutoff) {
            throw std::runtime_error(
                "brute force algorithm requires every periodic box width to "
                "be at least the cutoff; use the cell list instead"
            );
        }
    }

    auto matrix = box.matrix();
    auto inverse = matrix.inverse();
    auto cutoff2 = options.cutoff * options.cutoff;
    auto any_periodic = box.periodic(0) || box.periodic(1) || box.periodic(2);

    auto initial_capacity = std::max(capacity, raw_neighbors.length);
    auto neighbors = GrowableNeighborList{raw_neighbors, initial_capacity, options};
    neighbors.reset();

    auto want_shifts = options.return_shifts;
    auto want_distances = options.return_distances;
    auto want_vectors = options.return_vectors;
    auto full = options.full;

    // Direct writes with one capacity check per pair: the per-field set_*
    // members each re-check capacity through a non-inlined call, which
    // dominates dense O(n^2) fills.
    auto emit = [&](size_t first, size_t second, CellShift shift, Vector vector, double distance2) {
        auto index = neighbors.length();
        while (index + (full ? 2 : 1) > neighbors.capacity) {
            neighbors.grow();
        }
        auto& nb = neighbors.neighbors;
        nb.pairs[index][0] = first;
        nb.pairs[index][1] = second;
        double distance = 0.0;
        if (want_shifts) {
            nb.shifts[index][0] = shift[0];
            nb.shifts[index][1] = shift[1];
            nb.shifts[index][2] = shift[2];
        }
        if (want_distances) {
            distance = std::sqrt(distance2);
            nb.distances[index] = distance;
        }
        if (want_vectors) {
            nb.vectors[index][0] = vector[0];
            nb.vectors[index][1] = vector[1];
            nb.vectors[index][2] = vector[2];
        }
        nb.length = index + 1;

        if (full) {
            index += 1;
            nb.pairs[index][0] = second;
            nb.pairs[index][1] = first;
            if (want_shifts) {
                nb.shifts[index][0] = -shift[0];
                nb.shifts[index][1] = -shift[1];
                nb.shifts[index][2] = -shift[2];
            }
            if (want_distances) {
                nb.distances[index] = distance;
            }
            if (want_vectors) {
                nb.vectors[index][0] = -vector[0];
                nb.vectors[index][1] = -vector[1];
                nb.vectors[index][2] = -vector[2];
            }
            nb.length = index + 1;
        }
    };

    auto diagonal = matrix[0][1] == 0 && matrix[0][2] == 0 &&
                    matrix[1][0] == 0 && matrix[1][2] == 0 &&
                    matrix[2][0] == 0 && matrix[2][1] == 0;

    if (!any_periodic) {
        // Free boundaries: plain distance filter, no shifts at all.
        for (size_t i = 0; i < n_points; i++) {
            for (size_t j = i + 1; j < n_points; j++) {
                auto vector = points[j] - points[i];
                auto distance2 = vector.dot(vector);
                if (distance2 < cutoff2) {
                    if (visitor != nullptr && distance2 <= visit_cutoff2) {
                        visitor(visitor_data, i, j, vector[0], vector[1], vector[2], distance2);
                    }
                    emit(i, j, CellShift(), vector, distance2);
                }
            }
        }
    } else if (diagonal) {
        // Orthorhombic: fold each dimension directly, no fractional
        // round trip. inv == 0 turns the fold into a no-op for
        // non-periodic dimensions (floor(0.5) == 0).
        double width[3];
        double inv[3];
        for (size_t k = 0; k < 3; k++) {
            width[k] = matrix[k][k];
            inv[k] = (box.periodic(k) && width[k] != 0.0) ? 1.0 / width[k] : 0.0;
        }

        for (size_t i = 0; i < n_points; i++) {
            for (size_t j = i + 1; j < n_points; j++) {
                auto vector = points[j] - points[i];
                auto shift = CellShift();
                for (size_t k = 0; k < 3; k++) {
                    // round-to-nearest via truncate-cast: no libm floor call
                    // on -march targets without roundsd
                    auto t = vector[k] * inv[k];
                    auto s = -static_cast<double>(static_cast<int64_t>(
                        t + std::copysign(0.5, t)));
                    shift[k] = static_cast<int32_t>(s);
                    vector[k] += s * width[k];
                }

                auto distance2 = vector.dot(vector);
                if (distance2 < cutoff2) {
                    if (visitor != nullptr && distance2 <= visit_cutoff2) {
                        visitor(visitor_data, i, j, vector[0], vector[1], vector[2], distance2);
                    }
                    emit(i, j, shift, vector, distance2);
                }
            }
        }
    } else {
        for (size_t i = 0; i < n_points; i++) {
            for (size_t j = i + 1; j < n_points; j++) {
                auto vector = points[j] - points[i];

                auto shift = CellShift();
                auto fractional = vector * inverse;
                for (size_t k = 0; k < 3; k++) {
                    if (box.periodic(k)) {
                        shift[k] = static_cast<int32_t>(-std::floor(fractional[k] + 0.5));
                    }
                }
                if (shift[0] != 0 || shift[1] != 0 || shift[2] != 0) {
                    auto cartesian_shift = Vector{
                        static_cast<double>(shift[0]),
                        static_cast<double>(shift[1]),
                        static_cast<double>(shift[2]),
                    } * matrix;
                    vector = vector + cartesian_shift;
                }

                auto distance2 = vector.dot(vector);
                if (distance2 < cutoff2) {
                    if (visitor != nullptr && distance2 <= visit_cutoff2) {
                        visitor(visitor_data, i, j, vector[0], vector[1], vector[2], distance2);
                    }
                    emit(i, j, shift, vector, distance2);
                }
            }
        }
    }

    if (options.sorted) {
        neighbors.sort();
    }

    capacity = neighbors.capacity;
}

void vesin::cpu::cell_list_neighbors(
    const Vector* points,
    size_t n_points,
    BoundingBox box,
    VesinOptions options,
    VesinNeighborList& raw_neighbors,
    size_t& capacity
) {
    validate_algorithm(options);

    auto cell_list = CellList(std::move(box), options.cutoff);

    for (size_t i = 0; i < n_points; i++) {
        cell_list.add_point(i, points[i]);
    }

    auto cutoff2 = options.cutoff * options.cutoff;

    // the cell list creates too many pairs, we only need to keep the
    // one where the distance is actually below the cutoff
    auto initial_capacity = std::max(capacity, raw_neighbors.length);
    auto neighbors = GrowableNeighborList{raw_neighbors, initial_capacity, options};
    neighbors.reset();

    auto thread_count = std::max<size_t>(1, static_cast<size_t>(options.n_threads));
    thread_count = std::min(thread_count, cell_list.n_cells());
    auto& thread_pool = global_thread_pool();

    auto thread_locals = ThreadLocalNeighborLists(thread_count);

    auto local_thread_neighbors = std::vector<GrowableNeighborList>();
    local_thread_neighbors.reserve(thread_locals.raw.size());
    for (auto& local_raw : thread_locals.raw) {
        local_raw.device = {VesinCPU, 0};
        auto list = GrowableNeighborList{local_raw, 0, options};
        list.reset();
        local_thread_neighbors.emplace_back(list);
    }

    cell_list.foreach_pair(thread_pool, thread_count, [&](size_t first, size_t second, CellShift shift, size_t thread_id) {
        if (!options.full) {
            // filter out some pairs for half neighbor lists
            if (first > second) {
                return;
            }

            if (first == second) {
                // When creating pairs between a point and one of its periodic
                // images, the code generate multiple redundant pairs (e.g. with
                // shifts 0 1 1 and 0 -1 -1); and we want to only keep one of
                // these.
                if (shift[0] + shift[1] + shift[2] < 0) {
                    // drop shifts on the negative half-space
                    return;
                }

                if ((shift[0] + shift[1] + shift[2] == 0) && (shift[2] < 0 || (shift[2] == 0 && shift[1] < 0))) {
                    // drop shifts in the negative half plane or the negative
                    // shift[1] axis. See below for a graphical representation:
                    // we are keeping the shifts indicated with `O` and dropping
                    // the ones indicated with `X`
                    //
                    //  O O O │ O O O
                    //  O O O │ O O O
                    //  O O O │ O O O
                    // ─X─X─X─┼─O─O─O─
                    //  X X X │ X X X
                    //  X X X │ X X X
                    //  X X X │ X X X
                    return;
                }
            }
        }

        auto vector = points[second] - points[first] + shift.cartesian(box);
        auto distance2 = vector.dot(vector);

        if (distance2 < cutoff2) {
            GrowableNeighborList* local_list = nullptr;
            if (thread_id == 0) {
                local_list = &neighbors;
            } else {
                local_list = &local_thread_neighbors[thread_id - 1];
            }

            auto index = local_list->length();
            local_list->set_pair(index, first, second);

            if (options.return_shifts) {
                local_list->set_shift(index, shift);
            }

            if (options.return_distances) {
                local_list->set_distance(index, std::sqrt(distance2));
            }

            if (options.return_vectors) {
                local_list->set_vector(index, vector);
            }

            local_list->increment_length();
        }
    });

    // Threads 1 to n-1 have their neighbors in `thread_locals`; merge them
    // into the list created by thread 0 (`raw_neighbors`).
    auto final_length = neighbors.length();
    auto local_offsets = std::vector<size_t>(thread_locals.raw.size(), 0);
    for (size_t local_i = 0; local_i < thread_locals.raw.size(); local_i++) {
        local_offsets[local_i] = final_length;
        final_length += thread_locals.raw[local_i].length;
    }

    neighbors.grow(final_length);

    thread_pool.run(thread_count, thread_locals.raw.size(), [&](size_t local_i, size_t /*thread_id*/) {
        const auto& local_neighbors = thread_locals.raw[local_i];
        auto offset = local_offsets[local_i];

        for (size_t i = 0; i < local_neighbors.length; i++) {
            auto index = offset + i;
            neighbors.neighbors.pairs[index][0] = local_neighbors.pairs[i][0];
            neighbors.neighbors.pairs[index][1] = local_neighbors.pairs[i][1];

            if (options.return_shifts) {
                neighbors.neighbors.shifts[index][0] = local_neighbors.shifts[i][0];
                neighbors.neighbors.shifts[index][1] = local_neighbors.shifts[i][1];
                neighbors.neighbors.shifts[index][2] = local_neighbors.shifts[i][2];
            }

            if (options.return_distances) {
                neighbors.neighbors.distances[index] = local_neighbors.distances[i];
            }

            if (options.return_vectors) {
                neighbors.neighbors.vectors[index][0] = local_neighbors.vectors[i][0];
                neighbors.neighbors.vectors[index][1] = local_neighbors.vectors[i][1];
                neighbors.neighbors.vectors[index][2] = local_neighbors.vectors[i][2];
            }
        }
    });

    neighbors.neighbors.length = final_length;

    if (options.sorted) {
        neighbors.sort();
    }

    capacity = neighbors.capacity;
}

void vesin::cpu::neighbors(
    const Vector* points,
    size_t n_points,
    BoundingBox box,
    VesinOptions options,
    VesinNeighborList& raw_neighbors
) {
    validate_algorithm(options);

    auto& extra = extra_data(raw_neighbors);
    if (options.skin > 0.0) {
        if (extra.verlet_state == nullptr) {
            extra.verlet_state = new VerletList();
        }

        auto& state = *extra.verlet_state;
        state.set_options(options);

        if (state.needs_rebuild(points, n_points, box)) {
            state.rebuild(points, n_points, box);
        }

        state.filter(points, box, options, raw_neighbors, extra.capacity);
    } else {
        // Remove verlet state if it exists from a previous call
        delete extra.verlet_state;
        extra.verlet_state = nullptr;

        if (options.algorithm == VesinBruteForce) {
            brute_force_neighbors(points, n_points, std::move(box), options, raw_neighbors, extra.capacity);
        } else {
            cell_list_neighbors(points, n_points, std::move(box), options, raw_neighbors, extra.capacity);
        }
    }
}

/* ========================================================================== */

/// Maximal number of cells, we need to use this to prevent having too many
/// cells with a small bounding box and a large cutoff
#define MAX_NUMBER_OF_CELLS 1e5

/// Function to compute both quotient and remainder of the division of a by b.
/// This function follows Python convention, making sure the remainder have the
/// same sign as `b`.
static std::tuple<int32_t, int32_t> divmod(int32_t a, size_t b) {
    assert(b < (std::numeric_limits<int32_t>::max()));
    auto b_32 = static_cast<int32_t>(b);
    auto quotient = a / b_32;
    auto remainder = a % b_32;
    if (remainder < 0) {
        remainder += b_32;
        quotient -= 1;
    }
    return std::make_tuple(quotient, remainder);
}

/// Apply the `divmod` function to three components at the time
static std::tuple<std::array<int32_t, 3>, std::array<int32_t, 3>>
divmod(std::array<int32_t, 3> a, std::array<size_t, 3> b) {
    auto [qx, rx] = divmod(a[0], b[0]);
    auto [qy, ry] = divmod(a[1], b[1]);
    auto [qz, rz] = divmod(a[2], b[2]);
    return std::make_tuple(
        std::array<int32_t, 3>{qx, qy, qz},
        std::array<int32_t, 3>{rx, ry, rz}
    );
}

CellList::CellList(BoundingBox box, double cutoff):
    n_search_({0, 0, 0}),
    cells_shape_({0, 0, 0}),
    box_(std::move(box)) {
    auto distances_between_faces = box_.distances_between_faces();

    auto n_cells = Vector{
        std::clamp(std::trunc(distances_between_faces[0] / cutoff), 1.0, HUGE_VAL),
        std::clamp(std::trunc(distances_between_faces[1] / cutoff), 1.0, HUGE_VAL),
        std::clamp(std::trunc(distances_between_faces[2] / cutoff), 1.0, HUGE_VAL),
    };

    assert(std::isfinite(n_cells[0]) && std::isfinite(n_cells[1]) && std::isfinite(n_cells[2]));

    // limit memory consumption by ensuring we have less than `MAX_N_CELLS`
    // cells to look though
    auto n_cells_total = n_cells[0] * n_cells[1] * n_cells[2];
    if (n_cells_total > MAX_NUMBER_OF_CELLS) {
        // set the total number of cells close to MAX_N_CELLS, while keeping
        // roughly the ratio of cells in each direction
        auto ratio_x_y = n_cells[0] / n_cells[1];
        auto ratio_y_z = n_cells[1] / n_cells[2];

        n_cells[2] = std::trunc(std::cbrt(MAX_NUMBER_OF_CELLS / (ratio_x_y * ratio_y_z * ratio_y_z)));
        n_cells[1] = std::trunc(ratio_y_z * n_cells[2]);
        n_cells[0] = std::trunc(ratio_x_y * n_cells[1]);
    }

    // number of cells to search in each direction to make sure all possible
    // pairs below the cutoff are accounted for.
    this->n_search_ = std::array<int32_t, 3>{
        static_cast<int32_t>(std::ceil(cutoff * n_cells[0] / distances_between_faces[0])),
        static_cast<int32_t>(std::ceil(cutoff * n_cells[1] / distances_between_faces[1])),
        static_cast<int32_t>(std::ceil(cutoff * n_cells[2] / distances_between_faces[2])),
    };

    this->cells_shape_ = std::array<size_t, 3>{
        static_cast<size_t>(n_cells[0]),
        static_cast<size_t>(n_cells[1]),
        static_cast<size_t>(n_cells[2]),
    };

    for (size_t spatial = 0; spatial < 3; spatial++) {
        if (n_search_[spatial] < 1) {
            n_search_[spatial] = 1;
        }
    }

    // When there is mixed periodicity (some dimensions periodic, some not),
    // periodic shifts can bring atoms together across multiple cells in
    // non-periodic dimensions. Increase the search range in non-periodic
    // dimensions to cover the full cell extent in those dimensions.
    bool has_periodic = false;
    bool has_non_periodic = false;
    for (size_t d = 0; d < 3; d++) {
        if (box_.periodic(d)) {
            has_periodic = true;
        } else {
            has_non_periodic = true;
        }
    }

    if (has_periodic && has_non_periodic) {
        for (size_t spatial = 0; spatial < 3; spatial++) {
            if (!box_.periodic(spatial)) {
                n_search_[spatial] = std::max(
                    n_search_[spatial],
                    static_cast<int32_t>(cells_shape_[spatial]) - 1
                );
            }
        }
    }

    this->cells_.resize(cells_shape_[0] * cells_shape_[1] * cells_shape_[2]);
}

void CellList::add_point(size_t index, Vector position) {
    auto fractional = box_.cartesian_to_fractional(position);

    // find the cell in which this atom should go
    auto cell_index = std::array<int32_t, 3>{
        static_cast<int32_t>(std::floor(fractional[0] * static_cast<double>(cells_shape_[0]))),
        static_cast<int32_t>(std::floor(fractional[1] * static_cast<double>(cells_shape_[1]))),
        static_cast<int32_t>(std::floor(fractional[2] * static_cast<double>(cells_shape_[2]))),
    };

    // deal with pbc by wrapping the atom inside if it was outside of the cell
    CellShift shift;
    for (size_t spatial = 0; spatial < 3; spatial++) {
        auto result = divmod(cell_index[spatial], cells_shape_[spatial]);
        shift[spatial] = std::get<0>(result);
        cell_index[spatial] = std::get<1>(result);

        assert(box_.periodic(spatial) || shift[spatial] == 0);
    }

    this->get_cell(cell_index).emplace_back(Point{index, shift});
}

template <typename Function>
void CellList::foreach_pair(ThreadPool& thread_pool, size_t n_threads, Function callback) {
    auto total_cells = cells_shape_[0] * cells_shape_[1] * cells_shape_[2];

    thread_pool.run(n_threads, total_cells, [&](size_t linear_cell_i, size_t thread_id) {
        auto stride_x = cells_shape_[0];
        auto stride_xy = cells_shape_[0] * cells_shape_[1];

        auto cell_i_x = static_cast<int32_t>(linear_cell_i % stride_x);
        auto cell_i_y = static_cast<int32_t>((linear_cell_i / stride_x) % cells_shape_[1]);
        auto cell_i_z = static_cast<int32_t>(linear_cell_i / stride_xy);

        const auto& current_cell = this->get_cell({cell_i_x, cell_i_y, cell_i_z});

        // look through each neighboring cell
        for (int32_t delta_x = -n_search_[0]; delta_x <= n_search_[0]; delta_x++) {
            for (int32_t delta_y = -n_search_[1]; delta_y <= n_search_[1]; delta_y++) {
                for (int32_t delta_z = -n_search_[2]; delta_z <= n_search_[2]; delta_z++) {
                    // shift vector from one cell to the other
                    auto cell_shift = std::array<int32_t, 3>{0, 0, 0};
                    // index of the neighboring cell
                    auto neighbor_cell_i = std::array<int32_t, 3>{
                        cell_i_x + delta_x,
                        cell_i_y + delta_y,
                        cell_i_z + delta_z,
                    };

                    // only wrap (i.e. call divmod) cell indices in periodic dimensions,
                    // skip out-of-bounds cells in non-periodic ones
                    bool cell_is_valid = true;
                    for (int d = 0; d < 3; d++) {
                        if (box_.periodic(d)) {
                            auto [q, r] = divmod(neighbor_cell_i[d], cells_shape_[d]);
                            cell_shift[d] = q;
                            neighbor_cell_i[d] = r;
                        } else if (neighbor_cell_i[d] < 0 || neighbor_cell_i[d] >= static_cast<int32_t>(cells_shape_[d])) {
                            cell_is_valid = false;
                            break;
                        }
                    }

                    if (!cell_is_valid) {
                        continue;
                    }

                    for (const auto& atom_i : current_cell) {
                        for (const auto& atom_j : this->get_cell(neighbor_cell_i)) {
                            auto shift = CellShift(cell_shift) + atom_i.shift - atom_j.shift;
                            auto shift_is_zero = shift[0] == 0 && shift[1] == 0 && shift[2] == 0;

                            if ((shift[0] != 0 && !box_.periodic(0)) ||
                                (shift[1] != 0 && !box_.periodic(1)) ||
                                (shift[2] != 0 && !box_.periodic(2))) {
                                // do not create pairs crossing the periodic
                                // boundaries in a non-periodic box
                                continue;
                            }

                            if (atom_i.index == atom_j.index && shift_is_zero) {
                                // only create pairs with the same atom twice if the
                                // pair spans more than one bounding box
                                continue;
                            }

                            if constexpr (std::is_invocable_v<Function, size_t, size_t, CellShift, size_t>) {
                                callback(atom_i.index, atom_j.index, shift, thread_id);
                            } else {
                                callback(atom_i.index, atom_j.index, shift);
                            }
                        }
                    } // loop over atoms in current neighbor cells
                }
            }
        }
    });
}

CellList::Cell& CellList::get_cell(std::array<int32_t, 3> index) {
    size_t linear_index = (cells_shape_[0] * cells_shape_[1] * index[2]) + (cells_shape_[0] * index[1]) + index[0];
    return cells_[linear_index];
}

const CellList::Cell& CellList::get_cell(std::array<int32_t, 3> index) const {
    size_t linear_index = (cells_shape_[0] * cells_shape_[1] * index[2]) + (cells_shape_[0] * index[1]) + index[0];
    return cells_[linear_index];
}

/* ========================================================================== */

void GrowableNeighborList::set_pair(size_t index, size_t first, size_t second) {
    if (index >= this->capacity) {
        this->grow();
    }

    this->neighbors.pairs[index][0] = first;
    this->neighbors.pairs[index][1] = second;
}

void GrowableNeighborList::set_shift(size_t index, vesin::CellShift shift) {
    if (index >= this->capacity) {
        this->grow();
    }

    this->neighbors.shifts[index][0] = shift[0];
    this->neighbors.shifts[index][1] = shift[1];
    this->neighbors.shifts[index][2] = shift[2];
}

void GrowableNeighborList::set_distance(size_t index, double distance) {
    if (index >= this->capacity) {
        this->grow();
    }

    this->neighbors.distances[index] = distance;
}

void GrowableNeighborList::set_vector(size_t index, vesin::Vector vector) {
    if (index >= this->capacity) {
        this->grow();
    }

    this->neighbors.vectors[index][0] = vector[0];
    this->neighbors.vectors[index][1] = vector[1];
    this->neighbors.vectors[index][2] = vector[2];
}

template <typename scalar_t, size_t N>
using array_ptr = scalar_t (*)[N];

template <typename scalar_t, size_t N>
static array_ptr<scalar_t, N> alloc(array_ptr<scalar_t, N> ptr, size_t size, size_t new_size) {
    auto* new_ptr = reinterpret_cast<scalar_t(*)[N]>(std::realloc(ptr, new_size * sizeof(scalar_t[N])));

    if (new_ptr == nullptr) {
        return nullptr;
    }

#ifndef NDEBUG
    // initialize with a bit pattern that maps to NaN for double
    std::memset(new_ptr + size, 0b11111111, (new_size - size) * sizeof(scalar_t[N]));
#endif

    return new_ptr;
}

template <typename scalar_t>
static scalar_t* alloc(scalar_t* ptr, size_t size, size_t new_size) {
    auto* new_ptr = reinterpret_cast<scalar_t*>(std::realloc(ptr, new_size * sizeof(scalar_t)));

    if (new_ptr == nullptr) {
        return nullptr;
    }

#ifndef NDEBUG
    // initialize with a bit pattern that maps to NaN for double
    std::memset(new_ptr + size, 0b11111111, (new_size - size) * sizeof(scalar_t));
#endif

    return new_ptr;
}

void GrowableNeighborList::grow(size_t new_size) {
    assert(this->neighbors.device.type == VesinCPU);

    if (new_size == 0) {
        new_size = neighbors.length == 0 ? 1 : neighbors.length * 2;
    }

    auto* new_pairs = alloc<size_t, 2>(neighbors.pairs, neighbors.length, new_size);

    int32_t (*new_shifts)[3] = nullptr;
    if (options.return_shifts) {
        new_shifts = alloc<int32_t, 3>(neighbors.shifts, neighbors.length, new_size);
    }

    double* new_distances = nullptr;
    if (options.return_distances) {
        new_distances = alloc<double>(neighbors.distances, neighbors.length, new_size);
    }

    double (*new_vectors)[3] = nullptr;
    if (options.return_vectors) {
        new_vectors = alloc<double, 3>(neighbors.vectors, neighbors.length, new_size);
    }

    if (
        (new_pairs == nullptr) ||
        (options.return_shifts && new_shifts == nullptr) ||
        (options.return_distances && new_distances == nullptr) ||
        (options.return_vectors && new_vectors == nullptr)
    ) {
        std::free(new_pairs);
        std::free(new_shifts);
        std::free(new_distances);
        std::free(new_vectors);
        throw std::runtime_error("could not allocate memory for growing neighbor list");
    }

    this->neighbors.pairs = new_pairs;
    this->neighbors.shifts = new_shifts;
    this->neighbors.distances = new_distances;
    this->neighbors.vectors = new_vectors;

    this->capacity = new_size;
}

void GrowableNeighborList::reset() {
#ifndef NDEBUG
    auto size = this->neighbors.length;
    // set all allocated data to a bit pattern that maps to NaN for double
    std::memset(this->neighbors.pairs, 0b11111111, size * sizeof(size_t[2]));

    if (this->neighbors.shifts != nullptr) {
        std::memset(this->neighbors.shifts, 0b11111111, size * sizeof(int32_t[3]));
    }

    if (this->neighbors.distances != nullptr) {
        std::memset(this->neighbors.distances, 0b11111111, size * sizeof(double));
    }

    if (this->neighbors.vectors != nullptr) {
        std::memset(this->neighbors.vectors, 0b11111111, size * sizeof(double[3]));
    }
#endif

    // reset length (but keep the capacity where it's at)
    this->neighbors.length = 0;

    // allocate/deallocate pointers as required
    auto* shifts = this->neighbors.shifts;
    if (this->options.return_shifts && shifts == nullptr) {
        shifts = alloc<int32_t, 3>(shifts, 0, capacity);
    } else if (!this->options.return_shifts && shifts != nullptr) {
        std::free(shifts);
        shifts = nullptr;
    }

    auto* distances = this->neighbors.distances;
    if (this->options.return_distances && distances == nullptr) {
        distances = alloc<double>(distances, 0, capacity);
    } else if (!this->options.return_distances && distances != nullptr) {
        std::free(distances);
        distances = nullptr;
    }

    auto* vectors = this->neighbors.vectors;
    if (this->options.return_vectors && vectors == nullptr) {
        vectors = alloc<double, 3>(vectors, 0, capacity);
    } else if (!this->options.return_vectors && vectors != nullptr) {
        std::free(vectors);
        vectors = nullptr;
    }

    this->neighbors.shifts = shifts;
    this->neighbors.distances = distances;
    this->neighbors.vectors = vectors;
}

void GrowableNeighborList::sort() {
    if (this->length() == 0) {
        return;
    }

    // step 1: sort an array of indices, comparing the pairs at the indices
    auto indices = std::vector<int64_t>(this->length(), 0);
    std::iota(std::begin(indices), std::end(indices), 0);

    struct compare_pairs {
        compare_pairs(size_t (*pairs_)[2]):
            pairs(pairs_) {}

        bool operator()(int64_t a, int64_t b) const {
            return pairs[a][0] < pairs[b][0];
        }

        size_t (*pairs)[2];
    };

    std::sort(
        std::begin(indices),
        std::end(indices),
        compare_pairs(this->neighbors.pairs)
    );

    // step 2: move all data according to the sorted indices.
    auto* sorted_pairs = alloc<size_t, 2>(nullptr, 0, this->capacity);

    int32_t (*sorted_shifts)[3] = nullptr;
    if (options.return_shifts) {
        sorted_shifts = alloc<int32_t, 3>(nullptr, 0, this->capacity);
    }

    double* sorted_distances = nullptr;
    if (options.return_distances) {
        sorted_distances = alloc<double>(nullptr, 0, this->capacity);
    }

    double (*sorted_vectors)[3] = nullptr;
    if (options.return_vectors) {
        sorted_vectors = alloc<double, 3>(nullptr, 0, this->capacity);
    }

    if (
        (sorted_pairs == nullptr) ||
        (options.return_shifts && sorted_shifts == nullptr) ||
        (options.return_distances && sorted_distances == nullptr) ||
        (options.return_vectors && sorted_vectors == nullptr)
    ) {
        std::free(sorted_pairs);
        std::free(sorted_shifts);
        std::free(sorted_distances);
        std::free(sorted_vectors);
        throw std::runtime_error("could not allocate memory for sorting neighbor list");
    }

    for (size_t i = 0; i < this->neighbors.length; i++) {
        auto from = static_cast<size_t>(indices[i]);
        sorted_pairs[i][0] = this->neighbors.pairs[from][0];
        sorted_pairs[i][1] = this->neighbors.pairs[from][1];

        if (options.return_shifts) {
            sorted_shifts[i][0] = this->neighbors.shifts[from][0];
            sorted_shifts[i][1] = this->neighbors.shifts[from][1];
            sorted_shifts[i][2] = this->neighbors.shifts[from][2];
        }

        if (options.return_distances) {
            sorted_distances[i] = this->neighbors.distances[from];
        }

        if (options.return_vectors) {
            sorted_vectors[i][0] = this->neighbors.vectors[from][0];
            sorted_vectors[i][1] = this->neighbors.vectors[from][1];
            sorted_vectors[i][2] = this->neighbors.vectors[from][2];
        }
    }

    std::free(this->neighbors.pairs);
    this->neighbors.pairs = sorted_pairs;

    if (options.return_shifts) {
        std::free(this->neighbors.shifts);
        this->neighbors.shifts = sorted_shifts;
    }

    if (options.return_distances) {
        std::free(this->neighbors.distances);
        this->neighbors.distances = sorted_distances;
    }

    if (options.return_vectors) {
        std::free(this->neighbors.vectors);
        this->neighbors.vectors = sorted_vectors;
    }
}

void vesin::cpu::free_neighbors(VesinNeighborList& neighbors) {
    assert(neighbors.device.type == VesinCPU);

    if (neighbors.opaque != nullptr) {
        free_extra_data(neighbors);
    }

    std::free(neighbors.pairs);
    std::free(neighbors.shifts);
    std::free(neighbors.vectors);
    std::free(neighbors.distances);
}
#include <algorithm>
#include <cmath>
#include <cstring>


using namespace vesin;

// Sum of the two largest box-corner displacements, used to shrink the Verlet
// rebuild threshold when the box deforms (e.g. NPT).
static double corner_point_displacements(const Matrix& box, const Matrix& ref_box) {
    double delta1 = 0.0;
    double delta2 = 0.0;
    // The eight corners of the box are generated by the binary representation of the
    // index `i` (0-7). When a bit is 1, the corresponding box vector is included in the
    // corner position.
    for (size_t i = 0; i < 8; i++) {
        auto bits = Vector{
            static_cast<double>(i & 1),
            static_cast<double>((i >> 1) & 1),
            static_cast<double>((i >> 2) & 1),
        };
        auto displacement = bits * box - bits * ref_box;
        double delta = displacement.dot(displacement);
        if (delta > delta1) {
            delta2 = delta1;
            delta1 = delta;
        } else if (delta > delta2) {
            delta2 = delta;
        }
    }
    return std::sqrt(delta1) + std::sqrt(delta2);
}

cpu::VerletList::~VerletList() {
    this->clear_candidates();
}

void cpu::VerletList::clear_candidates() {
    if (candidates_.device.type == VesinCPU) {
        cpu::free_neighbors(candidates_);
    }

    candidates_ = VesinNeighborList();
    has_cache_ = false;
    ref_points_.clear();
}

void cpu::VerletList::set_options(VesinOptions options) {
    if (options_.cutoff != options.cutoff || options_.skin != options.skin || options_.full != options.full) {
        this->clear_candidates();
    }

    options_ = options;
    half_skin_sq_ = (options.skin / 2.0) * (options.skin / 2.0);
}

bool cpu::VerletList::needs_rebuild(
    const Vector* points,
    size_t n_points,
    const BoundingBox& box
) const {
    if (!has_cache_) {
        return true;
    }

    if (n_points != ref_points_.size()) {
        return true;
    }

    for (size_t d = 0; d < 3; d++) {
        if (box.periodic(d) != ref_periodic_[d]) {
            return true;
        }
    }
    // The corner displacement uses the full box matrix: non-periodic directions
    // carry no shift, so including them only over-estimates the bound (more
    // rebuilds, never fewer). This stays correct for any mix of periodicity.
    auto half_displacement = corner_point_displacements(box.matrix(), ref_matrix_) / 2.0;
    if (half_displacement * half_displacement > half_skin_sq_) {
        return true;
    }

    // Verlet-list reuse is valid while every atom stays within skin/2 of its
    // reference position: any pair inside cutoff is present in the cached
    // candidate list built at cutoff + skin. See Verlet, Phys. Rev. 159, 98-103
    // (1967), doi:10.1103/PhysRev.159.98, and Chialvo and Debenedetti, Comput.
    // Phys. Commun. 60, 215-224 (1990), doi:10.1016/0010-4655(90)90007-N.
    // One can also take the change of the box into account, see the
    // implementation of LAMMPS:
    // https://github.com/lammps/lammps/blob/3bfc12b02799eedf79d779d66fad8c4c60554084/src/neighbor.cpp#L2434-L2448
    auto half_threshold_sq = (sqrt(half_skin_sq_) - half_displacement) * (sqrt(half_skin_sq_) - half_displacement);
    for (size_t i = 0; i < n_points; i++) {
        auto displacement = points[i] - ref_points_[i];
        double displacement_sq = displacement.dot(displacement);
        if (displacement_sq > half_threshold_sq) {
            return true;
        }
    }

    return false;
}

void cpu::VerletList::rebuild(
    const Vector* points,
    size_t n_points,
    const BoundingBox& box
) {
    this->clear_candidates();

    auto build_options = VesinOptions();
    build_options.cutoff = options_.cutoff + options_.skin;
    build_options.full = options_.full;
    build_options.sorted = false;
    build_options.algorithm = VesinCellList;
    build_options.return_shifts = true;
    build_options.return_distances = false;
    build_options.return_vectors = false;
    build_options.skin = 0.0;

    candidates_.device = {VesinCPU, 0};

    auto periodic = std::array<bool, 3>{box.periodic(0), box.periodic(1), box.periodic(2)};
    auto candidate_box = BoundingBox(box.matrix(), periodic.data());
    candidate_box.make_bounding_for(reinterpret_cast<const double (*)[3]>(points), n_points);

    size_t candidate_capacity = 0;
    cpu::cell_list_neighbors(points, n_points, std::move(candidate_box), build_options, candidates_, candidate_capacity);

    ref_points_.resize(n_points);
    std::memcpy(ref_points_.data(), points, n_points * sizeof(Vector));
    ref_matrix_ = box.matrix();
    for (size_t d = 0; d < 3; d++) {
        ref_periodic_[d] = box.periodic(d);
    }

    has_cache_ = true;
}

void cpu::VerletList::filter(
    const Vector* points,
    const BoundingBox& box,
    VesinOptions options,
    VesinNeighborList& neighbors,
    size_t& output_capacity
) const {
    double cutoff_sq = options_.cutoff * options_.cutoff;

    auto initial_capacity = std::max(output_capacity, neighbors.length);

    auto growable = cpu::GrowableNeighborList{neighbors, initial_capacity, options};
    growable.reset();

    // The cached list is an over-complete Verlet candidate list. Each call
    // filters candidates with the exact cutoff and requested shift/vector outputs.
    for (size_t k = 0; k < candidates_.length; k++) {
        size_t i = candidates_.pairs[k][0];
        size_t j = candidates_.pairs[k][1];

        auto shift = CellShift(candidates_.shifts[k]);

        auto vec = points[j] - points[i] + shift.cartesian(box);
        double dist_sq = vec.dot(vec);

        if (dist_sq < cutoff_sq) {
            auto idx = growable.length();
            growable.set_pair(idx, i, j);

            if (options.return_shifts) {
                growable.set_shift(idx, shift);
            }

            if (options.return_distances) {
                growable.set_distance(idx, std::sqrt(dist_sq));
            }

            if (options.return_vectors) {
                growable.set_vector(idx, vec);
            }

            growable.increment_length();
        }
    }

    if (options.sorted) {
        growable.sort();
    }

    output_capacity = growable.capacity;
}
#include <stdexcept>

#ifndef VESIN_CUDA_HPP
#define VESIN_CUDA_HPP

#include <cstdint>

#ifndef VESIN_VERLET_CUDA_HPP
#define VESIN_VERLET_CUDA_HPP

#include <cstdint>


namespace vesin {
namespace cuda {

/// Verlet cache state. Candidates are generated at cutoff + skin and filtered
/// at the exact cutoff while the reference points remain valid.
class VerletCache {
public:
    VerletCache() = default;
    ~VerletCache();

    VerletCache(VerletCache&& other) noexcept;
    VerletCache& operator=(VerletCache&& other) noexcept;

    VerletCache(const VerletCache&) = delete;
    VerletCache& operator=(const VerletCache&) = delete;

    /// Run the verlet calculation, recomputing the cache as needed
    void run(
        const double (*points)[3],
        size_t n_points,
        const double box[3][3],
        const bool periodic[3],
        VesinOptions options,
        VesinNeighborList& neighbors
    );

private:
    // reference values for points/box/periodic, used to create the cache and
    // check its validity
    size_t ref_points_capacity_ = 0; // allocated capacity for d_ref_points
    double* d_ref_points_ = nullptr; // [n_ref_points * 3]
    size_t n_ref_points_ = 0;
    double ref_box_[9] = {0.0};
    bool ref_periodic_[3] = {false, false, false};

    VesinOptions options_;
    bool has_cache_ = false;
    VesinNeighborList candidates_;
    // cached allocation for rebuild_flag
    int32_t* d_rebuild_flag_ = nullptr;

    /// Free allocated buffers
    void free_buffers();

    /// Did the options changed since the cache was built?
    bool options_changed(VesinOptions options) const;
    /// Did the box size/shape change since the cache was built?
    bool box_size_changed(const double h_box[9]) const;
    /// Did the periodicity change since the cache was built?
    bool box_periodic_changed(const bool h_periodic[3]) const;
    /// Allocate the buffer for reference data
    void allocate_ref_buffers(size_t n_points);
    void rebuild_cache(const double (*points)[3], size_t n_points, const double box[3][3], const bool periodic[3], int32_t device_id, VesinOptions options);
    void filter_neighbors(const double (*d_points)[3], size_t n_points, const double d_box[3][3], VesinOptions options, VesinNeighborList& neighbors) const;
};

} // namespace cuda
} // namespace vesin

#endif // VESIN_CUDA_HPP

namespace vesin {
namespace cuda {

struct CudaNeighborListExtras;

#ifndef VESIN_CUDA_AT_LEAST_PAIRS_PER_POINT
/// Default value for the number of pairs per points in the CUDA implementation.
/// Unless `VESIN_CUDA_MAX_PAIRS_PER_POINT` is set in the environement, the
/// maximal number of pairs is `n_points *
/// max(VESIN_CUDA_AT_LEAST_PAIRS_PER_POINT, cutoff^3)`. This can be overriden
/// at compile time.
#define VESIN_CUDA_AT_LEAST_PAIRS_PER_POINT 128
#endif

/// @brief Buffers for cell list-based neighbor search
struct CellListBuffers {
    size_t points_capacity = 0; // Capacity for point-related arrays
    size_t cells_capacity = 0;  // Capacity for cell-related arrays

    // Per-particle arrays (device)
    int32_t* d_cell_indices = nullptr;    // [points_capacity] linear cell index per particle
    int32_t* d_particle_shifts = nullptr; // [points_capacity * 3] shift applied to wrap into cell

    // Per-cell arrays (device)
    int32_t* d_cell_counts = nullptr;  // [cells_capacity] number of particles in each cell
    int32_t* d_cell_starts = nullptr;  // [cells_capacity] starting index in sorted arrays
    int32_t* d_cell_offsets = nullptr; // [cells_capacity] working copy for scatter

    // Sorted particle data (device, for coalesced memory access)
    double* d_sorted_points = nullptr;        // [points_capacity * 3]
    int32_t* d_sorted_indices = nullptr;      // [points_capacity] original particle indices
    int32_t* d_sorted_shifts = nullptr;       // [points_capacity * 3] shifts for sorted particles
    int32_t* d_sorted_cell_indices = nullptr; // [points_capacity] cell indices in sorted order

    // Cell grid parameters (device, computed on device)
    int32_t* d_n_cells = nullptr;       // [3] number of cells in each direction
    int32_t* d_n_search = nullptr;      // [3] search range in each direction
    int32_t* d_n_cells_total = nullptr; // [1] total number of cells

    double* d_face_distances = nullptr; // [3] distances between faces of the box
    double* d_bounding_min = nullptr;   // [3] bottom of the bounding box

    void allocate(size_t n_points, size_t n_cells);

    CellListBuffers() = default;
    ~CellListBuffers();

    CellListBuffers(CellListBuffers&& other) noexcept;
    CellListBuffers& operator=(CellListBuffers&& other) noexcept;

    CellListBuffers(const CellListBuffers&) = delete;
    CellListBuffers& operator=(const CellListBuffers&) = delete;
};

struct SortBuffers {
    size_t capacity = 0;                  // Capacity for the buffers below (number of pairs)
    size_t (*d_pairs_tmp)[2] = nullptr;   // [capacity] temporary pair indices for sorting
    int32_t (*d_shifts_tmp)[3] = nullptr; // [capacity] temporary shifts for sorting
    double* d_distances_tmp = nullptr;    // [capacity] temporary distances for sorting
    double (*d_vectors_tmp)[3] = nullptr; // [capacity * 3] temporary vectors for sorting

    void allocate(size_t n, bool return_shifts, bool return_distances, bool return_vectors);

    SortBuffers() = default;
    ~SortBuffers();

    SortBuffers(SortBuffers&& other) noexcept;
    SortBuffers& operator=(SortBuffers&& other) noexcept;

    SortBuffers(const SortBuffers&) = delete;
    SortBuffers& operator=(const SortBuffers&) = delete;
};

struct CudaNeighborListExtras {
    size_t pairs_capacity = 0;           // Capacity for the pair buffers in the VesinNeighborList
    size_t* d_length_ptr = nullptr;      // GPU-side counter
    int32_t* d_cell_check_ptr = nullptr; // GPU-side status code for checking cell
    int32_t* d_overflow_flag = nullptr;  // GPU-side flag to detect overflow of pair buffers

    // Pinned host memory for async D2H copy
    size_t* pinned_length_ptr = nullptr;

    // Cell list buffers (allocated on demand for large systems)
    CellListBuffers cell_list;

    // Buffers for optimized brute force kernels
    double* d_box_diag = nullptr;     // [3] diagonal elements for orthogonal boxes
    double (*d_inv_box)[3] = nullptr; // [3][3] inverse box matrix for general boxes

    // Temporary buffers for on-device sorting
    SortBuffers sort_buffers;

    // Verlet cache state
    VerletCache verlet_cache;

    CudaNeighborListExtras() = default;
    ~CudaNeighborListExtras();

    CudaNeighborListExtras(CudaNeighborListExtras&& other) noexcept;
    CudaNeighborListExtras& operator=(CudaNeighborListExtras&& other) noexcept;

    CudaNeighborListExtras(const CudaNeighborListExtras&) = delete;
    CudaNeighborListExtras& operator=(const CudaNeighborListExtras&) = delete;
};

/// @brief Frees GPU memory associated with a VesinNeighborList.
///
/// This function should be called to release all CUDA-allocated memory
/// tied to the given neighbor list. It does not delete the structure itself,
/// only the device-side memory buffers.
///
/// @param neighbors Reference to the VesinNeighborList to clean up.
void free_neighbors(VesinNeighborList& neighbors);

/// @brief Computes the neighbor list on the GPU.
///
/// This function only works under Minimum Image Convention for now.
///
/// This function generates a neighbor list for a set of points within a
/// periodic simulation box using GPU acceleration. The output is stored in a
/// `VesinNeighborList` structure, which must be initialized for GPU usage.
///
/// @param points Pointer to an array of 3D points (shape: [n_points][3]).
/// @param n_points Number of points (atoms, particles, etc.).
/// @param box 3×3 matrix defining the bounding box of the system.
/// @param periodic Array of three booleans indicating periodicity in each dimension.
/// @param options Struct holding parameters such as cutoff, symmetry, etc.
/// @param neighbors Output neighbor list (device memory will be allocated as
/// needed).
void neighbors(
    const double (*points)[3],
    size_t n_points,
    const double box[3][3],
    const bool periodic[3],
    VesinOptions options,
    VesinNeighborList& neighbors
);

/// Get the `CudaNeighborListExtras` stored inside `VesinNeighborList`'s opaque pointer
inline CudaNeighborListExtras* get_cuda_extras(VesinNeighborList& neighbors) {
    if (neighbors.opaque == nullptr) {
        neighbors.opaque = new vesin::cuda::CudaNeighborListExtras();
    }
    return static_cast<vesin::cuda::CudaNeighborListExtras*>(neighbors.opaque);
}

/// Allocate output buffers in the `VesinNeighborList` according to the options
/// and the given number of pairs. If the current capacity is sufficient, this
/// is a no-op. Otherwise, existing buffers are freed and new ones are allocated
/// with the requested capacity.
void allocate_output_buffers(VesinNeighborList& neighbors, size_t n_pairs, VesinOptions options);

} // namespace cuda
} // namespace vesin

#endif // VESIN_CUDA_HPP

void vesin::cuda::free_neighbors(VesinNeighborList& neighbors) {
    throw std::runtime_error("CUDA neighbor list generation is not included in this build of vesin");
}

void vesin::cuda::neighbors(
    const double (*points)[3],
    size_t n_points,
    const double box[3][3],
    const bool periodic[3],
    VesinOptions options,
    VesinNeighborList& neighbors
) {
    throw std::runtime_error("CUDA neighbor list generation is not included in this build of vesin");
}
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <string>
#include <thread>


// used to store dynamically allocated error messages before giving a pointer
// to them back to the user
thread_local std::string LAST_ERROR;

static size_t resolve_n_threads(int32_t n_threads) {
    if (n_threads > 0) {
        return static_cast<size_t>(n_threads);
    }

    const char* omp_num_threads = std::getenv("OMP_NUM_THREADS");
    if (omp_num_threads != nullptr) {
        errno = 0;
        char* end = nullptr;
        auto parsed = std::strtol(omp_num_threads, &end, 10);
        if (errno == 0 && end != omp_num_threads && *end == '\0' && parsed > 0 && parsed <= INT32_MAX) {
            return static_cast<size_t>(parsed);
        }
    }

    auto hardware_threads = std::thread::hardware_concurrency();
    if (hardware_threads > 0) {
        return static_cast<size_t>(hardware_threads);
    }

    return 1;
}

extern "C" int vesin_neighbors(
    const double (*points)[3],
    size_t n_points,
    const double box[3][3],
    const bool periodic[3],
    VesinDevice device,
    VesinOptions options,
    VesinNeighborList* neighbors,
    const char** error_message
) {
    if (error_message == nullptr) {
        return EXIT_FAILURE;
    }

    if (points == nullptr) {
        *error_message = "`points` can not be a NULL pointer";
        return EXIT_FAILURE;
    }

    if (box == nullptr) {
        *error_message = "`cell` can not be a NULL pointer";
        return EXIT_FAILURE;
    }

    if (neighbors == nullptr) {
        *error_message = "`neighbors` can not be a NULL pointer";
        return EXIT_FAILURE;
    }

    if (!std::isfinite(options.cutoff) || options.cutoff <= 0) {
        *error_message = "cutoff must be a finite, positive number";
        return EXIT_FAILURE;
    }

    if (options.cutoff <= 1e-6) {
        *error_message = "cutoff is too small";
        return EXIT_FAILURE;
    }

    if (!std::isfinite(options.skin) || options.skin < 0.0) {
        *error_message = "skin must be a finite, non-negative number";
    }

    if (options.n_threads < 0) {
        *error_message = "n_threads must be zero or a positive integer";
        return EXIT_FAILURE;
    }

    if (neighbors->device.type != VesinUnknownDevice && neighbors->device.type != device.type) {
        *error_message = "`neighbors` device and data `device` do not match, free the neighbors first";
        return EXIT_FAILURE;
    }

    if (device.type == VesinUnknownDevice) {
        *error_message = "got an unknown device type";
        return EXIT_FAILURE;
    }

    if (neighbors->device.type == VesinUnknownDevice) {
        // initialize the device
        neighbors->device = device;
    } else if (neighbors->device.type != device.type) {
        *error_message = "`neighbors.device` and `device` do not match, free the neighbors first";
        return EXIT_FAILURE;
    }

    try {
        options.n_threads = static_cast<int32_t>(resolve_n_threads(options.n_threads));

        if (device.type == VesinCPU) {
            auto matrix = vesin::Matrix{{{
                {{box[0][0], box[0][1], box[0][2]}},
                {{box[1][0], box[1][1], box[1][2]}},
                {{box[2][0], box[2][1], box[2][2]}},
            }}};

            auto box = vesin::BoundingBox(matrix, periodic);
            box.make_bounding_for(points, n_points);

            vesin::cpu::neighbors(
                reinterpret_cast<const vesin::Vector*>(points),
                n_points,
                std::move(box),
                options,
                *neighbors
            );
        } else if (device.type == VesinCUDA) {
            vesin::cuda::neighbors(
                points,
                n_points,
                box,
                periodic,
                options,
                *neighbors
            );
        } else {
            throw std::runtime_error("unknown device " + std::to_string(device.type));
        }
    } catch (const std::bad_alloc&) {
        LAST_ERROR = "failed to allocate memory";
        *error_message = LAST_ERROR.c_str();
        return EXIT_FAILURE;
    } catch (const std::exception& e) {
        LAST_ERROR = e.what();
        *error_message = LAST_ERROR.c_str();
        return EXIT_FAILURE;
    } catch (...) {
        *error_message = "fatal error: unknown type thrown as exception";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

// eOn extension (upstream candidate); see vesin.h.
extern "C" int vesin_neighbors_visit(
    const double (*points)[3],
    size_t n_points,
    const double box[3][3],
    const bool periodic[3],
    VesinDevice device,
    VesinOptions options,
    double visit_cutoff,
    VesinPairVisitor visitor,
    void* user_data,
    VesinNeighborList* neighbors,
    const char** error_message
) {
    if (error_message == nullptr) {
        return EXIT_FAILURE;
    }

    if (points == nullptr || box == nullptr || neighbors == nullptr || visitor == nullptr) {
        *error_message = "`points`, `box`, `neighbors` and `visitor` can not be NULL pointers";
        return EXIT_FAILURE;
    }

    if (!std::isfinite(options.cutoff) || options.cutoff <= 1e-6) {
        *error_message = "cutoff must be a finite, positive number";
        return EXIT_FAILURE;
    }

    if (!std::isfinite(visit_cutoff) || visit_cutoff <= 0 || visit_cutoff > options.cutoff) {
        *error_message = "visit_cutoff must be positive and not larger than options.cutoff";
        return EXIT_FAILURE;
    }

    if (device.type != VesinCPU) {
        *error_message = "vesin_neighbors_visit is CPU-only";
        return EXIT_FAILURE;
    }

    if (neighbors->device.type == VesinUnknownDevice) {
        neighbors->device = device;
    } else if (neighbors->device.type != device.type) {
        *error_message = "`neighbors.device` and `device` do not match, free the neighbors first";
        return EXIT_FAILURE;
    }

    try {
        options.algorithm = VesinBruteForce;
        options.skin = 0.0; // the caller owns skin semantics in this mode

        auto matrix = vesin::Matrix{{{
            {{box[0][0], box[0][1], box[0][2]}},
            {{box[1][0], box[1][1], box[1][2]}},
            {{box[2][0], box[2][1], box[2][2]}},
        }}};

        auto bounding = vesin::BoundingBox(matrix, periodic);
        bounding.make_bounding_for(points, n_points);

        auto& extra = extra_data(*neighbors);
        delete extra.verlet_state;
        extra.verlet_state = nullptr;

        vesin::cpu::brute_force_neighbors(
            reinterpret_cast<const vesin::Vector*>(points),
            n_points,
            std::move(bounding),
            options,
            *neighbors,
            extra.capacity,
            visitor,
            user_data,
            visit_cutoff
        );
    } catch (const std::bad_alloc&) {
        LAST_ERROR = "failed to allocate memory";
        *error_message = LAST_ERROR.c_str();
        return EXIT_FAILURE;
    } catch (const std::exception& e) {
        LAST_ERROR = e.what();
        *error_message = LAST_ERROR.c_str();
        return EXIT_FAILURE;
    } catch (...) {
        *error_message = "fatal error: unknown type thrown as exception";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

extern "C" void vesin_free(VesinNeighborList* neighbors) {
    if (neighbors == nullptr) {
        return;
    }

    try {
        if (neighbors->device.type == VesinUnknownDevice) {
            // nothing to do
        } else if (neighbors->device.type == VesinCPU) {
            vesin::cpu::free_neighbors(*neighbors);
        } else if (neighbors->device.type == VesinCUDA) {
            vesin::cuda::free_neighbors(*neighbors);
        } else {
            throw std::runtime_error("unknown device " + std::to_string(neighbors->device.type) + " when freeing memory");
        }
    } catch (const std::exception& e) {
        std::cerr << "error in vesin_free: " << e.what() << std::endl;
        return;
    } catch (...) {
        std::cerr << "fatal error in vesin_free, unknown type thrown as exception" << std::endl;
        return;
    }

    std::memset(neighbors, 0, sizeof(VesinNeighborList));
}
