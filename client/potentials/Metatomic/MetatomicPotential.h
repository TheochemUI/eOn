#ifndef METATOMIC_POTENTIAL_H
#define METATOMIC_POTENTIAL_H

#include "../../Potential.h"

// Metatomic and torch headers
// These pragmas are included to suppress warnings from third-party libraries
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"

#include <torch/script.h>
#include <torch/version.h>
#include <torch/cuda.h>
#include <torch/mps.h>

#include "metatensor/torch.hpp"
#include "metatomic/torch.hpp"

#pragma GCC diagnostic pop


/**
 * @class MetatomicPotential
 * @brief A potential class that uses a metatomic model for energy and force calculations.
 *
 * This class loads a pre-trained atomistic model (via metatomic/PyTorch)
 * and uses it to compute the potential energy and corresponding forces on atoms.
 * It uses the `vesin` library to compute neighbor lists required by the model.
 *
 */
class MetatomicPotential : public Potential {
private:
    // --- Metatomic and Torch members ---
    torch::jit::Module model_;
    metatomic_torch::ModelCapabilities capabilities_;
    std::vector<metatomic_torch::NeighborListOptions> nl_requests_;
    metatomic_torch::ModelEvaluationOptions evaluations_options_;

    torch::ScalarType dtype_;
    torch::Device device_;
    bool check_consistency_;
    std::shared_ptr<Parameters> m_params;

    // --- Cached Tensors and Data for Performance ---
    torch::Tensor atomic_types_;
    std::vector<int> last_atomic_nrs_;

    /**
     * @brief Computes neighbor list using the vesin library.
     *
     * This function calls `vesin_neighbors` to build a neighbor list based on the
     * model's requirements (cutoff, full/half list) and converts it into a
     * metatensor::TensorBlock suitable for metatomic.
     *
     * @param request The neighbor list options requested by the model.
     * @param nAtoms The number of atoms in the system.
     * @param positions Pointer to the atomic positions array.
     * @param box Pointer to the simulation box matrix.
     * @return A metatensor_torch::TensorBlock containing the neighbor list.
     */
    metatensor_torch::TensorBlock computeNeighbors(
        metatomic_torch::NeighborListOptions request,
        long nAtoms,
        const double* positions,
        const double* box,
        bool periodic
    );

public:
    /**
     * @brief Constructor for the MetatomicPotential.
     *
     * @param params A shared pointer to the simulation parameters object.
     * This object should contain the settings needed for the metatomic model.
     */
    MetatomicPotential(std::shared_ptr<Parameters> params);

    /**
     * @brief Destructor.
     */
    ~MetatomicPotential() override = default;

    /**
     * @brief Calculates the energy and forces for a given atomic configuration.
     *
     * This is the core method of the potential. It takes the current atomic positions,
     * builds the necessary data structures for metatomic, executes the model to get the
     * potential energy, and uses PyTorch's autograd to compute forces.
     *
     * @param nAtoms Number of atoms.
     * @param positions Flat array of atomic positions (size nAtoms * 3).
     * @param atomicNrs Flat array of atomic numbers (size nAtoms), used for consistency checks.
     * @param forces Flat array to store the calculated forces (size nAtoms * 3).
     * @param energy Pointer to a double to store the calculated potential energy.
     * @param variance Pointer to a double to store the variance of the energy (currently unused, set to NULL).
     * @param box The simulation box vectors (3x3 matrix).
     */
    void force(long nAtoms, const double *positions, const int *atomicNrs,
               double *forces, double *energy, double *variance,
               const double *box) override;
};

#endif // METATOMIC_POTENTIAL_H
