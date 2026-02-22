#include "../MatrixHelpers.hpp"
#include "../Eigen.h"
#include "catch2/catch_amalgamated.hpp"

using namespace Catch::Matchers;

namespace tests {

TEST_CASE("Eigen Storage Order and Mapping", "[Core][Eigen]") {
    // 2 atoms, (1,2,3) and (4,5,6)
    // RowMajor memory: [1, 2, 3, 4, 5, 6]
    double data[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    SECTION("AtomMatrix Map (RowMajor)") {
        // AtomMatrix is typedef'd as RowMajor in Eigen.h
        AtomMatrix mapped = AtomMatrix::Map(data, 2, 3);

        REQUIRE(mapped.rows() == 2);
        REQUIRE(mapped.cols() == 3);

        // Check values
        REQUIRE(mapped(0, 0) == 1.0);
        REQUIRE(mapped(0, 1) == 2.0);
        REQUIRE(mapped(0, 2) == 3.0);
        REQUIRE(mapped(1, 0) == 4.0);
        REQUIRE(mapped(1, 1) == 5.0);
        REQUIRE(mapped(1, 2) == 6.0);
    }

    SECTION("RotationMatrix Map (RowMajor)") {
        double rotData[] = {1.0, 2.0, 3.0,
                            4.0, 5.0, 6.0,
                            7.0, 8.0, 9.0};
        RotationMatrix mapped = RotationMatrix::Map(rotData, 3, 3);

        REQUIRE(mapped.rows() == 3);
        REQUIRE(mapped.cols() == 3);

        REQUIRE(mapped(0, 0) == 1.0);
        REQUIRE(mapped(0, 1) == 2.0);
        REQUIRE(mapped(0, 2) == 3.0);
        REQUIRE(mapped(1, 0) == 4.0);
        REQUIRE(mapped(1, 1) == 5.0);
        REQUIRE(mapped(1, 2) == 6.0);
        REQUIRE(mapped(2, 0) == 7.0);
        REQUIRE(mapped(2, 1) == 8.0);
        REQUIRE(mapped(2, 2) == 9.0);
    }

    SECTION("MatrixXd Map (Default ColMajor) - DEMONSTRATING THE BUG") {
        // This is what caused the regression. MatrixXd defaults to ColMajor.
        // On RowMajor data [1, 2, 3, 4, 5, 6]:
        // Col 0: 1, 2
        // Col 1: 3, 4
        // Col 2: 5, 6
        // Result:
        // 1 3 5
        // 2 4 6
        Eigen::MatrixXd mapped = Eigen::MatrixXd::Map(data, 2, 3);

        REQUIRE(mapped(0, 0) == 1.0);
        REQUIRE(mapped(0, 1) == 3.0); // NOT 2.0
        REQUIRE(mapped(0, 2) == 5.0); // NOT 3.0
        REQUIRE(mapped(1, 0) == 2.0); // NOT 4.0
        REQUIRE(mapped(1, 1) == 4.0); // NOT 5.0
        REQUIRE(mapped(1, 2) == 6.0);
    }
}

} // namespace tests
