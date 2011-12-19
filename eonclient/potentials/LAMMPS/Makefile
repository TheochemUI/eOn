OBJECTS = LAMMPS_EON.o
MPI_OBJECTS = mpi.o
#CXX =		g++
CXXFLAGS += -I../..
#LINK =		g++

all: liblammps.a libfakempi.a

LIB := liblammps.a
$(LIB) : $(OBJECTS)
	ar cru $@ $^
	ranlib $@

MPI_LIB := libfakempi.a
$(MPI_LIB) : $(MPI_OBJECTS)
	ar cru $@ $^
	ranlib $@

clean :
	rm -f $(OBJECTS) $(LIB) $(MPI_OBJECTS) $(MPI_LIB)

.PHONY : clean
LAMMPS_EON.o: LAMMPS_EON.cpp LAMMPS_EON.h ../../Potential.h \
  ../../Parameters.h ../../Eigen.h ../../Eigen/Dense ../../Eigen/Core \
  ../../Eigen/src/Core/util/DisableStupidWarnings.h \
  ../../Eigen/src/Core/util/Macros.h \
  ../../Eigen/src/Core/util/Constants.h \
  ../../Eigen/src/Core/util/ForwardDeclarations.h \
  ../../Eigen/src/Core/util/Meta.h ../../Eigen/src/Core/util/XprHelper.h \
  ../../Eigen/src/Core/util/StaticAssert.h \
  ../../Eigen/src/Core/util/Memory.h ../../Eigen/src/Core/NumTraits.h \
  ../../Eigen/src/Core/MathFunctions.h \
  ../../Eigen/src/Core/GenericPacketMath.h \
  ../../Eigen/src/Core/arch/Default/Settings.h \
  ../../Eigen/src/Core/Functors.h ../../Eigen/src/Core/DenseCoeffsBase.h \
  ../../Eigen/src/Core/DenseBase.h \
  ../../Eigen/src/Core/../plugins/BlockMethods.h \
  ../../Eigen/src/Core/MatrixBase.h \
  ../../Eigen/src/Core/../plugins/CommonCwiseUnaryOps.h \
  ../../Eigen/src/Core/../plugins/CommonCwiseBinaryOps.h \
  ../../Eigen/src/Core/../plugins/MatrixCwiseUnaryOps.h \
  ../../Eigen/src/Core/../plugins/MatrixCwiseBinaryOps.h \
  ../../Eigen/src/Core/EigenBase.h ../../Eigen/src/Core/Assign.h \
  ../../Eigen/src/Core/util/BlasUtil.h \
  ../../Eigen/src/Core/DenseStorage.h ../../Eigen/src/Core/NestByValue.h \
  ../../Eigen/src/Core/ForceAlignedAccess.h \
  ../../Eigen/src/Core/ReturnByValue.h ../../Eigen/src/Core/NoAlias.h \
  ../../Eigen/src/Core/PlainObjectBase.h ../../Eigen/src/Core/Matrix.h \
  ../../Eigen/src/Core/Array.h ../../Eigen/src/Core/CwiseBinaryOp.h \
  ../../Eigen/src/Core/CwiseUnaryOp.h \
  ../../Eigen/src/Core/CwiseNullaryOp.h \
  ../../Eigen/src/Core/CwiseUnaryView.h \
  ../../Eigen/src/Core/SelfCwiseBinaryOp.h ../../Eigen/src/Core/Dot.h \
  ../../Eigen/src/Core/StableNorm.h ../../Eigen/src/Core/MapBase.h \
  ../../Eigen/src/Core/Stride.h ../../Eigen/src/Core/Map.h \
  ../../Eigen/src/Core/Block.h ../../Eigen/src/Core/VectorBlock.h \
  ../../Eigen/src/Core/Transpose.h ../../Eigen/src/Core/DiagonalMatrix.h \
  ../../Eigen/src/Core/Diagonal.h ../../Eigen/src/Core/DiagonalProduct.h \
  ../../Eigen/src/Core/PermutationMatrix.h \
  ../../Eigen/src/Core/Transpositions.h ../../Eigen/src/Core/Redux.h \
  ../../Eigen/src/Core/Visitor.h ../../Eigen/src/Core/Fuzzy.h \
  ../../Eigen/src/Core/IO.h ../../Eigen/src/Core/Swap.h \
  ../../Eigen/src/Core/CommaInitializer.h ../../Eigen/src/Core/Flagged.h \
  ../../Eigen/src/Core/ProductBase.h ../../Eigen/src/Core/Product.h \
  ../../Eigen/src/Core/TriangularMatrix.h \
  ../../Eigen/src/Core/SelfAdjointView.h \
  ../../Eigen/src/Core/SolveTriangular.h \
  ../../Eigen/src/Core/products/Parallelizer.h \
  ../../Eigen/src/Core/products/CoeffBasedProduct.h \
  ../../Eigen/src/Core/products/GeneralBlockPanelKernel.h \
  ../../Eigen/src/Core/products/GeneralMatrixVector.h \
  ../../Eigen/src/Core/products/GeneralMatrixMatrix.h \
  ../../Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h \
  ../../Eigen/src/Core/products/SelfadjointMatrixVector.h \
  ../../Eigen/src/Core/products/SelfadjointMatrixMatrix.h \
  ../../Eigen/src/Core/products/SelfadjointProduct.h \
  ../../Eigen/src/Core/products/SelfadjointRank2Update.h \
  ../../Eigen/src/Core/products/TriangularMatrixVector.h \
  ../../Eigen/src/Core/products/TriangularMatrixMatrix.h \
  ../../Eigen/src/Core/products/TriangularSolverMatrix.h \
  ../../Eigen/src/Core/products/TriangularSolverVector.h \
  ../../Eigen/src/Core/BandMatrix.h ../../Eigen/src/Core/BooleanRedux.h \
  ../../Eigen/src/Core/Select.h ../../Eigen/src/Core/VectorwiseOp.h \
  ../../Eigen/src/Core/Random.h ../../Eigen/src/Core/Replicate.h \
  ../../Eigen/src/Core/Reverse.h ../../Eigen/src/Core/ArrayBase.h \
  ../../Eigen/src/Core/../plugins/ArrayCwiseUnaryOps.h \
  ../../Eigen/src/Core/../plugins/ArrayCwiseBinaryOps.h \
  ../../Eigen/src/Core/ArrayWrapper.h \
  ../../Eigen/src/Core/GlobalFunctions.h \
  ../../Eigen/src/Core/util/ReenableStupidWarnings.h \
  ../../Eigen/Eigen2Support ../../Eigen/src/Eigen2Support/Macros.h \
  ../../Eigen/src/Eigen2Support/Memory.h \
  ../../Eigen/src/Eigen2Support/Meta.h \
  ../../Eigen/src/Eigen2Support/Lazy.h \
  ../../Eigen/src/Eigen2Support/Cwise.h \
  ../../Eigen/src/Eigen2Support/CwiseOperators.h \
  ../../Eigen/src/Eigen2Support/TriangularSolver.h \
  ../../Eigen/src/Eigen2Support/Block.h \
  ../../Eigen/src/Eigen2Support/VectorBlock.h \
  ../../Eigen/src/Eigen2Support/Minor.h \
  ../../Eigen/src/Eigen2Support/MathFunctions.h ../../Eigen/LU \
  ../../Eigen/src/misc/Solve.h ../../Eigen/src/misc/Kernel.h \
  ../../Eigen/src/misc/Image.h ../../Eigen/src/LU/FullPivLU.h \
  ../../Eigen/src/LU/PartialPivLU.h ../../Eigen/src/LU/Determinant.h \
  ../../Eigen/src/LU/Inverse.h ../../Eigen/src/Eigen2Support/LU.h \
  ../../Eigen/Cholesky ../../Eigen/src/Cholesky/LLT.h \
  ../../Eigen/src/Cholesky/LDLT.h ../../Eigen/QR ../../Eigen/Jacobi \
  ../../Eigen/src/Jacobi/Jacobi.h ../../Eigen/Householder \
  ../../Eigen/src/Householder/Householder.h \
  ../../Eigen/src/Householder/HouseholderSequence.h \
  ../../Eigen/src/Householder/BlockHouseholder.h \
  ../../Eigen/src/QR/HouseholderQR.h \
  ../../Eigen/src/QR/FullPivHouseholderQR.h \
  ../../Eigen/src/QR/ColPivHouseholderQR.h \
  ../../Eigen/src/Eigen2Support/QR.h ../../Eigen/Eigenvalues \
  ../../Eigen/src/Eigenvalues/Tridiagonalization.h \
  ../../Eigen/src/Eigenvalues/RealSchur.h \
  ../../Eigen/src/Eigenvalues/./EigenvaluesCommon.h \
  ../../Eigen/src/Eigenvalues/./HessenbergDecomposition.h \
  ../../Eigen/src/Eigenvalues/EigenSolver.h \
  ../../Eigen/src/Eigenvalues/./RealSchur.h \
  ../../Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h \
  ../../Eigen/src/Eigenvalues/./Tridiagonalization.h \
  ../../Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h \
  ../../Eigen/src/Eigenvalues/HessenbergDecomposition.h \
  ../../Eigen/src/Eigenvalues/ComplexSchur.h \
  ../../Eigen/src/Eigenvalues/ComplexEigenSolver.h \
  ../../Eigen/src/Eigenvalues/./ComplexSchur.h \
  ../../Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h ../../Eigen/SVD \
  ../../Eigen/src/SVD/JacobiSVD.h \
  ../../Eigen/src/SVD/UpperBidiagonalization.h \
  ../../Eigen/src/Eigen2Support/SVD.h ../../Eigen/Geometry \
  ../../Eigen/src/Geometry/OrthoMethods.h \
  ../../Eigen/src/Geometry/EulerAngles.h \
  ../../Eigen/src/Geometry/Homogeneous.h \
  ../../Eigen/src/Geometry/RotationBase.h \
  ../../Eigen/src/Geometry/Rotation2D.h \
  ../../Eigen/src/Geometry/Quaternion.h \
  ../../Eigen/src/Geometry/AngleAxis.h \
  ../../Eigen/src/Geometry/Transform.h \
  ../../Eigen/src/Geometry/Translation.h \
  ../../Eigen/src/Geometry/Scaling.h \
  ../../Eigen/src/Geometry/Hyperplane.h \
  ../../Eigen/src/Geometry/ParametrizedLine.h \
  ../../Eigen/src/Geometry/AlignedBox.h \
  ../../Eigen/src/Geometry/Umeyama.h \
  ../../Eigen/src/Eigen2Support/Geometry/All.h \
  ../../Eigen/src/Eigen2Support/Geometry/RotationBase.h \
  ../../Eigen/src/Eigen2Support/Geometry/Rotation2D.h \
  ../../Eigen/src/Eigen2Support/Geometry/Quaternion.h \
  ../../Eigen/src/Eigen2Support/Geometry/AngleAxis.h \
  ../../Eigen/src/Eigen2Support/Geometry/Transform.h \
  ../../Eigen/src/Eigen2Support/Geometry/Translation.h \
  ../../Eigen/src/Eigen2Support/Geometry/Scaling.h \
  ../../Eigen/src/Eigen2Support/Geometry/AlignedBox.h \
  ../../Eigen/src/Eigen2Support/Geometry/Hyperplane.h \
  ../../Eigen/src/Eigen2Support/Geometry/ParametrizedLine.h \
  ../../Eigen/Eigenvalues ../../Parameters.h library.h
