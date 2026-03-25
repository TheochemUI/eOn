# Build eonclient from a git tag
# Creates a worktree, builds with make, produces the binary

import os

rule build_eonclient:
    output:
        binary=f"build/{GIT_TAG}/client/eonclient"
    params:
        tag=GIT_TAG,
        eigen_include=os.environ.get("CONDA_PREFIX", "/usr") + "/include/eigen3"
    log:
        f"build/{GIT_TAG}/build.log"
    shell:
        """
        set -e
        BUILDDIR="build/{params.tag}"

        # Create worktree from tag
        git worktree add "$BUILDDIR" {params.tag} 2>/dev/null || true

        cd "$BUILDDIR/client"

        # Create QSC stub if missing
        if [ ! -f potentials/QSC/QSC.h ]; then
            mkdir -p potentials/QSC
            cat > potentials/QSC/QSC.h << 'QEOF'
#pragma once
#include "../../Potential.h"
class QSC : public Potential {{ public: void initialize() override {{}}
void force(long n, const double *r, const int *a, double *f,
double *e, const double *b) override {{ *e=0; for(long i=0;i<3*n;i++) f[i]=0; }}
}};
QEOF
            cat > potentials/QSC/Makefile << 'MEOF'
all: libQSC.a
libQSC.a:
	ar cru libQSC.a
	ranlib libQSC.a
clean:
	rm -f *.a
MEOF
        fi

        NO_FORTRAN=1 CXXFLAGS="-O2 -I{params.eigen_include}" \
            make -j4 eonclient > ../../build.log 2>&1
        """
