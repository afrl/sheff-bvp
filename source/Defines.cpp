#include "fea_analysis/Defines.h"

#include "sparse_linear_algebra/LinAlgFactory.h"
#include "fea_io/common/MPIUtils.h"

namespace fea_analysis
{
    void init(const std::string preferred_lin_alg_backend, int* argc, char*** argv)
    {
        utility::MPIUtils::init(argc, argv);
        fea_io::MPIUtils::init(argc, argv);

        auto backend = sparse_linear_algebra::LinAlgFactory::backend_from_str(preferred_lin_alg_backend);
        sparse_linear_algebra::LinAlgFactory::init(backend, argc, argv);
    }
}