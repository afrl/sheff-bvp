#include "bvp/Defines.h"

#include "fea_io/common/MPIUtils.h"
#include "sparse_linear_algebra/LinAlgFactory.h"
#include "utility/MPIUtils.h"

namespace bvp
{
/*
void init(const std::string preferred_lin_alg_backend, int* argc, char*** argv)
{
    utility::MPIUtils::init(argc, argv);
    fea_io::MPIUtils::init(argc, argv);

    auto backend =
      sparse_linear_algebra::LinAlgFactory::get().backend_from_str(preferred_lin_alg_backend);
    sparse_linear_algebra::LinAlgFactory::get().init(backend, argc, argv);
}
*/
}