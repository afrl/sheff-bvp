#pragma once

#include "flex_mesh/mesh/MeshPartition.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/algorithms/Selection.h"
#include "finite_element/FiniteElement.h"

namespace fea_analysis
{
	using Node = flex_mesh::Node;
	using Element = finite_element::FiniteElement;
	using Mesh = flex_mesh::MeshPartition<Node, Element>;
	using ElementSelection = flex_mesh::ElementSelection<Mesh, Element>;
	using NodeSelection = flex_mesh::NodeSelection<Mesh, Node>;

	/// Must be called ONLY once and BEFORE an analysis is reaches the setup phase.
	/// \param preferred_lin_alg_backend - a string referring to the type of backend
	///        to use for linear algebra. (Ex: MKL, PETSC, PARLAB)
	/// \param argc - count of commandline arguments to pass on to other modules.
	/// \param argc - commandline arguments to pass on to other modules.
	void init(const std::string preferred_lin_alg_backend, int* argc, char*** argv);
}