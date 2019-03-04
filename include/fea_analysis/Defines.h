#pragma once

#include "flex_mesh/mesh/Mesh.h"
#include "flex_mesh/node/Node.h"
#include "flex_mesh/algorithms/Search.h"

namespace fea_analysis
{
	using Node = flex_mesh::Node;
	using Element = finite_element::FiniteElement;
	using Mesh = flex_mesh::Mesh<Node, Element>;
	using ElementSelection = flex_mesh::ElementSelection<Mesh, Element>;
	using NodeSelection = flex_mesh::NodeSelection<Mesh, Node>;
}