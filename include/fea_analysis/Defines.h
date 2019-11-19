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
}