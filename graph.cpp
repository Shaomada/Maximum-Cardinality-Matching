#include "graph.hpp"

namespace ED
{

void Graph::add_edge (NodeId const node1_id, NodeId const node2_id)
{
   if (node1_id == node2_id) { throw "ED::Graph class does not support loops!"; }

   Node & node1 = _nodes[node1_id];
   node1.add_neighbor(node2_id);

   Node & node2 = _nodes[node2_id];
   node2.add_neighbor(node1_id);

   _num_edges += 1;
}

} // namespace ED
