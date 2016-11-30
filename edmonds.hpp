#ifndef EDMONDS_HPP
#define EDMONDS_HPP

#include "graph.hpp"

namespace ED // for Edmonds
{

/* NOTE Matchings are saved such that:
 *    for {v, w} matching edge, matching[v]=w
 *    for v uncovered, matching[v]=invalid_node_id
 */

/// For a matching in G, augments it in place until it's maximum.
void edmonds (Graph const &G, std::vector<NodeId> &matching);

} // namespace ED

#endif
