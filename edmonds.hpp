#ifndef EDMONDS_HPP
#define EDMONDS_HPP

#include "graph.hpp"

namespace ED // for Edmonds
{

void edmonds (Graph const &G, std::vector<NodeId> &matching);

} // namespace ED

#endif
