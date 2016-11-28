#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include "graph.hpp"
#include "edmonds.hpp"

template <typename Fn1, typename Fn2>
void parse(char *filename, Fn1 && fn1, Fn2 && fn2) {
   std::ifstream file(filename);
   std::string line;
   while (std::getline(file, line)) {
      std::stringstream s(line);
      std::string word;
      s >> word;
      if (!strcmp(word.c_str(), "p")) {
         ED::NodeId n;
         std::size_t m;
         s >> word >> n >> m;
         fn1(n, m);
      } else if (!strcmp(word.c_str(), "e")){
         ED::NodeId v;
         ED::NodeId w;
         s >> v >> w;
         fn2(v, w);
      }
   }
}

void print(ED::Graph const &G, std::vector<ED::NodeId> const &matching)
{
   ED::NodeId nr_covered_nodes = 0;
   for (ED::NodeId id = 0; id < G.num_nodes(); id++) {
      if (matching[id] != ED::invalid_node_id) {
         nr_covered_nodes++;
      }
   }
   std::cout << "p edge " << G.num_nodes() << " " << nr_covered_nodes/2 << std::endl;
   for (ED::NodeId id = 0; id < G.num_nodes(); id++) {
      if (matching[id] < id) {
         std::cout
            << "e "
            << ED::ED_id_to_dimacs_id(matching[id]) << " "
            << ED::ED_id_to_dimacs_id(id)
            << std::endl;
      }
      // std::cout << id << ": " << matching[id] << std::endl;
   }
}

int main (int argc, char **argv)
{
   if ( (argc == 3 && !strcmp(argv[1], "--graph"))
      ||(argc == 5 && !strcmp(argv[1], "--graph") && !strcmp(argv[3], "--hint"))) {
      ED::Graph *G;
      parse(
         argv[2],
         [&](ED::NodeId num_nodes, std::size_t) -> void {
            G = new ED::Graph(num_nodes);
         },
         [&](ED::NodeId v, ED::NodeId w) -> void {
            G->add_edge(ED::dimacs_id_to_ED_id(v), ED::dimacs_id_to_ED_id(w));
         }
      );
      std::vector<ED::NodeId> matching(G->num_nodes(), ED::invalid_node_id);
      if (argc == 5) {
         parse(
            argv[4],
            [&](ED::NodeId, std::size_t) -> void {},
            [&](ED::NodeId v, ED::NodeId w) -> void {
               matching[v] = w;
               matching[w] = v;
            }
         );
      }
      ED::edmonds(*G, matching);
      print(*G, matching);
      delete G;
   } else {
      std::cout
         << "Usage: "
         << argv[0]
         << " --graph file1.dmx"
         << " [--hint file2.dmx]"
         << std::endl;
   }
}
