#include "edmonds.hpp"

namespace ED // for Edmonds
{

typedef NodeId dist_t;
dist_t const invalid_dist = std::numeric_limits<dist_t>::max();
typedef NodeId lbl_t;
lbl_t const invalid_lbl = std::numeric_limits<lbl_t>::max();
struct LabelData {
   LabelData(NodeId id) : root(id), labeled_vertices(1, id) {};
   NodeId root;
   std::vector<NodeId> labeled_vertices;
};

void edmonds (Graph const &G, std::vector<NodeId> &matching)
{
   /* The following will be required for all try_augment calls.
    * For runtime reasons, we only initialise them once, then
    * use try_augment::clean to prepare them for the next call.
    */
   std::vector<int> deleted(G.num_nodes(), 0);
   std::vector<NodeId> prev(G.num_nodes(), invalid_node_id);
   std::vector<lbl_t> label(G.num_nodes(), invalid_lbl);
   std::vector<dist_t> d(G.num_nodes(), invalid_dist);
   std::vector<std::size_t> next_edge_idx(G.num_nodes(), 0);

   auto try_augment = [&](NodeId id) -> void {
      std::vector<NodeId> even_vertices(1, id);
      size_t next_vertex_idx = 0;

      std::vector<LabelData> label_data (1, LabelData(id));
      label[id] = 0;

      auto augment = [&](NodeId x, NodeId y) -> void {
         while (matching[x] != invalid_node_id) {
            NodeId z = matching[x];
            matching[y] = x;
            matching[x] = y;
            y = z;
            x = prev[z];
         }
         matching[y] = x;
         matching[x] = y;
      };

      auto grow = [&](NodeId x, NodeId y) -> void {
         prev[y] = x;
         label[y] = label_data.size();
         label_data.emplace_back(y);

         NodeId z = matching[y];
         even_vertices.push_back(z);
         label[z] = label_data.size();
         label_data.emplace_back(z);
      };

      auto contract = [&](NodeId x, NodeId y) -> void {
         std::vector<lbl_t> labels_found;
         auto merge_labels = [&](lbl_t lbl_root) -> void {
            auto size = [&](lbl_t lbl) -> std::size_t {
               return label_data[lbl].labeled_vertices.size();
            };

            lbl_t new_lbl = lbl_root;
            for (lbl_t &lbl : labels_found) {
               if (size(lbl) > size(new_lbl)) {
                  std::swap(new_lbl, lbl);
               }
            }

            label_data[new_lbl].root = label_data[lbl_root].root;
            for (size_t lbl : labels_found) {
               for (NodeId id : label_data[lbl].labeled_vertices) {
                  label[id] = new_lbl;
                  label_data[new_lbl].labeled_vertices.push_back(id);
               }
               label_data[lbl].labeled_vertices.clear();
            }
         };

         auto root_pseudonode = [&](NodeId id) -> NodeId {
            return label_data[label[id]].root;
         };

         if (prev[x] == invalid_node_id) prev[x] = y;
         if (prev[y] == invalid_node_id) prev[y] = x;

         while (label[x] != label[y]) {
            if (d[root_pseudonode(x)] > d[root_pseudonode(y)]) std::swap(x, y);
            NodeId z = matching[root_pseudonode(y)];
            y = prev[z];
            if (prev[y] == invalid_node_id) prev[y] = z;
            even_vertices.push_back(z);
            labels_found.push_back(label[y]);
            labels_found.push_back(label[z]);
         }
         merge_labels(label[x]);
      };

      auto clean = [&](bool augmented) -> void {
         for (LabelData const & data : label_data) {
            for (NodeId id : data.labeled_vertices) {
               if (augmented) {
                  prev[id] = invalid_node_id;
                  label[id] = invalid_lbl;
                  d[id] = invalid_dist;
                  next_edge_idx[id] = 0;
               } else {
                  deleted[id] = 1;
               }
            }
         }
      };

      while (next_vertex_idx < even_vertices.size()) {
         NodeId x = even_vertices[next_vertex_idx];
         if (next_edge_idx[x] == G.node(x).neighbors().size()) {
            next_vertex_idx++;
            continue;
         }
         NodeId y = G.node(x).neighbors()[next_edge_idx[x]++];
         if (deleted[y] || label[x] == label[y]) continue;

         if (label[y] == invalid_lbl) {
            if (matching[y] == invalid_node_id) {
               augment(x, y);
               clean(true);
               return;
            } else {
               grow(x, y);
            }
         } else if (d[y]%2 == 0) {
            contract(x, y);
         }
      }
      clean(false);
   };

   for (NodeId id = 0; id < G.num_nodes(); id++)
   {
      if (matching[id] == invalid_node_id) {
         try_augment(id);
      }
   }
}

}
