#include "edmonds.hpp"
#include <tuple>

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
   /// Instead of explicitly deleteing a vertex v, we set deleted[v]=1
   std::vector<int> deleted(G.num_nodes(), 0);

   /* NOTE we initialise the following vectors that will be used in each
    * call of try_augment here for runtime reasons:
    * prev, rep, d: We will only require the size of the vectors to be
    * G.num_nodes(), all entries used will be initialised during try_augment.
    * label, next_edge_idx: We require these vectors to be initialised with
    * invalid_node_id and 0 respectively at all not deleted vertices whenever
    * we call try_augment. To ensure this, at the end of each try_augment
    * call, we will call the lambda clean.
    */

   /** Let {v, w} be a non-matching edge that we can use to backtrack
    * matching-alternating within the current tree from w to the root
    * of the tree. Then we set prev[r]=v and rep[r]=w, where w is the
    * root of the pseudonode w currently belongs to
    */
   /// @{
   std::vector<NodeId> prev(G.num_nodes(), invalid_node_id);
   std::vector<NodeId> rep(G.num_nodes(), invalid_node_id);
   /// @}
   /** This does not cash the actual distance to the root. We only
    * will require d[v] even iff v belongs to an even pseudonode and
    * d[v] < d[w] if v and w roots of pseudonodes V, W with V between
    * the root of the tree and W.
    */
   std::vector<dist_t> d(G.num_nodes(), invalid_dist);
   /** We assign each vertex in the tree a label in order to check if
    * two vertices belong to the same pseudonode. This label is also
    * used as an index for the vector label_data in try_augment.
    */
   std::vector<lbl_t> label(G.num_nodes(), invalid_lbl);
   /// Used to iterate over the incident edges of vertices
   std::vector<std::size_t> next_edge_idx(G.num_nodes(), 0);

   auto try_augment = [&](NodeId id) -> void {
      // initialise the tree
      label[id] = 0;
      d[id] = 0;

      /// used with next_edge_idx to find an even vertex and an incident
      /// unprocessed edge or decide non exist in constant time
      /// @{
      std::vector<NodeId> even_vertices(1, id);
      size_t next_vertex_idx = 0;
      /// @}

      /// used to find the root of a label or all vertices with a label
      std::vector<LabelData> label_data (1, LabelData(id));

      auto augment = [&](NodeId x, NodeId y) -> void {
         /// keeps track of which edges we need to add to our matching
         /// @{
         std::vector< std::pair< NodeId, NodeId> > add(1, {x, y});
         std::size_t idx = 0;
         /// @}

         /// if v is covered, removes the covering edge {v, w} from matching,
         /// finds an edge we can use to re-cover w and saves it in add
         auto uncover = [&](NodeId v) -> void {
            if (matching[v] != invalid_node_id) {
               NodeId w = matching[v];
               matching[v] = invalid_node_id;
               matching[w] = invalid_node_id;
               add.emplace_back(prev[w], rep[w]);
            }
         };

         while (idx < add.size()) {
            std::tie(x, y) = add[idx++];
            uncover(x);
            uncover(y);
            matching[x] = y;
            matching[y] = x;
         }
      };

      auto grow = [&](NodeId x, NodeId y) -> void {
         // data stored for the edge {x, y}
         prev[y] = x;
         rep[y] = y;
         // data stored for the odd vertex x
         label[y] = label_data.size();
         d[y] = d[x]+1;
         label_data.emplace_back(y);
         /// data stored for the even vertex matching[x]
         NodeId z = matching[y];
         even_vertices.push_back(z);
         label[z] = label_data.size();
         d[z] = d[y]+1;
         label_data.emplace_back(z);
      };

      auto contract = [&](NodeId x, NodeId y) -> void {
         /// keeps track of labels during backtracking to merge them later
         std::vector<lbl_t> labels_found;

         /// merges all labels found while backtracking and the label lbl_root
         auto merge_labels = [&](lbl_t lbl_root) -> void {
            auto size = [&](lbl_t lbl) -> std::size_t {
               return label_data[lbl].labeled_vertices.size();
            };

            // choose new_lbl as the minimum size lbl
            lbl_t new_lbl = lbl_root;
            for (lbl_t &lbl : labels_found) {
               if (size(lbl) > size(new_lbl)) {
                  std::swap(new_lbl, lbl);
               }
            }

            // relabel vertices and updates label_data[root]
            label_data[new_lbl].root = label_data[lbl_root].root;
            for (size_t lbl : labels_found) {
               for (NodeId id : label_data[lbl].labeled_vertices) {
                  label[id] = new_lbl;
                  label_data[new_lbl].labeled_vertices.push_back(id);
               }
               // clearing needed only for runtime of clean
               label_data[lbl].labeled_vertices.clear();
            }
         };

         NodeId pred_x = y, root_x = label_data[label[x]].root;
         NodeId pred_y = x, root_y = label_data[label[y]].root;
         while (label[x] != label[y]) {
            if (d[root_x] < d[root_y]) {
               std::swap (x, y);
               std::swap (pred_x, pred_y);
               std::swap (root_x, root_y);
            }
            // The actual backtracking would be x=prev[matching[label_data[x].root]].
            // However, during backtracking we want to also safe the incoming edge of x,
            // safe both labels and update the previously odd vertex to even.
            prev[root_x] = pred_x;
            rep[root_x] = x;
            pred_x = matching[root_x];
            labels_found.push_back(label[x]);
            labels_found.push_back(label[pred_x]);
            d[pred_x] = d[x];
            even_vertices.push_back(pred_x);
            x = prev[pred_x];
            root_x = label_data[label[x]].root;
         }
         merge_labels(label[x]);
      };

      auto clean = [&](bool augmented) -> void {
         for (LabelData const & data : label_data) {
            for (NodeId id : data.labeled_vertices) {
               if (augmented) {
                  label[id] = invalid_lbl;
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
