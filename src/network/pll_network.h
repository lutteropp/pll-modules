/*
 Copyright (C) 2019 Sarah Lutteropp

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Contact: Diego Darriba <Diego.Darriba@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */
#ifndef PLL_NETWORK_H_
#define PLL_NETWORK_H_

#ifndef PLL_H_
#define PLL_H_
#include "pll.h"
#endif

#include "pll_tree.h"

/**
 * PLL Network utils module
 * Prefix: pll_network_, pll_unetwork_, pll_rnetwork_
 */
#define PLLMOD_NETWORK_DEFAULT_BRANCH_LENGTH 0.1

/* error codes (for this module, 3000-4000) ; B = 2^10+2^11*/
/* TBR errors (B + {2^2,2^1,2^0}) */
#define PLLMOD_NETWORK_ERROR_TBR_LEAF_BISECTION   3073 // B + {001}
#define PLLMOD_NETWORK_ERROR_TBR_OVERLAPPED_NODES 3074 // B + {010}
#define PLLMOD_NETWORK_ERROR_TBR_SAME_SUBNETWORK     3075 // B + {011}
#define PLLMOD_NETWORK_ERROR_TBR_MASK             3079 // B + {111}

/* NNI errors (B + {2^4,2^3}) */
#define PLLMOD_NETWORK_ERROR_NNI_INVALID_MOVE     3080 // B + {01...}
#define PLLMOD_NETWORK_ERROR_NNI_MASK             3096 // B + {11...}

/* SPR errors (B + {2^6,2^5}) */
#define PLLMOD_NETWORK_ERROR_SPR_INVALID_NODE     3104 // B + {01...}
#define PLLMOD_NETWORK_ERROR_SPR_MASK             3168 // B + {11...}

/* general errors (B + {2^8,2^7}) */
#define PLLMOD_NETWORK_ERROR_INTERCHANGE_LEAF     3200 // B + {01...}
#define PLLMOD_NETWORK_ERROR_INVALID_REARRAGE     3328 // B + {10...}
#define PLLMOD_NETWORK_ERROR_INVALID_NETWORK_SIZE    3456 // B + {10...}
#define PLLMOD_NETWORK_ERROR_INVALID_NETWORK         3584 // B + {10...}
#define PLLMOD_NETWORK_ERROR_INVALID_SPLIT        3712 // B + {10...}
#define PLLMOD_NETWORK_ERROR_EMPTY_SPLIT          3840 // B + {10...}
#define PLLMOD_NETWORK_ERROR_INVALID_THRESHOLD    3968 // B + {10...}
#define PLLMOD_NETWORK_ERROR_POLYPHYL_OUTGROUP    3970 // B + {10...}

#define PLLMOD_NETWORK_REARRANGE_SPR  0
#define PLLMOD_NETWORK_REARRANGE_NNI  1
#define PLLMOD_NETWORK_REARRANGE_TBR  2

#define PLLMOD_NETWORKINFO_PARTITION_ALL -1

typedef struct pll_displayed_tree_s
{
  double * branch_lengths;
  pll_operation_t * operations;
  unsigned int ops_count;
  unsigned int * pmatrix_indices;
  unsigned int matrix_count;
} pll_displayed_tree_t;

PLL_EXPORT int pllmod_rnetwork_tree_buildarrays(pll_rnetwork_t * network, uint64_t tree_number, pll_displayed_tree_t * result);

typedef struct pll_network_edge
{
  union
  {
    struct
    {
      pll_unetwork_node_t * parent;
      pll_unetwork_node_t * child;
    } unetwork;
    struct
    {
      pll_rnetwork_node_t * parent;
      pll_rnetwork_node_t * child;
    } rnetwork;
  } edge;
    double length;
} pll_network_edge_t;

typedef struct
{
  int rearrange_type;
  int rooted;
  double  likelihood;

  union {
    struct {
      void * prune_edge;
      void * regraft_edge;
      double prune_bl;        //! length of the pruned branch
      double prune_left_bl;   //! length of the removed branch when pruning
      double prune_right_bl;  //! length of the removed branch when pruning
      double regraft_bl;      //! length of the splitted branch when regrafting
    } SPR;
    struct {
      void * edge;
      double left_left_bl;
      double left_right_bl;
      double right_left_bl;
      double right_right_bl;
      double edge_bl;
      int type;
    } NNI;
    struct {
      void * bisect_edge;
      pll_network_edge_t reconn_edge;
      double bisect_left_bl;
      double bisect_right_bl;
      double reconn_parent_left_bl;
      double reconn_parent_right_bl;
      double reconn_child_left_bl;
      double reconn_child_right_bl;
    } TBR;
  };
} pll_network_rollback_t;

typedef struct networkinfo_edge
{
  unsigned int left_index;
  unsigned int right_index;
  unsigned int pmatrix_index;
  double brlen;
  double prob;
} pllmod_networkinfo_edge_t;

typedef struct networkinfo_topology
{
  unsigned int edge_count;
  unsigned int brlen_set_count;
  unsigned int root_index;
  pllmod_networkinfo_edge_t * edges;
  double ** branch_lengths;
} pllmod_networkinfo_topology_t;

typedef struct networkinfo
{
  // dimensions
  unsigned int partition_count;

  /* 0 = linked/shared, 1 = linked with scaler, 2 = unlinked */
  int brlen_linkage;
  double * linked_branch_lengths;

  pll_unetwork_node_t * root;
  pll_unetwork_t * network;

  unsigned int subnode_count;
  pll_unetwork_node_t ** subnodes;

  // partitions & partition-specific stuff
  pll_partition_t ** partitions;
  double * alphas;
  int * gamma_mode; /* discrete GAMMA rates computation mode (mean, median) */
  unsigned int ** param_indices;
  int ** subst_matrix_symmetries;
  double ** branch_lengths;
  double * brlen_scalers;
  double * partition_loglh;
  int * params_to_optimize;

  // partition that have been initialized (useful for parallelization)
  unsigned int init_partition_count;
  unsigned int * init_partition_idx;
  pll_partition_t ** init_partitions;

  /* tree topology constraint */
  unsigned int * constraint;

  /* precomputation buffers for derivatives (aka "sumtable") */
  double ** deriv_precomp;

  // invalidation flags
  char ** clv_valid;
  char ** pmatrix_valid;

  // buffers
  pll_unetwork_node_t ** travbuffer;
  unsigned int * matrix_indices;
  pll_operation_t * operations;

  // partition on which all operations should be performed
  int active_partition;

  // general-purpose counter
  unsigned int counter;

  // parallelization stuff
  void * parallel_context;
  void (*parallel_reduce_cb)(void *, double *, size_t, int);
} pllmod_networkinfo_t;

typedef struct
{
  unsigned int node_count;
  pll_unetwork_node_t ** nodes;

  unsigned int partition_count;
  unsigned int * partition_indices;

  pll_unetwork_t * network;
  double ** probs;
} pllmod_ancestral_network_t;

/* Additional utilities */
/* functions at pll_unetwork.c */

PLL_EXPORT int pllmod_unetwork_set_clv_minimal(pll_unetwork_node_t * root,
                                         unsigned int tip_count);

PLL_EXPORT int pllmod_unetwork_traverse_apply(pll_unetwork_node_t * root,
                                        int (*cb_pre_trav)(pll_unetwork_node_t *,
                                                           void *),
                                        int (*cb_in_trav)(pll_unetwork_node_t *,
                                                          void *),
                                        int (*cb_post_trav)(pll_unetwork_node_t *,
                                                            void *),
                                        void *data);

PLL_EXPORT int pllmod_unetwork_is_tip(const pll_unetwork_node_t * node);

PLL_EXPORT void pllmod_unetwork_set_length(pll_unetwork_node_t * edge,
                                     double length);

PLL_EXPORT void pllmod_unetwork_set_length_recursive(pll_unetwork_t * network,
                                                  double length,
                                                  int missing_only);


PLL_EXPORT void pllmod_unetwork_scale_branches(pll_unetwork_t * network,
                                            double branch_length_scaler);

PLL_EXPORT void pllmod_unetwork_scale_branches_all(pll_unetwork_node_t * root,
                                                double branch_length_scaler);

PLL_EXPORT void pllmod_unetwork_scale_subnetwork_branches(pll_unetwork_node_t * root,
                                                    double branch_length_scaler);

PLL_EXPORT pll_unetwork_t * pllmod_unetwork_resolve_multi(const pll_unetwork_t * multi_network,
                                                    unsigned int random_seed,
                                                    int * clv_index_map);

PLL_EXPORT int pllmod_unetwork_root_inplace(pll_unetwork_t * network);

PLL_EXPORT int pllmod_unetwork_outgroup_root(pll_unetwork_t * network,
                                          unsigned int * outgroup_tip_ids,
                                          unsigned int outgroup_size,
                                          int add_root_node);


PLL_EXPORT double pllmod_unetwork_compute_lk(pll_partition_t * partition,
                                       pll_unetwork_node_t * network,
                                       const unsigned int * params_indices,
                                       int update_pmatrices,
                                       int update_partials);

PLL_EXPORT int pllmod_rnetwork_traverse_apply(pll_rnode_t * root,
                                           int (*cb_pre_trav)(pll_rnode_t *,
                                                              void *),
                                           int (*cb_in_trav)(pll_rnode_t *,
                                                             void *),
                                           int (*cb_post_trav)(pll_rnode_t *,
                                                               void *),
                                           void *data);

/* networkinfo */

PLL_EXPORT pllmod_networkinfo_t * pllmod_networkinfo_create(pll_unetwork_t * network,
                                                     unsigned int tips,
                                                     unsigned int partitions,
                                                     int brlen_linkage);

PLL_EXPORT
int pllmod_networkinfo_set_parallel_context(pllmod_networkinfo_t * networkinfo,
                                         void * parallel_context,
                                         void (*parallel_reduce_cb)(void *,
                                                                    double *,
                                                                    size_t,
                                                                    int op));

PLL_EXPORT int pllmod_networkinfo_init_partition(pllmod_networkinfo_t * networkinfo,
                                           unsigned int partition_index,
                                           pll_partition_t * partition,
                                           int params_to_optimize,
                                           int gamma_mode,
                                           double alpha,
                                           const unsigned int * param_indices,
                                           const int * subst_matrix_symmetries);

PLL_EXPORT int pllmod_networkinfo_set_active_partition(pllmod_networkinfo_t * networkinfo,
                                                    int partition_index);

PLL_EXPORT int pllmod_networkinfo_set_root(pllmod_networkinfo_t * treeinfo,
                                        pll_unetwork_node_t * root);

PLL_EXPORT
int pllmod_networkinfo_get_branch_length_all(const pllmod_networkinfo_t * treeinfo,
                                          const pll_unetwork_node_t * edge,
                                          double * lengths);

PLL_EXPORT int pllmod_networkinfo_set_branch_length(pllmod_networkinfo_t * networkinfo,
                                                 pll_unetwork_node_t * edge,
                                                 double length);

PLL_EXPORT
int pllmod_networkinfo_set_branch_length_all(pllmod_networkinfo_t * networkinfo,
                                          pll_unetwork_node_t * edge,
                                          const double * lengths);

PLL_EXPORT
int pllmod_networkinfo_set_branch_length_partition(pllmod_networkinfo_t * networkinfo,
		                                        pll_unetwork_node_t * edge,
                                                int partition_index,
                                                double length);

PLL_EXPORT
pll_unetwork_t * pllmod_networkinfo_get_partition_network(const pllmod_networkinfo_t * networkinfo,
                                                 int partition_index);

PLL_EXPORT
pllmod_networkinfo_topology_t * pllmod_networkinfo_get_topology(const pllmod_networkinfo_t * networkinfo,
                                                          pllmod_networkinfo_topology_t * topol);

PLL_EXPORT
int pllmod_networkinfo_set_topology(pllmod_networkinfo_t * networkinfo,
                                 const pllmod_networkinfo_topology_t * topol);

PLL_EXPORT
int pllmod_networkinfo_destroy_topology(pllmod_networkinfo_topology_t * topol);

PLL_EXPORT int pllmod_betworkinfo_destroy_partition(pllmod_networkinfo_t * networkinfo,
                                                 unsigned int partition_index);

PLL_EXPORT void pllmod_networkinfo_destroy(pllmod_networkinfo_t * networkinfo);

PLL_EXPORT int pllmod_networkinfo_update_prob_matrices(pllmod_networkinfo_t * networkinfo,
                                                    int update_all);

PLL_EXPORT void pllmod_networkinfo_invalidate_all(pllmod_networkinfo_t * networkinfo);

PLL_EXPORT int pllmod_networkinfo_validate_clvs(pllmod_networkinfo_t * networkinfo,
		                                     pll_unetwork_node_t ** travbuffer,
                                             unsigned int travbuffer_size);

PLL_EXPORT void pllmod_networkinfo_invalidate_pmatrix(pllmod_networkinfo_t * networkinfo,
                                                   const pll_unetwork_node_t * edge);

PLL_EXPORT void pllmod_networkinfo_invalidate_clv(pllmod_networkinfo_t * networkinfo,
                                               const pll_unetwork_node_t * edge);

PLL_EXPORT double pllmod_networkinfo_compute_loglh(pllmod_networkinfo_t * networkinfo,
                                                int incremental, int update_pmatrices);

PLL_EXPORT double pllmod_networkinfo_compute_loglh_flex(pllmod_networkinfo_t * networkinfo,
                                                     int incremental,
                                                     int update_pmatrices);

PLL_EXPORT
int pllmod_networkinfo_scale_branches_all(pllmod_networkinfo_t * networkinfo, double scaler);

PLL_EXPORT
int pllmod_networkinfo_scale_branches_partition(pllmod_networkinfo_t * networkinfo,
                                             unsigned int partition_idx,
                                             double scaler);

PLL_EXPORT
int pllmod_networkinfo_normalize_brlen_scalers(pllmod_networkinfo_t * networkinfo);

PLL_EXPORT int pllmod_networkinfo_set_network(pllmod_networkinfo_t * networkinfo,
                                        pll_unetwork_t * network);

PLL_EXPORT int pllmod_networkinfo_set_constraint_clvmap(pllmod_networkinfo_t * networkinfo,
                                                     const int * clv_index_map);

PLL_EXPORT int pllmod_networkinfo_set_constraint_network(pllmod_networkinfo_t * networkinfo,
                                                   const pll_unetwork_t * cons_network);

PLL_EXPORT int pllmod_networkinfo_check_constraint(pllmod_networkinfo_t * networkinfo,
		                                          pll_unetwork_node_t * subnetwork,
		                                          pll_unetwork_node_t * regraft_edge);

/* Network construction */
/* functions at pll_rnetwork.c */

PLL_EXPORT pll_rnetwork_t * pllmod_rnetwork_create_random(unsigned int taxa_count,
                                                    const char * const* names,
                                                    unsigned int random_seed);

PLL_EXPORT int pllmod_rnetwork_extend_random(pll_rnetwork_t * network,
                                          unsigned int ext_taxa_count,
                                          const char * const* ext_names,
                                          unsigned int random_seed);

PLL_EXPORT
pll_rnetwork_t * pllmod_rnetwork_create_parsimony(unsigned int taxon_count,
                                            unsigned int seq_length,
                                            char * const * names,
                                            char * const * sequences,
                                            const unsigned int * site_weights,
                                            const pll_state_t * map,
                                            unsigned int states,
                                            unsigned int attributes,
                                            unsigned int random_seed,
                                            unsigned int * score);

pll_rnetwork_t * pllmod_rnetwork_create_parsimony_multipart(unsigned int taxon_count,
                                                      char * const * taxon_names,
                                                      unsigned int partition_count,
                                                      pll_partition_t * const * partitions,
                                                      unsigned int random_seed,
                                                      unsigned int * score);

PLL_EXPORT pll_rnetwork_t * pllmod_rnetwork_resolve_multi(const pll_rnetwork_t * multi_network,
                                                    unsigned int random_seed,
                                                    int * clv_index_map);

PLL_EXPORT int pllmod_rnetwork_is_tip(const pll_rnetwork_node_t * node);

PLL_EXPORT void pllmod_rnetwork_set_length(pll_rnetwork_node_t * edge,
                                     double length);

PLL_EXPORT void pllmod_rnetwork_set_length_recursive(pll_rnetwork_t * network,
                                                  double length,
                                                  int missing_only);

PLL_EXPORT int pllmod_rnetwork_outgroup_root(pll_rnetwork_t * network,
                                          unsigned int * outgroup_tip_ids,
                                          unsigned int outgroup_size,
                                          int add_root_node);

/* Discrete operations */
/* functions at unetwork_distances.c */

PLL_EXPORT unsigned int pllmod_unetwork_rf_distance(pll_unetwork_node_t * t1,
                                                 pll_unetwork_node_t * t2,
                                                 unsigned int tip_count,
												 unsigned int reticulation_count);

/* check that node ids and tip labels agree in both networks */
PLL_EXPORT int pllmod_unetwork_consistency_check(pll_unetwork_t * t1,
                                              pll_unetwork_t * t2);

/* if 2 different networks are parsed from newick node ids migh have been set
   in a different order, so this function sets node ids in t2 such that
   node ids and tip labels agree in both networks */
PLL_EXPORT int pllmod_unetwork_consistency_set(pll_unetwork_t * t1,
                                            pll_unetwork_t * t2);

PLL_EXPORT unsigned int pllmod_unetwork_split_rf_distance(pll_split_t * s1,
                                                       pll_split_t * s2,
                                                       unsigned int tip_count);

PLL_EXPORT pll_split_t * pllmod_unetwork_split_create(const pll_unetwork_node_t * network,
                                                   unsigned int tip_count,
												   unsigned int reticulation_count,
                                                   pll_unetwork_node_t ** split_to_node_map);

PLL_EXPORT pll_split_t pllmod_unetwork_split_from_tips(unsigned int * subnetwork_tip_ids,
                                                    unsigned int subnetwork_size,
                                                    unsigned int tip_count);

PLL_EXPORT void pllmod_unetwork_split_normalize_and_sort(pll_split_t * s,
                                                      unsigned int tip_count,
                                                      unsigned int n_splits,
                                                      int keep_first);

PLL_EXPORT void pllmod_unetwork_split_show(pll_split_t split,
                                        unsigned int tip_count);

PLL_EXPORT void pllmod_unetwork_split_destroy(pll_split_t * split_list);

PLL_EXPORT unsigned int pllmod_unetwork_split_lightside(pll_split_t split,
                                                     unsigned int tip_count);

PLL_EXPORT unsigned int pllmod_unetwork_split_hamming_distance(pll_split_t s1,
                                                            pll_split_t s2,
                                                            unsigned int tip_count);

PLL_EXPORT pll_split_t * pll_unetwork_split_newick_string(char * s,
                                                       unsigned int tip_count,
                                                       string_hashtable_t * names_hash);

PLL_EXPORT
bitv_hashtable_t * pllmod_unetwork_split_hashtable_create(unsigned int tip_count,
                                                       unsigned int slot_count);

PLL_EXPORT bitv_hash_entry_t *
pllmod_unetwork_split_hashtable_insert_single(bitv_hashtable_t * splits_hash,
                                           pll_split_t split,
                                           double support);

PLL_EXPORT bitv_hashtable_t *
pllmod_unetwork_split_hashtable_insert(bitv_hashtable_t * splits_hash,
                                    pll_split_t * splits,
                                    unsigned int tip_count,
                                    unsigned int split_count,
                                    const double * support,
                                    int update_only);

PLL_EXPORT bitv_hash_entry_t *
pllmod_unetwork_split_hashtable_lookup(bitv_hashtable_t * splits_hash,
                                    pll_split_t split,
                                    unsigned int tip_count);

PLL_EXPORT
void pllmod_unetwork_split_hashtable_destroy(bitv_hashtable_t * hash);


/* Topological rearrangements */
/* functions at pll_unetwork.c */

PLL_EXPORT int pllmod_unetwork_tbr(pll_unetwork_node_t * b_edge,
                                pll_network_edge_t * r_edge,
                                pll_network_rollback_t * rollback_info);

PLL_EXPORT int pllmod_unetwork_spr(pll_unetwork_node_t * p_edge,
                                pll_unetwork_node_t * r_edge,
                                pll_network_rollback_t * rollback_info);

/* type = {PLL_NNI_NEXT, PLL_NNI_NEXTNEXT} */
PLL_EXPORT int pllmod_unetwork_nni(pll_unetwork_node_t * edge,
                                int type,
                                pll_network_rollback_t * rollback_info);

PLL_EXPORT int pllmod_network_rollback(pll_network_rollback_t * rollback_info);

PLL_EXPORT pll_unetwork_node_t * pllmod_unetwork_serialize(pll_unetwork_node_t * network,
                                                unsigned int tip_count, unsigned int reticulation_count);

PLL_EXPORT pll_unetwork_t * pllmod_unetwork_expand(pll_unetwork_node_t * serialized_network,
                                             unsigned int tip_count, unsigned int reticulation_count);

PLL_EXPORT int pllmod_unetwork_draw_support(pll_unetwork_t * ref_network,
                                         const double * support,
                                         pll_unetwork_node_t ** node_map,
                                         char * (*cb_serialize)(double));

PLL_EXPORT int pllmod_unetwork_draw_support(pll_unetwork_t * ref_network,
                                         const double * support,
                                         pll_unetwork_node_t ** node_map,
                                         char * (*cb_serialize)(double));

/* functions at rnetwork_operations.c */

PLL_EXPORT int pllmod_rnetwork_spr(pll_rnetwork_node_t * p_node,
		pll_rnetwork_node_t * r_tree,
		pll_rnetwork_node_t ** root,
                                pll_network_rollback_t * rollback_info);

PLL_EXPORT int pllmod_rnetwork_get_sibling_pointers(pll_rnetwork_node_t * node,
		pll_rnetwork_node_t ***self,
		pll_rnetwork_node_t ***sister);

PLL_EXPORT pll_rnetwork_node_t * pllmod_rnetwork_prune(pll_rnetwork_node_t * node);

PLL_EXPORT int pllmod_rnetwork_regraft(pll_rnetwork_node_t * node,
                                    pll_rnetwork_node_t * network);

PLL_EXPORT int pllmod_rnetwork_connect_nodes(pll_rnetwork_node_t * parent,
                                          pll_rnetwork_node_t * child,
                                           double length, double prob);

/* Network construction */
/* functions at pll_unetwork.c */

PLL_EXPORT pll_unetwork_t * pllmod_unetwork_create_random(unsigned int taxa_count,
                                                    const char * const* names,
                                                    unsigned int random_seed);

PLL_EXPORT int pllmod_unetwork_extend_random(pll_unetwork_t * network,
                                          unsigned int ext_taxa_count,
                                          const char * const* ext_names,
                                          unsigned int random_seed);

PLL_EXPORT
pll_unetwork_t * pllmod_unetwork_create_parsimony(unsigned int taxon_count,
                                            unsigned int seq_length,
                                            char * const * names,
                                            char * const * sequences,
                                            const unsigned int * site_weights,
                                            const pll_state_t * map,
                                            unsigned int states,
                                            unsigned int attributes,
                                            unsigned int random_seed,
                                            unsigned int * score);

pll_unetwork_t * pllmod_unetwork_create_parsimony_multipart(unsigned int taxon_count,
                                                      char * const * taxon_names,
                                                      unsigned int partition_count,
                                                      pll_partition_t * const * partitions,
                                                      unsigned int random_seed,
                                                      unsigned int * score);

PLL_EXPORT pll_unetwork_t * pllmod_unetwork_resolve_multi(const pll_unetwork_t * multi_network,
                                                    unsigned int random_seed,
                                                    int * clv_index_map);

PLL_EXPORT int pllmod_unetwork_is_tip(const pll_unetwork_node_t * node);

PLL_EXPORT void pllmod_unetwork_set_length(pll_unetwork_node_t * edge,
                                     double length);

PLL_EXPORT void pllmod_unetwork_set_length_recursive(pll_unetwork_t * network,
                                                  double length,
                                                  int missing_only);

PLL_EXPORT int pllmod_unetwork_outgroup_root(pll_unetwork_t * network,
                                          unsigned int * outgroup_tip_ids,
                                          unsigned int outgroup_size,
                                          int add_root_node);

/* Additional utilities in pll_unetwork.c */
PLL_EXPORT int pllmod_unetwork_traverse_apply(pll_unetwork_node_t * root,
                                        int (*cb_pre_trav)(pll_unetwork_node_t *,
                                                           void *),
                                        int (*cb_in_trav)(pll_unetwork_node_t *,
                                                          void *),
                                        int (*cb_post_trav)(pll_unetwork_node_t *,
                                                            void *),
                                        void *data);


/* functions at unetwork_operations.c */

PLL_EXPORT int pllmod_unetwork_bisect(pll_unetwork_node_t * edge,
                                   pll_unetwork_node_t ** parent_subnetwork,
                                   pll_unetwork_node_t ** child_subnetwork);

PLL_EXPORT pll_network_edge_t pllmod_unetwork_reconnect(pll_network_edge_t * edge,
                                                  pll_unetwork_node_t * pruned_edge);

PLL_EXPORT pll_unetwork_node_t * pllmod_unetwork_prune(pll_unetwork_node_t * edge);

PLL_EXPORT int pllmod_unetwork_regraft(pll_unetwork_node_t * edge,
                                    pll_unetwork_node_t * network);

PLL_EXPORT int pllmod_unetwork_interchange(pll_unetwork_node_t * edge1,
                                        pll_unetwork_node_t * edge2);

PLL_EXPORT pll_unetwork_node_t * pllmod_unetwork_create_node(unsigned int clv_index,
                                                  int scaler_index,
                                                  char * label,
                                                  void * data);

PLL_EXPORT int pllmod_unetwork_connect_nodes(pll_unetwork_node_t * parent,
                                          pll_unetwork_node_t * child,
                                           double length, double prob);

PLL_EXPORT int pllmod_unetwork_nodes_at_node_dist(pll_unetwork_node_t * node,
                                            pll_unetwork_node_t ** outbuffer,
                                            unsigned int * node_count,
                                            unsigned int min_distance,
                                            unsigned int max_distance);

PLL_EXPORT int pllmod_unetwork_nodes_at_edge_dist(pll_unetwork_node_t * edge,
                                               pll_unetwork_node_t ** outbuffer,
                                               unsigned int * node_count,
                                               unsigned int min_distance,
                                               unsigned int max_distance);

#endif /* PLL_NETWORK_H_ */
