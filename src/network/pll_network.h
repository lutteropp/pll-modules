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

typedef struct pll_displayed_tree_s
{
  double * branch_lengths;
  pll_operation_t * operations;
  unsigned int ops_count;
  unsigned int * pmatrix_indices;
  unsigned int matrix_count;
} pll_displayed_tree_t;

PLL_EXPORT int pllmod_rnetwork_tree_buildarrays(pll_rnetwork_t * network, uint64_t tree_number, pll_displayed_tree_t * result);


typedef struct networkinfo_edge
{
  unsigned int left_index;
  unsigned int right_index;
  unsigned int pmatrix_index;
  double brlen;
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
  unsigned int tip_count;
  unsigned int partition_count;

  /* 0 = linked/shared, 1 = linked with scaler, 2 = unlinked */
  int brlen_linkage;
  double * linked_branch_lengths;

  pll_rnetwork_node_t * root;
  pll_rnetwork_t * network;

  unsigned int subnode_count;
  pll_rnetwork_node_t ** subnodes;

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
  pll_rnetwork_node_t ** travbuffer;
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


/* networkinfo */

PLL_EXPORT pllmod_networkinfo_t * pllmod_networkinfo_create(pll_rnetwork_node_t * root,
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
                                        pll_rnetwork_node_t * root);

PLL_EXPORT
int pllmod_networkinfo_get_branch_length_all(const pllmod_networkinfo_t * treeinfo,
                                          const pll_rnetwork_node_t * edge,
                                          double * lengths);

PLL_EXPORT int pllmod_networkinfo_set_branch_length(pllmod_networkinfo_t * networkinfo,
                                                 pll_rnetwork_node_t * edge,
                                                 double length);

PLL_EXPORT
int pllmod_networkinfo_set_branch_length_all(pllmod_networkinfo_t * networkinfo,
                                          pll_rnetwork_node_t * edge,
                                          const double * lengths);

PLL_EXPORT
int pllmod_networkinfo_set_branch_length_partition(pllmod_networkinfo_t * networkinfo,
                                                pll_rnetwork_node_t * edge,
                                                int partition_index,
                                                double length);

PLL_EXPORT
pll_rnetwork_t * pllmod_networkinfo_get_partition_tree(const pllmod_networkinfo_t * networkinfo,
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
                                             pll_rnetwork_node_t ** travbuffer,
                                             unsigned int travbuffer_size);

PLL_EXPORT void pllmod_networkinfo_invalidate_pmatrix(pllmod_networkinfo_t * networkinfo,
                                                   const pll_rnetwork_node_t * edge);

PLL_EXPORT void pllmod_networkinfo_invalidate_clv(pllmod_networkinfo_t * networkinfo,
                                               const pll_unode_t * edge);

PLL_EXPORT double pllmod_networkinfo_compute_loglh(pllmod_networkinfo_t * networkinfo,
                                                int incremental);

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
                                        pll_rnetwork_t * network);

PLL_EXPORT int pllmod_networkinfo_set_constraint_clvmap(pllmod_networkinfo_t * networkinfo,
                                                     const int * clv_index_map);

PLL_EXPORT int pllmod_networkinfo_set_constraint_network(pllmod_networkinfo_t * networkinfo,
                                                   const pll_rnetwork_t * cons_network);

PLL_EXPORT int pllmod_networkinfo_check_constraint(pllmod_networkinfo_t * networkinfo,
                                                pll_rnetwork_node_t * subtree,
                                                pll_rnetwork_node_t * regraft_edge);

#endif /* PLL_NETWORK_H_ */
