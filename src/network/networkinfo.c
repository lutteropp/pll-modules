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

    Contact: Alexey Kozlov <Alexey.Kozlov@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "pll_network.h"

#include "../pllmod_common.h"


PLL_EXPORT pllmod_networkinfo_t * pllmod_networkinfo_create(pll_rnetwork_node_t * root,
                                                     unsigned int tips,
                                                     unsigned int partitions,
                                                     int brlen_linkage)
{
  return PLL_FAILURE;
}

PLL_EXPORT
int pllmod_networkinfo_set_parallel_context(pllmod_networkinfo_t * networkinfo,
                                         void * parallel_context,
                                         void (*parallel_reduce_cb)(void *,
                                                                    double *,
                                                                    size_t,
                                                                    int op))
{
  return PLL_FAILURE;
}

PLL_EXPORT int pllmod_networkinfo_init_partition(pllmod_networkinfo_t * networkinfo,
                                           unsigned int partition_index,
                                           pll_partition_t * partition,
                                           int params_to_optimize,
                                           int gamma_mode,
                                           double alpha,
                                           const unsigned int * param_indices,
                                           const int * subst_matrix_symmetries)
{
  return PLL_FAILURE;
}

PLL_EXPORT int pllmod_networkinfo_set_active_partition(pllmod_networkinfo_t * networkinfo,
                                                    int partition_index)
{
  return PLL_FAILURE;
}

PLL_EXPORT int pllmod_networkinfo_set_root(pllmod_networkinfo_t * networkinfo,
                                        pll_rnetwork_node_t * root)
{
  return PLL_FAILURE;
}

PLL_EXPORT
int pllmod_networkinfo_get_branch_length_all(const pllmod_networkinfo_t * networkinfo,
                                          const pll_rnetwork_node_t * edge,
                                          double * lengths)
{
  return PLL_FAILURE;
}

PLL_EXPORT int pllmod_networkinfo_set_branch_length(pllmod_networkinfo_t * networkinfo,
                                                 pll_rnetwork_node_t * edge,
                                                 double length)
{
  return PLL_FAILURE;
}

PLL_EXPORT
int pllmod_networkinfo_set_branch_length_all(pllmod_networkinfo_t * networkinfo,
                                          pll_rnetwork_node_t * edge,
                                          const double * lengths)
{
  return PLL_FAILURE;
}

PLL_EXPORT
int pllmod_networkinfo_set_branch_length_partition(pllmod_networkinfo_t * networkinfo,
                                                pll_rnetwork_node_t * edge,
                                                int partition_index,
                                                double length)
{
  return PLL_FAILURE;
}

PLL_EXPORT
pll_rnetwork_t * pllmod_networkinfo_get_partition_network(const pllmod_networkinfo_t * networkinfo,
                                                 int partition_index)
{
  return PLL_FAILURE;
}

PLL_EXPORT
pllmod_networkinfo_topology_t * pllmod_networkinfo_get_topology(const pllmod_networkinfo_t * networkinfo,
                                                          pllmod_networkinfo_topology_t * topol)
{
  return PLL_FAILURE;
}

PLL_EXPORT
int pllmod_networkinfo_set_topology(pllmod_networkinfo_t * networkinfo,
                                 const pllmod_networkinfo_topology_t * topol)
{
  return PLL_FAILURE;
}

PLL_EXPORT
int pllmod_networkinfo_destroy_topology(pllmod_networkinfo_topology_t * topol)
{
  return PLL_FAILURE;
}

PLL_EXPORT int pllmod_betworkinfo_destroy_partition(pllmod_networkinfo_t * networkinfo,
                                                 unsigned int partition_index)
{
  return PLL_FAILURE;
}

PLL_EXPORT void pllmod_networkinfo_destroy(pllmod_networkinfo_t * networkinfo)
{
  return;
}

PLL_EXPORT int pllmod_networkinfo_update_prob_matrices(pllmod_networkinfo_t * networkinfo,
                                                    int update_all)
{
  return PLL_FAILURE;
}

PLL_EXPORT void pllmod_networkinfo_invalidate_all(pllmod_networkinfo_t * networkinfo)
{
  return;
}

PLL_EXPORT int pllmod_networkinfo_validate_clvs(pllmod_networkinfo_t * networkinfo,
                                             pll_rnetwork_node_t ** travbuffer,
                                             unsigned int travbuffer_size)
{
  return PLL_FAILURE;
}

PLL_EXPORT void pllmod_networkinfo_invalidate_pmatrix(pllmod_networkinfo_t * networkinfo,
                                                   const pll_rnetwork_node_t * edge)
{
  return;
}

PLL_EXPORT void pllmod_networkinfo_invalidate_clv(pllmod_networkinfo_t * networkinfo,
                                               const pll_unode_t * edge)
{
  return;
}

PLL_EXPORT double pllmod_networkinfo_compute_loglh(pllmod_networkinfo_t * networkinfo,
                                                int incremental)
{
  return PLL_FAILURE;
}

PLL_EXPORT double pllmod_networkinfo_compute_loglh_flex(pllmod_networkinfo_t * networkinfo,
                                                     int incremental,
                                                     int update_pmatrices)
{
  return PLL_FAILURE;
}

PLL_EXPORT
int pllmod_networkinfo_scale_branches_all(pllmod_networkinfo_t * networkinfo, double scaler)
{
  return PLL_FAILURE;
}

PLL_EXPORT
int pllmod_networkinfo_scale_branches_partition(pllmod_networkinfo_t * networkinfo,
                                             unsigned int partition_idx,
                                             double scaler)
{
  return PLL_FAILURE;
}

PLL_EXPORT
int pllmod_networkinfo_normalize_brlen_scalers(pllmod_networkinfo_t * networkinfo)
{
  return PLL_FAILURE;
}

PLL_EXPORT int pllmod_networkinfo_set_network(pllmod_networkinfo_t * networkinfo,
                                        pll_rnetwork_t * network)
{
  return PLL_FAILURE;
}

PLL_EXPORT int pllmod_networkinfo_set_constraint_clvmap(pllmod_networkinfo_t * networkinfo,
                                                     const int * clv_index_map)
{
  return PLL_FAILURE;
}

PLL_EXPORT int pllmod_networkinfo_set_constraint_network(pllmod_networkinfo_t * networkinfo,
                                                   const pll_rnetwork_t * cons_network)
{
  return PLL_FAILURE;
}

PLL_EXPORT int pllmod_networkinfo_check_constraint(pllmod_networkinfo_t * networkinfo,
                                                pll_rnetwork_node_t * subnetwork,
                                                pll_rnetwork_node_t * regraft_edge)
{
  return PLL_FAILURE;
}
