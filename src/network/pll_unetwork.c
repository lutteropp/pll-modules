/*
 Copyright (C) 2016 Diego Darriba, extended to networks 2019 by Sarah Lutteropp

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

 /**
  * @file pll_network.c
  *
  * @brief Operations on network structures
  *
  * @author Diego Darriba
  */

#include "pll_network.h"
#include "../tree/tree_hashtable.h"
#include "../pllmod_common.h"

#define UNIMPLEMENTED 0

static int rnetwork_rollback_tbr(pll_network_rollback_t * rollback_info);
static int rnetwork_rollback_spr(pll_network_rollback_t * rollback_info);
static int rnetwork_rollback_nni(pll_network_rollback_t * rollback_info);
static int unetwork_rollback_tbr(pll_network_rollback_t * rollback_info);
static int unetwork_rollback_spr(pll_network_rollback_t * rollback_info);
static int unetwork_rollback_nni(pll_network_rollback_t * rollback_info);
static int unetwork_find_node_in_subnetwork(pll_unetwork_node_t * root, pll_unetwork_node_t * node);
static int cb_update_matrices_partials(pll_unetwork_node_t * node, void *data);
static void shuffle_network_nodes(const pll_unetwork_t * network, unsigned int seed);
static void split_multi_node(pll_unetwork_t * network, pll_unetwork_node_t * first,
                             pll_unetwork_node_t * last, unsigned int degree);
static char * default_support_fmt(double support);

struct cb_params
{
  const unsigned int * params_indices;
  pll_partition_t * partition;
  int update_pmatrices;
  int update_partials;
};

/******************************************************************************/
/* Topological rearrangements */

/**
 * Performs one TBR move by applying a bisection and a reconnection.
 * The CLV, scaler and pmatrix indices are updated.
 *
 * @param[in] b_edge bisection point
 * @param[in] r_edge reconnection point
 * @param[out] rollback_info Rollback information for undoing this move.
 *                           If it is NULL, rollback information is ignored.
 *
 * @return PLL_SUCCESS if the move was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_unetwork_tbr(pll_unetwork_node_t * b_edge,
                                pll_network_edge_t * r_edge,
                                pll_network_rollback_t * rollback_info)
{
  pll_unetwork_node_t *parent, *child;

  /* validate if the move can be applied */

  /* 1. bisection point must not be a leaf branch */
  if (!(b_edge->next && b_edge->back->next))
  {
    pllmod_set_error(PLLMOD_NETWORK_ERROR_TBR_LEAF_BISECTION,
                     "attempting to bisect at a leaf node");
    return PLL_FAILURE;
  }

  /* 2. reconnection edges are different from bisection point */
  if (b_edge == r_edge->edge.unetwork.parent ||
      b_edge == r_edge->edge.unetwork.parent->back ||
      b_edge == r_edge->edge.unetwork.child ||
      b_edge == r_edge->edge.unetwork.child->back ||
      b_edge->back == r_edge->edge.unetwork.parent ||
      b_edge->back == r_edge->edge.unetwork.parent->back ||
      b_edge->back == r_edge->edge.unetwork.child ||
      b_edge->back == r_edge->edge.unetwork.child->back)
  {
    pllmod_set_error(PLLMOD_NETWORK_ERROR_TBR_OVERLAPPED_NODES,
                     "TBR nodes are overlapped");
    return PLL_FAILURE;
  }

  /* 3. reconnection edges must belong to different subnetworks rooted at b_edge
   *    and b_edge->back
   */
  if (!(unetwork_find_node_in_subnetwork(b_edge, r_edge->edge.unetwork.parent) &&
        unetwork_find_node_in_subnetwork(b_edge->back, r_edge->edge.unetwork.child)) &&
      !(unetwork_find_node_in_subnetwork(b_edge->back, r_edge->edge.unetwork.parent) &&
        unetwork_find_node_in_subnetwork(b_edge, r_edge->edge.unetwork.child)))
  {
    pllmod_set_error(PLLMOD_NETWORK_ERROR_TBR_SAME_SUBNETWORK,
                     "TBR reconnection in same subnetwork");
    return PLL_FAILURE;
  }

  /* save rollback information */
  if (rollback_info)
  {
    rollback_info->rearrange_type     = PLLMOD_NETWORK_REARRANGE_TBR;
    rollback_info->rooted             = 0;
    rollback_info->TBR.bisect_edge    = (void *) b_edge;
    rollback_info->TBR.reconn_edge.edge.unetwork.parent = b_edge->next->next;
    rollback_info->TBR.reconn_edge.edge.unetwork.child  = b_edge->back->next->next;
    rollback_info->TBR.reconn_edge.length = b_edge->length;

    rollback_info->TBR.bisect_left_bl = r_edge->edge.unetwork.parent->length;
    rollback_info->TBR.bisect_right_bl = r_edge->edge.unetwork.child->length;

    rollback_info->TBR.reconn_parent_left_bl  = b_edge->next->length;
    rollback_info->TBR.reconn_parent_right_bl = b_edge->next->next->length;
    rollback_info->TBR.reconn_child_left_bl   = b_edge->back->next->length;
    rollback_info->TBR.reconn_child_right_bl  = b_edge->back->next->next->length;
  }

  /* bisect at b_edge */
  pllmod_unetwork_bisect(b_edge, &parent, &child);

  /* reconnect at r_edge */
  pllmod_unetwork_reconnect(r_edge,
                      b_edge);

  return PLL_SUCCESS;
}



/**
 * Performs one SPR move
 * The CLV, scaler and pmatrix indices are updated.
 *
 * @param[in] p_edge Edge to be pruned
 * @param[in] r_edge Edge to be regrafted
 * @param[out] rollback_info Rollback information for undoing this move.
 *                           If it is NULL, rollback information is ignored.
 *
 * @return PLL_SUCCESS if the move was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_unetwork_spr(pll_unetwork_node_t * p_edge,
                                pll_unetwork_node_t * r_edge,
                                pll_network_rollback_t * rollback_info)
{
  int retval;

  if (pllmod_unetwork_is_tip(p_edge))
  {
    /* invalid move */
    pllmod_set_error(PLLMOD_NETWORK_ERROR_SPR_INVALID_NODE,
                     "Attempting to prune a leaf branch");
    return PLL_FAILURE;
  }

  /* save rollback information */
  if (rollback_info)
  {
    rollback_info->rearrange_type     = PLLMOD_NETWORK_REARRANGE_SPR;
    rollback_info->rooted             = 0;
    rollback_info->SPR.prune_edge     = (void *) p_edge;
    rollback_info->SPR.regraft_edge   = (void *) p_edge->next->back;
    rollback_info->SPR.prune_bl       = p_edge->length;
    rollback_info->SPR.prune_left_bl  = p_edge->next->length;
    rollback_info->SPR.prune_right_bl = p_edge->next->next->length;
    rollback_info->SPR.regraft_bl     = r_edge->length;
  }

  retval = pll_unetwork_spr(p_edge,
                         r_edge,
                         0, 0, 0);

  return retval;
}

/**
 * Performs one NNI move
 * The CLV, scaler and pmatrix indices are updated.
 *
 * @param[in] edge NNI interchange edge
 * @param[in] type move type: PLL_NNI_LEFT, PLL_NNI_RIGHT
 * @param[out] rollback_info Rollback information for undoing this move.
 *                           If it is NULL, rollback information is ignored.
 *
 * @return PLL_SUCCESS if the move was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_unetwork_nni(pll_unetwork_node_t * edge,
                                int type,
                                pll_network_rollback_t * rollback_info)
{
  /* validate preconditions */
  assert(edge && edge->back);

  if (!(type == PLL_UNETWORK_MOVE_NNI_LEFT || type == PLL_UNETWORK_MOVE_NNI_RIGHT))
  {
    /* invalid move */
    pllmod_set_error(PLLMOD_NETWORK_ERROR_NNI_INVALID_MOVE,
                     "Invalid NNI move type");
    return PLL_FAILURE;
  }
  if (pllmod_unetwork_is_tip(edge) || pllmod_unetwork_is_tip(edge->back))
  {
    /* invalid move */
    pllmod_set_error(PLLMOD_NETWORK_ERROR_INTERCHANGE_LEAF,
                     "Attempting to apply NNI on a leaf branch");
    return PLL_FAILURE;
  }

  /* save rollback information */
  if (rollback_info)
  {
    rollback_info->rearrange_type     = PLLMOD_NETWORK_REARRANGE_NNI;
    rollback_info->rooted             = 0;
    rollback_info->NNI.edge           = (void *) edge;
    rollback_info->NNI.type           = type;
    rollback_info->NNI.left_left_bl   = edge->next->length;
    rollback_info->NNI.left_right_bl  = edge->next->next->length;
    rollback_info->NNI.right_left_bl  = edge->back->next->length;
    rollback_info->NNI.right_right_bl = edge->back->next->next->length;
    rollback_info->NNI.edge_bl        = edge->length;
  }

  if (!pllmod_unetwork_nni(edge, type, 0))
    return PLL_FAILURE;

  return PLL_SUCCESS;
}

/**
 * Rollback the previous move
 * @param  rollback_info the rollback info returned by the previous move
 * @return PLL_SUCCESS if the rollback move was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_network_rollback(pll_network_rollback_t * rollback_info)
{
  int retval = PLL_FAILURE;
  switch (rollback_info->rearrange_type)
  {
    case PLLMOD_NETWORK_REARRANGE_TBR:
      {
        if (rollback_info->rooted)
          retval = rnetwork_rollback_tbr (rollback_info);
        else
          retval = unetwork_rollback_tbr (rollback_info);
      }
      break;
    case PLLMOD_NETWORK_REARRANGE_SPR:
      {
        if (rollback_info->rooted)
          retval = rnetwork_rollback_spr (rollback_info);
        else
          retval = unetwork_rollback_spr (rollback_info);
      }
      break;
    case PLLMOD_NETWORK_REARRANGE_NNI:
      {
        if (rollback_info->rooted)
          retval = rnetwork_rollback_nni (rollback_info);
        else
          retval = unetwork_rollback_nni (rollback_info);
      }
      break;
    default:
      /* unimplemented */
      assert(0);
      break;
  }
  return retval;
}



/******************************************************************************/
/* Tree construction */

PLL_EXPORT pll_unetwork_t * pllmod_unetwork_resolve_multi(const pll_unetwork_t * multi_network,
                                                    unsigned int random_seed,
                                                    int * clv_index_map)
{
  if (!multi_network)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Parameter multi_network is NULL.");
    return NULL;
  }

  if (multi_network->vroot->next &&
      multi_network->vroot->next->next == multi_network->vroot)
  {
    pllmod_set_error(PLLMOD_NETWORK_ERROR_INVALID_NETWORK,
                     "Unrooted network is expected but a rooted network was provided.");
    return NULL;
  }

  pll_unetwork_t * bin_network = pll_unetwork_clone(multi_network);

  unsigned int tip_count = bin_network->tip_count;
  unsigned int multi_node_count = bin_network->tip_count + bin_network->inner_tree_count + bin_network->reticulation_count;
  unsigned int bin_node_count = 2 * tip_count -2;

  // 1:1 CLV index mapping for existing nodes
  if (clv_index_map)
  {
    for (unsigned int i = 0; i < multi_node_count; ++i)
    {
      const unsigned int clv_id = bin_network->nodes[i]->clv_index;
      clv_index_map[clv_id] = clv_id;
    }
  }

  if (bin_network->binary)
    return bin_network;

  if (random_seed)
    shuffle_network_nodes(bin_network, random_seed);

  bin_network->nodes = (pll_unetwork_node_t **)realloc(bin_network->nodes,
                                            bin_node_count*sizeof(pll_unetwork_node_t *));

  // iterate over inner nodes, resolve multifurcations, and map new->old CLV indices
  unsigned int old_inner_count = bin_network->inner_tree_count + bin_network->reticulation_count;
  for (unsigned int i = tip_count; i < multi_node_count; ++i)
  {
    pll_unetwork_node_t * start = bin_network->nodes[i];
    pll_unetwork_node_t * end = NULL;
    pll_unetwork_node_t * snode = start;
    unsigned int degree = 0;
    do
    {
      end = snode;
      snode = snode->next;
      degree++;
    }
    while (snode && snode != start);

    split_multi_node(bin_network, start, end, degree);

    assert(bin_network->inner_tree_count + bin_network->reticulation_count == old_inner_count + degree-3);
    if (clv_index_map)
    {
      for (unsigned int j = old_inner_count; j < bin_network->inner_tree_count + bin_network->reticulation_count; ++j)
      {
        unsigned int new_clv_id = bin_network->nodes[tip_count + j]->clv_index;
        clv_index_map[new_clv_id] = start->clv_index;
      }
    }
    old_inner_count = bin_network->inner_tree_count + bin_network->reticulation_count;
  }

  assert(bin_network->inner_tree_count == bin_network->tip_count - 2 + bin_network->reticulation_count);

  bin_network->binary = 1;

  /* re-assign node indices such that:
   * (1) all 3 pll_unetwork_node's of an inner node have consecutive indices: (x, x+1, x+2)
   * (2) for any two random multifurcation resolutions R1 and R2 holds
   *      (x, x+1, x+2) in R1 iff (x, x+1, x+2) in R2
   */
  unsigned int max_node_index = tip_count;
  for (unsigned int i = tip_count; i < bin_node_count; ++i)
  {
    pll_unetwork_node_t * node = bin_network->nodes[i];
    node->node_index = max_node_index++;
    node->next->node_index = max_node_index++;
    node->next->next->node_index = max_node_index++;
  }
  assert(max_node_index == bin_network->tip_count + 3*(bin_network->inner_tree_count + bin_network->reticulation_count));

  return bin_network;
}

PLL_EXPORT int pllmod_unetwork_root_inplace(pll_unetwork_t * network)
{
  if (!network)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Empty network specified!");
    return PLL_FAILURE;
  }

  /* check if network is already rooted */
  if (network->vroot->next && network->vroot->next->next == network->vroot)
    return PLL_SUCCESS;

  pll_unetwork_node_t * root = network->vroot;
  pll_unetwork_node_t * root_back = root->back;
  pll_unetwork_node_t * root_left = (pll_unetwork_node_t *) calloc(1, sizeof(pll_unetwork_node_t));
  pll_unetwork_node_t * root_right = (pll_unetwork_node_t *) calloc(1, sizeof(pll_unetwork_node_t));
  root_left->next = root_right;
  root_right->next = root_left;
  double root_brlen = root->length / 2.;
  unsigned int last_clv_index = 0;
  int last_scaler_index = 0;
  unsigned int last_node_index = 0;
  unsigned int last_pmatrix_index = 0;
  unsigned int node_count = network->inner_tree_count + network->reticulation_count + network->tip_count;

  for (unsigned int i = 0; i < node_count; ++i)
  {
    const pll_unetwork_node_t * node  = network->nodes[i];
    last_clv_index = PLL_MAX(last_clv_index, node->clv_index);
    last_scaler_index = PLL_MAX(last_scaler_index, node->scaler_index);
    do
    {
      last_node_index = PLL_MAX(last_node_index, node->node_index);
      last_pmatrix_index = PLL_MAX(last_pmatrix_index, node->pmatrix_index);
      node = node->next;
    }
    while(node && node != network->nodes[i]);
  }

  root_left->clv_index = root_right->clv_index = ++last_clv_index;
  root_left->scaler_index = root_right->scaler_index = ++last_scaler_index;
  root_left->node_index = ++last_node_index;
  root_right->node_index = ++last_node_index;
  root_right->pmatrix_index = ++last_pmatrix_index;

  pllmod_unetwork_connect_nodes(root, root_left, root_brlen, 1.0);
  pllmod_unetwork_connect_nodes(root_right, root_back, root_brlen, 1.0);

  network->vroot = root_left;
  network->inner_tree_count++;
  network->edge_count++;
  node_count++;

  network->nodes = (pll_unetwork_node_t **) realloc(network->nodes,
                                         node_count*sizeof(pll_unetwork_node_t *));
  network->nodes[node_count-1] = root_left;

  return PLL_SUCCESS;
}

PLL_EXPORT int pllmod_unetwork_outgroup_root(pll_unetwork_t * network,
                                          unsigned int * outgroup_tip_ids,
                                          unsigned int outgroup_size,
                                          int add_root_node)
{
  pll_unetwork_node_t ** split_to_node_map = NULL;
  pll_split_t * network_splits = NULL;
  pll_unetwork_node_t * new_root = NULL;
  unsigned int tip_count;
  unsigned int split_count;

  if (!network || !outgroup_tip_ids || !outgroup_size)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
                     "Empty network and/or outgroup specified!");
    return PLL_FAILURE;
  }

  if (outgroup_size == 1)
  {
    // special case single-taxon outgroup: just find a tip by node_index
    for (unsigned int i = 0; i < network->tip_count; ++i)
    {
      const pll_unetwork_node_t * node  = network->nodes[i];
      if (node->node_index == outgroup_tip_ids[0])
      {
        new_root = node->back;
        break;
      }
    }
  }
  else
  {
    tip_count = network->tip_count;
    split_count = tip_count - 3;

    split_to_node_map = (pll_unetwork_node_t **) calloc(split_count,
                                                sizeof(pll_unetwork_node_t *));

    if (!split_to_node_map)
    {
      pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                       "Cannot allocate memory for split->node map!");
      return PLL_FAILURE;
    }

    network_splits = pllmod_unetwork_split_create(network->vroot, network->tip_count, network->reticulation_count,
                                            split_to_node_map);

    if (!network_splits)
    {
      assert(pll_errno);
      free(split_to_node_map);
      return PLL_FAILURE;
    }

    // create outgroup split
    pll_split_t outgroup_split = pllmod_unetwork_split_from_tips(outgroup_tip_ids,
                                                              outgroup_size,
                                                              tip_count);

    // check if this split is in the network
    size_t split_len = bitv_length(tip_count);
    for (unsigned int i = 0; i < split_count; ++i)
    {
      if (!bitv_compare(network_splits[i], outgroup_split, split_len))
      {
        new_root = split_to_node_map[i];
        break;
      }
    }

    pllmod_unetwork_split_destroy(network_splits);
    free(split_to_node_map);
    free(outgroup_split);
  }

  // set network->vroot to the outgroup split node
  if (new_root)
  {
    network->vroot = new_root;
    if (add_root_node)
      return pllmod_unetwork_root_inplace(network);
    else
      return PLL_SUCCESS;
  }
  else
  {
    pllmod_set_error(PLLMOD_NETWORK_ERROR_POLYPHYL_OUTGROUP,
                     "Outgroup is not monophyletic!");
    return PLL_FAILURE;
  }
}


int unetwork_insert_tips_random(pll_unetwork_node_t ** nodes, unsigned int taxa_count,
                             unsigned int start_tip, unsigned int random_seed)
{
  unsigned int i;
  unsigned int start_inner_count     = start_tip - 2;
  unsigned int start_branches        = 2 * start_tip - 3;
  unsigned int max_branches          = 2 * taxa_count - 3;
  unsigned int placed_branches_count = 0;
  unsigned int last_branch_id        = 0;

  pll_unetwork_node_t ** branches   = NULL;
  pll_random_state * rstate = NULL;

  branches = (pll_unetwork_node_t **) calloc(max_branches, sizeof(pll_unetwork_node_t *));

  if (!branches)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for branches!");
    return PLL_FAILURE;
  }

  rstate =  pll_random_create(random_seed);

  if (!rstate)
  {
    free(branches);
    return PLL_FAILURE;
  }

  // check pmatrix indices on tip branches
  for (i = 0; i < taxa_count; ++i)
    last_branch_id = PLL_MAX(last_branch_id, nodes[i]->pmatrix_index);

  for (i = taxa_count; i < taxa_count + start_inner_count; ++i)
  {
    pll_unetwork_node_t * snode = nodes[i];
    do
    {
      if (snode->clv_index > snode->back->clv_index)
      {
        branches[placed_branches_count++] = snode;
        last_branch_id = PLL_MAX(last_branch_id, snode->pmatrix_index);
      }
      snode = snode->next;
    }
    while (snode != nodes[i]);
  }
  assert(placed_branches_count == start_branches);

  for (i = start_tip; i < taxa_count; ++i)
  {
    /* take tips iteratively */
    pll_unetwork_node_t * next_tip = nodes[i];
    pll_unetwork_node_t * next_inner = nodes[taxa_count + i - 2];

    /* select random branch from the network */
    unsigned int rand_branch_id = pll_random_getint(rstate, placed_branches_count);
    pll_unetwork_node_t * next_branch = branches[rand_branch_id];

    /* connect tip to selected branch */
    if (next_branch->incoming)
    {
      pllmod_unetwork_connect_nodes(next_branch->back, next_inner,
                                   PLLMOD_NETWORK_DEFAULT_BRANCH_LENGTH, 1.0);
      pllmod_unetwork_connect_nodes(next_inner->next, next_branch,
                                   PLLMOD_NETWORK_DEFAULT_BRANCH_LENGTH, 1.0);
      pllmod_unetwork_connect_nodes(next_inner->next->next, next_tip,
                                   PLLMOD_NETWORK_DEFAULT_BRANCH_LENGTH, 1.0);

      assert(next_inner->incoming == 1);
      assert(next_inner->next->incoming == 0);
      assert(next_inner->next->next->incoming == 0);
    }
    else
    {
      pllmod_unetwork_connect_nodes(next_inner, next_branch->back,
    	                            PLLMOD_NETWORK_DEFAULT_BRANCH_LENGTH, 1.0);
      pllmod_unetwork_connect_nodes(next_branch, next_inner->next,
    	                            PLLMOD_NETWORK_DEFAULT_BRANCH_LENGTH, 1.0);
      pllmod_unetwork_connect_nodes(next_inner->next->next, next_tip,
    	                            PLLMOD_NETWORK_DEFAULT_BRANCH_LENGTH, 1.0);

      assert(next_inner->incoming == 0);
      assert(next_inner->next->incoming == 1);
      assert(next_inner->next->next->incoming == 0);
    }

    if (pllmod_unetwork_is_tip (next_inner->back))
    {
      next_inner->next->pmatrix_index = next_inner->next->back->pmatrix_index =
          ++last_branch_id;
    }
    else
    {
      next_inner->pmatrix_index = next_inner->back->pmatrix_index =
          ++last_branch_id;
    }

    /* store branches */
    branches[placed_branches_count++] = next_inner;
    branches[placed_branches_count++] = next_inner->next->next;
  }
  assert(placed_branches_count == max_branches);

  /* clean */
  free (branches);
  pll_random_destroy(rstate);

  return PLL_SUCCESS;
}

/**
 * Extend a network by inserting new taxa to randomly chosen branches
 */
PLL_EXPORT int pllmod_unetwork_extend_random(pll_unetwork_t * network,
                                          unsigned int ext_taxa_count,
                                          const char * const* ext_names,
                                          unsigned int random_seed)
{
  unsigned int old_taxa_count  = network->tip_count;
  unsigned int old_inner_tree_count = network->inner_tree_count;
  unsigned int old_reticulation_count = network->reticulation_count;
  unsigned int old_node_count  = old_taxa_count + old_inner_tree_count + old_reticulation_count;
  unsigned int new_taxa_count  = old_taxa_count + ext_taxa_count;
  unsigned int new_inner_tree_count = old_inner_tree_count + ext_taxa_count;
  unsigned int new_reticulation_count = old_reticulation_count;
  unsigned int new_node_count = new_taxa_count + new_inner_tree_count + new_reticulation_count;

  unsigned int last_clv_id     = 0;
  unsigned int last_pmatrix_id = 0;
  unsigned int last_node_id    = 0;
  int next_scaler_id  = 0;

  unsigned int i;
  int retval;

  pll_unetwork_node_t ** old_nodes = network->nodes;
  pll_unetwork_node_t ** new_nodes = (pll_unetwork_node_t **) calloc(new_node_count,
                                                     sizeof(pll_unetwork_node_t *));

  if (!new_nodes)
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for nodes!");
    return PLL_FAILURE;
  }

  // 1:1 mapping for old tips
  for (i = 0; i < old_taxa_count; ++i)
    new_nodes[i] = old_nodes[i];

  // copy old inner nodes and adjust clvs
  for (i = old_taxa_count; i < old_node_count; ++i)
  {
    unsigned int new_idx = i + ext_taxa_count;
    new_nodes[new_idx] = old_nodes[i];
    pll_unetwork_node_t * snode = new_nodes[new_idx];
    assert(snode->next);
    do
    {
      snode->clv_index += ext_taxa_count;
      snode->node_index += ext_taxa_count;
      last_clv_id = PLL_MAX(last_clv_id, snode->clv_index);
      last_node_id = PLL_MAX(last_node_id, snode->node_index);
      next_scaler_id = PLL_MAX(next_scaler_id, snode->scaler_index);
      last_pmatrix_id = PLL_MAX(last_pmatrix_id, snode->pmatrix_index);
      snode = snode->next;
    }
    while (snode != new_nodes[new_idx]);
  }

  // create new tip nodes
  for (i = old_taxa_count; i < new_taxa_count; ++i)
  {
    pll_unetwork_node_t * node = (pll_unetwork_node_t *) calloc(1, sizeof(pll_unetwork_node_t));
    node->clv_index = i;
    node->node_index = i;
    node->scaler_index = PLL_SCALE_BUFFER_NONE;
    node->pmatrix_index = ++last_pmatrix_id; // ????

    node->label = ext_names ? strdup(ext_names[i - old_taxa_count]) : NULL;

    new_nodes[i] = node;
  }

  // create new inner nodes
  for (i = old_node_count + ext_taxa_count; i < new_node_count; ++i)
  {
    pll_unetwork_node_t * node = pllmod_unetwork_create_node(++last_clv_id,
                                                  ++next_scaler_id,
                                                  NULL, NULL);

    node->node_index = ++last_node_id;
    node->next->node_index = ++last_node_id;
    node->next->next->node_index = ++last_node_id;

    new_nodes[i] = node;
  }

  retval = unetwork_insert_tips_random(new_nodes, new_taxa_count,
                                    old_taxa_count, random_seed);

  if (retval)
  {
    free(network->nodes);
    network->nodes = new_nodes;
    network->tip_count = new_taxa_count;
    network->inner_tree_count = new_inner_tree_count;
    network->reticulation_count = new_reticulation_count;
    network->edge_count += 2 * ext_taxa_count;
    return PLL_SUCCESS;
  }
  else
  {
    free(new_nodes);
    return PLL_FAILURE;
  }
}

/**
 * Creates a random topology with default branch lengths
 */
PLL_EXPORT pll_unetwork_t * pllmod_unetwork_create_random(unsigned int taxa_count,
                                                    const char * const* names,
                                                    unsigned int random_seed)
{
  /*
   * The algorithm works as follows:
   *    1. Build a minimal 3-tip network
   *    2. Select a branch at random
   *    3. Connect next tip to that branch
   *    4. Repeat 2 and 3 until no tips left
   */
  unsigned int i;
  unsigned int tip_node_count        = taxa_count;
  unsigned int inner_node_count      = taxa_count - 2;
  unsigned int node_count            = tip_node_count + inner_node_count;

  pll_unetwork_node_t ** nodes    = (pll_unetwork_node_t **) calloc(node_count,
                                                    sizeof(pll_unetwork_node_t *));

  pll_unetwork_node_t * network_root;

  pll_unetwork_t * wrapped_network;

  unsigned int node_id = 0;

  /* allocate tips */
  for (i=0; i<taxa_count; ++i)
  {
    nodes[i] = (pll_unetwork_node_t *)calloc(1, sizeof(pll_unetwork_node_t));
    nodes[i]->clv_index = i;
    nodes[i]->scaler_index = PLL_SCALE_BUFFER_NONE;
    nodes[i]->pmatrix_index = i;
    nodes[i]->node_index = node_id++;

    if (names)
    {
      nodes[i]->label = (char *) malloc( strlen(names[i]) + 1 );
      strcpy(nodes[i]->label, names[i]);
    }
    else
    {
      nodes[i]->label = NULL;
    }
  }

  /* allocate inner */
  for (i=taxa_count; i<node_count; ++i)
  {
    nodes[i] = pllmod_unetwork_create_node(i, (int)i, NULL, NULL);
    nodes[i]->scaler_index -= taxa_count;
    nodes[i]->next->scaler_index -= taxa_count;
    nodes[i]->next->next->scaler_index -= taxa_count;

    nodes[i]->node_index = node_id++;
    nodes[i]->next->node_index = node_id++;
    nodes[i]->next->next->node_index = node_id++;
  }
  assert(node_id == tip_node_count + inner_node_count * 3);

  /* set an inner node as return value */
  network_root = nodes[taxa_count];

  /* build minimal network with 3 tips and 1 inner node */
  pllmod_unetwork_connect_nodes(nodes[taxa_count], nodes[0],
                             PLLMOD_NETWORK_DEFAULT_BRANCH_LENGTH, 1.0);
  pllmod_unetwork_connect_nodes(nodes[taxa_count]->next, nodes[1],
                             PLLMOD_NETWORK_DEFAULT_BRANCH_LENGTH, 1.0);
  pllmod_unetwork_connect_nodes(nodes[taxa_count]->next->next, nodes[2],
                             PLLMOD_NETWORK_DEFAULT_BRANCH_LENGTH, 1.0);

  /* insert remaining taxa_count-3 tips into the network */
  unetwork_insert_tips_random(nodes, taxa_count, 3, random_seed);

  /* clean */
  free (nodes);

  wrapped_network = pll_unetwork_wrapnetwork(network_root, tip_node_count);
  return (wrapped_network);
}

/**
 * Creates a maximum parsimony topology using randomized stepwise-addition
 * algorithm. All branch lengths will be set to default.
 */
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
                                            unsigned int * score)
{
  size_t i;
  pll_unetwork_t * network = NULL;

  pll_partition_t * partition = pll_partition_create(taxon_count,
                                                     0,   /* number of CLVs */
                                                     states,
                                                     seq_length,
                                                     1,
                                                     1, /* pmatrix count */
                                                     1,  /* rate_cats */
                                                     0,  /* scale buffers */
                                                     attributes);

  if (!partition)
  {
    assert(pll_errno);
    return NULL;
  }

  /* set pattern weights and free the weights array */
  if (site_weights)
    pll_set_pattern_weights(partition, site_weights);

  /* find sequences in hash table and link them with the corresponding taxa */
  for (i = 0; i < taxon_count; ++i)
    pll_set_tip_states(partition, i, map, sequences[i]);

  network = pllmod_unetwork_create_parsimony_multipart(taxon_count,
                                                 names,
                                                 1,
                                                 &partition,
                                                 random_seed,
                                                 score);

  /* destroy all structures allocated for the concrete PLL partition instance */
  pll_partition_destroy(partition);

  return network;
}

/**
 * Creates a maximum parsimony topology using randomized stepwise-addition
 * algorithm. All branch lengths will be set to default.
 * This function can be used with partitioned alignments (e.g., combined DNA+AA data)
 */
PLL_EXPORT
pll_unetwork_t * pllmod_unetwork_create_parsimony_multipart(unsigned int taxon_count,
                                                      char * const * taxon_names,
                                                      unsigned int partition_count,
                                                      pll_partition_t * const * partitions,
                                                      unsigned int random_seed,
                                                      unsigned int * score)
{
  pll_unetwork_t * network = NULL;
  unsigned int i;

  pll_parsimony_t ** parsimony =
      (pll_parsimony_t **) calloc(partition_count, sizeof(pll_parsimony_t *));

  if (!parsimony)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return NULL;
  }

  for (i = 0; i < partition_count; ++i)
  {
    assert(taxon_count == partitions[i]->tips);
    parsimony[i] = pll_fastparsimony_init(partitions[i]);
    if (!parsimony[i])
    {
      assert(pll_errno);
      goto cleanup;
    }
  }

  network = pll_fastparsimony_stepwise_network(parsimony,
                                    taxon_names,
                                    score,
                                    partition_count,
                                    random_seed);

  if (network)
  {
    /* set default branch lengths */
    pllmod_unetwork_set_length_recursive(network,
                                      PLLMOD_NETWORK_DEFAULT_BRANCH_LENGTH,
                                      0);
  }
  else
    assert(pll_errno);

cleanup:
  /* destroy parsimony */
  for (i = 0; i < partition_count; ++i)
  {
    if (parsimony[i])
      pll_parsimony_destroy(parsimony[i]);
  }

  free(parsimony);

  return network;
}

/* static functions */

static int unetwork_find_node_in_subnetwork(pll_unetwork_node_t * root,
                                                 pll_unetwork_node_t * node)
{
  if (root == node)
  {
    return PLL_SUCCESS;
  }

  if (root->next)
  {
    if (root->next == node || root->next->next == node)
    {
      return PLL_SUCCESS;
    }

    return unetwork_find_node_in_subnetwork(root->next->back, node)
        || unetwork_find_node_in_subnetwork(root->next->next->back, node);
  }

  return PLL_FAILURE;
}



/******************************************************************************/
/* Additional utilities */

static int unetwork_traverse_apply(pll_unetwork_node_t * node,
                                int (*cb_pre_trav)(pll_unetwork_node_t *, void *),
                                int (*cb_in_trav)(pll_unetwork_node_t *, void *),
                                int (*cb_post_trav)(pll_unetwork_node_t *, void *),
                                void *data)
{
  int retval = 1;
  pll_unetwork_node_t * child_network = 0;

  if (cb_pre_trav && !cb_pre_trav(node,  data))
    return PLL_FAILURE;

  if (pllmod_unetwork_is_tip(node))
  {
    if (cb_in_trav)
      retval &= cb_in_trav(node, data);
    if (cb_post_trav)
      retval &= cb_post_trav(node, data);
    return retval;
  }

  child_network = node->next;
  while(child_network != node)
  {
    retval &= unetwork_traverse_apply(child_network->back,
                                   cb_pre_trav, cb_in_trav, cb_post_trav, data);

    if (cb_in_trav &&
        child_network->next != node &&
        !cb_in_trav(child_network, data))
      return PLL_FAILURE;

    child_network = child_network->next;
  }

  if (cb_post_trav)
    retval &= cb_post_trav(node,  data);

  return retval;
}

PLL_EXPORT int pllmod_unetwork_traverse_apply(pll_unetwork_node_t * root,
                                           int (*cb_pre_trav)(pll_unetwork_node_t *,
                                                               void *),
                                           int (*cb_in_trav)(pll_unetwork_node_t *,
                                                               void *),
                                           int (*cb_post_trav)(pll_unetwork_node_t *,
                                                               void *),
                                           void *data)
{
  int retval = 1;

  assert(root);

  if (pllmod_unetwork_is_tip(root)) return PLL_FAILURE;

  retval &= unetwork_traverse_apply(root->back,
                                 cb_pre_trav, cb_in_trav, cb_post_trav,
                                 data);
  retval &= unetwork_traverse_apply(root,
                                 cb_pre_trav, cb_in_trav, cb_post_trav,
                                 data);

  return retval;
}

PLL_EXPORT int pllmod_unetwork_is_tip(const pll_unetwork_node_t * node)
{
  return (node->next == NULL);
}

PLL_EXPORT void pllmod_unetwork_set_length(pll_unetwork_node_t * edge,
                                            double length)
{
  edge->length = edge->back->length = length;
}

PLL_EXPORT void pllmod_unetwork_set_length_recursive(pll_unetwork_t * network,
                                                  double length,
                                                  int missing_only)
{
  /* set branch lengths */
  unsigned int i;
  unsigned int tip_count = network->tip_count;
  unsigned int inner_count = network->inner_tree_count + network->reticulation_count;
  for (i = 0; i < tip_count + inner_count; ++i)
  {
    pll_unetwork_node_t * node = network->nodes[i];
    if (!node->length || !missing_only)
      pllmod_unetwork_set_length(node, length);
    if (node->next)
    {
      if (!node->next->length || !missing_only)
        pllmod_unetwork_set_length(node->next, length);
      if (!node->next->next->length || !missing_only)
        pllmod_unetwork_set_length(node->next->next, length);
    }
  }
}

PLL_EXPORT void pllmod_unetwork_scale_branches(pll_unetwork_t * network,
                                            double branch_length_scaler)
{
  /* scale branch lengths */
  unsigned int i;
  unsigned int tip_count = network->tip_count;
  unsigned int inner_count = network->inner_tree_count + network->reticulation_count;
  pll_unetwork_node_t ** nodes = network->nodes;
  for (i=0; i<tip_count; ++i)
  {
    nodes[i]->length *= branch_length_scaler;
  }
  for (i=tip_count; i<tip_count+inner_count; ++i)
  {
    nodes[i]->length *= branch_length_scaler;
    nodes[i]->next->length *= branch_length_scaler;
    nodes[i]->next->next->length *= branch_length_scaler;
  }
}

PLL_EXPORT void pllmod_unetwork_scale_branches_all(pll_unetwork_node_t * root,
                                                double branch_length_scaler)
{
  double root_length;

  /* scale all branches in a network */
  pllmod_unetwork_scale_subnetwork_branches(root, branch_length_scaler);
  root_length = root->length;
  pllmod_unetwork_scale_subnetwork_branches(root->back, branch_length_scaler);

  /* undo duplicated scaling */
  root->length = root->back->length = root_length;
}

PLL_EXPORT void pllmod_unetwork_scale_subnetwork_branches(pll_unetwork_node_t * root,
                                                    double branch_length_scaler)
{
 /* scale all branches in a subnetwork rooted at node */
  root->length *= branch_length_scaler;
  root->back->length *= branch_length_scaler;

  if (root->next)
  {
    pllmod_unetwork_scale_subnetwork_branches(root->next->back, branch_length_scaler);
    pllmod_unetwork_scale_subnetwork_branches(root->next->next->back, branch_length_scaler);
  }
}


/**
 * compute the likelihood on a unetwork structure
 * if update_pmatrices or update_partials are set, p-matrices and CLVs are
 * updated before computing the likelihood.
 */
PLL_EXPORT double pllmod_unetwork_compute_lk(pll_partition_t * partition,
                                       pll_unetwork_node_t * network,
                                       const unsigned int * params_indices,
                                       int update_pmatrices,
                                       int update_partials)
{
  struct cb_params parameters;
  assert (network);
  assert (network->pmatrix_index == network->back->pmatrix_index);

  parameters.partition      = partition;
  parameters.params_indices = params_indices;

  /* update pmatrices */
  if (update_pmatrices || update_partials)
  {
    parameters.update_pmatrices = update_pmatrices;
    parameters.update_partials  = update_partials;

    pllmod_unetwork_traverse_apply(network,
                                0,
                                0,
                                cb_update_matrices_partials,
                                (void *) &parameters);
  }

  double logl = pll_compute_edge_loglikelihood(partition,
                                              network->clv_index,
                                              network->scaler_index,
                                              network->back->clv_index,
                                              network->back->scaler_index,
                                              network->pmatrix_index,
                                              params_indices,
                                              NULL);
  return logl;
}

struct clv_set_data
{
  int * set_indices;
  unsigned int max_index;
  unsigned int tip_count;
};

static int cb_set_clv_minimal(pll_unetwork_node_t * node, void * data)
{
  unsigned int i, next_index;
  int index_found;
  struct clv_set_data * clv_data = (struct clv_set_data *)data;
  int * v = 0;

  if (!pllmod_unetwork_is_tip(node))
  {
    /* find next free position */
    v = clv_data->set_indices;
    next_index  = 0;
    index_found = 0;
    for (i=0; i<clv_data->max_index; ++i)
    {
      if (!v[i])
      {
        index_found = 1;
        next_index = i;
        v[i] = 1;
        break;
      }
    }
    assert(index_found);

    /* set clv index */
    node->clv_index =
      node->next->clv_index =
      node->next->next->clv_index =
       next_index + clv_data->tip_count;
    /* set scaler index */
    node->scaler_index =
       node->next->scaler_index =
       node->next->next->scaler_index =
        (int)(next_index + clv_data->tip_count);

    /* free indices from children */
    if (!pllmod_unetwork_is_tip(node->next->back))
    {
      v[node->next->back->clv_index - clv_data->tip_count] = 0;
    }
    if (!pllmod_unetwork_is_tip(node->next->next->back))
    {
      v[node->next->next->back->clv_index - clv_data->tip_count] = 0;
    }
  }

  /* continue */
  return 1;
}

PLL_EXPORT int pllmod_unetwork_set_clv_minimal(pll_unetwork_node_t * root,
                                            unsigned int tip_count)
{
  unsigned int clv_count = (unsigned int) ceil(log2(tip_count)) + 2;
  int * set_indices = (int *) calloc((size_t)clv_count, sizeof(int));
  struct clv_set_data data;
  data.set_indices = set_indices;
  data.max_index   = clv_count;
  data.tip_count   = tip_count;
  pllmod_unetwork_traverse_apply(root, 0, 0, cb_set_clv_minimal, (void *) &data);
  free(set_indices);

  return PLL_SUCCESS;
}

static int rnetwork_traverse_apply(pll_rnode_t * node,
                                int (*cb_pre_trav)(pll_rnode_t *, void *),
                                int (*cb_in_trav)(pll_rnode_t *, void *),
                                int (*cb_post_trav)(pll_rnode_t *, void *),
                                void *data)
{
  int retval = 1;

  if (cb_pre_trav && !cb_pre_trav(node,  data))
    return PLL_FAILURE;

  if (node->left)
  {
    retval &= rnetwork_traverse_apply(node->left,
                                   cb_pre_trav,
                                   cb_in_trav,
                                   cb_post_trav,
                                   data);

    if (cb_in_trav && !cb_in_trav(node,  data))
      return PLL_FAILURE;

    retval &= rnetwork_traverse_apply(node->right,
                                   cb_pre_trav,
                                   cb_in_trav,
                                   cb_post_trav,
                                   data);
  }

  if (cb_post_trav)
    retval &= cb_post_trav(node,  data);

  return retval;
}

PLL_EXPORT int pllmod_rnetwork_traverse_apply(pll_rnode_t * root,
                                           int (*cb_pre_trav)(pll_rnode_t *,
                                                               void *),
                                           int (*cb_in_trav)(pll_rnode_t *,
                                                               void *),
                                           int (*cb_post_trav)(pll_rnode_t *,
                                                               void *),
                                           void *data)
{
  int retval = 1;

  if (!root->left || !root->right) return PLL_FAILURE;

  retval &= rnetwork_traverse_apply(root,
                                 cb_pre_trav,
                                 cb_in_trav,
                                 cb_post_trav,
                                 data);

  return retval;
}

/* auxiliary structure for the callback function below */
struct serial_network_s {
  pll_unetwork_node_t * serialized_network;
  int node_count;
  int max_nodes;
};

/* callback function to fill the serialized network */
static int cb_serialize(pll_unetwork_node_t * network,
                        void * data)
{
  struct serial_network_s * list = (struct serial_network_s *) data;
  pll_unetwork_node_t * serialized_network = list->serialized_network;
  int cur_pos = list->node_count;

  assert(cur_pos < list->max_nodes);

  memcpy(&(serialized_network[cur_pos]), network, sizeof(pll_unetwork_node_t));
  serialized_network[cur_pos].data  = 0;
  serialized_network[cur_pos].label = 0;
  if (!pllmod_unetwork_is_tip(network))
  {
    /* set to arbitrary non-junk value */
    serialized_network[cur_pos].next = (pll_unetwork_node_t *) 1;
  }

  ++list->node_count;
  return 1;
}

//TODO: serialize/expand using a compressed format instead of pll_unetwork_node_t
PLL_EXPORT pll_unetwork_node_t * pllmod_unetwork_serialize(pll_unetwork_node_t * network,
                                                unsigned int tip_count, unsigned int reticulation_count)
{
  unsigned int node_count;
  pll_unetwork_node_t * serialized_network;
  struct serial_network_s data;

  node_count = 2*tip_count + reticulation_count - 2;

  /* allocate the serialized structure */
  serialized_network = (pll_unetwork_node_t *) malloc(node_count * sizeof (pll_unetwork_node_t));

  /* fill data for callback function */
  data.serialized_network = serialized_network;
  data.node_count      = 0;
  data.max_nodes       = node_count;

  /* if network is a tip, move to its back position */
  if (pllmod_unetwork_is_tip(network)) network = network->back;

  /* apply callback function to serialize */
  pllmod_unetwork_traverse_apply(network,
                              NULL,
                              NULL,
                              cb_serialize,
                              &data);

  if (data.node_count != data.max_nodes)
  {
    /* if the number of serialized nodes is not correct, return error */
    pllmod_set_error(PLLMOD_NETWORK_ERROR_INVALID_NETWORK,
                     "network structure ot tip_count are invalid");
    free(serialized_network);
    serialized_network = NULL;
  }

  return serialized_network;
}

PLL_EXPORT pll_unetwork_t * pllmod_unetwork_expand(pll_unetwork_node_t * serialized_network,
                                             unsigned int tip_count, unsigned int reticulation_count)
{
  unsigned int i, node_count, next_node_index;
  pll_unetwork_node_t ** network_stack;
  pll_unetwork_node_t * network;
  unsigned int network_stack_top;

  pllmod_reset_error();

  node_count  = 2*tip_count + reticulation_count - 2;

  /* allocate stack for at most 'n_tips' nodes */
  network_stack = (pll_unetwork_node_t **) malloc(tip_count * sizeof (pll_unetwork_node_t *));
  network_stack_top = 0;

  next_node_index = tip_count;

  /* read nodes */
  for (i=0; i<node_count; ++i)
  {
    pll_unetwork_node_t * t = 0;                   /* new node */
    pll_unetwork_node_t t_s = serialized_network[i];  /* serialized node */
    if (t_s.next)
    {
      /* build inner node and connect */
      pll_unetwork_node_t *t_cr, *t_r, *t_cl, *t_l;
      t = pllmod_unetwork_create_node(t_s.clv_index,
                                   t_s.scaler_index,
                                   0,  /* label */
                                   0); /* data */
      t_l = t->next;
      t_r = t->next->next;

      t->node_index = next_node_index++;
      t_r->node_index = next_node_index++;
      t_l->node_index = next_node_index++;

      /* pop and connect */
      t_cr = network_stack[--network_stack_top];
      t_r->back = t_cr; t_cr->back = t_r;
      t_cl = network_stack[--network_stack_top];
      t_l->back = t_cl; t_cl->back = t_l;

      /* set branch attributes */
      t->length = t_s.length;
      t->pmatrix_index = t_s.pmatrix_index;
      t_r->pmatrix_index = t_cr->pmatrix_index;
      t_r->length = t_cr->length;
      t_l->pmatrix_index = t_cl->pmatrix_index;
      t_l->length = t_cl->length;
    }
    else
    {
      t = (pll_unetwork_node_t *)calloc(1, sizeof(pll_unetwork_node_t));
      memcpy(t, &t_s, sizeof(pll_unetwork_node_t));
      assert(t->node_index < tip_count);
    }

    /* push */
    network_stack[network_stack_top++] = t;
  }

  /* root vertices must be in the stack */
  assert (network_stack_top == 2);
  assert (next_node_index == (4*tip_count - 6));

  network = network_stack[--network_stack_top];
  network->back = network_stack[--network_stack_top];
  network->back->back = network;

  if(network->pmatrix_index != network->back->pmatrix_index)
  {
    /* if pmatrix indices differ, connecting branch must be a tip */
    if(network->back->next)
    {
      pllmod_set_error(PLLMOD_NETWORK_ERROR_INVALID_NETWORK,
                       "pmatrix indices do not match in serialized network");
    }
    network->pmatrix_index = network->back->pmatrix_index;
  }

  if(network->length != network->back->length)
  {
    pllmod_set_error(PLLMOD_NETWORK_ERROR_INVALID_NETWORK,
                     "branch lengths do not matchin serialized network");
  }

  if (pll_errno)
  {
    pll_unetwork_graph_destroy(network, NULL);
    network = 0;
  }

  free(network_stack);

  return pll_unetwork_wrapnetwork(network, tip_count);
}

PLL_EXPORT int pllmod_unetwork_draw_support(pll_unetwork_t * ref_network,
                                         const double * support,
                                         pll_unetwork_node_t ** node_map,
                                         char * (*cb_serialize)(double) )
{
  if (!ref_network || !support)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Parameter is NULL!\n");
    return PLL_FAILURE;
  }

  unsigned int split_count = ref_network->edge_count - ref_network->tip_count;
  for (size_t i = 0; i < split_count; ++i)
  {
    pll_unetwork_node_t * node = node_map ? node_map[i] :
                                    ref_network->nodes[ref_network->tip_count + i];

    /* this has to be an inner node! */
    assert(node->next);

    if (node->label)
      free(node->label);

    node->label = node->next->label = node->next->next->label =
        cb_serialize ? cb_serialize(support[i]) : default_support_fmt(support[i]);
  }

  /* set support value for the root branch to both adjacent nodes */
  pll_unetwork_node_t * root = ref_network->vroot;
  if (!root->label && root->back->next)
  {
    assert(root->back->label);
    root->label = root->next->label = root->next->next->label = strdup(root->back->label);
  }
  else if (!root->back->label && root->next)
  {
    assert(root->label);
    root->back->label = root->back->next->label = root->back->next->next->label =
        strdup(root->label);
  }

  return PLL_SUCCESS;
}

/******************************************************************************/
/* Static functions */

static int rnetwork_rollback_tbr(pll_network_rollback_t * rollback_info)
{
  PLLMOD_UNUSED(rollback_info);
  assert(UNIMPLEMENTED);
  return PLL_FAILURE;
}

static int rnetwork_rollback_spr(pll_network_rollback_t * rollback_info)
{
  //TODO: Add preconditions

  pll_rnetwork_node_t * p = (pll_rnetwork_node_t *) rollback_info->SPR.prune_edge;
  pll_rnetwork_node_t * r = (pll_rnetwork_node_t *) rollback_info->SPR.regraft_edge;

  /* undo move */
  if (!pllmod_rnetwork_spr(p, r, 0, 0))
    return PLL_FAILURE;

  //TODO: set branch lengths
  //
  return PLL_SUCCESS;
}

static int rnetwork_rollback_nni(pll_network_rollback_t * rollback_info)
{
  PLLMOD_UNUSED(rollback_info);
  assert(UNIMPLEMENTED);
  return PLL_FAILURE;
}


static int unetwork_rollback_tbr(pll_network_rollback_t * rollback_info)
{
  assert(!rollback_info->rooted);
  assert(rollback_info->rearrange_type == PLLMOD_NETWORK_REARRANGE_TBR);

  pll_unetwork_node_t * p = (pll_unetwork_node_t *) rollback_info->TBR.bisect_edge;
  pll_unetwork_node_t * q = p->next->back;
  pll_unetwork_node_t * r = p->back->next->back;
  double reconn_length = rollback_info->TBR.reconn_edge.length;

  /* undo move */
    if (!pllmod_unetwork_tbr(p, &(rollback_info->TBR.reconn_edge), 0))
      return PLL_FAILURE;

  /* reset branches */
  pllmod_unetwork_set_length(p, reconn_length);
  pllmod_unetwork_set_length(q, rollback_info->TBR.bisect_left_bl);
  pllmod_unetwork_set_length(r, rollback_info->TBR.bisect_right_bl);
  pllmod_unetwork_set_length(p->next, rollback_info->TBR.reconn_parent_left_bl);
  pllmod_unetwork_set_length(p->next->next, rollback_info->TBR.reconn_parent_right_bl);
  pllmod_unetwork_set_length(p->back->next, rollback_info->TBR.reconn_child_left_bl);
  pllmod_unetwork_set_length(p->back->next->next, rollback_info->TBR.reconn_child_right_bl);

  return PLL_SUCCESS;
}

static int unetwork_rollback_spr(pll_network_rollback_t * rollback_info)
{
  assert(!rollback_info->rooted);
  assert(rollback_info->rearrange_type == PLLMOD_NETWORK_REARRANGE_SPR);

  pll_unetwork_node_t * p = (pll_unetwork_node_t *) rollback_info->SPR.prune_edge;
  pll_unetwork_node_t * r = (pll_unetwork_node_t *) rollback_info->SPR.regraft_edge;
  pll_unetwork_node_t * z1 =  p->next->back;
  pll_unetwork_node_t * z2 =  r->back;

  /* undo move */
  if (!pllmod_unetwork_spr(p, r, 0))
    return PLL_FAILURE;

  /* reset branches */
  pllmod_unetwork_set_length(z1, rollback_info->SPR.regraft_bl);
  pllmod_unetwork_set_length(p, rollback_info->SPR.prune_bl);
  pllmod_unetwork_set_length(r, rollback_info->SPR.prune_left_bl);
  pllmod_unetwork_set_length(z2, rollback_info->SPR.prune_right_bl);

  return PLL_SUCCESS;
}

static int unetwork_rollback_nni(pll_network_rollback_t * rollback_info)
{
  assert(!rollback_info->rooted);
  assert(rollback_info->rearrange_type == PLLMOD_NETWORK_REARRANGE_NNI);

  pll_unetwork_node_t * p = rollback_info->NNI.edge;
  pll_unetwork_node_t * q = p->back;

  /* undo move */
  if (!pllmod_unetwork_nni(p, rollback_info->NNI.type, 0))
      return PLL_FAILURE;

  /* reset branches */

  pllmod_unetwork_set_length(p, rollback_info->NNI.edge_bl);
  pllmod_unetwork_set_length(p->next, rollback_info->NNI.left_left_bl);
  pllmod_unetwork_set_length(p->next->next, rollback_info->NNI.left_right_bl);
  pllmod_unetwork_set_length(q->next, rollback_info->NNI.right_left_bl);
  pllmod_unetwork_set_length(q->next->next, rollback_info->NNI.right_right_bl);

  //assert(UNIMPLEMENTED);
  return PLL_SUCCESS;
}

/**
 * callback function for updating p-matrices and partials
 */
static int cb_update_matrices_partials(pll_unetwork_node_t * node, void *data)
{
  struct cb_params * st_data = (struct cb_params *) data;

  if (st_data->update_pmatrices)
  {
    unsigned int matrix_index = node->pmatrix_index;
    double branch_length = node->length;

    /* check integrity */
    assert(fabs(node->length - node->back->length) < 1e-8);
    assert(node->pmatrix_index == node->back->pmatrix_index);

    pll_update_prob_matrices (st_data->partition,
                              st_data->params_indices,
                              &matrix_index,
                              &branch_length,
                              1);
  }

  if (st_data->update_partials && !pllmod_unetwork_is_tip(node))
  {
    /* check integrity */
    assert(node->next->pmatrix_index == node->next->back->pmatrix_index);
    assert(node->next->next->pmatrix_index == node->next->next->back->pmatrix_index);

    pll_operation_t op;
    op.child1_clv_index    = node->next->back->clv_index;
    op.child1_scaler_index = node->next->back->scaler_index;
    op.child1_matrix_index = node->next->pmatrix_index;
    op.child2_clv_index    = node->next->next->back->clv_index;
    op.child2_scaler_index = node->next->next->back->scaler_index;
    op.child2_matrix_index = node->next->next->pmatrix_index;
    op.parent_clv_index    = node->clv_index;
    op.parent_scaler_index = node->scaler_index;

    pll_update_partials(st_data->partition, &op, 1);
  }

  return PLL_SUCCESS;
}

static void shuffle_network_nodes(const pll_unetwork_t * network, unsigned int seed)
{
  unsigned int node_count = network->tip_count + network->inner_tree_count + network->reticulation_count;
  pll_random_state * rstate =  pll_random_create(seed);
  pll_unetwork_node_t ** subnodes =  (pll_unetwork_node_t **) calloc(network->tip_count,
                                                     sizeof(pll_unetwork_node_t *));

  for (unsigned int i = network->tip_count; i < node_count; ++i)
  {
    pll_unetwork_node_t * node = network->nodes[i];
    unsigned int degree = 0;
    do
    {
      subnodes[degree] = node;
      degree++;
      node = node->next;
    }
    while (node != network->nodes[i]);

    // Fisher–Yates shuffle
    for (unsigned int j = degree-1; j > 0; --j)
    {
      unsigned int r = pll_random_getint(rstate, j+1);
      PLL_SWAP(subnodes[j], subnodes[r]);
    }

    // re-connect pll_unetwork_nodes in the new, shuffled order
    network->nodes[i] = node = subnodes[0];
    for (unsigned int j = 1; j < degree; ++j)
    {
      node->next = subnodes[j];
      node = node->next;
    }

    // close roundabout
    node->next = network->nodes[i];
  }

  pll_random_destroy(rstate);
  free(subnodes);
}


static void split_multi_node(pll_unetwork_t * network, pll_unetwork_node_t * first,
                             pll_unetwork_node_t * last, unsigned int degree)
{
  assert(last->next == first);
  if (degree > 3)
  {
    assert(first->next && first->next->next && first->next->next->next != first);

    // got a multifurcating node, split it in two
    unsigned int new_pmatrix_id = network->edge_count;
    unsigned int new_node_id = network->edge_count * 2;
    unsigned int new_clv_id = network->tip_count + network->inner_tree_count + network->reticulation_count;
    unsigned int new_scaler_id = network->inner_tree_count + network->reticulation_count;

    pll_unetwork_node_t * second = first->next;

    pll_unetwork_node_t * old_link = (pll_unetwork_node_t *) calloc(1, sizeof(pll_unetwork_node_t));
    pll_unetwork_node_t * new_link = (pll_unetwork_node_t *) calloc(1, sizeof(pll_unetwork_node_t));

    old_link->data = second->data;
    new_link->data = second->next->data;

    //close 'new' roundabout
    new_link->next = second->next;
    last->next = new_link;

    //close 'old' roundabout
    old_link->next = first;
    second->next = old_link;

    old_link->clv_index = second->clv_index;
    old_link->scaler_index = second->scaler_index;
    old_link->pmatrix_index = new_pmatrix_id;
    old_link->node_index = new_node_id;

    assert(new_link->next && new_link->next->next);

    new_link->pmatrix_index = new_pmatrix_id;
    new_link->node_index = new_node_id+1;

    new_link->clv_index = new_link->next->clv_index =
        new_link->next->next->clv_index = new_clv_id;
    new_link->scaler_index = new_link->next->scaler_index =
        new_link->next->next->scaler_index = new_scaler_id;

    //set backpointers old<->new
    pllmod_unetwork_connect_nodes(old_link, new_link, PLLMOD_NETWORK_DEFAULT_BRANCH_LENGTH, 1.0);

    network->nodes[network->inner_tree_count + network->reticulation_count + network->tip_count] = new_link;

    network->edge_count++;
    network->inner_tree_count++; // ?

    // split new node if needed
    split_multi_node(network, new_link, last, degree-1);
  }
}

static char * default_support_fmt(double support)
{
  char *sup_str;
  int size_alloced = asprintf(&sup_str, "%lf", support);

  return size_alloced >= 0 ? sup_str : NULL;
}

pll_unetwork_node_t * go_down_recursive(pll_unetwork_node_t * node, int * present)
{
  if (!node)
  {
	return NULL;
  }
  if (pll_unetwork_is_reticulation(node))
  {
	return go_down_recursive(pll_unetwork_get_reticulation_child(node), present);
  }
  if (present[node->clv_index])
  {
  	return node;
  }

  pll_unetwork_node_t * child1 = NULL;
  pll_unetwork_node_t * child2 = NULL;
  pll_unetwork_get_tree_children(node, &child1, &child2);

  // check left child and right child first
  if (present[child1->clv_index])
  {
	return child1;
  }
  if (present[child2->clv_index])
  {
	return child2;
  }

  if (child1->active) {
	  pll_unetwork_node_t * try1 = go_down_recursive(child1, present);
	  if (try1) {
		  return try1;
	  }
  }

  if (child2->active) {
  	  pll_unetwork_node_t * try2 = go_down_recursive(child2, present);
  	  if (try2) {
  		  return try2;
  	  }
  }

  return NULL; // we should never end up here though
}

double collect_branch_length_to_first_present_parent(pll_unetwork_node_t * node, int * present, uint64_t tree_number)
{
  if (pll_unetwork_is_reticulation(node))
  {
	return PLL_FAILURE; // this should not happen as reticulation nodes should not be in this post-order tree traversal.
  }
  double branch_sum = node->length;
  pll_unetwork_node_t * parent = pll_unetwork_get_active_parent(node);
  while (!present[parent->clv_index])
  {
	node = parent;
	branch_sum += node->length;
	parent = pll_unetwork_get_active_parent(node);
  }
  return branch_sum;
}

int cb_full_unetwork_traversal(pll_unetwork_node_t * node) {
	(void) node;
	return 1;
}

PLL_EXPORT int pllmod_unetwork_tree_buildarrays(pll_unetwork_t * network, uint64_t tree_number, pll_displayed_tree_t * result, unsigned int fake_clv_index, unsigned int fake_pmatrix_index) { // TODO: FIXME: replace path collapsing by using the fake node...
	unsigned int nodes_count = network->reticulation_count + network->inner_tree_count + network->tip_count;
	unsigned int inner_nodes_count = network->reticulation_count + network->inner_tree_count;
	unsigned int branch_count = network->inner_tree_count * 2 + network->reticulation_count;

	pll_unetwork_node_t ** trav_buffer = (pll_unetwork_node_t **) malloc(nodes_count * sizeof(pll_unetwork_node_t *));
	unsigned int trav_size;
	if (!pll_unetwork_tree_traverse(network, PLL_TREE_TRAVERSE_POSTORDER, cb_full_unetwork_traversal, trav_buffer, &trav_size, tree_number)) {
		return PLL_FAILURE;
	}

	unsigned int i;

    result->branch_lengths = (double *) malloc(branch_count * sizeof(double));
    result->operations = (pll_operation_t *) malloc(inner_nodes_count * sizeof(pll_operation_t));
    result->ops_count = 0;
    result->pmatrix_indices = (unsigned int *) malloc(branch_count * sizeof(unsigned int));
    result->matrix_count = 0;

    // TODO: Fill the operations array, note that we have to sum up the branch lengths that lie on a path...
    // (by going up through the parents, this gets a bit tricky when encountering reticulations because then we also have to figure out whether the edge belongs to the current tree)
    for (i = 0; i < trav_size; ++i)
    {
      pll_unetwork_node_t * node = trav_buffer[i];

      /* do not store the branch of the root, since it does not exist */
      if (i < trav_size-1)
      {
    	(result->branch_lengths)[result->matrix_count] = node->length;
        (result->pmatrix_indices)[result->matrix_count] = node->pmatrix_index;
        result->matrix_count = result->matrix_count + 1;
      }

      if (!pll_unetwork_is_leaf(node)) // inner tree node
      {
        result->operations[result->ops_count].parent_clv_index = node->clv_index;
        result->operations[result->ops_count].parent_scaler_index = node->scaler_index;

        // These values change! It could be that the child is not in the trav_buffer, in this case go down until a child has been found!
        // Keep in mind that only nodes with at most one non-dead child are thrown out from the trav_buffer, and dead nodes are thrown out, too

        pll_unetwork_node_t * child1 = NULL;
        pll_unetwork_node_t * child2 = NULL;
        if (pll_unetwork_is_inner_tree(node)) {
        	pll_unetwork_get_tree_children(node, &child1, &child2);
        }
        else {
        	child1 = pll_unetwork_get_reticulation_child(node);
        }

        if (child1 && child1->active) {
          result->operations[result->ops_count].child1_clv_index = child1->clv_index;
          result->operations[result->ops_count].child1_scaler_index = child1->scaler_index;
          result->operations[result->ops_count].child1_matrix_index = child1->pmatrix_index;
        } else {
          result->operations[result->ops_count].child1_clv_index = fake_clv_index;
          result->operations[result->ops_count].child1_scaler_index = -1;
          result->operations[result->ops_count].child1_matrix_index = fake_pmatrix_index;
        }

        if (child2 && child2->active) {
          result->operations[result->ops_count].child2_clv_index = child2->clv_index;
          result->operations[result->ops_count].child2_scaler_index = child2->scaler_index;
          result->operations[result->ops_count].child2_matrix_index = child2->pmatrix_index;
        } else {
          result->operations[result->ops_count].child2_clv_index = fake_clv_index;
          result->operations[result->ops_count].child2_scaler_index = -1;
          result->operations[result->ops_count].child2_matrix_index = fake_pmatrix_index;
        }

        result->ops_count = result->ops_count + 1;
      }
    }
    //free(present);
	return PLL_SUCCESS;
}
