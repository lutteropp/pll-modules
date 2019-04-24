/*
 Copyright (C) 2016 Diego Darriba

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
#include "pll_network.h"

#include "../pllmod_common.h"

/**
 * @file rnetwork_operations.c
 *
 * @brief Operations on rooted network structures
 *
 * @author Diego Darriba
 */

/* finds the self and sister pointers */
PLL_EXPORT int pllmod_rnetwork_get_sibling_pointers(pll_rnetwork_node_t * node,
                                                 pll_rnetwork_node_t ***self,
                                                 pll_rnetwork_node_t ***sister)
{
  if (!node->parent)
  {
    /* if there is no parent, there are no pointers */
    if (self) *self = NULL;
    if (sister) *sister = NULL;
  }
  else if (node->parent->left == node)
  {
    if (self) *self = &(node->parent->left);
    if (sister) *sister = &(node->parent->right);
  }
  else if (node->parent->right == node)
  {
    if (self) *self = &(node->parent->right);
    if (sister) *sister = &(node->parent->left);
  }
  else
  {
    /* `node` is not the left nor the right child of its parent */
    if (self) *self = NULL;
    if (sister) *sister = NULL;
    pllmod_set_error(PLLMOD_TREE_ERROR_INVALID_TREE,
                     "Tree is not consistent");
    return PLL_FAILURE;
  }

  return PLL_SUCCESS;
}

/**
 * @brief Prunes a subnetwork in a rooted network
 * @param node the node to prune
 * @return the new connected node, if the operation was applied correctly, NULL otherwise
 */
PLL_EXPORT pll_rnetwork_node_t * pllmod_rnetwork_prune(pll_rnetwork_node_t * node)
{
  pll_rnetwork_node_t **self_ptr, **sister_ptr, **parent_ptr;
  pll_rnetwork_node_t *connected_node = NULL;
  assert(node);

  if (!node->parent)
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_SPR_INVALID_NODE,
                     "Attempting to prune the root node");
    return NULL;
  }

  if (!pllmod_rnetwork_get_sibling_pointers(node,
                               &self_ptr,
                               &sister_ptr))
  {
    /* return and spread error */
    return NULL;
  }
  else
  {
    assert (self_ptr && sister_ptr);
  }

  /* connect adjacent subnetworks together */
  if (node->parent->parent)
  {
    /* connect parent->parent and sister */
    connected_node = node->parent->parent;
    if (!pllmod_rnetwork_get_sibling_pointers(node->parent,
                                 &parent_ptr,
                                 NULL))
    {
      /* return and spread error */
      return NULL;
    }
    else
    {
      assert (parent_ptr);
    }
    *parent_ptr = *sister_ptr;
    (*sister_ptr)->parent = node->parent->parent;

    /* disconnect pruned network */
    *sister_ptr = NULL;
    node->parent->parent = NULL;
  }
  else
  {
    /* re-root */
    connected_node = *sister_ptr;

    (*sister_ptr)->parent = NULL;

    /* disconnect pruned network */
    *sister_ptr = NULL;
  }

  return connected_node;
}

/**
 * @brief Regrafts a dettached subnetwork into a branch
 *
 * @param node the node to regraft
 * @param network the target branch
 *
 * @return PLL_SUCCESS if the operation was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_rnetwork_regraft(pll_rnetwork_node_t * node,
                                    pll_rnetwork_node_t * network)
{
  pll_rnetwork_node_t *parent_node;
  pll_rnetwork_node_t **edge_from_parent = 0,
              **edge_to_child    = 0;

  /* node must contain a dettached parent */
  if (!node->parent || node->parent->parent)
  {
    pllmod_set_error(PLLMOD_TREE_ERROR_SPR_INVALID_NODE,
                     "Attempting to regraft a node without dettached parent");
    return PLL_FAILURE;
  }

  /* `node` parent should contain a pointer to `node` */
  if (!pllmod_rnetwork_get_sibling_pointers(node,
                               NULL,
                               &edge_to_child))
  {
    /* return and spread error */
    assert(pll_errno);
    return PLL_FAILURE;
  }

  parent_node = node->parent;

  if (network->parent)
  {
    if (network->parent && !pllmod_rnetwork_get_sibling_pointers(network,
                                  &edge_from_parent,
                                  NULL))
    {
      /* return and spread error */
      assert(pll_errno);
      return PLL_FAILURE;
    }
    assert(edge_from_parent);

    /* set new parents */
    parent_node->parent = network->parent;

    /* set new children */
    *edge_from_parent = parent_node;
  }

  network->parent = parent_node;
  *edge_to_child = network;

  return PLL_SUCCESS;
}

/**
 * Performs one SPR move
 * The CLV, scaler and pmatrix indices are updated.
 *
 * @param[in] p_node Edge to be pruned
 * @param[in] r_network Edge to be regrafted
 * @param[in] root The network root (it might change)
 * @param[in,out] rollback_info Rollback information
 * @return PLL_SUCCESS if the move was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_rnetwork_spr(pll_rnetwork_node_t * p_node,
                                pll_rnetwork_node_t * r_network,
                                pll_rnetwork_node_t ** root,
                                pll_network_rollback_t * rollback_info)
{
  pll_rnetwork_node_t **self_ptr, **sister_ptr;

  if (!pllmod_rnetwork_get_sibling_pointers(p_node,
                               &self_ptr,
                               &sister_ptr))
  {
    /* return and spread error */
    assert(pll_errno);
    return PLL_FAILURE;
  }
  else
  {
    assert (self_ptr && sister_ptr);
  }

  /* save rollback information */
  if (rollback_info)
  {
    rollback_info->rearrange_type     = PLLMOD_TREE_REARRANGE_SPR;
    rollback_info->rooted             = 1;
    rollback_info->SPR.prune_edge     = (void *) p_node;
    rollback_info->SPR.regraft_edge   = (void *) *sister_ptr;
    //TODO: Set branch lengths
    // rollback_info->SPR.prune_bl       = p_edge->parent->length;
    // rollback_info->SPR.prune_left_bl  = (*sister)->length;
    // rollback_info->SPR.prune_right_bl = p_edge->->length;
    // rollback_info->SPR.regraft_bl     = r_network->length;
  }

  if (pllmod_rnetwork_prune(p_node) == NULL)
  {
    /* return and spread error */
    assert(pll_errno);
    return PLL_FAILURE;
  }

  if (!pllmod_rnetwork_regraft(p_node,
                            r_network))
  {
    /* return and spread error */
    assert(pll_errno);
    return PLL_FAILURE;
  }

  /* reset root in case it has changed */
  while (root && (*root)->parent) *root = (*root)->parent;

  return PLL_SUCCESS;
}

/* re-roots the network at the branch connecting `new_root` with its parent */
PLL_EXPORT int pllmod_rnetwork_reroot(pll_rnetwork_node_t * root,
                                   pll_rnetwork_node_t * new_root)
{
  return PLL_FAILURE;
}

static void rnetwork_nodes_at_node_dist_down(pll_rnetwork_node_t * root,
                                          pll_rnetwork_node_t ** outbuffer,
                                          unsigned int * node_count,
                                          int min_distance,
                                          int max_distance)
{
  if (max_distance < 0) return;

  if (min_distance < 0)
  {
    outbuffer[*node_count] = root;
    *node_count = *node_count + 1;
  }

  if (!(root->left && root->right)) return;

  rnetwork_nodes_at_node_dist_down(root->left,
                                outbuffer,
                                node_count,
                                min_distance-1,
                                max_distance-1);
  rnetwork_nodes_at_node_dist_down(root->right,
                                outbuffer,
                                node_count,
                                min_distance-1,
                                max_distance-1);
}

PLL_EXPORT int pllmod_rnetwork_nodes_at_node_dist(pll_rnetwork_node_t * root,
                                               pll_rnetwork_node_t ** outbuffer,
                                               unsigned int * node_count,
                                               int min_distance,
                                               int max_distance)
{
  pll_rnetwork_node_t * current_root = root;
  pll_rnetwork_node_t ** sister_ptr;

  if (max_distance < min_distance)
    {
      pllmod_set_error(PLLMOD_ERROR_INVALID_RANGE,
                 "Invalid distance range: %d..%d (max_distance < min_distance)",
                 min_distance, max_distance);
      return PLL_FAILURE;
    }

  *node_count = 0;

  while (current_root->parent)
  {
    if (!pllmod_rnetwork_get_sibling_pointers(current_root,
                                 NULL,
                                 &sister_ptr))
    {
      /* return and spread error */
      assert(pll_errno);
      return PLL_FAILURE;
    }

    --min_distance;
    --max_distance;

    current_root = current_root->parent;
    if (min_distance < 0)
    {
      outbuffer[*node_count] = current_root;
      *node_count = *node_count + 1;
    }
    rnetwork_nodes_at_node_dist_down(*sister_ptr,
                                  outbuffer,
                                  node_count,
                                  min_distance-1,
                                  max_distance-1);
  }


  return PLL_SUCCESS;
}

PLL_EXPORT int pllmod_rnetwork_connect_nodes(pll_rnetwork_node_t * parent, pll_rnetwork_node_t * child, double length, double prob) {
	if (!(parent && child))
		return PLL_FAILURE;

	if (parent->is_reticulation) {
		parent->child = child;
		if (child->is_reticulation) {
			if (!child->first_parent) {
				child->first_parent = parent;
			} else if (!child->second_parent) {
				child->second_parent = parent;
			} else {
				return PLL_FAILURE;
			}
		} else {
			child->parent = parent;
		}
		child->prob = prob;
	} else {
		if (!parent->left) {
			parent->left = child;
			if (child->is_reticulation) {
				if (!child->first_parent) {
					child->first_parent = parent;
				} else if (!child->second_parent) {
					child->second_parent = parent;
				} else {
					return PLL_FAILURE;
				}
			} else {
				child->parent = parent;
			}
			child->prob = prob;
		} else if (!parent->right) {
			parent->right = child;
			if (child->is_reticulation) {
				if (!child->first_parent) {
					child->first_parent = parent;
				} else if (!child->second_parent) {
					child->second_parent = parent;
				} else {
					return PLL_FAILURE;
				}
			} else {
				child->parent = parent;
			}
			child->prob = prob;
		} else {
			return PLL_FAILURE;
		}
	}

	/* PMatrix index is set to parent node */
	child->pmatrix_index = parent->pmatrix_index;

	return PLL_SUCCESS;
}
