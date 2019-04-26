/*
 Copyright (C) 2016 Diego Darriba, modified to networks in 2019 by Sarah Lutteropp

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
  * @file unetwork_operations.c
  *
  * @brief Operations on unrooted network structures
  *
  * @author Diego Darriba, modifications by Sarah Lutteropp
  */

#include "pll_network.h"

#include "../pllmod_common.h"

static void unetwork_nodes_at_dist(pll_unetwork_node_t * node,
                                pll_unetwork_node_t ** outbuffer,
                                unsigned int * index,
                                unsigned int min_distance,
                                unsigned int max_distance,
                                unsigned int depth);



/******************************************************************************/
/* Topological operations */

/**
 * @brief Bisects the network by removing one edge
 *
 * Removes the edge \p edge and frees the nodes defining that edge.
 * Reconnects the subnetworks at the sides of the edge (figure below).
 * The branch lengths of the new edges are the sum of the removed ones.
 * The join branch contains the pmatrix index of the parent edges
 * The removed pmatrix indices are returned in the field
 *     'additional_pmatrix_index' of both output subnetworks
 *
 * Returns the new parent and child edges, where parent is the closest to \p edge.
 *
 *   A            C              A        C
 *    \___edge___/       ---->   |        |
 *    /          \               |        |
 *   B            D              B        D
 *   A,B,C,D are subnetworks
 *
 * @param[in] edge            edge to remove
 * @param[out] parent_subnetwork edge corresponding to the 'edge' subnetwork
 * @param[out] child_subnetwork  edge corresponding to the 'edge->back' subnetwork
 * @return PLL_SUCCESS if OK
 */
PLL_EXPORT int pllmod_unetwork_bisect(pll_unetwork_node_t * edge,
                                   pll_unetwork_node_t ** parent_subnetwork,
                                   pll_unetwork_node_t ** child_subnetwork)
{
  assert(parent_subnetwork);
  assert(child_subnetwork);

  pll_unetwork_node_t * aux_network;

  if (!edge->next)
    return PLL_FAILURE;

  pll_unetwork_node_t * c_edge = edge->back;

  /* connect parent subnetwork */
  (*parent_subnetwork) = edge->next->back;
  aux_network = edge->next->next->back;

  pllmod_unetwork_connect_nodes(*parent_subnetwork,
                             aux_network,
                             (*parent_subnetwork)->length + aux_network->length, 1.0);

  edge->next->pmatrix_index = edge->next->next->pmatrix_index;

  /* connect child subnetwork */
  (*child_subnetwork) = c_edge->next->back;
  aux_network = c_edge->next->next->back;

  pllmod_unetwork_connect_nodes(*child_subnetwork,
                             aux_network,
                             (*child_subnetwork)->length + aux_network->length, 1.0);

  c_edge->next->pmatrix_index = c_edge->next->next->pmatrix_index;

  return PLL_SUCCESS;
}

/**
 * Reconnects two subnetworks by adding 2 new nodes and 1 edge.
 *
 * Adds 1 new edge connecting edges \p edge.parent and \p edge.child with
 * length \p edge.length.
 *
 *   A       C         A              C
 *   |       |  ---->   \            /
 *                       e1--edge--e2
 *   |       |          /            \
 *   B       D         B              D
 *   A,B,C,D are subnetworks
 *
 * @param edge                 new edge (edge structure)
 * @param pruned_edge          edge to prune, defined by a network node
 *
 * @return the new created edge
 */
PLL_EXPORT pll_network_edge_t pllmod_unetwork_reconnect(pll_network_edge_t * edge,
                                                  pll_unetwork_node_t * pruned_edge)
{
  /* create and connect 2 new nodes */
  pll_unetwork_node_t *parent_node, *child_node;
  assert(pruned_edge->back);

  parent_node = pruned_edge;
  child_node  = pruned_edge->back;
  assert(parent_node->back == child_node && child_node->back == parent_node);

  assert(!pllmod_unetwork_is_tip(parent_node));
  assert(!pllmod_unetwork_is_tip(child_node));

  pll_network_edge_t new_edge;
  new_edge.edge.unetwork.child = child_node;
  new_edge.length = edge->length;

  /* set length */
  pllmod_unetwork_set_length(parent_node, edge->length);

  /* reconnect parent close to edge.parent */
  pllmod_unetwork_connect_nodes(parent_node->next->next,
                             edge->edge.unetwork.parent->back,
                             edge->edge.unetwork.parent->back->length, edge->edge.unetwork.parent->back->prob);

  pllmod_unetwork_connect_nodes(edge->edge.unetwork.parent,
                             parent_node->next,
                             0, 1.0);

  /* reconnect child close to edge.child */
  pllmod_unetwork_connect_nodes(child_node->next->next,
                             edge->edge.unetwork.child->back,
                             edge->edge.unetwork.child->back->length, edge->edge.unetwork.child->back->prob);

  pllmod_unetwork_connect_nodes(edge->edge.unetwork.child,
                             child_node->next,
                             0, 1.0);

  return new_edge;
}

/**
 * @brief Prunes a subnetwork in an unrooted network
 *
 * Disconnecs an edge (e1) and connects the adjacent nodes. New branch (A-B)
 * length is set to the sum of lengths of previous branch (e1-A + e1-B)
 *
 *   A              C              A                   C
 *    \            /               |                  /
 *     e1--edge--e2        --->    |  +   e1--edge--e2
 *    /            \               |                  \
 *   B              D              B                   D
 *   A,B,C,D are subnetworks
 *
 *  Note that `edge` is disconnected after the operation
 *
 * @param edge the edge to prune
 * @return the new connected edge, if the operation was applied correctly
 */
PLL_EXPORT pll_unetwork_node_t * pllmod_unetwork_prune(pll_unetwork_node_t * edge)
{
  pll_unetwork_node_t *edge1, *edge2;

  assert(edge);
  if (!edge->next)
  {
    /* invalid node */
    pllmod_set_error(PLLMOD_NETWORK_ERROR_SPR_INVALID_NODE,
                     "Attempting to prune a tip node");
    return NULL;
  }

  /* connect adjacent subnetworks together */
  edge1 = edge->next->back;
  edge2 = edge->next->next->back;
  pllmod_unetwork_connect_nodes(edge1, edge2, edge1->length + edge2->length, 1.0);

  /* disconnect pruned edge */
  edge->next->back = edge->next->next->back = NULL;

  return edge1;
}

/**
 * @brief Regrafts an edge into a network
 *
 * Connects a disconnected edge (provided by `e2` in the graph below)
 * into a network
 *
 *  A                    C         A              C
 *   \                   |          \            /
 *    e1--edge--e2   +   |   --->    e1--edge--e2
 *   /                   |          /            \
 *  B                    D         B              D
 *   A,B,C,D are subnetworks
 *
 *  The length of the new branches (e2-C and e2-D) are set to half the length
 *  of the removed branch (C-D)
 *
 * @param edge the edge to regraft
 * @param network the network to connect `edge` to
 * @return PLL_SUCCESS if the operation was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_unetwork_regraft(pll_unetwork_node_t * edge,
                                    pll_unetwork_node_t * network)
{
  pll_unetwork_node_t *edge1, *edge2;
  double new_length;

  assert(edge && network);
  if (!edge->next)
  {
    /* invalid node */
    pllmod_set_error(PLLMOD_NETWORK_ERROR_SPR_INVALID_NODE,
                     "Attempting to regraft a tip node");
    return PLL_FAILURE;
  }
  if (edge->next->back || edge->next->next->back)
  {
    /* invalid node */
    pllmod_set_error(PLLMOD_NETWORK_ERROR_SPR_INVALID_NODE,
                     "Attempting to regraft a connected node");
    return PLL_FAILURE;
  }

  /* connect network with edge, splitting the branch designed by network */
  edge1      = network;
  edge2      = network->back;
  new_length = network->length/2;
  pllmod_unetwork_connect_nodes(edge1, edge->next,       new_length, 1.0);
  pllmod_unetwork_connect_nodes(edge->next->next, edge2, new_length, 1.0);

  return PLL_SUCCESS;
}

/**
 * @brief Interchanges 2 edges, represented by 2 internal nodes
 *
 * CLV and scaler indices, and labels are interchanged between nodes to match
 * the other 2 nodes in the triplet.
 *
 * @return PLL_SUCCESS if the operation was applied correctly,
 *         PLL_FAILURE otherwise (check pll_errmsg for details)
 */
PLL_EXPORT int pllmod_unetwork_interchange(pll_unetwork_node_t * node1,
                                        pll_unetwork_node_t * node2)
{
  pll_unetwork_node_t *next1 = node2->back;
  pll_unetwork_node_t *next2 = node1->back;

  pllmod_unetwork_connect_nodes(next1, node1, next1->length, next1->prob);
  pllmod_unetwork_connect_nodes(next2, node2, next2->length, next2->prob);

  return PLL_SUCCESS;
}

/**
 * @brief Creates a new circular node
 *
 *           n2
 *          / |
 *        n1  |
 *          \ |
 *           n3
 *
 * All parameters are shared among the nodes in the triplet
 *
 * @param clv_index    the clv_index
 * @param scaler_index the scaler index
 * @param label        the node label
 * @param data         the data pointer
 *
 * @return the new node
 */
PLL_EXPORT pll_unetwork_node_t * pllmod_unetwork_create_node(unsigned int clv_index,
                                                  int scaler_index,
                                                  char * label,
                                                  void * data)
{
  pll_unetwork_node_t * new_node = (pll_unetwork_node_t *)calloc(1, sizeof(pll_unetwork_node_t));
  new_node->next         = (pll_unetwork_node_t *)calloc(1, sizeof(pll_unetwork_node_t));
  new_node->next->next   = (pll_unetwork_node_t *)calloc(1, sizeof(pll_unetwork_node_t));
  if (!(new_node && new_node->next && new_node->next->next))
  {
    pllmod_set_error(PLL_ERROR_MEM_ALLOC,
                     "Cannot allocate memory for new node\n");
    return NULL;
  }

  new_node->next->next->next = new_node;
  new_node->label = label;
  new_node->next->label =
    new_node->next->next->label =
    new_node->label;
  new_node->next->data =
    new_node->next->next->data =
    new_node->data = data;
  new_node->next->length =
    new_node->next->next->length =
    new_node->length = 0;
  new_node->next->clv_index =
    new_node->next->next->clv_index =
    new_node->clv_index = clv_index;
  new_node->next->scaler_index =
    new_node->next->next->scaler_index =
    new_node->scaler_index = scaler_index;
  new_node->back =
    new_node->next->back =
    new_node->next->next->back = NULL;
  return new_node;
}

/**
 * @brief Connects 2 nodes and sets the pmatrix index and branch length
 *
 * connects `back` pointers of `child` and `parent`
 * pmatrix index for `child` is set to the one in `parent`
 *
 * @param[in,out] parent the parent node
 * @param[in,out] child  the child node
 * @param[in] length     the branch length
 *
 */
PLL_EXPORT int pllmod_unetwork_connect_nodes(pll_unetwork_node_t * parent,
                                          pll_unetwork_node_t * child,
                                          double length, double prob)
{
  if(!(parent && child))
    return PLL_FAILURE;

  parent->back = child;
  child->back = parent;
  pllmod_unetwork_set_length(parent, length);
  parent->prob = parent->back->prob = prob;
  parent->incoming = 0;
  child->incoming = 1;

  /* PMatrix index is set to parent node */
  child->pmatrix_index = parent->pmatrix_index;

  return PLL_SUCCESS;
}

/******************************************************************************/
/* Topological search */

/**
 * Returns the list of nodes at a distance between \p min_distance and
 * \p max_distance from a specified node
 *
 * @param[in] node the root node
 * @param[out] outbuffer the list of nodes. Outbuffer should be allocated
 * @param[out] node_count the number of nodes returned in \p outbuffer
 * @param[in] min_distance the minimum distance to check
 * @param[in] max_distance the maximum distance to check
 */
PLL_EXPORT int pllmod_unetwork_nodes_at_node_dist(pll_unetwork_node_t * node,
                                               pll_unetwork_node_t ** outbuffer,
                                               unsigned int * node_count,
                                               unsigned int min_distance,
                                               unsigned int max_distance)
{
  if (!node->next)
  {
    pllmod_set_error(PLLMOD_ERROR_INVALID_NODE_TYPE,
                     "Internal node expected, but tip node was provided");
    return PLL_FAILURE;
  }

  if (max_distance < min_distance)
    {
      pllmod_set_error(PLLMOD_ERROR_INVALID_RANGE,
                 "Invalid distance range: %d..%d (max_distance < min_distance)",
                 min_distance, max_distance);
      return PLL_FAILURE;
    }

  *node_count = 0;

  /* we will traverse an unrooted network in the following way

               1
             /
          --*
             \
               2
    */


  unetwork_nodes_at_dist(node, outbuffer, node_count, min_distance, max_distance, 0);

  return PLL_SUCCESS;
}

/**
 * Returns the list of nodes at a distance between \p min_distance and
 * \p max_distance from a specified edge
 *
 * @param[in] edge the root edge
 * @param[out] outbuffer the list of nodes. Outbuffer should be allocated
 * @param[out] node_count the number of nodes returned in \p outbuffer
 * @param[in] min_distance the minimum distance to check
 * @param[in] max_distance the maximum distance to check
 */

PLL_EXPORT int pllmod_unetwork_nodes_at_edge_dist(pll_unetwork_node_t * edge,
                                               pll_unetwork_node_t ** outbuffer,
                                               unsigned int * node_count,
                                               unsigned int min_distance,
                                               unsigned int max_distance)
{
  unsigned int depth = 0;

  if (!edge->next)
  {
    pllmod_set_error(PLLMOD_ERROR_INVALID_NODE_TYPE,
                     "Internal node expected, but tip node was provided");
    return PLL_FAILURE;
  }

  if (max_distance < min_distance)
    {
      pllmod_set_error(PLLMOD_ERROR_INVALID_RANGE,
                 "Invalid distance range: %d..%d (max_distance < min_distance)",
                 min_distance, max_distance);
      return PLL_FAILURE;
    }

  *node_count = 0;

  /* we will traverse an unrooted network in the following way

       3          1
        \        /
         * ---- *
        /        \
       4          2
   */

  unetwork_nodes_at_dist(edge->back, outbuffer, node_count,
                      min_distance, max_distance, depth+1);
  unetwork_nodes_at_dist(edge, outbuffer, node_count,
                      min_distance, max_distance, depth);

  return PLL_SUCCESS;
}


/******************************************************************************/
/* static functions */

static void unetwork_nodes_at_dist(pll_unetwork_node_t * node,
                                pll_unetwork_node_t ** outbuffer,
                                unsigned int * index,
                                unsigned int min_distance,
                                unsigned int max_distance,
                                unsigned int depth)
{
  if (depth >= min_distance && depth <= max_distance)
  {
    outbuffer[*index] = node;
    *index = *index + 1;
  }

  if (depth >= max_distance || !(node->next)) return;

  unetwork_nodes_at_dist(node->next->back, outbuffer, index,
                      min_distance, max_distance, depth+1);
  unetwork_nodes_at_dist(node->next->next->back, outbuffer, index,
                      min_distance, max_distance, depth+1);
}
