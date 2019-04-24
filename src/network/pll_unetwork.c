#include "pll_network.h"

#include "../pllmod_common.h"

int cb_full_network_traversal(pll_unetwork_node_t * node) {
	(void) node;
	return 1;
}

PLL_EXPORT pll_unetwork_t * pllmod_unetwork_create_random(unsigned int taxa_count,
                                                    const char * const* names,
                                                    unsigned int random_seed) {
	return PLL_FAILURE;
}

PLL_EXPORT int pllmod_unetwork_extend_random(pll_unetwork_t * network,
                                          unsigned int ext_taxa_count,
                                          const char * const* ext_names,
                                          unsigned int random_seed) {
	return PLL_FAILURE;
}

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
                                            unsigned int * score) {
	return PLL_FAILURE;
}

pll_unetwork_t * pllmod_unetwork_create_parsimony_multipart(unsigned int taxon_count,
                                                      char * const * taxon_names,
                                                      unsigned int partition_count,
                                                      pll_partition_t * const * partitions,
                                                      unsigned int random_seed,
                                                      unsigned int * score) {
	return PLL_FAILURE;
}

PLL_EXPORT pll_unetwork_t * pllmod_unetwork_resolve_multi(const pll_unetwork_t * multi_network,
                                                    unsigned int random_seed,
                                                    int * clv_index_map) {
	return PLL_FAILURE;
}

PLL_EXPORT int pllmod_unetwork_is_tip(const pll_unetwork_node_t * node)
{
  return node_is_leaf(node);
}

PLL_EXPORT void pllmod_unetwork_set_length(pll_unetwork_node_t * edge,
                                            double length)
{
  edge->length = edge->back->length = length;
}

PLL_EXPORT void pllmod_unetwork_set_prob(pll_unetwork_node_t * edge,
                                            double prob)
{
  edge->prob = edge->back->prob = prob;
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
  }
}

PLL_EXPORT int pllmod_unetwork_outgroup_root(pll_unetwork_t * network,
                                          unsigned int * outgroup_tip_ids,
                                          unsigned int outgroup_size,
                                          int add_root_node)
{
  return PLL_FAILURE;
}

static int unetwork_traverse_apply(pll_unetwork_node_t * node,
                                int (*cb_pre_trav)(pll_unetwork_node_t *, void *),
                                int (*cb_in_trav)(pll_unetwork_node_t *, void *),
                                int (*cb_post_trav)(pll_unetwork_node_t *, void *),
                                void *data)
{
  int retval = 1;

  if (!node->active)
	return retval;

  pll_unetwork_node_t * child_tree = 0;

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

  child_tree = node->next;
  while(child_tree != node)
  {
	if (child_tree->active) {
      retval &= unetwork_traverse_apply(child_tree->back,
                                     cb_pre_trav, cb_in_trav, cb_post_trav, data);

      if (cb_in_trav &&
          child_tree->next != node &&
          !cb_in_trav(child_tree, data))
        return PLL_FAILURE;
	}

    child_tree = child_tree->next;
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
