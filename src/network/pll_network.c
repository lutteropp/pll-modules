#include "pll_network.h"
#include "../pllmod_common.h"

int cb_full_traversal(pll_rnetwork_node_t * node) {
	(void) node;
	return 1;
}

pll_rnetwork_node_t * find_present_left_child(pll_rnetwork_node_t * node, int * present, uint64_t tree_number)
{
  pll_rnetwork_node_t* left = node->left;
  while (!present[left->idx])
  {
    // go down further
	if (left->is_reticulation)
	{
	  left = left->child;
	}
	else
	{ // inner tree node
	  // TODO: we need to go further down... but which child to take?
	}
  }
  return left;
}

pll_rnetwork_node_t * find_present_right_child(pll_rnetwork_node_t * node, int * present, uint64_t tree_number)
{
  pll_rnetwork_node_t * right = node->right;
  while (!present[right->idx])
  {
    // go down further
    if (right->is_reticulation)
	{
	  right = right->child;
	}
	else
	{ // inner tree node
	  // TODO: we need to go further down... but which child to take?
	}
  }
  return right;
}

double collect_branch_length_to_first_present_parent(pll_rnetwork_node_t * node, int * present)
{
  double branch_sum = node->length;
  while (!present[node->parent->idx])
  {
	node = node->parent;
	branch_sum += node->length;
  }
  return branch_sum;
}

PLL_EXPORT int pll_rnetwork_tree_buildarrays(pll_rnetwork_t * network, uint64_t tree_number, pll_displayed_tree_t * result) {
	unsigned int nodes_count = network->reticulation_count + network->inner_tree_count + network->tip_count;
	unsigned int inner_nodes_count = network->reticulation_count + network->inner_tree_count;
	unsigned int branch_count = network->inner_tree_count * 2 + network->reticulation_count;

	pll_rnetwork_node_t ** trav_buffer = (pll_rnetwork_node_t **) malloc(nodes_count * sizeof(pll_rnetwork_node_t *));
	unsigned int trav_size;
	if (!pll_rnetwork_tree_traverse(network, PLL_TREE_TRAVERSE_POSTORDER, cb_full_traversal, trav_buffer, &trav_size, tree_number)) {
		return PLL_FAILURE;
	}

	int * present = (int *) calloc(nodes_count, sizeof(int));
	unsigned int i;
	for (i = 0; i < trav_size; ++i)
	{
	  present[trav_buffer[i]->idx] = 1;
	}

    result->branch_lengths = (double *) malloc(branch_count * sizeof(double));
    result->operations = (pll_operation_t *) malloc(inner_nodes_count * sizeof(pll_operation_t));
    result->ops_count = 0;
    result->pmatrix_indices = (unsigned int *) malloc(branch_count * sizeof(unsigned int));
    result->matrix_count = 0;

    // TODO: Fill the operations array, note that we have to sum up the branch lengths that lie on a path...
    // (by going up through the parents, this gets a bit tricky when encountering reticulations because then we also have to figure out whether the edge belongs to the current tree)
    for (i = 0; i < trav_size; ++i)
    {
      pll_rnetwork_node_t * node = trav_buffer[i];

      if (node->is_reticulation)
      {
        free(present);
        return PLL_FAILURE; // because the reticulations should have been thrown out in the post-order tree-traversal already...
      }

      /* do not store the branch of the root, since it does not exist */
      if (i < trav_size-1)
      {
        *(result->branch_lengths)++ = collect_branch_length_to_first_present_parent(node, present);
        *(result->pmatrix_indices)++ = node->idx;
        result->matrix_count = result->matrix_count + 1;
      }

      if (node->left) // inner tree node
      {
        result->operations[result->ops_count].parent_clv_index = node->idx;
        result->operations[result->ops_count].parent_scaler_index = node->scaler_idx;

        // These values change! It could be that the child is not in the trav_buffer, in this case go down until a child has been found!
        // Keep in mind that only nodes with at most one non-dead child are thrown out from the trav_buffer, and dead nodes are thrown out, too

        pll_rnetwork_node_t * left = find_present_left_child(node, present, tree_number);
        pll_rnetwork_node_t * right = find_present_right_child(node, present, tree_number);

        result->operations[result->ops_count].child1_clv_index = left->idx;
        result->operations[result->ops_count].child1_scaler_index = left->scaler_idx;
        result->operations[result->ops_count].child1_matrix_index = left->idx;

        result->operations[result->ops_count].child2_clv_index = right->idx;
        result->operations[result->ops_count].child2_scaler_index = right->scaler_idx;
        result->operations[result->ops_count].child2_matrix_index = right->idx;

        result->ops_count = result->ops_count + 1;
      }
    }
    free(present);
	return PLL_SUCCESS;
}
