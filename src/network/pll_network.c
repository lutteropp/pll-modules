#include "pll_network.h"
#include "../pllmod_common.h"

int cb_full_traversal(pll_rnetwork_node_t * node) {
	(void) node;
	return 1;
}

PLL_EXPORT int pll_rnetwork_tree_buildarrays(pll_rnetwork_t * network, uint64_t tree_number, pll_displayed_tree_t * result) {
	unsigned int nodes_count = network->reticulation_count + network->inner_tree_count + network->tip_count;
	unsigned int inner_nodes_count = network->reticulation_count + network->inner_tree_count;
	unsigned int branch_count = network->inner_tree_count * 2 + network->reticulation_count;

	pll_rnetwork_node_t ** travbuffer = (pll_rnetwork_node_t **) malloc(nodes_count * sizeof(pll_rnetwork_node_t *));
	unsigned int trav_size;
	if (!pll_rnetwork_tree_traverse(network, PLL_TREE_TRAVERSE_POSTORDER, cb_full_traversal, travbuffer, &trav_size, tree_number)) {
		return PLL_FAILURE;
	}

    result->branch_lengths = (double *) malloc(branch_count * sizeof(double));
    result->operations = (pll_operation_t *) malloc(inner_nodes_count * sizeof(pll_operation_t));
    result->ops_count = 0;
    result->matrix_indices = (unsigned int *) malloc(branch_count * sizeof(unsigned int));
    result->matrix_count = 0;

    // TODO: Fill the operations array, note that we have to sum up the branch lengths that lie on a path...
    // (by going up through the parents, this gets a bit tricky when encountering reticulations because then we also have to figure out whether the edge belongs to the current tree)

	return PLL_SUCCESS;
}
