#include "pll_network.h"

#include "../pllmod_common.h"

int cb_full_network_traversal(pll_rnetwork_node_t * node) {
	(void) node;
	return 1;
}

pll_rnetwork_node_t * go_down_recursive(pll_rnetwork_node_t * node, int * present, uint64_t tree_number) {
	if (!node) {
		return NULL;
	}
	if (node->is_reticulation) {
		return go_down_recursive(node->child, present, tree_number);
	}
	if (present[node->idx]) {
		return node;
	}

	// check left child and right child first
	if (present[node->left->idx]) {
		return node->left;
	}
	if (present[node->right->idx]) {
		return node->right;
	}

	if (pll_rnetwork_can_go_tree(node, node->left, tree_number)) {
		pll_rnetwork_node_t * try_left = go_down_recursive(node->left, present, tree_number);
		if (try_left) {
			return try_left;
		}
	} else if (pll_rnetwork_can_go_tree(node, node->right, tree_number)) {
		return go_down_recursive(node->right, present, tree_number);
	}
	return NULL; // we should never end up here though
}

double collect_branch_length_to_first_present_parent(pll_rnetwork_node_t * node, int * present, uint64_t tree_number) {
	if (node->is_reticulation) {
		return PLL_FAILURE; // this should not happen as reticulation nodes should not be in the post-order tree traversal.
	}
	double branch_sum = node->length;
	pll_rnetwork_node_t* parent = node->parent;
	while (!present[parent->idx]) {
		node = parent;
		if (node->is_reticulation) {
			if (pll_rnetwork_can_go_tree(node->first_parent, node, tree_number)) {
				branch_sum += node->first_parent_length;
				parent = node->first_parent;
			} else {
				branch_sum += node->second_parent_length;
				parent = node->second_parent;
			}
		} else {
			branch_sum += node->length;
			parent = node->parent;
		}
	}
	return branch_sum;
}

PLL_EXPORT int pllmod_rnetwork_tree_build_operations(pll_rnetwork_t * network, uint64_t tree_number, pll_displayed_tree_t * result,
		unsigned int fake_clv_index, unsigned int fake_pmatrix_index) { // TODO: FIXME: replace path collapsing by using the fake node...
	unsigned int nodes_count = network->reticulation_count + network->inner_tree_count + network->tip_count;
	unsigned int inner_nodes_count = network->reticulation_count + network->inner_tree_count;
	unsigned int branch_count = network->inner_tree_count * 2 + network->reticulation_count;

	pll_rnetwork_node_t ** trav_buffer = (pll_rnetwork_node_t **) malloc(nodes_count * sizeof(pll_rnetwork_node_t *));
	unsigned int trav_size;
	if (!pll_rnetwork_tree_traverse(network, PLL_TREE_TRAVERSE_POSTORDER, cb_full_network_traversal, trav_buffer, &trav_size,
			tree_number)) {
		return PLL_FAILURE;
	}

	unsigned int i;

	result->operations = (pll_operation_t *) malloc(inner_nodes_count * sizeof(pll_operation_t));
	result->ops_count = 0;

	for (i = 0; i < trav_size; ++i) {
		pll_rnetwork_node_t * node = trav_buffer[i];

		if (node->left != NULL || node->right!= NULL || node->child != NULL) // inner tree node
				{
			result->operations[result->ops_count].parent_clv_index = node->clv_index;
			result->operations[result->ops_count].parent_scaler_index = node->scaler_index;

			pll_rnetwork_node_t * child1 = NULL;
			pll_rnetwork_node_t * child2 = NULL;
			if (!node->is_reticulation) {
				child1 = node->left;
				child2 = node->right;
				if (child1->is_reticulation) {
					if (!pll_rnetwork_can_go_tree(node, child1, tree_number)) {
						child1 = NULL;
					}
				}
				if (child2->is_reticulation) {
					if (!pll_rnetwork_can_go_tree(node, child2, tree_number)) {
						child2 = NULL;
					}
				}
			} else {
				child1 = node->child;
			}

			if (child1) {
				result->operations[result->ops_count].child1_clv_index = child1->clv_index;
				result->operations[result->ops_count].child1_scaler_index = child1->scaler_index;
				result->operations[result->ops_count].child1_matrix_index = child1->pmatrix_index;
			} else {
				result->operations[result->ops_count].child1_clv_index = fake_clv_index;
				result->operations[result->ops_count].child1_scaler_index = -1;
				result->operations[result->ops_count].child1_matrix_index = fake_pmatrix_index;
			}

			if (child2) {
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
	return PLL_SUCCESS;
}

PLL_EXPORT pll_rnetwork_t * pllmod_rnetwork_create_random(unsigned int taxa_count, const char * const * names, unsigned int random_seed) { // TODO: I think this is the problem! Why is this function body empty???
	return PLL_FAILURE;
}

PLL_EXPORT int pllmod_rnetwork_extend_random(pll_rnetwork_t * network, unsigned int ext_taxa_count, const char * const * ext_names,
		unsigned int random_seed) {
	return PLL_FAILURE;
}

PLL_EXPORT
pll_rnetwork_t * pllmod_rnetwork_create_parsimony(unsigned int taxon_count, unsigned int seq_length, char * const * names,
		char * const * sequences, const unsigned int * site_weights, const pll_state_t * map, unsigned int states, unsigned int attributes,
		unsigned int random_seed, unsigned int * score) {
	return PLL_FAILURE;
}

pll_rnetwork_t * pllmod_rnetwork_create_parsimony_multipart(unsigned int taxon_count, char * const * taxon_names,
		unsigned int partition_count, pll_partition_t * const * partitions, unsigned int random_seed, unsigned int * score) {
	return PLL_FAILURE;
}

PLL_EXPORT pll_rnetwork_t * pllmod_rnetwork_resolve_multi(const pll_rnetwork_t * multi_network, unsigned int random_seed,
		int * clv_index_map) {
	return PLL_FAILURE;
}

PLL_EXPORT int pllmod_rnetwork_is_tip(const pll_rnetwork_node_t * node) {
	return (node->left == NULL && node->right == NULL && node->child == NULL);
}

PLL_EXPORT void pllmod_rnetwork_set_length(pll_rnetwork_node_t * edge, double length) {
	edge->length = length;
}

PLL_EXPORT void pllmod_rnetwork_set_length_recursive(pll_rnetwork_t * network, double length, int missing_only) {
	/* set branch lengths */
	unsigned int i;
	unsigned int tip_count = network->tip_count;
	unsigned int inner_count = network->inner_tree_count + network->reticulation_count;
	for (i = 0; i < tip_count + inner_count; ++i) {
		pll_rnetwork_node_t * node = network->nodes[i];
		if (!node->length || !missing_only)
			pllmod_rnetwork_set_length(node, length);
	}
}

PLL_EXPORT int pllmod_rnetwork_outgroup_root(pll_rnetwork_t * network, unsigned int * outgroup_tip_ids, unsigned int outgroup_size,
		int add_root_node) {
	return PLL_FAILURE;
}
