#include "pll_network.h"

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
