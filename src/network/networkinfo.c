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

#include "../pllmod_common.h"
#include "pll_network.h"

static int networkinfo_check_network(pllmod_networkinfo_t * networkinfo,
                               pll_unetwork_t * network);
static int networkinfo_init_network(pllmod_networkinfo_t * networkinfo);

/* a callback function for performing a full traversal */
static int cb_full_traversal_network(pll_unetwork_node_t * node)
{
  PLLMOD_UNUSED(node);
  return PLL_SUCCESS;
}

/* a callback function for performing a partial traversal on invalid CLVs */
static int cb_partial_traversal_network(pll_unetwork_node_t * node)
{
  /* do not include tips */
  if (!node->next) return PLL_FAILURE;

  pllmod_networkinfo_t * networkinfo = (pllmod_networkinfo_t *) node->data;

  /* if clv is invalid, traverse the subnetwork to compute it */
  if (networkinfo->active_partition == PLLMOD_NETWORKINFO_PARTITION_ALL)
  {
    /* check if at least one per-partition CLV is invalid */
    for (unsigned int i = 0; i < networkinfo->init_partition_count; ++i)
    {
      unsigned int p = networkinfo->init_partition_idx[i];
      if (networkinfo->clv_valid[p][node->node_index] == 0)
        return PLL_SUCCESS;
    }

    /* CLVs for all partitions are valid -> skip subnetwork */
    return PLL_FAILURE;
  }
  else
    return (networkinfo->clv_valid[networkinfo->active_partition][node->node_index] == 0);
}

static int networkinfo_partition_active(pllmod_networkinfo_t * networkinfo,
                                     unsigned int partition_index)
{
  return (networkinfo->active_partition == PLLMOD_NETWORKINFO_PARTITION_ALL ||
		  networkinfo->active_partition == (int) partition_index);
}

PLL_EXPORT pllmod_networkinfo_t * pllmod_networkinfo_create(pll_unetwork_t * network, unsigned int tips, unsigned int partitions,
		int brlen_linkage) {
	/* create networkinfo instance */
	pllmod_networkinfo_t * networkinfo;

	if (!(networkinfo = (pllmod_networkinfo_t *) calloc(1, sizeof(pllmod_networkinfo_t)))) {
		pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory for networkinfo\n");
		return NULL;
	}

	/* save dimensions & options */
	networkinfo->partition_count = partitions;
	networkinfo->brlen_linkage = brlen_linkage;

	/* store it in networkinfo */
	networkinfo->network = network;

	if (!networkinfo->network || !networkinfo_check_network(networkinfo, networkinfo->network)) {
		assert(pll_errno);
		free(networkinfo);
		return NULL;
	}

	/* compute some derived dimensions */
	unsigned int inner_nodes_count = networkinfo->network->inner_tree_count + networkinfo->network->reticulation_count;
	unsigned int nodes_count = inner_nodes_count + tips;
	unsigned int branch_count = networkinfo->network->edge_count;
	networkinfo->subnode_count = tips + 3 * inner_nodes_count;

	/* allocate a buffer for storing pointers to nodes of the network in postorder
	 traversal */
	networkinfo->travbuffer = (pll_unetwork_node_t **) malloc(nodes_count * sizeof(pll_unetwork_node_t *));

	/* allocate a buffer for matrix indices */
	networkinfo->matrix_indices = (unsigned int *) malloc(branch_count * sizeof(unsigned int));

	/* allocate a buffer for operations (parent/child clv indices) */
	networkinfo->operations = (pll_operation_t *) malloc(inner_nodes_count * sizeof(pll_operation_t));

	networkinfo->subnodes = (pll_unetwork_node_t **) malloc(networkinfo->subnode_count * sizeof(pll_unetwork_node_t *));

	/* check memory allocation */
	if (!networkinfo->travbuffer || !networkinfo->matrix_indices || !networkinfo->operations || !networkinfo->subnodes) {
		pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory for networkinfo structures\n");
		return NULL;
	}

	/* allocate arrays for storing per-partition info */
	networkinfo->partitions = (pll_partition_t **) calloc(partitions, sizeof(pll_partition_t *));
	networkinfo->params_to_optimize = (int *) calloc(partitions, sizeof(int));
	networkinfo->alphas = (double *) calloc(partitions, sizeof(double));
	networkinfo->gamma_mode = (int *) calloc(partitions, sizeof(int));
	networkinfo->param_indices = (unsigned int **) calloc(partitions, sizeof(unsigned int*));
	networkinfo->subst_matrix_symmetries = (int **) calloc(partitions, sizeof(int*));
	networkinfo->branch_lengths = (double **) calloc(partitions, sizeof(double*));
	networkinfo->deriv_precomp = (double **) calloc(partitions, sizeof(double*));
	networkinfo->clv_valid = (char **) calloc(partitions, sizeof(char*));
	networkinfo->pmatrix_valid = (char **) calloc(partitions, sizeof(char*));
	networkinfo->partition_loglh = (double *) calloc(partitions, sizeof(double));

	networkinfo->init_partition_count = 0;
	networkinfo->init_partition_idx = (unsigned int *) calloc(partitions, sizeof(unsigned int));
	networkinfo->init_partitions = (pll_partition_t **) calloc(partitions, sizeof(pll_partition_t *));

	/* allocate array for storing linked/average branch lengths */
	networkinfo->linked_branch_lengths = (double *) malloc(branch_count * sizeof(double));

	/* allocate branch length scalers if needed */
	if (brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED)
		networkinfo->brlen_scalers = (double *) calloc(partitions, sizeof(double));
	else
		networkinfo->brlen_scalers = NULL;

	/* check memory allocation */
	if (!networkinfo->partitions || !networkinfo->alphas || !networkinfo->param_indices || !networkinfo->subst_matrix_symmetries
			|| !networkinfo->branch_lengths || !networkinfo->deriv_precomp || !networkinfo->clv_valid || !networkinfo->pmatrix_valid
			|| !networkinfo->linked_branch_lengths || !networkinfo->partition_loglh || !networkinfo->gamma_mode
			|| !networkinfo->init_partition_idx || (brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED && !networkinfo->brlen_scalers)) {
		pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory for networkinfo arrays\n");
		return NULL;
	}

	unsigned int p;
	for (p = 0; p < partitions; ++p) {
		/* use mean GAMMA rates per default */
		networkinfo->gamma_mode[p] = PLL_GAMMA_RATES_MEAN;

		/* allocate arrays for storing the per-partition branch lengths */
		if (brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED) {
			networkinfo->branch_lengths[p] = (double *) malloc(branch_count * sizeof(double));
		} else
			networkinfo->branch_lengths[p] = networkinfo->linked_branch_lengths;

		/* initialize all branch length scalers to 1 */
		if (networkinfo->brlen_scalers)
			networkinfo->brlen_scalers[p] = 1.;

		/* check memory allocation */
		if (!networkinfo->branch_lengths[p]) {
			pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory for arrays for partition %d\n", p);
			return NULL;
		}
	}

	/* by default, work with all partitions */
	networkinfo->active_partition = PLLMOD_NETWORKINFO_PARTITION_ALL;

	/* needs to be here since we use some of the arrays allocated above */
	if (!networkinfo_init_network(networkinfo)) {
		pllmod_networkinfo_destroy(networkinfo);
		assert(pll_errno);
		return NULL;
	}

	assert(networkinfo->network && networkinfo->network->tip_count == tips);

	return networkinfo;
}

PLL_EXPORT
int pllmod_networkinfo_set_parallel_context(pllmod_networkinfo_t * networkinfo, void * parallel_context,
		void (*parallel_reduce_cb)(void *, double *, size_t, int op)) {
	networkinfo->parallel_context = parallel_context;
	networkinfo->parallel_reduce_cb = parallel_reduce_cb;

	return PLL_SUCCESS;
}

PLL_EXPORT int pllmod_networkinfo_init_partition(pllmod_networkinfo_t * networkinfo, unsigned int partition_index,
		pll_partition_t * partition, int params_to_optimize, int gamma_mode, double alpha, const unsigned int * param_indices,
		const int * subst_matrix_symmetries) {
	if (!networkinfo) {
		pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Networkinfo structure is NULL\n");
		return PLL_FAILURE;
	} else if (partition_index >= networkinfo->partition_count) {
		pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Partition %d is out of bounds\n", partition_index);
		return PLL_FAILURE;
	} else if (networkinfo->partitions[partition_index]) {
		pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Partition %d is already initialized\n", partition_index);
		return PLL_FAILURE;
	}

	unsigned int local_partition_index = networkinfo->init_partition_count++;
	networkinfo->partitions[partition_index] = partition;
	networkinfo->init_partitions[local_partition_index] = partition;
	networkinfo->init_partition_idx[local_partition_index] = partition_index;
	networkinfo->params_to_optimize[partition_index] = params_to_optimize;
	networkinfo->gamma_mode[partition_index] = gamma_mode;
	networkinfo->alphas[partition_index] = alpha;

	/* compute some derived dimensions */
	unsigned int inner_nodes_count = networkinfo->network->tip_count - 2 + networkinfo->network->reticulation_count;
	unsigned int nodes_count = inner_nodes_count + networkinfo->network->tip_count;
	unsigned int branch_count = nodes_count - 1;
	unsigned int pmatrix_count = branch_count;
	unsigned int unetwork_count = inner_nodes_count * 3 + networkinfo->network->tip_count;

	/* allocate invalidation arrays */
	networkinfo->clv_valid[partition_index] = (char *) calloc(unetwork_count, sizeof(char));
	networkinfo->pmatrix_valid[partition_index] = (char *) calloc(pmatrix_count, sizeof(char));

	/* check memory allocation */
	if (!networkinfo->clv_valid[partition_index] || !networkinfo->pmatrix_valid[partition_index]) {
		pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory for parameter indices\n");
		return PLL_FAILURE;
	}

	/* allocate param_indices array and initialize it to all 0s,
	 * i.e. per default, all rate categories will use
	 * the same substitution matrix and same base frequencies */
	networkinfo->param_indices[partition_index] = (unsigned int *) calloc(partition->rate_cats, sizeof(unsigned int));

	/* check memory allocation */
	if (!networkinfo->param_indices[partition_index]) {
		pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory for parameter indices\n");
		return PLL_FAILURE;
	}

	/* if param_indices were provided, use them instead of default */
	if (param_indices)
		memcpy(networkinfo->param_indices[partition_index], param_indices, partition->rate_cats * sizeof(unsigned int));

	/* copy substitution rate matrix symmetries, if any */
	if (subst_matrix_symmetries) {
		const unsigned int symm_size = (partition->states * (partition->states - 1) / 2) * sizeof(int);
		networkinfo->subst_matrix_symmetries[partition_index] = (int *) malloc(symm_size);

		/* check memory allocation */
		if (!networkinfo->subst_matrix_symmetries[partition_index]) {
			pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory for substitution scheme\n");
			return PLL_FAILURE;
		}

		memcpy(networkinfo->subst_matrix_symmetries[partition_index], subst_matrix_symmetries, symm_size);
	} else
		networkinfo->subst_matrix_symmetries[partition_index] = NULL;

	/* allocate memory for derivative precomputation table */
	unsigned int sites_alloc = partition->sites;
	if (partition->attributes & PLL_ATTRIB_AB_FLAG)
		sites_alloc += partition->states;
	unsigned int precomp_size = sites_alloc * partition->rate_cats * partition->states_padded;

	networkinfo->deriv_precomp[partition_index] = (double *) pll_aligned_alloc(precomp_size * sizeof(double), partition->alignment);

	if (!networkinfo->deriv_precomp[partition_index]) {
		pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory for derivative buffers\n");
		return PLL_FAILURE;
	}

	memset(networkinfo->deriv_precomp[partition_index], 0, precomp_size * sizeof(double));

	return PLL_SUCCESS;
}

PLL_EXPORT int pllmod_networkinfo_set_active_partition(pllmod_networkinfo_t * networkinfo, int partition_index) {
	if (partition_index != PLLMOD_NETWORKINFO_PARTITION_ALL && partition_index >= (int) networkinfo->partition_count) {
		pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Partition %d is out of bounds\n", partition_index);
		return PLL_FAILURE;
	} else {
		networkinfo->active_partition = partition_index;
		return PLL_SUCCESS;
	}
}

PLL_EXPORT int pllmod_networkinfo_set_root(pllmod_networkinfo_t * networkinfo, pll_unetwork_node_t * root) {
	if (!networkinfo || !root || root->data != (void *) networkinfo) {
		pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Invalid root node!\n");
		return PLL_FAILURE;
	}

	/* root must be an inner node! */
	networkinfo->root = pllmod_unetwork_is_tip(root) ? root->back : root;
	networkinfo->network->vroot = networkinfo->root;

	return PLL_SUCCESS;
}

PLL_EXPORT
int pllmod_networkinfo_get_branch_length_all(const pllmod_networkinfo_t * networkinfo, const pll_unetwork_node_t * edge, double * lengths) {
	unsigned int pmatrix_index = edge->pmatrix_index;

	if (networkinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED) {
		for (unsigned int p = 0; p < networkinfo->partition_count; ++p) {
			if (networkinfo->partitions[p])
				*lengths++ = networkinfo->branch_lengths[p][pmatrix_index];
		}
	} else
		lengths[0] = networkinfo->branch_lengths[0][pmatrix_index];

	return PLL_SUCCESS;
}

PLL_EXPORT int pllmod_networkinfo_set_branch_length(pllmod_networkinfo_t * networkinfo, pll_unetwork_node_t * edge, double length) {
	assert(networkinfo->brlen_linkage != PLLMOD_COMMON_BRLEN_UNLINKED);
	return pllmod_networkinfo_set_branch_length_all(networkinfo, edge, &length);
}

PLL_EXPORT
int pllmod_networkinfo_set_branch_length_all(pllmod_networkinfo_t * networkinfo, pll_unetwork_node_t * edge, const double * lengths) {
	unsigned int pmatrix_index = edge->pmatrix_index;
	if (networkinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED) {
		for (unsigned int p = 0; p < networkinfo->partition_count; ++p) {
			if (networkinfo->partitions[p])
				networkinfo->branch_lengths[p][pmatrix_index] = *lengths++;
		}
	} else {
		networkinfo->branch_lengths[0][pmatrix_index] = lengths[0];
		pllmod_unetwork_set_length(edge, lengths[0]);
	}

#if 0
	pllmod_networkinfo_set_active_partition(networkinfo, PLLMOD_NETWORKINFO_PARTITION_ALL);

	/* invalidate p-matrices */
	pllmod_networkinfo_invalidate_pmatrix(networkinfo, edge);

	/* invalidate CLVs */
	if (edge->next)
	{
		pllmod_networkinfo_invalidate_clv(networkinfo, edge->next);
		pllmod_networkinfo_invalidate_clv(networkinfo, edge->next->next);
	}
	if (edge->back->next)
	{
		pllmod_networkinfo_invalidate_clv(networkinfo, edge->back->next);
		pllmod_networkinfo_invalidate_clv(networkinfo, edge->back->next->next);
	}
#endif

	return PLL_SUCCESS;
}

PLL_EXPORT
int pllmod_networkinfo_set_branch_length_partition(pllmod_networkinfo_t * networkinfo, pll_unetwork_node_t * edge, int partition_index,
		double length) {
	unsigned int pmatrix_index = edge->pmatrix_index;
	const int old_active_partition = networkinfo->active_partition;

	if (!pllmod_networkinfo_set_active_partition(networkinfo, partition_index))
		return PLL_FAILURE;

	if (partition_index != PLLMOD_NETWORKINFO_PARTITION_ALL)
		networkinfo->branch_lengths[partition_index][pmatrix_index] = length;
	else if (networkinfo->brlen_linkage != PLLMOD_COMMON_BRLEN_UNLINKED) {
		pllmod_unetwork_set_length(edge, length);
		networkinfo->branch_lengths[0][pmatrix_index] = length;
	} else {
		for (unsigned int p = 0; p < networkinfo->partition_count; ++p) {
			if (networkinfo->partitions[p])
				networkinfo->branch_lengths[p][pmatrix_index] = length;
		}
	}

	networkinfo->active_partition = old_active_partition;

#if 0
	/* invalidate p-matrices */
	pllmod_networkinfo_invalidate_pmatrix(networkinfo, edge);

	/* invalidate CLVs */
	if (edge->next)
	{
		pllmod_networkinfo_invalidate_clv(networkinfo, edge->next);
		pllmod_networkinfo_invalidate_clv(networkinfo, edge->next->next);
	}
	if (edge->back->next)
	{
		pllmod_networkinfo_invalidate_clv(networkinfo, edge->back->next);
		pllmod_networkinfo_invalidate_clv(networkinfo, edge->back->next->next);
	}
#endif

	return PLL_SUCCESS;
}

PLL_EXPORT
pll_unetwork_t * pllmod_networkinfo_get_partition_network(const pllmod_networkinfo_t * networkinfo, int partition_index) {
	pll_unetwork_t * pnetwork = pll_unetwork_clone(networkinfo->network);

	// set partition-specific branch lengths
	if (networkinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED) {
		unsigned int node_count = pnetwork->tip_count + pnetwork->inner_tree_count + pnetwork->reticulation_count;
		unsigned int edge_count = pnetwork->edge_count;
		const double * brlens = networkinfo->branch_lengths[partition_index];
		for (unsigned int i = 0; i < node_count; ++i) {
			pll_unetwork_node_t * snode = pnetwork->nodes[i];
			do {
				unsigned int pmat_idx = snode->pmatrix_index;
				if (pmat_idx > edge_count) {
					pllmod_set_error(PLLMOD_NETWORK_ERROR_INVALID_NETWORK, "p-matrix index out of bounds (%u). "
							"networkinfo structure requires that each branch "
							"is assigned a unique p-matrix index "
							"between 0 and branch_count-1\n", pmat_idx);
					return NULL;
				}
				snode->length = brlens[pmat_idx];
				snode = snode->next;
			} while (snode && snode != pnetwork->nodes[i]);
		}
	}

	return pnetwork;
}

PLL_EXPORT
pllmod_networkinfo_topology_t * pllmod_networkinfo_get_topology(const pllmod_networkinfo_t * networkinfo,
		pllmod_networkinfo_topology_t * topol) {
	unsigned int brlen_set_count = (networkinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED) ? networkinfo->init_partition_count : 0;

	if (!topol) {
		topol = (pllmod_networkinfo_topology_t *) calloc(1, sizeof(pllmod_networkinfo_topology_t));
		if (!topol) {
			pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory for topology structure\n");
			return PLL_FAILURE;
		}

		topol->edge_count = networkinfo->network->edge_count;
		topol->brlen_set_count = brlen_set_count;
		topol->root_index = networkinfo->root->node_index;
		topol->edges = (pllmod_networkinfo_edge_t *) calloc(topol->edge_count, sizeof(pllmod_networkinfo_edge_t));
		if (!topol->edges) {
			pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory for topology buffers\n");
			pllmod_networkinfo_destroy_topology(topol);
			return PLL_FAILURE;
		}

		if (topol->brlen_set_count) {
			topol->branch_lengths = (double **) calloc(topol->brlen_set_count, sizeof(double *));
			if (!topol->branch_lengths) {
				pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory for branch length buffers\n");
				pllmod_networkinfo_destroy_topology(topol);
				return PLL_FAILURE;
			}

			for (unsigned int i = 0; i < topol->brlen_set_count; ++i) {
				topol->branch_lengths[i] = (double *) calloc(topol->edge_count, sizeof(double));
				if (!topol->branch_lengths[i]) {
					pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Cannot allocate memory for branch length buffers\n");
					pllmod_networkinfo_destroy_topology(topol);
					return PLL_FAILURE;
				}
			}
		}
	} else if (networkinfo->network->edge_count != topol->edge_count || brlen_set_count != topol->brlen_set_count) {
		pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Incompatible networkinfo_topology structure!\n");
		return PLL_FAILURE;
	}

	// save topology as a list of edges
	unsigned int edge_num = 0;
	for (unsigned int i = 0; i < networkinfo->subnode_count; ++i) {
		const pll_unetwork_node_t * snode = networkinfo->subnodes[i];
		if (snode->node_index < snode->back->node_index) {
			topol->edges[edge_num].pmatrix_index = snode->pmatrix_index;
			topol->edges[edge_num].left_index = snode->node_index;
			topol->edges[edge_num].right_index = snode->back->node_index;
			topol->edges[edge_num].brlen = snode->length;
			edge_num++;
		}
	}
	assert(edge_num == topol->edge_count);

	// save unlinked branch lengths
	for (unsigned int i = 0; i < topol->brlen_set_count; ++i) {
		unsigned int p = networkinfo->init_partition_idx[i];
		memcpy(topol->branch_lengths[i], networkinfo->branch_lengths[p], topol->edge_count * sizeof(double));
	}

	return topol;
}

PLL_EXPORT
int pllmod_networkinfo_set_topology(pllmod_networkinfo_t * networkinfo, const pllmod_networkinfo_topology_t * topol) {
	unsigned int brlen_set_count = (networkinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED) ? networkinfo->init_partition_count : 0;

	if (!networkinfo || !topol) {
		pllmod_set_error(PLL_ERROR_PARAM_INVALID, "networkinfo and/or topology is empty!\n");
		return PLL_FAILURE;
	}
	if (networkinfo->network->edge_count != topol->edge_count) {
		pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Incompatible topology: edge count differs!\n");
		return PLL_FAILURE;
	}
	if (brlen_set_count != topol->brlen_set_count) {
		pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Incompatible topology: brlen set count differs!!\n");
		return PLL_FAILURE;
	}

	// re-connect branches and reset pmatrix indices
	for (unsigned int i = 0; i < topol->edge_count; ++i) {
		const pllmod_networkinfo_edge_t * edge = &topol->edges[i];
		pll_unetwork_node_t * left_node = networkinfo->subnodes[edge->left_index];
		pll_unetwork_node_t * right_node = networkinfo->subnodes[edge->right_index];
		pllmod_unetwork_connect_nodes(left_node, right_node, edge->brlen, edge->prob);
		left_node->pmatrix_index = right_node->pmatrix_index = edge->pmatrix_index;

		//    printf("load edge: %u %u %u %f\n", left_node->node_index, right_node->node_index,
		//           left_node->pmatrix_index, left_node->length);

		// linked/scaled brlens -> set a global value for all partitions
		if (!topol->branch_lengths)
			networkinfo->branch_lengths[0][edge->pmatrix_index] = edge->brlen;
	}

	// restore unlinked branch lengths
	if (topol->branch_lengths) {
		for (unsigned int i = 0; i < topol->brlen_set_count; ++i) {
			unsigned int p = networkinfo->init_partition_idx[i];
			memcpy(networkinfo->branch_lengths[p], topol->branch_lengths[i], topol->edge_count * sizeof(double));
		}
	}

	// restore virtual root
	networkinfo->root = networkinfo->network->vroot = networkinfo->subnodes[topol->root_index];
	assert(networkinfo->root->node_index == topol->root_index);

	pllmod_networkinfo_invalidate_all(networkinfo);

	return PLL_SUCCESS;
}

PLL_EXPORT
int pllmod_networkinfo_destroy_topology(pllmod_networkinfo_topology_t * topol) {
	if (topol) {
		if (topol->edges)
			free(topol->edges);
		if (topol->branch_lengths) {
			for (unsigned int i = 0; i < topol->brlen_set_count; ++i)
				free(topol->branch_lengths[i]);
			free(topol->branch_lengths);
		}
		free(topol);
	}
	return PLL_SUCCESS;
}

PLL_EXPORT int pllmod_networkinfo_destroy_partition(pllmod_networkinfo_t * networkinfo, unsigned int partition_index) {
	if (!networkinfo) {
		pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Networkinfo structure is NULL\n");
		return PLL_FAILURE;
	} else if (partition_index >= networkinfo->partition_count) {
		pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Partition %d is out of bounds\n", partition_index);
		return PLL_FAILURE;
	}

	if (networkinfo->clv_valid[partition_index]) {
		free(networkinfo->clv_valid[partition_index]);
		networkinfo->clv_valid[partition_index] = NULL;
	}

	if (networkinfo->pmatrix_valid[partition_index]) {
		free(networkinfo->pmatrix_valid[partition_index]);
		networkinfo->pmatrix_valid[partition_index] = NULL;
	}

	if (networkinfo->param_indices[partition_index]) {
		free(networkinfo->param_indices[partition_index]);
		networkinfo->param_indices[partition_index] = NULL;
	}

	if (networkinfo->subst_matrix_symmetries[partition_index]) {
		free(networkinfo->subst_matrix_symmetries[partition_index]);
		networkinfo->subst_matrix_symmetries[partition_index] = NULL;
	}

	if (networkinfo->deriv_precomp[partition_index]) {
		free(networkinfo->deriv_precomp[partition_index]);
		networkinfo->deriv_precomp[partition_index] = NULL;
	}

	return PLL_SUCCESS;
}

PLL_EXPORT void pllmod_networkinfo_destroy(pllmod_networkinfo_t * networkinfo) {
	if (!networkinfo)
		return;

	/* deallocate traversal buffer, branch lengths array, matrix indices
	 array and operations */
	free(networkinfo->travbuffer);
	free(networkinfo->matrix_indices);
	free(networkinfo->operations);
	free(networkinfo->subnodes);

	/* destroy all structures allocated for the concrete PLL partition instance */
	unsigned int p;
	for (p = 0; p < networkinfo->partition_count; ++p) {
		if (networkinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
			free(networkinfo->branch_lengths[p]);

		pllmod_networkinfo_destroy_partition(networkinfo, p);
	}

	if (networkinfo->subst_matrix_symmetries)
		free(networkinfo->subst_matrix_symmetries);

	if (networkinfo->constraint)
		free(networkinfo->constraint);

	/* free invalidation arrays */
	free(networkinfo->clv_valid);
	free(networkinfo->pmatrix_valid);

	free(networkinfo->linked_branch_lengths);

	/* free alpha and param_indices arrays */
	free(networkinfo->params_to_optimize);
	free(networkinfo->alphas);
	free(networkinfo->gamma_mode);
	free(networkinfo->param_indices);
	free(networkinfo->branch_lengths);
	free(networkinfo->partition_loglh);
	free(networkinfo->deriv_precomp);

	if (networkinfo->brlen_scalers)
		free(networkinfo->brlen_scalers);

	/* deallocate partition array */
	free(networkinfo->partitions);
	free(networkinfo->init_partitions);
	free(networkinfo->init_partition_idx);

	if (networkinfo->network) {
		free(networkinfo->network->nodes);
		free(networkinfo->network->reticulation_nodes);
		free(networkinfo->network);
	}

	/* finally, deallocate networkinfo object itself */
	free(networkinfo);
}

PLL_EXPORT int pllmod_networkinfo_update_prob_matrices(pllmod_networkinfo_t * networkinfo, int update_all) {
	unsigned int p, m;
	unsigned int updated = 0;
	unsigned int pmatrix_count = networkinfo->network->edge_count;

	for (p = 0; p < networkinfo->partition_count; ++p) {
		/* only selected partitioned will be affected */
		if (networkinfo_partition_active(networkinfo, p) && networkinfo->partitions[p]) {

			for (m = 0; m < pmatrix_count; ++m) {
				if (networkinfo->pmatrix_valid[p][m] && !update_all)
					continue;

				double p_brlen = networkinfo->branch_lengths[p][m];
				if (networkinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_SCALED)
					p_brlen *= networkinfo->brlen_scalers[p];

				int ret = pll_update_prob_matrices(networkinfo->partitions[p], networkinfo->param_indices[p], &m, &p_brlen, 1);

				if (!ret)
					return PLL_FAILURE;

				networkinfo->pmatrix_valid[p][m] = 1;
				updated++;
			}
		}
	}

	return PLL_SUCCESS;
}

PLL_EXPORT void pllmod_networkinfo_invalidate_all(pllmod_networkinfo_t * networkinfo) {
	unsigned int i, m;
	unsigned int clv_count = networkinfo->network->tip_count + (networkinfo->network->tip_count - 2 + networkinfo->network->reticulation_count) * 3;
	unsigned int pmatrix_count = networkinfo->network->edge_count;

	for (i = 0; i < networkinfo->init_partition_count; ++i) {
		unsigned int p = networkinfo->init_partition_idx[i];

		/* only selected partitioned will be affected */
		if (networkinfo_partition_active(networkinfo, p)) {
			for (m = 0; m < pmatrix_count; ++m)
				networkinfo->pmatrix_valid[p][m] = 0;

			for (m = 0; m < clv_count; ++m)
				networkinfo->clv_valid[p][m] = 0;
		}
	}
}

PLL_EXPORT int pllmod_networkinfo_validate_clvs(pllmod_networkinfo_t * networkinfo, pll_unetwork_node_t ** travbuffer,
		unsigned int travbuffer_size) {
	for (unsigned int i = 0; i < networkinfo->init_partition_count; ++i) {
		unsigned int p = networkinfo->init_partition_idx[i];

		/* only selected partitioned will be affected */
		if (networkinfo_partition_active(networkinfo, p)) {
			for (unsigned int j = 0; j < travbuffer_size; ++j) {
				const pll_unetwork_node_t * node = travbuffer[j];
				if (node->next) {
					networkinfo->clv_valid[p][node->node_index] = 1;

					/* since we have only 1 CLV vector per inner node,
					 * we must invalidate CLVs for other 2 directions */
					networkinfo->clv_valid[p][node->next->node_index] = 0;
					networkinfo->clv_valid[p][node->next->next->node_index] = 0;
				}
			}
		}
	}

	return PLL_SUCCESS;
}

PLL_EXPORT void pllmod_networkinfo_invalidate_pmatrix(pllmod_networkinfo_t * networkinfo, const pll_unetwork_node_t * edge) {
	for (unsigned int i = 0; i < networkinfo->init_partition_count; ++i) {
		unsigned int p = networkinfo->init_partition_idx[i];
		if (networkinfo->pmatrix_valid[p] && networkinfo_partition_active(networkinfo, p))
			networkinfo->pmatrix_valid[p][edge->pmatrix_index] = 0;
	}
}

PLL_EXPORT void pllmod_networkinfo_invalidate_clv(pllmod_networkinfo_t * networkinfo, const pll_unetwork_node_t * edge) {
	for (unsigned int i = 0; i < networkinfo->init_partition_count; ++i) {
		unsigned int p = networkinfo->init_partition_idx[i];
		if (networkinfo->clv_valid[p] && networkinfo_partition_active(networkinfo, p))
			networkinfo->clv_valid[p][edge->node_index] = 0;
	}
}

PLL_EXPORT double pllmod_networkinfo_compute_loglh(pllmod_networkinfo_t * networkinfo, int incremental, int update_pmatrices) { // TODO: This still needs to be adapted to networks!!!
	/* network root must be an inner node! */
	assert(!pllmod_unetwork_is_tip(networkinfo->root));

	unsigned int traversal_size;
	unsigned int ops_count;
	unsigned int i, p;

	const double LOGLH_NONE = (double) NAN;
	double total_loglh = 0.0;
	const int old_active_partition = networkinfo->active_partition;

	/* NOTE: in unlinked brlen mode, up-to-date brlens for partition p
	 * have to be prefetched to networkinfo->branch_lengths[p] !!! */
	int collect_brlen = (networkinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED ? 0 : 1);

	pllmod_networkinfo_set_active_partition(networkinfo, PLLMOD_NETWORKINFO_PARTITION_ALL);

	/* we need full traversal in 2 cases: 1) update p-matrices, 2) update all CLVs */
	if (!incremental || (update_pmatrices && collect_brlen)) {
		/* perform a FULL postorder traversal of the unrooted network */
		if (!pll_unetwork_tree_traverse(networkinfo->network,
		PLL_TREE_TRAVERSE_POSTORDER, cb_full_traversal_network, networkinfo->travbuffer, &traversal_size, 0))
			return LOGLH_NONE;
	}

	/* update p-matrices if asked for */
	if (update_pmatrices) {
		if (collect_brlen) {
			assert(traversal_size == networkinfo->network->tip_count * 2 - 2);
			for (i = 0; i < traversal_size; ++i) {
				pll_unetwork_node_t * node = networkinfo->travbuffer[i];
				networkinfo->branch_lengths[0][node->pmatrix_index] = node->length;
			}
		}

		pllmod_networkinfo_update_prob_matrices(networkinfo, !incremental);
	}

	if (incremental) {
		/* compute partial traversal and update only invalid CLVs */
		if (!pll_unetwork_tree_traverse(networkinfo->network,
		PLL_TREE_TRAVERSE_POSTORDER, cb_partial_traversal_network, networkinfo->travbuffer, &traversal_size, 0))
			return LOGLH_NONE;
	}

	/* create operations based on partial traversal obtained above */
	pll_unetwork_create_operations(networkinfo->travbuffer, traversal_size,
	NULL,
	NULL, networkinfo->operations,
	NULL, &ops_count);

	networkinfo->counter += ops_count;

	//  printf("Traversal size (%s): %u\n", incremental ? "part" : "full", ops_count);

	/* iterate over all partitions (we assume that traversal is the same) */
	for (p = 0; p < networkinfo->partition_count; ++p) {
		if (!networkinfo->partitions[p]) {
			/* this partition will be computed by another thread(s) */
			networkinfo->partition_loglh[p] = 0.0;
			continue;
		}

		/* all subsequent operation will affect current partition only */
		pllmod_networkinfo_set_active_partition(networkinfo, (int) p);

		/* use the operations array to compute all ops_count inner CLVs. Operations
		 will be carried out sequentially starting from operation 0 towards
		 ops_count-1 */
		pll_update_partials(networkinfo->partitions[p], networkinfo->operations, ops_count);

		pllmod_networkinfo_validate_clvs(networkinfo, networkinfo->travbuffer, traversal_size);

		/* compute the likelihood on an edge of the unrooted network by specifying
		 the CLV indices at the two end-point of the branch, the probability
		 matrix index for the concrete branch length, and the index of the model
		 of whose frequency vector is to be used */
		networkinfo->partition_loglh[p] = pll_compute_edge_loglikelihood(networkinfo->partitions[p], networkinfo->root->clv_index,
				networkinfo->root->scaler_index, networkinfo->root->back->clv_index, networkinfo->root->back->scaler_index,
				networkinfo->root->pmatrix_index, networkinfo->param_indices[p],
				NULL);
	}

	/* sum up likelihood from all threads */
	if (networkinfo->parallel_reduce_cb) {
		networkinfo->parallel_reduce_cb(networkinfo->parallel_context, networkinfo->partition_loglh, p,
		PLLMOD_COMMON_REDUCE_SUM);
	}

	/* accumulate loglh by summing up over all the partitions */
	for (p = 0; p < networkinfo->partition_count; ++p)
		total_loglh += networkinfo->partition_loglh[p];

	/* restore original active partition */
	pllmod_networkinfo_set_active_partition(networkinfo, old_active_partition);

	assert(total_loglh < 0.);

	return total_loglh;
}

PLL_EXPORT double pllmod_networkinfo_compute_loglh_flex(pllmod_networkinfo_t * networkinfo, int incremental, int update_pmatrices) {
	return pllmod_networkinfo_compute_loglh(networkinfo, incremental, update_pmatrices);
}

PLL_EXPORT
int pllmod_networkinfo_scale_branches_all(pllmod_networkinfo_t * networkinfo, double scaler) {
	unsigned int i, p;

	for (i = 0; i < networkinfo->subnode_count; ++i)
		networkinfo->subnodes[i]->length *= scaler;

	if (networkinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED) {
		for (p = 0; p < networkinfo->partition_count; ++p) {
			for (i = 0; i < networkinfo->network->edge_count; ++i)
				networkinfo->branch_lengths[p][i] *= scaler;
		}

	} else {
		for (i = 0; i < networkinfo->network->edge_count; ++i)
			networkinfo->branch_lengths[0][i] *= scaler;
	}

	return PLL_SUCCESS;
}

PLL_EXPORT
int pllmod_networkinfo_scale_branches_partition(pllmod_networkinfo_t * networkinfo, unsigned int partition_idx, double scaler) {
	unsigned int i;

	if (networkinfo->brlen_linkage != PLLMOD_COMMON_BRLEN_UNLINKED) {
		pllmod_set_error(PLLMOD_NETWORK_ERROR_INVALID_NETWORK, "Per-partition branch length scaling is supported "
				"in unlinked branch length mode only.\n");
		return PLL_FAILURE;
	}

	if (partition_idx >= networkinfo->partition_count) {
		pllmod_set_error(PLL_ERROR_PARAM_INVALID, "Partition ID out of bounds\n");
		return PLL_FAILURE;
	}

	for (i = 0; i < networkinfo->network->edge_count; ++i)
		networkinfo->branch_lengths[partition_idx][i] *= scaler;

	return PLL_SUCCESS;
}

PLL_EXPORT
int pllmod_networkinfo_normalize_brlen_scalers(pllmod_networkinfo_t * networkinfo) {
	double sum_scalers = 0.;
	double sum_sites = 0.;
	unsigned int p;

	for (p = 0; p < networkinfo->partition_count; ++p) {
		if (networkinfo->partitions[p]) {
			const double pat_sites = networkinfo->partitions[p]->pattern_weight_sum;
			sum_sites += pat_sites;
			sum_scalers += networkinfo->brlen_scalers[p] * pat_sites;
		}
	}

	/* sum up scalers and sites from all threads */
	if (networkinfo->parallel_reduce_cb) {
		// TODO: although this reduce is unlikely to become a bottleneck, it can
		// be avoided by storing _total_ pattern_weight_sum in networkinfo
		networkinfo->parallel_reduce_cb(networkinfo->parallel_context, &sum_scalers, 1,
		PLLMOD_COMMON_REDUCE_SUM);

		networkinfo->parallel_reduce_cb(networkinfo->parallel_context, &sum_sites, 1,
		PLLMOD_COMMON_REDUCE_SUM);
	}

	const double mean_rate = sum_scalers / sum_sites;
	pllmod_networkinfo_scale_branches_all(networkinfo, mean_rate);
	for (p = 0; p < networkinfo->partition_count; ++p) {
		if (networkinfo->partitions[p])
			networkinfo->brlen_scalers[p] /= mean_rate;
	}

	return PLL_SUCCESS;
}

static int networkinfo_check_network(pllmod_networkinfo_t * networkinfo,
                               pll_unetwork_t * network)
{
  if (!networkinfo || !network)
  {
    pllmod_set_error(PLL_ERROR_PARAM_INVALID,
                     "Parameter is NULL\n");
    return PLL_FAILURE;
  }

  if (networkinfo->network->tip_count != network->tip_count)
  {
    pllmod_set_error(PLLMOD_NETWORK_ERROR_INVALID_NETWORK_SIZE,
                     "Invalid network size. Got %d instead of %d\n",
					 network->tip_count, networkinfo->network->tip_count);
    return PLL_FAILURE;
  }

  if (!network->binary)
  {
    pllmod_set_error(PLLMOD_NETWORK_ERROR_INVALID_NETWORK,
                     "Binary network expected\n");
    return PLL_FAILURE;
  }

  return PLL_SUCCESS;
}

static int networkinfo_init_network(pllmod_networkinfo_t * networkinfo)
{
  pll_unetwork_t * network = networkinfo->network;
  unsigned int node_count = network->tip_count + network->inner_tree_count + network->reticulation_count;
  unsigned int edge_count = network->edge_count;

  /* save virtual root */
  networkinfo->root = network->vroot;

  /* 1. save back pointer to treeinfo stuct in in eacn node's _data_ field
   * 2. collect branch length from the tree and store in a separate array
   *    indexed by pmatrix_index */
  for (unsigned int i = 0; i < node_count; ++i)
  {
    pll_unetwork_node_t * snode = network->nodes[i];
    do
    {
      unsigned int pmat_idx = snode->pmatrix_index;
      if (pmat_idx > edge_count)
      {
        pllmod_set_error(PLLMOD_NETWORK_ERROR_INVALID_NETWORK,
                         "p-matrix index out of bounds (%u). "
                         "networkinfo structure require that each branch "
                         "branch is assigned a unique p-matrix index "
                         "between 0 and branch_count-1\n",
                         pmat_idx);
        return PLL_FAILURE;
      }
      networkinfo->subnodes[snode->node_index] = snode;
      networkinfo->branch_lengths[0][pmat_idx] = snode->length;
      snode->data = networkinfo;
      snode = snode->next;
    }
    while (snode && snode != network->nodes[i]);
  }

  /* in unlinked branch length mode, copy brlen to other partitions */
  if (networkinfo->brlen_linkage == PLLMOD_COMMON_BRLEN_UNLINKED)
  {
    for (unsigned int i = 1; i < networkinfo->partition_count; ++i)
    {
      // TODO: only save brlens for initialized partitions
//      if (networkinfo->partitions[i])
      {
        memcpy(networkinfo->branch_lengths[i], networkinfo->branch_lengths[0],
               edge_count * sizeof(double));
      }
    }
  }

  return PLL_SUCCESS;
}

PLL_EXPORT int pllmod_networkinfo_set_network(pllmod_networkinfo_t * networkinfo, pll_unetwork_t * network) {
	if (!networkinfo_check_network(networkinfo, network))
		return PLL_FAILURE;

	unsigned int node_count = network->tip_count + network->inner_tree_count + network->reticulation_count;

	/* make a shallow copy of pll_unetwork_t structure, that is,
	 * pll_unetwork_nodes will not be cloned! */
	pll_unetwork_node_t ** nodes = networkinfo->network->nodes;
	memcpy(networkinfo->network, network, sizeof(pll_unetwork_t));
	networkinfo->network->nodes = nodes;
	memcpy(networkinfo->network->nodes, network->nodes, node_count * sizeof(pll_unetwork_node_t *));

	return networkinfo_init_network(networkinfo);
}

PLL_EXPORT int pllmod_networkinfo_set_constraint_clvmap(pllmod_networkinfo_t * networkinfo, const int * clv_index_map) {
	const unsigned int tip_count = networkinfo->network->tip_count;
	const unsigned int inner_count = networkinfo->network->inner_tree_count + networkinfo->network->reticulation_count;
	const size_t cons_size = (tip_count + inner_count) * sizeof(unsigned int);

	if (!networkinfo->constraint) {
		networkinfo->constraint = malloc(cons_size);
		if (!networkinfo->constraint) {
			pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Can't allocate memory for topological constraint\n");
			return PLL_FAILURE;
		}
	}

	assert(networkinfo->constraint);

	for (unsigned int i = 0; i < tip_count + inner_count; ++i) {
		const pll_unetwork_node_t * node = networkinfo->network->nodes[i];
		const unsigned int cons_group_id = clv_index_map[node->clv_index] + 1;

		networkinfo->constraint[node->clv_index] = cons_group_id;
	}

	return PLL_SUCCESS;
}

PLL_EXPORT int pllmod_networkinfo_set_constraint_network(pllmod_networkinfo_t * networkinfo, const pll_unetwork_t * cons_network) {
	unsigned int node_count = cons_network->tip_count * 2 - 2;
	int * clv_index_map = (int *) calloc(node_count, sizeof(int));
	int retval;

	if (networkinfo->network->tip_count < cons_network->tip_count) {
		pllmod_set_error(PLLMOD_NETWORK_ERROR_INVALID_NETWORK_SIZE, "Invalid network size. Got %d instead of %d\n", cons_network->tip_count,
				networkinfo->network->tip_count);
		return PLL_FAILURE;
	}

	pll_unetwork_t * bin_cons_network = pllmod_unetwork_resolve_multi(cons_network, 0, clv_index_map);

	if (!bin_cons_network || !clv_index_map) {
		assert(pll_errno);
		return PLL_FAILURE;
	}

	// HACK; copy clv group ids from inner nodes to tips
	for (unsigned int i = 0; i < bin_cons_network->tip_count; ++i) {
		pll_unetwork_node_t * node = bin_cons_network->nodes[i];
		assert(!node->next);
		clv_index_map[node->clv_index] = clv_index_map[node->back->clv_index];
	}

	pll_unetwork_destroy(bin_cons_network, NULL);

	if (networkinfo->network->tip_count > cons_network->tip_count) {
		// non-comprehensive constraint
		unsigned int free_count = networkinfo->network->tip_count - cons_network->tip_count;
		unsigned int ext_node_count = networkinfo->network->tip_count * 2 - 2;
		int * ext_clv_index_map = (int *) calloc(ext_node_count, sizeof(int));

		for (unsigned int i = 0; i < ext_node_count; ++i) {
			unsigned int clv_id = networkinfo->network->nodes[i]->clv_index;
			if (clv_id < cons_network->tip_count) {
				// copy map for tips in the constraint network and adjust clv_ids
				ext_clv_index_map[clv_id] = clv_index_map[clv_id] + free_count;
			} else if (clv_id < networkinfo->network->tip_count) {
				// mark new tips as freely movable
				ext_clv_index_map[clv_id] = -1;
			} else if (clv_id < node_count + free_count) {
				// update clv_ids for old inner nodes by adding free tip count
				unsigned int old_clv_id = clv_id - free_count;
				ext_clv_index_map[clv_id] = clv_index_map[old_clv_id] + free_count;
			} else if (clv_id < ext_node_count) {
				// mark new inner nodes as freely movable
				ext_clv_index_map[clv_id] = -1;
			} else
				assert(0);
		}

		free(clv_index_map);
		clv_index_map = ext_clv_index_map;
	}

	retval = pllmod_networkinfo_set_constraint_clvmap(networkinfo, clv_index_map);

	free(clv_index_map);

	return retval;
}

static unsigned int find_cons_id_network(pll_unetwork_node_t * node, const unsigned int * constraint, unsigned int s) {
	unsigned int cons_group_id = constraint[node->clv_index];
	if (!node->next || cons_group_id > 0)
		return cons_group_id;
	else {
		unsigned int left_id = find_cons_id_network(node->next->back, constraint, s);
		unsigned int right_id = find_cons_id_network(node->next->next->back, constraint, s);
		if (left_id == right_id)
			return left_id;
		else if (!s)
			return PLL_MAX(left_id, right_id);
		else
			return (left_id == 0 || left_id == s) ? right_id : left_id;
	}
}

PLL_EXPORT int pllmod_networkinfo_check_constraint(pllmod_networkinfo_t * networkinfo, pll_unetwork_node_t * subnetwork,
		pll_unetwork_node_t * regraft_edge) {
	if (networkinfo->constraint) {
		int res;
		unsigned int s = networkinfo->constraint[subnetwork->clv_index];
		s = s ? s : find_cons_id_network(subnetwork->back, networkinfo->constraint, 0);

		if (s) {
			unsigned int r1 = find_cons_id_network(regraft_edge, networkinfo->constraint, s);
			unsigned int r2 = find_cons_id_network(regraft_edge->back, networkinfo->constraint, s);

			res = (s == r1 || s == r2) ? PLL_SUCCESS : PLL_FAILURE;
		} else
			res = PLL_SUCCESS;

		return res;
	} else
		return PLL_SUCCESS;
}

static pllmod_ancestral_network_t * pllmod_networkinfo_create_ancestral(const pllmod_networkinfo_t * networkinfo) {
	unsigned int i;

	pllmod_ancestral_network_t * ancestral = (pllmod_ancestral_network_t *) calloc(1, sizeof(pllmod_ancestral_network_t));

	if (!ancestral)
		return NULL;

	ancestral->node_count = networkinfo->network->inner_tree_count + networkinfo->network->reticulation_count;
	ancestral->partition_count = networkinfo->init_partition_count;

	ancestral->nodes = (pll_unetwork_node_t **) calloc(ancestral->node_count, sizeof(pll_unetwork_node_t *));
	ancestral->partition_indices = (unsigned int *) calloc(ancestral->partition_count, sizeof(unsigned int));

	ancestral->probs = (double **) calloc(ancestral->node_count, sizeof(double *));

	size_t vector_size = 0;
	for (i = 0; i < ancestral->partition_count; ++i) {
		const pll_partition_t * partition = networkinfo->init_partitions[i];
		vector_size += partition->sites * partition->states;
	}

	ancestral->network = pll_unetwork_clone(networkinfo->network);

	for (i = 0; i < ancestral->node_count; ++i)
		ancestral->probs[i] = (double *) calloc(vector_size, sizeof(double));

	return ancestral;
}

PLL_EXPORT void pllmod_networkinfo_destroy_ancestral(pllmod_ancestral_network_t * ancestral) {
	unsigned int i;

	if (!ancestral)
		return;

	free(ancestral->nodes);
	free(ancestral->partition_indices);

	for (i = 0; i < ancestral->node_count; ++i)
		free(ancestral->probs[i]);

	pll_unetwork_destroy(ancestral->network, NULL);
	free(ancestral->probs);
}

PLL_EXPORT
pllmod_ancestral_network_t * pllmod_networkinfo_compute_ancestral(pllmod_networkinfo_t * networkinfo) {
	unsigned int i, p;
	unsigned int traversal_size;
	unsigned int node_count = networkinfo->network->tip_count + networkinfo->network->inner_tree_count
			+ networkinfo->network->reticulation_count;

	pll_unetwork_node_t * old_root = networkinfo->root;

	pll_unetwork_node_t ** travbuffer = (pll_unetwork_node_t **) calloc(node_count, sizeof(pll_unetwork_node_t *));

	pllmod_ancestral_network_t * ancestral = pllmod_networkinfo_create_ancestral(networkinfo);

	if (!travbuffer || !ancestral) {
		if (travbuffer)
			free(travbuffer);
		if (ancestral)
			pllmod_networkinfo_destroy_ancestral(ancestral);
		pllmod_set_error(PLL_ERROR_MEM_ALLOC, "Can't allocate memory for ancestral probabilities\n");
		return NULL;
	}

	/* perform a FULL postorder traversal of the unrooted network */
	if (!pll_unetwork_tree_traverse(ancestral->network,
	PLL_TREE_TRAVERSE_POSTORDER, cb_full_traversal_network, travbuffer, &traversal_size, 0)) {
		free(travbuffer);
		pllmod_networkinfo_destroy_ancestral(ancestral);
		return NULL;
	}

	assert(traversal_size == node_count);

#ifdef DEBUG
	printf("\nTraversal: ");
	for (i = 0; i < traversal_size; ++i)
	printf("%u ", travbuffer[i]->clv_index);
	printf("\n\n");
#endif

	/* collect internal nodes from the traversal */
	unsigned int n = 0;
	for (i = 0; i < traversal_size; ++i) {
		pll_unetwork_node_t * node = travbuffer[i];

		if (pllmod_unetwork_is_tip(node))
			continue;

		/* generate inner node label if it is empty */
		if (!node->label) {
			char buf[100];
			snprintf(buf, 100, "Node%u", n + 1);
			node->label = node->next->label = node->next->next->label = strdup(buf);
		}

		ancestral->nodes[n] = node;
		++n;
	}

	assert(n == ancestral->node_count);

	free(travbuffer);

	/* now compute ancestral state probs for every internal node */
	for (i = 0; i < ancestral->node_count; ++i) {
		pll_unetwork_node_t * node = ancestral->nodes[i];
		pll_unetwork_node_t * networkinfo_node = networkinfo->subnodes[node->node_index];
		double * ancp = ancestral->probs[i];

		networkinfo->root = networkinfo_node;
		pllmod_networkinfo_compute_loglh(networkinfo, 1, 1);

		for (p = 0; p < networkinfo->init_partition_count; ++p) {
			pll_partition_t * partition = networkinfo->init_partitions[p];
			size_t part_span = partition->sites * partition->states;
			unsigned int pidx = networkinfo->init_partition_idx[p];

			ancestral->partition_indices[p] = pidx;

			if (!pll_compute_node_ancestral(partition, node->clv_index, node->scaler_index, node->back->clv_index, node->back->scaler_index,
					node->pmatrix_index, networkinfo->param_indices[pidx], ancp)) {
				pllmod_networkinfo_destroy_ancestral(ancestral);
				return NULL;
			}

			ancp += part_span;
		}
	}

	networkinfo->root = old_root;

	return ancestral;
}
