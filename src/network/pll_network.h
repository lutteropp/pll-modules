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

 Contact: Diego Darriba <Diego.Darriba@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */
#ifndef PLL_NETWORK_H_
#define PLL_NETWORK_H_

#ifndef PLL_H_
#define PLL_H_
#include "pll.h"
#endif

typedef struct pll_displayed_tree_s
{
  double * branch_lengths;
  pll_operation_t * operations;
  unsigned int ops_count;
  unsigned int * matrix_indices;
  unsigned int matrix_count;
} pll_displayed_tree_t;

PLL_EXPORT int pll_rnetwork_tree_buildarrays(pll_rnetwork_t * network, uint64_t tree_number, pll_displayed_tree_t * result);

#endif /* PLL_NETWORK_H_ */
