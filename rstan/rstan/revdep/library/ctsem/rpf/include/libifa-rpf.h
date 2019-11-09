/*
  Copyright 2012-2017 Joshua Nathaniel Pritikin and contributors

  libifa-rpf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _LIBIFA_RPF_
#define _LIBIFA_RPF_

enum RPF_ISpec {
  RPF_ISpecID,
  RPF_ISpecOutcomes,
  RPF_ISpecDims,
  RPF_ISpecCount
};

#define RPF_ISpecFacts RPF_ISpecDims

typedef int (*rpf_numSpec_t)(const double *spec);
typedef int (*rpf_numParam_t)(const double *spec);
typedef void (*rpf_paramInfo_t)(const double *spec, const int param,
				const char **type, double *upper, double *lower);
typedef void (*rpf_prob_t)(const double *spec,
			   const double *param, const double *th,
			   double *out);
typedef void (*rpf_dLL1_t)(const double *spec,
			   const double *param,
			   const double *where,
			   const double *weight, double *out);
typedef void (*rpf_dLL2_t)(const double *spec, const double *param, double *out);
typedef void (*rpf_rescale_t)(const double *spec, double *param, const int *paramMask,
			      const double *mean, const double *choleskyCov);
typedef void (*rpf_dTheta_t)(const double *spec, const double *param,
			     const double *where, const double *dir,
			     double *grad, double *hess);

struct rpf {
  const char name[10];
  rpf_numSpec_t numSpec;
  rpf_numParam_t numParam;
  rpf_paramInfo_t paramInfo;
  rpf_prob_t prob;
  rpf_prob_t logprob;
  rpf_dLL1_t dLL1;
  rpf_dLL2_t dLL2;
  rpf_dTheta_t dTheta;
  rpf_rescale_t rescale;
};

/* R_GetCCallable */
#define LIBIFA_RPF_API_VERSION 420832   /* this is a random number */

typedef void (*get_librpf_t)(int version, int *numModels, const struct rpf **model);

#endif
