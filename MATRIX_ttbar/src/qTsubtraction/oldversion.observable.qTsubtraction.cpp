#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"
#include "../qTsubtraction/header/more.h"

void observable_set::determine_psp_weight_VT(phasespace_set & psi){
  Logger logger("observable_set::determine_psp_weight_VT");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  this_psp_weight = integrand;
  this_psp_weight2 = pow(this_psp_weight, 2.);
  step_sum_weight += this_psp_weight;
  step_sum_weight2 += this_psp_weight2;
  
  if (switch_CV != 0){
    for (int c = 0; c < n_qTcut; c++){
      for (int s = 0; s < n_scales_CV; s++){
	if (c <= cut_ps[0]){
	  this_psp_weight_CV[c][s] = var_rel_alpha_S_CV[s] * psi_ps_factor * VA_b_ME2 * QT_virt_CV[s];
	  this_psp_weight2_CV[c][s] = pow(this_psp_weight_CV[c][s], 2.);
	  step_sum_weight_CV[c][s] += this_psp_weight_CV[c][s];
	  step_sum_weight2_CV[c][s] += this_psp_weight2_CV[c][s];
	}
      }
    }
  }
  
  if (n_moments != 0){
    for (int i_m = 0; i_m < moment.size(); i_m++){
      for (int ia = 0; ia < moment[i_m].size(); ia++){
	this_psp_moment[i_m] = this_psp_weight * moment[i_m][ia];
	this_psp_moment2[i_m] = pow(this_psp_moment[i_m], 2.);
	step_sum_moment[i_m] += this_psp_moment[i_m];
	step_sum_moment2[i_m] += this_psp_moment2[i_m];
	if (switch_CV != 0){
	  for (int c = 0; c < n_qTcut; c++){
	    for (int s = 0; s < n_scales_CV; s++){
	      if (c <= cut_ps[0]){
		this_psp_moment_CV[i_m][c][s] = this_psp_weight_CV[c][s] * moment[i_m][ia];
		this_psp_moment2_CV[i_m][c][s] = pow(this_psp_moment_CV[i_m][c][s], 2.);
		step_sum_moment_CV[i_m][c][s] += this_psp_moment_CV[i_m][c][s];
		step_sum_moment2_CV[i_m][c][s] += this_psp_moment2_CV[i_m][c][s];
	      }
	    }
	  }
	}
      }
    }
  }
  
  //  exit(0);

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
