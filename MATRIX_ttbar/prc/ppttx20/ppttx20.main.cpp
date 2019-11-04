#include "ppttx20.header.cxx"
#include <header/pptt40.amplitude.doublevirtual.h>
int main(int argc, char *argv[]){
  cout << "BEGIN" << endl;

  munich munich(argc, argv, "pp-tt~+X");

  if (munich.subprocess != ""){
    // generic functions - relevant for all subcontributions
    munich.generic.particles = & ppttx20_particles;
    munich.generic.cuts = & ppttx20_cuts;
    munich.generic.cuts_test = & ppttx20_cuts_test;
    munich.generic.calculate_dynamic_scale = & ppttx20_calculate_dynamic_scale;
    munich.generic.calculate_dynamic_scale_TSV = & ppttx20_calculate_dynamic_scale_TSV;
    munich.generic.moments = & ppttx20_moments;
    if (munich.csi.type_contribution == "born" || 
        munich.csi.type_contribution == "VA" || 
        munich.csi.type_contribution == "CA" || 
        munich.csi.type_contribution == "VT" || 
        munich.csi.type_contribution == "CT" || 
        munich.csi.type_contribution == "VT2" || 
        munich.csi.type_contribution == "CT2" || 
        munich.csi.type_contribution == "NLL_LO" || 
        munich.csi.type_contribution == "NLL_NLO" || 
        munich.csi.type_contribution == "NNLL_LO" || 
        munich.csi.type_contribution == "NNLL_NLO" || 
        munich.csi.type_contribution == "NNLL_NNLO"){
      // generic functions - relevant for all subcontributions located on born phase-space in order to determine subprocess information
      munich.generic.determination_subprocess_psp = & ppttx20_determination_subprocess_born;
      munich.generic.combination_subprocess_psp = & ppttx20_combination_subprocess_born;
      // generic functions - relevant for all subcontributions located on born phase-space in order to perform phase-space generation
      munich.generic.optimize_minv_psp = & ppttx20_optimize_minv_born;
      munich.generic.determination_MCchannels_psp = & ppttx20_determination_MCchannels_born;
      munich.generic.ax_psp_psp = & ppttx20_ax_psp_born;
      munich.generic.ac_tau_psp_psp = & ppttx20_ac_tau_psp_born;
      munich.generic.ac_psp_psp = & ppttx20_ac_psp_born;
      munich.generic.ag_psp_psp = & ppttx20_ag_psp_born;
      if (munich.csi.type_contribution == "CA"){munich.generic.phasespacepoint_psp = & ppttx20_phasespacepoint_collinear;}
      else {munich.generic.phasespacepoint_psp = & ppttx20_phasespacepoint_born;}
      if (munich.csi.type_contribution == "VT2" || 
          munich.csi.type_contribution == "VJ2" || 
          munich.csi.type_contribution == "NNLL_NNLO"){
        munich.generic.calculate_H2 = & pptt40_calculate_H2;
      }
    }

    else if (munich.csi.type_contribution == "RA" || 
             munich.csi.type_contribution == "RT" || 
             munich.csi.type_contribution == "RVA" || 
             munich.csi.type_contribution == "RCA"){
      // generic functions - relevant for all subcontributions located on real phase-space in order to determine subprocess information
      munich.generic.determination_subprocess_psp = & ppttx20_determination_subprocess_real;
      munich.generic.combination_subprocess_psp = & ppttx20_combination_subprocess_real;
      // generic functions - relevant for all subcontributions located on real phase-space in order to perform phase-space generation
      munich.generic.determination_MCchannels_psp = & ppttx20_determination_MCchannels_real;
      munich.generic.optimize_minv_psp = & ppttx20_optimize_minv_real;
      munich.generic.ac_tau_psp_psp = & ppttx20_ac_tau_psp_real;
      munich.generic.ax_psp_psp = & ppttx20_ax_psp_real;
      munich.generic.ac_psp_psp = & ppttx20_ac_psp_real;
      munich.generic.ag_psp_psp = & ppttx20_ag_psp_real;
      munich.generic.calculate_dynamic_scale_RA = & ppttx20_calculate_dynamic_scale_RA;
      munich.generic.determination_no_subprocess_dipole = & ppttx20_determination_no_subprocess_born;
      munich.generic.determination_subprocess_dipole = & ppttx20_determination_subprocess_born;
      munich.generic.determination_MCchannels_dipole = & ppttx20_determination_MCchannels_born;
      munich.generic.ac_tau_psp_dipole = & ppttx20_ac_tau_psp_born;
      munich.generic.ax_psp_dipole = & ppttx20_ax_psp_born;
      munich.generic.ac_psp_dipole = & ppttx20_ac_psp_born;
      munich.generic.ag_psp_dipole = & ppttx20_ag_psp_born;
      if (munich.csi.type_contribution == "RCA" || munich.csi.type_contribution == "RCJ"){munich.generic.phasespacepoint_psp = & ppttx20_phasespacepoint_realcollinear;}
      else {munich.generic.phasespacepoint_psp = & ppttx20_phasespacepoint_real;}
    }

    else if (munich.csi.type_contribution == "RRA"){
      // generic functions - relevant for all subcontributions located on doublereal phase-space in order to determine subprocess information
      munich.generic.determination_subprocess_psp = & ppttx20_determination_subprocess_doublereal;
      munich.generic.combination_subprocess_psp = & ppttx20_combination_subprocess_doublereal;
      // generic functions - relevant for all subcontributions located on doublereal phase-space in order to perform phase-space generation
      munich.generic.determination_MCchannels_psp = & ppttx20_determination_MCchannels_doublereal;
      munich.generic.ax_psp_psp = & ppttx20_ax_psp_doublereal;
      munich.generic.ac_tau_psp_psp = & ppttx20_ac_tau_psp_doublereal;
      munich.generic.ac_psp_psp = & ppttx20_ac_psp_doublereal;
      munich.generic.ag_psp_psp = & ppttx20_ag_psp_doublereal;
      munich.generic.optimize_minv_psp = & ppttx20_optimize_minv_doublereal;
      munich.generic.calculate_dynamic_scale_RA = & ppttx20_calculate_dynamic_scale_RA;
      munich.generic.determination_no_subprocess_dipole = & ppttx20_determination_no_subprocess_real;
      munich.generic.determination_subprocess_dipole = & ppttx20_determination_subprocess_real;
      munich.generic.determination_MCchannels_dipole = & ppttx20_determination_MCchannels_real;
      munich.generic.ac_tau_psp_dipole = & ppttx20_ac_tau_psp_real;
      munich.generic.ax_psp_dipole = & ppttx20_ax_psp_real;
      munich.generic.ac_psp_dipole = & ppttx20_ac_psp_real;
      munich.generic.ag_psp_dipole = & ppttx20_ag_psp_real;
      munich.generic.phasespacepoint_psp = & ppttx20_phasespacepoint_doublereal;

      munich.generic.determination_no_subprocess_doubledipole = & ppttx20_determination_no_subprocess_born;
      munich.generic.determination_subprocess_doubledipole = & ppttx20_determination_subprocess_born;
      munich.generic.determination_MCchannels_doubledipole = & ppttx20_determination_MCchannels_born;
      munich.generic.ac_tau_psp_doubledipole = & ppttx20_ac_tau_psp_born;
      munich.generic.ax_psp_doubledipole = & ppttx20_ax_psp_born;
      munich.generic.ac_psp_doubledipole = & ppttx20_ac_psp_born;
      munich.generic.ag_psp_doubledipole = & ppttx20_ag_psp_born;
    }

    munich.run_initialization();
    munich.run_integration();

    cout << "END " << munich.csi.type_contribution << " " << munich.csi.type_correction << endl;
    return 0;
  }
  else {
    munich.generic.list_subprocess_born = & ppttx20_list_subprocess_born;
    munich.generic.list_subprocess_V_QCD = & ppttx20_list_subprocess_V_QCD;
    munich.generic.list_subprocess_C_QCD = & ppttx20_list_subprocess_C_QCD;
    munich.generic.list_subprocess_R_QCD = & ppttx20_list_subprocess_R_QCD;
    munich.generic.list_subprocess_V2_QCD = & ppttx20_list_subprocess_V2_QCD;
    munich.generic.list_subprocess_C2_QCD = & ppttx20_list_subprocess_C2_QCD;
    munich.generic.list_subprocess_RV_QCD = & ppttx20_list_subprocess_RV_QCD;
    munich.generic.list_subprocess_RC_QCD = & ppttx20_list_subprocess_RC_QCD;
    munich.generic.list_subprocess_RR_QCD = & ppttx20_list_subprocess_RR_QCD;

    munich.get_summary();

    cout << "END RESULT/DISTRIBUTION" << endl;
    return 0;
  }
}
