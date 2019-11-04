#include "../include/classes.cxx"

// old version: only needed to restore CT version --- will be removed  soon.
void phasespace_set::initialization_phasespace_IS_CT_from_CX(){
  Logger logger("phasespace_set::initialization_phasespace_IS_CT_from_CX");

  QT_g_z1 = 0.;
  QT_g_z2 = 0.;

  QT_zx_pdf.resize(3, 0.);
  QT_xz_pdf.resize(3, 0.);
  zz_pdf.resize(3, 0.);
  z_pdf.resize(3, 0.);

  QT_random_z1 = QT_random_z[1];
  QT_random_z2 = QT_random_z[2];
  //    QT_random_qt = QT_random_qt;

  logger << LOG_DEBUG << "finished" << endl;
}

// old version: will be replaced by ...CX and ...QT soon.
void phasespace_set::initialization_phasespace_IS_CT(){
  Logger logger("phasespace_set::initialization_phasespace_IS_CT");
  logger << LOG_DEBUG << "started" << endl;

  QT_qt2 = 0.;
  QT_jacqt2 = 0.;
  QT_random_qt2 = 0.;
  QT_random_IS = 0.;
  QT_g_z1 = 0.;
  QT_g_z2 = 0.;
  QT_eps = 0;

  QT_zx_pdf.resize(3, 0.);
  QT_xz_pdf.resize(3, 0.);
  //  zx_pdf.resize(3, 0.);
  //  xz_pdf.resize(3, 0.);
  zz_pdf.resize(3, 0.);
  z_pdf.resize(3, 0.);
  // zall_pdf;

  logger << LOG_DEBUG << "n_z1z2_events = " << n_z1z2_events << endl;
  logger << LOG_DEBUG << "switch_IS_z1z2 = " << switch_IS_z1z2 << endl;
  logger << LOG_DEBUG << "n_z1z2_steps = " << n_z1z2_steps << endl;
  logger << LOG_DEBUG << "n_z1z2_bins = " << n_z1z2_bins << endl;

  QT_random_z1 = new randomvariable(n_z1z2_events, switch_IS_z1z2, n_z1z2_steps, n_z1z2_bins, "z1");
  QT_random_z2 = new randomvariable(n_z1z2_events, switch_IS_z1z2, n_z1z2_steps, n_z1z2_bins, "z2");
  QT_random_qt = new randomvariable(n_z1z2_events, switch_IS_z1z2, n_z1z2_steps, n_z1z2_bins, "qt2");

  random_manager.register_variable(QT_random_z1);
  random_manager.register_variable(QT_random_z2);
  random_manager.register_variable(QT_random_qt);
  //  random_manager.register_variable(QT_random_z1, true);
  //  random_manager.register_variable(QT_random_z2, true);
  //  random_manager.register_variable(QT_random_qt, true);
  
  logger << LOG_DEBUG << "finished" << endl;
}


// old version: will be replaced by ...CX  soon.
void phasespace_set::initialization_phasespace_IS_VT(){
  Logger logger("phasespace_set::initialization_phasespace_IS_VT");
  logger << LOG_DEBUG << "started" << endl;

  QT_qt2 = 0.;
  QT_jacqt2 = 0.;
  QT_random_qt2 = 0.;
  QT_random_IS = 0.;
  QT_g_z1 = 0.;
  QT_g_z2 = 0.;
  QT_eps = 0;

  QT_zx_pdf.resize(3, 0.);
  QT_xz_pdf.resize(3, 0.);
  //  zx_pdf.resize(3, 0.);
  //  xz_pdf.resize(3, 0.);
  zz_pdf.resize(3, 0.);
  z_pdf.resize(3, 0.);
  // zall_pdf;

  logger << LOG_DEBUG << "n_z1z2_events = " << n_z1z2_events << endl;
  logger << LOG_DEBUG << "switch_IS_z1z2 = " << switch_IS_z1z2 << endl;
  logger << LOG_DEBUG << "n_z1z2_steps = " << n_z1z2_steps << endl;
  logger << LOG_DEBUG << "n_z1z2_bins = " << n_z1z2_bins << endl;

  
  QT_random_z1 = new randomvariable(n_z1z2_events, switch_IS_z1z2, n_z1z2_steps, n_z1z2_bins, "z1");
  QT_random_z2 = new randomvariable(n_z1z2_events, switch_IS_z1z2, n_z1z2_steps, n_z1z2_bins, "z2");
  
  random_manager.register_variable(QT_random_z1);
  random_manager.register_variable(QT_random_z2);
  //  random_manager.register_variable(QT_random_z1, true);
  //  random_manager.register_variable(QT_random_z2, true);
  
  //  if (do_resummation) {
  //    QT_random_qt = new randomvariable(n_z1z2_events, switch_IS_z1z2, n_z1z2_steps, n_z1z2_bins, "qt2");
  //    random_manager.register_variable(QT_random_qt, true);
  //  }
 
  logger << LOG_DEBUG << "finished" << endl;
}



