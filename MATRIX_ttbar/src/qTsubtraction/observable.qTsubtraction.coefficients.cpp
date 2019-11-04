#include "../include/classes.cxx"
#include "../include/definitions.phasespace.set.cxx"

#define no_zx ncollinear[ncollinear[i_c].no_endpoint[y_e]].no_pdf
#define no_xx ncollinear[ncollinear[i_c].no_endpoint[0]].no_pdf
#define x_e ncollinear[i_c].x_a
#define y_e ncollinear[i_c].x_b

void observable_set::determine_splitting_tH1F(phasespace_set & psi){
  static Logger logger("observable_set::determine_splitting_tH1F");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tH1F_contribution[i_l].clear();}
  //  tH1F: prefactors of splittings, to be multiplied with pdf factors
  for (int i_c = 1; i_c < ncollinear.size(); i_c++){
    logger << LOG_DEBUG_VERBOSE << "ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << endl;
    if (x_e == 0){logger << LOG_FATAL << "Should not happen at 1st order!!!" << endl; exit(1);}

    // type_splitting_full (leg-wise):
    // x_e - y_e     ---  ---

    //   1 -   0   +   0 -  1   // hard process with g from g -> g (g) splitting
    //   2 -   0   +   0 -  2   // hard process with q from q -> q (g) splitting
    //   3 -   0   +   0 -  3   // hard process with g from q -> g (q) splitting
    //   4 -   0   +   0 -  4   // hard process with q from g -> q (qx) splitting       
    
    if (ncollinear[i_c].type_splitting_full[x_e] == 0){logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = 0 -> no splitting type defined !!!" << endl; exit(1);}
    else if (ncollinear[i_c].type_splitting_full[x_e] == 1){ // hard process with g from g -> g (g) splitting
      coll_tH1F_contribution[no_zx].push_back(3 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] + Pggreg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      // -log( x1 ) * (pdf_factor_z1x2 * Pggreg(z1)
      // -log( x1 ) * (pdf_factor_z1x2)*3/(1-z1)
      coll_tH1F_contribution[no_xx].push_back(-psi_z_coll[x_e] * 3 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] - 3 * D0int(x_pdf[x_e]) + beta0);
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_tH1F_contribution[no_xx = " << no_xx << "] += " << coll_tH1F_contribution[no_xx][coll_tH1F_contribution[no_xx].size() - 1] << endl;
      // -log( x1 ) * (- pdf_factor_x1x2 * z1 * 3 / (1 - z1)) - 3 * D0int(x1) * pdf_factor_x1x2
      // + beta0 * pdf_factor_x1x2 (half of this (because added for each leg)!)
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 2){ // hard process with q from q -> q (g) splitting
      coll_tH1F_contribution[no_zx].push_back(Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      // - g_z1 * (pdf_factor_z1x2) * Pqq(z1)  // from H1F
      coll_tH1F_contribution[no_xx].push_back((-psi_z_coll[x_e] * Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - Pqqint(x_pdf[x_e])));
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1F_contribution[no_xx = " << no_xx << "] += " << coll_tH1F_contribution[no_xx][coll_tH1F_contribution[no_xx].size() - 1] << endl;
      // - g_z1 * (- pdf_factor_x1x2 * z1) * Pqq(z1) - pdf_factor_x1x2 * Pqqint(x1)  // from H1F
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 3){ // hard process with g from q -> g (q) splitting
      coll_tH1F_contribution[no_zx].push_back(Pgq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "q->g(q):  coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      // -log( x1 ) * (pdf_factor_qx2 * Pgq(z1) )
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 4){ // hard process with q from g -> q (qx) splitting
      coll_tH1F_contribution[no_zx].push_back(Pqg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->q(qx): coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      // - g_z1 * Pqg(z1) * (pdf_factor_gx2)  // from H1F
    }
    else {logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << " -> selected splitting type not allowed !!!" << endl; exit(1);}
  }
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    coll_tH1F[i_l] = accumulate(coll_tH1F_contribution[i_l].begin(), coll_tH1F_contribution[i_l].end(), 0.);
    logger << LOG_DEBUG_VERBOSE << "coll_tH1F[" << i_l << "] = " << coll_tH1F[i_l] << endl;
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::determine_splitting_tH1F_tH1(phasespace_set & psi){
  static Logger logger("observable_set::determine_splitting_tH1F_tH1");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tH1F_contribution[i_l].clear();}
  //  tH1F: prefactors of splittings, to be multiplied with pdf factors
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tH1_contribution[i_l].clear();}
  //  tH1: prefactors of splittings, to be multiplied with pdf factors

  for (int i_c = 1; i_c < ncollinear.size(); i_c++){
    logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << "   x_e = " << x_e << endl;
    if (x_e == 0){logger << LOG_DEBUG_VERBOSE << "ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << " does not contribute here." << endl; continue;}
    logger << LOG_DEBUG_POINT << "ncollinear[" << i_c << "]: " << "type_splitting_full[x_e = " << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << "   type_splitting_full[y_e = " << y_e << "] = " << ncollinear[i_c].type_splitting_full[y_e] << endl;
    
    // type_splitting_full (leg-wise):
    // x_e - y_e     ---  ---

    //  14 -   0   +   0 - 14   // hard process with q from g -> g (g)  -> q (qx) splitting (-> 4  in sig11 calculation)
    // (42 -   0   +   0 - 42   // hard process with q from g -> q (q)  -> q (g)  splitting (part of  14 -  0   +   0 - 14))
    //  11 -   0   +   0 - 11   // hard process with g from g -> g (g)  -> g (g)  splitting (-> 1  in sig11 calculation)
    // (43 -   0   +   0 - 43   // hard process with g from g -> q (qx) -> g (Q)  splitting (part of  11 -  0   +   0 - 11))
    //  23 -   0   +   0 - 23   // hard process with g from q -> q (g)  -> g (Q)  splitting (-> 3  in sig11 calculation)
    // (31 -   0   +   0 - 31   // hard process with g from q -> g (q)  -> g (g)  splitting (part of  23 -  0   +   0 - 23))
    //  22 -   0   +   0 - 22   // hard process with q from q -> q (g)  -> q (g)  splitting (-> 2  in sig11 calculation)

    //  34 -   0   +   0 - 34   // hard process with q from q -> g (q)  -> q (qx) splitting (----  no equivalent)
    //  54 -   0   +   0 - 54   // hard process with q from qx -> ...   -> q      splitting (non-singlet)    (----  no equivalent)
    //   2 -   2                // hard process with (leg 1) q from q -> q (g)  splitting and (leg 2) q from q -> q (g)  splitting
    //   2 -   4   +   4 -  2   // hard process with (leg 1) q from q -> q (g)  splitting and (leg 2) g from q -> g (q)  splitting
    //   4 -   4                // hard process with (leg 1) g from q -> g (q)  splitting and (leg 2) g from q -> g (q)  splitting
    //   1 -   1                // hard process with (leg 1) g from g -> g (g)  splitting and (leg 2) g from g -> g (g)  splitting
    //   1 -   3   +   3 -  1   // hard process with (leg 1) g from g -> g (g)  splitting and (leg 2) g from q -> g (q)  splitting
    //   3 -   3                // hard process with (leg 1) g from q -> g (q)  splitting and (leg 2) g from q -> g (q)  splitting

    if (ncollinear[i_c].type_splitting_full[x_e] == 0){logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = 0 -> no splitting type defined !!!" << endl; exit(1);}
    else if ((ncollinear[i_c].type_splitting_full[x_e] == 11 && ncollinear[i_c].type_splitting_full[y_e] == 0) ||
	     (ncollinear[i_c].type_splitting_full[x_e] == 1 && ncollinear[i_c].type_splitting_full[y_e] == 0)){ // hard process with q from g -> g (g) splitting
      coll_tH1F_contribution[no_zx].push_back(3 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] + Pggreg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      coll_tH1F_contribution[no_xx].push_back(-psi_z_coll[x_e] * 3 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] - 3 * D0int(x_pdf[x_e]) + beta0);
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_tH1F_contribution[no_xx = " << no_xx << "] += " << coll_tH1F_contribution[no_xx][coll_tH1F_contribution[no_xx].size() - 1] << endl;
      coll_tH1_contribution[no_xx].push_back(.5 * QT_H1_delta);
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_tH1_contribution[no_xx = " << no_xx << "] += " << coll_tH1_contribution[no_xx][coll_tH1_contribution[no_xx].size() - 1] << endl;
    }
    else if ((ncollinear[i_c].type_splitting_full[x_e] == 23 && ncollinear[i_c].type_splitting_full[y_e] == 0) ||
	     (ncollinear[i_c].type_splitting_full[x_e] == 3 && ncollinear[i_c].type_splitting_full[y_e] == 0)){ // hard process with g from q -> g (q) splitting
      coll_tH1F_contribution[no_zx].push_back(Pgq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "q->g(q):  coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      coll_tH1_contribution[no_zx].push_back(Cgq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "q->g(q):  coll_tH1_contribution[no_zx = " << no_zx << "] += " << coll_tH1_contribution[no_zx][coll_tH1_contribution[no_zx].size() - 1] << endl;
    }
    else if ((ncollinear[i_c].type_splitting_full[x_e] == 22 && ncollinear[i_c].type_splitting_full[y_e] == 0) ||
	     (ncollinear[i_c].type_splitting_full[x_e] == 2 && ncollinear[i_c].type_splitting_full[y_e] == 0)){ // hard process with q from q -> q (g) splitting
      coll_tH1F_contribution[no_zx].push_back(Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      coll_tH1F_contribution[no_xx].push_back((-psi_z_coll[x_e] * Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - Pqqint(x_pdf[x_e])));
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1F_contribution[no_xx = " << no_xx << "] += " << coll_tH1F_contribution[no_xx][coll_tH1F_contribution[no_xx].size() - 1] << endl;
      coll_tH1_contribution[no_zx].push_back(Cqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1_contribution[no_zx = " << no_zx << "] += " << coll_tH1_contribution[no_zx][coll_tH1_contribution[no_zx].size() - 1] << endl;
      coll_tH1_contribution[no_xx].push_back(.5 * QT_H1_delta);
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1_contribution[no_xx = " << no_xx << "] += " << coll_tH1_contribution[no_xx][coll_tH1_contribution[no_xx].size() - 1] << endl;
    }
    else if ((ncollinear[i_c].type_splitting_full[x_e] == 14 && ncollinear[i_c].type_splitting_full[y_e] == 0) ||
	     (ncollinear[i_c].type_splitting_full[x_e] == 4 && ncollinear[i_c].type_splitting_full[y_e] == 0)){ // hard process with q from g -> q (qx) splitting
      coll_tH1F_contribution[no_zx].push_back(Pqg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->q(qx): coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      coll_tH1_contribution[no_zx].push_back(Cqg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->q(qx): coll_tH1_contribution[no_zx = " << no_zx << "] += " << coll_tH1_contribution[no_zx][coll_tH1_contribution[no_zx].size() - 1] << endl;
    }
    else {
      //  logger << LOG_DEBUG_VERBOSE << "ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << " does not contribute here." << endl; 
      continue;
    }
  }
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    coll_tH1F[i_l] = accumulate(coll_tH1F_contribution[i_l].begin(), coll_tH1F_contribution[i_l].end(), 0.);
    coll_tH1[i_l] = accumulate(coll_tH1_contribution[i_l].begin(), coll_tH1_contribution[i_l].end(), 0.);
    logger << LOG_DEBUG_VERBOSE << "coll_tH1F[" << i_l << "] = " << coll_tH1F[i_l] << endl;
    logger << LOG_DEBUG_VERBOSE << "coll_tH1[" << i_l << "] = " << coll_tH1[i_l] << endl;
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_splitting_tH1F_tH1_without_H1_delta(phasespace_set & psi){
  static Logger logger("observable_set::determine_splitting_tH1F_tH1_without_H1_delta");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tH1F_contribution[i_l].clear();}
  //  tH1F: prefactors of splittings, to be multiplied with pdf factors
  //  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tH1_without_H1_delta_contribution[i_l].clear();}
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tH1_contribution[i_l].clear();}
  //  tH1: prefactors of splittings, to be multiplied with pdf factors

  for (int i_c = 1; i_c < ncollinear.size(); i_c++){
    logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << "   x_e = " << x_e << endl;
    if (x_e == 0){logger << LOG_DEBUG_VERBOSE << "ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << " does not contribute here." << endl; continue;}

    // type_splitting_full (leg-wise):
    // x_e - y_e     ---  ---

    //  14 -   0   +   0 - 14   // hard process with q from g -> g (g)  -> q (qx) splitting (-> 4  in sig11 calculation)
    // (42 -   0   +   0 - 42   // hard process with q from g -> q (q)  -> q (g)  splitting (part of  14 -  0   +   0 - 14))
    //  11 -   0   +   0 - 11   // hard process with g from g -> g (g)  -> g (g)  splitting (-> 1  in sig11 calculation)
    // (43 -   0   +   0 - 43   // hard process with g from g -> q (qx) -> g (Q)  splitting (part of  11 -  0   +   0 - 11))
    //  23 -   0   +   0 - 23   // hard process with g from q -> q (g)  -> g (Q)  splitting (-> 3  in sig11 calculation)
    // (31 -   0   +   0 - 31   // hard process with g from q -> g (q)  -> g (g)  splitting (part of  23 -  0   +   0 - 23))
    //  22 -   0   +   0 - 22   // hard process with q from q -> q (g)  -> q (g)  splitting (-> 2  in sig11 calculation)

    //  34 -   0   +   0 - 34   // hard process with q from q -> g (q)  -> q (qx) splitting (----  no equivalent)
    //  54 -   0   +   0 - 54   // hard process with q from qx -> ...   -> q      splitting (non-singlet)    (----  no equivalent)
    //   2 -   2                // hard process with (leg 1) q from q -> q (g)  splitting and (leg 2) q from q -> q (g)  splitting
    //   2 -   4   +   4 -  2   // hard process with (leg 1) q from q -> q (g)  splitting and (leg 2) g from q -> g (q)  splitting
    //   4 -   4                // hard process with (leg 1) g from q -> g (q)  splitting and (leg 2) g from q -> g (q)  splitting
    //   1 -   1                // hard process with (leg 1) g from g -> g (g)  splitting and (leg 2) g from g -> g (g)  splitting
    //   1 -   3   +   3 -  1   // hard process with (leg 1) g from g -> g (g)  splitting and (leg 2) g from q -> g (q)  splitting
    //   3 -   3                // hard process with (leg 1) g from q -> g (q)  splitting and (leg 2) g from q -> g (q)  splitting

    if (ncollinear[i_c].type_splitting_full[x_e] == 0){logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = 0 -> no splitting type defined !!!" << endl; exit(1);}
    else if ((ncollinear[i_c].type_splitting_full[x_e] == 11 && ncollinear[i_c].type_splitting_full[y_e] == 0) ||
	     (ncollinear[i_c].type_splitting_full[x_e] == 1 && ncollinear[i_c].type_splitting_full[y_e] == 0)){ // hard process with q from g -> g (g) splitting
      coll_tH1F_contribution[no_zx].push_back(3 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] + Pggreg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      coll_tH1F_contribution[no_xx].push_back(-psi_z_coll[x_e] * 3 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] - 3 * D0int(x_pdf[x_e]) + beta0);
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_tH1F_contribution[no_xx = " << no_xx << "] += " << coll_tH1F_contribution[no_xx][coll_tH1F_contribution[no_xx].size() - 1] << endl;
      ///      coll_tH1_contribution[no_xx].push_back(.5 * QT_H1_delta);
      //      coll_tH1_without_H1_delta_contribution[no_xx].push_back(0.); // to be removed later !!!
      coll_tH1_contribution[no_xx].push_back(0.); // to be removed later !!!
      //      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_tH1_without_H1_delta_contribution[no_xx = " << no_xx << "] += " << coll_tH1_without_H1_delta_contribution[no_xx][coll_tH1_without_H1_delta_contribution[no_xx].size() - 1] << endl;
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_tH1_contribution[no_xx = " << no_xx << "] += " << coll_tH1_contribution[no_xx][coll_tH1_contribution[no_xx].size() - 1] << endl;
    }
    else if ((ncollinear[i_c].type_splitting_full[x_e] == 23 && ncollinear[i_c].type_splitting_full[y_e] == 0) ||
	     (ncollinear[i_c].type_splitting_full[x_e] == 3 && ncollinear[i_c].type_splitting_full[y_e] == 0)){ // hard process with g from q -> g (q) splitting
      coll_tH1F_contribution[no_zx].push_back(Pgq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "q->g(q):  coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      //      coll_tH1_without_H1_delta_contribution[no_zx].push_back(Cgq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tH1_contribution[no_zx].push_back(Cgq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      //      logger << LOG_DEBUG_VERBOSE << "q->g(q):  coll_tH1_without_H1_delta_contribution[no_zx = " << no_zx << "] += " << coll_tH1_without_H1_delta_contribution[no_zx][coll_tH1_without_H1_delta_contribution[no_zx].size() - 1] << endl;
      logger << LOG_DEBUG_VERBOSE << "q->g(q):  coll_tH1_contribution[no_zx = " << no_zx << "] += " << coll_tH1_contribution[no_zx][coll_tH1_contribution[no_zx].size() - 1] << endl;
    }
    else if ((ncollinear[i_c].type_splitting_full[x_e] == 22 && ncollinear[i_c].type_splitting_full[y_e] == 0) ||
	     (ncollinear[i_c].type_splitting_full[x_e] == 2 && ncollinear[i_c].type_splitting_full[y_e] == 0)){ // hard process with q from q -> q (g) splitting
      coll_tH1F_contribution[no_zx].push_back(Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      coll_tH1F_contribution[no_xx].push_back((-psi_z_coll[x_e] * Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e] - Pqqint(x_pdf[x_e])));
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1F_contribution[no_xx = " << no_xx << "] += " << coll_tH1F_contribution[no_xx][coll_tH1F_contribution[no_xx].size() - 1] << endl;
      //      coll_tH1_without_H1_delta_contribution[no_zx].push_back(Cqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tH1_contribution[no_zx].push_back(Cqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      //      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1_without_H1_delta_contribution[no_zx = " << no_zx << "] += " << coll_tH1_without_H1_delta_contribution[no_zx][coll_tH1_without_H1_delta_contribution[no_zx].size() - 1] << endl;
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1_contribution[no_zx = " << no_zx << "] += " << coll_tH1_contribution[no_zx][coll_tH1_contribution[no_zx].size() - 1] << endl;
      ///      coll_tH1_contribution[no_xx].push_back(.5 * QT_H1_delta);
      //      coll_tH1_without_H1_delta_contribution[no_xx].push_back(0.); // to be removed later !!!
      coll_tH1_contribution[no_xx].push_back(0.); // to be removed later !!!
      //      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1_without_H1_delta_contribution[no_xx = " << no_xx << "] += " << coll_tH1_without_H1_delta_contribution[no_xx][coll_tH1_without_H1_delta_contribution[no_xx].size() - 1] << endl;
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1_contribution[no_xx = " << no_xx << "] += " << coll_tH1_contribution[no_xx][coll_tH1_contribution[no_xx].size() - 1] << endl;
    }
    else if ((ncollinear[i_c].type_splitting_full[x_e] == 14 && ncollinear[i_c].type_splitting_full[y_e] == 0) ||
	     (ncollinear[i_c].type_splitting_full[x_e] == 4 && ncollinear[i_c].type_splitting_full[y_e] == 0)){ // hard process with q from g -> q (qx) splitting
      coll_tH1F_contribution[no_zx].push_back(Pqg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->q(qx): coll_tH1F_contribution[no_zx = " << no_zx << "] += " << coll_tH1F_contribution[no_zx][coll_tH1F_contribution[no_zx].size() - 1] << endl;
      //      coll_tH1_without_H1_delta_contribution[no_zx].push_back(Cqg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tH1_contribution[no_zx].push_back(Cqg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      //      logger << LOG_DEBUG_VERBOSE << "g->q(qx): coll_tH1_without_H1_delta_contribution[no_zx = " << no_zx << "] += " << coll_tH1_without_H1_delta_contribution[no_zx][coll_tH1_without_H1_delta_contribution[no_zx].size() - 1] << endl;
      logger << LOG_DEBUG_VERBOSE << "g->q(qx): coll_tH1_contribution[no_zx = " << no_zx << "] += " << coll_tH1_contribution[no_zx][coll_tH1_contribution[no_zx].size() - 1] << endl;
    }
    else {
      //  logger << LOG_DEBUG_VERBOSE << "ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << " does not contribute here." << endl; 
      continue;
    }
  }
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    coll_tH1F[i_l] = accumulate(coll_tH1F_contribution[i_l].begin(), coll_tH1F_contribution[i_l].end(), 0.);
    //    coll_tH1[i_l] = accumulate(coll_tH1_without_H1_delta_contribution[i_l].begin(), coll_tH1_without_H1_delta_contribution[i_l].end(), 0.);
    coll_tH1[i_l] = accumulate(coll_tH1_contribution[i_l].begin(), coll_tH1_contribution[i_l].end(), 0.);
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}



void observable_set::determine_splitting_tH1_only_H1_delta(phasespace_set & psi){
  static Logger logger("observable_set::determine_splitting_tH1_only_H1_delta");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tH1_only_H1_delta_contribution[i_l].clear();}
  //  tH1: prefactors of splittings, to be multiplied with pdf factors

  for (int i_c = 1; i_c < ncollinear.size(); i_c++){
    logger << LOG_DEBUG_VERBOSE << "i_c = " << i_c << "   x_e = " << x_e << endl;
    if (x_e == 0){logger << LOG_DEBUG_VERBOSE << "ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << " does not contribute here." << endl; continue;}

    // type_splitting_full (leg-wise):
    // x_e - y_e     ---  ---

    //  14 -   0   +   0 - 14   // hard process with q from g -> g (g)  -> q (qx) splitting (-> 4  in sig11 calculation)
    // (42 -   0   +   0 - 42   // hard process with q from g -> q (q)  -> q (g)  splitting (part of  14 -  0   +   0 - 14))
    //  11 -   0   +   0 - 11   // hard process with g from g -> g (g)  -> g (g)  splitting (-> 1  in sig11 calculation)
    // (43 -   0   +   0 - 43   // hard process with g from g -> q (qx) -> g (Q)  splitting (part of  11 -  0   +   0 - 11))
    //  23 -   0   +   0 - 23   // hard process with g from q -> q (g)  -> g (Q)  splitting (-> 3  in sig11 calculation)
    // (31 -   0   +   0 - 31   // hard process with g from q -> g (q)  -> g (g)  splitting (part of  23 -  0   +   0 - 23))
    //  22 -   0   +   0 - 22   // hard process with q from q -> q (g)  -> q (g)  splitting (-> 2  in sig11 calculation)

    //  34 -   0   +   0 - 34   // hard process with q from q -> g (q)  -> q (qx) splitting (----  no equivalent)
    //  54 -   0   +   0 - 54   // hard process with q from qx -> ...   -> q      splitting (non-singlet)    (----  no equivalent)
    //   2 -   2                // hard process with (leg 1) q from q -> q (g)  splitting and (leg 2) q from q -> q (g)  splitting
    //   2 -   4   +   4 -  2   // hard process with (leg 1) q from q -> q (g)  splitting and (leg 2) g from q -> g (q)  splitting
    //   4 -   4                // hard process with (leg 1) g from q -> g (q)  splitting and (leg 2) g from q -> g (q)  splitting
    //   1 -   1                // hard process with (leg 1) g from g -> g (g)  splitting and (leg 2) g from g -> g (g)  splitting
    //   1 -   3   +   3 -  1   // hard process with (leg 1) g from g -> g (g)  splitting and (leg 2) g from q -> g (q)  splitting
    //   3 -   3                // hard process with (leg 1) g from q -> g (q)  splitting and (leg 2) g from q -> g (q)  splitting

    if (ncollinear[i_c].type_splitting_full[x_e] == 0){logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = 0 -> no splitting type defined !!!" << endl; exit(1);}
    else if ((ncollinear[i_c].type_splitting_full[x_e] == 11 && ncollinear[i_c].type_splitting_full[y_e] == 0) ||
	     (ncollinear[i_c].type_splitting_full[x_e] == 1 && ncollinear[i_c].type_splitting_full[y_e] == 0)){ // hard process with q from g -> g (g) splitting
      coll_tH1_only_H1_delta_contribution[no_xx].push_back(.5 * QT_H1_delta);
      logger << LOG_DEBUG_VERBOSE << "g->g(g):  coll_tH1_only_H1_delta_contribution[no_xx = " << no_xx << "] += " << coll_tH1_only_H1_delta_contribution[no_xx][coll_tH1_only_H1_delta_contribution[no_xx].size() - 1] << endl;
    }
    else if ((ncollinear[i_c].type_splitting_full[x_e] == 22 && ncollinear[i_c].type_splitting_full[y_e] == 0) ||
	     (ncollinear[i_c].type_splitting_full[x_e] == 2 && ncollinear[i_c].type_splitting_full[y_e] == 0)){ // hard process with q from q -> q (g) splitting
      coll_tH1_only_H1_delta_contribution[no_xx].push_back(.5 * QT_H1_delta);
      logger << LOG_DEBUG_VERBOSE << "q->q(g):  coll_tH1_only_H1_delta_contribution[no_xx = " << no_xx << "] += " << coll_tH1_only_H1_delta_contribution[no_xx][coll_tH1_only_H1_delta_contribution[no_xx].size() - 1] << endl;
    }
    else {
      //  logger << LOG_DEBUG_VERBOSE << "ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << " does not contribute here." << endl; 
      continue;
    }
  }
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    coll_tH1_only_H1_delta[i_l] = accumulate(coll_tH1_only_H1_delta_contribution[i_l].begin(), coll_tH1_only_H1_delta_contribution[i_l].end(), 0.);
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

#undef no_zx 
#undef no_xx
#undef x_e
#undef y_e



#define no_zz ncollinear[ncollinear[i_c].no_endpoint[3]].no_pdf
#define no_xz ncollinear[ncollinear[i_c].no_endpoint[x_e]].no_pdf
#define no_zx ncollinear[ncollinear[i_c].no_endpoint[y_e]].no_pdf
#define no_xx ncollinear[ncollinear[i_c].no_endpoint[0]].no_pdf
#define x_e ncollinear[i_c].x_a
#define y_e ncollinear[i_c].x_b
void observable_set::determine_splitting_tgaga_tcga_tgamma2(phasespace_set & psi){
  static Logger logger("observable_set::determine_splitting_tgaga_tcga_tgamma2");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tgaga_contribution[i_l].clear();}
  //  tgaga: prefactors of splittings, to be multiplied with pdf factors
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tcga_contribution[i_l].clear();}
  //  tcga: prefactors of splittings, to be multiplied with pdf factors
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tgamma2_contribution[i_l].clear();}
  //  tgamma2: prefactors of splittings, to be multiplied with pdf factors

  for (int i_c = 1; i_c < ncollinear.size(); i_c++){

    // type_splitting_full (leg-wise):
    // x_e - y_e     ---  ---

    //  14 -   0   +   0 - 14   // hard process with q from g -> g (g)  -> q (qx) splitting (-> 4  in sig11 calculation)
    // (42 -   0   +   0 - 42   // hard process with q from g -> q (q)  -> q (g)  splitting (part of  14 -  0   +   0 - 14))
    //  11 -   0   +   0 - 11   // hard process with g from g -> g (g)  -> g (g)  splitting (-> 1  in sig11 calculation)
    // (43 -   0   +   0 - 43   // hard process with g from g -> q (qx) -> g (Q)  splitting (part of  11 -  0   +   0 - 11))
    //  23 -   0   +   0 - 23   // hard process with g from q -> q (g)  -> g (Q)  splitting (-> 3  in sig11 calculation)
    // (31 -   0   +   0 - 31   // hard process with g from q -> g (q)  -> g (g)  splitting (part of  23 -  0   +   0 - 23))
    //  22 -   0   +   0 - 22   // hard process with q from q -> q (g)  -> q (g)  splitting (-> 2  in sig11 calculation)

    //  34 -   0   +   0 - 34   // hard process with q from q -> g (q)  -> q (qx) splitting (----  no equivalent)
    //  54 -   0   +   0 - 54   // hard process with q from qx -> ...   -> q      splitting (non-singlet)    (----  no equivalent)
    //   2 -   2                // hard process with (leg 1) q from q -> q (g)  splitting and (leg 2) q from q -> q (g)  splitting
    //   2 -   4   +   4 -  2   // hard process with (leg 1) q from q -> q (g)  splitting and (leg 2) g from q -> g (q)  splitting
    //   4 -   4                // hard process with (leg 1) g from q -> g (q)  splitting and (leg 2) g from q -> g (q)  splitting
    //   1 -   1                // hard process with (leg 1) g from g -> g (g)  splitting and (leg 2) g from g -> g (g)  splitting
    //   1 -   3   +   3 -  1   // hard process with (leg 1) g from g -> g (g)  splitting and (leg 2) g from q -> g (q)  splitting
    //   3 -   3                // hard process with (leg 1) g from q -> g (q)  splitting and (leg 2) g from q -> g (q)  splitting

    if (ncollinear[i_c].type_splitting_full[x_e] == 0){logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = 0 -> no splitting type defined !!!" << endl; exit(1);}


    //    tgaga
    // 11 - 0
    //   -  log( x1 )*( (pdf_factor_z1x2-pdf_factor_x1x2* z1 )*(D0gggg/(1-z1)+D1gggg* log(1-z1)/(1-z1)) + pdf_factor_z1x2*(Pggggreg(z1, beta0)+Pgqqg(z1, nf)) ) + (Deltagggg-D0gggg*D0int( x1 )-D1gggg*D1int( x1 ))*pdf_factor_x1x2;
    // second leg (automatically included)   -  log( x2 )*( (pdf_factor_x1z2-pdf_factor_x1x2* z2 )*(D0gggg/(1-z2)+D1gggg* log(1-z2)/(1-z2)) + pdf_factor_x1z2*(Pggggreg(z2, beta0)+Pgqqg(z2, nf)) ) + (Deltagggg-D0gggg*D0int( x2 )-D1gggg*D1int( x2 ))*pdf_factor_x1x2;
    // 23 - 0
    //   -  log( x1 ) * (Pgqqq(z1)+Pgggq(z1, beta0)) * pdf_factor_qx2;
    // second leg (automatically included)   -  log( x2 ) * (Pgqqq(z2)+Pgggq(z2, beta0)) * pdf_factor_x1q;
    //   double tx1 = 3.0*log(x1)*z1/(1-z1) - 3*D0int(x1) + beta0,  tz1 = -log(x1)*( 3.0/(1-z1) + Pggreg(z1) ) ;
    //   double tx2 = 3.0*log(x2)*z2/(1-z2) - 3*D0int(x2) + beta0,  tz2 = -log(x2)*( 3.0/(1-z2) + Pggreg(z2) ) ;
    // 1 - 1
    //   2 * (pdf_factor_x1x2*tx1*tx2 + pdf_factor_x1z2*tx1*tz2 + pdf_factor_z1x2*tz1*tx2 + pdf_factor_z1z2*tz1*tz2);
    // 1 - 3
    //   -2 * ( pdf_factor_x1q*tx1 + pdf_factor_z1q*tz1 )*log(x2)*Pgq(z2);
    // second leg (automatically included)      -2 * log(x1)*Pgq(z1)*( pdf_factor_qx2*tx2 + pdf_factor_qz2*tz2 );
    // 3 - 3
    //   2 * log(x1)*Pgq(z1)*log(x2)*Pgq(z2)*pdf_factor_qq;
    

    //   tcga
    // 11 - 0
    //   tcga = tcga - CgqPqg(z1, nf)* log( x1 )     * pdf_factor_z1x2;
    // second leg (automatically included)   tcga = tcga - CgqPqg(z2, nf)* log( x2 )     * pdf_factor_x1z2;
    // 23 - 0
    //   tcga = tcga - CgqPqq(z1)* log( x1 )*pdf_factor_qx2;
    // second leg (automatically included)   tcga = tcga - CgqPqq(z2)* log( x2 )*pdf_factor_x1q;
    // 11 - 0
    //   double diffg10_diffc20 = ( pdf_factor_z1x2*tz1 + pdf_factor_x1x2*tx1 )*C1ggdelta(A_F) ;
    // 3 - 3
    //   double diffg1f_diffc2f = pdf_factor_qq * log(x1)*log(x2)*Pgq(z1)*Cgq(z2) ;
    // 1 - 3
    //   double diffg10_diffc2f = -log(x2)*Cgq(z2)*( pdf_factor_z1q*tz1 + pdf_factor_x1q*tx1 ) ;
    // 23 - 0
    //   double diffg1f_diffc20 = -log(x1)*Pgq(z1)*C1ggdelta(A_F)*pdf_factor_qx2;
    // 0 - 11
    // second leg (automatically included)   double diffc10_diffg20 = C1ggdelta(A_F)*( pdf_factor_x1z2*tz2 + pdf_factor_x1x2*tx2 ) ;
    // 3 - 3
    //   double diffc1f_diffg2f = pdf_factor_qq*log(x1)*log(x2)*Cgq(z1)*Pgq(z2) ;
    // 3 - 1
    // second leg (automatically included)   double diffc1f_diffg20 = -log(x1)*Cgq(z1)*(pdf_factor_qz2*tz2 + pdf_factor_qx2*tx2) ;
    // 0 - 23 
    // second leg (automatically included)   double diffc10_diffg2f = -C1ggdelta(A_F)*log(x2)*Pgq(z2)*pdf_factor_x1q ;
    //   tcga = tcga + diffg10_diffc20 + diffg1f_diffc2f + diffg10_diffc2f + diffg1f_diffc20;
    //   tcga = tcga + diffc10_diffg20 + diffc1f_diffg2f + diffc1f_diffg20 + diffc10_diffg2f;


    //   tgamma2
    //  gamma2: diagonal part
    //   First leg
    // 11 - 0
    //    tgamma2 = tgamma2 - ( 1.5  *Kappa*(  log( x1 )*(pdf_factor_z1x2-pdf_factor_x1x2* z1       )/(1-z1) + D0int( x1 )*pdf_factor_x1x2) +  log( x1 )*P2gg(z1, nf)*pdf_factor_z1x2 );  
    //   Second leg
    // second leg (automatically included)   tgamma2 = tgamma2 - ( 1.5  *Kappa*(  log( x2 )*(pdf_factor_x1z2-pdf_factor_x1x2* z2       )/(1-z2) + D0int( x2 )*pdf_factor_x1x2) +  log( x2 )*P2gg(z2, nf)*pdf_factor_x1z2 );
    //  gamma2: qg channel
    //   First leg
    // 23 - 0
    //   tgamma2 = tgamma2 -  log( x1 )*P2gq(z1, nf)* pdf_factor_qx2;
    //  Second leg
    // second leg (automatically included)   tgamma2 = tgamma2 -  log( x2 )*P2gq(z2, nf)* pdf_factor_x1q;


    else if (ncollinear[i_c].type_splitting_full[x_e] == 11 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with g from g -> g (g) -> g (g)  or  g -> Q (Qx) -> g (Q) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with g from g -> g (g) -> g (g)  or  g -> Q (Qx) -> g (Q) splitting" << endl;
      double temp_Dxgggg = (D0gggg / (1. - psi_z_coll[x_e]) + D1gggg * log(1. - psi_z_coll[x_e]) / (1. - psi_z_coll[x_e]));
      coll_tgaga_contribution[no_zx].push_back(temp_Dxgggg / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tgaga  [no_zx = " << setw(2) << no_zx << "] += " << coll_tgaga_contribution[no_zx][coll_tgaga_contribution[no_zx].size() - 1] << endl;
      coll_tgaga_contribution[no_zx].push_back((Pggggreg(psi_z_coll[x_e], beta0) + Pgqqg(psi_z_coll[x_e], N_f)) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(gg/qqx):  tgaga  [no_zx = " << setw(2) << no_zx << "] += " << coll_tgaga_contribution[no_zx][coll_tgaga_contribution[no_zx].size() - 1] << endl;
      coll_tgaga_contribution[no_xx].push_back(-psi_z_coll[x_e] * temp_Dxgggg / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tgaga  [no_xx = " << setw(2) << no_xx << "] += " << coll_tgaga_contribution[no_xx][coll_tgaga_contribution[no_xx].size() - 1] << endl;
      coll_tgaga_contribution[no_xx].push_back(Deltagggg - D0gggg * D0int(psi_x_pdf[x_e]) - D1gggg * D1int(psi_x_pdf[x_e]));
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tgaga  [no_xx = " << setw(2) << no_xx << "] += " << coll_tgaga_contribution[no_xx][coll_tgaga_contribution[no_xx].size() - 1] << endl;

      coll_tcga_contribution[no_zx].push_back(CgqPqg(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(qqx):     tcga   [no_zx = " << setw(2) << no_zx << "] += " << coll_tcga_contribution[no_zx][coll_tcga_contribution[no_zx].size() - 1] << endl;

      double temp_C1ggdelta_A_F = C1ggdelta(A_F);
      double tx1 = -3 * psi_z_coll[x_e] / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] - 3 * D0int(psi_x_pdf[x_e]) + beta0;
      double tz1 = (3 / (1. - psi_z_coll[x_e]) + Pggreg(psi_z_coll[x_e])) / psi_g_z_coll[x_e];
      coll_tcga_contribution[no_zx].push_back(tz1 * temp_C1ggdelta_A_F); // scale dependent
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tcga   [no_zx = " << setw(2) << no_zx << "] += " << coll_tcga_contribution[no_zx][coll_tcga_contribution[no_zx].size() - 1] << endl;
      coll_tcga_contribution[no_xx].push_back(tx1 * temp_C1ggdelta_A_F); // scale dependent
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tcga   [no_xx = " << setw(2) << no_xx << "] += " << coll_tcga_contribution[no_xx][coll_tcga_contribution[no_xx].size() - 1] << endl;

      coll_tgamma2_contribution[no_zx].push_back(1.5 * Kappa / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] );
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tgamma2[no_zx = " << setw(2) << no_zx << "] += " << coll_tgamma2_contribution[no_zx][coll_tgamma2_contribution[no_zx].size() - 1] << endl;
      coll_tgamma2_contribution[no_zx].push_back(P2gg(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tgamma2[no_zx = " << setw(2) << no_zx << "] += " << coll_tgamma2_contribution[no_zx][coll_tgamma2_contribution[no_zx].size() - 1] << endl;
      coll_tgamma2_contribution[no_xx].push_back(-1.5 * Kappa * psi_z_coll[x_e] / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tgamma2[no_xx = " << setw(2) << no_xx << "] += " << coll_tgamma2_contribution[no_xx][coll_tgamma2_contribution[no_xx].size() - 1] << endl;
      coll_tgamma2_contribution[no_xx].push_back(-1.5 * Kappa * D0int(psi_x_pdf[x_e]));
      logger << LOG_DEBUG_VERBOSE << "g->g(gg):      tgamma2[no_xx = " << setw(2) << no_xx << "] += " << coll_tgamma2_contribution[no_xx][coll_tgamma2_contribution[no_xx].size() - 1] << endl;
    }

    else if (ncollinear[i_c].type_splitting_full[x_e] == 23 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with g from Q -> Q (g) -> g (g)  or  Q -> g (Q) -> g (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with g from Q -> Q (g) -> g (g)  or  Q -> g (Q) -> g (g) splitting" << endl;
      coll_tgaga_contribution[no_zx].push_back((Pgqqq(psi_z_coll[x_e]) + Pgggq(psi_z_coll[x_e], beta0)) / psi_g_z_coll[x_e]);
 
      coll_tcga_contribution[no_zx].push_back(CgqPqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      // reactivate this term: !!!
      coll_tcga_contribution[no_zx].push_back(Pgq(psi_z_coll[x_e]) * C1ggdelta(A_F) / psi_g_z_coll[x_e]); // scale dependent
      // until here.

      coll_tgamma2_contribution[no_zx].push_back(P2gq(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
    }

    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 1 && ncollinear[i_c].type_splitting_full[y_e] == 1){ 
      // hard process with both legs x_e and y_e with g from g -> g (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << " with g from g -> g (g) splitting" << endl;

      double tx1 = -3 * psi_z_coll[x_e] / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] - 3 * D0int(psi_x_pdf[x_e]) + beta0;
      double tz1 = (3 / (1. - psi_z_coll[x_e]) + Pggreg(psi_z_coll[x_e])) / psi_g_z_coll[x_e];
      double tx2 = -3 * psi_z_coll[y_e] / (1. - psi_z_coll[y_e]) / psi_g_z_coll[y_e] - 3 * D0int(psi_x_pdf[y_e]) + beta0;
      double tz2 = (3 / (1. - psi_z_coll[y_e]) + Pggreg(psi_z_coll[y_e])) / psi_g_z_coll[y_e];

      coll_tgaga_contribution[no_zz].push_back(2 * tz1 * tz2);
      coll_tgaga_contribution[no_zx].push_back(2 * tz1 * tx2);
      coll_tgaga_contribution[no_xz].push_back(2 * tx1 * tz2);
      coll_tgaga_contribution[no_xx].push_back(2 * tx1 * tx2);
    }

    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 3 && ncollinear[i_c].type_splitting_full[y_e] == 3){
      // hard process with both legs x_e and y_e with g from q -> g (q) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << "with g from q -> g (q) splitting" << endl;
      coll_tgaga_contribution[no_zz].push_back(2 * Pgq(psi_z_coll[x_e]) * Pgq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);

      coll_tcga_contribution[no_zz].push_back(Pgq(psi_z_coll[x_e]) * Cgq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      coll_tcga_contribution[no_zz].push_back(Cgq(psi_z_coll[x_e]) * Pgq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);

      /*
      coll_tgamma2_contribution[no_zz].push_back();
      */
    }

    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 1 && ncollinear[i_c].type_splitting_full[y_e] == 3){
      // hard process with leg x_e with g from g -> g (g) splitting and leg y_e with g from q -> g (q) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with leg " << x_e << " with g from g -> g (g) splitting and leg " << y_e << " with g from q -> g (q) splitting" << endl;
      double tx1 = -3 * psi_z_coll[x_e] / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] - 3 * D0int(psi_x_pdf[x_e]) + beta0;
      double tz1 = (3 / (1. - psi_z_coll[x_e]) + Pggreg(psi_z_coll[x_e])) / psi_g_z_coll[x_e];
      coll_tgaga_contribution[no_zz].push_back(-2 * tz1 * log(psi_x_pdf[y_e]) * Pgq(psi_z_coll[y_e]));
      coll_tgaga_contribution[no_xz].push_back(-2 * tx1 * log(psi_x_pdf[y_e]) * Pgq(psi_z_coll[y_e]));

      coll_tcga_contribution[no_zz].push_back(Cgq(psi_z_coll[y_e]) * tz1 / psi_g_z_coll[y_e]);
      coll_tcga_contribution[no_xz].push_back(Cgq(psi_z_coll[y_e]) * tx1 / psi_g_z_coll[y_e]);

      /*
      coll_tgamma2_contribution[no_zz].push_back();
      coll_tgamma2_contribution[no_xz].push_back();
      */
    }

    else if (ncollinear[i_c].type_splitting_full[x_e] == 14 && ncollinear[i_c].type_splitting_full[y_e] == 0){ 
      // hard process with q from g -> g (g) -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from g -> g (g) -> q (qx) splitting" << endl;
      coll_tgaga_contribution[no_zx].push_back(1. / psi_g_z_coll[x_e] * (Pqqqg(psi_z_coll[x_e]) + Pqggg(psi_z_coll[x_e], beta0)));
      coll_tcga_contribution[no_zx].push_back((CqqPqg(psi_z_coll[x_e]) + CqgPgg(psi_z_coll[x_e], beta0)) / psi_g_z_coll[x_e]);
      coll_tgamma2_contribution[no_zx].push_back(P2qg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 34 && ncollinear[i_c].type_splitting_full[y_e] == 0){ 
      // hard process with q from Q -> g (Q) -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from Q -> g (Q) -> q (qx) splitting" << endl;
      coll_tgaga_contribution[no_zx].push_back(1. / psi_g_z_coll[x_e] * Pqggq(psi_z_coll[x_e]));
      coll_tcga_contribution[no_zx].push_back(CqgPgq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tgamma2_contribution[no_zx].push_back(P2qqS(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 22 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with q from q -> q (g) -> q (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from q -> q (g) -> q (g) splitting" << endl;
      double temp_Dqqqq = (D0qqqq / (1. - psi_z_coll[x_e]) + D1qqqq * log(1. - psi_z_coll[x_e]) / (1. - psi_z_coll[x_e])) / psi_g_z_coll[x_e];
      coll_tgaga_contribution[no_zx].push_back(temp_Dqqqq);
      coll_tgaga_contribution[no_zx].push_back(Pqqqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tgaga_contribution[no_xx].push_back(-psi_z_coll[x_e] * temp_Dqqqq);
      coll_tgaga_contribution[no_xx].push_back((Deltaqqqq - D0qqqq * D0int(psi_x_pdf[x_e]) - D1qqqq * D1int(psi_x_pdf[x_e])));
      coll_tcga_contribution[no_zx].push_back(CqqPqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tgamma2_contribution[no_zx].push_back(P2qqV(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      coll_tgamma2_contribution[no_zx].push_back(2.0 / 3 * Kappa * 1. / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tgamma2_contribution[no_xx].push_back(-2.0 / 3 * Kappa * (psi_z_coll[x_e] / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e] + D0int(psi_x_pdf[x_e])));
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 54 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with q from qx -> ... -> q  splitting (non-singlet)
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from qx -> ... -> q  splitting (non-singlet)" << endl;
      coll_tgamma2_contribution[no_zx].push_back(P2qqbV(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
    }
    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 2 && ncollinear[i_c].type_splitting_full[y_e] == 2){ 
      // hard process with both legs x_e and y_e with q from q -> q (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << " with q from q -> q (g) splitting" << endl;
      double temp_Pqq_zz = 2 * Pqq(psi_z_coll[x_e]) * Pqq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e];
      coll_tgaga_contribution[no_zz].push_back(temp_Pqq_zz);
      coll_tgaga_contribution[no_xz].push_back(-2 * Pqq(psi_z_coll[y_e]) * Pqqint(psi_x_pdf[x_e]) / psi_g_z_coll[y_e]);
      coll_tgaga_contribution[no_xz].push_back(-temp_Pqq_zz * psi_z_coll[x_e]);
      coll_tgaga_contribution[no_zx].push_back(-2 * Pqq(psi_z_coll[x_e]) * Pqqint(psi_x_pdf[y_e]) / psi_g_z_coll[x_e]);
      coll_tgaga_contribution[no_zx].push_back(-temp_Pqq_zz * psi_z_coll[y_e]);
      coll_tgaga_contribution[no_xx].push_back(2 * Pqqint(psi_x_pdf[x_e]) * Pqqint(psi_x_pdf[y_e]));
      coll_tgaga_contribution[no_xx].push_back(2 * Pqq(psi_z_coll[x_e]) * Pqqint(psi_x_pdf[y_e]) / psi_g_z_coll[x_e] * psi_z_coll[x_e]);
      coll_tgaga_contribution[no_xx].push_back(2 * Pqq(psi_z_coll[y_e]) * Pqqint(psi_x_pdf[x_e]) / psi_g_z_coll[y_e] * psi_z_coll[y_e] );
      coll_tgaga_contribution[no_xx].push_back(temp_Pqq_zz * psi_z_coll[x_e] * psi_z_coll[y_e]);
      double temp_C1qqdelta_A_F = C1qqdelta(A_F);
      //      double temp_C1qqdelta_A_F = C1qqdelta(value_A_F[sd][ss]);
      coll_tcga_contribution[no_zz].push_back((Cqq(psi_z_coll[x_e]) * Pqq(psi_z_coll[y_e]) + Pqq(psi_z_coll[x_e]) * Cqq(psi_z_coll[y_e])) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      coll_tcga_contribution[no_xz].push_back(temp_C1qqdelta_A_F * Pqq(psi_z_coll[y_e]) / psi_g_z_coll[y_e]);
      coll_tcga_contribution[no_xz].push_back(-Cqq(psi_z_coll[y_e]) / psi_g_z_coll[y_e] * (Pqq(psi_z_coll[x_e]) * psi_z_coll[x_e] / psi_g_z_coll[x_e] + Pqqint(psi_x_pdf[x_e])));
      coll_tcga_contribution[no_zx].push_back(Pqq(psi_z_coll[x_e]) * temp_C1qqdelta_A_F / psi_g_z_coll[x_e]);
      coll_tcga_contribution[no_zx].push_back(-Cqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e] * (Pqq(psi_z_coll[y_e]) * psi_z_coll[y_e] / psi_g_z_coll[y_e] + Pqqint(psi_x_pdf[y_e])));
      coll_tcga_contribution[no_xx].push_back(-temp_C1qqdelta_A_F * (psi_z_coll[x_e] * Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e] + Pqqint(psi_x_pdf[x_e])));
      coll_tcga_contribution[no_xx].push_back(-temp_C1qqdelta_A_F * (psi_z_coll[y_e] * Pqq(psi_z_coll[y_e]) / psi_g_z_coll[y_e] + Pqqint(psi_x_pdf[y_e])));


    }
    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 4 && ncollinear[i_c].type_splitting_full[y_e] == 4){
      // hard process with both legs x_e and y_e with q from g -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << "with q from g -> q (qx) splitting" << endl;
      coll_tgaga_contribution[no_zz].push_back(2 * Pqg(psi_z_coll[x_e]) * Pqg(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      coll_tcga_contribution[no_zz].push_back((Pqg(psi_z_coll[x_e]) * Cqg(psi_z_coll[y_e]) + Cqg(psi_z_coll[x_e]) * Pqg(psi_z_coll[y_e])) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
    }
    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 2 && ncollinear[i_c].type_splitting_full[y_e] == 4){
      // hard process with leg x_e with q from q -> q (g) splitting and leg y_e with q from g -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with leg " << x_e << " with q from q -> q (g) splitting and leg " << y_e << " with q from g -> q (qx) splitting" << endl;
      coll_tgaga_contribution[no_zz].push_back(2 * Pqg(psi_z_coll[y_e]) / psi_g_z_coll[y_e] * Pqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tgaga_contribution[no_xz].push_back(-2 * Pqg(psi_z_coll[y_e]) / psi_g_z_coll[y_e] * (Pqq(psi_z_coll[x_e]) * psi_z_coll[x_e] / psi_g_z_coll[x_e] + Pqqint(psi_x_pdf[x_e])));
      coll_tcga_contribution[no_zz].push_back((Pqq(psi_z_coll[x_e]) * Cqg(psi_z_coll[y_e]) + Cqq(psi_z_coll[x_e]) * Pqg(psi_z_coll[y_e])) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      coll_tcga_contribution[no_xz].push_back(-(Pqq(psi_z_coll[x_e]) * psi_z_coll[x_e] / psi_g_z_coll[x_e] + Pqqint(psi_x_pdf[x_e])) * Cqg(psi_z_coll[y_e]) / psi_g_z_coll[y_e]);
      //      coll_tcga_contribution[no_xz].push_back(C1qqdelta(value_A_F[sd][ss]) * Pqg(psi_z_coll[y_e]) / psi_g_z_coll[y_e]);
      coll_tcga_contribution[no_xz].push_back(C1qqdelta(A_F) * Pqg(psi_z_coll[y_e]) / psi_g_z_coll[y_e]);
    }
    else {logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << "   ncollinear[" << i_c << "].type_splitting_full[" << y_e << "] = " << ncollinear[i_c].type_splitting_full[y_e] << " -> selected splitting type not allowed !!!" << endl; exit(1);}
  }
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    for (int i_c = 0; i_c < coll_tcga_contribution[i_l].size(); i_c++){
      logger << LOG_DEBUG_VERBOSE << "tcga[" << i_l << "][" << i_c << "] = " << coll_tcga_contribution[i_l][i_c] << endl;
    }

    coll_tgaga[i_l] = accumulate(coll_tgaga_contribution[i_l].begin(), coll_tgaga_contribution[i_l].end(), 0.);
    coll_tcga[i_l] = accumulate(coll_tcga_contribution[i_l].begin(), coll_tcga_contribution[i_l].end(), 0.);
    coll_tgamma2[i_l] = accumulate(coll_tgamma2_contribution[i_l].begin(), coll_tgamma2_contribution[i_l].end(), 0.);
  }
  
  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}

void observable_set::determine_splitting_tH2(phasespace_set & psi){
  static Logger logger("observable_set::determine_splitting_tH2");
  logger << LOG_DEBUG_VERBOSE << "started" << endl;

  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){coll_tH2_contribution[i_l].clear();}

  ///  static double H2qqD0 = -404. / 27 + (56. * N_f) / 81 + 14. * zeta3;

  //  coll: prefactors of splittings, to be multiplied with pdf factors
  //  coll[0] -> tH2
  for (int i_c = 1; i_c < ncollinear.size(); i_c++){
    // type_splitting_full (leg-wise):
    // x_e y_e      --- ---

    // 14 -  0   +   0 - 14   // hard process with q from g -> g (g) -> q (qx) splitting             (-> 4  in sig11 calculation)
    //(42 -  0   +   0 - 42   // hard process with q from g -> q (qx) -> q (g) splitting             (part of  14 -  0   +   0 - 14))

    // 34 -  0   +   0 - 34   // hard process with q from Q -> g (Q) -> q (qx) splitting             (----  no equivalent
    // 54 -  0   +   0 - 54   // hard process with q from qx -> ... -> q  splitting (non-singlet)    (----  no equivalent

    // 11 -  0   +   0 - 11   // hard process with g from g -> g (g) -> g (g) splitting              (-> 1  in sig11 calculation)
    //(43 -  0   +   0 - 43   // hard process with g from g -> Q (Qx) -> g (Q) splitting             (part of  11 -  0   +   0 - 11))

    // 23 -  0   +   0 - 23   // hard process with g from Q -> Q (g) -> g (Q) splitting              (-> 3  in sig11 calculation)
    //(31 -  0   +   0 - 31   // hard process with g from Q -> g (Q) -> g (g) splitting              (part of  23 -  0   +   0 - 23))

    // 22 -  0   +   0 - 22   // hard process with q from q -> q (g) -> q (g) splitting              (-> 2  in sig11 calculation)

    //  2 -  2                // hard process with (leg 1) q from q -> q (g) splitting and (leg 2) q from q -> q (g) splitting
    //  2 -  4   +   4 -  2   // hard process with (leg 1) q from q -> q (g) splitting and (leg 2) g from q -> g (q) splitting
    //  4 -  4                // hard process with (leg 1) g from q -> g (q) splitting and (leg 2) g from q -> g (q) splitting





    //  gg channel  
    // 11 - 0 (with a factor .5) or 1 - 1 (without !!!)
    // analogously to qqbar channel: 1 - 1
    //  tH2 = tH2 + H2_delta*pdf_factor_x1x2;
    // 11 - 0
    //  tH2 = tH2 - 0.5*( log(x1)*(pdf_factor_z1x2-pdf_factor_x1x2*z1)*H2ggD0/(1-z1) + H2ggD0*D0int(x1)*pdf_factor_x1x2 );
    // second leg (automatically included)  tH2 = tH2 - 0.5*( log(x2)*(pdf_factor_x1z2-pdf_factor_x1x2*z2)*H2ggD0/(1-z2) + H2ggD0*D0int(x2)*pdf_factor_x1x2 );
    // 11 - 0
    //  tH2 = tH2 - log(x1)*pdf_factor_z1x2*C2ggreg(z1,nf); 
    // second leg (automatically included)  tH2 = tH2 - log(x2)*pdf_factor_x1z2*C2ggreg(z2,nf);
    // 1 - 1
    //  tH2 = tH2 + log(x1)*log(x2)*Ggg(z1)*Ggg(z2)*pdf_factor_z1z2;
    //  qg channel  
    // 23 - 0
    //  tH2 = tH2 - log(x1)*pdf_factor_qx2*C2gq(z1,nf);
    // second leg (automatically included)  tH2 = tH2 - log(x2)*pdf_factor_x1q*C2gq(z2,nf);
    //  tH2 = tH2 - log(x1)*H1_delta*pdf_factor_qx2*Cgq(z1);          // doesn't appear explicitly in HHNLO
    // second leg (automatically included)  tH2 = tH2 - log(x2)*H1_delta*pdf_factor_x1q*Cgq(z2);          // doesn't appear explicitly in HHNLO 
    // 1 - 3 (both needed)
    //  tH2 = tH2 + log(x1)*log(x2)*Ggg(z1)*Ggq(z2)*pdf_factor_z1q;
    //  tH2 = tH2 + log(x1)*log(x2)*Ggq(z1)*Ggg(z2)*pdf_factor_qz2;
    //  qq channel (C1C1 and G1G1)
    // 3 - 3
    //  tH2 = tH2 + log(x1)*log(x2)*Cgq(z1)*Cgq(z2)*pdf_factor_qq;    
    //  tH2 = tH2 + log(x1)*log(x2)*Ggq(z1)*Ggq(z2)*pdf_factor_qq;


    //    cout << ncollinear[i_c].type_splitting_full[x_e] << "   " << ncollinear[i_c].type_splitting_full[y_e] << "   no_xx = " << no_xx << "   no_zx = " << no_zx << "   no_xz = " << no_xz << "   no_zz = " << no_zz << endl;   


    if (ncollinear[i_c].type_splitting_full[x_e] == 0){logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = 0 -> no splitting type defined !!!" << endl; exit(1);}
    else if (ncollinear[i_c].type_splitting_full[x_e] == 11 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with g from g -> g (g) -> g (g)  or  g -> Q (Qx) -> g (Q) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with g from g -> g (g) -> g (g)  or  g -> Q (Qx) -> g (Q) splitting" << endl;
      coll_tH2_contribution[no_zx].push_back(.5 * H2ggD0 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tH2_contribution[no_zx].push_back(C2ggreg(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      coll_tH2_contribution[no_xx].push_back(-.5 * psi_z_coll[x_e] * H2ggD0 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tH2_contribution[no_xx].push_back(-.5* H2ggD0 * D0int(psi_x_pdf[x_e]));
    }

    else if (ncollinear[i_c].type_splitting_full[x_e] == 23 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with g from Q -> Q (g) -> g (g)  or  Q -> g (Q) -> g (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with g from Q -> Q (g) -> g (g)  or  Q -> g (Q) -> g (g) splitting" << endl;
      coll_tH2_contribution[no_zx].push_back(C2gq(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      coll_tH2_contribution[no_zx].push_back(QT_H1_delta * Cgq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
    }

    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 1 && ncollinear[i_c].type_splitting_full[y_e] == 1){ 
      // hard process with both legs x_e and y_e with g from g -> g (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << " with g from g -> g (g) splitting" << endl;

      if (QT_finalstate_massive_coloured){
	// !!! works only for ttx at present, until general implementation available...
 	coll_tH2_contribution[no_zz].push_back(QT_H0_doublespinflip * Ggg(psi_z_coll[x_e]) * Ggg(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
	coll_tH2_contribution[no_xz].push_back(QT_H0_DG * Ggg(psi_z_coll[y_e]) / psi_g_z_coll[y_e]);
	coll_tH2_contribution[no_zx].push_back(QT_H0_DG * Ggg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      }
      else {
	coll_tH2_contribution[no_zz].push_back(Ggg(psi_z_coll[x_e]) * Ggg(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      }
      
      coll_tH2_contribution[no_xx].push_back(QT_H2_delta);
    }

    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 3 && ncollinear[i_c].type_splitting_full[y_e] == 3){
      // hard process with both legs x_e and y_e with g from q -> g (q) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << "with g from q -> g (q) splitting" << endl;
      coll_tH2_contribution[no_zz].push_back(Cgq(psi_z_coll[x_e]) * Cgq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);

      if (QT_finalstate_massive_coloured){
	// !!! works only for ttx at present, until general implementation available...
	coll_tH2_contribution[no_zz].push_back(QT_H0_doublespinflip * Ggq(psi_z_coll[x_e]) * Ggq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      }
      else {
	coll_tH2_contribution[no_zz].push_back(Ggq(psi_z_coll[x_e]) * Ggq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      }
    }

    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 1 && ncollinear[i_c].type_splitting_full[y_e] == 3){
      // hard process with leg x_e with g from g -> g (g) splitting and leg y_e with g from q -> g (q) splitting
      // first leg is Xgg, second is Xgq (called twice!) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // no_zx may not be used here (it's automatically included from no_xz by the inverted call, and "no_zx" points would not be assigned the correct PDFs !!!).
      
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with leg " << x_e << " with g from g -> g (g) splitting and leg " << y_e << " with g from q -> g (q) splitting" << endl;

      if (QT_finalstate_massive_coloured){
	// !!! works only for ttx at present, until general implementation available...
	coll_tH2_contribution[no_zz].push_back(QT_H0_doublespinflip * Ggg(psi_z_coll[x_e]) * Ggq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
	coll_tH2_contribution[no_xz].push_back(QT_H0_DG * Ggq(psi_z_coll[y_e]) / psi_g_z_coll[y_e]);
      }
      else {
	coll_tH2_contribution[no_zz].push_back(Ggg(psi_z_coll[x_e]) * Ggq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      }
      
    }

    else if (ncollinear[i_c].type_splitting_full[x_e] == 14 && ncollinear[i_c].type_splitting_full[y_e] == 0){ 
      // hard process with q from g -> g (g) -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from g -> g (g) -> q (qx) splitting" << endl;
      coll_tH2_contribution[no_zx].push_back(C2qg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      // no corresponding term (23-0) proportional to QT_H1_delta !!!
      // could be the one that is added outside of tH2 ???
      // new term, shifted here from outside this function:
      coll_tH2_contribution[no_zx].push_back(QT_H1_delta * Cqg(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 34 && ncollinear[i_c].type_splitting_full[y_e] == 0){ 
      // hard process with q from Q -> g (Q) -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from Q -> g (Q) -> q (qx) splitting" << endl;
      // this contribution should be subtracted again for the specified quarks q and qx, which are treated in 22 and 54, respectively !!!
      coll_tH2_contribution[no_zx].push_back(C2qqp(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 22 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with q from q -> q (g) -> q (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from q -> q (g) -> q (g) splitting" << endl;
      coll_tH2_contribution[no_zx].push_back(.5 * H2qqD0 / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tH2_contribution[no_zx].push_back(C2qqreg(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      // this contribution should be subtracted again for the specified quarks q and qx, which are treated in 22 and 54, respectively (from 34)!!!
      coll_tH2_contribution[no_zx].push_back(-C2qqp(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      coll_tH2_contribution[no_xx].push_back(-.5 * H2qqD0 * psi_z_coll[x_e] / (1. - psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
      coll_tH2_contribution[no_xx].push_back(-.5 * H2qqD0 * D0int(psi_x_pdf[x_e]));
      // new term, shifted here from outside this function:
      coll_tH2_contribution[no_zx].push_back(QT_H1_delta * Cqq(psi_z_coll[x_e]) / psi_g_z_coll[x_e]);
    }
    else if (ncollinear[i_c].type_splitting_full[x_e] == 54 && ncollinear[i_c].type_splitting_full[y_e] == 0){
      // hard process with q from qx -> ... -> q  splitting (non-singlet)
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with q from qx -> ... -> q  splitting (non-singlet)" << endl;
      coll_tH2_contribution[no_zx].push_back(C2qqb(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
      // this contribution should be subtracted again for the specified quarks q and qx, which are treated in 22 and 54, respectively (from 34)!!!
      coll_tH2_contribution[no_zx].push_back(-C2qqp(psi_z_coll[x_e], N_f) / psi_g_z_coll[x_e]);
    }
    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 2 && ncollinear[i_c].type_splitting_full[y_e] == 2){ 
      // hard process with both legs x_e and y_e with q from q -> q (g) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << " with q from q -> q (g) splitting" << endl;
      coll_tH2_contribution[no_zz].push_back(Cqq(psi_z_coll[x_e]) * Cqq(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
      coll_tH2_contribution[no_xx].push_back(QT_H2_delta);
    }
    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 4 && ncollinear[i_c].type_splitting_full[y_e] == 4){
      // hard process with both legs x_e and y_e with q from g -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with both legs " << x_e << " and " << y_e << "with q from g -> q (qx) splitting" << endl;
      coll_tH2_contribution[no_zz].push_back(Cqg(psi_z_coll[x_e]) * Cqg(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
    }
    else if (ncollinear[i_c].leg_emission == 0 && ncollinear[i_c].type_splitting_full[x_e] == 2 && ncollinear[i_c].type_splitting_full[y_e] == 4){
      // hard process with leg x_e with q from q -> q (g) splitting and leg y_e with q from g -> q (qx) splitting
      logger << LOG_DEBUG_VERBOSE << setw(19) << ncollinear[i_c].name << "   [" << setw(2) << i_c << "]   hard process with leg " << x_e << " with q from q -> q (g) splitting and leg " << y_e << " with q from g -> q (qx) splitting" << endl;
      coll_tH2_contribution[no_zz].push_back(Cqq(psi_z_coll[x_e]) * Cqg(psi_z_coll[y_e]) / psi_g_z_coll[x_e] / psi_g_z_coll[y_e]);
    }
    else {logger << LOG_FATAL << "Should not happen: ncollinear[" << i_c << "].type_splitting_full[" << x_e << "] = " << ncollinear[i_c].type_splitting_full[x_e] << "   ncollinear[" << i_c << "].type_splitting_full[" << y_e << "] = " << ncollinear[i_c].type_splitting_full[y_e] << " -> selected splitting type not allowed !!!" << endl; exit(1);}
  }
  for (int i_l = 0; i_l < list_combination_pdf.size(); i_l++){
    coll_tH2[i_l] = accumulate(coll_tH2_contribution[i_l].begin(), coll_tH2_contribution[i_l].end(), 0.);
  }

  logger << LOG_DEBUG_VERBOSE << "finished" << endl;
}
#undef no_zz 
#undef no_xz
#undef no_zx 
#undef no_xx
#undef x_e
#undef y_e


