    /*
    if ((oset.name_process[0] == 'g' && oset.name_process[1] == 'g') || (oset.name_process[0] != 'g' && oset.name_process[1] != 'g')){
//    if (oset.name_process[0] == 'g' && oset.name_process[1] == 'g'){
//      logger << LOG_FATAL << "VT or CT2 term for gg channel is not implemented yet !!!" << endl;
//
//      exit(1);
//    }
//    if (oset.name_process[0] != 'g' && oset.name_process[1] != 'g'){
//      logger << LOG_FATAL << "VT or CT2 term for qqx channel is not implemented yet !!!" << endl;

      logger << LOG_FATAL << "VT or CT2 term for gg and qqx channel is not implemented yet !!!" << endl;
      //void ol_tree_colbasis_dim(int id, int* ncolb, int* colelemsz, int* nheltot); 

      int ncolb = 0;
      int colelemsz = 0;
      int nheltot = 0;

      ol_tree_colbasis_dim(1, &ncolb, &colelemsz, &nheltot); 

      cout << "ncolb = " << ncolb << endl;
      cout << "colelemsz = " << colelemsz << endl;
      cout << "nheltot = " << nheltot << endl;

      //      int basis[ncolb][colelemsz];
      //      int needed[ncolb][ncolb];

      int *basis;
      int *needed;
      basis = new int[ncolb*colelemsz];
      needed = new int[ncolb*ncolb];
      int cpp_basis[ncolb][colelemsz];
      int cpp_needed[ncolb][ncolb];
      // int basis[ncolb][colelemsz];
      // needed[ncolb][ncolb];
      ol_tree_colbasis(1, basis, needed);

 
      cout << endl;
      for (int i_c = 0; i_c < ncolb; i_c++){
	for (int i_i = 0; i_i < colelemsz; i_i++){
	  cpp_basis[i_c][i_i] = basis[i_c * colelemsz + i_i];
	  cout << "basis[" << i_c << "][" << i_i << "] = " << setw(5) << cpp_basis[i_c][i_i] << setw(5) << basis[i_c * colelemsz + i_i] << endl;
	}
	cout << endl;
      }

      // qqx -> ttx
      // 3 1 0 2 4 -> T^{i1}_31 * T^{i1}_24
      // 3 4 0 2 1 -> T^{i1}_34 * T^{i1}_21

      // gg -> ttx
      // 1 2 3 4 0 -> [T^{1}T^{2}]_34
      // 2 1 3 4 0 -> [T^{2}T^{1}]_34

      for (int i_c = 0; i_c < ncolb; i_c++){
	for (int j_c = 0; j_c < ncolb; j_c++){
	  cpp_needed[i_c][j_c] = needed[i_c * ncolb + j_c];
	  cout << "needed[" << i_c << "][" << j_c << "] = " << setw(5) << needed[i_c * ncolb + j_c] << setw(5) << cpp_needed[i_c][j_c] << endl;
	}
      }

      double *amp;
      amp = new double[nheltot*2*ncolb];
      double_complex cpp_amp[nheltot][ncolb];

      double cpp_amp2[nheltot][ncolb];

      // double amp[nheltot][2*ncolb];
      int nhel = 0;

      cout << "nheltot = " << nheltot << endl;

      ol_evaluate_tree_colvect(1, P, amp, &nhel); 

      cout << "nhel = " << nhel << endl;

      double amp2_col0 = 0.;
      for (int i_h = 0; i_h < nheltot; i_h++){
	for (int i_c = 0; i_c < ncolb; i_c++){
	  //	  cpp_amp[i_h][i_c] = double_complex(amp[2 * (i_h * ncolb + i_c)], amp[2 * (i_h * ncolb + i_c) + 1]);
	  //	  cout << "amp[" << i_h << "][" << i_c << "] = " << setprecision(15) << setw(23) << amp[2 * (i_h * ncolb + i_c)] << " + " << setprecision(15) << setw(23) << amp[2 * (i_h * ncolb + i_c) + 1] << "i   " << setprecision(15) << setw(23) << cpp_amp[i_h][i_c] << endl;
	  cpp_amp[i_h][i_c] = double_complex(amp[2 * i_h + 2 * i_c * nheltot], amp[2 * i_h + 2 * i_c * nheltot + 1]);
	  cpp_amp2[i_h][i_c] = pow(abs(cpp_amp[i_h][i_c]), 2);
	  cout << "amp[" << i_h << "][" << i_c << "] = " << setprecision(15) << setw(23) << amp[2 * i_h + 2 * i_c * nheltot] << " + " << setprecision(15) << setw(23) << amp[2 * i_h + 2 * i_c * nheltot + 1] << "i   " << setprecision(15) << setw(23) << cpp_amp[i_h][i_c] << endl;
	}
	  amp2_col0 += cpp_amp2[i_h][0];
      }
      cout << "amp2_col0 = " << amp2_col0 << endl;

      double b_ME2 = 0.;
      ol_evaluate_tree(1, P, &b_ME2);

      cout << "amp2_col0 / b_ME2 = " << amp2_col0 / b_ME2 << endl;

      
      for (int i_h = 0; i_h < nheltot*2*ncolb; i_h++){
	cout << "amp[" << i_h << "] = " << setprecision(15) << setw(23) << amp[i_h] << endl;
      }
      
      // for comparison
      static int n_cc = (osi_n_particle + 2) * (osi_n_particle + 1) / 2;
      static double ewcc;
      static double *M2cc;
      M2cc = new double[n_cc];
      ol_evaluate_cc(1, P, &osi_ME2, M2cc, &ewcc);

      cout << "amp2_col0 / M2cc[0] = " << amp2_col0 / M2cc[0] << endl;

      for (int i_c = 0; i_c < osi_QT_correlationoperator.size(); i_c++){
	logger << LOG_INFO << "QT_ME2_cf[" << i_c << "] = " << osi_QT_ME2_cf[i_c] << "   charge_factor = " << osi_QT_correlationoperator[i_c].charge_factor << "   no_BLHA_entry = " << osi_QT_correlationoperator[i_c].no_BLHA_entry << endl;
      } 

      for (int i_c = 0; i_c < n_cc; i_c++){
	logger << LOG_INFO << "M2cc[" << i_c << "] = " << M2cc[i_c] << endl;
      }

      delete [] M2cc;




      exit(1);
    }
*/
