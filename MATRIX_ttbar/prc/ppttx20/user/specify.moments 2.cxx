if (osi_n_moments == 4){
  for (int i_a = 0; i_a < osi_n_ps; i_a++){
    osi_moment[0][i_a] = 1.;
    double mu_dynamic = 1.;
    double log_mu_dynamic = log(mu_dynamic);
    osi_moment[1][i_a] = mu_dynamic;
    osi_moment[2][i_a] = pow(mu_dynamic, 2);
    if (log_mu_dynamic == log_mu_dynamic){
      osi_moment[3][i_a] = log_mu_dynamic;
    }
    else {
      cout << "log_mu_dynamic = nan" << endl;
      osi_moment[3][i_a] = 0.;
    }
  }
}
else {
}
