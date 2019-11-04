{

if (sd == 1){
  fourvector Q;
  for (int i_p = 3; i_p < 5; i_p++){Q = Q + osi_p_parton[i_a][i_p];}
  temp_mu_central = Q.m();
}
else if (sd == 2){
  // H_T
  double msq  = PARTICLE("top")[0].m2;
  double pTsq1 = PARTICLE("top")[0].pT2;
  double pTsq2 = PARTICLE("atop")[0].pT2;
  temp_mu_central = sqrt(msq + pTsq1) + sqrt(msq + pTsq2);
}

}

