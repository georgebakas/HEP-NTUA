{
  if (sd == 1){
    // invariant mass m of the ttx system
    temp_mu_central = (PARTICLE("top")[0].momentum + PARTICLE("atop")[0].momentum).m();
  }
  else if (sd == 2){
    // sum of transverse masses of t and tx - HT
    double ET_t = PARTICLE("top")[0].ET;
    double ET_tx = PARTICLE("atop")[0].ET;
    temp_mu_central = ET_t + ET_tx;
  }
  else if (sd == 3){
    // transverse mass of the ttx system
    double m  = (PARTICLE("top")[0].momentum + PARTICLE("atop")[0].momentum).m();
    double pT = (PARTICLE("top")[0].momentum + PARTICLE("atop")[0].momentum).pT();
    temp_mu_central = sqrt(pow(m, 2) + pow(pT, 2));
  }
  else if (sd == 4){
    // geometric average of transverse masses of t and tx - ET
    double ET_t = PARTICLE("top")[0].ET;
    double ET_tx = PARTICLE("atop")[0].ET;
    temp_mu_central = sqrt(ET_t * ET_tx);
  }
  else if (sd == 5){
    // transverse mass of t - mT_t
    temp_mu_central = PARTICLE("top")[0].ET;
  }
  else if (sd == 6){
    // transverse mass of tx - mT_tx
    temp_mu_central = PARTICLE("atop")[0].ET;
  }
  else if (sd == 7){
    // transverse mass of the leading t/tx - mT_t1
    temp_mu_central = PARTICLE("tjet")[0].ET;
  }
  else if (sd == 8){
    // transverse mass of the subleading t/tx - mT_t2
    temp_mu_central = PARTICLE("tjet")[1].ET;
  }
}
