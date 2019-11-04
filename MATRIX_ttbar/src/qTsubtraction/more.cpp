#ifdef MORE

#include "../include/classes.cxx"
#include "header/more.h"

void ConvertToFortran(char* fstring, std::size_t fstring_len,
                      const char* cstring)
{

    std::size_t inlen = std::strlen(cstring);
    std::size_t cpylen = std::min(inlen, fstring_len);

    if (inlen > fstring_len)
    {
        std::cout << "ConvertToFortran: illegal sizes (" << inlen << ", " << fstring_len << ")" << std::endl;
        exit(1);
    }

    std::copy(cstring, cstring + cpylen, fstring);
    std::fill(fstring + cpylen, fstring + fstring_len, ' ');
}

void initMore(int resorder, int pertorder, observable_set & oset) {
  resorder_.resorder=resorder;
  resorder_.pertorder=pertorder;

  resumchan_.pdf_content_modify=-1; //FIXME oset.pdf_content_modify;

  lopowers_.pcf=0;
  lopowers_.kcf=0;

  flaginitres_.doinitres=1;

  char test[64]="MSTW2008nnlo90cl.LHgrid               "; //FIXME
  for (int i=0; i<oset.LHAPDFname.size(); i++) {
                test[i]=oset.LHAPDFname.c_str()[i];
  }
  for (int i=oset.LHAPDFname.size(); i<63; i++) {
          test[i]=' ';
  }
  char fstring[64];
  ConvertToFortran(fstring,sizeof test,test);

  for (int i=0; i<64; i++) {
          pdfhres_.name[i]=fstring[i];
  }

  pdfhres_.member=1;
  pdfhres_.ih1=1;
  pdfhres_.ih2=1;

  nf_.nf = oset.N_f;

  pdfset_();

  somescale_.startscale=oset.var_mu_fact;
}


// boosts system to the specified transverse momentum. the azimuth is generated uniformly
void performQTBoost(double qt2,double Q,double y,vector<vector<vector<particle> > > &p) {
  // boost final state system
  fourvector q_lab;
//  double Q = sqrt(xbs_all[0][0]);
//  double y = log(sqrt(x_pdf[1]/x_pdf[2]));

//  q_lab = osi_particle_event[24][0][0].momentum()+osi_particle_event[24][0][1].momentum();
  double pt = 0.5*Q*(exp(y)+exp(-y));
  double pz = 0.5*Q*(exp(y)-exp(-y));


//  cout << "q_lab = " << q_lab << endl;
  q_lab = fourvector(pt,0,0,pz);
//  cout << "q_lab = " << q_lab << ", y = " << q_lab.rapidity() << ", eta = " << q_lab.eta() << endl;

//  cout << "pz=" << pz << ", " << q_lab.x3() << "; pt=" << pt << ", " << q_lab.x0() << endl;
  double mT = sqrt(Q*Q+qt2);

  // boost to the correct z
  double B = pt*pz/(pt*pt+mT*mT/Q/Q*pz*pz);
  double C = (pz*pz-mT*mT/Q/Q*pz*pz)/(pt*pt+mT*mT/Q/Q*pz*pz);
  double beta_z = B-sqrt(B*B-C);

  if (pz<0) {
    beta_z = B+sqrt(B*B-C);
  }

//  cout << beta_z << endl;
//  if (beta_z<0) {
//    beta_z = B+sqrt(B*B-C);
//  }
  q_lab = q_lab.zboost(beta_z);
//  cout << "correctness: " << 0.5*mT*(exp(y)-exp(-y)) << ", " << (1.0/sqrt((1-beta_z*beta_z)))*(pz-beta_z*pt) << endl;

//  cout << beta_z << endl;
//  cout << q_lab << "; " << q_lab.m() << endl;
//  cout << "pz=" << q_lab.x3() << ", " << 0.5*mT*(exp(y)-exp(-y)) << "; pt=" << q_lab.x0() << ", " << 0.5*mT*(exp(y)+exp(-y)) << endl;

  // boost to the correct qT
  double beta_x = -sqrt(qt2/(q_lab.x0()*q_lab.x0()+qt2));

  if (q_lab.x0()<0)
    beta_x *= -1.0;

  q_lab = q_lab.xboost(beta_x);
//  cout << "pz=" << q_lab.x3() << ", " << 0.5*mT*(exp(y)-exp(-y)) << "; pt=" << q_lab.x0() << ", " << 0.5*mT*(exp(y)+exp(-y)) << endl;

  // generate angle for qT
  double qT_phi = 2*M_PI*double(rand())/RAND_MAX;
//  cout << qT_phi << endl;
  q_lab = q_lab.rot(qT_phi,1);
//  cout << "q_lab = " << q_lab << "; " << q_lab.pT() << endl;

  // modify all final state momenta according to the prescription above
//  cout << "beta_z, beta_x, qT_phi = " << beta_z << ", " << beta_x << ", " << qT_phi << endl;
  for (int j1=0; j1<p.size(); j1++) {
    for (int i_a=0; i_a<p[j1].size(); i_a++) {
      for (int i_p=0; i_p<p[j1][i_a].size(); i_p++) {
          // event vector might have unused entries; skip these
          if (p[j1][i_a][i_p].momentum == physconst::nullvector) {
            continue;
          }
          //cout << j1 << ", " << i_a << ", " << i_p << endl;
          //cout << p[j1][i_a][i_p].momentum << endl;

          q_lab = p[j1][i_a][i_p].momentum;

          q_lab = q_lab.zboost(beta_z);
          q_lab = q_lab.xboost(beta_x);
          q_lab = q_lab.rot(qT_phi,1);

          p[j1][i_a][i_p] = particle(q_lab);
          //cout << p[j1][i_a][i_p].momentum << endl;
      }
    }
  }
}

// boosts system to the specified transverse momentum. the azimuth is generated uniformly
void performQTBoost_cms(double qt2,double Q,double y,vector<vector<fourvector> > &p) {
  // if we are boosting fourvectors in the cms system, we first have to transform to the lab system
  // and then transform back at the very end
  double boost_cms=0;
  boost_cms = (exp(y)-exp(-y))/(exp(y)+exp(-y));
  //cout << "boost_cms = " << boost_cms << endl;
//  boost = (x_pdf[1] - x_pdf[2]) / (x_pdf[1] + x_pdf[2]);


  // boost final state system
  fourvector q_lab;
//  double Q = sqrt(xbs_all[0][0]);
//  double y = log(sqrt(x_pdf[1]/x_pdf[2]));

//  q_lab = osi_particle_event[24][0][0].momentum()+osi_particle_event[24][0][1].momentum();
  double pt = 0.5*Q*(exp(y)+exp(-y));
  double pz = 0.5*Q*(exp(y)-exp(-y));


//  cout << "q_lab = " << q_lab << endl;
  q_lab = fourvector(pt,0,0,pz);
//  cout << "q_lab = " << q_lab << ", y = " << q_lab.rapidity() << ", eta = " << q_lab.eta() << endl;

//  cout << "pz=" << pz << ", " << q_lab.x3() << "; pt=" << pt << ", " << q_lab.x0() << endl;
  double mT = sqrt(Q*Q+qt2);

  // boost to the correct z
  double B = pt*pz/(pt*pt+mT*mT/Q/Q*pz*pz);
  double C = (pz*pz-mT*mT/Q/Q*pz*pz)/(pt*pt+mT*mT/Q/Q*pz*pz);
  double beta_z = B-sqrt(B*B-C);

  if (pz<0) {
    beta_z = B+sqrt(B*B-C);
  }

//  cout << beta_z << endl;
//  if (beta_z<0) {
//    beta_z = B+sqrt(B*B-C);
//  }
  q_lab = q_lab.zboost(beta_z);
//  cout << "correctness: " << 0.5*mT*(exp(y)-exp(-y)) << ", " << (1.0/sqrt((1-beta_z*beta_z)))*(pz-beta_z*pt) << endl;

//  cout << beta_z << endl;
//  cout << q_lab << "; " << q_lab.m() << endl;
//  cout << "pz=" << q_lab.x3() << ", " << 0.5*mT*(exp(y)-exp(-y)) << "; pt=" << q_lab.x0() << ", " << 0.5*mT*(exp(y)+exp(-y)) << endl;

  // boost to the correct qT
  double beta_x = -sqrt(qt2/(q_lab.x0()*q_lab.x0()+qt2));

  if (q_lab.x0()<0)
    beta_x *= -1.0;

  q_lab = q_lab.xboost(beta_x);
//  cout << "pz=" << q_lab.x3() << ", " << 0.5*mT*(exp(y)-exp(-y)) << "; pt=" << q_lab.x0() << ", " << 0.5*mT*(exp(y)+exp(-y)) << endl;

  // generate angle for qT
  double qT_phi = 2*M_PI*double(rand())/RAND_MAX;
//  cout << qT_phi << endl;
  q_lab = q_lab.rot(qT_phi,1);
//  cout << "q_lab = " << q_lab << "; " << q_lab.pT() << endl;

  // modify all final state momenta according to the prescription above
//  cout << "beta_z, beta_x, qT_phi = " << beta_z << ", " << beta_x << ", " << qT_phi << endl;
    for (int i_a=0; i_a<p.size(); i_a++) {
      for (int i_p=0; i_p<p[i_a].size(); i_p++) {
          // event vector might have unused entries; skip these
          if (p[i_a][i_p]== physconst::nullvector) {
            continue;
          }
          //cout << i_a << ", " << i_p << endl;
          //cout << p[i_a][i_p]<< endl;

          q_lab = p[i_a][i_p];

          q_lab = q_lab.zboost(-boost_cms);

          q_lab = q_lab.zboost(beta_z);
          q_lab = q_lab.xboost(beta_x);
          q_lab = q_lab.rot(qT_phi,1);

          q_lab = q_lab.zboost(boost_cms);

          p[i_a][i_p] = q_lab;
          //cout << p[i_a][i_p] << endl;
      }
  }
}


#endif
