#pragma once

#include <vector>

using namespace std;

#include "../../classes/header/fourvector.h"
#include "../../classes/header/particle.h"

extern "C" {
	extern struct{
		int verbose_mode;
	} handling_params_;
}

extern "C" {
	extern struct{
		unsigned char name[63];
		int member,ih1,ih2;
	} pdfhres_;
}

extern "C" {
	extern struct{
		int resorder,pertorder,approxNNLL;
	} resorder_;
}

extern "C" {
	extern struct{
		int gg_initiated,channel,pdf_content_modify;
	} resumchan_;
}

extern "C" {
	extern struct{
		int nf;
	} nf_;
}

extern "C" {
	extern struct{
		double pcf,kcf;
	} lopowers_;
}

extern "C" {
	extern struct{
		int doinitres;
	} flaginitres_;
}

extern "C" {
	extern struct{
		double qt,xtau,q,q2,mur2,muf2,qres2,startscale;
	} somescale_;
}

extern "C" {
	extern struct{
		double md,mu,ms,mc,mbot,mtop,mel,mmu,mtau,higgsmass,hwidth,wmass,wwidth,zmass,zwidth,width,tauwidth,mtausq,mcsq,mbsq;
	} masses_;
}

extern "C" {
	void pdfset_();
}

extern "C" {
	void resumm_(double *resummed,int *gg_initiated,int *channel,double *Msyst,double *ysyst,double *pTsyst,double *Qres,double *muF,double *muR,double *sroot,double *bornfactor,double *H1_in,double *H2_in);
}

extern "C" {
	void init_const_and_coeff_();
}

void ConvertToFortran(char* fstring, std::size_t fstring_len, const char* cstring);
void performQTBoost(double qt2,double Q,double y,vector<vector<vector<particle> > > &p);
void performQTBoost_cms(double qt2,double Q,double y,vector<vector<fourvector> > &p);
