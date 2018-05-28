
/*
*  vmnevofit.h
*  
*  Created by Duncan MacGregor 
*  University of Edinburgh 2017
*  Released under MIT license, see https://opensource.org/licenses/MIT
*
*/


#ifndef VMNEVOFIT_H
#define VMNEVOFIT_H


#include "evofitbasic.h"
#include "vmnmodel.h"


class EvoFitVMN_CPU : public wxThread
{
public:
	int celldex;
	int numcells;
	float *chromepop;
	int paramcount;
	int threadcount;
	int blocksize;
	int runtime;
	float *Ints;
	float *ISIs;
	float *Temp;
	float *SpikeCounts;

	DiagBox *diagbox;

	EvoFitVMN_CPU(int numcells, int index, float *chromepop, int paramcount, int threadcount, int blocksize, int runtime, float *Ints, float *ISIs, float *Temp, float *SpikeCounts);

	void SpikeGen();
	virtual void *Entry();
};


class EvoFitVMN : public EvoFit
{
public:
	VMNModel *mod;
	EvoFitBox *fitbox;
	DiagBox *diagbox;

	vector<EvoChrome> *chromepop;
	vector<EvoChrome> *chromeresult;
	EvoChrome *nextgen;
	float *chromearray;
	EvoChrome *chrome;

	EvoChrome *fitchrome;
	int chromeparams;

	int popsize;
	int parentrange;
	int generations;
	//int numparams;
	int numruns;
	double mutateprob;
	double dualfit;
	bool burstmode;
	TextFile ofp;
	bool genmon;
	bool genlysis;
	bool seedgen;
	long evoseed;
	int netmode;
	int cpumode;
	int runmode;
	int diagnostic;

	int blocksize;
	int gpuparams;
	float runtime;
	int *threaddata;
	bool multirun;

	// CPU mode
	EvoFitVMN_CPU *chromethread[1000];

	EvoFitVMN(VMNModel *, EvoFitBox *);
	//~EvoFitNet();

	void InitPop();
	void Evaluate(int start, int pop, double dual=0);
	void Evolve();

	virtual void *Entry();
};


// CUDA GPU function - takes a chrome population and runs the model with each chrome generating spike times and a quad binned ISI histogram

extern void EvoFitVMN_GPU(int netmode, int numcells, float *chromepop, int gpuparams, int popsize, int blocksize, float runtime, float *Ints, float *ISIs, float *SpikeCounts);

//extern void EvoFitVMN_CPU(int netmode, int numcells, float *chromepop, int gpuparams, int popsize, int blocksize, float runtime, float *Ints, float *ISIs, float *SpikeCounts);


#endif