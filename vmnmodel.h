/*
 *  vmnmodel.h
 *
 *  Created by Duncan MacGregor
 *  University of Edinburgh 2018
 *  Released under MIT license, see https://opensource.org/licenses/MIT
 *
 */


#ifndef VMNMOD_H
#define VMNMOD_H


#include "wx/wx.h"
#include <hypomodel.h>
#include "evofitbasic.h"


enum {
	ID_vmnflag = 3000,
	ID_datneuro,
	ID_gpunetmode,
	ID_cpumode,
	ID_autosum,
	ID_L1L2
};



class VMNModel;


class VMNNeuron : public NeuroDat{
public:
	double memtau, tauHAP, tauAHP, tauDAP;
	double kHAP, kAHP, kDAP;
	double th0, vrest;
	double input, inputdensity;
	double dend0e, dend0i;
	double dend1e, dend1i;
	double ae, ai, aesyn;
	double esynL1;
	
	double v, th;
	double ttime, neurotime;
	double tHAP, tAHP, tDAP; 
	double epspt0, ipspt0, epspt1, ipspt1;
	double synv, inputv;
	double esynsum, isynsum;
	double tB;
	
	double vrestsdgen;
	double kHAPsdgen;
	double tauHAPsdgen;
	double inputsdgen;
	double esyn, esynsdgen;
	
	datdouble inputrec;
	datdouble raterec;
	
	// Network connections
	int network[50];
	int enetwork[1000];
	int inetwork[1000];
	double eweight[1000];
	double iweight[1000];
	double esynqueue[100];
	double isynqueue[100];
	int econnect;
	int iconnect;
	int connect;
	int type;

	std::vector<double> mag;   // spike triggered signal magnitude, tB at corresponding spike time

	VMNNeuron();
};


class SignalBox: public ParamBox
{
public:
	VMNModel *mod;
	
	SignalBox(VMNModel *mod, const wxString& title, const wxPoint& pos, const wxSize& size);
	void OnRun(wxCommandEvent& event);
};


class VMNNeuroBox: public ParamBox
{
public:
	VMNModel *mod;

	VMNNeuroBox(VMNModel *mod, const wxString& title, const wxPoint& pos, const wxSize& size);
	
	void OnBox(wxCommandEvent& event);
	void PanelFull();
	void PanelBasic();
};


class VMNNetBox: public ParamBox
{
public:
	VMNNeuron *neurons;
	SpikeDat *netdat;
	VMNModel *mod;
	
	int vmhneurons;
	int layer1, layer2, layer3;
	int neuroindex;
	int sumflag;
	
	wxStaticText *spikes;
	wxStaticText *esyn;
	wxStaticText *vrest;
	wxStaticText *kHAP;
	wxStaticText *tauHAP;
	wxTextCtrl *datneuron;
	wxStaticText *freqL1, *freqL2, *freqL3;
	
	wxSpinButton *datspin;
	wxCheckBox *netcheck;
	wxCheckBox *cellcheck;
	wxCheckBox *seedcheck;
	wxCheckBox *unicheck;
	wxTextCtrl *storetag;
	//wxCheckBox *synccheck;        // moved to parent ParamBox
	wxToolBar *toolbar;
	wxPanel *toolpanel;
	
	VMNNetBox(VMNModel *, const wxString& title, const wxPoint& pos, const wxSize& size);
	
	void ModData(VMNNeuron *);
	void OnNeuroNext(wxSpinEvent& event);
	void OnNeuroPrev(wxSpinEvent& event);
	void OnSpin(wxSpinEvent& event);
	void OnRun(wxCommandEvent& event);
	void OnSum(wxCommandEvent& event);
	void OffSum();
	void NeuroData();
	void OnStore(wxCommandEvent& event);
	void OnLoad(wxCommandEvent& event);
	void OnParamStore(wxCommandEvent& event);
	void OnParamLoad(wxCommandEvent& event);
	void OnOutput(wxCommandEvent& event);
	void OnSynQueue(wxCommandEvent& event);
	void OnParamMenu(wxCommandEvent& event);
	void OnModelMenu(wxCommandEvent& event);
	void OnControlsMenu(wxCommandEvent& event);
	void netcalcvmh(SpikeDat *);
	void netcalcvmn(SpikeDat *, int layer=1);
	void SetNeuroCount();
	void NetworkSum();      // network spiking analysis, generate summed spike rates and ISI analysis 

	void PanelFull();
	void PanelBasic();
	
	virtual wxToolBar* OnCreateToolBar(long style, wxWindowID id, const wxString &name); 	
};


class VMNProtoBox : public ParamBox
{
public:
	VMNModel *mod;
	wxStaticText *currentinput;
	wxStaticText *currentpulse;
	wxStaticText *currentrange;
	wxStaticText *status;
	
	VMNProtoBox(VMNModel *model, const wxString& title, const wxPoint& pos, const wxSize& size);
	
	void OnRun(wxCommandEvent& event);
};


class VMNDat: public ModDat
{
public:
	VMNDat(int size);

	int runtime;
	datdouble input1;
	datdouble input2;
	datdouble noisig;
};


class VMNNeuroDat: public ModDat
{
public:
	VMNNeuroDat(int size);

	datdouble V;
	datdouble V2;
	datdouble HAP;
	datdouble AHP;
	datdouble netsyn;
	datdouble inpsyn;
	datdouble memsyn;
	datdouble memsyn2;
	datdouble sigsyn;
	datdouble sigsyn2;
	datdouble facB;
};


class VMHspikegen;


class VMNMod : public ModThread
{
public:
	double runtime;
	int datsample;
	int vmhneurons;
	int vmhL1;
	int vmhL2;
	int vmhL3;
	int parallel;
	int synchcount;
	int synchcall, spikecall;
	int synchrecord[10000];
	int synchOK;
	int synchflag[10];
	int done;
	int goflag[10];
	int queuelength;
	int synqueue;
	static const int numtypes = 3;

	datint *activity;
	int *active;
	FILE *ofp;
	wxString text;

	bool vmndiag;

	int spikecount[1000];
	wxMessageQueue <int> *synapse;
	
	VMNModel *mod;
	VMNNeuron *vmhneuron;
	VMNNeuroBox *neurobox;
	VMNNetBox *netbox;
	VMNProtoBox *protobox;
	ScaleBox *scalebox;
	
	// Model Flags
	bool revpots;
	bool iratioflag;

	// Protocol Flags
	bool rampflag;
	int prototype;

	// Analysis Flags
	int autosum;

	//Protocol Parameters
	int rampstart[numtypes], rampstop[numtypes];
	int rampbase[numtypes], rampinit[numtypes];
	double rampstep[numtypes];
	double rampinput[numtypes];

	int pulsestart[numtypes][5], pulsestop[numtypes][5];
	int pulsebase[numtypes][5], pulseinit[numtypes][5];
	double pulsestep[numtypes][5];
	double pulsetau[numtypes][5];
	double pulsemag[numtypes][5];
	double pulseinput[numtypes];

	int pulseprotocount;

	int rangestart, rangestop, rangestep;

	// Network Flags
	int cellgen;
	int netgen;
	
	// IGF parameters
	double numspikes;
	double hstep;
	double vthre[numtypes], vrest[numtypes];
	double ae[numtypes], ai[numtypes], ve, vi;
	double pspmag;
	double ire, iratio;
	double psphalflife;
	double kHAP[numtypes], kAHP[numtypes], kDAP[numtypes];
	double tauHAP[numtypes], tauAHP[numtypes], tauDAP[numtypes];
	double psptau;
	double emax;                // New March 2015
	
	// VMH Net parameters
	double vmhconnect;
	double vmhmaxcon;
	double vmhinput[numtypes];
	double esynweightL1;
	double esynweightL12;
	double esynweightL2;
	double esynweightL21;
	double isynweightL1;
	double isynweightL12;
	double isynweightL2;
	double isynweightL21;
	double maxsyn;
	double syndelay;
	double inputcycle;
	double absref; 
	double waveamp;
	double esynL1;
	double isynL1;
	double esynL12;
	double isynL12;
	double esynL2;
	double isynL2;
	double esynL21;
	double isynL21;
	double esyntrans;
	double isyntrans;
	double syntau;
	double synmag;
	double syndelrange;
	double esynsd;

	// 3rd Layer                                          17/7/14
	double esynL3, isynL3;
	double esynL23, isynL23;
	double esynweightL3, isynweightL3;
	double esynweightL23, isynweightL23;
	
	double esyn, isyn, vmhinput1, synweight;
	double aesyn, aisyn;
	
	double vrestsd[numtypes];
	double kHAPsd[numtypes];
	double tauHAPsd[numtypes];
	double inputsd[numtypes];

	// Noise Signal
	double noimean, noitau, noiamp;
	bool noisemode, signalmode;

	// Synaptic Wave Signal
	double synwaveamp;
	double synwavecycle;
	double synwaveshift;
	int siglayer;
	double sigIratio;

	// Output Signal Facilitation
	double kB;
	double halflifeB, tauB;
	
	// Model flags
	int unigen;
	int seedgen;
	long modseed;
	
	double input0, input1;
	double memtau;
	double th0;
	double re0, re1, ri0, ri1;
	double epsph, ipsph;
	
	int network[1000][50];
	int connect[1000];
	
	double dendinput[1000][2];
	
	double v_rest;
	
	VMNMod(VMNModel *);
	~VMNMod();
	virtual void *Entry();
	
	void spikegen(int, int, int *);
	void initialise();
	void networkgen();	
	void networkdisp();
	void networkgen2();	
	void networkdisp2();
	void inputoutput(int);
};




class VMNModel : public NeuroMod
{
public:
	int basicmode;

	VMNDat *vmndata;
	VMNNetBox *netbox;
	VMNNeuroBox *neurobox;
	SignalBox *signalbox;
	VMNProtoBox *protobox;
	//OutBox *outbox;
	//CellBox *cellbox;

	EvoFitBox *fitbox;
	EvoChrome *fitchrome;
	SpikeFitDat *spikefitdata;
	SpikeDat *fitboxdata;

	VMNNeuron *vmhneuron;
	SpikeDat *currvmn, *netdat;
	SpikeDat *netdat1, *netdat2, *netdat3;
	AnaDat *analysisdata;
	datdouble sigsim;
	datdouble raterec;
	datdouble raterec2;
	VMNNeuroDat *neurodata;

	datdouble rangeref;
	datdouble rangedata[3];
	double rangedata2[1000];

	int numtypes, prototypes;
	int datsample;
	int rangecount;

	VMNModel(int, wxString, HypoMain *);
	~VMNModel();

	void RunModel();
	void GraphData();
	void GSwitch(GraphDisp *gpos, ParamStore *gflags);
	void StoreClear();
	void SigSim(SpikeDat *);
	void RangeOut(int layers = 1);
	void Output();
	void EvoInit();
	void EvoRun();
};


#endif



