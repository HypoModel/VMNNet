/*
*  vmnmod.cpp
*  HypoModel
*
*  Created by Duncan MacGregor
*  University of Edinburgh 2018
*  Released under MIT license, see https://opensource.org/licenses/MIT
*
*
*/


#include "vmnmodel.h"


VMNMod::VMNMod(VMNModel *vmnmodel)
	: ModThread(vmnmodel->modbox, vmnmodel->mainwin) 
{
	mod = vmnmodel;
	datsample = mod->datsample;
	neurobox = mod->neurobox;
	netbox = mod->netbox;
	mainwin = mod->mainwin;
	scalebox = mod->mainwin->scalebox;
	protobox = mod->protobox;

	vmhneuron = mod->vmhneuron;

	vmndiag = false;
	queuelength = 100;

	active = new int[1000];
}


VMNMod::~VMNMod()
{
	delete []active;
}


void *VMNMod::Entry()
{
	int i, spiketotal = 0;
	int inprate, runcount, layers;

	initialise();

	layers = 1;
	if(vmhL2) layers = 2;
	if(vmhL3) layers = 3;

	if(prototype == range) {
		runcount = 0;
		for(inprate = rangestart; inprate<=rangestop; inprate+=rangestep) {
			protobox->currentrange->SetLabel(text.Format("%d", inprate));
			vmhinput[0] = inprate;
			spikegen(0, vmhneurons, active);
			if(vmhL1) netbox->netcalcvmn(mod->netdat1, 1);
			if(vmhL2) netbox->netcalcvmn(mod->netdat2, 2);
			if(vmhL3) netbox->netcalcvmn(mod->netdat3, 3);
			netbox->NeuroData();
			mod->rangeref[runcount] = inprate;
			if(vmhL1) mod->rangedata[0][runcount] = mod->netdat1->freq;
			if(vmhL2) mod->rangedata[1][runcount] = mod->netdat2->freq;
			if(vmhL3) mod->rangedata[2][runcount] = mod->netdat3->freq;
			if(vmhL1) netbox->freqL1->SetLabel(text.Format("%.2f", mod->netdat1->freq));
			if(vmhL2) netbox->freqL2->SetLabel(text.Format("%.2f", mod->netdat2->freq));
			if(vmhL3) netbox->freqL3->SetLabel(text.Format("%.2f", mod->netdat3->freq));
			runcount++;
		}
		mod->rangecount = runcount;
		mod->RangeOut(layers);
	}
	else {
		for(i=0; i<vmhneurons; i++) spikecount[i] = 0;
		spikegen(0, vmhneurons, active);
		if(autosum) netbox->NetworkSum();
		netbox->NeuroData();
	}

	return NULL;
}


void VMNMod::initialise()
{
	int i, celltype;
	FILE *ofp;

	//IGF parameters
	ParamStore *igfparams = neurobox->GetParams();
	numspikes = (*igfparams)["numspikes"];
	hstep = (*igfparams)["hstep"];
	vthre[0] = (*igfparams)["vthre"];
	vrest[0] = (*igfparams)["vrest"];
	//ve = (*igfparams)["ve"];
	//vi = (*igfparams)["vi"];
	pspmag = (*igfparams)["pspmag"];
	iratio = (*igfparams)["iratio"];
	ire = (*igfparams)["ire"];
	psphalflife = (*igfparams)["halflife"];
	kHAP[0] = (*igfparams)["kHAP"];	
	tauHAP[0] = 1 / (log((double)2) / (*igfparams)["halflifeHAP"]);
	kAHP[0] = (*igfparams)["kAHP"];
	tauAHP[0] = 1 / (log((double)2) / (*igfparams)["halflifeAHP"]);
	kDAP[0] = (*igfparams)["kDAP"];
	tauDAP[0] = 1 / (log((double)2) / (*igfparams)["halflifeDAP"]);
	psptau = log((double)2) / psphalflife;
	emax = (*igfparams)["emax"];
	absref = (*igfparams)["absref"];

	// VMH Network parameters
	ParamStore *netparams = netbox->GetParams();
	vmhL1 = (*netparams)["neuronsL1"];
	vmhL2 = (*netparams)["neuronsL2"];
	vmhL3 = (*netparams)["neuronsL3"];

	//vmhconnect = (*netparams)["vmhconnect"];
	vmhinput[0] = (*netparams)["vmhinput1"];

	esynweightL1 = (*netparams)["synweightL1"];
	isynweightL1 = esynweightL1;
	esynweightL2 = (*netparams)["synweightL2"];
	isynweightL2 = esynweightL2;
	esynweightL12 = (*netparams)["synweightL12"];
	isynweightL12 = esynweightL12;
	esynweightL21 = (*netparams)["synweightL21"];
	isynweightL21 = esynweightL21;
	esynweightL3 = (*netparams)["synweightL3"];
	esynweightL23 = (*netparams)["synweightL23"];

	
	inputcycle = (*netparams)["inputcycle"];
	waveamp = (*netparams)["waveamp"];
	syndelay = (*netparams)["syndelay"];
	syndelrange = (*netparams)["syndelrange"];

	esynL1 = (*netparams)["esynL1"];
	isynL1 = (*netparams)["isynL1"];
	esynsd = (*netparams)["esynsd"];
	esynL2 = (*netparams)["esynL2"];
	isynL2 = 0;
	esynL2sd = (*netparams)["esynL2sd"];
	esynL12 = (*netparams)["esynL12"];
	isynL12 = 0;
	esynL12sd = (*netparams)["esynL12sd"];
	esynL21 = (*netparams)["esynL21"];
	isynL21 = 0;
	esynL21sd = (*netparams)["esynL21sd"];

	vrestsd[0] = (*netparams)["vrestsd"];
	kHAPsd[0] = (*netparams)["kHAPsd"];
	tauHAPsd[0] = (*netparams)["tauHAPsd"];
	maxsyn = (*netparams)["maxsyn"];
	esyntrans = (*netparams)["esyntrans"];
	inputsd[0] = (*netparams)["inputsd"];

	vmhinput[1] = (*netparams)["vmhinput2"];
	vrest[1] = (*netparams)["vrest2"];
	kHAP[1] = (*netparams)["kHAP2"];
	tauHAP[1] = 1 / (log((double)2) / (*netparams)["halflifeHAP2"]);
	kDAP[1] = (*netparams)["kDAP2"];
	tauDAP[1] = 1 / (log((double)2) / (*netparams)["halflifeDAP2"]);
	tauHAPsd[1] = (*netparams)["tauHAP2sd"];

	syntau = 1 / (log((double)2) / (*netparams)["synhl"]);
	synmag = (*netparams)["synmag"];

	vmhinput[2] = (*netparams)["vmhinput3"];
	vrest[2] = (*netparams)["vrest3"];
	kHAP[2] = (*netparams)["kHAP3"];
	tauHAP[2] = 1 / (log((double)2) / (*netparams)["halflifeHAP3"]);
	esynL3 = (*netparams)["esynL3"];
	esynL23 = (*netparams)["esynL23"];

	vthre[2] = vthre[0];
	kHAPsd[2] = 0;
	tauHAPsd[2] = 0;
	vrestsd[2] = 0;
	inputsd[2] = 0;
	kAHP[2] = kAHP[0];
	tauAHP[2] = tauAHP[0];
	kDAP[2] = 0;
	tauDAP[2] = 1 / (log((double)2) / (*netparams)["halflifeDAP2"]);

	vthre[1] = vthre[0];
	kAHP[1] = kAHP[0];
	tauAHP[1] = tauAHP[0];
	kHAPsd[1] = 0;
	//tauHAPsd[1] = 0;
	vrestsd[1] = 0;
	inputsd[1] = inputsd[0];

	if(mod->basicmode) {
		emax = 0;
		maxsyn = 0;
		vmhL3 = 0;
		vrestsd[0] = 0;
		kHAPsd[0] = 0;
		tauHAPsd[0] = 0;
		tauHAPsd[1] = 0;
		inputsd[0] = 0;
		inputcycle = 0;
	}

	if(mod->revisionmode) {
		emax = (*igfparams)["emax"];
	}


	vmhneurons = vmhL1 + vmhL2 + vmhL3;

	// Flags
	netgen = (*netbox->modflags)["netgen"];
	cellgen = (*netbox->modflags)["cellgen"];
	cellgen2 = (*netbox->modflags)["cellgen"];
	unigen = (*netbox->modflags)["unigen"];
	seedgen = (*netbox->modflags)["seedgen"];

	synqueue = (*netbox->modflags)["synqueue"];
	fixeddelay = (*netbox->modflags)["fixeddelay"];

	revpots = (*neurobox->modflags)["revpots"];
	iratioflag = (*neurobox->modflags)["Iratioflag"];

	autosum = (*netbox->conflags)["autosum"];
	modseed = (*netparams)["modseed"];

	// Signal Parameters  - input signal and output signal simulation
	ParamStore *sigparams = mod->signalbox->GetParams();
	
	// Noise Parameters
	noimean = (*sigparams)["noimean"];
	noitau = (*sigparams)["noitau"];
	noiamp = (*sigparams)["noiamp"];

	noisemode = (*mod->signalbox->modflags)["noiseflag"];
	siglayer = (*mod->signalbox->modflags)["layerflag"];
	signalmode = noisemode;    // signalmode replacing noisemode

	sigIratio = (*sigparams)["sigIratio"];

	// Wave Parameters
	synwaveamp = (*sigparams)["synwaveamp"];
	synwavecycle = (*sigparams)["synwavecycle"];
	synwaveshift = (*sigparams)["synwaveshift"];

	// Facilitation Parameters
	kB = (*sigparams)["kB"];
	halflifeB = (*sigparams)["halflifeB"];
	tauB = 1 / (log((double)2) / halflifeB);

	if(mod->basicmode) {
		synwaveamp = 0;
		kB = 0;
	}

	mod->diagbox->Write(text.Format("Signal kB %.4f halflifeB %.2f tauB %.4f\n", kB, halflifeB, tauB));


	// Protocol Parameters

	wxString tag[2];
	prototype = (*mod->modeflags)["prototype"];

	ParamStore *protoparams = protobox->GetParams();

	for(i=0; i<mod->prototypes; i++) {
		tag[i].Printf("%d", i);
		rampbase[i] = (*protoparams)["rampbase" + tag[i]];
		rampstart[i] = (*protoparams)["rampstart" + tag[i]];
		rampstop[i] = (*protoparams)["rampstop" + tag[i]];
		rampinit[i] = (*protoparams)["rampinit" + tag[i]];
		rampstep[i] = (*protoparams)["rampstep" + tag[i]];

		mainwin->diagbox->Write(text.Format("ramp proto %d  base %d  step %.4f\n", i, rampbase[i], rampstep[i]));
	}

	pulseprotocount = 2;
	celltype = 0;

	for(i=0; i<pulseprotocount; i++) {
		tag[i].Printf("%d", i);
		pulsebase[celltype][i] = (*protoparams)["pulsebase" + tag[i]];
		pulsestart[celltype][i] = (*protoparams)["pulsestart" + tag[i]];
		pulsestop[celltype][i] = (*protoparams)["pulsestop" + tag[i]];
		pulseinit[celltype][i] = (*protoparams)["pulseinit" + tag[i]];
		pulsetau[celltype][i] = 1 / (log((double)2) / ((*protoparams)["pulsehl" + tag[i]]) / 1000);
	}

	rangestart = (*protoparams)["rangestart"];
	rangestop = (*protoparams)["rangestop"];
	rangestep = (*protoparams)["rangestep"];

	// Synch Recording

	for(i=0; i<10000; i++) synchrecord[i] = 0;

	// Network generation

	if(seedgen) modseed = (unsigned)(time(NULL));
	netbox->paramset->GetCon("modseed")->SetValue(modseed);
	init_mrand(modseed);
	
	if(netgen) networkgen2();
	if(vmndiag) networkdisp2();

	// Log Parameters

	if(vmndiag) {
		ofp = fopen("vmhnetparams.txt", "w");
		fprintf(ofp, "Model run parameters\n\n");
		fprintf(ofp, "Generating spikes on %d neurons, %.0fms\n\n", vmhneurons, numspikes*1000);
		fprintf(ofp, "HAP k = %.2f  tau = %.4f\n", kHAP[0], tauHAP[0]);
		fprintf(ofp, "AHP k = %.2f  tau = %.4f\n", kAHP[0], tauAHP[0]);
		fprintf(ofp, "DAP k = %.2f  tau = %.4f\n", kDAP[0], tauDAP[0]);
		fprintf(ofp, "EPSP rate = %.2f   ISP ratio = %.2f\n", ire, iratio);
		fprintf(ofp, "vrest = %.2f, th0 = %.2f\n", vrest[0], vthre[0]); 
		fprintf(ofp, "ae %.2f  ai %.2f  ve %.2f  vi %.2f\n\n", ae, ai, ve, vi);
		fprintf(ofp, "vmhconnect %.2f synweight %.2f\n", vmhconnect, esynweightL2);
		fprintf(ofp, "esynL2 %.2f  esynweightL2 %.2f\n", esynL2, esynweightL2);
		fprintf(ofp, "esynL12 %.2f  esynweightL12 %.2f\n", esynL12, esynweightL12);
		fprintf(ofp, "esynL21 %.2f  esynweightL21 %.2f\n", esynL21, esynweightL21);
		fprintf(ofp, "absref %.2f %.2f synqueue %d\n", absref, syndelay, synqueue);
		fprintf(ofp, "netgen %d  cellgen %d\n", netgen, cellgen);
		fclose(ofp);
	}
}

// Outdated diagnostic function
void VMNMod::inputoutput(int n)
{
	int i;
	FILE *ofp;

	ofp = fopen("vmhinputtrace.txt", "w");
	fprintf(ofp, "VMH network neuron %d\n\n", n);

	//for(i=0; i<10000; i++) {
	//	fprintf(ofp, "%d ms : %.2f mV\n", i, vmhneuron[n].inputrec[i]);
	//}
	fclose(ofp);
}


void VMNMod::networkdisp()
{
	int i, b, con;
	FILE *ofp;

	ofp = fopen("vmhnet.txt", "w");

	fprintf(ofp, "VMH Network\n\n");
	fprintf(ofp, "%d neurons, %.0f synapses\n\n", vmhneurons, vmhconnect);
	
	for(i=0; i<vmhneurons; i++) {
		for(con = 0; con<connect[i]; con++) fprintf(ofp, "%d ", network[i][con]);
		fprintf(ofp, "\n");
	}		
	fclose(ofp);
}


void VMNMod::networkdisp2()
{
	int i, c;
	FILE *ofp;

	ofp = fopen("vmhnet2.txt", "w");

	fprintf(ofp, "VMH Network\n\n");
	fprintf(ofp, "%d neurons, esyn %.2f isyn %.2f\n\n", vmhneurons, esyn, isyn);

	for(i=0; i<vmhneurons; i++) {
		fprintf(ofp, "N%d Econ: ", i);
		for(c = 0; c<vmhneuron[i].econnect; c++) fprintf(ofp, "%d ", vmhneuron[i].enetwork[c]);
		fprintf(ofp, " Icon: ");
		for(c = 0; c<vmhneuron[i].iconnect; c++) fprintf(ofp, "%d ", vmhneuron[i].inetwork[c]);
		fprintf(ofp, "\n");
	}

	fprintf(ofp, "\n\n");
	for(i=0; i<vmhneurons; i++) {
		fprintf(ofp, "N%d Eweight: ", i);
		for(c = 0; c<vmhneuron[i].econnect; c++) fprintf(ofp, "%.1f ", vmhneuron[i].eweight[c]);
		fprintf(ofp, "\n");
	}

	fclose(ofp);
}


void VMNMod::networkgen()           // Out Of Use
{
	int i, con, syn, n;
	double d;
	int pregrid[1000];

	for(i=0; i<vmhneurons; i++) {
		for(con=0; con<vmhmaxcon; con++) network[i][con] = -1;
		connect[i] = vmhconnect;
		for(n=0; n<vmhneurons; n++) pregrid[n] = 0;
		pregrid[i] = 1;
		for(con=0; con<connect[i]; con++) {
			d = mrand01();
			syn = floor(d * vmhneurons);
			while(pregrid[syn] > 0) {
				d = mrand01();
				syn = floor(d * vmhneurons);
			}
			pregrid[syn]++;
			network[i][con] = syn;
		}
	}	
}


void VMNMod::networkgen2()
{
	int i, j;
	double d;
	int startL1, endL1;
	int startL2, endL2;
	int startL3, endL3;
	wxString text;

	startL1 = 0;
	endL1 = vmhL1;
	startL2 = vmhL1;
	endL2 = vmhL1 + vmhL2;
	startL3 = vmhL1 + vmhL2;
	endL3 = vmhL1 + vmhL2 + vmhL3;

	// Layer1

	for(i=0; i<vmhL1; i++) {
		vmhneuron[i].econnect = 0;
		vmhneuron[i].iconnect = 0;

		if(cellgen) vmhneuron[i].esynsdgen = gaussian(0, 1); 
		else vmhneuron[i].esynsdgen = 0;
		vmhneuron[i].esynL1 = esynL1 + vmhneuron[i].esynsdgen * esynsd;

		//mod->diagbox->Write(text.Format("netgen neuron %d, esyn = %.2f\n", i, vmhneuron[i].esyn));

		for(j=0; j<vmhL1; j++) {
			d = mrand01();
			if(d <= vmhneuron[i].esynL1 && i != j) {
				vmhneuron[i].eweight[vmhneuron[i].econnect] = esynweightL1;
				vmhneuron[i].edelay[vmhneuron[i].econnect] = (syndelrange + 1) * mrand01();   // new November 2018 - fixed connection delay
				vmhneuron[i].enetwork[vmhneuron[i].econnect++] = j;

			}			
			d = mrand01();
			if(d <= isynL1 && i != j) {
				vmhneuron[i].iweight[vmhneuron[i].iconnect] = isynweightL1;
				vmhneuron[i].inetwork[vmhneuron[i].iconnect++] = j;
			}
		}

		if(cellgen2) vmhneuron[i].esynL21sdgen = gaussian(0, 1); 
		else vmhneuron[i].esynL21sdgen = 0;
		vmhneuron[i].esynL21 = esynL21 + vmhneuron[i].esynL21sdgen * esynL21sd;

		for(j=vmhL1; j<vmhL1+vmhL2; j++) {
			d = mrand01();
			if(d <= vmhneuron[i].esynL21 && i != j) {
				vmhneuron[i].eweight[vmhneuron[i].econnect] = esynweightL21;
				vmhneuron[i].edelay[vmhneuron[i].econnect] = (syndelrange + 1) * mrand01();   // new November 2018 - fixed connection delay
				vmhneuron[i].enetwork[vmhneuron[i].econnect++] = j;
			}
			d = mrand01();
			if(d <= isynL21 && i != j) {
				vmhneuron[i].iweight[vmhneuron[i].iconnect] = isynweightL21;
				vmhneuron[i].inetwork[vmhneuron[i].iconnect++] = j;
			}
		}
	}


	// Layer 2

	for(i=vmhL1; i<vmhL1+vmhL2; i++) {
		vmhneuron[i].econnect = 0;
		vmhneuron[i].iconnect = 0;

		if(cellgen2) vmhneuron[i].esynL12sdgen = gaussian(0, 1); 
		else vmhneuron[i].esynL12sdgen = 0;
		vmhneuron[i].esynL12 = esynL12 + vmhneuron[i].esynL12sdgen * esynL12sd;

		for(j=0; j<vmhL1; j++) {
			d = mrand01();
			if(d <= vmhneuron[i].esynL12 && i != j) {
				vmhneuron[i].eweight[vmhneuron[i].econnect] = esynweightL12;
				vmhneuron[i].edelay[vmhneuron[i].econnect] = (syndelrange + 1) * mrand01();   // new November 2018 - fixed connection delay
				vmhneuron[i].enetwork[vmhneuron[i].econnect++] = j;
			}
			d = mrand01();
			if(d <= isynL12 && i != j) {
				vmhneuron[i].iweight[vmhneuron[i].iconnect] = isynweightL12;
				vmhneuron[i].inetwork[vmhneuron[i].iconnect++] = j;
			}
		}

		if(cellgen2) vmhneuron[i].esynL2sdgen = gaussian(0, 1); 
		else vmhneuron[i].esynL2sdgen = 0;
		vmhneuron[i].esynL2 = esynL2 + vmhneuron[i].esynL2sdgen * esynL2sd;

		for(j=vmhL1; j<vmhL1+vmhL2; j++) {
			d = mrand01();
			if(d <= vmhneuron[i].esynL2 && i != j) {
				vmhneuron[i].eweight[vmhneuron[i].econnect] = esynweightL2;
				vmhneuron[i].edelay[vmhneuron[i].econnect] = (syndelrange + 1) * mrand01();   // new November 2018 - fixed connection delay
				vmhneuron[i].enetwork[vmhneuron[i].econnect++] = j;
			}
			d = mrand01();
			if(d <= isynL2 && i != j) {
				vmhneuron[i].iweight[vmhneuron[i].iconnect] = isynweightL2;
				vmhneuron[i].inetwork[vmhneuron[i].iconnect++] = j;
			}
		}
	}


	// Layer 3

	for(i=startL3; i<endL3; i++) {
		vmhneuron[i].econnect = 0;
		vmhneuron[i].iconnect = 0;

		for(j=startL2; j<endL2; j++) {
			d = mrand01();
			if(d <= esynL23 && i != j) {
				vmhneuron[i].eweight[vmhneuron[i].econnect] = esynweightL23;
				vmhneuron[i].enetwork[vmhneuron[i].econnect++] = j;
			}
			d = mrand01();
			if(d <= isynL23 && i != j) {
				vmhneuron[i].iweight[vmhneuron[i].iconnect] = isynweightL23;
				vmhneuron[i].inetwork[vmhneuron[i].iconnect++] = j;
			}
		}

		for(j=startL3; j<endL3; j++) {
			d = mrand01();
			if(d <= esynL3 && i != j) {
				vmhneuron[i].eweight[vmhneuron[i].econnect] = esynweightL3;
				vmhneuron[i].enetwork[vmhneuron[i].econnect++] = j;
			}
			d = mrand01();
			if(d <= isynL3 && i != j) {
				vmhneuron[i].iweight[vmhneuron[i].iconnect] = isynweightL3;
				vmhneuron[i].inetwork[vmhneuron[i].iconnect++] = j;
			}
		}
	}

}


void VMNMod::spikegen(int nstart, int nstop, int *activity)
{			
	int i, j, tmax = 0, n, tstep, c;
	int nepsp0, nipsp0, nepsp1, nipsp1;
	int numconnect, con;
	int modstep;
	int modeltime;
	int syndel;
	int celltype, pulseproto;
	int maxcells = 500;
	int maxspikes = 100000;

	int s[1000];
	double maxtime, nettime;
	double maxdap = 50;

	double ire1 = ire;
	double iratio1 = iratio;
	double dend0e, dend0i, dend1e, dend1i;
	double connectsum;
	double dendsum0, dendsum1, endosum;
	double esynsum, isynsum;
	double esyninput, isyninput, sinput;
	double suckstim;
	double waveinput, noiseinput;
	double noisig, wavesig, synsig;
	double pi = 3.14159265;
	double esynh;
	double synrand;

	float **weights = new float*[maxcells];
	for(i=0; i<maxcells; i++) weights[i] = new float[maxcells];

	float **esynqueue = new float*[maxcells];
	for(i=0; i<maxcells; i++) esynqueue[i] = new float[queuelength];

	int **enetwork = new int*[maxcells];
	for(i=0; i<maxcells; i++) enetwork[i] = new int[maxcells];

	float **delays = new float*[maxcells];
	for(i=0; i<maxcells; i++) delays[i] = new float[maxcells];

	FILE *ofp;
	FILE *tofp, *inputofp;
	if(vmndiag) tofp = fopen("vmhtest.txt", "w");
	if(vmndiag) inputofp = fopen("vmninput.txt", "w");

	/// Param copy for testing

	th0 = vthre[0];
	v_rest = vrest[0];

	maxtime = numspikes * 1000;

	memtau = 1 / (log((double)2) / psphalflife);

	dend0e = 100; //vmhinput[0] / 1000;       // out of use - see vmhneuron.dend0e
	dend0i = dend0e * iratio;                 //

	dend1e = 100; //vmhinput[1] / 1000;       // out of use
	dend1i = dend1e * iratio1;                //

	if(vmndiag) {
		ofp = fopen("vmhdatanet.txt", "w");
		fprintf(ofp, "Model run parameters\n\n");
		fprintf(ofp, "Generating spikes on %d neurons, %.0fms\n\n", vmhneurons, maxtime);
		fprintf(ofp, "Type 0 %d neurons, Type 1 %d neurons\n", vmhL1, vmhL2);
		fprintf(ofp, "Type 0 HAP k = %.2f  tau = %.4f\n", kHAP[0], tauHAP[0]);
		fprintf(ofp, "Type 1 HAP k = %.2f  tau = %.4f\n", kHAP[1], tauHAP[1]);
		fprintf(ofp, "AHP k = %.2f  tau = %.4f\n", kAHP[0], tauAHP[0]);
		fprintf(ofp, "DAP k = %.2f  tau = %.4f\n", kDAP[0], tauDAP[0]);
		fprintf(ofp, "memtau = %.2f\n", memtau);
		fprintf(ofp, "EPSP rate = %.2f   ISP ratio = %.2f\n", ire, iratio);
		//fprintf(ofp, "revpots = %d  newVth = %d   memtau = %.4f\n", revpots, newVth, memtau);
		//fprintf(ofp, "dend0e = %.2f, dend0i = %.2f, dend1e = %.2f, dend1i = %.2f\n", dend0e, dend0i, dend1e, dend1i);
		fprintf(ofp, "vrest0 = %.2f, thresh0 = %.2f\n", vrest[0], vthre[0]); 
		fprintf(ofp, "vrest1 = %.2f, thresh1 = %.2f\n", vrest[1], vthre[1]); 
		fprintf(ofp, "input0 = %.2f, input1 = %.2f\n", vmhinput[0], vmhinput[1]); 
		//fprintf(ofp, "ae %.2f  ai %.2f  ve %.2f  vi %.2f\n\n", ae, ai, ve, vi);
		fprintf(ofp, "ve %.2f  vi %.2f\n\n", ve, vi);
		fprintf(ofp, "vmhconnect %.2f synweight %.2f\n", vmhconnect, esynweightL1);
		fprintf(ofp, "absref %.2f syndelay %.2f\n", absref, syndelay);
		fprintf(ofp, "noimean %.2f  noitau %.2f  noiamp %.2f\n", noimean, noitau, noiamp);
		fprintf(ofp, "synwaveamp %.2f  synwavecycle %.2f  synwaveshift %.2f\n", synwaveamp, synwavecycle, synwaveshift);
		fclose(ofp);
	}


	// Protocol
	double rampstep1ms[2]; 

	rampstep1ms[0] = (double)rampstep[0] / 1000;
	rampstep1ms[1] = (double)rampstep[1] / 1000;

	celltype = 0;
	pulsemag[celltype][0] = pulseinit[celltype][0];
	pulsemag[celltype][1] = pulseinit[celltype][1];


	// initial values

	for(i=0; i<vmhneurons; i++) {

		// Select and set neuron type
		if(i < vmhL1) {
			celltype = 0; 
			vmhneuron[i].type = 0;
		}
		if(i >= vmhL1 && i < vmhL1 + vmhL2) {
			celltype = 1;
			vmhneuron[i].type = 1;
		}
		if(i >= vmhL1 + vmhL2) {
			celltype = 2;
			vmhneuron[i].type = 2;
		}

		double min_tauHAP = 5;

		// Random parameter generation
		if(cellgen) {
			vmhneuron[i].vrestsdgen = gaussian(0, 1);
			if(unigen) vmhneuron[i].vrestsdgen = mrand01();
			vmhneuron[i].kHAPsdgen = gaussian(0, 1);
			vmhneuron[i].tauHAPsdgen = gaussian(0, 1);
			vmhneuron[i].inputsdgen = gaussian(0, 1);
		}

		vmhneuron[i].vrest = vrest[celltype] + vmhneuron[i].vrestsdgen * vrestsd[celltype];
		
		if(unigen) vmhneuron[i].vrest = (vrest[celltype] - vrestsd[celltype]) + (vrestsd[celltype] * 2 * vmhneuron[i].vrestsdgen);

		vmhneuron[i].kHAP = kHAP[celltype] + vmhneuron[i].kHAPsdgen * kHAPsd[celltype];
		vmhneuron[i].tauHAP = tauHAP[celltype] + vmhneuron[i].tauHAPsdgen * tauHAPsd[celltype];
		if(vmhneuron[i].tauHAP < min_tauHAP) vmhneuron[i].tauHAP = min_tauHAP;      // new April 2019

		//vmhneuron[i].input = vmhinput[type] + vmhneuron[i].inputsdgen * inputsd[type];
		vmhneuron[i].inputdensity = 1 + vmhneuron[i].inputsdgen * inputsd[celltype];            // new Feb 2014

		vmhneuron[i].input = vmhinput[celltype] * vmhneuron[i].inputdensity;                    // new Feb 2014
		vmhneuron[i].dend0e = vmhneuron[i].input / 1000;         // currently only used for epspt initiation
		vmhneuron[i].dend0i = vmhneuron[i].dend0e * iratio;

		vmhneuron[i].ae = pspmag / (ve - vmhneuron[i].vrest);
		vmhneuron[i].ai = pspmag / (vmhneuron[i].vrest - vi);
		vmhneuron[i].aesyn = synmag / (ve - vmhneuron[i].vrest);
		vmhneuron[i].inputv = 0;
		vmhneuron[i].synv = 0;

		//vmhneuron[i].vrest = gaussian(vrest, vrestsd);
		//vmhneuron[i].vrest = uniform(vrest, vrestsd);
		vmhneuron[i].th0 = gaussian(vthre[celltype], 0);
		//vmhneuron[i].kHAP = gaussian(kHAP, kHAPsd);
		vmhneuron[i].kAHP = gaussian(kAHP[celltype], 0);
		vmhneuron[i].kDAP = gaussian(kDAP[celltype], 0);
		//vmhneuron[i].tauHAP = gaussian(tauHAP, tauHAPsd);
		vmhneuron[i].tauAHP = gaussian(tauAHP[celltype], 0);
		vmhneuron[i].tauDAP = gaussian(tauDAP[celltype], 0);

		vmhneuron[i].ttime = 0;

		vmhneuron[i].v = vmhneuron[i].vrest;
		vmhneuron[i].tHAP = vmhneuron[i].kHAP;
		vmhneuron[i].tAHP = vmhneuron[i].kAHP;
		vmhneuron[i].tDAP = vmhneuron[i].kDAP;
		vmhneuron[i].epspt0 = -log(1 - mrand01()) / dend0e;
		vmhneuron[i].ipspt0 = -log(1 - mrand01()) / dend0i;
		vmhneuron[i].epspt1 = -log(1 - mrand01()) / dend1e;
		vmhneuron[i].ipspt1 = -log(1 - mrand01()) / dend1i;
		s[i] = 0;
		vmhneuron[i].neurotime = 0;
		activity[i] = 0;

		vmhneuron[i].tB = 0;
		
		if(celltype == 1) {
			vmhneuron[i].inputrec.data.resize(110000);
			vmhneuron[i].inputrec.max = 100000;
		}
		if(celltype == 0) {
			vmhneuron[i].inputrec.data.resize(110000);
			vmhneuron[i].inputrec.max = 100000;
		}

		for(c=0; c<vmhneuron[i].econnect; c++) {            // copy to local arrays for speed
			weights[i][c] = vmhneuron[i].eweight[c];
			delays[i][c] = vmhneuron[i].edelay[c];
			enetwork[i][c] = vmhneuron[i].enetwork[c];
		}	
		for(j=0; j<queuelength; j++) esynqueue[i][j] = 0;
	}

	noisig = noimean;

	modeltime = 0;
	nettime = 0;
	modstep = 0;

	if(vmndiag) fprintf(tofp, "Start spike gen, run for %.0fms\n", maxtime);

	i = 0;
	for(i=0; i<vmhneurons; i++) {
		vmhneuron[i].th = vmhneuron[i].th0 + vmhneuron[i].tHAP + vmhneuron[i].tAHP - vmhneuron[i].tDAP;
		if(vmndiag) {
			fprintf(tofp, "neuron %d  time %.1f  v %.2f  vrest %.2f  vthresh %.2f\n\n", i, nettime, vmhneuron[i].v, vmhneuron[i].vrest, vmhneuron[i].th); 
			fprintf(tofp, "inputv %.2f  synv %.2f\n", vmhneuron[i].inputv, vmhneuron[i].synv);
			/*fprintf(tofp, "neuron %d  time %.1f  input0 %.2f  input1 %.2f  sinput %.2f  v %.2f  thresh %.2f\n", 
					i, nettime, input0, input1, sinput, vmhneuron[i].v, vmhneuron[i].th); 

			//fprintf(tofp, "dendinput %.2f  sinput %.2f\n", dendinput[i][1], sinput);
			fprintf(tofp, "esyninput %.2f  esynqueue %.2f isyninput %.2f\n", esyninput, esynqueue[i][0], isyninput);
			fprintf(tofp, "inputv %.2f  synv %.2f\n", vmhneuron[i].inputv, vmhneuron[i].synv);
			fprintf(tofp, "nepsp %d  epsph %.2f  nipsp %d  ipsph %.2f\n", nepsp0, epsph, nipsp0, ipsph);
			fprintf(tofp, "dend1e %.2f  dend1i %.2f\n\n", vmhneuron[i].dend1e, vmhneuron[i].dend1i);*/
		}
		//fprintf(tofp, "dendinput %.2f  sinput %.2f\n", dendinput[i][1], sinput);
		//fprintf(tofp, "sinput %.2f\n", sinput);
		//fprintf(tofp, "nepsp %d  epsph %.2f  nipsp %d  ipsph %.2f\n\n", nepsp0, epsph, nipsp0, ipsph);
	}


	////////////////////////

	mod->diagbox->Write(text.Format("\nprototype %d\n", prototype)); 

	for(pulseproto=0; pulseproto<pulseprotocount; pulseproto++) {
		mod->diagbox->Write(text.Format("Type 0 Proto %d Start %d Stop %d Step %.2f\n", pulseproto, pulsestart[0][pulseproto], pulsestop[0][pulseproto], pulsemag[0][pulseproto]));
	}


	netbox->countmark = 0;

	while(nettime<maxtime) {

		modstep++;
		nettime = nettime + hstep;
		if((int)nettime%1000 == 0) netbox->SetCount(floor(nettime)/maxtime*100);
		//if(int(runtime / 10) % (int)tmax == 0) parambox->SetCount(runtime / 10 / tmax);

		//fprintf(tofp, "time %.2f\n", nettime);

		// Protocol /////////////////////////


		// Noise Input         new 9/3/17

		noisig = noisig + hstep * ((noimean - noisig) / noitau) + noiamp * sqrt(hstep) * gaussian(0, 1); 

		// Wave Signal (synaptic)          new 13/3/18

		wavesig = 0;
		if(synwavecycle) wavesig = synwaveamp * 0.5 * (1 + sin((nettime + synwaveshift) * 2 * pi / synwavecycle)); 

		synsig = noisig + wavesig;


		// double noisediff = (noimean - noise) * noidecay + noiamp * sqrt(stepsize) * gaussian(0, 1);

		// Wave Input

		waveinput = 0;
		if(inputcycle) waveinput = waveamp * 0.5 * (1 + sin(nettime*2*pi/inputcycle));

		// Ramp Input

		//mod->diagbox->Write(text.Format("\n proto types %d\n", mod->numtypes)); 

		for(celltype=0; celltype<mod->prototypes; celltype++) {
			rampinput[celltype] = 0;

			//mod->diagbox->Write(text.Format("\nprototype %d\n", prototype)); 

			if(prototype == ramp) {
				//mod->diagbox->Write("set ramp\n"); 
				if(nettime < rampstart[celltype]*1000) rampinput[celltype] = rampbase[celltype];
				if(nettime >= rampstart[celltype]*1000 && nettime < rampstop[celltype]*1000) 
					rampinput[celltype] = rampinit[celltype] + (nettime - rampstart[celltype]*1000) * rampstep1ms[celltype] * hstep;
				if(nettime >= rampstop[celltype]*1000) rampinput[celltype] = rampbase[celltype];
				vmhinput[celltype] = rampinput[celltype];
				if(vmhinput[celltype] < 0) vmhinput[celltype] = 0;
			}
		}

		
		for(celltype=0; celltype<mod->prototypes; celltype++) pulseinput[celltype] = 0;

		celltype = 0;
		if(prototype == pulse) {
			for(pulseproto=0; pulseproto<pulseprotocount; pulseproto++) {
				//if(nettime < pulsestart[celltype][pulseproto]*1000) pulseinput[celltype] = pulsebase[celltype][pulseproto];
				if(nettime >= pulsestart[celltype][pulseproto] * 1000 && nettime < pulsestop[celltype][pulseproto] * 1000) { 
					pulseinput[celltype] = pulsemag[celltype][pulseproto];
					if(pulsetau[celltype][pulseproto] > 1) pulsemag[celltype][pulseproto] = pulsemag[celltype][pulseproto] - hstep * pulsemag[celltype][pulseproto] / pulsetau[celltype][pulseproto]; 
				}
				if(nettime >= pulsestop[celltype][pulseproto] * 1000) pulseinput[celltype] = 0;
				//vmhinput[celltype] = pulseinput[celltype];
				//if(pulseinput[celltype] < 0) pulseinput[celltype] = 0;
			}
			if(modstep%1000 == 0 && pulseinput[celltype] > 0) mod->diagbox->Write(text.Format("Type 0 Pulse Input %.2f\n", pulseinput[celltype]));
		}

		//fprintf(tofp, "Inputs: neuron %d vmhconnect = %.2f synweight = %.2f dendsum1 = %.2f\n", i, vmhconnect, synweight, dendsum1);

		// Monitor

		if(modstep%mod->datsample == 0 && modstep<mod->storesize*datsample) {
			mod->vmndata->input1[modstep/datsample] = vmhinput[0];
			mod->vmndata->input2[modstep/datsample] = vmhinput[1];
			mod->vmndata->noisig[modstep/datsample] = synsig;
		}


		// Network input

		for(i=0; i<vmhneurons; i++) {

			vmhneuron[i].ttime = vmhneuron[i].ttime + hstep;
			vmhneuron[i].neurotime = vmhneuron[i].neurotime + hstep;

			esynsum = 0;

			if(!synqueue) {
				for(c=0; c<vmhneuron[i].econnect; c++) 
					if(activity[enetwork[i][c]] == 1) {
						esynsum = esynsum + weights[i][c];
					}
					vmhneuron[i].esynsum = esynsum;
			}
			else {
				for(c=0; c<vmhneuron[i].econnect; c++) 
					if(activity[enetwork[i][c]] == 1) {
						synrand = mrand01();
						if(esyntrans >= synrand) { 
							if(!fixeddelay) syndel = (syndelay - 1) + (syndelrange + 1) * (synrand * (1/esyntrans));    // scaled use of synrand (max value = esyntrans) allows second use for random delay
							else syndel = (syndelay - 1) + delays[i][c];
							syndel = (double)syndel/hstep;
							esynqueue[i][syndel] = esynqueue[i][syndel] + weights[i][c];	
							if(vmndiag && i == 0) fprintf(tofp, "time %.1f connect %d syndel %d\n", nettime, c, syndel);
						}
					}
			}
			/*else {
				for(c=0; c<vmhneuron[i].econnect; c++) 
					if(activity[enetwork[i][c]] == 1) {
						synrand = mrand01();
						if(esyntrans >= synrand) { 
							if(!fixeddelay) syndel = (syndelay - 1) + (syndelrange + 1) * (synrand * (1/esyntrans));    // scaled use of synrand (max value = esyntrans) allows second use for random delay
							else syndel = (syndelay - 1) + delays[i][c];
							esynqueue[i][syndel] = esynqueue[i][syndel] + weights[i][c];	
							if(vmndiag && i == 0) fprintf(tofp, "time %.1f connect %d syndel %d\n", nettime, c, syndel);
						}
					}
			}*/
		}

		// Spiking model

		for(i=0; i<vmhneurons; i++) {

			// Signal input

			celltype = vmhneuron[i].type;
			vmhneuron[i].dend0e = (vmhneuron[i].inputdensity * (vmhinput[celltype] + pulseinput[celltype])) / 1000;
			if(iratioflag) vmhneuron[i].dend0i = vmhneuron[i].dend0e * iratio;

			if(signalmode) {   
				//mod->diagbox->Write("noise mode on\n");// independent noise signal   February 2018
				if(celltype == siglayer) vmhneuron[i].dend1e = synsig / 1000; 
				else vmhneuron[i].dend1e = 0;
				vmhneuron[i].dend1i = vmhneuron[i].dend1e * sigIratio; //* iratio1;
				//mod->diagbox->Write("dend1e set\n");
			}

			if(i == 0) mod->raterec[modstep] = vmhneuron[i].dend0e * 1000;
			if(i == vmhL1) mod->raterec2[modstep] = vmhneuron[i].dend0e * 1000;

			if(activity[i]) activity[i] = activity[i] - 1;

			if(synqueue) {
				vmhneuron[i].esynsum = esynqueue[i][0];
				for(j=0; j<queuelength-1; j++) esynqueue[i][j] = esynqueue[i][j+1];
				esynqueue[i][queuelength-1] = 0;
			}

			isynsum = 0;
		
			esyninput = vmhneuron[i].esynsum;
			//isyninput = vmhneuron[i].isynsum;
			isyninput = 0;

			//neudex = i;

			if((int)modstep%1000 == 0) vmhneuron[i].raterec[modstep/1000] = vmhneuron[i].dend0e * 1000; 

			nepsp0 = 0;
			if(vmhneuron[i].dend0e > 0) {
				while(vmhneuron[i].epspt0 < hstep) {
					nepsp0++;
					vmhneuron[i].epspt0 = -log(1 - mrand01()) / vmhneuron[i].dend0e + vmhneuron[i].epspt0;
				}
				vmhneuron[i].epspt0 = vmhneuron[i].epspt0 - hstep;
			}

			nipsp0 = 0;
			if(vmhneuron[i].dend0i > 0) {
				while(vmhneuron[i].ipspt0 < hstep) {
					nipsp0++;
					vmhneuron[i].ipspt0 = -log(1 - mrand01()) / vmhneuron[i].dend0i + vmhneuron[i].ipspt0;
				}
				vmhneuron[i].ipspt0 = vmhneuron[i].ipspt0 - hstep;
			}

			nepsp1 = 0;
			nipsp1 = 0;
			
			if(signalmode) {
				if(nettime < 10 && i == 0) mod->diagbox->Write(text.Format("noise syn dend1e %.2f\n", vmhneuron[i].dend1e));
				if(vmhneuron[i].dend1e > 0) {
					while(vmhneuron[i].epspt1 < hstep) {
						nepsp1++;
						vmhneuron[i].epspt1 = -log(1 - mrand01()) / vmhneuron[i].dend1e + vmhneuron[i].epspt1;
					}
					vmhneuron[i].epspt1 = vmhneuron[i].epspt1 - hstep;
				}

				if(vmhneuron[i].dend1i > 0) {
					while(vmhneuron[i].ipspt1 < hstep) {
						nipsp1++;
						vmhneuron[i].ipspt1 = -log(1 - mrand01()) / vmhneuron[i].dend1i + vmhneuron[i].ipspt1;
					}
					vmhneuron[i].ipspt1 = vmhneuron[i].ipspt1 - hstep;
				}
			}

			//revpots = true;

			//fprintf(tofp, "psps generated\n");
			if(revpots) {
				epsph = vmhneuron[i].ae * (ve - vmhneuron[i].v);
				ipsph = vmhneuron[i].ai * (vmhneuron[i].v - vi);
				synh = vmhneuron[i].aesyn * (ve - vmhneuron[i].v);
			}
			else {
				epsph = pspmag;
				ipsph = pspmag;
				synh = synmag;
			}

			input0 = epsph * nepsp0 - ipsph * nipsp0; 
			input1 = epsph * nepsp1 - ipsph * nipsp1;
		
			sinput = synh * esyninput - synmag * isyninput;

			if(maxsyn && sinput > maxsyn) sinput = maxsyn;
			vmhneuron[i].synv = vmhneuron[i].synv - hstep * vmhneuron[i].synv / syntau + sinput;
			if(emax && vmhneuron[i].synv > emax) vmhneuron[i].synv = emax;

			vmhneuron[i].inputv = vmhneuron[i].inputv - hstep * vmhneuron[i].inputv / memtau + input0 + input1 + waveinput;
			vmhneuron[i].v = vmhneuron[i].vrest + vmhneuron[i].inputv + vmhneuron[i].synv;

			if(nettime < 100000 && i < vmhL1) vmhneuron[i].inputrec[modstep] = vmhneuron[i].inputv;

			vmhneuron[i].tHAP = vmhneuron[i].tHAP - hstep * vmhneuron[i].tHAP / vmhneuron[i].tauHAP;

			vmhneuron[i].tAHP = vmhneuron[i].tAHP - hstep * vmhneuron[i].tAHP / vmhneuron[i].tauAHP;  

			vmhneuron[i].tDAP = vmhneuron[i].tDAP - hstep * vmhneuron[i].tDAP / vmhneuron[i].tauDAP;  

			vmhneuron[i].th = vmhneuron[i].th0 + vmhneuron[i].tHAP + vmhneuron[i].tAHP - vmhneuron[i].tDAP;

			vmhneuron[i].tB = vmhneuron[i].tB - hstep * vmhneuron[i].tB / tauB;  

			if(vmndiag && nettime < 1000 && i == 0) {
				fprintf(tofp, "neuron %d  time %.1f  input0 %.2f  input1 %.2f  waveinput %.2f  sinput %.2f  v %.2f  thresh %.2f\n", 
					i, nettime, input0, input1, waveinput, sinput, vmhneuron[i].v, vmhneuron[i].th); 

				//fprintf(tofp, "dendinput %.2f  sinput %.2f\n", dendinput[i][1], sinput);
				fprintf(tofp, "esyninput %.2f  esynqueue %.2f isyninput %.2f\n", esyninput, esynqueue[i][0], isyninput);
				fprintf(tofp, "inputv %.2f  synv %.2f\n", vmhneuron[i].inputv, vmhneuron[i].synv);
				fprintf(tofp, "nepsp %d  epsph %.2f  nipsp %d  ipsph %.2f\n", nepsp0, epsph, nipsp0, ipsph);
				fprintf(tofp, "noisig %.2f  wavesig %.2f  synsig %.2f  dend1e %.2f  dend1i %.2f\n\n", noisig, wavesig, synsig, vmhneuron[i].dend1e, vmhneuron[i].dend1i);
			}

			if(vmndiag && nettime < 5000 && i == 0) {
				fprintf(inputofp, "%.4f %.4f %.4f %.4f\n", nettime, vmhneuron[i].inputv, input0, waveinput);
			}

			// Monitor L1 Neuron 0
			if(nettime < 500000 && i == 0) {
				mod->neurodata->V[nettime] = vmhneuron[i].v - vmhneuron[i].vrest;
				mod->neurodata->HAP[nettime] = vmhneuron[i].tHAP;
				mod->neurodata->AHP[nettime] = vmhneuron[i].tAHP;
				//mod->neurodata->netsyn[nettime] = vmhneuron[i].synv;
				mod->neurodata->netsyn[nettime] = synsig; // input1; //noisig;
				//mod->neurodata->netsyn[nettime] = waveinput;
				mod->neurodata->inpsyn[nettime] = input0;
				mod->neurodata->memsyn[nettime] = vmhneuron[i].inputv;
				mod->neurodata->facB[nettime] = vmhneuron[i].tB;
			}

			// Monitor L2 Neuron 0
			if(nettime < 500000 && i == vmhL1) {
				mod->neurodata->V2[nettime] = vmhneuron[i].v - vmhneuron[i].vrest;
				mod->neurodata->memsyn2[nettime] = vmhneuron[i].inputv + vmhneuron[i].synv;
				mod->neurodata->sigsyn2[nettime] = vmhneuron[i].inputv; // input1;
			}


			if(vmhneuron[i].v >= vmhneuron[i].th && vmhneuron[i].ttime >= absref) {          // Spike Fired

				if(s[i] < maxspikes) {
					vmhneuron[i].times[s[i]] = nettime;
					vmhneuron[i].mag[s[i]] = vmhneuron[i].tB;
					s[i]++;	
				}

				//if(i == 0) mod->diagbox->Write(text.Format("N%d spike %d  time %.2f  int %.2f\n", i, s[i], nettime/1000, vmhneuron[i].ttime));
				if(synqueue) activity[i] = 1;
				else activity[i] = syndelay * (double)1/hstep;

				vmhneuron[i].tHAP = vmhneuron[i].tHAP + vmhneuron[i].kHAP;
				vmhneuron[i].tAHP = vmhneuron[i].tAHP + vmhneuron[i].kAHP;
				vmhneuron[i].tDAP = vmhneuron[i].tDAP + vmhneuron[i].kDAP;
				if(vmhneuron[i].tDAP > maxdap) vmhneuron[i].tDAP = maxdap;

				vmhneuron[i].tB = vmhneuron[i].tB + kB;

				vmhneuron[i].ttime = 0;
			}
		}
	}

	for(i=0; i<vmhneurons; i++) {
		vmhneuron[i].spikecount = s[i];
		if(vmndiag) fprintf(tofp, "neuron %d fired %d spikes at %.2fHz, time %.2f\n", i, s[i], s[i]/(maxtime/1000), vmhneuron[i].neurotime);
	}	

	if(vmndiag) fclose(tofp);
	if(vmndiag) fclose(inputofp);
	//networkdisp();
	//inputoutput(0);

	for(i=0; i<maxcells; i++) delete[] weights[i];
	delete[] weights;

	for(i=0; i<maxcells; i++) delete[] delays[i];
	delete[] delays;

	for(i=0; i<maxcells; i++) delete[] esynqueue[i];
	delete[] esynqueue;

	for(i=0; i<maxcells; i++) delete[] enetwork[i];
	delete[] enetwork;
}

