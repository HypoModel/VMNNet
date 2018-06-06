

/*
*  vmnevofit.cpp
*  
*  Created by Duncan MacGregor and Tom Clayton 
*  University of Edinburgh 2017
*  Released under MIT license, see https://opensource.org/licenses/MIT
*
*/

#include "vmnmodel.h"
#include "vmnevofit.h"
#include "hypodef.h"


EvoFitVMN::EvoFitVMN(VMNModel *model, EvoFitBox *fbox)
	: EvoFit((Model *)model, fbox)
{
	mod = model;
	fitbox = fbox;
	diagbox = mod->mainwin->diagbox;
}


void EvoFitVMN::Evaluate(int start, int popcount, double dualfit)
{
	int i, j;
	wxString text;


	ParamStore *netparams = mod->netbox->GetParams();

	int numcells = (*netparams)["neuronsL1"];

	ParamStore *params = fitchrome->parambox->GetParams();

	diagbox->Write("Net Evaluate, generate chrome array...\n");

	for(i=0; i<popcount; i++) {
		chrome = &(*chromepop)[i + start];
		chromearray[0 * popsize + i] = chrome->Val("ire");  
		chromearray[1 * popsize + i] = chrome->Val("ire") * chrome->Val("iratio");
		chromearray[2 * popsize + i] = chrome->Val("pspmag");
		chromearray[3 * popsize + i] = log((double)2) / chrome->Val("halflifeSyn");
		chromearray[4 * popsize + i] = chrome->Val("Vthresh") - (*params)["vrest"];
		chromearray[5 * popsize + i] = chrome->Val("kHAP");
		chromearray[6 * popsize + i] = log((double)2) / chrome->Val("halflifeHAP");
		chromearray[7 * popsize + i] = chrome->Val("kDAP");
		chromearray[8 * popsize + i] = log((double)2) / chrome->Val("halflifeDAP");
		chromearray[9 * popsize + i] = chrome->Val("kAHP");
		chromearray[10 * popsize + i] = log((double)2) / chrome->Val("halflifeAHP");
		chromearray[11 * popsize + i] = chrome->Val("esynL1");  
		chromearray[12 * popsize + i] = chrome->Val("esynweight");  
		chromearray[13 * popsize + i] = chrome->Val("syndelay");  
		chromearray[14 * popsize + i] = chrome->Val("syndelrange");  
	}

	diagbox->Write(text.Format("esynL1 %.2f\n", chrome->Val("esynL1")));
	diagbox->Write(text.Format("Vthresh %.4f\n", chrome->Val("Vthresh")));
	diagbox->Write(text.Format("Vrest %.4f\n", (*params)["vrest"]));
	diagbox->Write(text.Format("halflifeHAP %.4f\n", chrome->Val("halflifeHAP")));

	diagbox->Write("Evaluate, Chrome Array OK\n");
	
	//for(i=0; i<gpuparams; i++) diagbox->Write(text.Format("param %d : %.4f\n", i, chromearray[i * popsize]));

	if(!runmode) {
		diagbox->Write("Running CPU Fit " + text.Format("%.0f steps\n", runtime));
		float *Temp = new float[popcount * 512];
		for(i=0; i<popcount; i++) {
			chromethread[i] = new EvoFitVMN_CPU(numcells, i, chromearray, gpuparams, popcount, 1, runtime*1000, 
				fitbox->spikefitdata->Ints, fitbox->spikefitdata->ISIs, Temp, fitbox->spikefitdata->SpikeCounts);
			chromethread[i]->diagbox = diagbox;
			chromethread[i]->Create();
		}
		for(i=0; i<popcount; i++) chromethread[i]->Run();
		for(i=0; i<popcount; i++) chromethread[i]->Wait();
		for(i=0; i<popcount; i++) delete chromethread[i];	
		delete[] Temp;
		diagbox->Write(text.Format("CPU Fit OK, POP %d to %d\n", start, popcount + start));
	}
	else {
		diagbox->Write("Running GPU Fit " + text.Format("%.0f steps\n", runtime));
		EvoFitVMN_GPU(runmode, numcells, chromearray, gpuparams, popcount, blocksize, runtime*1000, 
			fitbox->spikefitdata->Ints, fitbox->spikefitdata->ISIs, fitbox->spikefitdata->SpikeCounts);
		diagbox->Write(text.Format("GPU Fit OK, POP %d to %d\n", start, popcount + start));
	}

	for(i=0; i<popcount; i++) {
		for(j=0; j<512; j++) {
			fitbox->evodata->histquadsm[j] = fitbox->spikefitdata->ISIs[j + i * 512];
		}
		if(netmode) fitbox->spikefitdata->spikecounts[i] = fitbox->spikefitdata->SpikeCounts[i * 128];
		else fitbox->spikefitdata->spikecounts[i] = fitbox->spikefitdata->SpikeCounts[i];

		// IoD non-parallel code
		//fitbox->evodata->spikecount = fitbox->spikefitdata->SpikeCounts[i * 128];	
	    fitbox->evodata->spikecount = fitbox->spikefitdata->spikecounts[i];	
		if(fitbox->evodata->spikecount > 16384) fitbox->evodata->spikecount = 16384;
		for(j=0; j<fitbox->evodata->spikecount; j++) fitbox->evodata->times[j] = fitbox->spikefitdata->Ints[j+32*512*i];

		fitbox->evodata->FitScoreOxy(fitbox->expdata, (*chromepop)[i + start].fitdata, fitbox->fitset, fitbox->fitconset);
		(*chromepop)[i + start].fitness = (*chromepop)[i + start].fitdata->score;
		(*chromepop)[i + start].fithead = (*chromepop)[i + start].fitdata->RMSFirstNBins;
		(*chromepop)[i + start].fittail = (*chromepop)[i + start].fitdata->RMSBinRange;
		(*chromepop)[i + start].fithaz = (*chromepop)[i + start].fitdata->RMSHaz;
		(*chromepop)[i + start].fitIoD = (*chromepop)[i + start].fitdata->RMSIoD;
		(*chromepop)[i + start].index = i;

		//if(i%32 == 0) diagbox->Write(text.Format("Chrome %d  Spikes %.0f  Score %.2f\n", i, fitbox->spikefitdata->SpikeCounts[i * 128], (*chromepop)[i + start].fitness)); 
		if(diagnostic && i%8 == 0) { 
			diagbox->Write(text.Format("Chrome %d  Score %.2f\n", i, (*chromepop)[i + start].fitness)); 
			diagbox->Write("Spikes ");
			for(j=0; j<numcells; j++) diagbox->Write(text.Format("%.0f  ", fitbox->spikefitdata->SpikeCounts[j + i * 128])); 
			diagbox->Write("\n");
		}
		else if(i%32 == 0) diagbox->Write(text.Format("Chrome %d  Spikes %.0f  Score %.2f\n", i, fitbox->spikefitdata->SpikeCounts[i], (*chromepop)[i + start].fitness)); 
	}

	diagbox->Write(text.Format("Evaluate OK\n"));
}


void EvoFitVMN::InitPop()
{
	int i, j;
	int size;
	wxString text, tag;
	
	fitchrome = fitbox->fitchrome;
	chromeparams = fitchrome->numparams;

	size = popsize + parentrange;
	if(fitbox->chromepopinit) chromepop->clear();
	fitbox->chromepopinit = true;

	fitchrome->Diagnostic();

	// Population chrome vector
	for(i=0; i<size; i++) {
		chromepop->push_back(EvoChrome());
		for(j=0; j<chromeparams; j++) {
			tag = fitchrome->tags[j];
			(*chromepop)[i].diagbox = mod->diagbox;
			//(*chromepop)[i].AddParam(tag, fitchrome->cons[j], &fitchrome->GetParam(tag)); 
			(*chromepop)[i].AddParam(tag, fitchrome->cons[j], &fitchrome->params[j]); 
		}
	}

	// Multirun result chrome vector
	if(multirun) {
		chromeresult->clear();	
		for(i=0; i<numruns; i++) {
			chromeresult->push_back(EvoChrome());
			for(j=0; j<chromeparams; j++) {
				tag = fitchrome->tags[j];
				(*chromeresult)[i].diagbox = mod->diagbox;
				(*chromeresult)[i].AddParam(tag, fitchrome->cons[j], &fitchrome->GetParam(tag)); 
			}
		}
	}
}


void *EvoFitVMN::Entry()
{
	int i, j;
	wxString text;
	int runindex;

	genmon = false;
	gpuparams = 15;
	runtime = 1000;

	// Read GA parameters

	ParamStore *evoparams = fitbox->GetParams();
	ParamStore *evoflags = fitbox->modflags;

	runtime = (float)(*evoparams)["evosteps"];
	popsize = (int)(*evoparams)["evopop"];
	parentrange = (int)(*evoparams)["parentrange"];
	blocksize = (int)(*evoparams)["blocksize"];
	generations = (int)(*evoparams)["evogens"];
	mutateprob = (*evoparams)["evomutprob"];
	numruns = (int)(*evoparams)["numruns"];

	multirun = (*fitbox->modflags)["multirun"];
	genlysis = (*fitbox->modflags)["genlysis"];
	
	netmode = (*mod->netbox->modflags)["GPUnetmode"];
	cpumode = (*mod->netbox->modflags)["cpumode"];
	diagnostic = (*mod->netbox->modflags)["GPUdiagnostic"];

	if(cpumode) runmode = 0;
	else if(netmode) runmode = 2;
	else runmode = 1;

	evoseed = (*evoparams)["evoseed"];
	seedgen = (*evoflags)["seedgen"];
	if(seedgen) {
		evoseed = (unsigned)(time(NULL));
		fitbox->paramset->GetCon("evoseed")->SetValue(evoseed);
	}
	init_mrand(evoseed);

	chromepop = &(fitbox->chromepop);
	chromeresult = &(fitbox->chromeresult);

	// Allocate or clear population storage
	InitPop();

	// Allocate chrome array for GPU transfer
	chromearray = new float[popsize * gpuparams];

	// Allocate fit data storage arrays for GPU output
	if(fitbox->spikefitdata->chromecount < popsize) {
		fitbox->spikefitdata->DeAllocate();
		fitbox->spikefitdata->Ints = new float[popsize * 512 * 32];
		fitbox->spikefitdata->ISIs = new float[popsize * 512];
		//fitbox->spikefitdata->IntraFreq = new float[popsize];
		//fitbox->spikefitdata->ExtraFreq = new float[popsize];
		fitbox->spikefitdata->SpikeCounts = new float[popsize * 128];
		fitbox->spikefitdata->chromecount = popsize;
	}

	threaddata = new int[popsize * 4];

	mod->mainwin->SetStatus("Fit Storage OK");

	if(!multirun) Evolve();          
	else {
		for(i=0; i<numruns; i++) {
			fitbox->runstatus->SetLabel(text.Format("%d", i+1));
			Evolve();
			(*chromeresult)[i] = (*chromepop)[0];
			(*chromeresult)[i].index = i;
			diagbox->Write(text.Format("\nRun %d Best chrome fitness %.2f\n\n", i, (*chromeresult)[i].fitness));
		}	
		TextFile resofp;
		wxString filetag = fitchrome->parambox->paramstoretag->GetValue();
		resofp.New(text.Format("chrome results %s.txt", filetag));
		quicksort(chromeresult, 0, numruns-1);
		for(i=0; i<numruns; i++) (*chromeresult)[i].Output(&resofp, 1); 
		resofp.Close();
	}

	// Analysis Output
	if(genlysis) {
		TextFile genofp;
		wxString outline;

		genofp.New(text.Format("genlysis results.txt"));
		for(i=1; i<=generations; i++) {
			outline.Printf("Gen\t%d\t", i);
			for(j=0; j<chromeparams; j++) 
				if((*chromepop)[0].params[j].adapt) 	
					outline += text.Format("P%d:\t%.6f\t", j, fitbox->paramsd[i][j] / fitbox->parammean[i][j]);	
			outline += text.Format("Fit\t%.6f\t%.6f\t", fitbox->fitmean[i], fitbox->fitsd[i]);
			genofp.WriteLine(outline);
		}
		genofp.Close();
	}

	diagbox->Write("GA Complete OK\n");

	fitbox->ChromeData();

	return NULL;
}


void EvoFitVMN::Evolve()
{
	wxString text;
	int i, j;

	int gen, pA, pB;
	bool orient;
	EvoChrome temp, parentA, parentB;
	EvoParam newparam;
	short crossA, crossB;
	double muteA, muteB, offset;

	temp = (*chromepop)[0];

	bool diagfile = false;
	bool gpudiag = false;


	// Generate initial population
	for(i=0; i<popsize; i++) {
		(*chromepop)[i].index = i;
		(*chromepop)[i].Initialise();	
	}

	diagbox->Write("Population Generation OK\n");

	Evaluate(0, popsize, dualfit);

	quicksort(chromepop, 0, popsize-1);

	diagbox->Write("Sort OK\n");

	if(gpudiag) {
		FILE *dofp;
		dofp = fopen("gpu-diagnostic.txt", "w");
		for(i=0; i<popsize; i++) {
			fprintf(dofp, "pop %d data index %d tid %d blocksize %d bid %d  syn %.4f\n\n", i, threaddata[0 + i*4], threaddata[1 + i*4], 
				threaddata[2 + i*4], threaddata[3 + i*4], fitbox->spikefitdata->SpikeCounts[i]);
		}
		fclose(dofp);
	}

	mod->mainwin->SetStatus("GA Running");

	// Diagnostic file

	if(diagfile) {
		ofp.New("gatest.txt");
		ofp.WriteLine(text.Format("Running GA, %d parameters", chromeparams));
		ofp.WriteLine("");
		ofp.WriteLine("");
	}

	diagbox->Write(text.Format("Running GA, %d parameters, %d parents, %d gens\n\n\n", chromeparams, parentrange, generations));


	for(gen=1; gen<=generations; gen++) {
		if(diagfile) {
			ofp.WriteLine("");
			ofp.WriteLine(text.Format("generation %d", gen));
		}
		diagbox->Write(text.Format("Generation %d...\n", gen));

		// Generate new generation

		for(i=0; i<popsize; i++) {
			pA = (int)(mrand01() * parentrange); 
			pB = (int)(mrand01() * parentrange); 
			parentA = (*chromepop)[pA]; 
			parentB = (*chromepop)[pB]; 
			crossA = ((int)(mrand01() * (chromeparams - 3)) + 1);
			crossB = ((int)(mrand01() * (chromeparams - 2)) + crossA);
			if(mrand01() > 0.5) orient = false; else orient = true;

			for(j=0; j<chromeparams; j++) {
				newparam = temp.params[j];
				if(temp.params[j].adapt) {
					if(mrand01() < mutateprob) {
						newparam.Generate();            // Mutate
						if(diagfile && genmon) ofp.WriteLine("mutate");
					}
					else {                           // Crossover
						if(diagfile && genmon) ofp.WriteLine("crossover");
						if(j < crossA || j > crossB) {
							if(orient) newparam = parentA.params[j];
							else newparam = parentB.params[j];
						}
						else {
							if(orient) newparam = parentB.params[j];
							else newparam = parentA.params[j];
						}
						muteA = mrand01();
						muteB = mrand01();
						offset = (muteA - muteB) * 0.5;
						offset = offset * (parentA.params[j].value - parentB.params[j].value);
						newparam.value += offset;
						if(newparam.value < newparam.min) newparam.value = newparam.min;
						if(newparam.value > newparam.max) newparam.value = newparam.max;
					}
				}
				(*chromepop)[i + parentrange].params[j] = newparam;
			}
			(*chromepop)[i + parentrange].index = popsize * gen + i;
		}

		diagbox->Write("Generation OK, Evaluate...\n");

		Evaluate(parentrange, popsize, dualfit);

		// Insert new generation
		quicksort(chromepop, 0, popsize + parentrange - 1);
	}

	if(diagfile) ofp.Close();
}
