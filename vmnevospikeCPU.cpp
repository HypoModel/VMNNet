
/*
*  evospikegen.cu
*  
*  Created by Tom Clayton and Duncan MacGregor.
*  University of Edinburgh 2017
*  Released under MIT license, see https://opensource.org/licenses/MIT
*
*/


// CPU version of CUDA spike model code, used to run GA with network model
//
// Currently not fully working


#include "vmnevofit.h"


EvoFitVMN_CPU::EvoFitVMN_CPU(int num, int index, float *pop, int pcount, int thcount, int block, int rtime, float *pInts, float *pISIs, float *pTemp, float *pSpikeCounts)
	: wxThread(wxTHREAD_JOINABLE)
{
	celldex = index;
	numcells = num;
	chromepop = pop;
	paramcount = pcount;
	threadcount = thcount;
	blocksize = block;
	runtime = rtime;
	Ints = pInts;
	ISIs = pISIs;
	Temp = pTemp;
	SpikeCounts = pSpikeCounts;
}


void *EvoFitVMN_CPU::Entry()
{
	wxString text;

	if(!celldex) diagbox->Write(text.Format("CPU thread %d running\n", celldex));

	SpikeGen();

	if(!celldex) diagbox->Write(text.Format("CPU thread %d SpikeGen OK, cell 0 %.0f spikes\n", celldex, SpikeCounts[0 + celldex*128]));

	return NULL;
}


void EvoFitVMN_CPU::SpikeGen()
{
	int i, j;
	wxString text;

	const unsigned int tid = celldex;            // fix this  9/11/17  DONE
	const unsigned int num_threads = 1;
	const unsigned int bid = 0;

	int d_iChromeCount = threadcount;
	float *d_fChrome = chromepop;
	int d_iRunTime = runtime;
	float d_fCutOffTop = 7052;
	float d_fInitBinWidth = 1;
	float d_fBinIncAmount = 0.05;
	float *d_oISI = ISIs;
	float *d_oInts = Ints;
	float *d_oTemp = Temp;
	float *d_oSpikeCounts = SpikeCounts;

	//FILE *ofp;
	//ofp = fopen("netcpu.txt", "w");
	//fprintf(ofp, "netcpu diagnostics\n\n");

	// Read neuron chrome
	float P_ESynRate				=	d_fChrome[0*d_iChromeCount + (tid + bid*num_threads)];
	float P_ISynRate				=	d_fChrome[1*d_iChromeCount + (tid + bid*num_threads)];
	float P_PSPMag					=	d_fChrome[2*d_iChromeCount + (tid + bid*num_threads)];
	float P_PSPDecay				=	d_fChrome[3*d_iChromeCount + (tid + bid*num_threads)];
	float P_RestToThreshold			=	d_fChrome[4*d_iChromeCount + (tid + bid*num_threads)];

	// Post Spike Potentials
	float P_kHAP					=	d_fChrome[5*d_iChromeCount + (tid + bid*num_threads)];
	float P_tauHAP					=	d_fChrome[6*d_iChromeCount + (tid + bid*num_threads)];
	float P_kDAP					=	d_fChrome[7*d_iChromeCount + (tid + bid*num_threads)];
	float P_tauDAP					=	d_fChrome[8*d_iChromeCount + (tid + bid*num_threads)];
	float P_kAHP					=	d_fChrome[9*d_iChromeCount + (tid + bid*num_threads)];
	float P_tauAHP					=	d_fChrome[10*d_iChromeCount + (tid + bid*num_threads)];

	// Network parameters
	float esynL1					=   d_fChrome[11*d_iChromeCount + (tid + bid*num_threads)];
	float esynweight				=   d_fChrome[12*d_iChromeCount + (tid + bid*num_threads)];
	float syndelay					=   d_fChrome[13*d_iChromeCount + (tid + bid*num_threads)];
	float syndelrange			    =   d_fChrome[14*d_iChromeCount + (tid + bid*num_threads)];

	int absref = 2;


	if(!celldex) diagbox->Write(text.Format("CPU SpikeGen thread %d runtime %d numcells %d esynweight %.2f\n", tid, d_iRunTime, numcells, esynweight));


	/////////////////////////////////////////////////////////////////////////////////////////////////////
	//Pre Calculate the stuff that is used to find the correct bin.  Most of this is just the parts of
	//the quadratic solver equation (-b + sqrt(b*b -4ac))/2a.
	float b = d_fInitBinWidth - d_fBinIncAmount/2;
	float BSquared = b*b;
	float FourA = 2*d_fBinIncAmount;
	float TwoA = d_fBinIncAmount;


	float nepsp, nipsp;
	float ExpectedNoE = P_ESynRate / 1000.0f;
	float ExpectedNoI = P_ISynRate / 1000.0f;

	int id = tid + bid * num_threads;
	//curandState randstate = state[id];



	// Tom Random Poisson Code
	//Create look-up table for poisson Distribution
	float InputE[8], InputI[8];
	float Factorial = 1;
	float CumulativeE = 0;
	float CumulativeI = 0;
	for(i=0; i<8; i++)
	{
		Factorial = 1;
		for(j=1; j<i+1; j++)
		{
			Factorial *= j;
		}

		CumulativeE += (pow(ExpectedNoE, i) * pow(2.718281828f, -ExpectedNoE)) / Factorial;
		InputE[(int)i] = CumulativeE;

		CumulativeI += (pow(ExpectedNoI, i) * pow(2.718281828f, -ExpectedNoI)) / Factorial;
		InputI[(int)i] = CumulativeI;
	}

	unsigned int RandExcit = 234248430 - tid + bid*num_threads;
	unsigned int RandInhib = 936753243 - tid + bid*num_threads;

	//numcells = 1;

	// Network Data
	//float esynL1 = 0.6;
	float esyntrans = 0.5;
	//float syndelay = 5;
	//float syndelrange = 0;
	//float esynweight = 1;

	unsigned char econnect[50];
	unsigned char enetwork[50][50];
	float esynqueue[50][20];
	unsigned char activity[50];
	
	float esynsum[50];
	int Spikes[50];


	/////////////////////////////////////////////////////////////////////////////////////////////////////
	//Create the Variables (HAP, DAP, AHP, LastDAP, Drift, Membrane Potential, Inhibition, Balance Point)


	//Input
	float V_Input[50];

	//Post Spike Potentials
	float V_HAP[50];
	float V_AHP[50];
	float V_DAP[50];

	// Initialise neuron variables
	for(i=0; i<numcells; i++) { 
		V_Input[i] = 0;
		V_HAP[i] = 0;
		V_DAP[i] = 0;
		V_AHP[i] = 0;
		activity[i] = 0;
		Spikes[i] = 0;
	}

	//Spikes
	//float Spike	= 0;
	float LastSpike	= 0;

	//Max average of 500 spikes per second
	int MaxSpikes = ((int)(d_iRunTime/1000)) * 500;

	//__syncthreads();


	//clear ISI memory locations
	for(i = 0; i < 512; i++) d_oISI[i + tid*512 + bid*num_threads*512] = 0.0f;

	//__syncthreads();


	float epspt = 0;
    float ipspt = 0;

	// Network Generation
	for(i=0; i<numcells; i++) {
		econnect[i] = 0;
		for(j=0; j<numcells; j++) {
			float d = mrand01(); // curand_uniform(&randstate);
			if(d <= esynL1 && i != j) enetwork[i][econnect[i]++] = j;
		}
		for(j=0; j<20; j++) esynqueue[i][j] = 0;       // queue max length fixed at 20
	}


	////////////////////////////////////////////////////////////
	// Model Loop
	for(float step=0; step<d_iRunTime; step++) {	
	
		// Network Input
		for(i=0; i<numcells; i++) {
			// Add network activity to input queue
			for(int c=0; c<econnect[i]; c++) 
				if(activity[enetwork[i][c]] == 1) {
					float synrand = mrand01(); //curand_uniform(&randstate);
					if(esyntrans >= synrand) { 
						float syndel = (syndelay - 1) + (syndelrange + 1) * (synrand * (1/esyntrans));
						esynqueue[i][(int)syndel] = esynqueue[i][(int)syndel] + esynweight;	
					}
				}

			//if(activity[i]) activity[i] = activity[i] - 1;

			// Read and shift input queue
			esynsum[i] = esynqueue[i][0];
			for(j=0; j<20-1; j++) esynqueue[i][j] = esynqueue[i][j+1];
			esynqueue[i][20-1] = 0;
		}
		

		// Membrane Activity and Spiking
		for(i=0; i<numcells; i++) {

			// Tom Poisson Code
			RandExcit = (22695477 * RandExcit + 1) & 0xFFFFFFFF;
			RandInhib = (22695477 * RandInhib + 1) & 0xFFFFFFFF;

			float fRandExcit = RandExcit / 4294967295.0f;
			float fRandInhib = RandInhib / 4294967295.0f;

			float Total = 0;
			for(int p=0; p<8; p++)
			{
				if (fRandExcit > InputE[p]) Total += P_PSPMag;
				if (fRandInhib > InputI[p]) Total -= P_PSPMag;
			}
		

			// CUDA Poisson Generator

			//nepsp = curand_poisson(&randstate, (double)ExpectedNoE);
			//nipsp = curand_poisson(&randstate, (double)ExpectedNoI);
		

			// CUDA Uniform Generator
			/*
			nepsp = 0;
			if(ExpectedNoE > 0) {
				while(epspt < 1) {
					nepsp++;
					epspt = -log(1 - curand_uniform(&randstate)) / ExpectedNoE + epspt;
				}
				epspt = epspt - 1;
			}

			nipsp = 0;
			if(ExpectedNoI > 0) {
				while(ipspt < 1) {
					nipsp++;
					ipspt = -log(1 - curand_uniform(&randstate)) / ExpectedNoI + ipspt;
				}
				ipspt = ipspt - 1;
			}
		

			float Total = P_PSPMag * nepsp - P_PSPMag * nipsp;*/

			float NetInput = P_PSPMag * esynsum[i];
			//float NetInput = 0;
		
			V_Input[i] = V_Input[i] - V_Input[i] * P_PSPDecay + Total + NetInput;


			//Calculate the Exponential Decays; Decay to the membrane rest

			V_HAP[i] = V_HAP[i] - V_HAP[i] * P_tauHAP + P_kHAP * activity[i];
			V_DAP[i] = V_DAP[i] - V_DAP[i] * P_tauDAP + P_kDAP * activity[i];
			V_AHP[i] = V_AHP[i] - V_AHP[i] * P_tauAHP + P_kAHP * activity[i];


			//Check to see if it has fired
			if((V_Input[i] - V_HAP[i] + V_DAP[i] - V_AHP[i]) > P_RestToThreshold && (step - LastSpike) >= absref) {
				//Spike = 1;
				activity[i] = 1;

				if(i == 0) {
					//If the model fires record and timestamp the event
					//float Diff = (step - LastSpike) * 0.001f;
					float Diff = step - LastSpike;

					if(Spikes[i] < 16384) d_oInts[Spikes[i] + tid*16384 + bid*num_threads*16384] = step;           
			
					if((Diff > 2) && (Diff < d_fCutOffTop)) {
						//This equation defines an ever increasing bin width
						//Even though it is quantised, it will increase by BinSpacing
						//every bin.
						float fBinNo = (-b + sqrt(BSquared + FourA * Diff)) / TwoA;
						int iBinNo = fBinNo;
						//Increment the relevent bin within Memory
						if((iBinNo < 512) & (iBinNo > 2)) d_oISI[iBinNo + tid*512 + bid*num_threads*512]++;
					}
					LastSpike = step;	
				}
				Spikes[i]++;

			}
			else {
				//Spike = 0;
				activity[i] = 0;
			}		
		}
		if(Spikes[0] > MaxSpikes) break;
	}

	////////////////////   End of model loop

	//__syncthreads();


	// Spike time analysis - generate quad binned ISI histogram

	//reset the temp
	for(i=0; i<512; i++) {
		d_oTemp[i + tid*512 + bid*num_threads*512] = 0.0f;
	} 

	//Smooth and Scale ISI Histogram
	float NumBinnedEvents = 0;
	for(i=0; i<512; i++) {
		float Div = 0.0f;

		for(j=(i-2); j<(i+3); j++) {	
			if((j > -1)&(j < 512)) {
				d_oTemp[i + tid*512 + bid*num_threads*512] += d_oISI[j + tid*512 + bid*num_threads*512];
				Div++;
			}
		}
		Div = d_oTemp[i + tid*512 + bid*num_threads*512] / Div;
		NumBinnedEvents += Div;
		d_oTemp[i + tid*512 + bid*num_threads*512] = Div;
	}

	//__syncthreads();

	if(NumBinnedEvents > 0.0f) {
		for(i = 0; i < 512; i++)
			d_oISI[i + tid*512 + bid*num_threads*512] = d_oTemp[i + tid*512 + bid*num_threads*512]/NumBinnedEvents;
	}
	else {
		for(i = 0; i < 512; i++) d_oISI[i + tid*512 + bid*num_threads*512] = 0.0f;
	}

	
	//__syncthreads();

	//If there are too many events... clean the Exp
	if(Spikes[0] > MaxSpikes) {
		for(i=0; i<512; i++)
			d_oISI[i + tid*512 + bid*num_threads*512] = 0.0f;
	}

	for(i=0; i<numcells; i++) d_oSpikeCounts[i + tid*128 + bid*num_threads*128] = Spikes[i];
	
	//__syncthreads();	
	
	//fprintf(ofp, "\nesyn %d  isyn %d\n", esyncount/1000, isyncount/1000);
	//fprintf(ofp, "\nSpikes %d\n", Spikes[0]);
	//fprintf(ofp, "\nSpikes %.0f\n", d_oMeanExtraBurstNoise[tid + bid*num_threads]);
	//for(i=0; i<10; i++) fprintf(ofp, "Spike %d Index %d Time %.2f\n", i, i + tid*16384 + bid*num_threads*16384, d_oInts[i + tid*16384 + bid*num_threads*16384]);
	//fclose(ofp);
}
