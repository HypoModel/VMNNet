

/*
*  evospikegen.cu
*  
*  Created by Tom Clayton and Duncan MacGregor.
*  University of Edinburgh 2016
*  Released under MIT license, see https://opensource.org/licenses/MIT
*
*/

#include <cuda_runtime.h>
#include <cuda.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <curand.h>
#include <curand_kernel.h>
#include <mersenne.h>


#define BLOCK_DIM 16


__global__ void spikenetGPU(
	float *d_fChrome,
	int	d_iChromeCount,
	float d_fRunTime, 
	float d_fCutOffTop,
	float d_fInitBinWidth,
	float d_fBinIncAmount,
	float *d_oInts,
	float *d_oISI,
	float *d_oTemp,
	float *d_oInputVars,
	float *d_oMeanIntraBurstNoise,
	float *d_oMeanExtraBurstNoise,
	int *d_iThreadData) 
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// access thread id + No thread stuff.  This is used for addressing
	const unsigned int tid = threadIdx.x;
	const unsigned int num_threads = blockDim.x;
	const unsigned int bid = blockIdx.x;


	/////////////////////////////////////////////////////////////////////////////////////////////////////
	//Copy Control Parameters (CP) that are in perpetual use to local memory

	// Synaptic Input (IPSPs and EPSPs)
	
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


	d_iThreadData[0 + 4*tid + bid*num_threads*4] = (int)P_ESynRate;
	d_iThreadData[1 + 4*tid + bid*num_threads*4] = tid;
	d_iThreadData[2 + 4*tid + bid*num_threads*4] = num_threads;
	d_iThreadData[3 + 4*tid + bid*num_threads*4] = bid;

	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	//Pre Calculate the stuff that is used to find the correct bin.  Most of this is just the parts of
	//the quadratic solver equation (-b + sqrt(b*b -4ac))/2a.
	float b = d_fInitBinWidth - d_fBinIncAmount/2;
	float BSquared = b*b;
	float FourA = 2*d_fBinIncAmount;
	float TwoA = d_fBinIncAmount;

	//Create look-up table for poisson Distribution
	float InputE[8], InputI[8];
	float ExpectedNoE = P_ESynRate / 1000.0f;
	float ExpectedNoI = P_ISynRate / 1000.0f;
	float Factorial = 1;
	float CumulativeE = 0;
	float CumulativeI = 0;
	for(float i=0; i<8; i++)
	{
		Factorial = 1;
		for(int j=1; j<i+1; j++)
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


	/////////////////////////////////////////////////////////////////////////////////////////////////////
	//Create the Variables (HAP, DAP, AHP, LastDAP, Drift, Membrane Potential, Inhibition, Balance Point)

	//Input
	float V_Input = 0;

	//Post Spike Potentials
	float V_HAP	= 0;
	float V_DAP = 0;
	float V_AHP = 0;

	//Spikes
	float Spike	= 0;
	float LastSpike	= 0;

	//Max average of 500 spikes per second
	int MaxSpikes = ((int)(d_fRunTime/1000)) * 500;

	__syncthreads();

	int Spikes = 0;

	//clear ISI memory locations
	for(int i = 0; i < 512; i++)
		d_oISI[i + tid*512 + bid*num_threads*512] = 0.0f;

	//Input
	d_oInputVars[0 + tid*4 + bid*num_threads*4] = 0.0f; 
	d_oInputVars[1 + tid*4 + bid*num_threads*4] = 0.0f; 
	d_oInputVars[2 + tid*4 + bid*num_threads*4] = 0.0f; 
	d_oInputVars[3 + tid*4 + bid*num_threads*4] = 0.0f; 

	__syncthreads();

	////////////////////////////////////////////////////////////
	//Model Starts Here
	for(float step=0; step < d_fRunTime; step++)
	{	
		RandExcit = (22695477 * RandExcit + 1) & 0xFFFFFFFF;
		RandInhib = (22695477 * RandInhib + 1) & 0xFFFFFFFF;

		float fRandExcit = RandExcit / 4294967295.0f;
		float fRandInhib = RandInhib / 4294967295.0f;

		float Total = 0;
		for(int i=0; i<8; i++)
		{
			if (fRandExcit > InputE[i]) Total += P_PSPMag;
			if (fRandInhib > InputI[i]) Total -= P_PSPMag;
		}

		V_Input = V_Input - V_Input * P_PSPDecay + Total;


		//Calculate the Exponential Decays; Decay to the membrane rest

		V_HAP = V_HAP - V_HAP * P_tauHAP + P_kHAP * Spike;
		V_DAP = V_DAP - V_DAP * P_tauDAP + P_kDAP * Spike;
		V_AHP = V_AHP - V_AHP * P_tauAHP + P_kAHP * Spike;


		//Check to see if it has fired
		if((V_Input - V_HAP + V_DAP - V_AHP) > P_RestToThreshold) {
			Spike = 1;

			//If the model fires record and timestamp the event
			//float Diff = (step - LastSpike) * 0.001f;
			float Diff = step - LastSpike;

			if(Spikes < 16384) d_oInts[Spikes + tid*16384 + bid*num_threads*16384] = step;           
			
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
			Spikes++;
		}
		else Spike = 0;

		if(Spikes > MaxSpikes) break;
	}

	////////////////////   End of model loop

	__syncthreads();


	// Spike time analysis - generate quad binned ISI histogram

	//reset the temp
	for(int i=0; i<512; i++) {
		d_oTemp[i + tid*512 + bid*num_threads*512] = 0.0f;
	} 

	//Smooth and Scale ISI Histogram
	float NumBinnedEvents = 0;
	for(int i=0; i<512; i++) {
		float Div = 0.0f;

		for (int j = (i-2); j < (i+3); j++) {	
			if((j > -1)&(j < 512)) {
				d_oTemp[i + tid*512 + bid*num_threads*512] += d_oISI[j + tid*512 + bid*num_threads*512];
				Div++;
			}
		}
		Div = d_oTemp[i + tid*512 + bid*num_threads*512] / Div;
		NumBinnedEvents += Div;
		d_oTemp[i + tid*512 + bid*num_threads*512] = Div;
	}

	__syncthreads();

	if(NumBinnedEvents > 0.0f) {
		for(int i = 0; i < 512; i++)
			d_oISI[i + tid*512 + bid*num_threads*512] = d_oTemp[i + tid*512 + bid*num_threads*512]/NumBinnedEvents;
	}
	else {
		for(int i = 0; i < 512; i++) d_oISI[i + tid*512 + bid*num_threads*512] = 0.0f;
	}

	//Save the Noise Value;
	d_oMeanIntraBurstNoise[tid + bid*num_threads] = d_oInputVars[0 + tid*4 + bid*num_threads*4]/d_oInputVars[1 + tid*4 + bid*num_threads*4];
	d_oMeanExtraBurstNoise[tid + bid*num_threads] = d_oInputVars[2 + tid*4 + bid*num_threads*4]/d_oInputVars[3 + tid*4 + bid*num_threads*4];

	__syncthreads();

	//If there are too many events... clean the Exp
	if(Spikes > MaxSpikes) {
		for(int i=0; i<512; i++)
			d_oISI[i + tid*512 + bid*num_threads*512] = 0.0f;
	}


	d_oMeanExtraBurstNoise[tid + bid*num_threads] = Spikes;
	//d_oMeanExtraBurstNoise[tid + bid*num_threads] = P_SynRate;
	__syncthreads();	
}


__global__ void spikevmnnetGPU(
	int numcells,
	float *d_fChrome,
	int d_iChromeCount,
	int d_iRunTime,
	float d_fCutOffTop,
	float d_fInitBinWidth,
	float d_fBinIncAmount,
	float *d_oInts,
	float *d_oISI,
	float *d_oTemp,
	float *d_oSpikeCounts,
	curandState *state)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int num_threads = blockDim.x;
	const unsigned int bid = blockIdx.x;

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

	int id = threadIdx.x + blockIdx.x * blockDim.x;
	curandState randstate = state[id];



	// Tom Random Poisson Code
	//Create look-up table for poisson Distribution
	float InputE[8], InputI[8];
	float Factorial = 1;
	float CumulativeE = 0;
	float CumulativeI = 0;
	for(float i=0; i<8; i++)
	{
		Factorial = 1;
		for(int j=1; j<i+1; j++)
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

	int i;
	//numcells = 1;

	// Network Data
	//float esynL1 = 0.6;
	float esyntrans = 0.5;
	//float syndelay = 5;
	//float syndelrange = 0;
	//float esynweight = 1;

	unsigned char econnect[50];
	unsigned char enetwork[50][50];
	unsigned char esynqueue[50][20];
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

	__syncthreads();


	//clear ISI memory locations
	for(i = 0; i < 512; i++) d_oISI[i + tid*512 + bid*num_threads*512] = 0.0f;

	__syncthreads();


	float epspt = 0;
    float ipspt = 0;

	// Network Generation
	for(i=0; i<numcells; i++) {
		econnect[i] = 0;
		for(float j=0; j<numcells; j++) {
			float d = curand_uniform(&randstate);
			if(d <= esynL1 && i != j) enetwork[i][econnect[i]++] = j;
		}
		for(int j=0; j<20; j++) esynqueue[i][j] = 0;       // queue max length fixed at 20
	}


	////////////////////////////////////////////////////////////
	// Model Loop
	for(float step=0; step < d_iRunTime; step++) {	
	
		// Network Input
		for(i=0; i<numcells; i++) {
			// Add network activity to input queue
			for(int c=0; c<econnect[i]; c++) 
				if(activity[enetwork[i][c]] == 1) {
					float synrand = curand_uniform(&randstate);
					if(esyntrans >= synrand) { 
						float syndel = (syndelay - 1) + (syndelrange + 1) * (synrand * (1/esyntrans));
						esynqueue[i][(int)syndel] = esynqueue[i][(int)syndel] + esynweight;	
					}
				}

			//if(activity[i]) activity[i] = activity[i] - 1;

			// Read and shift input queue
			esynsum[i] = esynqueue[i][0];
			for(int j=0; j<20-1; j++) esynqueue[i][j] = esynqueue[i][j+1];
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
			if((V_Input[i] - V_HAP[i] + V_DAP[i] - V_AHP[i]) > P_RestToThreshold) {
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

	__syncthreads();


	// Spike time analysis - generate quad binned ISI histogram

	//reset the temp
	for(int i=0; i<512; i++) {
		d_oTemp[i + tid*512 + bid*num_threads*512] = 0.0f;
	} 

	//Smooth and Scale ISI Histogram
	float NumBinnedEvents = 0;
	for(int i=0; i<512; i++) {
		float Div = 0.0f;

		for (int j = (i-2); j < (i+3); j++) {	
			if((j > -1)&(j < 512)) {
				d_oTemp[i + tid*512 + bid*num_threads*512] += d_oISI[j + tid*512 + bid*num_threads*512];
				Div++;
			}
		}
		Div = d_oTemp[i + tid*512 + bid*num_threads*512] / Div;
		NumBinnedEvents += Div;
		d_oTemp[i + tid*512 + bid*num_threads*512] = Div;
	}

	__syncthreads();

	if(NumBinnedEvents > 0.0f) {
		for(int i = 0; i < 512; i++)
			d_oISI[i + tid*512 + bid*num_threads*512] = d_oTemp[i + tid*512 + bid*num_threads*512]/NumBinnedEvents;
	}
	else {
		for(int i = 0; i < 512; i++) d_oISI[i + tid*512 + bid*num_threads*512] = 0.0f;
	}

	
	__syncthreads();

	//If there are too many events... clean the Exp
	if(Spikes[0] > MaxSpikes) {
		for(int i=0; i<512; i++)
			d_oISI[i + tid*512 + bid*num_threads*512] = 0.0f;
	}


	for(i=0; i<numcells; i++) d_oSpikeCounts[i + tid*128 + bid*num_threads*128] = Spikes[i];
	
	__syncthreads();	



	/*
	// curand test code

	curandState_t state;

	curand_init(0, 0, 0, &state);

	for(int i=0; i<10; i++) d_oInts[i + tid*16384 + bid*num_threads*16384] = curand(&state)*(1.0/UINT_MAX);
	d_oSpikeCounts[tid + bid*num_threads] = 10;

	__syncthreads();	
	*/
}


__global__ void spikegenvmnGPU(
	float *d_fChrome,
	int d_iChromeCount,
	int d_iRunTime,
	float d_fCutOffTop,
	float d_fInitBinWidth,
	float d_fBinIncAmount,
	float *d_oInts,
	float *d_oISI,
	float *d_oTemp,
	float *d_oSpikeCounts,
	curandState *state)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int num_threads = blockDim.x;
	const unsigned int bid = blockIdx.x;

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

	int id = threadIdx.x + blockIdx.x * blockDim.x;
	curandState randstate = state[id];


	// Tom Random Poisson Code
	//Create look-up table for poisson Distribution
	float InputE[8], InputI[8];
	float Factorial = 1;
	float CumulativeE = 0;
	float CumulativeI = 0;
	for(float i=0; i<8; i++)
	{
		Factorial = 1;
		for(int j=1; j<i+1; j++)
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

	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	//Create the Variables (HAP, DAP, AHP, LastDAP, Drift, Membrane Potential, Inhibition, Balance Point)


	//Input
	float V_Input = 0;

	//Post Spike Potentials
	float V_HAP = 0;
	float V_AHP = 0;
	float V_DAP = 0;

	//Spikes
	float Spike	= 0;
	float LastSpike	= 0;

	//Max average of 500 spikes per second
	int MaxSpikes = ((int)(d_iRunTime/1000)) * 500;

	__syncthreads();

	int Spikes = 0;

	//clear ISI memory locations
	for(int i=0; i<512; i++) d_oISI[i + tid*512 + bid*num_threads*512] = 0.0f;

	__syncthreads();


	////////////////////////////////////////////////////////////
	//Model Starts Here
	for(float step=0; step < d_iRunTime; step++) {

		// Membrane Activity and Spiking
		
		// Tom Poisson Code
		RandExcit = (22695477 * RandExcit + 1) & 0xFFFFFFFF;
		RandInhib = (22695477 * RandInhib + 1) & 0xFFFFFFFF;

		float fRandExcit = RandExcit / 4294967295.0f;
		float fRandInhib = RandInhib / 4294967295.0f;

		float Total = 0;
		for(int p=0; p<8; p++) {
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

		//float NetInput = P_PSPMag * esynsum[i];
		float NetInput = 0;
		
		V_Input = V_Input - V_Input * P_PSPDecay + Total;


		//Calculate the Exponential Decays; Decay to the membrane rest

		V_HAP = V_HAP - V_HAP * P_tauHAP + P_kHAP * Spike;
		V_DAP = V_DAP - V_DAP * P_tauDAP + P_kDAP * Spike;
		V_AHP = V_AHP - V_AHP * P_tauAHP + P_kAHP * Spike;


		//Check to see if it has fired
		if((V_Input - V_HAP + V_DAP - V_AHP) > P_RestToThreshold) {
			Spike = 1;
			
			//If the model fires record and timestamp the event
			//float Diff = (step - LastSpike) * 0.001f;
			float Diff = step - LastSpike;

			if(Spikes < 16384) d_oInts[Spikes + tid*16384 + bid*num_threads*16384] = step;           
			
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
			Spikes++;
		}
		else Spike = 0;

	    if(Spikes > MaxSpikes) break;
	}

	////////////////////   End of model loop

	__syncthreads();


	// Spike time analysis - generate quad binned ISI histogram

	//reset the temp
	for(int i=0; i<512; i++) {
		d_oTemp[i + tid*512 + bid*num_threads*512] = 0.0f;
	} 

	//Smooth and Scale ISI Histogram
	float NumBinnedEvents = 0;
	for(int i=0; i<512; i++) {
		float Div = 0.0f;

		for (int j = (i-2); j < (i+3); j++) {	
			if((j > -1)&(j < 512)) {
				d_oTemp[i + tid*512 + bid*num_threads*512] += d_oISI[j + tid*512 + bid*num_threads*512];
				Div++;
			}
		}
		Div = d_oTemp[i + tid*512 + bid*num_threads*512] / Div;
		NumBinnedEvents += Div;
		d_oTemp[i + tid*512 + bid*num_threads*512] = Div;
	}

	__syncthreads();

	if(NumBinnedEvents > 0.0f) {
		for(int i = 0; i < 512; i++)
			d_oISI[i + tid*512 + bid*num_threads*512] = d_oTemp[i + tid*512 + bid*num_threads*512]/NumBinnedEvents;
	}
	else {
		for(int i = 0; i < 512; i++) d_oISI[i + tid*512 + bid*num_threads*512] = 0.0f;
	}

	
	__syncthreads();

	//If there are too many events... clean the Exp
	if(Spikes > MaxSpikes) {
		for(int i=0; i<512; i++)
			d_oISI[i + tid*512 + bid*num_threads*512] = 0.0f;
	}


	d_oSpikeCounts[tid + bid*num_threads] = Spikes;
	
	__syncthreads();	



	/*
	// curand test code

	curandState_t state;

	curand_init(0, 0, 0, &state);

	for(int i=0; i<10; i++) d_oInts[i + tid*16384 + bid*num_threads*16384] = curand(&state)*(1.0/UINT_MAX);
	d_oSpikeCounts[tid + bid*num_threads] = 10;

	__syncthreads();	
	*/
}


__global__ void randinitGPU(unsigned long seed, curandState *state)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
    /* Each thread gets same seed, a different sequence 
       number, no offset */
    curand_init(seed, id, 0, &state[id]);
}


void spikevmnnetCPU(
	int numcells,
	float *d_fChrome,
	int d_iChromeCount,
	int d_iRunTime,
	float d_fCutOffTop,
	float d_fInitBinWidth,
	float d_fBinIncAmount,
	float *d_oInts,
	float *d_oISI,
	float *d_oTemp,
	float *d_oSpikeCounts)
{
	//const unsigned int tid = threadIdx.x;
	//const unsigned int num_threads = blockDim.x;
	//const unsigned int bid = blockIdx.x;

	const unsigned int tid = 0;
	const unsigned int num_threads = 1;
	const unsigned int bid = 0;

	int i, j;

	FILE *ofp;
	ofp = fopen("netcpu.txt", "w");
	fprintf(ofp, "netcpu diagnostics\n\n");

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
	unsigned char esynqueue[50][20];
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
	for(float step=0; step < d_iRunTime; step++) {	
	
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
			if((V_Input[i] - V_HAP[i] + V_DAP[i] - V_AHP[i]) > P_RestToThreshold) {
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
	fprintf(ofp, "\nSpikes %d\n", Spikes[0]);
	//fprintf(ofp, "\nSpikes %.0f\n", d_oMeanExtraBurstNoise[tid + bid*num_threads]);
	for(i=0; i<10; i++) fprintf(ofp, "Spike %d Index %d Time %.2f\n", i, i + tid*16384 + bid*num_threads*16384, d_oInts[i + tid*16384 + bid*num_threads*16384]);
	fclose(ofp);
}


void EvoFitVMN_GPU(int runmode, int numcells, float *chromepop, int paramcount, int threadcount, int blocksize, float runtime, float *Ints, float *ISIs, float *SpikeCounts)
{
	float *d_chromepop;
	float *d_Ints, *d_SpikeCounts;
	float *d_ISIs, *d_Temp;
	unsigned long spikeseed;

	curandState *state;

    cudaSetDevice(0);

	cudaMallocManaged((void **) &d_chromepop, threadcount * paramcount * sizeof(float));
	cudaMallocManaged((void **) &d_Ints, threadcount * 512 * 32 * sizeof(float));
	cudaMallocManaged((void **) &d_ISIs, threadcount * 512 * sizeof(float));
	cudaMallocManaged((void **) &d_Temp, threadcount * 512 * sizeof(float));
	cudaMallocManaged((void **) &d_SpikeCounts, threadcount * 128 * sizeof(float));
	cudaMemcpy(d_chromepop, chromepop, threadcount * paramcount * sizeof(float), cudaMemcpyHostToDevice);
	
	int blockSize = blocksize;
	int blocks = threadcount / blockSize;

	dim3 grid(1, 1, 1);
	dim3 threads(threadcount, 1, 1);

	float *Temp = new float[threadcount * 512];

	FILE *ofp;
	ofp = fopen("netgpu.txt", "w");
	fprintf(ofp, "netgpu diagnostics\n\n");

	// Allocate space for Poisson RNG states on device
    cudaMallocManaged((void **) &state, threadcount * sizeof(curandState));
	// Initialise PRNGs
	spikeseed = (unsigned long)time(NULL);
	randinitGPU<<< blocks, blocksize>>>(spikeseed, state);

	if(runmode == 0) {
		fprintf(ofp, "running CPU net mode\n");
		spikevmnnetCPU(numcells, chromepop, threadcount, runtime, 7052, 1, 0.05, Ints, ISIs, Temp, SpikeCounts);
	}
	else if(runmode == 1) spikegenvmnGPU<<< blocks, blockSize >>>(d_chromepop, threadcount, runtime, 7052, 1, 0.05, d_Ints, d_ISIs, d_Temp, d_SpikeCounts, state);
	else {
		fprintf(ofp, "running gpu net mode\n");
		spikevmnnetGPU<<< blocks, blockSize >>>(numcells, d_chromepop, threadcount, runtime, 7052, 1, 0.05, d_Ints, d_ISIs, d_Temp, d_SpikeCounts, state);
	}
	
	if(runmode) {
		cudaMemcpy(Ints, d_Ints, threadcount * 512 * 32 * sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(ISIs, d_ISIs, threadcount * 512 * sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(SpikeCounts, d_SpikeCounts, threadcount * 128 * sizeof(float), cudaMemcpyDeviceToHost);
	}
	
	cudaFree(d_chromepop);
	cudaFree(d_Ints);
	cudaFree(d_ISIs);
	cudaFree(d_Temp);
	cudaFree(d_SpikeCounts);

	delete[] Temp;

	fclose(ofp);
}