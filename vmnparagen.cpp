


#include "vmnmodel.h"


DEFINE_EVENT_TYPE(wxEVT_SPIKE)
DEFINE_EVENT_TYPE(wxEVT_SYNCH)


VMHspikegen::VMHspikegen(int nstart, int nstop, datint *act, VMNMod *p_vmhmod, int id)
: wxThread(wxTHREAD_JOINABLE)
{
	threadID = id;
	nfirst = nstart;
	nlast = nstop;
	vmhmod = p_vmhmod;
	vmhneuron = vmhmod->vmhneuron;
	actmutex = vmhmod->actmutex;
	activity = act;
	vmhbox = vmhmod->netbox;
	synapse = vmhmod->synapse;
	
	numspikes = vmhmod->numspikes;
	hstep = vmhmod->hstep;
	vthre = vmhmod->vthre[0];
	vrest = vmhmod->vrest[0];
	//ae = vmhmod->ae;
	//ai = vmhmod->ai;
	ve = vmhmod->ve;
	vi = vmhmod->vi;
	iratio = vmhmod->iratio;
	ire = vmhmod->ire;
	halflife = vmhmod->halflife;
	kHAP = vmhmod->kHAP[0];
	tauHAP = vmhmod->tauHAP[0];
	kAHP = vmhmod->kAHP[0];
	tauAHP = vmhmod->tauAHP[0];
	kDAP = vmhmod->kDAP[0];
	tauDAP = vmhmod->tauDAP[0];
	gamma = log((double)2) / halflife;
	
	// VMH Network parameters
	
	vmhneurons = vmhmod->vmhneurons;
	vmhconnect = vmhmod->vmhconnect;
	vmhinput1 = vmhmod->vmhinput1;
	synweight = vmhmod->synweight;
	absref = vmhmod->absref;
	inputcycle = vmhmod->inputcycle;
	waveamp = vmhmod->waveamp;
	syndelay = vmhmod->syndelay;
	esyn = vmhmod->esyn;
	isyn = vmhmod->isyn;
	vrestsd = vmhmod->vrestsd[0];
	kHAPsd = vmhmod->kHAPsd[0];
	tauHAPsd = vmhmod->tauHAPsd[0];
	maxsyn = vmhmod->maxsyn;
}


void *VMHspikegen::Entry()
{
	//activity[0] = 0;
	spikegen();
	return NULL;
}


void VMHspikegen::spikegen()
{		
	// October 2009
	// March 2010
	
	int i, j, tmax = 0, n, tstep, c;
	int nepsp0, nipsp0, nepsp1, nipsp1;
	int numconnect, con;
	int modstep;
	int modeltime;
	int active[1000];
	int esynsum[1000];
	int postok;
	int synch = 0;
	
	int s[1000];
	double maxtime, nettime;
	double maxdap = 50;
	
	double ire1 = ire;
	double iratio1 = iratio;
	double dend0e, dend0i, dend1e, dend1i;
	double connectsum;
	double dendsum0, dendsum1, endosum;
	double esyninput, isyninput, sinput;
	double suckstim;
	double waveinput;
	double pi = 3.14159265;
	
	wxCommandEvent synchevt(wxEVT_SYNCH, wxID_ANY);
	wxCommandEvent spikeevt(wxEVT_SPIKE, wxID_ANY);
	synchevt.SetInt(threadID);
	
	FILE *ofp;
	FILE *tofp[10];
	if(threadID == 0) tofp[threadID] = fopen("vmhtest0.txt", "w");
	if(threadID == 1) tofp[threadID] = fopen("vmhtest1.txt", "w");
	
	/// Param copy for testing
	
	th0 = vthre;
	//v_rest = vrest;
	
	maxtime = numspikes * 1000;
	
	//for(i=0; i<=Nspikes; i++)
	//	netisi[i] = 0;
	
	memtau = 1 / (log((double)2) / halflife);
	
	dend0e = vmhinput1 / 1000;
	dend0i = dend0e * iratio;
	
	dend1e = vmhinput1 / 1000;
	dend1i = dend1e * iratio1;
	
	
	// initial values
	
	for(i=nfirst; i<=nlast; i++) {
		vmhneuron[i].vrest = gaussian(vrest, vrestsd);
		//vmhneuron[i].vrest = uniform(vrest, vrestsd);
		vmhneuron[i].th0 = gaussian(vthre, 0);
		vmhneuron[i].kHAP = gaussian(kHAP, kHAPsd);
		vmhneuron[i].kAHP = gaussian(kAHP, 0);
		vmhneuron[i].kDAP = gaussian(kDAP, 0);
		vmhneuron[i].tauHAP = gaussian(tauHAP, tauHAPsd);
		vmhneuron[i].tauAHP = gaussian(tauAHP, 0);
		vmhneuron[i].tauDAP = gaussian(tauDAP, 0);
		
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
		//actmutex->Lock();
		(*activity)[i] = 0;
		//actmutex->Unlock();
		active[i] = 0;
		//vmhneuron[i].inputrec[0] = vrest;
	}
	
	modeltime = 0;
	nettime = 0;
	hstep = 1;
	modstep = 0;
	int gocheck;
	
	//vmhmod->synchmutex[threadID]->Lock();
	
	
	
	fprintf(tofp[threadID], "Start spike gen, run for %.0fms  Thread %d\n", maxtime, threadID);
	fflush(tofp[threadID]);
	
	////////////////////////
	
	while(nettime<maxtime) {
		
		//vmhmod->AddPendingEvent(synchevt);
		//synapse->Receive(postok);
		//vmhmod->synchmutex[threadID]->Unlock();
		//fprintf(tofp[threadID], "%.2f ms Thread %d Call synch\n", nettime, threadID);
		//fflush(tofp[threadID]);
		//vmhmod->goflag[threadID] = 0;
		//vmhmod->OnSynch(synchevt);                                                      // change for ModThread
		//fprintf(tofp[threadID], "%.2f ms Thread %d Wait\n", nettime, threadID);
		//fflush(tofp[threadID]);
		//vmhmod->synchmutex[threadID]->Lock();
		//gocheck = 0;
		//while(gocheck == 0) {
		//	vmhmod->synchmutex[threadID]->Lock();
		//	gocheck = vmhmod->goflag[threadID];
		//	vmhmod->synchmutex[threadID]->Unlock();
		//}
		//fprintf(tofp[threadID], "%.2f ms Thread %d Check\n", nettime, threadID);
		//fflush(tofp[threadID]);
		//Sleep(1);
		//}
		//fprintf(tofp[threadID], "%.2f ms Thread %d Go\n", nettime, threadID);
		//fflush(tofp[threadID]);
		//fprintf(tofp, "Try Lock  ");
		//fflush(tofp);
		//vmhmod->synchsemaphore->Post();
		//vmhmod->synchsemaphore->Wait();
		
		//vmhmod->synchmutex[0]->Lock();
		//vmhmod->synchcond[0]->Signal();
		//vmhmod->synchmutex[0]->Unlock();
		vmhmod->synchcond[threadID]->Wait();
		
		//fprintf(tofp, "Signal OK\n");
		//fflush(tofp);
		//synch = vmhmod->synchOK;
		//synch = 1;
		
		//actmutex->Lock();
		//vmhmod->synchcond[threadID]->Wait();
		//while(synch) synch = vmhmod->synchflag[threadID]; 
		//synapse->Receive(postok);
		
		modstep++;
		nettime = nettime + hstep;
		if((int)nettime%1000 == 0) vmhbox->SetCount(floor(nettime)/maxtime*100);
		
		
		//fprintf(tofp, "time %.2f\n", nettime);
		
		waveinput = 0;
		if(inputcycle) waveinput = waveamp * 0.5 * (1 + sin(nettime*2*pi/inputcycle));
		
		/*
		 dend0e = vmhinput1 / 1000;
		 dend0i = dend0e * iratio;
		 
		 dend1e = vmhinput1 / 1000;
		 dend1i = dend1e * iratio1;
		 
		 for(i=0; i<vmhneurons; i++) {
		 dendsum1 = 0;
		 for(con=0; con<vmhconnect; con++) 
		 if(activity[network[i][con]] == 1) dendsum1++;	
		 dendinput[i][1] = synweight * dendsum1;
		 }
		 */
		
		//fprintf(tofp, "Inputs: neuron %d vmhconnect = %.2f synweight = %.2f dendsum1 = %.2f\n", i, vmhconnect, synweight, dendsum1);
		
		for(i=nfirst; i<=nlast; i++) esynsum[i] = 0;
		if(esyn > 0) {
			actmutex->Lock(); 
			for(i=nfirst; i<=nlast; i++) {
				if((*activity)[i]) (*activity)[i] = (*activity)[i] - 1;
				esynsum[i] = 0;
				for(c=0; c<vmhneuron[i].econnect; c++) 
					if((*activity)[vmhneuron[i].enetwork[c]] == 1) esynsum[i]++;
			}
			actmutex->Unlock(); 
		}
		
		
		for(i=nfirst; i<=nlast; i++) {
			
			//fprintf(tofp, "neuron %d\n", i);
			vmhneuron[i].ttime = vmhneuron[i].ttime + hstep;
			vmhneuron[i].neurotime = vmhneuron[i].neurotime + hstep;
			
			/*
			 if(active[i]) {
			 active[i]--;
			 if(active[i] == 0) {
			 //wxCommandEvent evt(wxEVT_SPIKE, wxID_ANY);
			 spikeevt.SetInt(i); 
			 vmhmod->AddPendingEvent(spikeevt);
			 }
			 }
			 //wxCommandEvent evt(wxEVT_SYNCH, wxID_ANY);
			 //evt.SetInt(i); 
			 //vmhmod->AddPendingEvent(synchevt);
			 */
			/*
			 if(active) {
			 active--;
			 actmutex->Lock(); 
			 if(active == 1) for(c=0; c<vmhneuron[i].econnect; c++) (*activity)[vmhneuron[i].enetwork[c]]++;
			 if(active == 0) for(c=0; c<vmhneuron[i].econnect; c++) (*activity)[vmhneuron[i].enetwork[c]]--;
			 actmutex->Unlock();   
			 }
			 
			 
			 if(esyn == 0) {
			 esynsum = 0;
			 isynsum = 0;
			 }
			 else {//actsynch(i);
			 */
			
			//esynsum = esynsum[i];
			isynsum = 0;
			
			//actmutex->Lock();                                   // **Activity Synch Section**
			//wxMutexLocker lock(*actmutex);
			//if((*activity)[i]) (*activity)[i] = (*activity)[i] - 1;
			//esynsum = (*activity)[i];
			
			//for(c=0; c<vmhneuron[i].econnect; c++) 
			//	if((*activity)[vmhneuron[i].enetwork[c]] == 1) esynsum++;
			
			
			//for(c=0; c<vmhneuron[i].iconnect; c++) 
			//	if((*activity)[vmhneuron[i].inetwork[c]] == 1) isynsum++;
			//actmutex->Unlock();                                      // **Synch End**
			
			//actsynch(i);
			//}
			
			esyninput = (synweight / vmhneurons) * esynsum[i];
			isyninput = (synweight / vmhneurons) * isynsum;
			
			re0 = dend0e;
			ri0 = dend0i;
			re1 = dend1e;
			ri1 = dend1i;
			
			//neudex = i;
			
			nepsp0 = 0;
			while(vmhneuron[i].epspt0 < hstep) {
				nepsp0++;
				vmhneuron[i].epspt0 = -log(1 - mrand01()) / dend0e + vmhneuron[i].epspt0;
			}
			vmhneuron[i].epspt0 = vmhneuron[i].epspt0 - hstep;
			
			nipsp0 = 0;
			while(vmhneuron[i].ipspt0 < hstep) {
				nipsp0++;
				vmhneuron[i].ipspt0 = -log(1 - mrand01()) / dend0i + vmhneuron[i].ipspt0;
			}
			vmhneuron[i].ipspt0 = vmhneuron[i].ipspt0 - hstep;
			
			/*
			 nepsp1 = 0;
			 while(vmhneuron[i].epspt1 < hstep) {
			 nepsp1++;
			 vmhneuron[i].epspt1 = -log(1 - mrand01()) / dend1e + vmhneuron[i].epspt1;
			 //randcall++;
			 }
			 vmhneuron[i].epspt1 = vmhneuron[i].epspt1 - hstep;
			 
			 nipsp1 = 0;
			 while(vmhneuron[i].ipspt1 < hstep) {
			 nipsp1++;
			 vmhneuron[i].ipspt1 = -log(1 - mrand01()) / dend1i + vmhneuron[i].ipspt1;
			 }
			 vmhneuron[i].ipspt1 = vmhneuron[i].ipspt1 - hstep;
			 */
			
			//fprintf(tofp, "psps generated\n");
			
			epsph = ae * (ve - vmhneuron[i].v);
			ipsph = ai * (vmhneuron[i].v - vi);
			
			input0 = epsph * nepsp0 - ipsph * nipsp0 ; 
			//input1 = epsph * nepsp1 - ipsph * nipsp1 ; 
			//input1 = epsph * dendinput[i][1];
			input1 = 0;
			sinput = epsph * esyninput - ipsph * isyninput;
			
			if(maxsyn && sinput > maxsyn) sinput = maxsyn;
			
			vmhneuron[i].v = vmhneuron[i].v - (vmhneuron[i].v - vmhneuron[i].vrest) / memtau * hstep + input0 + input1 + sinput + waveinput;
			
			//if(nettime < 100000) vmhneuron[i].inputrec[modstep] = input0;
			
			vmhneuron[i].tHAP = vmhneuron[i].tHAP - hstep * vmhneuron[i].tHAP / vmhneuron[i].tauHAP;
			
			vmhneuron[i].tAHP = vmhneuron[i].tAHP - hstep * vmhneuron[i].tAHP / vmhneuron[i].tauAHP;  
			
			vmhneuron[i].tDAP = vmhneuron[i].tDAP - hstep * vmhneuron[i].tDAP / vmhneuron[i].tauDAP;  
			
			vmhneuron[i].th = vmhneuron[i].th0 + vmhneuron[i].tHAP + vmhneuron[i].tAHP - vmhneuron[i].tDAP;
			
			//if(nettime < 10000 && i<10) {
			//	vmhneudat[i].inputrec[modstep] = sinput; //v[i];
			//	//vmhneuron[i].inputrec[modstep] = waveinput;
			//	vmhneudat[i].threshrec[modstep] = vmhneuron[i].th;
			//}
			
			//if(nettime < 100) {
			//	fprintf(tofp, "time %.1f  input0 %.2f  input1 %.2f  v %.2f  thresh %.2f\n", 
			//		nettime, input0, input1, vmhneuron[i].v, vmhneuron[i].th); 
			//	fprintf(tofp, "dendinput %.2f  sinput %.2f\n", dendinput[i][1], sinput);
			//	fprintf(tofp, "nepsp %d  epsph %.2f  nipsp %d  ipsph %.2f\n\n", nepsp0, epsph, nipsp0, ipsph);
			//}
			
			if(vmhneuron[i].v >= vmhneuron[i].th && vmhneuron[i].ttime >= absref) {
				//fprintf(tofp, "spike fired\n");
				//netisi[i][s[i]] = vmhneuron[i].ttime;
				//vmhneuron[i].isis[s[i]] = vmhneuron[i].ttime;
				vmhneuron[i].times[s[i]] = nettime/1000;
				actmutex->Lock();
				(*activity)[i] = syndelay;
				actmutex->Unlock();
				//active[i] = syndelay;
				s[i]++;
				
				vmhneuron[i].tHAP = vmhneuron[i].tHAP + vmhneuron[i].kHAP;
				vmhneuron[i].tAHP = vmhneuron[i].tAHP + vmhneuron[i].kAHP;
				vmhneuron[i].tDAP = vmhneuron[i].tDAP + vmhneuron[i].kDAP;
				if(vmhneuron[i].tDAP > maxdap) vmhneuron[i].tDAP = maxdap;
				
				vmhneuron[i].ttime = 0;
			}
		}
	}
	
	for(i=nfirst; i<=nlast; i++) {
		//neurospikenum[i] = s[i];
		vmhneuron[i].count = s[i];
		//fprintf(tofp, "neuron %d fired %d spikes at %.2fHz, time %.2f\n",
		//	i, s[i], s[i]/(maxtime/1000), vmhneuron[i].neurotime);
	}	
	actmutex->Lock();
	vmhmod->done++;
	actmutex->Unlock();
	
	fclose(tofp[threadID]);
	//networkdisp();
	//inputoutput(0);
}



void VMHspikegen::actsynch(int i)
{
	int c;
	//wxMutexLocker lock(*actmutex);
	actmutex->Lock();  
	
	if((*activity)[i]) (*activity)[i] = (*activity)[i] - 1;
	
	esynsum = 0;
	for(c=0; c<vmhneuron[i].econnect; c++) if((*activity)[vmhneuron[i].enetwork[c]] == 1) esynsum++;
	
	isynsum = 0;
	for(c=0; c<vmhneuron[i].iconnect; c++) if((*activity)[vmhneuron[i].inetwork[c]] == 1) isynsum++;
	
	actmutex->Unlock();  
}