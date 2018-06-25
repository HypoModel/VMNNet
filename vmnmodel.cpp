
/*
*  vmnmodel.cpp
*  HypoModel
*
*  Created by Duncan MacGregor
*  University of Edinburgh 2016
*
*/


#include "vmnmodel.h"
#include "vmnevofit.h"


VMNDat::VMNDat(int size)
{
	input1.setsize(size);
	input2.setsize(size);
	noisig.setsize(size);
}


VMNNeuron::VMNNeuron()
{
	inputrec.setsize(100000);
	raterec.setsize(10000);

	mag.resize(maxspikes);
}


VMNNeuroDat::VMNNeuroDat(int size)
{
	V.setsize(size);
	V2.setsize(size);
	HAP.setsize(size);
	AHP.setsize(size);
	netsyn.setsize(size);
	inpsyn.setsize(size);
	memsyn.setsize(size);
	memsyn2.setsize(size);
	sigsyn.setsize(size);
	sigsyn2.setsize(size);
	facB.setsize(size);
}


VMNModel::VMNModel(int type, wxString name, HypoMain *main)
	: NeuroMod(type, name, main)
{
	int i;

	if(mainwin->user || mainwin->basic) basicmode = 1;

	evodata = NULL;
	spikefitdata = NULL;
	fitboxdata = NULL;

	datsample = mainwin->datsample;
	storesize = 500000;

	mainwin->SetMinSize(wxSize(470, 690));
	mainwin->xstretch = 0;
	mainwin->mod = this;
	xmin = 0;
	path = "VMN";
	oldhist = false;

	vmhneuron = new VMNNeuron[200];   // 500
	numtypes = 3;
	prototypes = 2;

	currvmn = new SpikeDat();
	currvmn->diagbox = diagbox;

	netdat = new SpikeDat();
	netdat1 = new SpikeDat();
	netdat2 = new SpikeDat();
	netdat3 = new SpikeDat();
	vmndata = new VMNDat(storesize);
	analysisdata = new AnaDat();
	sigsim.setsize(1000100);
	raterec.setsize(1000100);
	raterec2.setsize(1000100);
	neurodata = new VMNNeuroDat(storesize);

	for(i=0; i<3; i++) rangedata[i].setsize(1000);
	rangeref.setsize(1000);

	currvmn->mainwin = mainwin;

	netbox = new VMNNetBox(this, "VMN Network", wxPoint(320, 0), wxSize(470, 610));
	signalbox = new SignalBox(this, "Signal Box", wxPoint(0, 300), wxSize(400, 500));
	protobox = new VMNProtoBox(this, "Protocol", wxPoint(0, 0), wxSize(320, 500));
	outbox = new OutBox(this, "Data Output", wxPoint(0, 0), wxSize(320, 500), 2000, 200);
	//cellbox = new CellBox(this, "Cell Data", wxPoint(0, 0), wxSize(320, 500));
	neurobox = new VMNNeuroBox(this, "VMN Neuron", wxPoint(0, 0), wxSize(320, 700));
	diagbox = mainwin->diagbox;
	modbox = netbox;

	(*toolflags)["spikebox"] = 1;		// Select universal model tools
	mainwin->ToolLoad();				// Load universal model tools

	EvoInit();

	//fitchrome->Output("chrometext.txt");
	
	modtools.AddBox(neurobox);
	modtools.AddBox(netbox);
	modtools.AddBox(signalbox, true);
	modtools.AddBox(protobox, true);
	modtools.AddBox(outbox, true);
	//modtools.AddBox(cellbox, true);
	modtools.AddBox(diagbox, false, true);
	modtools.AddBox(fitbox, true);

	neurobox->SetPanel(ID_EvoFit, fitbox); 
	neurobox->canclose = false;
	netbox->canclose = false;
	//conbox->SetWindowStyleFlag(wxFRAME_FLOAT_ON_PARENT | wxFRAME_TOOL_WINDOW | wxCAPTION | wxRESIZE_BORDER);

	ModLoad();
	for(i=0; i<modtools.numtools; i++) {
		modtools.box[i]->ReSize();
		modtools.box[i]->Show(modtools.box[i]->visible);
	}
	modbox->ParamLoad("default");
	graphload = true;

	graphbase->initpath = path;
	GraphData();
	EvoGraphs();

	graphbase->GetGraph("vmnhaz5ms")->ylabelplaces  = 2;
	graphbase->GetGraph("exphaz5ms")->ylabelplaces  = 2;

	gsync = 0;

	modbox->Show(true);
}


void VMNModel::EvoInit()
{
	// Evo Fit Data Storage
	spikefitdata = new SpikeFitDat();
	fitboxdata = new SpikeDat();
	fitboxdata->burstdata = new BurstDat();
	evodata = new SpikeDat();
	evodata->burstdata = new BurstDat();
	evodata->diagbox = diagbox;

	// Evo Fit Chrome initialise
	fitchrome = new EvoChrome();
	fitchrome->diagbox = diagbox;
	fitchrome->parambox = neurobox;
	fitchrome->AddParam("ire", netbox->GetCon("vmhinput1"), 200, 1000, true, "Syn Rate");
	fitchrome->AddParam("iratio", neurobox->GetCon("iratio"), 0, 1, false, "Iratio");
	fitchrome->AddParam("pspmag", neurobox->GetCon("pspmag"), 0.1, 10, false, "PSP mag");
	fitchrome->AddParam("halflifeSyn", neurobox->GetCon("halflife"), 2, 20, false, "Syn HL");
	fitchrome->AddParam("Vthresh", neurobox->GetCon("vthre"), -60, -40, false, "V Thresh");
	fitchrome->AddParam("kHAP", neurobox->GetCon("kHAP"), 10, 100, false, "HAP k");
	fitchrome->AddParam("halflifeHAP", neurobox->GetCon("halflifeHAP"), 2, 20, false, "HAP HL");
	fitchrome->AddParam("kDAP", neurobox->GetCon("kDAP"), 0, 5, false, "DAP k");
	fitchrome->AddParam("halflifeDAP", neurobox->GetCon("halflifeDAP"), 50, 500, false, "DAP HL");
	fitchrome->AddParam("kAHP", neurobox->GetCon("kAHP"), 0, 0.001, false, "AHP k");
	fitchrome->AddParam("halflifeAHP", neurobox->GetCon("halflifeAHP"), 100, 20000, false, "AHP HL");
	if(!basicmode) {
		fitchrome->AddParam("esynL1", netbox->GetCon("esynL1"), 0, 1, false, "esynL1");
		fitchrome->AddParam("esynweight", netbox->GetCon("synweightL1"), 0, 2, false, "weightL1");
		fitchrome->AddParam("syndelay", netbox->GetCon("syndelay"), 0, 10, false, "syndelay");
		fitchrome->AddParam("syndelrange", netbox->GetCon("syndelrange"), 0, 20, false, "sdrange");
	}
	fitchrome->Output("chrometest.txt");

	fitchrome->Diagnostic();
	
	fitbox = new EvoFitBox((Model *)this, fitchrome, "Evo Spike Fit", wxPoint(320, 455), wxSize(320, 430));
	fitbox->expdata = mainwin->expdata;
	fitbox->moddata = currvmn;
	fitbox->loaddata = fitboxdata; 
	fitbox->evodata = evodata;
	fitbox->spikefitdata = spikefitdata;
	fitbox->burstbox = mainwin->burstbox;
}


void VMNModel::EvoRun()
{
	fitbox->evothread = new EvoFitVMN(this, fitbox);
	fitbox->evothread->Create();
	fitbox->evothread->Run();
}


void VMNModel::Output()
{
	int i, j;
	TextFile outfile;
	wxString filename, filetag, text;
	wxString sectag, outdir;

	outdir = mainwin->outpath + "/output";

	for(i=0; i<netbox->vmhneurons; i++) {
		outbox->textgrid->SetCell(0, i, text.Format("m%d", i));
		diagbox->Write(text.Format("grid row %d : %d spikes\n", i, vmhneuron[i].spikecount));
		for(j=0; j<vmhneuron[i].spikecount; j++) {
			outbox->textgrid->SetCell(j+1, i, text.Format("%.5f", vmhneuron[i].times[j] / 1000));
		}
	}
}


void VMNModel::RangeOut(int datacount)
{
	int i, j;
	wxString text, tag, path, outline;
	TextFile outfile;

	path = GetPath() + "\\Output\\";
	if(!wxDirExists(path)) wxMkdir(path);

	tag = netbox->paramstoretag->GetValue();

	outfile.New(path + "rangedata-" + tag + ".txt");
	for(i=0; i<rangecount; i++) {
		outline.Printf("%.0f", rangeref[i]);
		for(j=0; j<datacount; j++) outline += text.Format("\t%.4f", rangedata[j][i]);
		outfile.WriteLine(outline);
		//outfile.WriteLine(text.Format("%.0f %.4f %.4f %.4f", rangeref[i], rangedata[0][i], rangedata[1][i], rangedata[2][i]));
	}
	outfile.Close();
}


void VMNModel::SigSim(SpikeDat *sdat)
{
	int i;
	double memtau;
	double v, epsph, nepsp;
	double syninput;
	double vrest = -62;
	double ae = 0.04;
	double ve = 38;
	double halflife;
	double pspheight;
	int revpots;

	ParamStore *sigparams = signalbox->GetParams();
	halflife = (*sigparams)["psphalflife"];
	vrest = (*sigparams)["vrest"];
	pspheight = (*sigparams)["pspheight"];
	double kB = (*sigparams)["kB"];

	revpots = 0;
	//revpots = (*signalbox->modflags)["revpots"];

	memtau = 1 / (log((double)2) / halflife);
	v = vrest;
	//memtau = 3.6067;

	for(i=0; i<1000000; i++) {
		nepsp = sdat->srate1.data[i];
		if(revpots) epsph = ae * (ve - v);
		else epsph = pspheight;
		//if(kB  && nepsp) epsph = sdat->mag[i];
		syninput = epsph * nepsp;
		//v = v - (v - vrest) * 0.1 + syninput;
		v = v - (v - vrest) / memtau + syninput;
		sdat->synsim[i] = v;
	}
}


void VMNModel::StoreClear()
{
	int i;

	for(i=0; i<storesize; i++) {
		vmndata->input1[i] = 0;
		vmndata->input2[i] = 0;
	}
}


void VMNModel::RunModel()
{
	if(mainwin->diagnostic) mainwin->SetStatusText("VMN Model Run");
	modthread = new VMNMod(this);
	modthread->Create();
	modthread->Run();
}


VMNModel::~VMNModel()
{
	
	delete vmndata;
	delete[] vmhneuron;
	delete currvmn; 
	//delete netdat;
	delete netdat1; 
	delete netdat2;
	delete netdat3;
	delete analysisdata; 
	delete neurodata;

	if(spikefitdata) delete spikefitdata;
	if(fitboxdata) delete fitboxdata;
	//if(evodata) delete evodata;
	//delete fitchrome;
}


void VMNModel::GraphData()
{
	GraphSet *graphset;

	// Data graphs
	//
	// GraphDat(data pointer, xfrom, xto, yfrom, yto, label string, graph type, bin size, colour)
	// ----------------------------------------------------------------------------------
	
	currvmn->GraphSet(graphbase, "VMN ", blue, 1, "vmn");
	netdat1->GraphSet(graphbase, "Net L1 ", blue, 1, "net1");
	netdat2->GraphSet(graphbase, "Net L2 ", green, 1, "net2");

	graphset = graphbase->NewSet("VMN Spikes", "vmnspikes");
	graphset->AddFlag("timeres", 1);
	graphset->AddFlag("nettog", 10);
	graphset->AddFlag("rateres", 100);
	graphset->Add(graphbase->tagindex["vmnrate1s"], 0);
	graphset->Add(graphbase->tagindex["vmnspikes1ms"], 1);
	graphset->Add(graphbase->tagindex["net1rate1s"], 10);
	graphset->Add(graphbase->tagindex["net1spikes1ms"], 11);
	graphset->Add("vmnrate10s", 100);
	if(diagbox) diagbox->textbox->AppendText(graphset->Display());

	graphset = graphbase->NewSet("VMN Intervals", "vmnintervals");
	graphset->AddFlag("nettog", 100);
	graphset->AddFlag("hazmode1", 10);
	graphset->AddFlag("binrestog1", 1);
	graphset->AddFlag("normtog", 1000);
	graphset->Add("vmnhist1ms", 0);
	graphset->Add("vmnhaz1ms", 10);
	graphset->Add("vmnhaz1ms", 1010);
	graphset->Add("vmnhist5ms", 1);
	graphset->Add("vmnhaz5ms", 11);
	graphset->Add("vmnhaz5ms", 1011);
	graphset->Add("net1hist1ms", 100);
	graphset->Add("net1haz1ms", 110);
	graphset->Add("net1hist5ms", 101);
	graphset->Add("net1haz5ms", 111);
	graphset->Add("vmnnormhist1ms", 1000);
	graphset->Add("vmnnormhist5ms", 1001);
	if(diagbox) diagbox->textbox->AppendText(graphset->Display());

	graphset = graphbase->NewSet("L2 Spikes", "l2spikes");
	graphset->AddFlag("timeres", 1);
	graphset->Add("net2rate1s", 0);
	graphset->Add("net2spikes1ms", 1);
	if(diagbox) diagbox->textbox->AppendText(graphset->Display());

	graphset = graphbase->NewSet("L2 Intervals", "l2intervals");
	graphset->AddFlag("hazmode1", 10);
	graphset->AddFlag("binrestog1", 1);
	graphset->Add("net2hist1ms", 0);
	graphset->Add("net2haz1ms", 10);
	graphset->Add("net2hist5ms", 1);
	graphset->Add("net2haz5ms", 11);
	if(diagbox) diagbox->textbox->AppendText(graphset->Display());

	graphbase->Add(GraphDat(&netdat1->synsim, 0, 1, -80, 20, "L1 Syn Sim", 4, 0.001, green), "l1synsim");
	graphbase->Add(GraphDat(&netdat2->synsim, 0, 1, -80, 20, "L2 Syn Sim", 4, 0.001, lightblue), "l2synsim");
	graphbase->Add(GraphDat(&currvmn->synsim, 0, 1, -80, 20, "Net Input Sim", 4, 0.001, green), "netinputsim");
	graphbase->Add(GraphDat(&currvmn->netinputrec, 0, 1000, -80, 20, "Net Input", 4, 0.001, green), "netinput");
	graphbase->Add(GraphDat(&currvmn->raterec, 0, 1000, 0, 500, "Input Rate", 4, 1, green), "inputrate");
	graphbase->Add(GraphDat(&analysisdata->autocorr, 0, 1000, 0, 1000, "Auto Correlation", 1, 10, green), "autocorr");
	graphbase->Add(GraphDat(&raterec, 0, 1000, 0, 1000, "Input Rate L1", 4, 0.001, blue), "inputrate1");
	graphbase->Add(GraphDat(&raterec2, 0, 1000, 0, 1000, "Input Rate L2", 4, 0.001, blue), "inputrate2");

	graphbase->Add(GraphDat(&vmndata->input1, 0, 500, 0, 500, "Input L1", 4, 1, lightblue, 1000/datsample), "input1");
	graphbase->Add(GraphDat(&vmndata->input2, 0, 500, 0, 500, "Input L2", 4, 1, lightblue, 1000/datsample), "input2");
	graphbase->Add(GraphDat(&vmndata->noisig, 0, 500, 0, 500, "Noise Signal", 4, 1, lightred, 1000/datsample), "noisig");

	graphbase->Add(GraphDat(&neurodata->V, 0, 1000, 0, 100, "Neuron V", 4, 0.001, green), "neurov");
	graphbase->Add(GraphDat(&neurodata->HAP, 0, 1000, 0, 100, "Neuron HAP", 4, 0.001, green), "neurohap");
	graphbase->Add(GraphDat(&neurodata->AHP, 0, 1000, 0, 100, "Neuron AHP", 4, 0.001, green), "neuroahp");
	graphbase->Add(GraphDat(&neurodata->netsyn, 0, 1000, 0, 100, "Neuron NetSyn", 5, 0.001, green), "neuronetsyn");
	graphbase->Add(GraphDat(&neurodata->inpsyn, 0, 1000, 0, 100, "Neuron InpSyn", 5, 0.001, green), "neuroinpsyn");
	graphbase->Add(GraphDat(&neurodata->memsyn, 0, 1000, 0, 100, "Neuron MemSyn L1", 5, 0.001, green), "neuromemsyn");
	graphbase->Add(GraphDat(&neurodata->memsyn2, 0, 1000, 0, 100, "Neuron MemSyn L2", 5, 0.001, blue), "neuromemsyn2");
	graphbase->Add(GraphDat(&neurodata->sigsyn2, 0, 1000, 0, 100, "Neuron SigSyn L2", 5, 0.001, red), "neurosigsyn2");
	graphbase->Add(GraphDat(&neurodata->facB, 0, 1000, 0, 100, "Neuron FacB", 5, 0.001, green), "neurofacb");

	graphbase->Add(GraphDat(&rangedata[0], 0, 1000, 0, 100, "L1 Range Freq", 2, 1, green), "l1range");
	graphbase->Add(GraphDat(&rangedata[1], 0, 1000, 0, 100, "L2 Range Freq", 2, 1, green), "l2range");
	graphbase->Add(GraphDat(&rangedata[2], 0, 1000, 0, 100, "L3 Range Freq", 2, 1, green), "l3range");

	graphbase->GetGraph("l1range")->gdatax = &rangeref;
	graphbase->GetGraph("l2range")->gdatax = &rangeref;
	graphbase->GetGraph("l3range")->gdatax = &rangeref;

	graphbase->Add(GraphDat(&currvmn->IoDdata, 0, 70, 0, 2, "VMN IoD", 9, 1, lightgreen), "iodvmn");
	graphbase->GetGraph("iodvmn")->gdatax = &currvmn->IoDdataX;
	graphbase->GetGraph("iodvmn")->xcount = 7;  
	graphbase->GetGraph("iodvmn")->synchx = false; 
	graphbase->GetGraph("iodvmn")->barshift  = 20;
	graphbase->GetGraph("iodvmn")->ylabelplaces  = 0;

	graphbase->Add(GraphDat(&currvmn->meanV, 0, 100, 0, 300, "Mean Spike Form", 4, 1, lightblue), "meanV");

	gcodes[0] = "vmnspikes";
	gcodes[1] = "vmnintervals";
	gcodes[2] = "l2spikes";
	gcodes[3] = "autocorr";
	gcodes[4] = "neurov";
	gcodes[5] = "neuroinpsyn";

	gcount = 6;
	gsmode = 1;
}


void VMNModel::GSwitch(GraphDisp *gpos, ParamStore *gflags)
{
	int i, gdex;
	GraphSet *graphset;
	wxString text;

	if(diagbox) diagbox->textbox->AppendText("\n");

	for(i=0; i<gcount; i++) {
		graphset = graphbase->GetSet(gcodes[i]);
		gdex = graphset->GetPlot(gflags);
		if(diagbox) diagbox->textbox->AppendText(text.Format("gpos %d   gcode %s   set %s   plot %d   modesum %d   sdex %d\n", 
			i, gcodes[i], graphset->tag, gdex, graphset->modesum, graphset->sdex));
		gpos[i].Front((*graphbase)[gdex]);
		gpos[i].sdex = graphset->sdex;
	}
}

