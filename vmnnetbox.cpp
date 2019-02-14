/*
*  vmnnetbox.cpp
*  HypoModel
*
*  Created by Duncan MacGregor
*  University of Edinburgh 2018
*  Released under MIT license, see https://opensource.org/licenses/MIT
*
*/

#include "vmnmodel.h"


// Network control box
VMNNetBox::VMNNetBox(VMNModel *vmnmodel, const wxString& title, const wxPoint& pos, const wxSize& size)
: ParamBox(vmnmodel, title, pos, size, "VMNNET")
{
	column = 0;
	labelwidth = 60;

	wxSize toolsize;
	wxColour toolback, textcolour;
	boxname = "VMN";
	mod = vmnmodel;
	mainwin = mod->mainwin;

	neuroindex = 0;
	vmhneurons = 50;                 // 100
	//buttonheight = 23;
	sumflag = 0;
	toolsize.x = 1;
	toolsize.y = 1;
	toolback.Set("#ddddff");       // 8080ff light blue
	textcolour.Set("#000000");

	neurons = mod->vmhneuron;
	netdat = mod->netdat;
	wxFont boxfont(8, wxFONTFAMILY_SWISS, wxNORMAL, wxNORMAL, false, "Tahoma");
	panel->SetFont(boxfont);

	InitMenu();

	if(mod->basicmode) PanelBasic();
	else PanelFull();
}


void VMNNetBox::PanelFull() 
{
	SetConFlag(ID_autosum, "autosum", "Auto Net Sum", 1);

	SetModFlag(ID_synqueue, "synqueue", "Syn Queue", 1); 
	SetModFlag(ID_gpunetmode, "GPUnetmode", "Net Fit Mode", 0); 
	SetModFlag(ID_cpumode, "cpumode", "CPU Fit Mode", 0); 
	SetModFlag(ID_Diag, "GPUdiagnostic", "GPU Diagnostic", 0); 
	SetModFlag(ID_formfilter, "formfilter", "Mean Spike Filter", 0); 

	paramset->AddCon("neuronsL1", "Neurons L1", vmhneurons, 1, 0, labelwidth);
	//paramset->AddCon("vmhconnect", "Synapse", 0, 1, 0, labelwidth);
	paramset->AddCon("vmhinput1", "Input L1", 300, 1, 1, labelwidth);
	paramset->AddCon("synweightL1", "Synaptic L1", 1.0, 1, 2, labelwidth);
	paramset->AddCon("absref", "Abs Ref", 5, 1, 1, labelwidth);
	paramset->AddCon("inputcycle", "Inp Cycle", 0, 1, 1, labelwidth);
	paramset->AddCon("waveamp", "Inp Amp", 0, 1, 1, labelwidth);
	paramset->AddCon("syndelay", "Syn Delay", 3, 1, 1, labelwidth);
	paramset->AddCon("syndelrange", "Syn Range", 0, 1, 1, labelwidth);
	paramset->AddCon("esynL1", "E Syn L1", 0.7, 0.01, 2, labelwidth);
	paramset->AddCon("esynsd", "E Syn SD", 0, 0.01, 2, labelwidth);
	paramset->AddCon("isynL1", "I Syn L1", 0, 0.01, 2, labelwidth);
	paramset->AddCon("vrestsd", "Vrest SD", 0, 0.1, 2, labelwidth);
	paramset->AddCon("kHAPsd", "kHAP SD", 0, 0.1, 2, labelwidth);
	paramset->AddCon("tauHAPsd", "tauHAP SD", 0, 0.1, 2, labelwidth);
	paramset->AddCon("maxsyn", "Max Syn", 0, 1, 0, labelwidth);
	paramset->AddCon("esyntrans", "E Syn Trans", 0.5, 0.01, 2, labelwidth);
	paramset->AddCon("inputsd", "Input SD", 0, 0.01, 2, labelwidth);

	paramset->AddCon("neuronsL2", "Neurons L2", 0, 1, 0, labelwidth);
	paramset->AddCon("esynL2", "E Syn L2", 0.5, 0.01, 2, labelwidth);
	paramset->AddCon("synweightL2", "Synaptic L2", 1.0, 1, 1, labelwidth);
	paramset->AddCon("esynL12", "E Syn L12", 0.35, 0.01, 2, labelwidth);
	paramset->AddCon("synweightL12", "Synap L12", 1, 0.1, 1, labelwidth);
	paramset->AddCon("esynL21", "E Syn L21", 0.5, 0.01, 2, labelwidth);
	paramset->AddCon("synweightL21", "Synap L21", 1.0, 0.1, 1, labelwidth);
	paramset->AddCon("vmhinput2", "Input L2", 100, 1, 1, labelwidth);
	paramset->AddCon("vrest2", "V Rest L2", -62, 0.1, 2, labelwidth);
	paramset->AddCon("kHAP2", "HAP k L2", 60, 0.1, 2, labelwidth);
	paramset->AddCon("halflifeHAP2", "HAP HL L2", 5, 0.1, 2, labelwidth);
	paramset->AddCon("kDAP2", "L2 DAP k", 0, 0.1, 2, labelwidth);
	paramset->AddCon("halflifeDAP2", "L2 DAP HL", 500, 10, 2, labelwidth);
	paramset->AddCon("synhl", "Syn HL", 2.5, 0.1, 2, labelwidth);
	paramset->AddCon("synmag", "Syn Mag", 4, 0.1, 2, labelwidth);

	paramset->AddCon("neuronsL3", "Neurons L3", 0, 1, 0, labelwidth);
	paramset->AddCon("esynL3", "E Syn L3", 0.5, 0.01, 2, labelwidth);
	paramset->AddCon("synweightL3", "Synaptic L3", 1.0, 1, 1, labelwidth);
	paramset->AddCon("esynL23", "E Syn L23", 0.35, 0.01, 2, labelwidth);
	paramset->AddCon("synweightL23", "Synap L23", 1, 0.1, 1, labelwidth);
	paramset->AddCon("vmhinput3", "Input L3", 100, 1, 1, labelwidth);
	paramset->AddCon("vrest3", "V Rest L3", -62, 0.1, 2, labelwidth);
	paramset->AddCon("kHAP3", "HAP k L3", 60, 0.1, 2, labelwidth);
	paramset->AddCon("halflifeHAP3", "HAP HL L3", 5, 0.1, 2, labelwidth);

	ParamLayout(2);

	int datwidth = 50;

	runcount = new wxStaticText(panel, -1, wxT("---"), wxDefaultPosition, wxSize(40, -1), wxALIGN_CENTRE|wxBORDER_RAISED|wxST_NO_AUTORESIZE);

	mean = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);
	freq = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);
	sd = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);
	spikes = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);
	esyn = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);
	vrest = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);
	kHAP = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);
	tauHAP = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);
	input = new wxStaticText(panel, -1, "0", wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);

	wxGridSizer *datagrid = new wxGridSizer(2, 5, 5);
	datagrid->Add(new wxStaticText(panel, -1, "Spikes"), 0, wxALIGN_CENTRE);
	datagrid->Add(spikes);
	datagrid->Add(new wxStaticText(panel, -1, "Freq"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	datagrid->Add(freq);
	datagrid->Add(new wxStaticText(panel, -1, "Mean"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	datagrid->Add(mean);
	datagrid->Add(new wxStaticText(panel, -1, "Std Dev"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	datagrid->Add(sd);
	datagrid->AddSpacer(2);
	datagrid->AddSpacer(2);
	datagrid->Add(new wxStaticText(panel, -1, "E Syn"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	datagrid->Add(esyn);
	datagrid->Add(new wxStaticText(panel, -1, "Vrest"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	datagrid->Add(vrest);
	datagrid->Add(new wxStaticText(panel, -1, "kHAP"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	datagrid->Add(kHAP);
	datagrid->Add(new wxStaticText(panel, -1, "tauHAP"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	datagrid->Add(tauHAP);
	datagrid->Add(new wxStaticText(panel, -1, "Input"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	datagrid->Add(input);

	datneuron = new wxTextCtrl(panel, ID_datneuro, "---", wxDefaultPosition, wxSize(50, -1), wxALIGN_LEFT|wxBORDER_SUNKEN|wxST_NO_AUTORESIZE|wxTE_PROCESS_ENTER);
	datspin = new wxSpinButton(panel, wxID_ANY, wxDefaultPosition, wxSize(40, 17), wxSP_HORIZONTAL|wxSP_ARROW_KEYS);
	wxButton *sumbutton = new wxButton(panel, ID_sum, "Sum", wxDefaultPosition, wxSize(30, 17));
	wxBoxSizer *datbox = new wxBoxSizer(wxHORIZONTAL);
	datbox->Add(datspin, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);
	datbox->AddSpacer(5);
	datbox->Add(sumbutton, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);

	wxBoxSizer *neurobox = new wxBoxSizer(wxHORIZONTAL);
	neurobox->Add(new wxStaticText(panel, wxID_ANY, "Neuron"), 1, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	neurobox->Add(datneuron, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);

	wxStaticBoxSizer *databox = new wxStaticBoxSizer(wxVERTICAL, panel, "");
	databox->AddSpacer(2);
	databox->Add(neurobox, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);
	databox->AddSpacer(5);
	databox->Add(datbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	databox->AddSpacer(5);
	databox->Add(datagrid, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);

	wxBoxSizer *runbox = new wxBoxSizer(wxHORIZONTAL);
	runbox->Add(new wxStaticText(panel, -1, "Run :"), 0, wxALIGN_CENTRE_VERTICAL|wxRIGHT, 5);
	runbox->Add(runcount);
	runbutton = new wxButton(panel, ID_Run, "Run Model", wxDefaultPosition, wxSize(60, buttonheight));
	//wxButton *prunbutton = new wxButton(panel, ID_pararun, "Para Run", wxDefaultPosition, wxSize(60, buttonheight));
	buttonbox = new wxBoxSizer(wxHORIZONTAL);
	buttonbox->Add(runbutton);
	buttonbox->AddSpacer(5);
	//buttonbox->Add(prunbutton);

	wxBoxSizer *seedbox = new wxBoxSizer(wxHORIZONTAL);
	paramset->AddNum("modseed", "", 0, 0, 0, 80);
	paramset->GetCon("modseed")->SetMinMax(0, 1000000000000);
	seedbox->Add(paramset->GetCon("modseed"), 0, wxALIGN_CENTRE_VERTICAL);
	seedcheck = SetModCheck(ID_seedcheck, "seedgen", "Seed", true);
	seedbox->Add(seedcheck, 0, wxALIGN_CENTRE_VERTICAL|wxLEFT, 5);

	wxStaticBoxSizer *netbox = new wxStaticBoxSizer(wxVERTICAL, panel, "Network Generation");
	netcheck = new wxCheckBox(panel, ID_netcheck, "Network");
	netcheck->SetValue(true);
	cellcheck = new wxCheckBox(panel, ID_cellcheck, "Cells");
	cellcheck->SetValue(true);
	unicheck = new wxCheckBox(panel, ID_cellcheck, "Uniform");
	unicheck->SetValue(false);
	storetag = new wxTextCtrl(panel, wxID_ANY, "l1n10", wxDefaultPosition, wxSize(100, -1));
	wxButton *storebutton = new wxButton(panel, ID_netstore, "Store", wxDefaultPosition, wxSize(40, buttonheight));
	wxButton *loadbutton = new wxButton(panel, ID_netload, "Load", wxDefaultPosition, wxSize(40, buttonheight));

	wxBoxSizer *checks = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer *checks2 = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer *netbuttons = new wxBoxSizer(wxHORIZONTAL);
	checks->Add(netcheck, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5); 
	checks->Add(cellcheck, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5); 
	checks2->Add(unicheck, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5); 
	netbox->Add(checks, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 3);
	netbox->Add(checks2, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 3);
	netbox->AddSpacer(5);
	netbox->Add(storetag, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 2);
	netbuttons->Add(storebutton, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 2); 
	netbuttons->Add(loadbutton, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 2); 
	netbox->Add(netbuttons, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 2);

	//wxBoxSizer *paramfilebox = StoreBox("ref2b");
	//synccheck = new wxCheckBox(panel, wxID_ANY, "Sync");
	//synccheck->SetValue(true);
	//paramfilebox->Add(synccheck, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 2);

	wxBoxSizer *storebox = StoreBoxSync();

	wxBoxSizer *rightbox = new wxBoxSizer(wxVERTICAL);
	rightbox->Add(netbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);
	rightbox->AddStretchSpacer(10);
	rightbox->Add(buttonbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);
	rightbox->Add(runbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);
	rightbox->AddStretchSpacer(10);
	rightbox->Add(seedbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);
	rightbox->AddStretchSpacer(10);
	rightbox->Add(databox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);

	wxBoxSizer *leftbox = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer *footbox = new wxBoxSizer(wxHORIZONTAL);
	leftbox->Add(parambox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	footbox->Add(storebox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	AddButton(ID_Compare, "Comp", 40, footbox);
	
	/*
	popfreq = NumPanel(50, wxALIGN_RIGHT, "0");
	ratemean = NumPanel(50, wxALIGN_RIGHT, "0");

	wxGridSizer *popdatagrid = new wxFlexGridSizer(2, 3, 3);
	popdatagrid->Add(TextLabel("Freq"), 0, wxALIGN_CENTRE);
	popdatagrid->Add(popfreq);
	popdatagrid->Add(TextLabel("SD"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	popdatagrid->Add(popsd);
	popdatagrid->Add(TextLabel("Resp"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	popdatagrid->Add(ratemean);

	wxStaticBoxSizer *popdatabox = new wxStaticBoxSizer(wxVERTICAL, panel, "Population");
	popdatabox->AddSpacer(5);
	popdatabox->Add(popdatagrid, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 2);*/

	freqL1 = NumPanel(50, wxALIGN_RIGHT, "0");
	freqL2 = NumPanel(50, wxALIGN_RIGHT, "0");
	freqL3 = NumPanel(50, wxALIGN_RIGHT, "0");

	wxGridSizer *laydatagrid = new wxFlexGridSizer(2, 3, 3);
	laydatagrid->Add(TextLabel("Freq L1"), 0, wxALIGN_CENTRE);
	laydatagrid->Add(freqL1);
	laydatagrid->Add(TextLabel("Freq L2"), 0, wxALIGN_CENTRE);
	laydatagrid->Add(freqL2);
	laydatagrid->Add(TextLabel("Freq L3"), 0, wxALIGN_CENTRE);
	laydatagrid->Add(freqL3);

	wxStaticBoxSizer *laydatabox = new wxStaticBoxSizer(wxVERTICAL, panel, "Layer Data");
	laydatabox->AddSpacer(5);
	laydatabox->Add(laydatagrid, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 2);

	rightbox->AddStretchSpacer(5);
	rightbox->Add(laydatabox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);
	rightbox->AddSpacer(5);

	leftbox->Add(footbox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);

	wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);
	hbox->Add(leftbox, 0, wxALL, 0);
	hbox->AddStretchSpacer();
	hbox->Add(rightbox, 0, wxALIGN_CENTRE_VERTICAL|wxALL, 0);

	status = StatusBar();
	wxBoxSizer *statusbox = new wxBoxSizer(wxHORIZONTAL);
	statusbox->Add(status, 1, wxEXPAND);

	mainbox->AddSpacer(5);
	mainbox->Add(hbox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	mainbox->AddStretchSpacer(5);
	mainbox->Add(statusbox, 0, wxEXPAND);

	panel->Layout();

	Connect(ID_Run, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnRun));
	Connect(ID_pararun, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnRun));
	Connect(ID_sum, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnSum));
	Connect(ID_netstore, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnStore));
	Connect(ID_netload, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnLoad));
	Connect(ID_output, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnOutput));
	Connect(ID_paramstore, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnParamStore));
	Connect(ID_paramload, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnParamLoad));
	Connect(ID_Compare, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnParamLoad));
	Connect(ID_controls, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnControlsMenu));
	Connect(ID_model, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnModelMenu));

	Connect(wxEVT_COMMAND_TEXT_ENTER, wxCommandEventHandler(VMNNetBox::OnRun));
	Connect(wxEVT_SCROLL_LINEUP, wxSpinEventHandler(VMNNetBox::OnNeuroNext));
	Connect(wxEVT_SCROLL_LINEDOWN, wxSpinEventHandler(VMNNetBox::OnNeuroPrev));
	Connect(wxEVT_SCROLL_THUMBTRACK, wxSpinEventHandler(VMNNetBox::OnSpin));
}


void VMNNetBox::PanelBasic()
{
	SetConFlag(ID_autosum, "autosum", "Auto Net Sum", 1);

	SetModFlag(ID_synqueue, "synqueue", "Syn Queue", 1); 
	SetModFlag(ID_gpunetmode, "GPUnetmode", "Net Fit Mode", 0); 
	SetModFlag(ID_cpumode, "cpumode", "CPU Fit Mode", 0); 
	SetModFlag(ID_Diag, "GPUdiagnostic", "GPU Diagnostic", 0); 
	SetModFlag(ID_formfilter, "formfilter", "Mean Spike Filter", 0); 
	SetModFlag(ID_fixeddelay, "fixeddelay", "Fixed Connect Delay", 0);

	paramset->AddCon("neuronsL1", "Neurons L1", vmhneurons, 1, 0, labelwidth);
	//paramset->AddCon("vmhconnect", "Synapse", 0, 1, 0, labelwidth);
	paramset->AddCon("vmhinput1", "Input L1", 300, 1, 1, labelwidth);
	paramset->AddCon("synweightL1", "Synaptic L1", 1.0, 1, 2, labelwidth);
	//paramset->AddCon("inputcycle", "Inp Cycle", 0, 1, 1, labelwidth);
	//paramset->AddCon("waveamp", "Inp Amp", 0, 1, 1, labelwidth);
	paramset->AddCon("syndelay", "Syn Delay", 3, 1, 1, labelwidth);
	paramset->AddCon("syndelrange", "Syn Range", 0, 1, 1, labelwidth);
	paramset->AddCon("esynL1", "E Syn L1", 0.7, 0.01, 2, labelwidth);
	//paramset->AddCon("esynsd", "E Syn SD", 0, 0.01, 2, labelwidth);
	//paramset->AddCon("isynL1", "I Syn L1", 0, 0.01, 2, labelwidth);
	//paramset->AddCon("vrestsd", "Vrest SD", 0, 0.1, 2, labelwidth);
	if(mod->revisionmode) {
		paramset->AddCon("kHAPsd", "kHAP SD", 0, 0.1, 2, labelwidth);
		paramset->AddCon("tauHAPsd", "tauHAP SD", 0, 0.1, 2, labelwidth);
	}
	//paramset->AddCon("maxsyn", "Max Syn", 0, 1, 0, labelwidth);
	paramset->AddCon("esyntrans", "E Syn Trans", 0.5, 0.01, 2, labelwidth);
	//paramset->AddCon("inputsd", "Input SD", 0, 0.01, 2, labelwidth);

	paramset->AddCon("neuronsL2", "Neurons L2", 10, 1, 0, labelwidth);
	if(mod->revisionmode) {
		paramset->AddCon("esynL2", "E Syn L2", 0, 0.01, 2, labelwidth);
		paramset->AddCon("synweightL2", "Synaptic L2", 1.0, 1, 1, labelwidth);
	}
	paramset->AddCon("esynL12", "E Syn L12", 0.35, 0.01, 2, labelwidth);
	paramset->AddCon("synweightL12", "Synap L12", 1, 0.1, 1, labelwidth);
	if(mod->revisionmode) {
		paramset->AddCon("esynL21", "E Syn L21", 0, 0.01, 2, labelwidth);
		paramset->AddCon("synweightL21", "Synap L21", 1.0, 0.1, 1, labelwidth);
	}
	paramset->AddCon("vmhinput2", "Input L2", 100, 1, 1, labelwidth);
	paramset->AddCon("vrest2", "V Rest L2", -62, 0.1, 2, labelwidth);
	paramset->AddCon("kHAP2", "HAP k L2", 60, 0.1, 2, labelwidth);
	paramset->AddCon("halflifeHAP2", "HAP HL L2", 5, 0.1, 2, labelwidth);
	paramset->AddCon("kDAP2", "DAP k L2", 0, 0.1, 2, labelwidth);
	paramset->AddCon("halflifeDAP2", "DAP HL L2", 500, 10, 2, labelwidth);
	paramset->AddCon("synhl", "Syn HL", 2.5, 0.1, 2, labelwidth);
	paramset->AddCon("synmag", "Syn Mag", 4, 0.1, 2, labelwidth);

	//paramset->AddCon("neuronsL3", "Neurons L3", 10, 1, 0, labelwidth);
	//paramset->AddCon("esynL3", "E Syn L3", 0.5, 0.01, 2, labelwidth);
	//paramset->AddCon("synweightL3", "Synaptic L3", 1.0, 1, 1, labelwidth);
	//paramset->AddCon("esynL23", "E Syn L23", 0.35, 0.01, 2, labelwidth);
	//paramset->AddCon("synweightL23", "Synap L23", 1, 0.1, 1, labelwidth);
	//paramset->AddCon("vmhinput3", "Input L3", 100, 1, 1, labelwidth);
	//paramset->AddCon("vrest3", "V Rest L3", -62, 0.1, 2, labelwidth);
	//paramset->AddCon("kHAP3", "HAP k L3", 60, 0.1, 2, labelwidth);
	//paramset->AddCon("halflifeHAP3", "HAP HL L3", 5, 0.1, 2, labelwidth);

	ParamLayout(2);

	int datwidth = 50;

	runcount = new wxStaticText(panel, -1, wxT("---"), wxDefaultPosition, wxSize(40, -1), wxALIGN_CENTRE|wxBORDER_RAISED|wxST_NO_AUTORESIZE);

	mean = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);
	freq = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);
	sd = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);
	spikes = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);
	//esyn = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);
	//vrest = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);
	//kHAP = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);
	//tauHAP = new wxStaticText(panel, -1, wxT("0"), wxDefaultPosition, wxSize(datwidth, -1), wxALIGN_RIGHT|wxBORDER_RAISED|wxST_NO_AUTORESIZE);

	wxGridSizer *datagrid = new wxGridSizer(2, 5, 5);
	datagrid->Add(new wxStaticText(panel, -1, "Spikes"), 0, wxALIGN_CENTRE);
	datagrid->Add(spikes);
	datagrid->Add(new wxStaticText(panel, -1, "Freq"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	datagrid->Add(freq);
	datagrid->Add(new wxStaticText(panel, -1, "Mean"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	datagrid->Add(mean);
	datagrid->Add(new wxStaticText(panel, -1, "Std Dev"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	datagrid->Add(sd);
	//datagrid->AddSpacer(5);
	//datagrid->AddSpacer(5);
	//datagrid->Add(new wxStaticText(panel, -1, "E Syn"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	//datagrid->Add(esyn);
	//datagrid->Add(new wxStaticText(panel, -1, "Vrest"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	//datagrid->Add(vrest);
	//datagrid->Add(new wxStaticText(panel, -1, "kHAP"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	//datagrid->Add(kHAP);
	//datagrid->Add(new wxStaticText(panel, -1, "tauHAP"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	//datagrid->Add(tauHAP);

	datneuron = new wxTextCtrl(panel, ID_datneuro, "---", wxDefaultPosition, wxSize(50, -1), wxALIGN_LEFT|wxBORDER_SUNKEN|wxST_NO_AUTORESIZE|wxTE_PROCESS_ENTER);
	datspin = new wxSpinButton(panel, wxID_ANY, wxDefaultPosition, wxSize(40, 17), wxSP_HORIZONTAL|wxSP_ARROW_KEYS);
	wxButton *sumbutton = new wxButton(panel, ID_sum, "Sum", wxDefaultPosition, wxSize(30, 17));
	wxBoxSizer *datbox = new wxBoxSizer(wxHORIZONTAL);
	datbox->Add(datspin, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);
	datbox->AddSpacer(5);
	datbox->Add(sumbutton, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);

	wxBoxSizer *neurobox = new wxBoxSizer(wxHORIZONTAL);
	neurobox->Add(new wxStaticText(panel, wxID_ANY, "Neuron"), 1, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	neurobox->Add(datneuron, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);

	wxStaticBoxSizer *databox = new wxStaticBoxSizer(wxVERTICAL, panel, "");
	databox->AddSpacer(2);
	databox->Add(neurobox, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);
	databox->AddSpacer(5);
	databox->Add(datbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	databox->AddSpacer(5);
	databox->Add(datagrid, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);

	wxBoxSizer *runbox = new wxBoxSizer(wxHORIZONTAL);
	runbox->Add(new wxStaticText(panel, -1, "Run :"), 0, wxALIGN_CENTRE_VERTICAL|wxRIGHT, 5);
	runbox->Add(runcount);
	runbutton = new wxButton(panel, ID_Run, "Run Model", wxDefaultPosition, wxSize(60, buttonheight));
	//wxButton *prunbutton = new wxButton(panel, ID_pararun, "Para Run", wxDefaultPosition, wxSize(60, buttonheight));
	buttonbox = new wxBoxSizer(wxHORIZONTAL);
	buttonbox->Add(runbutton);
	buttonbox->AddSpacer(5);
	//buttonbox->Add(prunbutton);

	wxBoxSizer *seedbox = new wxBoxSizer(wxHORIZONTAL);
	paramset->AddNum("modseed", "", 0, 0, 0, 80);
	paramset->GetCon("modseed")->SetMinMax(0, 1000000000000);
	seedbox->Add(paramset->GetCon("modseed"), 0, wxALIGN_CENTRE_VERTICAL);
	seedcheck = SetModCheck(ID_seedcheck, "seedgen", "Seed", true);
	seedbox->Add(seedcheck, 0, wxALIGN_CENTRE_VERTICAL|wxLEFT, 5);

	/*
	wxStaticBoxSizer *netbox = new wxStaticBoxSizer(wxVERTICAL, panel, "Network Generation");
	netcheck = new wxCheckBox(panel, ID_netcheck, "Network");
	netcheck->SetValue(true);
	cellcheck = new wxCheckBox(panel, ID_cellcheck, "Cells");
	cellcheck->SetValue(true);
	unicheck = new wxCheckBox(panel, ID_cellcheck, "Uniform");
	unicheck->SetValue(false);
	storetag = new wxTextCtrl(panel, wxID_ANY, "l1n10", wxDefaultPosition, wxSize(100, -1));
	wxButton *storebutton = new wxButton(panel, ID_netstore, "Store", wxDefaultPosition, wxSize(40, buttonheight));
	wxButton *loadbutton = new wxButton(panel, ID_netload, "Load", wxDefaultPosition, wxSize(40, buttonheight));

	wxBoxSizer *checks = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer *checks2 = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer *netbuttons = new wxBoxSizer(wxHORIZONTAL);
	checks->Add(netcheck, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5); 
	checks->Add(cellcheck, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5); 
	checks2->Add(unicheck, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5); 
	netbox->Add(checks, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 3);
	netbox->Add(checks2, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 3);
	netbox->AddSpacer(5);
	netbox->Add(storetag, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 2);
	netbuttons->Add(storebutton, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 2); 
	netbuttons->Add(loadbutton, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 2); 
	netbox->Add(netbuttons, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 2);*/

	//wxBoxSizer *paramfilebox = StoreBox("ref2b");
	//synccheck = new wxCheckBox(panel, wxID_ANY, "Sync");
	//synccheck->SetValue(true);
	//paramfilebox->Add(synccheck, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 2);

	wxBoxSizer *storebox = StoreBoxSync();

	wxBoxSizer *rightbox = new wxBoxSizer(wxVERTICAL);
	//rightbox->Add(netbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);
	//rightbox->AddStretchSpacer(10);
	rightbox->Add(buttonbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);
	rightbox->Add(runbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);
	rightbox->AddStretchSpacer(10);
	rightbox->Add(seedbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);
	rightbox->AddStretchSpacer(10);
	rightbox->Add(databox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);

	wxBoxSizer *leftbox = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer *footbox = new wxBoxSizer(wxHORIZONTAL);
	leftbox->Add(parambox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	footbox->Add(storebox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	AddButton(ID_Compare, "Comp", 40, footbox);
	
	/*
	popfreq = NumPanel(50, wxALIGN_RIGHT, "0");
	ratemean = NumPanel(50, wxALIGN_RIGHT, "0");

	wxGridSizer *popdatagrid = new wxFlexGridSizer(2, 3, 3);
	popdatagrid->Add(TextLabel("Freq"), 0, wxALIGN_CENTRE);
	popdatagrid->Add(popfreq);
	popdatagrid->Add(TextLabel("SD"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	popdatagrid->Add(popsd);
	popdatagrid->Add(TextLabel("Resp"), 0, wxALIGN_CENTRE|wxST_NO_AUTORESIZE);
	popdatagrid->Add(ratemean);

	wxStaticBoxSizer *popdatabox = new wxStaticBoxSizer(wxVERTICAL, panel, "Population");
	popdatabox->AddSpacer(5);
	popdatabox->Add(popdatagrid, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 2);*/

	freqL1 = NumPanel(50, wxALIGN_RIGHT, "0");
	freqL2 = NumPanel(50, wxALIGN_RIGHT, "0");
	//freqL3 = NumPanel(50, wxALIGN_RIGHT, "0");

	wxGridSizer *laydatagrid = new wxFlexGridSizer(2, 3, 3);
	laydatagrid->Add(TextLabel("Freq L1"), 0, wxALIGN_CENTRE);
	laydatagrid->Add(freqL1);
	laydatagrid->Add(TextLabel("Freq L2"), 0, wxALIGN_CENTRE);
	laydatagrid->Add(freqL2);
	//laydatagrid->Add(TextLabel("Freq L3"), 0, wxALIGN_CENTRE);
	//laydatagrid->Add(freqL3);

	wxStaticBoxSizer *laydatabox = new wxStaticBoxSizer(wxVERTICAL, panel, "Layer Data");
	laydatabox->AddSpacer(5);
	laydatabox->Add(laydatagrid, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 2);

	rightbox->AddStretchSpacer(5);
	rightbox->Add(laydatabox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);
	rightbox->AddSpacer(5);

	leftbox->Add(footbox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);

	wxBoxSizer *hbox = new wxBoxSizer(wxHORIZONTAL);
	hbox->Add(leftbox, 0, wxALL, 0);
	hbox->AddStretchSpacer();
	hbox->Add(rightbox, 0, wxALIGN_CENTRE_VERTICAL|wxALL, 0);

	status = StatusBar();
	wxBoxSizer *statusbox = new wxBoxSizer(wxHORIZONTAL);
	statusbox->Add(status, 1, wxEXPAND);

	mainbox->AddSpacer(5);
	mainbox->Add(hbox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	mainbox->AddStretchSpacer(5);
	mainbox->Add(statusbox, 0, wxEXPAND);

	panel->Layout();

	Connect(ID_Run, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnRun));
	//Connect(ID_pararun, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnRun));
	Connect(ID_sum, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnSum));
	//Connect(ID_netstore, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnStore));
	//Connect(ID_netload, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnLoad));
	//Connect(ID_output, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnOutput));
	Connect(ID_paramstore, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnParamStore));
	Connect(ID_paramload, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnParamLoad));
	Connect(ID_Compare, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnParamLoad));
	Connect(ID_controls, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnControlsMenu));
	Connect(ID_model, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNetBox::OnModelMenu));

	Connect(wxEVT_COMMAND_TEXT_ENTER, wxCommandEventHandler(VMNNetBox::OnRun));
	Connect(wxEVT_SCROLL_LINEUP, wxSpinEventHandler(VMNNetBox::OnNeuroNext));
	Connect(wxEVT_SCROLL_LINEDOWN, wxSpinEventHandler(VMNNetBox::OnNeuroPrev));
	Connect(wxEVT_SCROLL_THUMBTRACK, wxSpinEventHandler(VMNNetBox::OnSpin));
}


// Not in use
wxToolBar* VMNNetBox::OnCreateToolBar(long style, wxWindowID id, const wxString &name)
{
	//return ParamBox::OnCreateToolBar(style, id, name);
	//return new NewToolBar(this, id, wxDefaultPosition, wxDefaultSize, style, name);
	return NULL;
}


void VMNNetBox::OnModelMenu(wxCommandEvent& event)
{
	int id = event.GetId();
	wxWindow *pos = FindWindowById(id, toolpanel);
	wxPoint point = pos->GetPosition();
	wxSize size = pos->GetSize();
	PopupMenu(menuModel, point.x, point.y + size.y);
}


void VMNNetBox::OnControlsMenu(wxCommandEvent& event)
{
	int id = event.GetId();
	wxWindow *pos = FindWindowById(id, toolpanel);
	wxPoint point = pos->GetPosition();
	wxSize size = pos->GetSize();
	PopupMenu(menuControls, point.x, point.y + size.y);
}


void VMNNetBox::OnSynQueue(wxCommandEvent& WXUNUSED(event))                 // Out of use, replaced by automatic modflag
{
	if((*modflags)["synqueue"]) {
		(*modflags)["synqueue"] = 0;
		menuModel->Check(ID_synqueue, 0);
	}
	else {	
		(*modflags)["synqueue"] = 1;
		menuModel->Check(ID_synqueue, 1);
	}
}


void VMNNetBox::ModData(VMNNeuron *data)
{

	if(data->netflag) snum = "sum";
	else snum = numstring(neuroindex, 0);
	datneuron->SetLabel(snum);
	snum.Printf("%d", data->spikecount);
	spikes->SetLabel(snum);
	snum.Printf("%.2f", data->freq);
	freq->SetLabel(snum);
	snum.Printf("%.1f", data->meanisi);
	mean->SetLabel(snum);
	snum.Printf("%.2f", data->isivar);
	sd->SetLabel(snum);

	if(!mod->basicmode) {
		snum.Printf("%.2f", data->esynL1);
		esyn->SetLabel(snum);
		snum.Printf("%.2f", data->vrest);
		vrest->SetLabel(snum);
		snum.Printf("%.1f", data->kHAP);
		kHAP->SetLabel(snum);
		snum.Printf("%.2f", data->tauHAP);
		tauHAP->SetLabel(snum);
		snum.Printf("%.2f", data->input);
		input->SetLabel(snum);
	}
}


void VMNNetBox::NetworkSum()
{
	netcalcvmn(mod->netdat1, 1);
	netcalcvmn(mod->netdat2, 2);
	netcalcvmn(mod->netdat3, 3);

	if(!mod->basicmode) {
		freqL1->SetLabel(text.Format("%.2f", mod->netdat1->freq));
		freqL2->SetLabel(text.Format("%.2f", mod->netdat2->freq));
		freqL3->SetLabel(text.Format("%.2f", mod->netdat3->freq));
	}

	freqL1->SetLabel(text.Format("%.2f", mod->netdat1->freq));
	freqL2->SetLabel(text.Format("%.2f", mod->netdat2->freq));
}


void VMNNetBox::OnSum(wxCommandEvent& WXUNUSED(event))
{
	wxString text; 

	if(sumflag) {
		OffSum();
		NeuroData();
	}
	else {
		sumflag = 1;
		ModData(&(neurons[vmhneurons]));	
	}	

	NetworkSum();
	mainwin->scalebox->GraphUpdate();
}


void VMNNetBox::OffSum()
{
	sumflag = 0;
	mainwin->scalebox->NetSwitch(0);
}


void VMNNetBox::NeuroData()
{
	mod->currvmn->neurocalc(&(neurons[neuroindex]));
	//mod->analysisdata->autocalc(mod->currvmh);
	mod->currvmn->id = neuroindex;
	mainwin->scalebox->inputflag = 0;
	mainwin->scalebox->vmhflag = 1;

	ParamStore *sigparams = mod->signalbox->GetParams();
	int timerange = (*sigparams)["timerange"];

	mod->currvmn->MeanSpikeForm(mod->neurodata->V, timerange, (*modflags)["formfilter"]);

	ModData(&(neurons[neuroindex]));
	mainwin->scalebox->GraphUpdate();
}


void VMNNetBox::OnNeuroPrev(wxSpinEvent& WXUNUSED(event))
{
	if(sumflag) OffSum();
	else {
		if(neuroindex > 0) neuroindex--;
		else neuroindex = vmhneurons-1;
	}
	NeuroData();
}


void VMNNetBox::OnNeuroNext(wxSpinEvent& WXUNUSED(event))
{
	if(sumflag) OffSum();
	else {
		if(neuroindex < vmhneurons-1) neuroindex++;
		else neuroindex = 0;
	}
	NeuroData();
}


void VMNNetBox::OnSpin(wxSpinEvent& event)
{
	if(autorun) OnRun(event);
}


void VMNNetBox::SetNeuroCount()
{
	GetParams();

	if(mod->basicmode) {
		(*modflags)["netgen"] = 1;
		(*modflags)["cellgen"] = 1;
		(*modflags)["unigen"] = 1;

		layer1 = (*modparams)["neuronsL1"];
		layer2 = (*modparams)["neuronsL2"];
		layer3 = 0;
	}
	else {
		(*modflags)["netgen"] = netcheck->GetValue();
		(*modflags)["cellgen"] = cellcheck->GetValue();
		(*modflags)["unigen"] = unicheck->GetValue();

		layer1 = (*modparams)["neuronsL1"];
		layer2 = (*modparams)["neuronsL2"];
		layer3 = (*modparams)["neuronsL3"];
	}

	vmhneurons = layer1 + layer2 + layer3;
}


void VMNNetBox::OnRun(wxCommandEvent& event)
{
	int id = event.GetId();
	int parallel = 0;
	long data;

	// Enter pressed for neuron selection
	if(id == ID_datneuro) {
		datneuron->GetValue().ToLong(&data);
		if(data >= 0 && data < vmhneurons) {
			neuroindex = data;
			NeuroData();
		}
		return;
	}

	if(id == ID_pararun) parallel = 1;
	countmark = 0;

	SetNeuroCount();

	(*mod->modeflags)["prototype"] = none;

	mod->RunModel();
}


void VMNNetBox::netcalcvmn(SpikeDat *spikedata, int layer)
{
	int i, step, stepmax;
	int vmnL1, vmnL2, vmnL3;
	int cellfrom, cellto, numcells;
	wxString text;

	double tB, kB, halflifeB, tauB;

	vmnL1 = (*modparams)["neuronsL1"];
	vmnL2 = (*modparams)["neuronsL2"];
	vmnL3 = (*modparams)["neuronsL3"];
	stepmax = (*mod->neurobox->modparams)["numspikes"];

	if(mod->basicmode) vmnL3 = 0;

	if(layer == 2) {
		cellfrom = vmnL1;
		cellto = vmnL1 + vmnL2;
		numcells = vmnL2;
	}
	if(layer == 1) {
		cellfrom = 0;
		cellto = vmnL1;
		numcells = vmnL1;
	}
	if(layer == 3) {
		cellfrom = vmnL1 + vmnL2;
		cellto = vmnL1 + vmnL2 + vmnL3;
		numcells = vmnL3;
	}

	neurons[vmhneurons].meanisi = 0;
	neurons[vmhneurons].count = 0;
	neurons[vmhneurons].netflag = 1;

	mainwin->diagbox->Write(text.Format("\nnetcalc  maxspikes %d\n", stepmax));


	/*ParamStore *sigparams = mod->signalbox->GetParams();
	halflifeB = (*sigparams)["halflifeB"];
	kB = (*sigparams)["kB"];

	tauB = log((double)2) / halflifeB;*/

	// fix stepmax

	for(i=0; i<=stepmax; i++) spikedata->srate[i] = 0;
	for(i=0; i<1000000; i++) spikedata->srate1[i] = 0;
	for(i=0; i<10000; i++) spikedata->hist1[i] = 0;
	for(i=0; i<10000; i++) spikedata->hist5[i] = 0;
	for(i=0; i<10000; i++) spikedata->haz1[i] = 0;
	for(i=0; i<10000; i++) spikedata->haz5[i] = 0;

	for(i=cellfrom; i<cellto; i++) {
		netdat->neurocalc(&(neurons[i]));
		neurons[vmhneurons].meanisi = neurons[vmhneurons].meanisi + neurons[i].meanisi/numcells;
		neurons[vmhneurons].count = neurons[vmhneurons].count + neurons[i].count;	
		for(step=0; step<=stepmax; step++) spikedata->srate[step] += netdat->srate[step];

		//tB = 0;
		for(step=0; step<1000000; step++) {
			spikedata->srate1[step] += netdat->srate1[step];
		}

		for(step=0; step<10000;step++) {
			spikedata->hist1[step] += netdat->hist1[step];
			spikedata->hist5[step] += netdat->hist5[step];
			spikedata->haz1[step] += netdat->haz1[step];
			spikedata->haz5[step] += netdat->haz5[step];
		}
	}
	if(neurons[vmhneurons].meanisi > 0) neurons[vmhneurons].freq = 1000/neurons[vmhneurons].meanisi;
	else neurons[vmhneurons].freq = 0;
	spikedata->freq = neurons[vmhneurons].freq;
	mod->SigSim(spikedata);
}


void VMNNetBox::netcalcvmh(SpikeDat *currvmh)
{
	int i, step, stepmax;
	int vmhL1;

	vmhL1 = (*modparams)["neuronsL1"];

	stepmax = (*modparams)["numspikes"];

	neurons[vmhneurons].meanisi = 0;
	neurons[vmhneurons].count = 0;
	neurons[vmhneurons].netflag = 1;

	for(i=0; i<=stepmax; i++) currvmh->srate.data[i] = 0;
	for(i=0; i<1000000; i++) currvmh->srate1.data[i] = 0;
	for(i=0; i<10000; i++) currvmh->hist1.data[i] = 0;
	for(i=0; i<10000; i++) currvmh->hist5.data[i] = 0;
	for(i=0; i<10000; i++) currvmh->haz1.data[i] = 0;
	for(i=0; i<10000; i++) currvmh->haz5.data[i] = 0;

	for(i=0; i<vmhL1; i++) {
		netdat->neurocalc(&(neurons[i]));
		neurons[vmhneurons].meanisi = neurons[vmhneurons].meanisi + neurons[i].meanisi/vmhneurons;
		neurons[vmhneurons].count = neurons[vmhneurons].count + neurons[i].count;	
		for(step=0; step<=stepmax; step++) currvmh->srate.data[step] = currvmh->srate.data[step] + netdat->srate.data[step];
		for(step=0; step<1000000; step++) currvmh->srate1.data[step] = currvmh->srate1.data[step] + netdat->srate1.data[step];
		for(step=0; step<10000;step++) {
			currvmh->hist1.data[step] = currvmh->hist1.data[step] + netdat->hist1.data[step];
			currvmh->hist5.data[step] = currvmh->hist5.data[step] + netdat->hist5.data[step];
			currvmh->haz1.data[step] = currvmh->haz1.data[step] + netdat->haz1.data[step];
			currvmh->haz5.data[step] = currvmh->haz5.data[step] + netdat->haz5.data[step];
		}
	}
	neurons[vmhneurons].freq = 1000/neurons[vmhneurons].meanisi;
	currvmh->inputsim((*modparams)["halflife"]);
}


void VMNNetBox::OnParamLoad(wxCommandEvent& event)
{
	wxString filetag, filename;

	if(synccheck->GetValue()) {
		filetag = paramstoretag->GetValue();
		mod->neurobox->paramstoretag->SetValue(filetag);
		mod->neurobox->OnParamLoad(event);
		mod->signalbox->paramstoretag->SetValue(filetag);
		mod->signalbox->OnParamLoad(event);
		mod->protobox->paramstoretag->SetValue(filetag);
		mod->protobox->OnParamLoad(event);
		mod->fitbox->paramstoretag->SetValue(filetag);
		mod->fitbox->OnParamLoad(event);
	}
	ParamBox::OnParamLoad(event);	
}


void VMNNetBox::OnParamStore(wxCommandEvent& event)
{
	wxString filetag, filename;

	if(synccheck->GetValue()) {
		filetag = paramstoretag->GetValue();
		mod->neurobox->paramstoretag->SetValue(filetag);
		mod->neurobox->OnParamStore(event);
		mod->signalbox->paramstoretag->SetValue(filetag);
		mod->signalbox->OnParamStore(event);
		mod->protobox->paramstoretag->SetValue(filetag);
		mod->protobox->OnParamStore(event);
		mod->fitbox->paramstoretag->SetValue(filetag);
		mod->fitbox->OnParamStore(event);
	}
	ParamBox::OnParamStore(event);	
}

// Store generated network cells and connections
void VMNNetBox::OnStore(wxCommandEvent& event)
{
	int i, c;
	wxString filetag, filename, filepath;
	wxString outline;

	filepath = mod->GetPath() + "/Params";

	// Cell data file
	filetag = storetag->GetValue();
	filename = filepath + "/" + filetag + "-cellgen.dat";

	wxTextFile cellfile(filename);
	if(!cellfile.Exists()) cellfile.Create();
	cellfile.Open();
	cellfile.Clear();
	outline.Printf("VMH Network: cell data");
	cellfile.AddLine(outline);
	cellfile.AddLine("");
	outline.Printf("%d neurons", vmhneurons);
	cellfile.AddLine(outline);
	cellfile.AddLine("");

	for(i=0; i<vmhneurons; i++) {
		outline.Printf("N%d: vrestsd %.4f khapsd %.4f tauhapsd %.4f", 
			i, neurons[i].vrestsdgen, neurons[i].kHAPsdgen, neurons[i].tauHAPsdgen);
		cellfile.AddLine(outline);
	}
	cellfile.Write();
	cellfile.Close();

	// Network data file
	filename = filepath + "/" + filetag + "-netgen.dat";
	wxTextFile netfile(filename);
	if(!netfile.Exists()) netfile.Create();
	netfile.Open();
	netfile.Clear();

	outline.Printf("VMH Network: network data");
	netfile.AddLine(outline);
	netfile.AddLine("");
	outline.Printf("%d neurons", vmhneurons);
	netfile.AddLine(outline);
	netfile.AddLine("");

	for(i=0; i<vmhneurons; i++) {
		outline.Printf("N%d: econnect %d Econ:", i, neurons[i].econnect);
		for(c = 0; c<neurons[i].econnect; c++) {
			snum.Printf("%d ", neurons[i].enetwork[c]);
			outline = outline + snum;
		}
		snum.Printf(" iconnect %d Icon:", neurons[i].iconnect);
		outline = outline + snum;
		for(c = 0; c<neurons[i].iconnect; c++) {
			snum.Printf("%d ", neurons[i].inetwork[c]);
			outline = outline + snum;
		}	
		netfile.AddLine(outline);
	}

	netfile.AddLine("");

	for(i=0; i<vmhneurons; i++) {
		outline.Printf("N%d: Eweight:", i);
		for(c = 0; c<neurons[i].econnect; c++) {
			snum.Printf("%.2f ", neurons[i].eweight[c]);
			outline = outline + snum;
		}
		snum.Printf(" Iweight:", neurons[i].iconnect);
		outline = outline + snum;
		for(c = 0; c<neurons[i].iconnect; c++) {
			snum.Printf("%.2f ", neurons[i].iweight[c]);
			outline = outline + snum;
		}	
		netfile.AddLine(outline);
	}
	netfile.Write();
	netfile.Close();
}


void VMNNetBox::OnLoad(wxCommandEvent& event)
{
	int i, c;
	long numcells;
	long econnect, iconnect;
	long econ, icon;
	long eweight, iweight;

	wxString filetag, filename, filepath;
	wxString readline, sdat, condat;
	wxString weightdat;

	filepath = mod->GetPath() + "/Params";

	// Cell data file
	filetag = storetag->GetValue();
	filename = filepath + "/" + filetag + "-cellgen.dat";

	wxTextFile cellfile(filename);

	if(!cellfile.Exists()) {
		storetag->SetValue("Not found");
		return;
	}
	cellfile.Open();
	readline = cellfile.GetLine(2);
	sdat = readline.BeforeFirst(' ');
	sdat.ToLong(&numcells);
	snum.Printf("%d loaded", numcells);
	storetag->SetValue(snum);

	for(i=0; i<numcells; i++) {
		readline = cellfile.GetLine(i+4);
		readline = readline.AfterFirst('d');
		sdat = readline.BeforeFirst('k');
		sdat.Trim();
		sdat.ToDouble(&neurons[i].vrestsdgen);
		readline = readline.AfterFirst('d');
		sdat = readline.BeforeFirst('t');
		sdat.Trim();
		sdat.ToDouble(&neurons[i].kHAPsdgen);
		readline = readline.AfterFirst('d');
		sdat = readline.Trim();
		sdat.ToDouble(&neurons[i].tauHAPsdgen);
		neurons[i].vrest = -62 + 10 * neurons[i].vrestsdgen;
	}

	snum.Printf("%d cells read", numcells);
	storetag->SetValue(snum);

	ModData(&(neurons[neuroindex]));

	// Net data file
	filename = filepath + "/" + filetag + "-netgen.dat";

	wxTextFile netfile(filename);

	if(!netfile.Exists()) {
		storetag->SetValue("Net not found");
		return;
	}
	netfile.Open();
	readline = netfile.GetLine(2);
	sdat = readline.BeforeFirst(' ');
	sdat.ToLong(&numcells);
	snum.Printf("net load %d cells", numcells);
	storetag->SetValue(snum);

	for(i=0; i<numcells; i++) {
		readline = netfile.GetLine(i+4);
		readline = readline.AfterFirst('t');
		sdat = readline.BeforeFirst('E');
		sdat.Trim();
		sdat.ToLong(&econnect);
		neurons[i].econnect = econnect;
		readline = readline.AfterFirst(':');
		condat = readline.BeforeFirst('i');

		for(c=0; c<econnect; c++) {
			sdat = condat.BeforeFirst(' ');
			sdat.ToLong(&econ);
			neurons[i].enetwork[c] = econ;
			condat = condat.AfterFirst(' ');
		}

		readline = readline.AfterFirst('t');
		sdat = readline.BeforeFirst('I');
		sdat.Trim();
		sdat.ToLong(&iconnect);
		neurons[i].iconnect = iconnect;
		readline = readline.AfterFirst(':');
		condat = readline;
		for(c=0; c<iconnect; c++) {
			sdat = condat.BeforeFirst(' ');
			sdat.ToLong(&icon);
			neurons[i].inetwork[c] = icon;
			condat = condat.AfterFirst(' ');
		}
	}	

	// Read synaptic weights

	for(i=0; i<numcells; i++) {
		readline = netfile.GetLine(i+numcells+5);
		readline = readline.AfterFirst(':');
		readline = readline.AfterFirst(':');
		weightdat = readline.BeforeFirst('I');

		//mainwin->diagbox->Write("weightdat " + weightdat + "\n");

		for(c=0; c<neurons[i].econnect; c++) {
			sdat = weightdat.BeforeFirst(' ');
			sdat.ToDouble(&neurons[i].eweight[c]);
			weightdat = weightdat.AfterFirst(' ');
		}
	
		readline = readline.AfterFirst(':');
		weightdat = readline;
		for(c=0; c<neurons[i].iconnect; c++) {
			sdat = weightdat.BeforeFirst(' ');
			sdat.ToDouble(&neurons[i].iweight[c]);
			weightdat = weightdat.AfterFirst(' ');
		}
	}	

	cellfile.Close();
	netfile.Close();	
}


void VMNNetBox::OnOutput(wxCommandEvent& event)
{
	wxString snum;
	wxString outdir = mainwin->outpath + "/output";

	mod->diagbox->Write(snum.Format("Data output to path: %s", outdir));

	snum.Printf("vmncell%d", neuroindex);
	if(sumflag) snum.Printf("vmncellnet");
	mod->currvmn->output(snum, outdir);

	mod->netdat1->output("L1", outdir);
	mod->netdat2->output("L2", outdir);
	mod->netdat3->output("L3", outdir);

	mod->diagbox->Write(" OK\n");
}