/*
*  vmnpanels.cpp
*  HypoModel
*
*  Created by Duncan MacGregor
*  University of Edinburgh 2016
*
*/

#include "vmnmodel.h"


// Input signal protocol box
VMNProtoBox::VMNProtoBox(VMNModel *model, const wxString& title, const wxPoint& pos, const wxSize& size)
	: ParamBox(model, title, pos, size, "PROTO", 1)
{
	int pnum, inpnum, rampnum;
	int pulsenum0, pulsenum1, rampnum0, rampnum1;
	int numwidth;
	boxname = "PROTO";
	mod = model;
	wxString tag;

	long notestyle = wxAUI_NB_TOP | wxAUI_NB_TAB_SPLIT | wxAUI_NB_TAB_MOVE | wxAUI_NB_SCROLL_BUTTONS;
	wxAuiNotebook *tabpanel = new wxAuiNotebook(this, wxID_ANY, wxDefaultPosition, wxDefaultSize, notestyle);

	pnum = 0;
	status = mod->netbox->status;


	// Range Panel

	ToolPanel *rangepanel = new ToolPanel(this, tabpanel);
	rangepanel->SetFont(boxfont);
	wxBoxSizer *rangesizer = new wxBoxSizer(wxVERTICAL);
	rangepanel->SetSizer(rangesizer);

	activepanel = rangepanel;
	paramset->panel = activepanel;

	labelwidth = 50;
	numwidth = 45;

	paramset->AddNum("rangestart", "Start", 200, 0, labelwidth, numwidth); 
	paramset->AddNum("rangestop", "Stop", 200, 0, labelwidth, numwidth); 
	paramset->AddNum("rangestep", "Step", 300, 0, labelwidth, numwidth); 

	wxStaticBoxSizer *rangebox0 = new wxStaticBoxSizer(wxVERTICAL, rangepanel, "Range Input L1");
	for(pnum=pnum; pnum<paramset->numparams; pnum++) {
		rangebox0->Add(paramset->con[pnum], 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxRIGHT|wxLEFT, 5);
	}

	wxBoxSizer *rangerunbox = new wxBoxSizer(wxHORIZONTAL);
	rangerunbox->Add(TextLabel("Input"), 1, wxALIGN_CENTRE);
	currentrange = NumPanel(40, wxALIGN_CENTRE, "---"); 
	rangerunbox->AddSpacer(10);
	rangerunbox->Add(currentrange, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);
	rangebox0->AddSpacer(10);
	rangebox0->Add(rangerunbox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxRIGHT|wxLEFT, 5);

	rangebox0->AddSpacer(10);
	AddButton(ID_Range, "Run", 50, rangebox0);

	wxBoxSizer *rangebox = new wxBoxSizer(wxHORIZONTAL);
	rangebox->Add(rangebox0, 0, wxALL, 5);
	rangesizer->AddSpacer(10);
	rangesizer->Add(rangebox, 1, wxALL, 0);
	rangepanel->Layout();


	// Pulse Panel 

	ToolPanel *pulsepanel = new ToolPanel(this, tabpanel);
	pulsepanel->SetFont(boxfont);
	wxBoxSizer *pulsesizer = new wxBoxSizer(wxVERTICAL);
	pulsepanel->SetSizer(pulsesizer);

	activepanel = pulsepanel;
	paramset->panel = activepanel;

	labelwidth = 50;
	numwidth = 45;

	paramset->AddNum("pulsebase0", "Base", 200, 0, labelwidth, numwidth); 
	paramset->AddNum("pulsestart0", "Start", 200, 0, labelwidth, numwidth); 
	paramset->AddNum("pulsestop0", "Stop", 300, 0, labelwidth, numwidth); 
	paramset->AddNum("pulseinit0", "Initial", 200, 0, labelwidth, numwidth); 
	paramset->AddNum("pulsehl0", "Halflife", 0.1, 2, labelwidth, numwidth); 
	//paramset->SetMinMax("rampstep0", -1000, 1000);
	pulsenum0 = paramset->numparams;

	paramset->GetCon("pulseinit0")->SetMinMax(-10000, 10000);

	wxStaticBoxSizer *pulsebox0 = new wxStaticBoxSizer(wxVERTICAL, pulsepanel, "Pulse Input L1");
	for(pnum=pnum; pnum<paramset->numparams; pnum++) {
		pulsebox0->Add(paramset->con[pnum], 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxRIGHT|wxLEFT, 5);
	}
	wxBoxSizer *pulserunbox = new wxBoxSizer(wxHORIZONTAL);

	pulserunbox->Add(TextLabel("Input"), 1, wxALIGN_CENTRE);
	currentpulse = NumPanel(40, wxALIGN_CENTRE, "---"); 
	pulserunbox->AddSpacer(10);
	pulserunbox->Add(currentpulse, 1, wxALIGN_CENTRE_VERTICAL);
	pulsebox0->AddSpacer(10);
	pulsebox0->Add(pulserunbox, 1, wxRIGHT|wxLEFT, 5);
	pulsebox0->AddSpacer(10);
	AddButton(ID_Pulse, "Run", 50, pulsebox0);

	paramset->AddNum("pulsebase1", "Base", 200, 0, labelwidth, numwidth); 
	paramset->AddNum("pulsestart1", "Start", 200, 0, labelwidth, numwidth); 
	paramset->AddNum("pulsestop1", "Stop", 300, 0, labelwidth, numwidth); 
	paramset->AddNum("pulseinit1", "Initial", 200, 0, labelwidth, numwidth); 
	paramset->AddNum("pulsehl1", "Halflife", 0.1, 2, labelwidth, numwidth); 

	paramset->GetCon("pulseinit1")->SetMinMax(-10000, 10000);
	//paramset->SetMinMax("rampstep0", -1000, 1000);

	wxStaticBoxSizer *pulsebox1 = new wxStaticBoxSizer(wxVERTICAL, pulsepanel, "Pulse Input L2");
	for(pnum=pnum; pnum<paramset->numparams; pnum++) {
		pulsebox1->Add(paramset->con[pnum], 1, wxALIGN_CENTRE_HORIZONTAL|wxRIGHT|wxLEFT, 5);
	}

	pulsebox1->AddSpacer(10);
	AddButton(ID_Pulse, "Run", 50, pulsebox1);


	wxBoxSizer *pulsebox = new wxBoxSizer(wxHORIZONTAL);
	pulsebox->Add(pulsebox0, 0, wxALL, 5);
	pulsebox->AddStretchSpacer();
	pulsebox->Add(pulsebox1, 0, wxALL, 5);
	pulsesizer->AddSpacer(10);
	pulsesizer->Add(pulsebox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	pulsepanel->Layout();


	// Ramp Panel
	ToolPanel *ramppanel = new ToolPanel(this, tabpanel);
	ramppanel->SetFont(boxfont);
	wxBoxSizer *rampsizer = new wxBoxSizer(wxVERTICAL);
	ramppanel->SetSizer(rampsizer);

	activepanel = ramppanel;
	paramset->panel = activepanel;

	labelwidth = 50;
	numwidth = 45;

	paramset->AddNum("rampbase0", "Base", 200, 0, labelwidth, numwidth); 
	paramset->AddNum("rampstart0", "Start", 200, 0, labelwidth, numwidth); 
	paramset->AddNum("rampstop0", "Stop", 300, 0, labelwidth, numwidth); 
	paramset->AddNum("rampinit0", "Initial", 200, 0, labelwidth, numwidth); 
	paramset->AddNum("rampstep0", "1s Step", 0.1, 2, labelwidth, numwidth); 
	paramset->SetMinMax("rampstep0", -1000, 1000);
	rampnum0 = paramset->numparams;

	wxStaticBoxSizer *rampbox0 = new wxStaticBoxSizer(wxVERTICAL, ramppanel, "Input Ramp L1");
	for(pnum=pnum; pnum<paramset->numparams; pnum++) {
		rampbox0->Add(paramset->con[pnum], 1, wxALIGN_CENTRE_HORIZONTAL|wxRIGHT|wxLEFT, 5);
	}
	wxBoxSizer *inputrunbox = new wxBoxSizer(wxHORIZONTAL);

	inputrunbox->Add(TextLabel("Input"), 1, wxALIGN_CENTRE);
	currentinput = NumPanel(40, wxALIGN_CENTRE, "---"); 
	inputrunbox->AddSpacer(10);
	inputrunbox->Add(currentinput, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL);
	rampbox0->AddSpacer(10);
	rampbox0->Add(inputrunbox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxRIGHT|wxLEFT, 5);
	rampbox0->AddSpacer(10);
	AddButton(ID_Ramp, "Run", 50, rampbox0);

	paramset->AddNum("rampbase1", "Base", 200, 0, labelwidth, numwidth); 
	paramset->AddNum("rampstart1", "Start", 200, 0, labelwidth, numwidth); 
	paramset->AddNum("rampstop1", "Stop", 300, 0, labelwidth, numwidth); 
	paramset->AddNum("rampinit1", "Initial", 200, 0, labelwidth, numwidth); 
	paramset->AddNum("rampstep1", "1s Step", 0.1, 2, labelwidth, numwidth); 

	paramset->SetMinMax("rampstep1", -1000, 1000);
	rampnum1 = paramset->numparams;

	wxStaticBoxSizer *rampbox1 = new wxStaticBoxSizer(wxVERTICAL, ramppanel, "Input Ramp L2");
	for(pnum=pnum; pnum<paramset->numparams; pnum++) {
		rampbox1->Add(paramset->con[pnum], 1, wxALIGN_CENTRE_HORIZONTAL|wxRIGHT|wxLEFT, 5);
	}
	rampbox1->AddSpacer(10);
	AddButton(ID_Ramp, "Run", 50, rampbox1);

	wxBoxSizer *rampbox = new wxBoxSizer(wxHORIZONTAL);
	rampbox->Add(rampbox0, 0, wxALL, 5);
	rampbox->AddStretchSpacer();
	rampbox->Add(rampbox1, 0, wxALL, 5);
	rampsizer->AddSpacer(10);
	rampsizer->Add(rampbox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALL, 0);
	ramppanel->Layout();


	//////////////////////////////////////////////////
	// Main Structure

	tabpanel->Freeze();
	tabpanel->AddPage(ramppanel, "Ramp" , false);
	tabpanel->AddPage(pulsepanel, "Pulse" , true);
	tabpanel->AddPage(rangepanel, "Range" , false);
	tabpanel->Thaw();

	ToolPanel *storepanel = new ToolPanel(this, wxDefaultPosition, wxDefaultSize);
	wxBoxSizer *storesizer = new wxBoxSizer(wxVERTICAL);
	storepanel->SetSizer(storesizer);

	activepanel = storepanel;
	wxBoxSizer *paramfilebox = StoreBox("ns1", storepanel);

	storesizer->Add(paramfilebox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);	
	storesizer->Layout();

	winman->AddPane(tabpanel, wxAuiPaneInfo().Name("tabpane").CentrePane().PaneBorder(false));
	winman->AddPane(storepanel, wxAuiPaneInfo().Name("storepane").Bottom().CaptionVisible(false).BestSize(-1, 30).PaneBorder(false));
	winman->Update();

	Connect(ID_Ramp, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNProtoBox::OnRun));
	Connect(ID_Pulse, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNProtoBox::OnRun));
	Connect(ID_Range, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNProtoBox::OnRun));
}


void VMNProtoBox::OnRun(wxCommandEvent& event)
{
	if(event.GetId() == ID_Ramp) (*(mod->modeflags))["prototype"] = ramp;
	if(event.GetId() == ID_Pulse) (*(mod->modeflags))["prototype"] = pulse;
	if(event.GetId() == ID_Range) (*(mod->modeflags))["prototype"] = range;

	mod->netbox->SetNeuroCount();

	mod->RunModel();
} 

// Control box for testing signal processing
SignalBox::SignalBox(VMNModel *vmnmodel, const wxString& title, const wxPoint& pos, const wxSize& size)
	: ParamBox(vmnmodel, title, pos, size)
{
	column = 0;
	boxname = "Signal";
	mod = vmnmodel;

	InitMenu();

	//SetModFlag(ID_revpots, "revpots", "Reversal Potentials", 0); 
	SetModFlag(ID_noise, "noiseflag", "Noise Signal", 0); 
	SetModFlag(ID_L1L2, "layerflag", "Layer L2", 0); 


	// Parameter controls
	//
	// AddCon(tag string, display string, initial value, click increment, decimal places)
	// ----------------------------------------------------------------------------------

	paramset->AddCon("vrest", "V Rest", -62, 0.1, 2);
	paramset->AddCon("pspheight", "PSP mag", 4, 0.1, 2);
	paramset->AddCon("psphalflife", "PSP halflife", 7.5, 0.01, 3);

	paramset->AddCon("noimean", "Noise Mean", 300, 1, 2);
	paramset->AddCon("noitau", "Noise Tau", 1000, 1, 2);
	paramset->AddCon("noiamp", "Noise Amp", 1, 0.1, 2);
	
	paramset->AddCon("sigIratio", "Iratio", 0, 0.1, 2);

	if(!mod->basicmode) {
		paramset->AddCon("synwaveamp", "Wave Amp", 0, 1, 2);
		paramset->AddCon("synwavecycle", "Wave Cycle", 1000, 1, 2);
		paramset->AddCon("synwaveshift", "Wave Shift", 0, 0.1, 2);

		paramset->AddCon("kB", "kB", 0, 1, 4);
		paramset->AddCon("halflifeB", "halflifeB", 100, 1, 2);
		paramset->AddCon("timerange", "timerange", 100, 1, 0);
	}

	ParamLayout(2);

	// ----------------------------------------------------------------------------------

	defbutt = 0;
	wxBoxSizer *runbox = RunBox();

	wxSizer *paramfilebox = StoreBox("test1");

	mainbox->AddSpacer(5);
	mainbox->Add(parambox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	mainbox->AddStretchSpacer();
	mainbox->Add(runbox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);	
	mainbox->AddSpacer(5);
	mainbox->Add(paramfilebox, 0, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);	
	mainbox->AddSpacer(2);

	panel->Layout();
}


void SignalBox::OnRun(wxCommandEvent& event)
{
	mod->SigSim(mod->netdat1);
	mod->SigSim(mod->netdat2);

	ParamStore *sigparams = GetParams();
	int timerange = (*sigparams)["timerange"];
	if((*modflags)["layerflag"] == 0) mod->currvmn->MeanSpikeForm(mod->neurodata->V, timerange, (*mod->netbox->modflags)["formfilter"]);
	if((*modflags)["layerflag"] == 1) mod->currvmn->MeanSpikeForm(mod->neurodata->V2, timerange, (*mod->netbox->modflags)["formfilter"]);

	mainwin->scalebox->GraphUpdate();
}


// Individual neuron model parameter control box
VMNNeuroBox::VMNNeuroBox(VMNModel *vmnmodel, const wxString& title, const wxPoint& pos, const wxSize& size)
	: ParamBox(vmnmodel, title, pos, size)
{
	column = 0;
	mod = vmnmodel;
	boxname = "VMNNeuro";
	buttonwidth = 50;

	InitMenu();
	
	if(mod->basicmode) PanelBasic();
	else PanelFull();
}


void VMNNeuroBox::PanelFull()
{
	SetModFlag(ID_revpots, "revpots", "Reversal Potentials", 1); 
	SetModFlag(ID_DAPcap, "DAPcapflag", "DAP cap", 0); 
	SetModFlag(ID_DAP2, "DAP2flag", "Old DAP2", 0); 
	SetModFlag(ID_vsyn, "vsynflag", "V Syn", 1); 
	SetModFlag(ID_runtime, "runtimeflag", "Runtime", 1); 
	SetModFlag(ID_DAP, "DAPflag", "IKleak DAP", 0); 
	SetModFlag(ID_vDAP, "vDAPflag", "IKleak vDAP", 0); 
	SetModFlag(ID_Dyno, "Dynoflag", "IKleak Dyno", 0); 
	SetModFlag(ID_vIKleak, "vIKleakflag", "IKleak V", 0); 
	SetModFlag(ID_Iratio, "Iratioflag", "Use Iratio", 1); 
	SetModFlag(ID_artspikes, "artspikesflag", "Artificial Spikes", 0); 

	paramset->AddCon("numspikes", "Num Spikes", 1000, 1, 0);
	paramset->AddCon("hstep", "h Step", 1, 0.1, 1);
	paramset->AddCon("vthre", "V Threshold", -50, 0.1, 2);
	paramset->AddCon("vrest", "V Rest", -62, 0.1, 2);
	paramset->AddCon("absref", "Abs Ref", 5, 1, 1);
	paramset->AddCon("pspmag", "PSP mag", 4, 0.1, 2);
	paramset->AddCon("ve", "EPSP ve", 20, 1, 2);
	paramset->AddCon("vi", "IPSP vi", -110, 1, 2);
	paramset->AddCon("iratio", "IPSP ratio", 1, 0.1, 2);
	paramset->AddCon("ire", "EPSP freq", 300, 10, 1);
	paramset->AddCon("halflife", "half-life", 7.5, 0.1, 2);
	paramset->AddCon("emax", "EPSP max", 0, 0.1, 2);

	paramset->AddCon("kHAP", "HAP k", 60, 1, 2);
	paramset->AddCon("halflifeHAP", "HAP HL", 7, 0.1, 2);
	paramset->AddCon("kAHP", "AHP k", 0.0, 1, 2);
	paramset->AddCon("halflifeAHP", "AHP HL", 350, 1, 2);
	paramset->AddCon("kDAP", "DAP k", 0, 1, 2);
	paramset->AddCon("halflifeDAP", "DAP HL", 150, 1, 2);
	paramset->AddCon("iri", "IPSP freq", 300, 10, 1);

	SetVBox(2);
	if(!column) column = paramset->numparams / 2 + paramset->numparams % 2;

	for(i=0; i<column; i++) {
		vbox[0]->Add(paramset->con[i], 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxRIGHT|wxLEFT, 5);
		vbox[0]->AddSpacer(5);
	}

	for(i=column; i<paramset->numparams; i++) {
		vbox[1]->Add(paramset->con[i], 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxRIGHT|wxLEFT, 5);
		vbox[1]->AddSpacer(5);
	}

	parambox->Add(vbox[0], 0);
	parambox->Add(vbox[1], 0);

	runbutton = new wxButton(panel, ID_Run, "RUN", wxDefaultPosition, wxSize(70, buttonheight));  //50
	runbutton->SetFont(confont);
	resetbutton = new wxButton(panel, ID_Reset, "RESET", wxDefaultPosition, wxSize(70, buttonheight));  //50
	resetbutton->SetFont(confont);
	wxBoxSizer *runbox = new wxBoxSizer(wxHORIZONTAL);
	runbox->Add(runbutton);
	runbox->AddSpacer(20);
	runbox->Add(resetbutton);

	wxSizer *paramfilebox = StoreBox("n0fitb8");

	buttonbox = new wxBoxSizer(wxHORIZONTAL);

	AddButton(ID_Signal, "SIGNAL", buttonwidth, buttonbox);
	buttonbox->AddSpacer(5);
	buttonbox->AddStretchSpacer();

	AddButton(ID_Protocol, "PROTO", buttonwidth, buttonbox);
	buttonbox->AddSpacer(5);
	buttonbox->AddStretchSpacer();

	AddButton(ID_Output, "OUTPUT", 55, buttonbox);
	SetPanel(ID_Output, mod->outbox); 

	AddButton(ID_EvoFit, "FIT", 55, buttonbox);
	SetPanel(ID_EvoFit, mod->fitbox); 

	mainbox->AddSpacer(5);
	mainbox->Add(parambox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	mainbox->AddStretchSpacer(5);
	mainbox->Add(runbox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);	
	mainbox->AddStretchSpacer(5);
	mainbox->Add(paramfilebox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);	
	mainbox->AddStretchSpacer(5);
	mainbox->Add(buttonbox, 1, wxALIGN_CENTRE_HORIZONTAL | wxALIGN_CENTRE_VERTICAL | wxALL, 5);
	mainbox->AddSpacer(2);

	panel->Layout();

	Connect(ID_paramstore, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNeuroBox::OnParamStore));
	Connect(ID_paramload, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNeuroBox::OnParamLoad));
	Connect(ID_Signal, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNeuroBox::OnBox));
	Connect(ID_Protocol, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNeuroBox::OnBox));
	//Connect(ID_EvoFit, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNeuroBox::OnBox));
}


void VMNNeuroBox::PanelBasic()
{
	//SetModFlag(ID_revpots, "revpots", "Reversal Potentials", 1); 
	//SetModFlag(ID_DAPcap, "DAPcapflag", "DAP cap", 0); 
	//SetModFlag(ID_DAP2, "DAP2flag", "Old DAP2", 0); 
	SetModFlag(ID_vsyn, "vsynflag", "V Syn", 1); 
	SetModFlag(ID_runtime, "runtimeflag", "Runtime", 1); 
	//SetModFlag(ID_DAP, "DAPflag", "IKleak DAP", 0); 
	//SetModFlag(ID_vDAP, "vDAPflag", "IKleak vDAP", 0); 
	//SetModFlag(ID_Dyno, "Dynoflag", "IKleak Dyno", 0); 
	//SetModFlag(ID_vIKleak, "vIKleakflag", "IKleak V", 0); 
	SetModFlag(ID_Iratio, "Iratioflag", "Use Iratio", 1); 
	//SetModFlag(ID_artspikes, "artspikesflag", "Artificial Spikes", 0); 

	paramset->AddCon("numspikes", "Num Spikes", 1000, 1, 0);
	paramset->AddCon("hstep", "h Step", 1, 0.1, 1);
	paramset->AddCon("vthre", "V Threshold", -50, 0.1, 2);
	paramset->AddCon("vrest", "V Rest", -62, 0.1, 2);
	paramset->AddCon("absref", "Abs Ref", 2, 1, 1);
	paramset->AddCon("pspmag", "PSP mag", 4, 0.1, 2);
	//paramset->AddCon("ve", "EPSP ve", 20, 1, 2);
	//paramset->AddCon("vi", "IPSP vi", -110, 1, 2);
	paramset->AddCon("iratio", "IPSP ratio", 1, 0.1, 2);
	//paramset->AddCon("ire", "EPSP freq", 300, 10, 1);
	paramset->AddCon("halflife", "PSP HL", 7.5, 0.1, 2);
	//paramset->AddCon("emax", "EPSP max", 0, 0.1, 2);

	paramset->AddCon("kHAP", "HAP k", 60, 1, 2);
	paramset->AddCon("halflifeHAP", "HAP HL", 7, 0.1, 2);
	paramset->AddCon("kAHP", "AHP k", 0.0, 1, 2);
	paramset->AddCon("halflifeAHP", "AHP HL", 350, 1, 2);
	paramset->AddCon("kDAP", "DAP k", 0, 1, 2);
	paramset->AddCon("halflifeDAP", "DAP HL", 150, 1, 2);
	//paramset->AddCon("iri", "IPSP freq", 300, 10, 1);

	ParamLayout(2);

	runbutton = new wxButton(panel, ID_Run, "RUN", wxDefaultPosition, wxSize(70, buttonheight));  //50
	runbutton->SetFont(confont);
	resetbutton = new wxButton(panel, ID_Reset, "RESET", wxDefaultPosition, wxSize(70, buttonheight));  //50
	resetbutton->SetFont(confont);
	wxBoxSizer *runbox = new wxBoxSizer(wxHORIZONTAL);
	runbox->Add(runbutton);
	runbox->AddSpacer(20);
	runbox->Add(resetbutton);

	wxSizer *paramfilebox = StoreBox("n0fitb8");

	buttonbox = new wxBoxSizer(wxHORIZONTAL);

	AddButton(ID_Signal, "SIGNAL", buttonwidth, buttonbox);
	buttonbox->AddSpacer(5);
	buttonbox->AddStretchSpacer();

	AddButton(ID_Protocol, "PROTO", buttonwidth, buttonbox);
	buttonbox->AddSpacer(5);
	buttonbox->AddStretchSpacer();

	AddButton(ID_Output, "OUTPUT", 55, buttonbox);
	SetPanel(ID_Output, mod->outbox); 
	buttonbox->AddSpacer(5);
	buttonbox->AddStretchSpacer();

	AddButton(ID_EvoFit, "FIT", 55, buttonbox);

	mainbox->AddSpacer(5);
	mainbox->Add(parambox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);
	mainbox->AddStretchSpacer(5);
	mainbox->Add(runbox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 5);	
	mainbox->AddStretchSpacer(5);
	mainbox->Add(paramfilebox, 1, wxALIGN_CENTRE_HORIZONTAL|wxALIGN_CENTRE_VERTICAL|wxALL, 0);	
	mainbox->AddStretchSpacer(5);
	mainbox->Add(buttonbox, 1, wxALIGN_CENTRE_HORIZONTAL | wxALIGN_CENTRE_VERTICAL | wxALL, 5);
	mainbox->AddSpacer(2);

	panel->Layout();

	Connect(ID_paramstore, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNeuroBox::OnParamStore));
	Connect(ID_paramload, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNeuroBox::OnParamLoad));
	Connect(ID_Signal, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNeuroBox::OnBox));
	Connect(ID_Protocol, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(VMNNeuroBox::OnBox));
}


void VMNNeuroBox::OnBox(wxCommandEvent& event)
{
	int id = event.GetId();
	ToolBox *box;

	if(id == ID_Signal) box = mod->signalbox;
	if(id == ID_Protocol) box = mod->protobox;

	if(box->IsShown()) box->Show(false);
	else box->Show(true);
}
