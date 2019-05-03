#include "mapping.h"
#include <TH1.h>
#include <TF1.h>
#include <TH2D.h>
#include <TFitResultPtr.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <fstream>
#include <map>
#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <vector>
#include "TVirtualFFT.h"
using namespace std;

 Double_t xpar[128][2];
 Double_t ypar[48][2];
 
 Double_t xpar_xia[128][2];
 Double_t ypar_xia[48][2];
 Double_t sl=2.56038,cept=-33.8167;

void mapping::Readpar(){
  ifstream in("EF_equal.dat");
  if(!in.is_open()){
	 cout<<"cannot open file !"<<endl;
	 for(Int_t i=0;i<128;i++){
		xpar[i][0]=1;
		xpar[i][0]=0;
	}
	 for(Int_t i=0;i<48;i++){
		ypar[i][0]=1;
		ypar[i][0]=0;
	}
	 return;
  }
  Double_t xp,yp;
  Int_t ch;
  while(1){
    in>>ch>>xp>>yp;
	if(!in.good()) break;
	if(ch<128) {
	  xpar[ch][0]=xp;
	  xpar[ch][1]=yp;
	} else{
	  ypar[ch-128][0]=xp;
	  ypar[ch-128][1]=yp;
	}
  }
}

void mapping::Readpar_xia(){
  TFile *file = new TFile("map00315xia.root");
  fChain = (TTree*)file->Get("tree");
  Double_t tpar[2];
  memset(tpar,0,sizeof(tpar));
  TGraph *grb[128];
  TGraph *grf[48];
  Int_t cb[128];
  Int_t cf[48];
  memset(cb,0,sizeof(cb));
  memset(cf,0,sizeof(cf));
  TF1 *f = new TF1("f","pol1",0,140000);
  for(Int_t i=0;i<128;i++)
  {
  if(i<48)grf[i] = new TGraph;
  grb[i] = new TGraph;
  }
  
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("er",&er,&b_er);
   fChain->SetBranchAddress("rte",&rte,&b_rte);
   fChain->SetBranchAddress("np",&np,&b_np);
   fChain->SetBranchAddress("det",&det,&b_det);
   fChain->SetBranchAddress("str",&str,&b_str);
   fChain->SetBranchAddress("fr",&fr,&b_fr);

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
   if(fr==2&&np==1&&det==2)
   {
   grb[str]->SetPoint(cb[str],rte[0],er);
   cb[str]++;
   }
   
   if(fr==1&&np==1&&det==2)
   {
   grf[str]->SetPoint(cf[str],rte[0],er);
   cf[str]++;
   }
  // cout<<er<<"  "<<rte[0]<<"  "<<det<<"  "<<str<<"  "<<fr<<endl;
}
  Int_t ch=0;
  ofstream out("cali_xia.dat");

  for(Int_t i=0;i<128;i++)
  {
  grb[i]->Fit(f,"1Q");
  f->GetParameters(&tpar[0]);
  if(tpar[1]<0){tpar[1]=1;tpar[0]=0;}
  if(tpar[1]==0&&tpar[0]==0)continue;
  out<<ch<<"  "<<tpar[1]<<"  "<<tpar[0]<<endl;
  delete grb[i];
  ch++;
  }
  for(Int_t i=0;i<48;i++)
  {
  grf[i]->Fit(f,"1Q");
  f->GetParameters(&tpar[0]);
  if(tpar[1]<0){tpar[1]=1;tpar[0]=0;}
  if(tpar[1]==0&&tpar[0]==0)continue;
  out<<ch<<"  "<<tpar[1]<<"  "<<tpar[0]<<endl;
  delete grf[i];
  ch++;
  }
  delete f;
  //read parameters
  ifstream in("cali_xia.dat");
  if(!in.is_open()){
	 cout<<"cannot open file !"<<endl;
	 for(Int_t i=0;i<128;i++){
		xpar_xia[i][0]=1;
		xpar_xia[i][0]=0;
	}
	 for(Int_t i=0;i<48;i++){
		ypar_xia[i][0]=1;
		ypar_xia[i][0]=0;
	}
	 return;
  }
  Double_t xp,yp;
  while(1){
    in>>ch>>xp>>yp;
	if(!in.good()) break;
	if(ch<128) {
	  xpar_xia[ch][0]=xp;
	  xpar_xia[ch][1]=yp;
	} else{
	  ypar_xia[ch-128][0]=xp;
	  ypar_xia[ch-128][1]=yp;
	}
  }
  
}

 TH1D *hx[128],*hnx,*htx,*hehx;
 TH1D *hy[48],*hny,*hty,*hehy;
 Int_t nx[128],ny[48];

void mapping::SetHist()
{
  TString sx,sy;
  Int_t nx[128],ny[48];
  for(Int_t i=0;i<128;i++) {
    sx.Form("hx%03i",i);
    hx[i]=new TH1D(sx.Data(),sx.Data(),740,0,740);
    nx[i]=0;
    if(i<48) {
      sy.Form("hy%03i",i);
      hy[i]=new TH1D(sy.Data(),sy.Data(),740,0,740);
      ny[i]=0;
    } 
  }
    hnx=new TH1D("hnx","hnx",128,0,128);
    hny=new TH1D("hny","hny",48,0,48);
    htx=new TH1D("htx","htx",128,0,128);
    hty=new TH1D("hty","hty",48,0,48);
    hehx=new TH1D("hehx","hehx",128,0,128);
    hehy=new TH1D("hehy","hehy",48,0,48);
}

TF1 *f1;
TF1 *f2;
TF1 *f3;
TGraph *gr = new TGraph;
TCanvas *c1 = new TCanvas;
TGraph *gra = new TGraph;
void mapping::FillHist()
{
  if(det!=2) return;
  if(tsp[0]<125||tsp[0]>137) return;
  if(np!=1) return;
/*  delete c1;
  c1=new TCanvas;
  delete gra;
  gra=new TGraph;
       for(int i=0;i<750;i++)
        gra->SetPoint(i,i,dwave[i]);
        
        gra->Draw();
        c1->SaveAs("c1.eps");
        sleep(5);*/

  if(fr==2) {
    Int_t x=str;
    if(nx[x]<1000 && e>4000 && e<10000) {//choose a peak.5200-5140//4800-5800
    //align for tsp[0]==132, superpulse range from 0-740;
     Int_t tx0 = tsp[0];
     for(Int_t i=0;i<50;i++)
       hx[x]->Fill(i,0);
     
     for(Int_t i=132;i>0;i--){
       hx[x]->Fill(i,dwave[tx0]*2000./er); //since the tx0
       tx0--;
       if(tx0<1) break;
       }
     
     tx0 = tsp[0]+1;
     for(Int_t i=133;i<740;i++){
       hx[x]->Fill(i,dwave[tx0]*2000./er);
       tx0++;
       if(tx0>749) break;
       }
       nx[x]++;
   }
  }
  
  if(fr==1) {
    Int_t y=str;
    if(ny[y]<1000 && e>4000 && e<10000) {
     Int_t ty0 = tsp[0];
     for(Int_t i=0;i<50;i++)
       hy[y]->Fill(i,0);
     
     for(Int_t i=132;i>0;i--){
       hy[y]->Fill(i,dwave[ty0]*2000./er);
       ty0--;
       if(ty0<1)break;
       }
     
     ty0 = tsp[0]+1;//repeat fill increase the bin value
     for(Int_t i=133;i<740;i++) {
       hy[y]->Fill(i,dwave[ty0]*2000./er);
       ty0++;
       if(ty0>749) break;
       }
       ny[y]++;
    }
  }
}

void mapping::WriteHist(Int_t runnum)
{
  char opfh[124];
  number=runnum;
  sprintf(opfh,"./hist%05d.root",number);
  TFile *opf1=new TFile(opfh,"RECREATE");
  for(Int_t i=0;i<128;i++) {
    if(nx[i]<1)continue;
    hx[i]->Scale(1./nx[i]);
    hx[i]->Write();
    hnx->Fill(i,nx[i]);
  //  cout<<"i,nx[i]:"<<i<<" "<<nx[i]<<endl;
    for(Int_t j=0;j<750;j++) dwave[j]=hx[i]->GetBinContent(j+1);
    ResetOpt();
    det=2;
    fr=2;
    str=i;
    PuDet();
    htx->Fill(i,tsp[0]);
    hehx->Fill(i,eh[0]);
  //  if(np==1) cout<<fr<<" x "<<str<<" "<<np<<" "<<tsp[0]<<" "<<eh[0]<<endl;
   // if(np==2) cout<<fr<<" x "<<str<<" "<<np<<" "<<tsp[0]<<" "<<eh[0]<<" "<<tsp[1]<<" "<<eh[2]<<endl;  
  }

  for(Int_t i=0;i<48;i++) {
   if(ny[i]<1)continue;
   hy[i]->Scale(1./ny[i]);
   hy[i]->Write();
   hny->Fill(i,ny[i]);
    for(Int_t j=0;j<750;j++) dwave[j]=hy[i]->GetBinContent(j+1);
    ResetOpt();
    det=2;
    fr=1;
    str=i;
    PuDet();
    hty->Fill(i,tsp[0]);
    hehy->Fill(i,eh[0]);
  //  if(np==1) cout<<fr<<" y "<<str<<" "<<np<<" "<<tsp[0]<<" "<<eh[0]<<endl;
 //   if(np==2) cout<<fr<<" y "<<str<<" "<<np<<" "<<tsp[0]<<" "<<eh[0]<<" "<<tsp[1]<<" "<<eh[1]<<endl;  
  // cout<<"i,ny[i]:"<<i<<" "<<ny[i]<<endl;
   }
  hnx->Write();
  hny->Write();
  htx->Write();
  hty->Write();
  hehx->Write();
  hehy->Write();
  opf1->Close();
}

   TRandom3 *ru;
void mapping::Loop(TTree *ipt,Int_t bdn,TTree *opt)
{
   Init(ipt);
   boardn=bdn;
   ru=new TRandom3();

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   Int_t be = 2705;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //if(Data.TriggerTimeTag!=3981594868210)continue;
      // if (Cut(ientry) < 0) continue;
      if(det==4)continue;
	  ResetOpt();
	  if(Mapping()){
	    opt->Fill();
            //if(jentry*1000/nentries%50==0) cout<<" processing "<<Int_t(jentry*1000./nentries)/10.<<"%"<<endl;
#ifdef HIST
            FillHist();//fill histogram
#endif
	    nevt++;
	   // if(jentry%1000==0)cout<<"Processing..."<<"  "<<jentry<<"  Total events:"<<nentries<<endl;
	   }
   }
}

mapping::mapping(TTree *opt) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if(opt==0) {
    cout<<"no output tree!"<<endl;
	exit(1);
  }
  BranchOpt(opt);
  nevt=0;
}

mapping::~mapping()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mapping::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t mapping::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void mapping::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // poInt_ters of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch poInt_ters
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("BankNo", &BankNo, &b_BankNo);
   fChain->SetBranchAddress("ChannelNo", &ChannelNo, &b_ChannelNo);
   fChain->SetBranchAddress("SysTime", &SysTime, &b_SysTime);
   fChain->SetBranchAddress("Data", &Data, &b_Data);
   fChain->SetBranchAddress("Sample", Sample, &b_Sample);
   Notify();
}

Bool_t mapping::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mapping::Show(Long64_t entry)
{
// PrInt_t contents of entry.
// If entry is not specified, prInt_t current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mapping::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void mapping::ResetOpt(){
  np=0;
  trig0=125;
  ene_thr=1;//1MeV
  det=0;
  fr=0;
  str=-1;
  er=0;
  flu=0;
  e=0;
  st=false;
  thr1=0;
  thr2=0;
  min=0;
  memset(ratio,0,sizeof(ratio));//trapezoid calibrate
  memset(te,0,sizeof(te));
  memset(rte,0,sizeof(rte));
  memset(cate,0,sizeof(cate));
  memset(slowfilter,0,sizeof(slowfilter));
  memset(trap,0,sizeof(trap));
  memset(trap1,0,sizeof(trap1));
  memset(fd,0,sizeof(fd));
  memset(sd,0,sizeof(sd));
  memset(msd,0,sizeof(msd));
  memset(tsp,0,sizeof(tsp));
  memset(msd1,0,sizeof(msd1));
  memset(tsp1,0,sizeof(tsp1));
  memset(fft,0,sizeof(fft));
  memset(fftre,0,sizeof(fftre));
  memset(fftim,0,sizeof(fftim));
  memset(fftb,0,sizeof(fftb));
  memset(fftp,0,sizeof(fftp));
  //showcanvas = new TCanvas();
  //offlinemultigraph = new TMultiGraph();

}

void mapping::ffttr(){
  Int_t n=750;
  TH1D *hst = new TH1D("hst","hst",n+1,0,750);
  for(Int_t i=0;i<750;i++) hst->SetBinContent(i+1,dwave[i]);
  TH1 *ht = 0;
  TVirtualFFT::SetTransform(0);
  ht=hst->FFT(ht, "MAG");
  Int_t n1=ht->GetEntries();
  ht->GetBinContent(100);
  
  //DC and Nyquist 
  Double_t re, im;
  //That's the way to get the current transform object:
  TVirtualFFT *fft1 = TVirtualFFT::GetCurrentTransform();
  //Use the following method to get just one point of the output
  fft1->GetPointComplex(0, re, im);
  fft1->GetPointComplex(n/2+1, re, im);
  nyq=re;
  
  Double_t *re_full = new Double_t[n];
  Double_t *im_full = new Double_t[n];
  fft1->GetPointsComplex(re_full,im_full);
  //frequency filter
  /*
  for(int i=0;i<n;i++){
    if(i>10) re_full[i]=0;
    if(i>10) im_full[i]=0;
  }*/
  
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
  fft_back->SetPointsComplex(re_full,im_full);
  fft_back->Transform();
  TH1 *hb = 0;
  hb = TH1::TransformHisto(fft_back,hb,"Re");
  
  //using array to get image and real part
  Double_t in[750];
  for(Int_t i=0;i<n;i++) in[i] = dwave[i];
  TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &n, "R2C ES K");
  if (!fft_own) return;
  fft_own->SetPoints(in);
  fft_own->Transform();
  //Copy all the output points
  fft_own->GetPoints(in);
  TH1 *hr = 0;
  hr=TH1::TransformHisto(fft_own, hr, "RE");
  TH1 *hm = 0;
  hm = TH1::TransformHisto(fft_own, hm, "IM");
  Double_t wm=0, ftbm=0;
  for(int i=0;i<750;i++){
  if(hb->GetBinContent(i+1)>ftbm) ftbm=hb->GetBinContent(i+1);
  if(dwave[i]>wm) wm=dwave[i];
  }
  TH1 *hp = 0;
  hp = hst->FFT(hp, "PH");
  
  //DC and Nyquist
  for(int i=0;i<750;i++){
  fft[i] = ht->GetBinContent(i+1);
  fftb[i] = hb->GetBinContent(i+1)*wm/ftbm;
  fftre[i] = hr->GetBinContent(i+1);
  fftim[i] = hm->GetBinContent(i+1);
  fftp[i] = hp->GetBinContent(i+1);
  //cout<<fft[i]<<"  "<<sqrt(fftre[i]*fftre[i]+fftim[i]*fftim[i])<<"  "<<fftre[i]<<"  "<<fftim[i]<<endl;
  }
  delete ht;
  delete hb;
  delete hst;
  delete fft1;
  delete re_full;
  delete im_full;
  delete fft_back;
  delete fft_own;
  delete hr;
  delete hm;
  delete hp;
}

Bool_t mapping::BranchOpt(TTree *opt){
  if(opt==0) return 0;
  opt->Branch("det",&det,"det/I");
  opt->Branch("fr",&fr,"fr/I");
  opt->Branch("str",&str,"str/I");
  opt->Branch("ts",&ts,"ts/l");
  opt->Branch("er",&er,"er/D");
  opt->Branch("e",&e,"e/D");
  opt->Branch("np",&np,"np/I");
  opt->Branch("tsp",&tsp,"tsp[np]/l");
  opt->Branch("eh",&eh,"eh[np]/D");
  opt->Branch("te",&te,"te[np]/D");
  opt->Branch("rte",&rte,"rte[np]/D");
  opt->Branch("cate",&cate,"cate[np]/D");
  opt->Branch("msd",&msd,"msd[np]/D");
  opt->Branch("flu",&flu,"flu/D");
  opt->Branch("slowfilter",slowfilter,"slowfilter[1500]/D");
#ifdef TEST
  opt->Branch("a2",&a2,"a2/D");
  opt->Branch("off2",&off2,"off2/I");
#endif
  char varp[120];
  opt->Branch("dwave",&dwave,"dwave[750]/D");
  opt->Branch("twave",&twave,"twave[750]/D");
  
  sprintf(varp,"rwave[%d]/S",NSAMPLE_V1724);
  opt->Branch("rwave",&rwave,varp);
  opt->Branch("fd",&fd,"fd[750]/D");
  opt->Branch("sd",&sd,"sd[750]/D");
  opt->Branch("trap",&trap,"trap[750]/D");
  opt->Branch("trap1",&trap1,"trap1[750]/D");
  sprintf(varp,"dp1[%d]/S",NSAMPLE_V1724);
  opt->Branch("dp1",&Data.Dprobe1,varp);
  sprintf(varp,"dp2[%d]/S",NSAMPLE_V1724);
  opt->Branch("dp2",&Data.Dprobe2,varp);
  opt->Branch("pu",&npu,"npu/I");
  opt->Branch("sample",sample,"sample[750]/S");
  opt->Branch("Sample",Sample,"Sample[1500]/S");
  opt->Branch("slope",&slope,"slope/D");
  opt->Branch("thr1",&thr1,"thr1/D");
  opt->Branch("thr2",&thr2,"thr2/D");
  opt->Branch("min",&min,"min/D");
  
  opt->Branch("bd",&boardn,"bd/I");
  opt->Branch("ch",&ChannelNo,"ch/b");
  opt->Branch("nevt",&nevt,"nevt/I");
  opt->Branch("slopechi2",&slopechi2,"slopechi2/D");
  opt->Branch("st",&st,"st/B");
  opt->Branch("eh1",&eh1,"eh1[np]/D");
  opt->Branch("bslope",&bslope,"bslope/D");
  opt->Branch("gfl",&gfl,"gfl/I");
  opt->Branch("fft",&fft,"fft[750]/D");
  opt->Branch("fftre",&fftre,"fftre[750]/D");
  opt->Branch("fftim",&fftim,"fftim[750]/D");    
  opt->Branch("fftb",&fftb,"fftb[750]/D"); 
  opt->Branch("fftp",&fftp,"fftp[750]/D"); 
  opt->Branch("nyq",&nyq,"nyq/D");
  //opt->Branch("dfl1",&dfl1,"dfl1/D");
  //opt->Branch("dhl",&dhl,"dhl/D");
  return 1;
}

Bool_t mapping::Mapping(){
 
  if(boardn>=0&&boardn<16){
    det=2;//dssd
    fr=2;
    str=boardn*8+ChannelNo;
  }
  else if(boardn>=16&&boardn<22) {
    det=2;
    fr=1;
    str=(boardn-16)*8+ChannelNo;
  }
  else if(boardn==22&&ChannelNo==5){
    det=1;//mwpc
  }
  else if(boardn==24&&ChannelNo>=0&&ChannelNo<5){
    //return 0; 
    det=4;//gamma detector
    str=ChannelNo;
  }
  else if(boardn==22&&ChannelNo<3&&ChannelNo>=0){
    det=3;//veto
    str=ChannelNo;
  }
  else if(boardn==23&&ChannelNo==0){
    det=2;
    fr=2;
    str=45;
  }
  else {
    return 0;
  }
  
  //if(ChannelNo!=0)return 0;
  ts=Data.TriggerTimeTag;
  er=Data.Energy;
  
  //if(er<=0) return 0;
  if(det==2){
  if(fr==1) e=er*ypar[str][0]+ypar[str][1];
  if(fr==2) e=er*xpar[str][0]+xpar[str][1];
  
  e=e*sl+cept;
  //cout<<te<<"  "<<e<<endl;
  }
  else e=er;
  //if(det!=2&&er>0)cout<<e<<"  "<<det<<endl;
  memcpy(rwave,Data.Probe2,sizeof(Short_t)*NSAMPLE_V1724);
  memcpy(wave,Data.Probe2,sizeof(Short_t)*NSAMPLE_V1724);
  memcpy(dp1,Data.Dprobe1,sizeof(Short_t)*NSAMPLE_V1724);
  memcpy(dp2,Data.Dprobe2,sizeof(Short_t)*NSAMPLE_V1724);
  
  for(Int_t i=0;i<750;i++){//for even sample
  sample[i] = i;
  Sample[i] = i;
  Sample[i+1] = i+1;
  dwave[i]=wave[i*2];
  }


  if(det==2 ) {
    if(!BaseLineCorrection()) return 0;    
    ffttr();//fft transform
#ifdef TEST
    DPulseGen();//for test
#endif

       if(!PuDet()) return 0;
       if(np>0) pulsedet();
       // if(np==0)cout<<det<<"  "<<e<<"  "<<np<<endl;
      //  if(!Pileup())return 0;
    
    //TrapezoidXia(40,20,300);
 }

  return 1;
}


void mapping::DPulseGen()
{
  Int_t t0=261;
  a2=ru->Uniform(0.1,0.9);
  off2=ru->Uniform(10,200);
  Double_t twave[750];
  for(Int_t i=0;i<750;i++) {
    if(i>off2) twave[i]=a2*dwave[i-off2];
    else twave[i]=0;
  }
  
  for(Int_t i=0;i<750;i++) dwave[i]+=twave[i];
}

bool judge(const pair<Double_t,int> a, const pair<Double_t,int> b)
{
  return a.first>b.first;
}

bool judge_t(const pair<Double_t,int> a, const pair<Double_t,int> b)
{
  return a.second<b.second;
}

bool compare(int a, int b)
{
  return a<b;
}

Bool_t mapping::PuDet()
{
   //1. no negative signal in pulse
   //2. maximum in flat part should be lower than minimum in peak area.
   //fr=1,str0-31
      Double_t rat=AmpRatio();  //avoid the difference of strips
   //   if(min<-50) return 0;//avoid distortion
     /* for(int i=0;i<750;i++)
        gra->SetPoint(i,i,dwave[i]);*/
     
      for(int i=0;i<750;i++) twave[i]=1.*dwave[i];//amplitude normalization
   
      //fast filter with small raise time
      Int_t n=1, gg=1;
      Int_t ll=n;
      gfl=0;
      
      Double_t s0=0,s1=0,tmax=0;
      for(int i=0; i<ll+1; i++) s0 += twave[i];
      for(int i=ll+gg; i<2*ll+gg+1; i++) s1 +=twave[i];
      Int_t gf=0;
      for(int i=0;i<750;i++) {
        if(i>650&&i<700&&dwave[i]<40) gf++;
	//fd,sd
	if(i>2) {
          fd[i]=(twave[i]-twave[i-2]);
	  sd[i]=fd[i]*fd[i]-fd[i-1]*fd[i-1];
	  // cout<<i<<" "<<sd[i]<<endl;
	}
        //trapz shape
	if(i>=2*ll+gg+1 && i<750) {
	  s0 += twave[i-ll-gg]-twave[i-2*ll-gg-1];
	  s1 += twave[i]-twave[i-ll-1];
	  int k=i-ll-gg/2;
	  trap[k]=(s1-s0)/ll;
	  if(trap[k]<0)trap[k]=0;
	  if(tmax<trap[k]) tmax=trap[k];
	}
      }
      if(gf>10) gfl=1;
      /*for(int i=100;i<750;i++)
          sd[i]=trap[i]*trap[i]-trap[i-1]*trap[i-1];*/
          
     //iterator element
      vector<pair<Double_t,int> > vpeak;
      for(int i=n+1;i<750-n-1;i++) {
	if(trap[i]<=0) continue;
	int bp=0;
	for(int j=1;j<3;j++)
	  if(trap[i]>=trap[i+j] && trap[i]>trap[i-j])
	    bp++;
	    
	if(bp==2) vpeak.push_back(make_pair(trap[i],i));
      }
      
      //select largest amplitude of trapezoid 
      sort(vpeak.begin(),vpeak.end(),judge);
      
      if(vpeak.size()>4)
        for(int i=4;i<6;i++)
          thr1+=vpeak[i].first/2;
        else thr1=40.;
      
      if(vpeak.size()>3)
        vpeak.erase(vpeak.begin()+3,vpeak.end());//remove redundant element
        
      np=0;
      // reduce fluctuation
      for(int i=0;i<int(vpeak.size());i++)
        if(vpeak[i].first>2.5*thr1){//trap:20,15
          tsp[np]=vpeak[i].second;
          np++;
          }

      sort(tsp,tsp+np,compare);//timestamp sort
      
      for(int i=0;i<np;i++)
        msd[i]=trap[tsp[i]];
      //if(np>1) for(int i=0;i<np;i++)cout<<tsp[i]<<"  "<<msd[i]<<endl;
      //if(np==0) cout<<vpeak[0].first<<"  "<<vpeak[1].first<<"  "<<vpeak[2].first<<"  "<<thr1<<"  "<<vpeak.size()<<endl;
      
      vpeak.clear();
      // another parameters
      n=10, gg=1;
      ll=n;
      // initialize trap
      vector<pair<Double_t,int> > vpeak1;
      s0=0,s1=0,tmax=0;
      for(int i=0; i<ll+1; i++) s0 += twave[i];
      for(int i=ll+gg; i<2*ll+gg+1; i++) s1 +=twave[i];
      
      for(int i=0;i<750;i++) {
        //trapz shape
	if(i>=2*ll+gg+1 && i<750) {
	  s0 += twave[i-ll-gg]-twave[i-2*ll-gg-1];
	  s1 += twave[i]-twave[i-ll-1];
	  int k=i-ll-gg/2;
	  trap1[k]=(s1-s0)/ll;
	  if(trap1[k]<0) trap1[k]=0;
	  if(tmax<trap1[k]) tmax=trap1[k];
	}
      }
     //iterator element
      for(int i=n+1;i<750-n-1;i++) {
	if(trap1[i]<=0) continue;
        int bp=0;
	for(int j=1;j<3;j++)
	  if(trap1[i]>=trap1[i+j] && trap1[i]>trap1[i-j])
	    bp++;
	if(bp==2) vpeak1.push_back(make_pair(trap1[i],i));
      }
      
       //select largest amplitude of trap
      sort(vpeak1.begin(),vpeak1.end(),judge);
      
      if(vpeak1.size()>3)
        thr2=vpeak1[3].first;
        else thr2=3.;
      
      if(vpeak1.size()>3)
        vpeak1.erase(vpeak1.begin()+3,vpeak1.end());//remove redundant element 
      
      np1=0;
      for(int i=0;i<int(vpeak1.size());i++) {
        if(fr==2&&vpeak1[i].first<3.0||vpeak1[i].first/vpeak1[i+1].first<1.5)//threshold
          vpeak1.erase(vpeak1.begin()+i,vpeak1.end());
        if(fr==1&&vpeak1[i].first<3.0||vpeak1[i].first/vpeak1[i+1].first<1.5)//threshold
          vpeak1.erase(vpeak1.begin()+i,vpeak1.end());
        if(vpeak1[i+1].first==0) vpeak1.erase(vpeak1.begin()+i,vpeak1.end());
      }
      sort(vpeak1.begin(),vpeak1.end(),judge_t);//timestamp sort
      
      for(int i=0;i<int(vpeak1.size());i++) {
        msd1[i]=vpeak1[i].first;
        tsp1[i]=vpeak1[i].second;
        np1++;
      }
      // cout<<np1<<endl;
      vpeak1.clear();
      
      //combine the results
      for(Int_t i=0;i<np;i++)
        for(Int_t j=0;j<np1;j++){
          if(i+1<np&&tsp1[j]>tsp[i]&&tsp1[j]<tsp[i+1]){
            if(tsp1[j]-tsp[i]<4||tsp[i+1]-tsp1[j]<4) continue; //tsp1[j] the same as tsp[i] or tsp[i+1]
            for(Int_t k=np;k>i;k--){
              tsp[k]=tsp[k-1];
              msd[k]=msd[k-1];
              }
              tsp[i+1]=tsp1[j];
              msd[i+1]=msd1[j];
              i++;
              np++;
          }
          else if(i+1>=np&&tsp1[j]>tsp[i]){
            if(tsp1[j]-tsp[i]<4) continue; //tsp1[j] the same as tsp[i]
            tsp[i+1]=tsp1[j];
            msd[i+1]=msd1[j];
            i++;
            np++;
            }
          }
        //fluctuation estimation
        Int_t ta=tsp[np-1]+1,tb=740;
        flu=10.;
        
        if(ta<730)
        if(np>1){
        Int_t mm=2;
        Int_t dt=tb-ta;
        for(int i=1;i<np;i++) 
         if(dt<tsp[i]-tsp[i-1]) {
           ta=tsp[i-1]+mm;
           tb=ta+5;
           dt=tsp[i]-tsp[i-1];
          }
         }
          else if(np==1){
              ta=tsp[0]+2;
              tb=ta+5;
              }
              else{ 
         	ta=150;
         	tb=ta+5;
         	}
           else{
             ta=734;
             tb=739;
           }
         	
        for(int i=ta;i<tb;i++)
          flu+=abs(dwave[i+1]-dwave[i])/(tb-ta);
    //amplitude estimation and filter
     double h0[10],h1[10],hend=0;
     for(int i=745;i<750;i++) hend+=dwave[i]/5.;
     h0[np]=hend;
     Int_t anp=np;
     np=0;
    for(int i=0;i<anp;i++) {
      h0[i]=0;
      h1[i]=0;
      double af,bf;
      af=0.;
      bf=0.;
      int r=tsp[i]+5;
      if(tsp[i]>742) r=749;
      for(int j=tsp[i]-4;j<tsp[i];j++)
        af+=dwave[j]/4.;
      for(int j=tsp[i]+1;j<r;j++)
        bf+=dwave[j]/(r-tsp[i]-1);
        
      double mfl,nfl,mfl1,nfl1;
      int rl,rh;
      mfl=dwave[tsp[i]+3];
      nfl=dwave[tsp[i]+3];
      mfl1=dwave[tsp[i]+8];
      nfl1=dwave[tsp[i]+8];
      
      rl=tsp[i]+1;
      if(tsp[i]+8<749)
        rh=tsp[i]+8;
        else rh=749;
      if(i>0)
        for(int j=rl;j<rh;j++){
          if(mfl<dwave[j]) mfl=dwave[j];
          if(nfl>dwave[j]) nfl=dwave[j];
          }
      rl=tsp[i]+3;
      if(tsp[i]+15<749)
        rh=tsp[i]+15;
        else rh=749;
      if(i>0)
        for(int j=rl;j<rh;j++){
          if(mfl1<dwave[j]) mfl1=dwave[j];
          if(nfl1>dwave[j]) nfl1=dwave[j];
        }
      dfl=mfl-nfl;
      dfl1=mfl1-nfl1;  
      double mhl,nhl;
      mhl=dwave[tsp[i]-3];
      nhl=dwave[tsp[i]-3];
      
      if(i>0)
        for(int j=tsp[i]-7;j<tsp[i]-1;j++){
          if(mhl<dwave[j]) mhl=dwave[j];
          if(nhl>dwave[j]) nhl=dwave[j];
          }
      dhl=mhl-nhl;
      int n0 = tsp[i]-2;
      int n1 = tsp[i]+4;
      
      if (i>0&&np>0) n0 = (tsp[np-1]+2 > tsp[i]-2) ? (tsp[np-1]+1) : (tsp[i]-3);
      if (i>0 && i+1<anp) n1 = (tsp[i+1]-1 < tsp[i]+4) ? (tsp[i+1]-1) : (tsp[i]+4);
      for(int j=n0;j<tsp[i]-1;j++) h0[i] += dwave[j] / (tsp[i]-1-n0);
      for(int j=tsp[i]+2;j<n1;j++) h1[i] += dwave[j] / (n1-tsp[i]-2);
      
      if(i==0&&anp>1&&tsp[anp-1]-tsp[i]==3) h1[i]=dwave[tsp[i]+1];
      
       if(h0[i]<0) h0[i]=0;
       
       eh[i]=h1[i]-h0[i];
       
       //cout<<i<<"  "<<flu<<"  "<<eh[i]<<"  "<<tsp[i]<<"  "<<bf-af<<"  "<<mhl-nfl<<"  "<<mfl-nfl<<"  "<<mfl1-nfl1<<"  "<<mhl-nhl<<"  "<<eh[i]/eh[0]<<"  "<<msd[i]<<"  "<<thr1<<"  "<<dwave[tsp[i]+1]<<endl;

       if(eh[i]<1.2*flu) continue; //w186 0.5*flu;1.2
       
       if(np==0){//the first raise
         if(dwave[tsp[i]]>400.&&dwave[tsp[i]]<4.*flu) continue;
         if(dwave[tsp[i]]<=400.&&dwave[tsp[i]+1]<2.*flu) continue;//very little amplitude pulses
         }
       if(np>0){//the raise follwed the first        
           if(154>tsp[i]){
             if(eh[i]<2.5*flu) continue;
             if(eh[i]/eh[0]<0.02) continue;
             if(eh[i]<8.*flu&&(msd[i]<4.*thr1||mhl-nfl>3.*flu||mfl-nfl>3.*flu)) continue;//for big fluctuation after the raise
             }
           else if(tsp[i]>=154){
             if(dwave[tsp[i]+4]-dwave[tsp[i]-1]<0.3*flu) continue;     
                   
             if(dwave[tsp[i]+2]-dwave[tsp[i]-1]>4.*flu){
               tsp[np]=tsp[i];
               eh[np]=eh[i];
               msd[np]=msd[i];
               np++;  
               //if(np>1) cout<<np<<"  "<<eh[i]<<"  "<<tsp[i]<<"  "<<msd[i]<<"  "<<thr1<<"  "<<flu<<"  "<<mfl-nfl<<"  "<<mhl-nhl<<"  "<<eh[i]/eh[i-1]<<"  "<<ts<<"  "<<abs(mhl-nfl)<<endl;
               continue;             
               }
             if(((tsp[i+1]==0||tsp[i+1]-tsp[i]>=10)&&tsp[i]-tsp[np-1]>=10)&&(mfl-nfl>2.1*flu||mhl-nhl>1.6*flu)) continue;            
             if(((tsp[i+1]==0||tsp[i+1]-tsp[i]>=15)&&tsp[i]-tsp[np-1]>=10)&&mfl1-nfl1>1.6*flu) continue;
             if(tsp[i]-tsp[np-1]<10&&(mhl-dwave[tsp[i]-2]>4.*flu||eh[i]<4.*flu)) continue;
             if(msd[i]<0.6*thr1) continue;
          // cout<<tsp[i]<<"  "<<eh[i]<<"  "<<tsp[i+1]<<"  "<<tsp[i]<<"  "<<tsp[i+1]-tsp[i]<<"  "<<mfl1<<"  "<<nfl1<<"  "<<mfl1-nfl1<<"  "<<flu<<endl;
           }
         }
       
       tsp[np]=tsp[i];
       eh[np]=eh[i];
       msd[np]=msd[i];
       np++;
         //if(np>1) cout<<np<<"  "<<eh[i]<<"  "<<tsp[i]<<"  "<<msd[i]<<"  "<<thr1<<"  "<<flu<<"  "<<mfl-nfl<<"  "<<mhl-nhl<<"  "<<eh[i]/eh[i-1]<<"  "<<ts<<"  "<<abs(mhl-nfl)<<endl;
    }
       if(np>=0&&np<5)
         return 1;
       else
         return 0;
}

Double_t mapping::AmpRatio()
{

//from map00575.root
//tree->Draw("eh[0]:str>>(128,0,128,1500,0.1500)","np==1&&fr==2","colz")
//fr==2, reference 1100.
//fr==1, reference 250. x 4.4=1100
// er is not correct for estimating height ratio among differnt strips.
 /* if(det==2&&fr==2) {
    if(str>15&&str<32) return 1100./820.;
    if(str>111&&str<128) return 1100/250.;
    return 1.;
  }

  if(det==2&&fr==1) {
    if(str==0 || str>31&&str<48) return 250./185.*4.4;
    return 4.4;
  }*/
  //Amplify for different strips, 
  
  for(Int_t i=120;i<140;i++)
    if(dwave[i]<min) min=dwave[i];
  
  if(det==2&&fr==2) {
    if(str>15&&str<32) return 55./14.5;
    if(str>111&&str<128) return 55./17.;
    return 55./19.5;
  }

  if(det==2&&fr==1) {
    if(str==0 || str>31&&str<48) return 55./13.;
    return 55./17.5;
  }
  //average amplitude normalie to 3000.
  /*
  Double_t ave=0;
  for(Int_t i=150;i<750;i++)
    ave+=dwave[i]/600.;*/
  
  //return 3000./ave;
}
  TGraph *gra1 = new TGraph;

//perform the baseline correction and determined trigger positon from dp1
bool mapping::BaseLineCorrection()
{
     // for(Int_t i=0;i<750;i++) cout<<i*2<<" "<<dwave[i*2]<<endl;
     Int_t i0=trig0-10;
     Double_t ave=0;

     for(Int_t j=0;j<i0;j++) ave+= dwave[j];
     ave /=i0;
     npu=0;
     for(Int_t i=0;i<750;i++) {
       dwave[i]-=ave;
       if(dp1[i]>0) npu++;
    //   if(i>trig0+50 && dwave[i]<thr()) return 0; //abnormal fluctuation due to noise.
     }
     for(Int_t i=0;i<trig0-10;i++)
     dwave[i]=0;//for superpulse fitting
     
     return 1;
}

void mapping::pulsedet()
{
//slope calculation
 Double_t tmp[2];
 Int_t count;

 for(Int_t i=0;i<749;i++)
 {
     count=0;
     gr->SetPoint(i,i,dwave[i]);
     
     if(i>140)
     if(dwave[i-1]==dwave[i]) {
       for(Int_t k=0;k<30;k++)
       if(dwave[i-1]==dwave[i-1+k])
       count++;
     }
     
     if(count==30)st=true;
     }
     TF1 *f = new TF1("f","expo",tsp[np-1]+5,749);
 //    f->SetParameters(0,3e-5);
//     cout<<st<<"  "<<count<<endl;
     gr->Draw();
     gr->Fit(f,"RQ");

     f->GetParameters(&tmp[0]);
     slopechi2=f->GetChisquare()/f->GetNDF();
     //f->GetNdf();
     slope=-tmp[1];
     //cout<<slope<<endl;
     gr->Clear();
     for(Int_t i=0;i<120;i++)
       gra1->SetPoint(i,i,rwave[i]);
     TF1 *ff = new TF1("ff","expo",0,120);
     gra1->Fit(ff,"RQ");
     ff->GetParameters(&tmp[0]);
     bslope=-tmp[1];
     gra1->Clear();
     delete f;
     delete ff;
}

void mapping::TrapezoidXia(Int_t SlowLen,Int_t SlowGap,double preamptau)
{
  unsigned int esum0[32768], esum1[32768], esum2[32768];
  unsigned int offset, x, y;
  double deltaT;
  Int_t RcdTraceLength=750;
  double b1, c0, c1, c2;
  unsigned int bsum0, bsum1, bsum2;
  double baseline;

  doublesample = new double[RcdTraceLength];
  doubledata = new double[RcdTraceLength];
  doubleslowfilter = new double[RcdTraceLength];
  doublebline = new double[RcdTraceLength];
  
//Get raw data for file
  for(Int_t i=0;i<RcdTraceLength;i++)
  {
  doubledata[i]=(double)Data.Probe2[i*2];
  doublesample[i]=(double)i;
  }
  
//computer slow filter coefficients
  deltaT = 1.0/((double)50);//presetting
  b1 = exp(-1.0*deltaT/preamptau);
  c0 = -(1.0 - b1) * pow(b1, (double)SlowLen) * 4.0 / (1.0 - pow(b1, (double)SlowLen));
  c1 = (1.0 - b1) * 4.0;
  c2 = (1.0 - b1) * 4.0 / (1.0 - pow(b1, (double)SlowLen));
  
 // cout<<SlowLen<<"  "<<SlowGap<<"  "<<SlowFilterRange<<"  "<<preamptau<<"  "<<b1<<"  "<<c0<<"  "<<c1<<"  "<<c2<<endl;
 // Compute baseline
  bsum0 = 0;
  for(y=0; y<SlowLen; y++)
    {
      bsum0 += doubledata[y];
    }
  bsum1 = 0;
  for(y=SlowLen; y<SlowLen+SlowGap; y++)
    {
      bsum1 += doubledata[y];
    }
  bsum2 = 0;
  for(y=(SlowLen+SlowGap); y<(2*SlowLen+SlowGap); y++)
    {
      bsum2 += doubledata[y];
    }
  baseline = c0 * (double)bsum0 + c1 * (double)bsum1 + c2 * (double)bsum2;

  // Compute slow filter response
  offset = 2*SlowLen + SlowGap - 1;
  for(x=offset; x<RcdTraceLength; x++)
    {
      esum0[x] = 0;
      for(y=(x-offset); y<(x-offset+SlowLen); y++)
	{
	  esum0[x] += doubledata[y];
	}
      esum1[x] = 0;
      for(y=(x-offset+SlowLen); y<(x-offset+SlowLen+SlowGap); y++)
	{
	  esum1[x] += doubledata[y];
	}
      esum2[x] = 0;
      for(y=(x-offset+SlowLen+SlowGap); y<(x-offset+2*SlowLen+SlowGap); y++)
	{
	  esum2[x] += doubledata[y];
	}
      doubleslowfilter[x] = c0 * (double)esum0[x] + c1 * (double)esum1[x] + c2 * (double)esum2[x] - baseline;
    }

  for(x=0; x<offset; x++)
    {
      doubleslowfilter[x] = doubleslowfilter[offset];
    }
  
  //baseline and slowfilter for comparision
  for(Int_t i=0;i<RcdTraceLength;i++)
  {
    doublebline[i]=doubleslowfilter[offset];
    slowfilter[i] = doubleslowfilter[i];
    }
  
  //get flat top value 
  Norm();
  double tmpe;
  for(Int_t i=0;i<np;i++)
  {
   // tmpe=0;
   // for(Int_t j=tsp[i]+SlowLen;j<tsp[i]+SlowLen+SlowGap;j++)
   //   if(slowfilter[j]>tmpe)tmpe=slowfilter[j];
   
   for(int j=tsp[i]+SlowLen+2;j<tsp[i]+SlowLen+SlowGap-2;j++)
    rte[i]+=slowfilter[j];
    
    rte[i]/=SlowGap-4;
 
    if(fr==1){
      te[i]=rte[i]*ypar_xia[str][0]+ypar_xia[str][1];
     // cout<<fr<<"  "<<str<<"  "<<ypar_xia[str][0]<<"  "<<ypar[str][1]<<endl;
      }
    if(fr==2){
      te[i]=rte[i]*xpar_xia[str][0]+xpar_xia[str][1];
     // cout<<fr<<"  "<<str<<"  "<<xpar_xia[str][0]<<"  "<<xpar[str][1]<<endl;
      }
    //te[i]=ratio[1]*rte[i]+ratio[0];
    
    if(fr==1)cate[i]=te[i]*ypar[str][0]+ypar[str][1];
    if(fr==2)cate[i]=te[i]*xpar[str][0]+xpar[str][1];
    
    cate[i]=cate[i]*sl+cept;//
  }
  
  //Draw filter 
  /*rawdata = new TGraph(RcdTraceLength,doublesample,doubledata);
  
  gslowfilter = new TGraph(RcdTraceLength,doublesample,doubleslowfilter);
  
  bline = new TGraph(RcdTraceLength,doublesample,doublebline);
  
  showcanvas->cd();
  showcanvas->Clear();
  gslowfilter->SetLineColor(3);
  rawdata->SetLineColor(4);
  bline->SetLineColor(6);
  offlinemultigraph->Clear();
  offlinemultigraph->Add(rawdata);
  offlinemultigraph->Add(gslowfilter);
  offlinemultigraph->Add(bline);
  offlinemultigraph->Draw("AL");
  cout<<tsp[0]<<endl;
  showcanvas->SaveAs("c1.eps");
  sleep(10);*/
  
  
  delete doublesample;
  delete doubledata;
  delete doubleslowfilter;
  delete doublebline;
}


 void mapping::Norm()
 {
 if(fr==1)
 {
   ratio[0]=16.44;
   ratio[1]=1.77;
   }
   else if(fr==2&&str>111&&str<128)
   {
   ratio[0]=-1.154;
   ratio[1]=1.791;
   }
   else
   {
   ratio[0]=1.673;
   ratio[1]=0.4457;
   }
   
   }

void mapping::AveraDif()
{
  double h0[10],h1[10];
  for(Int_t i=0;i<np;i++)
  {
    h0[i]=0;
    h1[i]=0;
    
    Int_t n0 = tsp[i]-9;
    Int_t n1 = tsp[i]+10;
    if(i==0&&tsp[i+1]-tsp[i]<17) n1 = tsp[i+1]-2;
    
    if(i>0&&tsp[i]-tsp[i-1]<17) n0 = tsp[i-1]+4;
    
    for(Int_t j=n0;j<tsp[i]-2;j++) h0[i] += dwave[j] / Int_t(tsp[i]-2-n0);
    for(Int_t j=tsp[i]+4;j<n1;j++) h1[i] += dwave[j] / Int_t(n1-tsp[i]-4);
    if(i==0) h0[i]=0;
    eh1[i]=h1[i]-h0[i];
    }

}
