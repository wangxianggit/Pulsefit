//define HIST in mapping.h
#include "mapping.h"
#include <iostream>

using namespace std;

int main(int argc, char** argv){
        int runnum;
	if(argc!=2 ) {
	        cout<<"USAGE1 : "<<argv[0]<<" #[runnum] for all boards"<<endl;
		return -1;
	}
	char ipfn[124];
	char opfn[124];
	runnum=atoi(argv[1]);
	//runnum1=atoi(argv[2]);
	sprintf(ipfn,"/data/d1/SHANS_201605_223/data20180805/run%05d_wave.root",runnum);
	TFile *ipf=new TFile(ipfn);
	if(!ipf->IsOpen()) return 1;
        sprintf(opfn,"/data/d1/SHANS_201605_223/mapping/map%05d.root",runnum); 
	//sprintf(opfn,"./map%05d.root",runnum); 
	TFile *opf=new TFile(opfn,"RECREATE");
	TTree *opt=new TTree("tree","tree");
	mapping *mp=new mapping(opt);
        mp->Readpar();
        mp->Readpar_xia();
#ifdef HIST
        mp->SetHist();//create histogram for superpulse
#endif
	for(int i=0;i<27;i++){//27
		char tn[124];
		sprintf(tn,"DGTZ%02d",i);
		TTree *ipt=(TTree*)ipf->Get(tn);
		if(ipt==0) {
			cout<<"Cannot find tree: "<<tn<<endl;
			continue;
		}
		opf->cd();
		mp->Loop(ipt,i,opt);
		cout<<"Board N: "<<i<<" done! run"<<runnum<<endl;
	}
	opt->Write();

#ifdef HIST
        mp->WriteHist(runnum);
#endif
      	opf->Close();
	return 0;
}
