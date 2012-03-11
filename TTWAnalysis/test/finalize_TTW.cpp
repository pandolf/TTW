
#include "Ntp1Finalizer_TTW.h"
#include "TMath.h"
#include <iostream>







int main( int argc, char* argv[] ) {

  if( argc!=3 && argc!=4 && argc!=5 ) {
    std::cout << "USAGE: ./finalize_TTW [dataset] [selectionType] [bTaggerType=\"SSVHE\"] [leptType=\"ALL\"]" <<std::endl;
    return 13;
  }


  std::string dataset(argv[1]);
  std::string selectionType(argv[2]);

  std::string bTaggerType="SSVHE";
  if( argc>3 ) {
    std::string bTaggerType_str(argv[3]);
    bTaggerType = bTaggerType_str;
  }

  std::string leptType="ALL";
  if( argc>4 ) {
    std::string leptType_str(argv[4]);
    leptType = leptType_str;
  }



  Ntp1Finalizer_TTW* nf = new Ntp1Finalizer_TTW( dataset, selectionType, bTaggerType, leptType );


  if( dataset=="DATA_HR11_v2" ) {
   
    nf->addFile("DoubleMu_Run2011A_FULL"); //first muons! important!
    nf->addFile("DoubleMu_Run2011B_v2"); //first muons! important!
    nf->addFile("SingleMu_Run2011A_FULL"); //first muons! important!
    nf->addFile("SingleMu_Run2011B_v2"); //first muons! important!
    nf->addFile("DoubleElectron_Run2011A_FULL");
    nf->addFile("DoubleElectron_Run2011B_v2");

  } else if( dataset=="VV_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1" ) {

    nf->addFile("ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1");
    nf->addFile("WZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1");
    nf->addFile("WW_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1");

  } else {
  
    nf->addFile( dataset );

  }

  nf->finalize();


  return 0;

}

