#include "THcLADGEMModule.h"
#include "THcGlobals.h"
#include "THcLADGEM.h"
#include "THcLADSpectrometer.h"
#include "THcParmList.h"
#include "THaRunBase.h"

using namespace std;

ClassImp(THcLADGEMModule)

    //____________________________________________________________________________________
    THcLADGEMModule::THcLADGEMModule(const char *name, const char *description, Int_t imod, THaDetectorBase *parent)
    : THaSubDetector(name, description, parent) {

  fModuleID = imod;

  // Default values of MPD map parameters
  fN_APV25_CHAN    = 128; // number of channel per apv
  fN_MPD_TIME_SAMP = 6;   // mpd time sample
  fMPDMAP_ROW_SIZE = 9;

  fNClus = 0;

  fMAX2DHITS = 10000;

  fIsMC = false;
}

//____________________________________________________________________________________
THcLADGEMModule::~THcLADGEMModule() {
  // Default desctructor
  Clear();
  if (fIsSetup)
    RemoveVariables();
}

//____________________________________________________________________________________
void THcLADGEMModule::Clear(Option_t *opt) {

  //  cout << "THcLADGEMModule::Clear" << endl;
  THaSubDetector::Clear(opt);

  fNstrips_hit         = 0;
  fNstrips_hitU        = 0;
  fNstrips_hitV        = 0;
  fNstrips_hitU_neg    = 0;
  fNstrips_hitV_neg    = 0;
  fNdecoded_ADCsamples = 0;
  fIsDecoded           = false;

  // numbers of strips passing basic zero suppression thresholds and timing cuts:
  fNstrips_keep  = 0;
  fNstrips_keepU = 0;
  fNstrips_keepV = 0;

  // numbers of strips passing basic zero suppression thresholds, timing cuts, and higher max. sample and strip sum
  // thresholds for
  //  local max:
  fNstrips_keep_lmax  = 0;
  fNstrips_keep_lmaxU = 0;
  fNstrips_keep_lmaxV = 0;

  fTrigTime = 0.0;

  fCM_online.assign(fN_MPD_TIME_SAMP, 0.0);

  fxcmin.clear();
  fxcmax.clear();
  fycmin.clear();
  fycmax.clear();

  fClustersU.clear();
  fClustersV.clear();

  fBadStrips.clear();

  fNClus = 0;

  fN2Dhits = 0;
}

//____________________________________________________________________________________
THaAnalysisObject::EStatus THcLADGEMModule::Init(const TDatime &date) {

  //  cout << "THcLADGEMModule::Init()" << endl;
  EStatus status;
  if ((status = THaSubDetector::Init(date)))
    return fStatus = status;

  fParent = GetParent();

  TString appname = (static_cast<THaDetector *>(GetParent()) )->GetApparatus()->GetName();
  appname.ReplaceAll(".","_");
  appname += "_";
  TString detname = GetParent()->GetName();
  detname.Prepend(appname);
  detname.ReplaceAll(".","_");
  detname += "_";
  detname += GetName();
  fPedHistosInitialized = kFALSE;
  fPedDiagHistosInitialized = kFALSE;
  fCommonModePlotsInitialized = false;
  std::cout << "THcLADGEMModule::Init() - Initializing histograms for " << detname.Data() << std::endl;

  if( (fPedestalMode || fMakeCommonModePlots) && !fPedHistosInitialized ){ //make pedestal histograms:

    //Procedure:
    // 1. Analyze pedestal data with no pedestals supplied from the database.
    // 2. Extract "coarse" strip offsets from raw ADC histograms (i.e., common-mode means)
    // 3. Extract "fine" strip offsets from common-mode-subtracted histograms 
    // 3. Add "coarse" strip offsets (common-mode means) back into "common-mode subtracted" ADC 

    //U strips:
    hpedrmsU_distribution = new TH1D( TString::Format( "hpedrmsU_distribution_%s", detname.Data() ), "Pedestal RMS distribution, U strips; Ped. RMS", 200, 0, 100 );
    hpedmeanU_distribution = new TH1D( TString::Format( "hpedmeanU_distribution_%s", detname.Data() ), "Pedestal mean distribution, U strips; Ped. mean", 250, -250, 250 );
    hpedrmsU_by_strip = new TH1D( TString::Format( "hpedrmsU_by_strip_%s", detname.Data() ), "Pedestal rms U by strip; strip; ped. RMS", fNstripsU, -0.5, fNstripsU - 0.5 );
    hpedmeanU_by_strip = new TH1D( TString::Format( "hpedmeanU_by_strip_%s", detname.Data() ), "Pedestal mean U by strip; strip; ped. mean", fNstripsU, -0.5, fNstripsU - 0.5 );

    //V strips:
    hpedrmsV_distribution = new TH1D( TString::Format( "hpedrmsV_distribution_%s", detname.Data() ), "Pedestal RMS distribution, V strips; Ped. RMS", 200, 0, 100 );
    hpedmeanV_distribution = new TH1D( TString::Format( "hpedmeanV_distribution_%s", detname.Data() ), "Pedestal mean distribution, V strips; Ped. mean", 250, -250, 250 );
    hpedrmsV_by_strip = new TH1D( TString::Format( "hpedrmsV_by_strip_%s", detname.Data() ), "Pedestal rms V by strip; strip; Ped. RMS", fNstripsV, -0.5, fNstripsV - 0.5 );
    hpedmeanV_by_strip = new TH1D( TString::Format( "hpedmeanV_by_strip_%s", detname.Data() ), "Pedestal mean V by strip; strip; Ped. mean", fNstripsV, -0.5, fNstripsV - 0.5 );

    if( !fPedestalMode ){ //fill the above histograms with the pedestal info loaded from the database:
      for( UInt_t istrip = 0; istrip<fNstripsU; istrip++ ){
	      hpedmeanU_distribution->Fill( fPedestalU[istrip] );
	      hpedrmsU_distribution->Fill( fPedRMSU[istrip] );
	      hpedmeanU_by_strip->SetBinContent( istrip+1, fPedestalU[istrip] );
	      hpedrmsU_by_strip->SetBinContent( istrip+1, fPedRMSU[istrip] );	
      }

      for( UInt_t istrip=0; istrip<fNstripsV; istrip++ ){
      	hpedmeanV_distribution->Fill( fPedestalV[istrip] );
    	  hpedrmsV_distribution->Fill( fPedRMSV[istrip] );
      	hpedmeanV_by_strip->SetBinContent( istrip+1, fPedestalV[istrip] );
      	hpedrmsV_by_strip->SetBinContent( istrip+1, fPedRMSV[istrip] );
      }
      
    }
    hrawADCs_by_stripU = new TH2D( TString::Format( "hrawADCs_by_stripU_%s", detname.Data() ), "Raw ADCs by U strip number, no corrections; U/X strip; Raw ADC",
				   fNstripsU, -0.5, fNstripsU-0.5,
				   512, -0.5, 4095.5 );
    hrawADCs_by_stripV = new TH2D( TString::Format( "hrawADCs_by_stripV_%s", detname.Data() ), "Raw ADCs by V strip number, no corrections; V/Y strip; Raw ADC",
				   fNstripsV, -0.5, fNstripsV-0.5,
				   512, -0.5, 4095.5 );

    hrawADCs_by_stripU_nopedsub = new TH2D( TString::Format( "hrawADCs_by_stripU_nopedsub_%s", detname.Data() ), "Raw ADCs by U strip number, no corrections; U/X strip; Raw ADC",
				   fNstripsU, -0.5, fNstripsU-0.5,
				   512, -0.5, 4095.5 );
    hrawADCs_by_stripV_nopedsub = new TH2D( TString::Format( "hrawADCs_by_stripV_nopedsub_%s", detname.Data() ), "Raw ADCs by V strip number, no corrections; V/Y strip; Raw ADC",
				   fNstripsV, -0.5, fNstripsV-0.5,
				   512, -0.5, 4095.5 );
    
    
    hcommonmode_subtracted_ADCs_by_stripU = new TH2D( TString::Format( "hpedestalU_%s", detname.Data() ), "ADCs by U strip number, w/common mode correction, no ped. subtraction; U/X strip; ADC - Common-mode",
						      fNstripsU, -0.5, fNstripsU-0.5,
						      500, -500.0, 5000.0 );
    hcommonmode_subtracted_ADCs_by_stripV = new TH2D( TString::Format( "hpedestalV_%s", detname.Data() ), "ADCs by V strip number, w/common mode correction, no ped. subtraction; V/Y strip; ADC - Common-mode",
						      fNstripsV, -0.5, fNstripsV-0.5,
						      500, -500.0, 5000.0 );

    hpedestal_subtracted_ADCs_by_stripU = new TH2D( TString::Format( "hADCpedsubU_%s", detname.Data() ), "Pedestal and common-mode subtracted ADCs by U strip number; U/X strip; ADC - Common-mode - pedestal",
						    fNstripsU, -0.5, fNstripsU-0.5,
						    500,-500.,4500. );
    hpedestal_subtracted_ADCs_by_stripV = new TH2D( TString::Format( "hADCpedsubV_%s", detname.Data() ), "Pedestal and common-mode subtracted ADCs by V strip number; V/Y strip; ADC - Common-mode - pedestal",
						    fNstripsV, -0.5, fNstripsV-0.5,
						    500,-500.,4500. );

    hpedestal_subtracted_rawADCs_by_stripU = new TH2D( TString::Format( "hrawADCpedsubU_%s", detname.Data() ), "ADCs by U strip, ped-subtracted, no common-mode correction; U/X strip; ADC - pedestal",
						       fNstripsU, -0.5, fNstripsU-0.5,
						       500,-500.,4500. );
    hpedestal_subtracted_rawADCs_by_stripV = new TH2D( TString::Format( "hrawADCpedsubV_%s", detname.Data() ), "ADCs by V strip, ped-subtracted, no common-mode correction; V/Y strip; ADC - pedestal",
						       fNstripsV, -0.5, fNstripsV-0.5,
						       500,-500.,4500. );


  
  
    hpedestal_subtracted_rawADCsU = new TH1D( TString::Format( "hrawADCpedsubU_allstrips_%s", detname.Data() ), "distribution of ped-subtracted U strip ADCs w/o common-mode correction; ADC - pedestal", 1250, -500.,4500. );
    hpedestal_subtracted_rawADCsV = new TH1D( TString::Format( "hrawADCpedsubV_allstrips_%s", detname.Data() ), "distribution of ped-subtracted V strip ADCs w/o common-mode correction; ADC - pedestal", 1250, -500.,4500. );

    fPedHistosInitialized = true;
  }
  if( (fPedestalMode || fMakeCommonModePlots) && !fPedDiagHistosInitialized ){

    //Let's make versions of these for "all" events and "regular" full readout events only:
    hpedestal_subtracted_ADCsU = new TH1D( TString::Format( "hADCpedsubU_allstrips_%s", detname.Data() ), "full readout: ped-subtracted U strip ADCs w/common-mode correction; ADC - Common-mode - pedestal", 1250, -500.,4500. );
    hpedestal_subtracted_ADCsV = new TH1D( TString::Format( "hADCpedsubV_allstrips_%s", detname.Data() ), "full readout: ped-subtracted V strip ADCs w/common-mode correction; ADC - Common-mode - pedestal", 1250, -500.,4500. );

    hpedestal_subtracted_ADCsU_goodCM = new TH1D( TString::Format( "hADCpedsubU_allstrips_%s_goodCM", detname.Data() ), "full readout: ped-subtracted U strip ADCs w/common-mode correction; ADC - Common-mode - pedestal", 1250, -500.,4500. );
    hpedestal_subtracted_ADCsV_goodCM = new TH1D( TString::Format( "hADCpedsubV_allstrips_%s_goodCM", detname.Data() ), "full readout: ped-subtracted V strip ADCs w/common-mode correction; ADC - Common-mode - pedestal", 1250, -500.,4500. );

    
    UInt_t nAPVs_U = fNstripsU/fN_APV25_CHAN;
    hcommonmode_mean_by_APV_U = new TH2D( TString::Format( "hCommonModeMean_by_APV_U_%s", detname.Data() ), "full readout events: common-mode means for U strips; APV card; Common-mode", nAPVs_U, -0.5, nAPVs_U-0.5, 1024, -0.5, 4095.5 );

    hcommonmode_online_by_APV_U = new TH2D( TString::Format( "hCommonModeOnline_by_APV_U_%s", detname.Data() ), "Online ZS events: common-mode means for U strips; APV card; online common-mode", nAPVs_U, -0.5, nAPVs_U-0.5, 1024, -0.5, 4095.5 );
    
    hrawADCs_by_APV_U = new TH2D( TString::Format( "hrawADCs_by_APV_U_%s", detname.Data() ), "U strip Full readout events; APV; raw ADC", nAPVs_U, -0.5, nAPVs_U-0.5, 1250,-500.,4500. );

    hrawADCs_by_APV_U_goodCM = new TH2D( TString::Format( "hrawADCs_by_APV_U_%s_goodCM", detname.Data() ), "U strip Full readout events; APV; raw ADC", nAPVs_U, -0.5, nAPVs_U-0.5, 1250,-500.,4500. );

    hrawADCs_by_APV_U_nopedsub = new TH2D( TString::Format( "hrawADCs_by_APV_U_nopedsub_%s", detname.Data() ), "U strip Full readout events; APV; raw ADC", nAPVs_U, -0.5, nAPVs_U-0.5, 1250,-500.,4500. );

    hrawADCs_by_APV_U_nopedsub_goodCM = new TH2D( TString::Format( "hrawADCs_by_APV_U_nopedsub_%s_goodCM", detname.Data() ), "U strip Full readout events; APV; raw ADC", nAPVs_U, -0.5, nAPVs_U-0.5, 1250,-500.,4500. );
    
    
    hADCs_by_APV_U = new TH2D( TString::Format( "hADCs_by_APV_U_%s", detname.Data() ), "U strip Full readout events; APV; corrected ADC", nAPVs_U, -0.5, nAPVs_U-0.5, 1250,-500.,4500. );
    
    hCM_OR_by_APV_U = new TH2D( TString::Format( "hCM_OR_by_APV_U_%s", detname.Data() ), "U strip all events; APV; CM out of range?", nAPVs_U, -0.5, nAPVs_U-0.5, 2, -0.5, 1.5 );


    hADCs_by_APV_U_goodCM = new TH2D( TString::Format( "hADCs_by_APV_U_%s_goodCM", detname.Data() ), "U strip Full readout events; APV; corrected ADC", nAPVs_U, -0.5, nAPVs_U-0.5, 1250,-500.,4500. );
   
    
    
    UInt_t nAPVs_V = fNstripsV/fN_APV25_CHAN;
    hcommonmode_mean_by_APV_V = new TH2D( TString::Format( "hCommonModeMean_by_APV_V_%s", detname.Data() ), "distribution of common-mode means for V strip pedestal data; APV card; Common-mode", nAPVs_V, -0.5, nAPVs_V-0.5, 1024, -0.5, 4095.5 );

    hcommonmode_online_by_APV_V = new TH2D( TString::Format( "hCommonModeOnline_by_APV_V_%s", detname.Data() ), "Online ZS events: common-mode means for V strips; APV card; online common-mode", nAPVs_V, -0.5, nAPVs_V-0.5, 1024, -0.5, 4095.5 );
    
    hrawADCs_by_APV_V = new TH2D( TString::Format( "hrawADCs_by_APV_V_%s", detname.Data() ), "V strip Full readout events; APV; raw ADC", nAPVs_V, -0.5, nAPVs_V-0.5, 1250,-500.,4500. );
    
    hrawADCs_by_APV_V_goodCM = new TH2D( TString::Format( "hrawADCs_by_APV_V_%s_goodCM", detname.Data() ), "V strip Full readout events; APV; raw ADC", nAPVs_V, -0.5, nAPVs_V-0.5, 1250,-500.,4500. );

    hrawADCs_by_APV_V_nopedsub = new TH2D( TString::Format( "hrawADCs_by_APV_V_nopedsub_%s", detname.Data() ), "V strip Full readout events; APV; raw ADC", nAPVs_V, -0.5, nAPVs_V-0.5, 1250,-500.,4500. );
    
    hrawADCs_by_APV_V_nopedsub_goodCM = new TH2D( TString::Format( "hrawADCs_by_APV_V_nopedsub_%s_goodCM", detname.Data() ), "V strip Full readout events; APV; raw ADC", nAPVs_V, -0.5, nAPVs_V-0.5, 1250,-500.,4500. );
    
    hADCs_by_APV_V = new TH2D( TString::Format( "hADCs_by_APV_V_%s", detname.Data() ), "V strip Full readout events; APV; corrected ADC", nAPVs_V, -0.5, nAPVs_V-0.5, 1250,-500.,4500. );

    hADCs_by_APV_V_goodCM = new TH2D( TString::Format( "hADCs_by_APV_V_%s_goodCM", detname.Data() ), "V strip Full readout events; APV; corrected ADC", nAPVs_V, -0.5, nAPVs_V-0.5, 1250,-500.,4500. );
    
    hCM_OR_by_APV_V = new TH2D( TString::Format( "hCM_OR_by_APV_V_%s", detname.Data() ), "V strip all events; APV; CM out of range?", nAPVs_V, -0.5, nAPVs_V-0.5, 2, -0.5, 1.5 );

    fPedDiagHistosInitialized = true;
  }

  if( fMakeCommonModePlots && !fCommonModePlotsInitialized ){

    //std::cout << "SBSGEMModule::Begin: making common-mode histograms for module " << GetName() << std::endl;
    //U strips:
   
    fCommonModeDistU = new TH2D( TString::Format( "hcommonmodeU_%s", detname.Data() ), "U/X strips common-mode (user); APV card; Common-mode - Common-mode mean (user)", fNAPVs_U, -0.5, fNAPVs_U-0.5, 500, -200.0,500.0 );

    fCommonModeDistU_Histo = new TH2D( TString::Format( "hcommonmodeU_histo_%s", detname.Data() ), "U/X strips common-mode (user); APV card; Common-mode - Common-mode mean (Histogram)", fNAPVs_U, -0.5, fNAPVs_U-0.5, 500, -200.0,500.0 );
    
    
    fCommonModeDistU_Sorting = new TH2D( TString::Format( "hcommonmodeU_sorting_%s", detname.Data() ), "U/X strips common-mode (Sorting); APV card; Common-mode (Sorting) - Common-mode mean", fNAPVs_U, -0.5, fNAPVs_U-0.5, 500, -200.0,500.0 );
    
    
    fCommonModeDistU_Danning = new TH2D( TString::Format( "hcommonmodeU_danning_%s", detname.Data() ), "U/X strips common-mode (Danning); APV card; Common-mode (Danning) - Common-mode mean", fNAPVs_U, -0.5, fNAPVs_U-0.5, 500, -200.0,500.0 );
    
   
    fCommonModeDiffU = new TH2D( TString::Format( "hcommonmodeU_diff_%s", detname.Data() ), "U/X strips (all events); APV card; Common-mode (User) - Common-mode (Danning online)", fNAPVs_U, -0.5, fNAPVs_U-0.5, 200, -200.0, 500.0 );

    fCommonModeDiffU_Uncorrected = new TH2D( TString::Format( "hcommonmodeU_diff_uncorrected_%s", detname.Data() ), "U/X strips (uncorrected events only); APV card; CM (User) - CM (Danning online)", fNAPVs_U, -0.5, fNAPVs_U-0.5, 200, -200.0, 500.0 );
    
    fCommonModeCorrectionU = new TH2D( TString::Format( "hcommonmodeU_corr_%s", detname.Data() ), "U/X strips; APV card; applied CM correction", fNAPVs_U, -0.5, fNAPVs_U-0.5, 200, -200.0, 200.0 );
    
    //fCommonModeNotCorrectionU = new TH2D( TString::Format( "hcommonmodeU_notcorr_%s", detname.Data() ), "U/X strips; APV card; CM Sorting - CM Danning without Correction", fNAPVs_U, -0.5, fNAPVs_U-0.5, 250, -100.0, 100.0 );
    fCommonModeResidualBiasU = new TH2D( TString::Format( "hcommonmodeU_bias_%s", detname.Data() ), "U/X strips; APV card; residual bias (corr. - true)", fNAPVs_U, -0.5, fNAPVs_U-0.5, 200, -200, 200 );

    fCommonModeResidualBiasU_corrected = new TH2D( TString::Format( "hcommonmodeU_bias_corrected_%s", detname.Data() ), "U/X strips; APV card; residual bias (corr. - true)", fNAPVs_U, -0.5, fNAPVs_U-0.5, 200, -200, 200 );

    fCommonModeResidualBias_vs_OccupancyU = new TH2D( TString::Format( "hcommonmodeU_bias_vs_occupancy_%s", detname.Data() ), "U/X strips; APV in-range occupancy; residual bias (corr. - true)", 100, 0.0, 1.0, 200, -200, 200 );

    //V strips:
    fCommonModeDistV = new TH2D( TString::Format( "hcommonmodeV_%s", detname.Data() ), "V/Y strips common-mode (user); APV card; Common-mode - Common-mode mean (user)", fNAPVs_V, -0.5, fNAPVs_V-0.5, 500, -200.0,500.0 );

    fCommonModeDistV_Histo = new TH2D( TString::Format( "hcommonmodeV_histo_%s", detname.Data() ), "V/Y strips common-mode (user); APV card; Common-mode - Common-mode mean (Histogram)", fNAPVs_V, -0.5, fNAPVs_V-0.5, 500, -200.0,500.0 );
    
    fCommonModeDistV_Sorting = new TH2D( TString::Format( "hcommonmodeV_sorting_%s", detname.Data() ), "V/Y strips common-mode (Sorting); APV card; Common-mode (Sorting) - Common-mode mean", fNAPVs_V, -0.5, fNAPVs_V-0.5, 500, -200.0,500.0 );
    
    
    fCommonModeDistV_Danning = new TH2D( TString::Format( "hcommonmodeV_danning_%s", detname.Data() ), "V/Y strips common-mode (Danning); APV card; Common-mode (Danning) - Common-mode mean", fNAPVs_V, -0.5, fNAPVs_V-0.5, 500, -200.0,500.0 );
    
    
    fCommonModeDiffV = new TH2D( TString::Format( "hcommonmodeV_diff_%s", detname.Data() ), "V/Y strips (all events); APV card; Common-mode (User) - Common-mode (Danning online)", fNAPVs_V, -0.5, fNAPVs_V-0.5, 200, -200.0,500.0 );

    fCommonModeDiffV_Uncorrected = new TH2D( TString::Format( "hcommonmodeV_diff_uncorrected_%s", detname.Data() ), "V/Y strips (uncorrected events only); APV card; CM (User) - CM (Danning online)", fNAPVs_V, -0.5, fNAPVs_V-0.5, 200, -200.0, 500.0 );
    
    fCommonModeCorrectionV = new TH2D( TString::Format( "hcommonmodeV_corr_%s", detname.Data() ), "V/Y strips; APV card; applied CM correction", fNAPVs_V, -0.5, fNAPVs_V-0.5, 200, -200.0, 200.0 );

    fCommonModeResidualBiasV = new TH2D( TString::Format( "hcommonmodeV_bias_%s", detname.Data() ), "V/Y strips; APV card; residual bias (corr. - true)", fNAPVs_V, -0.5, fNAPVs_V-0.5, 200, -200, 200 );

    fCommonModeResidualBiasV_corrected = new TH2D( TString::Format( "hcommonmodeV_bias_corrected_%s", detname.Data() ), "V/Y strips; APV card; residual bias (corr. - true)", fNAPVs_V, -0.5, fNAPVs_V-0.5, 200, -200, 200 );

    fCommonModeResidualBias_vs_OccupancyV = new TH2D( TString::Format( "hcommonmodeV_bias_vs_occupancy_%s", detname.Data() ), "V/Y strips; APV in-range occupancy; residual bias (corr. - true)", 100, 0.0, 1.0, 200, -200, 200 );
    
    //fCommonModeNotCorrectionV = new TH2D( TString::Format( "hcommonmodeV_notcorr_%s", detname.Data() ), "V/Y strips; APV card; CM Sorting - CM Danning without Correction", fNAPVs_V, -0.5, fNAPVs_V-0.5, 250, -100.0, 100.0 );

    // fCommonModeDistU->Print();
    // fCommonModeDistU_Sorting->Print();
    // fCommonModeDistU_Danning->Print();
    // fCommonModeDiffU->Print();

    // fCommonModeDistV->Print();
    // fCommonModeDistV_Sorting->Print();
    // fCommonModeDistV_Danning->Print();
    // fCommonModeDiffV->Print();
    fCommonModePlotsInitialized = true;
  }


  return fStatus = kOK;
}

Int_t   THcLADGEMModule::Begin( THaRunBase* r){ 
  fZeroSuppress = fZeroSuppress && !fPedestalMode; 
  fCODA_BUILD_ALL_SAMPLES = r->GetDAQInfo("VTP_MPDRO_BUILD_ALL_SAMPLES") == "1";
  fCODA_CM_ENABLED = r->GetDAQInfo("VTP_MPDRO_ENABLE_CM") == "1";

  if(r->GetDAQInfo("VTP_MPDRO_ENABLE_CM").empty()){
    fCODA_BUILD_ALL_SAMPLES = -1;
    fCODA_CM_ENABLED = -1;
  }
  std::cout << "THcLADGEMModule::Begin: CODA DAQ info: VTP_MPDRO_BUILD_ALL_SAMPLES = " << r->GetDAQInfo("VTP_MPDRO_BUILD_ALL_SAMPLES") 
      << ", VTP_MPDRO_ENABLE_CM = " << r->GetDAQInfo("VTP_MPDRO_ENABLE_CM") << std::endl;
  return 0;
}
//____________________________________________________________________________________
Int_t THcLADGEMModule::ReadDatabase(const TDatime &date) {
  //  cout << "THcLADGEMModule::ReadDatabase" << endl;

  // Define default values
  fZeroSuppress    = kTRUE;
  fZeroSuppressRMS = 5.0; // threshold in units of RMS:

  fmin_strip_per_clust = 1;
  fmax_strip_per_clust = -1;
  fMinTimeSamp        = 0;
  fMaxTimeSamp        = fN_MPD_TIME_SAMP - 1;
  fNegSignalStudy      = kFALSE;

  fPedestalMode                = kFALSE;
  fSubtractPedBeforeCommonMode = false; // only affects the pedestal-mode analysis
  fOnlineZeroSuppression       = kFALSE;

  fCommonModeFlag    = 1; //
  fCommonModeOnlFlag = 3; // 3 = Danning method during GMn, 4 = Danning method during GEn
  // Default: discard highest and lowest 28 strips for "sorting method" common-mode calculation:
  fCommonModeNstripRejectHigh = 28;
  fCommonModeNstripRejectLow  = 28;
  fCommonModeNumIterations    = 3;
  fCommonModeMinStripsInRange = 10;
  fMakeCommonModePlots        = false;
  fCommonModePlotsInitialized = false;
  fCommonModePlots_DBoverride = false;

  fPedSubFlag = 1; // default to online ped subtraction, as that is the mode we will run in most of the time

  fTrigTime              = 0.0;
  fMaxTrigTimeCorrection = 25.0;
  fTrigTimeSlope         = 1.0; // GEM time vs trig time correlation

  fSamplePeriod   = 24.0;   // nanoseconds:
  fSigma_hitshape = 0.0004; // 0.4 mm; controls cluster-splitting algorithm

  // Default clustering parameters:
  fThresholdSample     = 50.0;
  fThresholdStripSum   = 250.0;
  fThresholdClusterSum = 500.0;

  fThresholdSampleDeconv      = 50.0;
  fThresholdDeconvADCMaxCombo = 75.0;
  fThresholdClusterSumDeconv  = 150.0;

  fADCasymCut    = 0.8;
  fTimeCutUVdiff = 30.0;
  //  fCorrCoeffCut = -1.1;
  //  fCorrCoeffCutDeconv = -1.1;

  fADCasymSigma   = 0.06;
  fADCratioSigma  = 0.1;
  fTimeCutUVsigma = 3.0; // ns

  fTimeCutUVdiffDeconv  = 40.0; // ns
  fTimeCutUVsigmaDeconv = 7.0;

  fTimeCutUVdiffFit  = 30.0;
  fTimeCutUVsigmaFit = 3.0;

  // default pitch = 0.0004 for all modules offset = zero:
  fUStripPitch  = 0.0004;
  fVStripPitch  = 0.0004;
  fUStripOffset = 0.0;
  fVStripOffset = 0.0;

  fChan_CM_flags       = 640; // default to 640 (so as not to step on up to 40 MPDs per VTP crate):
  fChan_TimeStamp_low  = 641;
  fChan_TimeStamp_high = 642;
  fChan_MPD_EventCount = 643;
  fChan_MPD_Debug      = 644;

  UInt_t MAXNSAMP_PER_APV = fN_APV25_CHAN * fN_MPD_TIME_SAMP;
  // arrays to hold raw data from one APV card:
  fStripAPV.resize(MAXNSAMP_PER_APV);
  fRawStripAPV.resize(MAXNSAMP_PER_APV);
  fRawADC_APV.resize(MAXNSAMP_PER_APV);
  fPedSubADC_APV.resize(MAXNSAMP_PER_APV);
  fCommonModeSubtractedADC_APV.resize(MAXNSAMP_PER_APV);
  fCM_online.resize(fN_MPD_TIME_SAMP);

  fRMS_ConversionFactor = sqrt(fN_MPD_TIME_SAMP); //=2.45

  fAPVmapping = LADGEM::kUVA_XY; // default to UVA X/Y style APV mapping, but require this in the database::

  fModuleGain             = 1.0;
  fCommonModeRange_nsigma = 5.0;
  fSuppressFirstLast      = 0; // suppress strips peaking in first or last time sample by default:
  fStripTau           = 56.0;  // ns, default value. Eventually load this from DB. This is not actually used as of yet.
  fUseStripTimingCuts = 0;
  fUseTSchi2cut       = false;

  for (int axis = 0; axis < 2; axis++) {

    fStripMaxTcut_central[axis] = 87.0; // ns
    fStripMaxTcut_width[axis]   = 4.5;  // sigmas
    fStripMaxTcut_sigma[axis]   = 7.0;  // ns, for purpose of "hit quality chi2" calculation

    fStripMaxTcut_central_deconv[axis] = 50.0; // ns
    fStripMaxTcut_width_deconv[axis]   = 4.5;  // number of sigmas
    fStripMaxTcut_sigma_deconv[axis]   = 15.0; // ns

    fStripMaxTcut_central_fit[axis] = 20.0; // ns
    fStripMaxTcut_width_fit[axis]   = 4.5;  // sigmas
    fStripMaxTcut_sigma_fit[axis]   = 10.0; // ns

    fHitTimeMean[axis]        = 87.0;
    fHitTimeSigma[axis]       = 7.0;
    fHitTimeMeanDeconv[axis]  = 50.0;
    fHitTimeSigmaDeconv[axis] = 15.0;
    fHitTimeMeanFit[axis]     = 20.0;
    fHitTimeSigmaFit[axis]    = 10.0;
  }

  fSigmaHitTimeAverageCorrected = 5.0;  // ns
  fStripAddTcut_width           = 50.0; // this one we keep in ns
  fStripAddCorrCoeffCut         = 0.25;
  fStripTSchi2Cut               = 10.0; // not yet clear what is a good value for this.

  // fN_MPD_TIME_SAMP default value set to 6 at the beginning
  fGoodStrip_TSfrac_mean.resize(fN_MPD_TIME_SAMP);
  fGoodStrip_TSfrac_sigma.resize(fN_MPD_TIME_SAMP);

  // Define some defaults for these:
  double fracmean_default[6]  = {0.055, 0.135, 0.203, 0.224, 0.208, 0.174};
  double fracsigma_default[6] = {0.034, 0.039, 0.019, 0.021, 0.032, 0.035};

  for (int isamp = 0; isamp < fN_MPD_TIME_SAMP; isamp++) {
    if (isamp < fN_MPD_TIME_SAMP && isamp < 6) {
      fGoodStrip_TSfrac_mean[isamp]  = fracmean_default[isamp];
      fGoodStrip_TSfrac_sigma[isamp] = fracsigma_default[isamp];
    }
  }

  fMeasureCommonMode         = true;
  fNeventsCommonModeLookBack = 10;

  fCorrectCommonMode          = false;
  fCorrectCommonModeMinStrips = 20;
  fCorrectCommonMode_Nsigma   = 5.0;

  fCommonModeBinWidth_Nsigma  = 1.0; // Bin width +/- 1 sigma by default
  fCommonModeScanRange_Nsigma = 4.0; // Scan window +/- 4 sigma
  fCommonModeStepSize_Nsigma  = 0.2; // sigma/5 for step size:

  fClusteringFlag    = 0; //"standard" clustering based on sum of six time samples on a strip
  fDeconvolutionFlag = 0; // Default should be zero

  // Default cut parameters related clustering
  fThresh_2ndMax_nsigma   = 3.0;
  fThresh_2ndMax_fraction = 0.15;

  fMaxNeighborsU_totalcharge = 4;
  fMaxNeighborsV_totalcharge = 4;
  fMaxNeighborsU_hitpos      = 3;
  fMaxNeighborsV_hitpos      = 3;

  //  std::vector<Double_t> rawpedu,rawpedv;
  //  std::vector<Double_t> rawrmsu,rawrmsv;

  const DBRequest list[] = {{"lgem_chanmap_file", &fChanMapFileName, kString}, {"is_mc", &fIsMC, kInt, 0, 1}, {0}};
  gHcParms->LoadParmValues((DBRequest *)&list, "");
  if (fIsMC) {
    fChanMapFileName = "MAPS/LAD/DETEC/GEM/lgem_chan_mc.map"; // Todo: make this a parameter, but this works for now
  }

  string prefix = "lgem_";
  prefix += GetName(); // e.g. prefix: lgem_m0

  for (int i = 0; i < 3; i++)
    fCenter[i] = 0.; // init

  for (int i = 0; i < 2; i++)
    fGEMAngle[i] = 0.; // init
  for (int i = 0; i < 3; i++)
    fRotation[i] = 0.; // init
  

  // Geometry for each module
  const DBRequest list1[] = {{"_layer", &fLayer, kInt, 0, 1},
                             {"_apvmap", &fAPVmapping, kInt, 0, 1},
                             {"_nstripsu", &fNstripsU, kInt, 0, 1},
                             {"_nstripsv", &fNstripsV, kInt, 0, 1},
                             {"_uangle", &fUAngle, kDouble, 0, 1}, // mandatory: Angle of "U" strips wrt X axis
                             {"_vangle", &fVAngle, kDouble, 0, 1}, // mandatory: Angle of "V" strips wrt X axis
                             {"_position", fCenter, kDouble, static_cast<UInt_t>(3)},
                             {"_angle", fGEMAngle, kDouble, static_cast<UInt_t>(2)}, // angle between centerline to the beam line(rotate along y, then z)
                             {"_rotation", fRotation, kDouble, static_cast<UInt_t>(3)}, // rotation around the center of the detector along x, y, and z-axis
                             {0}};
  gHcParms->LoadParmValues((DBRequest *)&list1, prefix.c_str());
  // common for all layers -- move to THcLADGEM::ReadDatabase?
  const DBRequest list2[] = {{"lgem_zerosuppress", &fZeroSuppress, kInt, 0, 1}, // optional
                             {"lgem_zerosuppress_nsigma", &fZeroSuppressRMS, kDouble, 0, 1},
                             {"lgem_do_neg_signal_study", &fNegSignalStudy, kInt, 0, 1},
                             {"lgem_pedestal_mode", &fPedestalMode, kInt, 0, 1},
                             {"lgem_onlinezerosuppress", &fOnlineZeroSuppression, kInt, 0, 1},
                             {"lgem_commonmode_flag", &fCommonModeFlag, kInt, 0, 1},
                             {"lgem_commonmode_online_flag", &fCommonModeOnlFlag, kInt, 0, 1},
                             {"lgem_commonmode_nstriphi", &fCommonModeNstripRejectHigh, kInt, 0, 1},
                             {"lgem_commonmode_nstriplo", &fCommonModeNstripRejectLow, kInt, 0, 1},
                             {"lgem_commonmode_niter", &fCommonModeNumIterations, kInt, 0, 1},
                             {"lgem_commonmode_range_nsigma", &fCommonModeRange_nsigma, kDouble, 0, 1},
                             {"lgem_commonmode_minstrips", &fCommonModeMinStripsInRange, kInt, 0, 1},
                             {"lgem_cmplots_flag", &fMakeCommonModePlots, kInt, 0, 1},
                             {"lgem_pedsub_online", &fPedSubFlag, kInt, 0, 1},
                             {"lgem_max2Dhits", &fMAX2DHITS, kInt, 0, 1},
                             {"lgem_maxtrigtime_correction", &fMaxTrigTimeCorrection, kDouble, 0, 1},
                             {"lgem_trigtime_slope", &fTrigTimeSlope, kDouble, 0, 1},
                             {"lgem_sample_period", &fSamplePeriod, kDouble, 0, 1},
                             {"lgem_sigmahitshape", &fSigma_hitshape, kDouble, 0, 1},
                             {"lgem_threshold_sample", &fThresholdSample, kDouble, 0, 1},
                             {"lgem_threshold_stripsum", &fThresholdStripSum, kDouble, 0, 1},
                             {"lgem_threshold_clustersum", &fThresholdClusterSum, kDouble, 0, 1},
                             {"lgem_threshold_sample_deconv", &fThresholdSampleDeconv, kDouble, 0, 1},
                             {"lgem_threshold_maxcombo_deconv", &fThresholdDeconvADCMaxCombo, kDouble, 0, 1},
                             {"lgem_threshold_clustersum_deconv", &fThresholdClusterSumDeconv, kDouble, 0, 1},
                             {"lgem_adcasym_cut", &fADCasymCut, kDouble, 0, 1},
                             {"lgem_adcasym_sigma", &fADCasymSigma, kDouble, 0, 1},
                             {"lgem_deltat_cut", &fTimeCutUVdiff, kDouble, 0, 1},
                             {"lgem_deltat_sigma", &fTimeCutUVsigma, kDouble, 0, 1},
                             {"lgem_deltat_cut_deconv", &fTimeCutUVdiffDeconv, kDouble, 0, 1},
                             {"lgem_deltat_sigma_deconv", &fTimeCutUVsigmaDeconv, kDouble, 0, 1},
                             {"lgem_deltat_cut_fit", &fTimeCutUVdiffFit, kDouble, 0, 1},
                             {"lgem_deltat_sigma_fit", &fTimeCutUVsigmaFit, kDouble, 0, 1},
                             {"lgem_adcratio_sigma", &fADCratioSigma, kDouble, 0, 1},
                             {"lgem_upitch", &fUStripPitch, kDouble, 0, 1},
                             {"lgem_vpitch", &fVStripPitch, kDouble, 0, 1},
                             {"lgem_chan_cm_flags", &fChan_CM_flags, kInt, 0, 1},
                             {"lgem_chan_timestamp_low", &fChan_TimeStamp_low, kInt, 0, 1},
                             {"lgem_chan_timestamp_high", &fChan_TimeStamp_high, kInt, 0, 1},
                             {"lgem_chan_event_count", &fChan_MPD_EventCount, kInt, 0, 1},
                             {"lgem_chan_mpd_debug", &fChan_MPD_Debug, kInt, 0, 1},
                             {"lgem_modulegain", &fModuleGain, kDouble, 0, 1},
                             {"lgem_suppressfirstlast", &fSuppressFirstLast, kInt, 0, 1},
                             {"lgem_deconvolution_tau", &fStripTau, kDouble, 0, 1},
                             {"lgem_usestriptimingcut", &fUseStripTimingCuts, kInt, 0, 1},
                             {"lgem_useTSchi2cut", &fUseTSchi2cut, kInt, 0, 1},
                             {"lgem_sigma_tcorr", &fSigmaHitTimeAverageCorrected, kDouble, 0, 1},
                             {"lgem_addstrip_tcut", &fStripAddTcut_width, kDouble, 0, 1},
                             {"lgem_addstrip_ccor_cut", &fStripAddCorrCoeffCut, kDouble, 0, 1},
                             {"lgem_striptschi2_cut", &fStripTSchi2Cut, kDouble, 0, 1},
                             {"lgem_measure_common_mode", &fMeasureCommonMode, kInt, 0, 1},
                             {"lgem_commonmode_nevents_lookback", &fNeventsCommonModeLookBack, kUInt, 0, 1},
                             {"lgem_correct_common_mode", &fCorrectCommonMode, kInt, 0, 1},
                             {"lgem_correct_common_mode_minstrips", &fCorrectCommonModeMinStrips, kUInt, 0, 1},
                             {"lgem_correct_common_mode_nsigma", &fCorrectCommonMode_Nsigma, kDouble, 0, 1},
                             {"lgem_commonmode_binwidth_nsigma", &fCommonModeBinWidth_Nsigma, kDouble, 0, 1},
                             {"lgem_commonmode_scanrange_nsigma", &fCommonModeScanRange_Nsigma, kDouble, 0, 1},
                             {"lgem_commonmode_stepsize_nsigma", &fCommonModeStepSize_Nsigma, kDouble, 0, 1},
                             {"lgem_clustering_flag", &fClusteringFlag, kInt, 0, 1},
                             {"lgem_deconvolution_flag", &fDeconvolutionFlag, kInt, 0, 1},
                             {"lgem_peakprominence_minsigma", &fThresh_2ndMax_nsigma, kDouble, 0, 1},
                             {"lgem_peakprominence_minfraction", &fThresh_2ndMax_fraction, kDouble, 0, 1},
                             {"lgem_maxnu_charge", &fMaxNeighborsU_totalcharge, kInt, 0, 1},
                             {"lgem_maxnv_charge", &fMaxNeighborsV_totalcharge, kInt, 0, 1},
                             {"lgem_maxnu_pos", &fMaxNeighborsU_hitpos, kInt, 0, 1},
                             {"lgem_maxnv_pos", &fMaxNeighborsV_hitpos, kInt, 0, 1},
                             {"lgem_min_strip_per_clust", &fmin_strip_per_clust, kInt, 0, 1},
                             {"lgem_max_strip_per_clust", &fmax_strip_per_clust, kInt, 0, 1},
                             {"lgem_min_time_samp", &fMinTimeSamp, kInt, 0, 1},
                             {"lgem_max_time_samp", &fMaxTimeSamp, kInt, 0, 1},
                             {0}};
  gHcParms->LoadParmValues((DBRequest *)&list2, "");

  InitAPVMAP();

  // prevent the user from defining something silly for the common-mode stuff:
  fCommonModeNstripRejectLow  = std::min(50, std::max(0, fCommonModeNstripRejectLow));
  fCommonModeNstripRejectHigh = std::min(50, std::max(0, fCommonModeNstripRejectHigh));
  fCommonModeNumIterations    = std::min(10, std::max(2, fCommonModeNumIterations));
  fCommonModeMinStripsInRange = std::min(fN_APV25_CHAN - 25, std::max(1, fCommonModeMinStripsInRange));

  double x           = fSamplePeriod / fStripTau;
  fDeconv_weights[0] = exp(x - 1.0) / x;     //~1.32
  fDeconv_weights[1] = -2.0 * exp(-1.0) / x; //~ -1.72
  fDeconv_weights[2] = exp(-1.0 - x) / x;    // 0.56

  if (fIsMC) {
    fCommonModeFlag        = -1;
    fPedestalMode          = false;
    fOnlineZeroSuppression = true;
    fAPVmapping            = LADGEM::kMC;
  }

  fPxU = cos(fUAngle * TMath::DegToRad());
  fPyU = sin(fUAngle * TMath::DegToRad());
  fPxV = cos(fVAngle * TMath::DegToRad());
  fPyV = sin(fVAngle * TMath::DegToRad());

  // FIXME: THcParmList has no support for vector type input
  // Read channel map data from a separate functiion using Database::LoadDatabase from analyzer for now
  // This should be updated by adding the support for future use
  GetChannelMap(prefix.c_str(), date);

  fNAPVs_U = std::max(fNAPVs_U, fNstripsU / fN_APV25_CHAN);
  fNAPVs_V = std::max(fNAPVs_V, fNstripsV / fN_APV25_CHAN);

  // Initialize all pedestals to zero, RMS values to default:
  fPedestalU.clear();
  fPedestalU.resize(fNstripsU);

  fPedRMSU.clear();
  fPedRMSU.resize(fNstripsU);

  for (UInt_t istrip = 0; istrip < fNstripsU; istrip++) {
    fPedestalU[istrip] = 0.0;
    fPedRMSU[istrip]   = 10.0; // placeholder to be replaced by value from database
    // Add Reading Pedestal values
  }

  fPedestalV.clear();
  fPedestalV.resize(fNstripsV);

  fPedRMSV.clear();
  fPedRMSV.resize(fNstripsV);

  for (UInt_t istrip = 0; istrip < fNstripsV; istrip++) {
    fPedestalV[istrip] = 0.0;
    fPedRMSV[istrip]   = 10.0;
  }

  // //resize all the "decoded strip" arrays to their maximum possible values for this module:
  UInt_t nstripsmax = fNstripsU + fNstripsV;

  fStrip.resize(nstripsmax);
  fAxis.resize(nstripsmax);
  fADCsamples.resize(nstripsmax);
  fRawADCsamples.resize(nstripsmax);
  fADCsamples_deconv.resize(nstripsmax);

  // The lines below are problematic and unnecessary
  for (unsigned int istrip = 0; istrip < nstripsmax; istrip++) {
    fADCsamples[istrip].resize(fN_MPD_TIME_SAMP);
    fRawADCsamples[istrip].resize(fN_MPD_TIME_SAMP);
    fADCsamples_deconv[istrip].resize(fN_MPD_TIME_SAMP);
  }

  fADCsums.resize(nstripsmax);
  fADCsumsDeconv.resize(nstripsmax);
  fStripADCavg.resize(nstripsmax);
  fStripIsU.resize(nstripsmax);
  fStripIsV.resize(nstripsmax);
  fStripOnTrack.resize(nstripsmax);
  fStripIsNeg.resize(nstripsmax);
  fStripIsNegU.resize(nstripsmax);
  fStripIsNegV.resize(nstripsmax);
  fStripIsNegOnTrack.resize(nstripsmax);
  fStripIsNegOnTrackU.resize(nstripsmax);
  fStripIsNegOnTrackV.resize(nstripsmax);
  fStripRaw.resize(nstripsmax);
  fStripEvent.resize(nstripsmax);
  fStripCrate.resize(nstripsmax);
  fStripMPD.resize(nstripsmax);
  fStripADC_ID.resize(nstripsmax);
  fStripTrackIndex.resize(nstripsmax);
  fKeepStrip.resize(nstripsmax);
  fMaxSamp.resize(nstripsmax);
  fMaxSampDeconv.resize(nstripsmax);
  fMaxSampDeconvCombo.resize(nstripsmax);
  fADCmax.resize(nstripsmax);
  fADCmaxDeconv.resize(nstripsmax);
  fADCmaxDeconvCombo.resize(nstripsmax);
  fTmean.resize(nstripsmax);
  fTmeanDeconv.resize(nstripsmax);
  fTsigma.resize(nstripsmax);
  fStripTdiff.resize(nstripsmax);
  fStripTSchi2.resize(nstripsmax);
  fStripTSprob.resize(nstripsmax);
  fStripCorrCoeff.resize(nstripsmax);
  fStripTfit.resize(nstripsmax);
  fTcorr.resize(nstripsmax);
  // Storing these by individual strip is redundant but convenient:
  fStrip_ENABLE_CM.resize(nstripsmax);
  fStrip_CM_GOOD.resize(nstripsmax);
  fStrip_BUILD_ALL_SAMPLES.resize(nstripsmax);

  fADCsamples1D.resize(nstripsmax * fN_MPD_TIME_SAMP);
  fRawADCsamples1D.resize(nstripsmax * fN_MPD_TIME_SAMP);
  fADCsamplesDeconv1D.resize(nstripsmax * fN_MPD_TIME_SAMP);

  // default all common-mode mean and RMS values to 0 and 10 respectively if they were
  //  NOT loaded from the DB and/or they are loaded with the wrong size:
  if (fCommonModeMeanU.size() != fNAPVs_U) {
    fCommonModeMeanU.resize(fNAPVs_U);
    for (unsigned int iAPV = 0; iAPV < fNAPVs_U; iAPV++) {
      fCommonModeMeanU[iAPV] = 0.0;
    }
  }

  if (fCommonModeRMSU.size() != fNAPVs_U) {
    fCommonModeRMSU.resize(fNAPVs_U);
    for (unsigned int iAPV = 0; iAPV < fNAPVs_U; iAPV++) {
      fCommonModeRMSU[iAPV] = 10.0;
    }
  }

  // default all common-mode mean and RMS values to 0 and 10 respectively if they were
  //  NOT loaded from the DB and/or they were loaded with the wrong size:
  if (fCommonModeMeanV.size() != fNAPVs_V) {
    fCommonModeMeanV.resize(fNAPVs_V);
    for (unsigned int iAPV = 0; iAPV < fNAPVs_V; iAPV++) {
      fCommonModeMeanV[iAPV] = 0.0;
    }
  }

  if (fCommonModeRMSV.size() != fNAPVs_V) {
    fCommonModeRMSV.resize(fNAPVs_V);
    for (unsigned int iAPV = 0; iAPV < fNAPVs_V; iAPV++) {
      fCommonModeRMSV[iAPV] = 10.0;
    }
  }

  // Initialize default "CM correction bias" values to zero if they were not loaded from the DB:
  if (fCMbiasU.size() != fNAPVs_U) {
    fCMbiasU.resize(fNAPVs_U);
    for (unsigned int iAPV = 0; iAPV < fNAPVs_U; iAPV++) {
      fCMbiasU[iAPV] = 0.0;
    }
  }

  if (fCMbiasV.size() != fNAPVs_V) {
    fCMbiasV.resize(fNAPVs_V);
    for (unsigned int iAPV = 0; iAPV < fNAPVs_V; iAPV++) {
      fCMbiasV[iAPV] = 0.0;
    }
  }

  // default all gains to 1 if they were not loaded from the DB and/or if they were loaded with the
  // wrong size:
  if (fUgain.size() != fNAPVs_U) {
    fUgain.resize(fNAPVs_U);
    for (unsigned int iAPV = 0; iAPV < fNAPVs_U; iAPV++) {
      fUgain[iAPV] = 1.0;
    }
  }

  // Multiply in Module gain:
  for (unsigned int iAPV = 0; iAPV < fNAPVs_U; iAPV++) {
    fUgain[iAPV] *= fModuleGain;
  }

  if (fVgain.size() != fNAPVs_V) {
    fVgain.resize(fNAPVs_V);
    for (unsigned int iAPV = 0; iAPV < fNAPVs_V; iAPV++) {
      fVgain[iAPV] = 1.0;
    }
  }

  for (unsigned int iAPV = 0; iAPV < fNAPVs_V; iAPV++) {
    fVgain[iAPV] *= fModuleGain;
  }
  std::cout<< "THcLADGEMModule: fPedestalMode = " << fPedestalMode << std::endl;
  if (fPedestalMode) {
    fZeroSuppress          = false;
    fOnlineZeroSuppression = false;
    // fPedSubFlag = 0;
    std::cout << "THcLADGEMModule: Pedestal mode ON, zero suppression OFF" << std::endl;
  }

  // Hard coded. Should fix later
  // axis: 0=x, 1=y
  int n_badstrips        = 0;
  int bad_pos_lo[]       = {50, 0, -40, 30, 27};
  int bad_pos_hi[]       = {55, 5, -33, 35, 30};
  int bad_strips_modid[] = {0, 0, 1, 1, 1};
  int bad_strips_axis[]  = {0, 1, 0, 0, 1};
  fBadStrips.resize(n_badstrips);
  for (int i = 0; i < n_badstrips; i++) {
    fBadStrips[i].pos_lo    = bad_pos_lo[i];
    fBadStrips[i].pos_hi    = bad_pos_hi[i];
    fBadStrips[i].module_id = bad_strips_modid[i];
    fBadStrips[i].axis      = bad_strips_axis[i];
  }

  return 0;
}

//____________________________________________________________________________________
Int_t THcLADGEMModule::GetChannelMap(const char *prefix, const TDatime &date) {

  // Keept it separate for now from ReadDatabase
  // as a reminder to merge the format of channel map into Hall C type parameters
  // or update parsing class to inclue vector type -- SP FIXME

  Int_t status;

  static const char *const here = "THcLADGEMModule::ReadDatabase";
  FILE *file                    = Podd::OpenDBFile(fChanMapFileName.c_str(), date, here, "r", 0);
  if (!file)
    return kFileError;

  const DBRequest request[] = {{"_chanmap", &fChanMapData, kIntV, 0, 0, 0}, {0}};

  status = LoadDB(file, date, request, prefix, 1);
  if (status != 0) {
    cout << "THcLADGEMModule::GetChannelMap Unable to read channel map data" << endl;
    fclose(file);
    return status;
  }

  fNAPVs_U = 0;
  fNAPVs_V = 0;

  // Read MPD channel map
  fMPDmap.clear();
  Int_t nentry = fChanMapData.size() / fMPDMAP_ROW_SIZE;

  fCommonModeResultContainer_by_APV.resize(nentry);
  fCommonModeRollingAverage_by_APV.resize(nentry);
  fCommonModeRollingRMS_by_APV.resize(nentry);
  fNeventsRollingAverage_by_APV.resize(nentry);

  fCMbiasResultContainer_by_APV.resize(nentry);
  fCommonModeOnlineBiasRollingAverage_by_APV.resize(nentry);
  fCommonModeOnlineBiasRollingRMS_by_APV.resize(nentry);
  fNeventsOnlineBias_by_APV.resize(nentry);

  for (Int_t mapline = 0; mapline < nentry; mapline++) {
    mpdmap_t thisdata;

    thisdata.crate  = fChanMapData[0 + mapline * fMPDMAP_ROW_SIZE];
    thisdata.slot   = fChanMapData[1 + mapline * fMPDMAP_ROW_SIZE];
    thisdata.mpd_id = fChanMapData[2 + mapline * fMPDMAP_ROW_SIZE];
    thisdata.gem_id = fChanMapData[3 + mapline * fMPDMAP_ROW_SIZE];
    thisdata.adc_id = fChanMapData[4 + mapline * fMPDMAP_ROW_SIZE];
    thisdata.i2c    = fChanMapData[5 + mapline * fMPDMAP_ROW_SIZE];
    thisdata.pos    = fChanMapData[6 + mapline * fMPDMAP_ROW_SIZE];
    thisdata.invert = fChanMapData[7 + mapline * fMPDMAP_ROW_SIZE];
    thisdata.axis   = fChanMapData[8 + mapline * fMPDMAP_ROW_SIZE];
    thisdata.index  = mapline;

    if (thisdata.axis == LADGEM::kUaxis) {
      fNAPVs_U++;
    } else {
      fNAPVs_V++;
    }

    fMPDmap.push_back(thisdata);

    fEventCount_by_APV.push_back(0);
    fT0_by_APV.push_back(0);
    fTcoarse_by_APV.push_back(0);
    fTfine_by_APV.push_back(0);
    fTimeStamp_ns_by_APV.push_back(0);

    // CommonMode related stuff
    fCommonModeResultContainer_by_APV[mapline].resize(fNeventsCommonModeLookBack * fN_MPD_TIME_SAMP);
    fCommonModeRollingAverage_by_APV[mapline] = 0.0;
    fCommonModeRollingRMS_by_APV[mapline]     = 10.0;
    fNeventsRollingAverage_by_APV[mapline]    = 0; // Really will be the number of time samples = 6 * number of events

    fCMbiasResultContainer_by_APV[mapline].resize(fNeventsCommonModeLookBack * fN_MPD_TIME_SAMP);
    fCommonModeOnlineBiasRollingAverage_by_APV[mapline] = 0.0;
    fCommonModeOnlineBiasRollingRMS_by_APV[mapline]     = 10.0;
    fNeventsOnlineBias_by_APV[mapline]                  = 0;
  } // Loop over mpd map data

  return 0;
}

//____________________________________________________________________________________
Int_t THcLADGEMModule::DefineVariables(EMode mode) {
  //  cout << "THcLADGEMModule::DefineVariables" << endl;

  if (mode == kDefine && fIsSetup)
    return kOK;
  fIsSetup = (mode == kDefine);

  // Strip level variables
  VarDef vars[] = {
      {"strip.nstripsfired", "Number of strips fired", kUInt, 0, &fNstrips_hit},
      {"strip.nstripsfired_pos", "Number of strips fired pos", kUInt, 0, &fNstrips_hit_pos},
      {"strip.nstripsfired_neg", "Number of strips fired neg", kUInt, 0, &fNstrips_hit_neg},
      {"strip.crate", "strip crate number", kUInt, 0, &(fStripCrate[0]), &fNstrips_hit},
      {"strip.mpd", "strip mpd number", kUInt, 0, &(fStripMPD[0]), &fNstrips_hit},
      {"strip.adc_id", "strip adc channel number", kUInt, 0, &(fStripADC_ID[0]), &fNstrips_hit},
      {"strip.istrip", "strip index", kUInt, 0, &(fStrip[0]), &fNstrips_hit},
      {"strip.IsU", "U strip?", kUInt, 0, &(fStripIsU[0]), &fNstrips_hit},
      {"strip.IsV", "V strip?", kUInt, 0, &(fStripIsV[0]), &fNstrips_hit},
      {"strip.ADCsamples", "ADC samples (index = isamp+Nsamples*istrip)", kDouble, 0, &(fADCsamples1D[0]),
       &fNdecoded_ADCsamples},
      {"strip.rawADCsamples", "raw ADC samples (no baseline subtraction)", kInt, 0, &(fRawADCsamples1D[0]),
       &fNdecoded_ADCsamples},
      {"strip.DeconvADCsamples", "Deconvoluted ADC samples (index = isamp+Nsamples*istrip)", kDouble, 0,
       &(fADCsamplesDeconv1D[0]), &fNdecoded_ADCsamples},
      {"strip.ADCsum", "Sum of ADC samples on a strip", kDouble, 0, &(fADCsums[0]), &fNstrips_hit},
      {"strip.DeconvADCsum", "Sum of deconvoluted ADC samples on a strip", kDouble, 0, &(fADCsumsDeconv[0]),
       &fNstrips_hit},
      {"strip.isampmax", "sample in which max ADC occurred on a strip", kUInt, 0, &(fMaxSamp[0]), &fNstrips_hit},
      {"strip.isampmaxDeconv", "sample in which max deconvoluted ADC occurred", kUInt, 0, &(fMaxSampDeconv[0]),
       &fNstrips_hit},
      {"strip.isampmaxDeconvCombo", "first of max. pair of deconvoluted samples", kUInt, 0, &(fMaxSampDeconvCombo[0]),
       &fNstrips_hit},
      {"strip.ADCmax", "Value of max ADC sample on a strip", kDouble, 0, &(fADCmax[0]), &fNstrips_hit},
      {"strip.DeconvADCmax", "Value of max deconvoluted ADC sample", kDouble, 0, &(fADCmaxDeconv[0]), &fNstrips_hit},
      {"strip.DeconvADCmaxCombo", "max sum of two adjacent deconv. samples", kDouble, 0, &(fADCmaxDeconvCombo[0]),
       &fNstrips_hit},
      {"strip.Tmean", "ADC-weighted mean strip time", kDouble, 0, &(fTmean[0]), &fNstrips_hit},
      {"strip.Tsigma", "ADC-weighted rms strip time", kDouble, 0, &(fTsigma[0]), &fNstrips_hit},
      {"strip.TmeanDeconv", "ADC-weighted mean deconvoluted strip time", kDouble, 0, &(fTmeanDeconv[0]), &fNstrips_hit},
      {"strip.Tcorr", "Corrected strip time", kDouble, 0, &(fTcorr[0]), &fNstrips_hit},
      {"strip.Tfit", "Fitted strip time", kDouble, 0, &(fStripTfit[0]), &fNstrips_hit},
      {"strip.Tdiff", "time diff. wrt max strip in cluster (or perhaps cluster tmean)", kDouble, 0, &(fStripTdiff[0]),
       &fNstrips_hit},
      {"strip.TSchi2", "chi2 of strip pulse shape (time samples) wrt average good strip pulse shape", kDouble, 0,
       &(fStripTSchi2[0]), &fNstrips_hit},
      {"strip.TSprob", "p-Value wrt average good strip pulse shape", kDouble, 0, &(fStripTSprob[0]), &fNstrips_hit},
      {"strip.CorrCoeff", "Correlation coefficient of strip wrt max strip on cluster (or perhaps cluster tmean)",
       kDouble, 0, &(fStripCorrCoeff[0]), &fNstrips_hit},
      {"strip.ADCavg", "average of ADC samples on a strip", kDouble, 0, &(fStripADCavg[0]), &fNstrips_hit},
      {nullptr}};

  DefineVarsFromList(vars, mode);

  // Add variable for number of 2D hits
  VarDef vars2D[] = {{"n2Dhits", "Number of 2D hits", kUInt, 0, &fN2Dhits}, {nullptr}};

  return DefineVarsFromList(vars2D, mode);
}

//____________________________________________________________________________________
Int_t THcLADGEMModule::Decode(const THaEvData &evdata) {

  // Hit counter
  fNstrips_hit     = 0;
  fNstrips_hit_neg = 0;
  fNstrips_hit_pos = 0;

  fNstrips_hitU     = 0;
  fNstrips_hitV     = 0;
  fNstrips_hitU_neg = 0;
  fNstrips_hitV_neg = 0;

  vector<UInt_t> &Strip                     = fStripAPV;
  vector<UInt_t> &rawStrip                  = fRawStripAPV;
  vector<Int_t> &rawADC                     = fRawADC_APV;
  vector<Double_t> &pedsubADC               = fPedSubADC_APV; // ped-subtracted, not necessarily common-mode subtracted
  vector<Double_t> &commonModeSubtractedADC = fCommonModeSubtractedADC_APV;

  int apvcounter  = 0;
  bool firstevent = true;
  //  UInt_t FirstEvCnt = 0;

  for (std::vector<mpdmap_t>::iterator it = fMPDmap.begin(); it != fMPDmap.end(); it++) {
    Int_t effChan = it->mpd_id << 4 | it->adc_id; // left-shift mpd id by 4 bits and take the bitwise OR with ADC_id to
                                                  // uniquely identify the APV card.

    // FIXME: take it out of the loop. take care of it later -- SP
    UInt_t nhits_timestamp_low  = evdata.GetNumHits(it->crate, it->slot, fChan_TimeStamp_low);
    UInt_t nhits_timestamp_high = evdata.GetNumHits(it->crate, it->slot, fChan_TimeStamp_high);
    UInt_t nhits_event_count    = evdata.GetNumHits(it->crate, it->slot, fChan_MPD_EventCount);

    if (nhits_timestamp_low > 0 && nhits_timestamp_high == nhits_timestamp_low &&
        nhits_event_count == nhits_timestamp_low) {
      for (unsigned int ihit = 0; ihit < nhits_timestamp_low; ihit++) {
        unsigned int fiber = evdata.GetRawData(it->crate, it->slot, fChan_TimeStamp_low, ihit);
        if (fiber == it->mpd_id) { // this is the channel we want:
          UInt_t Tlow  = evdata.GetData(it->crate, it->slot, fChan_TimeStamp_low, ihit);
          UInt_t Thigh = evdata.GetData(it->crate, it->slot, fChan_TimeStamp_high, ihit);
          UInt_t EvCnt = evdata.GetData(it->crate, it->slot, fChan_MPD_EventCount, ihit);

          if (firstevent) {
            // FirstEvCnt = EvCnt;
            firstevent = false;
          }

          fEventCount_by_APV[apvcounter] = EvCnt;

          // Fine time stamp is in the first 8 bits of Tlow;
          fTfine_by_APV[apvcounter] = Tlow & 0x000000FF;

          /*
          if( fMakeEventInfoPlots && fEventInfoPlotsInitialized ){
            hMPD_EventCount_Alignment->Fill( EvCnt - FirstEvCnt );
            hMPD_EventCount_Alignment_vs_Fiber->Fill( fiber, EvCnt - FirstEvCnt );

            hMPD_FineTimeStamp_vs_Fiber->Fill( fiber, fTfine_by_APV[apvcounter] * 4.0 );
          }
          */

          Long64_t Tcoarse = Thigh << 16 | (Tlow << 8);
          double Tc        = double(Tcoarse);

          if (EvCnt == 0)
            fT0_by_APV[apvcounter] = Tc;

          // T ref is the coarse time stamp of the reference APV (the first one, in this case)
          if (apvcounter == 0)
            fTref_coarse = Tc - fT0_by_APV[apvcounter];

          // This SHOULD make fTcoarse_by_APV the Tcoarse RELATIVE to the
          // "reference" APV
          fTcoarse_by_APV[apvcounter] = Tc - fT0_by_APV[apvcounter] - fTref_coarse;

          // FIXME: hard-coded constants from SBS-offline
          // We probably don't want to hard-code 24 ns and 4 ns here for the units of
          // Tcoarse and Tfine, but this should be fine for initial checkout of decoding:
          fTimeStamp_ns_by_APV[apvcounter] = 24.0 * fTcoarse_by_APV[apvcounter] + 4.0 * (fTfine_by_APV[apvcounter] % 6);

          break;
        }
      }
    }

    // Common-mode flags
    Bool_t CM_ENABLED        = fCommonModeFlag != 0 && fCommonModeFlag != 1 && !fPedestalMode;
    Bool_t BUILD_ALL_SAMPLES = !fOnlineZeroSuppression;
    Bool_t CM_OUT_OF_RANGE   = false;

    // These are supposed to be obtained from DAQ info
    // Leave it disabled for now until we know the parameter names
    if (fCODA_BUILD_ALL_SAMPLES != -1) {
      BUILD_ALL_SAMPLES = fCODA_BUILD_ALL_SAMPLES;
      fPedSubFlag       = (fCODA_BUILD_ALL_SAMPLES == 0);
    }
    if (fCODA_CM_ENABLED != -1)
      CM_ENABLED = fCODA_CM_ENABLED;

    UInt_t cm_flags      = 4 * CM_OUT_OF_RANGE + 2 * CM_ENABLED + BUILD_ALL_SAMPLES;
    UInt_t nhits_cm_flag = evdata.GetNumHits(it->crate, it->slot, fChan_CM_flags);

    if (nhits_cm_flag > 0) {

      // If applicable, find the common-mode/zero-suppression settings loaded from the raw data for this APV:
      // In principle in the SSP/VTP event format, there should be exactly one "hit" per APV in this "channel":
      for (unsigned int ihit = 0; ihit < nhits_cm_flag; ihit++) {
        int chan_temp = evdata.GetRawData(it->crate, it->slot, fChan_CM_flags, ihit);
        if (chan_temp == effChan) { // assume that this is only filled once per MPD per event, and exit the loop when we
                                    // find this MPD:
          // std::cout << "Before decoding cm flags, CM_ENABLED, BUILD_ALL_SAMPLES = " << CM_ENABLED << ", "
          // 	    << BUILD_ALL_SAMPLES << std::endl;
          cm_flags = evdata.GetData(it->crate, it->slot, fChan_CM_flags, ihit);
          // cm_flags_found = true;
          break;
        }
      }
    }

    // The proper logic of common-mode calculation/subtraction and zero suppression is as follows:
    //  1. If CM_ENABLED is true, we never calculate the common-mode ourselves, it has already been subtracted from the
    //  data:
    //  2. If BUILD_ALL_SAMPLES is false, then online zero suppression is enabled. We can, in addition, apply our own
    //  higher thresholds if we want:
    //  3. If CM_ENABLED is true, the pedestal has also been subtracted, so we don't subtract it again.
    //  4. If CM_ENABLED is false, we need to subtract the pedestals (maybe) AND calculate and subtract the common-mode:
    //  5. If BUILD_ALL_SAMPLES is false then CM_ENABLED had better be true!
    //  6. If CM_OUT_OF_RANGE is true then BUILD_ALL_SAMPLES must be true and CM_ENABLED
    //     must be false!

    CM_OUT_OF_RANGE   = cm_flags / 4;
    CM_ENABLED        = cm_flags / 2;
    BUILD_ALL_SAMPLES = cm_flags % 2;

    if (!BUILD_ALL_SAMPLES && !CM_ENABLED) { // This should never happen: skip this APV card
      continue;
    }
    if (CM_OUT_OF_RANGE && !BUILD_ALL_SAMPLES) { // this should also never happen: skip this APV card:
      continue;
    }

    if (CM_OUT_OF_RANGE && CM_ENABLED) { // force CM_ENABLED to false if CM_OUT_OF_RANGE is true;
      // This should have already been done online:
      CM_ENABLED = false;
    }

    if (fIsMC) {
      CM_ENABLED        = true;
      BUILD_ALL_SAMPLES = false;
    }

    // Decode MPD debug header
    UInt_t nhits_MPD_debug = 0;

    UInt_t CMcalc[fN_MPD_TIME_SAMP];
    Int_t CMcalc_signed[fN_MPD_TIME_SAMP];

    if (CM_ENABLED) { // try to decode MPD debug headers and see if the results make any sense:
      nhits_MPD_debug = evdata.GetNumHits(it->crate, it->slot, fChan_MPD_Debug);

      if (nhits_MPD_debug > 0) { // we expect to get three words per APV card:
        UInt_t wcount = 0;

        UInt_t MPDdebugwords[3];

        for (unsigned int ihit = 0; ihit < nhits_MPD_debug; ihit++) {
          int chan_temp = evdata.GetRawData(it->crate, it->slot, fChan_MPD_Debug, ihit);
          int word_temp = evdata.GetData(it->crate, it->slot, fChan_MPD_Debug, ihit);
          if (chan_temp == effChan && wcount < 3) {
            MPDdebugwords[wcount++] = word_temp;
          }
          if (wcount == 3)
            break; // if we found all 3 MPD debug words for this channel, exit the loop
        }
        if (wcount == 3) { // Then let's decode the debug headers:

          for (unsigned int iw = 0; iw < 3; iw++) {
            CMcalc[2 * iw] = (MPDdebugwords[iw] & 0xFFF) | ((MPDdebugwords[iw] & 0x1000) ? 0xFFFFF000 : 0x0);
            CMcalc[2 * iw + 1] =
                ((MPDdebugwords[iw] >> 13) & 0xFFF) | (((MPDdebugwords[iw] >> 13) & 0x1000) ? 0xFFFFF000 : 0x0);
            CMcalc_signed[2 * iw]     = Int_t(CMcalc[2 * iw]);
            CMcalc_signed[2 * iw + 1] = Int_t(CMcalc[2 * iw + 1]);

            fCM_online[2 * iw]     = double(CMcalc_signed[2 * iw]);
            fCM_online[2 * iw + 1] = double(CMcalc_signed[2 * iw + 1]);

            //Thesee are filled six times per APV per event:
	          if( it->axis == LADGEM::kUaxis ){
	            hcommonmode_online_by_APV_U->Fill( it->pos, fCM_online[2*iw] );
	            hcommonmode_online_by_APV_U->Fill( it->pos, fCM_online[2*iw+1] );
	          } else {
	            hcommonmode_online_by_APV_V->Fill( it->pos, fCM_online[2*iw] );
	            hcommonmode_online_by_APV_V->Fill( it->pos, fCM_online[2*iw+1] );
	          }
          }
        }
      }
    } // CM_ENABLED

    LADGEM::GEMaxis_t axis = it->axis == 0 ? LADGEM::kUaxis : LADGEM::kVaxis;

    Int_t nsamp = evdata.GetNumHits(it->crate, it->slot, effChan);
    if (nsamp > 0) {
      // Temporary variable to store the number of hits above negative saturation threshold
      // by time sample. If we have a full readout event and we fail to calculate a good
      // common-mode value in any time sample, we might want to ignore that APV's data? 
      
      //These are filled once per APV card per event
      // as long as the APV card has non-zero hits:
      int flag_cm_or = CM_OUT_OF_RANGE ? 1 : 0;
      if (fPedDiagHistosInitialized){
        if( axis == LADGEM::kUaxis ){
      	  hCM_OR_by_APV_U->Fill( it->pos, flag_cm_or );
        } else {
	        hCM_OR_by_APV_V->Fill( it->pos, flag_cm_or );
        }
      }
      Int_t nstrips = nsamp / fN_MPD_TIME_SAMP; // number of strips fired on this APV card (should be exactly 128 if
                                                // online zero suppression is NOT used):
      bool fullreadout = !CM_ENABLED && BUILD_ALL_SAMPLES && nstrips == fN_APV25_CHAN;

      double commonMode[fN_MPD_TIME_SAMP];
      double CommonModeCorrection[fN_MPD_TIME_SAMP]; // possible correction to apply, initialize to zero:
      for (int isamp = 0; isamp < fN_MPD_TIME_SAMP; isamp++) {
        commonMode[isamp]           = 0.0;
        CommonModeCorrection[isamp] = 0.0;
      }

      // First loop over the hits: populate strip, raw strip, raw ADC, ped sub ADC and common-mode-subtracted aDC:
      for (int iraw = 0; iraw < nsamp; iraw++) { // NOTE: iraw = isamp + fN_MPD_TIME_SAMP * istrip
        int strip             = evdata.GetRawData(it->crate, it->slot, effChan, iraw);
        UInt_t decoded_rawADC = evdata.GetData(it->crate, it->slot, effChan, iraw);

        int isamp = iraw % fN_MPD_TIME_SAMP;

        Int_t ADC = Int_t(decoded_rawADC);

        rawStrip[iraw] = strip;
        Strip[iraw]    = GetStripNumber(strip, it->axis, it->pos, it->invert);
        rawADC[iraw]   = ADC;

        // Note: this prints out all 128 channels for 6 samples for each APV
        // cout << "nsamp, iraw, rawstrip, strip, ADC: " << nsamp << " " << iraw << " "<< strip << " " << Strip[iraw] <<
        // " " << ADC << endl;

        double ped = (axis == LADGEM::kUaxis) ? fPedestalU[Strip[iraw]] : fPedestalV[Strip[iraw]];

        // If pedestal subtraction was done online, don't do it again:
        // In pedestal mode, the DAQ should NOT have subtracted the pedestals,
        // but even if it did, it shouldn't affect the pedestal analysis to first order whether we subtract the
        // pedestals or not:
        if (fPedSubFlag != 0 && !fPedestalMode && !fIsMC)
          ped = 0.0;
        if (CM_ENABLED && !fIsMC) {
          ped = 0.0; // If this is true then the pedestal was DEFINITELY always calculated online:
          // rawADC[iraw] += fCM_online[isamp] (not yet sure if we want to add back the online-CM) ;
        }

        pedsubADC[iraw]               = double(ADC) - ped;
        commonModeSubtractedADC[iraw] = double(ADC) - ped;

        if (fPedestalMode) {
          // do simple common-mode calculation involving the simple average of all 128 (ped-subtracted) ADC
          // values

          if (fSubtractPedBeforeCommonMode) {
            commonMode[isamp] += pedsubADC[iraw] / double(fN_APV25_CHAN);
          } else {
            commonMode[isamp] += rawADC[iraw] / double(fN_APV25_CHAN);
          }
        } // Pedestal mode

      } // sample loop

      if (fullreadout) {                              // then we need to calculate the common-mode:
        if (fMakeCommonModePlots || !fPedestalMode) { // calculate both ways:

          for (int isamp = 0; isamp < fN_MPD_TIME_SAMP; isamp++) {

            if (fMakeCommonModePlots) {
              double cm_danning        = GetCommonMode(isamp, 1, *it);
              double cm_histo          = GetCommonMode(isamp, 2, *it);
              double cm_sorting        = GetCommonMode(isamp, 0, *it);
              double cm_danning_online = GetCommonMode(isamp, fCommonModeOnlFlag, *it);

              fCM_online[isamp] = cm_danning_online;

              if (!fPedestalMode) {
                switch (fCommonModeFlag) {
                case 2:
                  commonMode[isamp] = cm_histo;
                  break;
                case 1:
                default:
                  commonMode[isamp] = cm_danning;
                  break;
                case 0:
                  commonMode[isamp] = cm_sorting;
                  break;
                }
              }

              double cm_mean;
              UInt_t iAPV = it->pos;

              //if (!CM_OUT_OF_RANGE || fPedestalMode) {
                if (axis == LADGEM::kUaxis) {
                  cm_mean = fCommonModeMeanU[iAPV];
                  /*std::cout << "APV, CMmean, CM_danning, CM_histo, CM_sorting, CM_online: " << iAPV << ", "
                            << cm_mean << ", " << cm_danning << ", " << cm_histo << ", " << cm_sorting << ", "
                            << cm_danning_online << std::endl;
                  */
                  
                  fCommonModeDistU->Fill( iAPV, commonMode[isamp] - cm_mean );
                  fCommonModeDistU_Histo->Fill( iAPV, cm_histo - cm_mean );
                  fCommonModeDistU_Sorting->Fill( iAPV, cm_sorting - cm_mean );
                  fCommonModeDistU_Danning->Fill( iAPV, cm_danning - cm_mean );
                  fCommonModeDiffU->Fill( iAPV, commonMode[isamp] - cm_danning_online );
                  

                } else {
                  cm_mean = fCommonModeMeanV[iAPV];

                  
                  fCommonModeDistV->Fill( iAPV, commonMode[isamp] - cm_mean );
                  fCommonModeDistV_Histo->Fill( iAPV, cm_histo - cm_mean );
                  fCommonModeDistV_Sorting->Fill( iAPV, cm_sorting - cm_mean );
                  fCommonModeDistV_Danning->Fill( iAPV, cm_danning - cm_mean );
                  fCommonModeDiffV->Fill( iAPV, commonMode[isamp] - cm_danning_online );
                  
                }
              //}

            } else if (!fPedestalMode) { // if not doing diagnostic plots, just calculate whichever way the user wanted:
              commonMode[isamp] = GetCommonMode(isamp, fCommonModeFlag, *it);

              if (fCorrectCommonMode) {
                // always calculate online CM if doing corrections:
                fCM_online[isamp] = GetCommonMode(isamp, fCommonModeOnlFlag, *it);
              }
            }

            if (!CM_OUT_OF_RANGE) {
              UpdateRollingAverage(apvcounter, commonMode[isamp], fCommonModeResultContainer_by_APV,
                                   fCommonModeRollingAverage_by_APV, fCommonModeRollingRMS_by_APV,
                                   fNeventsRollingAverage_by_APV);
            }
          } // loop over time samples

          if (fCorrectCommonMode) { // For full readout events we are mainly interested in monitoring the "bias" of the
                                    // ONLINE calculation,
            // NOT correcting the offline calculation
            for (int isamp = 0; isamp < fN_MPD_TIME_SAMP; isamp++) {
              // Are we passing sufficient information for this purpose? Let's see:
              UInt_t ngoodhits  = 0;
              double Correction = GetCommonModeCorrection(isamp, *it, ngoodhits, fN_APV25_CHAN, true);

              double bias = fCM_online[isamp] - Correction - commonMode[isamp];

              if (Correction != 0. && !CM_OUT_OF_RANGE) {
                UpdateRollingAverage(apvcounter, bias, fCMbiasResultContainer_by_APV,
                                     fCommonModeOnlineBiasRollingAverage_by_APV, fCommonModeOnlineBiasRollingRMS_by_APV,
                                     fNeventsOnlineBias_by_APV);
              }

              UInt_t iAPV = it->pos;

              
              // Disable histogram stuff for now
              if( fMakeCommonModePlots ){
                if( Correction != 0. ){

                  double CMbiasDB = ( it->axis == LADGEM::kUaxis ) ? fCMbiasU[iAPV] : fCMbiasV[iAPV];
                  double CMbias = CMbiasDB;

                  if( fNeventsOnlineBias_by_APV[apvcounter] >= std::min( UInt_t(100), std::max(UInt_t(10),
              fN_MPD_TIME_SAMP * fNeventsCommonModeLookBack) ) ){ CMbias =
              fCommonModeOnlineBiasRollingAverage_by_APV[apvcounter];
                  }

                  if( it->axis == LADGEM::kUaxis ){
                    fCommonModeCorrectionU->Fill( iAPV, -Correction );
                    fCommonModeResidualBiasU->Fill( iAPV, fCM_online[isamp] -Correction - commonMode[isamp] );
                    fCommonModeResidualBias_vs_OccupancyU->Fill( double(ngoodhits)/double(fN_APV25_CHAN),
              fCM_online[isamp] - Correction - commonMode[isamp] ); fCommonModeResidualBiasU_corrected->Fill( iAPV,
              fCM_online[isamp] -Correction - commonMode[isamp]
              - 2.0*CMbias*(1.0-double(ngoodhits)/double(fN_APV25_CHAN)) ); } else { fCommonModeCorrectionV->Fill( iAPV,
              -Correction ); fCommonModeResidualBiasV->Fill( iAPV, fCM_online[isamp] - Correction - commonMode[isamp] );
                    fCommonModeResidualBias_vs_OccupancyV->Fill( double(ngoodhits)/double(fN_APV25_CHAN),
              fCM_online[isamp] - Correction - commonMode[isamp] ); fCommonModeResidualBiasV_corrected->Fill( iAPV,
              fCM_online[isamp] -Correction - commonMode[isamp]
              - 2.0*CMbias*(1.0-double(ngoodhits)/double(fN_APV25_CHAN)) );
                  }
                } else {
                  if( it->axis == LADGEM::kUaxis ){
                    fCommonModeDiffU_Uncorrected->Fill( iAPV, commonMode[isamp] - fCM_online[isamp] );
                  } else {
                    fCommonModeDiffV_Uncorrected->Fill( iAPV, commonMode[isamp] - fCM_online[isamp] );
                  }
                }
              }// fMakeCommonModePlots
              
            }
          } // fCorrectCommonMode

        } // check if conditions are satisfied to require offline common-mode calculation
      } // End check !CM_ENABLED && BUILD_ALL_SAMPLES

      if (CM_ENABLED && fCorrectCommonMode) {
        // Under certain conditions we want to attempt to correct the ADC values for all strips on an APV card using
        // either the rolling average over a certain number of previous events, or the CM mean from the database. There
        // are two conditions that must be satisfied to attempt correcting the ADC values for an event: 1) The online
        // calculated CM must be more than some number of std. deviations below the "expected" CM according to the
        // rolling average 2) The number of strips with ADC values within some number of std. deviations of the
        // "expected" CM must exceed some threshold to allow us to
        //    to obtain a new estimate of the "true" common-mode for that event
        // If both 1) and 2) are satisfied, then we will attempt a new common-mode calculation using the strips that
        // passed zero suppression using the online common-mode calculation.

        for (int isamp = 0; isamp < fN_MPD_TIME_SAMP; isamp++) {

          UInt_t ngood     = 0;
          UInt_t nhitstemp = UInt_t(nstrips);

          // std::cout << "Attempting common-mode correction for online zero-suppressed event sample " << isamp <<
          // "...";

          CommonModeCorrection[isamp] = GetCommonModeCorrection(isamp, *it, ngood, nhitstemp);

          if (CommonModeCorrection[isamp] != 0.0) { // if we are applying a correction, correct it for bias:
            UInt_t iAPV = it->pos;

            double CMbiasDB = (it->axis == LADGEM::kUaxis) ? fCMbiasU[iAPV] : fCMbiasV[iAPV];

            double CMbias = CMbiasDB;

            if (fNeventsOnlineBias_by_APV[apvcounter] >=
                std::min(UInt_t(100), std::max(UInt_t(10), fN_MPD_TIME_SAMP * fNeventsCommonModeLookBack))) {
              CMbias = fCommonModeOnlineBiasRollingAverage_by_APV[apvcounter];
            }

            // bias is DEFINED as Online common-mode MINUS correction MINUS "true" common-mode:
            //"correction" is DEFINED as Online common-mode MINUS "corrected common-mode" and is to be ADDED to the ADC
            // values:

            // bias = online CM - (online CM - corrected CM) - true CM = corrected CM - true CM
            // --> true CM = corrected CM - bias
            // corrected ADC = ADC + online CM - true CM = uncorrected ADC + [online CM - (corrected CM - bias)]
            // = uncorrected ADC + [correction + bias]
            // [...] = correction to be ADDED to ADC
            // --> corrected correction = correction + bias
            CommonModeCorrection[isamp] += 2.0 * CMbias * (1.0 - double(ngood) / double(fN_APV25_CHAN));
          }
        }
      } // Here we are done with common mode stuff

      // Last loop over all the strips and samples in the data and populate/calculate global variables that are passed
      // to track-finding:
      for (Int_t istrip = 0; istrip < nstrips; ++istrip) {

        // Temporary vector to hold ped-subtracted ADC samples for this strip:
        std::vector<double> ADCtemp(fN_MPD_TIME_SAMP);
        std::vector<int> rawADCtemp(fN_MPD_TIME_SAMP);
        std::vector<Double_t> DeconvADCtemp(fN_MPD_TIME_SAMP, 0.0);

        // sums over time samples
        double ADCsum_temp = 0.0;
        double maxADC      = 0.0;
        double minADC      = 10000.0; // Negative pulse
        Int_t iSampMax     = -1;
        Int_t iSampMin     = -1; // Negative pulse

        // crude timing calculations:
        double Tsum  = 0.0;
        double T2sum = 0.0;

        // grab decoded strip number directly:
        int strip = Strip[fN_MPD_TIME_SAMP * istrip];

        //"pedtemp" is only used to fill pedestal histograms as of now:
        // double pedtemp = ( axis == LADGEM::kUaxis ) ? fPedestalU[strip] : fPedestalV[strip];
        // if( fPedSubFlag != 0 && !fIsMC && !fPedestalMode ) pedtemp = 0.0;

        double rmstemp  = (axis == LADGEM::kUaxis) ? fPedRMSU[strip] : fPedRMSV[strip];
        double gaintemp = (axis == LADGEM::kUaxis) ? fUgain[strip / fN_APV25_CHAN] : fVgain[strip / fN_APV25_CHAN];

        // Now loop over the time samples:
        for (Int_t adc_samp = 0; adc_samp < fN_MPD_TIME_SAMP; adc_samp++) {

          int iraw = adc_samp + fN_MPD_TIME_SAMP * istrip;

          // We need to subtract the common-mode if it was calculated offline:
          if (!CM_ENABLED && BUILD_ALL_SAMPLES && nstrips == fN_APV25_CHAN) {
            commonModeSubtractedADC[iraw] = pedsubADC[iraw] - commonMode[adc_samp];
          }

          if (fCorrectCommonMode)
            commonModeSubtractedADC[iraw] += CommonModeCorrection[adc_samp];

          Int_t RawADC         = rawADC[iraw]; // this value has no corrections applied:
          rawADCtemp[adc_samp] = RawADC;       // raw only

          // The following value already has pedestal and common-mode subtracted (if applicable):
          double ADCvalue = commonModeSubtractedADC[iraw]; // zero-suppress BEFORE we apply gain correction

          ADCtemp[adc_samp] = ADCvalue; // common-mode AND pedestal subtracted

          ADCsum_temp += ADCvalue;

          if (iSampMax < 0 || ADCvalue > maxADC) {
            maxADC   = ADCvalue;
            iSampMax = adc_samp;
          }

          ///// Used for negative pulse study
          if (iSampMin < 0 || ADCvalue < minADC) {
            minADC   = ADCvalue;
            iSampMin = adc_samp;
          }

          //	  cout << "nstrip,correctCM: " << nstrips << " " << fCorrectCommonMode << " " << pedsubADC[iraw] << " "
          //<< commonMode[adc_samp] << endl; 	  cout << "istrip, iraw, samp, rawADC, ADC: " << istrip << " " << iraw
          //<< " "
          //<< adc_samp << " " << RawADC << " " << ADCvalue << endl;

          // for crude strip timing, just take simple time bins at the center of each sample (we'll worry about trigger
          // time words later):
          double Tsamp = fSamplePeriod * (adc_samp + 0.5);

          Tsum += Tsamp * ADCvalue;
          T2sum += pow(Tsamp, 2) * ADCvalue;

        } // loop over time sample

        CalcDeconvolutedSamples(ADCtemp, DeconvADCtemp);

        // Disable histogram stuff for now
        
        if( (fPedestalMode || fMakeCommonModePlots) && !CM_ENABLED ){
          int iAPV = strip/fN_APV25_CHAN;
          double pedtemp = ( axis == LADGEM::kUaxis ) ? fPedestalU[strip] : fPedestalV[strip];
          if( axis == LADGEM::kUaxis ){
            for( int isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){

              hrawADCs_by_stripU->Fill( strip, rawADCtemp[isamp] );
              hpedestal_subtracted_ADCs_by_stripU->Fill( strip, ADCtemp[isamp] ); //common-mode AND ped-subtracted
              hcommonmode_subtracted_ADCs_by_stripU->Fill( strip, ADCtemp[isamp] + pedtemp ); //common-mode subtraction only, no ped: 
              hpedestal_subtracted_rawADCs_by_stripU->Fill( strip, ADCtemp[isamp] + commonMode[isamp] );
        //pedestal subtraction only, no common-mode 
              hpedestal_subtracted_rawADCsU->Fill( ADCtemp[isamp] +commonMode[isamp] ); //1D distribution of ped-subtracted ADCs w/o common-mode subtraction
              hpedestal_subtracted_ADCsU->Fill( ADCtemp[isamp] ); //1D distribution of ped-and-common-mode subtracted ADCs  
              hcommonmode_mean_by_APV_U->Fill( iAPV, commonMode[isamp] );

              //if( iSampMax != 0 ){
              //  hdeconv_ADCsU->Fill( DeconvADCtemp[isamp] );
             // }
            }
          } else {
            for( int isamp=0; isamp<fN_MPD_TIME_SAMP; isamp++ ){

              hrawADCs_by_stripV->Fill( strip, rawADCtemp[isamp] );
              hpedestal_subtracted_ADCs_by_stripV->Fill( strip, ADCtemp[isamp] );
              hcommonmode_subtracted_ADCs_by_stripV->Fill( strip, ADCtemp[isamp] + pedtemp );
              hpedestal_subtracted_rawADCs_by_stripV->Fill( strip, ADCtemp[isamp] + commonMode[isamp] );
              hpedestal_subtracted_rawADCsV->Fill( ADCtemp[isamp] + commonMode[isamp] );
              hpedestal_subtracted_ADCsV->Fill( ADCtemp[isamp] );
              hcommonmode_mean_by_APV_V->Fill( iAPV, commonMode[isamp] );

              //if( iSampMax != 0 ){
              //  hdeconv_ADCsV->Fill( DeconvADCtemp[isamp] );
              //}

            }
          }
        }
        

        // the ROOTgui multicrate uses a threshold on the AVERAGE ADC sample (not the MAX). To be consistent
        // with how the "Hit" root files are produced, let's use the same threshold;
        // this amounts to using a higher effective threshold than cutting on the max ADC sample would have been:
        // fThresholdStripSum is in many respects redundant with fZeroSuppressRMS
        if (!fZeroSuppress ||
            (ADCsum_temp / double(fN_MPD_TIME_SAMP) > fZeroSuppressRMS * rmstemp)) { // Default threshold is 5-sigma!
          // Increment hit count and populate decoded data structures:
          // threshold on the average ADC

          // Slight reorganization: compute Tmean and Tsigma before applying gain correction:
          //(since these sums were computed using the uncorrected ADC samples)
          double Tmean_temp  = Tsum / ADCsum_temp;
          double Tsigma_temp = sqrt(T2sum / ADCsum_temp - pow(Tmean_temp, 2));

          // NOW apply gain correction:
          // Don't apply any gain correction if we are doing pedestal mode analysis:
          if (!fPedestalMode) { // only apply gain correction if we aren't in pedestal-mode:
            for (Int_t isamp = 0; isamp < fN_MPD_TIME_SAMP; isamp++) {
              ADCtemp[isamp] *= gaintemp;
              DeconvADCtemp[isamp] *= gaintemp;
            }
            maxADC *= gaintemp;
            ADCsum_temp *= gaintemp;
          }

          fStrip[fNstrips_hit]    = strip;
          fAxis[fNstrips_hit]     = axis;
          fStripRaw[fNstrips_hit] = rawStrip[fN_MPD_TIME_SAMP * istrip];

          fKeepStrip[fNstrips_hit] = true;

          fMaxSamp[fNstrips_hit] = iSampMax;

          // fSuppressFirstLast:
          // 0 = allow peaking in first or last sample
          // 1 = suppress peaking in first and last sample
          // -1 = suppress peaking in first sample only (or other negative number)
          // -2 = suppress peaking in last sample only:
          if (fDeconvolutionFlag ==
              0) { // if "deconvolution flag" is non-zero, then set "keep strip" based on deconvoluted variables
            if (fSuppressFirstLast != 0) {
              bool peakfirst = iSampMax == 0;
              bool peaklast  = iSampMax + 1 == fN_MPD_TIME_SAMP;

              if (peakfirst) {
                if (fSuppressFirstLast > 0 || fSuppressFirstLast != -2) {
                  fKeepStrip[fNstrips_hit] = false;
                }
              } else if (peaklast) {
                if (fSuppressFirstLast > 0 || fSuppressFirstLast == -2) {
                  fKeepStrip[fNstrips_hit] = false;
                }
              }
            }
          }

          double ADCsum_deconv = 0.0;
          double maxdeconv     = 0.0;
          int imaxdeconv       = 0;

          double Tsum_deconv = 0.0;

          double maxcombo = 0.0;

          double combotemp = 0.0;
          int imaxcombo    = 0;

          // we don't need to repeat the calculation of deconvoluted samples. That was moved to its own helper method:

          for (Int_t isamp = 0; isamp < fN_MPD_TIME_SAMP; isamp++) {
            fADCsamples[fNstrips_hit][isamp]    = ADCtemp[isamp];
            fRawADCsamples[fNstrips_hit][isamp] = rawADCtemp[isamp];

            fADCsamples_deconv[fNstrips_hit][isamp] = DeconvADCtemp[isamp];

            ADCsum_deconv += DeconvADCtemp[isamp];

            if (isamp == 0) {
              combotemp = DeconvADCtemp[isamp];
              maxcombo  = combotemp;
              imaxcombo = isamp;
            } else {
              combotemp = DeconvADCtemp[isamp] + DeconvADCtemp[isamp - 1];
              if (combotemp > maxcombo) {
                maxcombo  = combotemp;
                imaxcombo = isamp;
              }
            }
            if (isamp == 5) {
              combotemp = DeconvADCtemp[isamp];
              if (combotemp > maxcombo) {
                maxcombo  = combotemp;
                imaxcombo = fN_MPD_TIME_SAMP;
              }
            }

            if (isamp == 0 || DeconvADCtemp[isamp] > maxdeconv) {
              imaxdeconv = isamp;
              maxdeconv  = DeconvADCtemp[isamp];
            }

            Tsum_deconv += (fSamplePeriod * (isamp + 0.5)) * DeconvADCtemp[isamp];

            fADCsamples1D[isamp + fN_MPD_TIME_SAMP * fNstrips_hit]       = ADCtemp[isamp];
            fRawADCsamples1D[isamp + fN_MPD_TIME_SAMP * fNstrips_hit]    = rawADCtemp[isamp];
            fADCsamplesDeconv1D[isamp + fN_MPD_TIME_SAMP * fNstrips_hit] = fADCsamples_deconv[fNstrips_hit][isamp];

            // Disable histogram stuff for now
            /*
            if( fKeepStrip[fNstrips_hit] && hADCfrac_vs_timesample_allstrips != NULL ){
              hADCfrac_vs_timesample_allstrips->Fill( isamp, ADCtemp[isamp]/ADCsum_temp );
            }
            */
          }

          fADCsumsDeconv[fNstrips_hit] = ADCsum_deconv;

          fStripTrackIndex[fNstrips_hit] = -1;
          fStripOnTrack[fNstrips_hit]    = 0;

          // This block is only used for the negative signal studies
          fStripIsNeg[fNstrips_hit]         = 0;
          fStripIsNegU[fNstrips_hit]        = 0;
          fStripIsNegV[fNstrips_hit]        = 0;
          fStripIsNegOnTrack[fNstrips_hit]  = 0;
          fStripIsNegOnTrackU[fNstrips_hit] = 0;
          fStripIsNegOnTrackV[fNstrips_hit] = 0;

          // These are used for saving numbers to a text file for event displays
          fStripEvent[fNstrips_hit]  = evdata.GetEvNum();
          fStripCrate[fNstrips_hit]  = it->crate;
          fStripMPD[fNstrips_hit]    = it->mpd_id;
          fStripADC_ID[fNstrips_hit] = it->adc_id;

          fADCmax[fNstrips_hit]            = maxADC;
          fADCmaxDeconv[fNstrips_hit]      = maxdeconv;
          fADCmaxDeconvCombo[fNstrips_hit] = maxcombo;

          fMaxSampDeconv[fNstrips_hit]      = imaxdeconv;
          fMaxSampDeconvCombo[fNstrips_hit] = imaxcombo;

          if (fDeconvolutionFlag != 0) {
            // rmstemp is the rms of the average of six time samples. To get individual sample noise, we take:
            double sigma_1sample = rmstemp * fRMS_ConversionFactor;

            // 5 * 8 * sqrt(6) ~= 100

            if (maxcombo <= fZeroSuppressRMS * sigma_1sample) {
              fKeepStrip[fNstrips_hit] = false;
            }

            if (fSuppressFirstLast != 0 && imaxcombo == 0) {
              fKeepStrip[fNstrips_hit] = false;
            }
          }

          fTmeanDeconv[fNstrips_hit] = Tsum_deconv / ADCsum_deconv - fTrigTimeSlope * fTrigTime;

          fTmean[fNstrips_hit] = Tmean_temp - fTrigTimeSlope * fTrigTime;

          fTsigma[fNstrips_hit] = Tsigma_temp;

          fStripTSchi2[fNstrips_hit] = StripTSchi2(fNstrips_hit);
          fStripTSprob[fNstrips_hit] = TMath::Prob(fStripTSchi2[fNstrips_hit], 6);

          if (fUseTSchi2cut && fStripTSchi2[fNstrips_hit] > fStripTSchi2Cut) {
            fKeepStrip[fNstrips_hit] = false;
          }

          fStripTdiff[fNstrips_hit]     = -1000.; // This will become meaningful only at the clustering stage
          fStripCorrCoeff[fNstrips_hit] = -1000.; // This will become meaningful only at the clustering stage
          fTcorr[fNstrips_hit]          = fTmean[fNstrips_hit];

          fStripTfit[fNstrips_hit] = FitStripTime(fNstrips_hit, rmstemp * 2.45);

          fStrip_ENABLE_CM[fNstrips_hit]         = CM_ENABLED;
          fStrip_CM_GOOD[fNstrips_hit]           = !CM_OUT_OF_RANGE;
          fStrip_BUILD_ALL_SAMPLES[fNstrips_hit] = BUILD_ALL_SAMPLES;

          fADCsums[fNstrips_hit] = ADCsum_temp;

          fStripADCavg[fNstrips_hit] = ADCsum_temp / double(fN_MPD_TIME_SAMP);

          UInt_t isU = (axis == LADGEM::kUaxis) ? 1 : 0;
          UInt_t isV = (axis == LADGEM::kVaxis) ? 1 : 0;

          fStripIsU[fNstrips_hit] = isU;
          fStripIsV[fNstrips_hit] = isV;

          fNstrips_hitU += isU;
          fNstrips_hitV += isV;

          if (fKeepStrip[fNstrips_hit]) {
            fNstrips_keep++;
            fNstrips_keepU += isU;
            fNstrips_keepV += isV;
            if (fADCmax[fNstrips_hit] >= fThresholdSample && fADCsums[fNstrips_hit] >= fThresholdStripSum) {
              fNstrips_keep_lmax++;
              fNstrips_keep_lmaxU += isU;
              fNstrips_keep_lmaxV += isV;
            }
          }

          fNstrips_hit++;
          fNstrips_hit_pos++;

        } // check if passed zero suppression cuts

        // FIXME: Decide we need to keep it. If so, try to do it in one loop instead of repeating it -- SP
        /////// Negative pulse study, This is an exact copy of the loop above but instead stores the negative ADC info.
        ///"Keep" is set
        /////// to false regardless so these strips will not be used for any of the normal clustering and tracking
        /// algorithms. They
        /////// are differentiated from positive strips by fStripIsNeg.
        if (ADCsum_temp / double(fN_MPD_TIME_SAMP) < -1.0 * fZeroSuppressRMS * rmstemp && BUILD_ALL_SAMPLES &&
            !CM_ENABLED && fNegSignalStudy) {

          // Increment hit count and populate decoded data structures:

          // threshold on the average ADC

          // Slight reorganization: compute Tmean and Tsigma before applying gain correction:
          //(since these sums were computed using the uncorrected ADC samples)
          double Tmean_temp  = Tsum / ADCsum_temp;
          double Tsigma_temp = sqrt(-1.0 * T2sum / ADCsum_temp - pow(Tmean_temp, 2));

          // NOW apply gain correction:
          // Don't apply any gain correction if we are doing pedestal mode analysis:
          if (!fPedestalMode) { // only apply gain correction if we aren't in pedestal-mode:
            for (Int_t isamp = 0; isamp < fN_MPD_TIME_SAMP; isamp++) {
              ADCtemp[isamp] *= gaintemp;
            }
            minADC *= gaintemp; // Use min ADC instead of max for negative signals
            ADCsum_temp *= gaintemp;
          }

          fStrip[fNstrips_hit]    = strip;
          fAxis[fNstrips_hit]     = axis;
          fStripRaw[fNstrips_hit] = rawStrip[fN_MPD_TIME_SAMP * istrip];

          fMaxSamp[fNstrips_hit] = iSampMin;

          // fSuppressFirstLast:
          // 0 = allow peaking in first or last sample
          // 1 = suppress peaking in first and last sample
          // -1 = suppress peaking in first sample only (or other negative number)
          // -2 = suppress peaking in last sample only:
          if (fSuppressFirstLast != 0) {
            bool peakfirst = iSampMin == 0;
            bool peaklast  = iSampMin + 1 == fN_MPD_TIME_SAMP;

            if (peakfirst) {
              if (fSuppressFirstLast > 0 || fSuppressFirstLast != -2) {
                fKeepStrip[fNstrips_hit] = false;
              }
            } else if (peaklast) {
              if (fSuppressFirstLast > 0 || fSuppressFirstLast == -2) {
                fKeepStrip[fNstrips_hit] = false;
              }
            }
          }

          for (Int_t isamp = 0; isamp < fN_MPD_TIME_SAMP; isamp++) {
            fADCsamples[fNstrips_hit][isamp]    = ADCtemp[isamp];
            fRawADCsamples[fNstrips_hit][isamp] = rawADCtemp[isamp];

            fADCsamples1D[isamp + fN_MPD_TIME_SAMP * fNstrips_hit]    = ADCtemp[isamp];
            fRawADCsamples1D[isamp + fN_MPD_TIME_SAMP * fNstrips_hit] = rawADCtemp[isamp];

            // Disable histogram stuff for now
            /*
            if( fKeepStrip[fNstrips_hit] && hADCfrac_vs_timesample_allstrips != NULL ){
              hADCfrac_vs_timesample_allstrips->Fill( isamp, ADCtemp[isamp]/ADCsum_temp );
            }
            */
          }

          fStripTrackIndex[fNstrips_hit] = -1;
          fStripOnTrack[fNstrips_hit]    = 0;

          // Variables used for negative signal studies
          fStripIsNeg[fNstrips_hit]         = 1;
          fStripIsNegU[fNstrips_hit]        = (axis == LADGEM::kUaxis) ? 1 : 0;
          fStripIsNegV[fNstrips_hit]        = (axis == LADGEM::kVaxis) ? 1 : 0;
          fStripIsNegOnTrack[fNstrips_hit]  = 0;
          fStripIsNegOnTrackU[fNstrips_hit] = 0;
          fStripIsNegOnTrackV[fNstrips_hit] = 0;

          // These are used for saving numbers to a text file for event displays
          fStripEvent[fNstrips_hit]  = evdata.GetEvNum();
          fStripCrate[fNstrips_hit]  = it->crate;
          fStripMPD[fNstrips_hit]    = it->mpd_id;
          fStripADC_ID[fNstrips_hit] = it->adc_id;

          // Variables used for negative signal studies

          fADCmax[fNstrips_hit] = minADC; // Use minADC instead for negative strips

          fTmean[fNstrips_hit] = Tmean_temp;

          fTsigma[fNstrips_hit] = Tsigma_temp;

          fStripTSchi2[fNstrips_hit] = StripTSchi2(fNstrips_hit);

          if (fUseTSchi2cut && fStripTSchi2[fNstrips_hit] > fStripTSchi2Cut) {
            fKeepStrip[fNstrips_hit] = false;
          }

          fStripTdiff[fNstrips_hit]     = -1000.; // This will become meaningful only at the clustering stage
          fStripCorrCoeff[fNstrips_hit] = -1000.; // This will become meaningful only at the clustering stage
          fTcorr[fNstrips_hit]          = fTmean[fNstrips_hit];

          fStripTfit[fNstrips_hit] = FitStripTime(fNstrips_hit, rmstemp * 2.45);

          fStrip_ENABLE_CM[fNstrips_hit]         = CM_ENABLED;
          fStrip_CM_GOOD[fNstrips_hit]           = !CM_OUT_OF_RANGE;
          fStrip_BUILD_ALL_SAMPLES[fNstrips_hit] = BUILD_ALL_SAMPLES;

          fADCsums[fNstrips_hit] = ADCsum_temp;

          fStripADCavg[fNstrips_hit] = ADCsum_temp / double(fN_MPD_TIME_SAMP);

          UInt_t isU = (axis == LADGEM::kUaxis) ? 1 : 0;
          UInt_t isV = (axis == LADGEM::kVaxis) ? 1 : 0;

          fStripIsU[fNstrips_hit] = isU;
          fStripIsV[fNstrips_hit] = isV;

          fNstrips_hitU_neg += isU;
          fNstrips_hitV_neg += isV;

          if (fKeepStrip[fNstrips_hit]) {
            fNstrips_keep++;
            fNstrips_keepU += isU;
            fNstrips_keepV += isV;
            if (fADCmax[fNstrips_hit] >= fThresholdSample && fADCsums[fNstrips_hit] >= fThresholdStripSum) {
              fNstrips_keep_lmax++;
              fNstrips_keep_lmaxU += isU;
              fNstrips_keep_lmaxV += isV;
            }
          }

          fNstrips_hit++;
          fNstrips_hit_neg++;

        } // end loop over negative zero suppression for full readout events

      } // end loop over strips on this APV card
    } // end if( nsamp > 0 )

    apvcounter++;
  } // end loop on decode map entries for this module

  fNdecoded_ADCsamples = fNstrips_hit * fN_MPD_TIME_SAMP;

  fIsDecoded = true;

  return 0;
}

//____________________________________________________________________________________
Int_t THcLADGEMModule::CoarseProcess(TClonesArray &tracks) {
  //  cout << "THcLADGEMModule::CoarseProcess" << endl;
  // Find 1D clusters for each axis
  FindClusters1D(LADGEM::kUaxis); // +input ucenter, 0.5*(umax-umin) for u strips
  FindClusters1D(LADGEM::kVaxis); // +input ucenter, 0.5*(umax-umin) for v strips

  // Find 2D hits
  Find2DHits();

  return 0;
}

//____________________________________________________________________________________
Int_t THcLADGEMModule::FineProcess(TClonesArray &tracks) {
  // cout << "THcLADGEMModule::FineProcess" << endl;

  return 0;
}

//____________________________________________________________________________________
void THcLADGEMModule::FindClusters1D(LADGEM::GEMaxis_t axis) {

  UShort_t maxsep      = (axis == LADGEM::kUaxis) ? fMaxNeighborsU_totalcharge : fMaxNeighborsV_totalcharge;
  UShort_t maxsepcoord = (axis == LADGEM::kUaxis) ? fMaxNeighborsU_hitpos : fMaxNeighborsV_hitpos;
  UInt_t Nstrips       = (axis == LADGEM::kUaxis) ? fNstripsU : fNstripsV;
  Double_t pitch       = (axis == LADGEM::kUaxis) ? fUStripPitch : fVStripPitch;
  Double_t offset      = (axis == LADGEM::kUaxis) ? fUStripOffset : fVStripOffset;

  // Temporary containers to store strip hit information to iterate
  std::set<UShort_t> striplist;        // sorted list of strips for 1D clustering
  std::map<UShort_t, UInt_t> hitindex; // key = strip ID, mapped value = index in decoded hit array, needed to access
                                       // the other information efficiently:
  std::map<UShort_t, Double_t> pedrms_strip;
  std::map<UShort_t, Double_t>
      ADC_strip; // These are the (configuration-dependent) quantities we use for clustering. They depend on the values
                 // of fClusteringFlag and fSuppressFirstLast and fDeconvolution_flag
  std::map<UShort_t, Double_t> ADC_maxsamp;  //
  std::map<UShort_t, Double_t> Tmean_strip;  // strip mean time with first and/or last samples removed (if applicable)
  std::map<UShort_t, Double_t> Tfit_strip;   // strip "fit" time
  std::map<UShort_t, Double_t> Tsigma_strip; // strip rms time with first and/or last samples removed (if applicable)

  std::set<UShort_t> striplist_neg; // same as above but for negative strips
  std::map<UShort_t, UInt_t> hitindex_neg;
  std::map<UShort_t, Double_t> pedrms_strip_neg;

  for (int ihit = 0; ihit < fNstrips_hit; ihit++) {
    if (fAxis[ihit] == axis) {

      bool newstrip = (striplist.insert(fStrip[ihit])).second;

      if (newstrip) { // should always be true:
        hitindex[fStrip[ihit]] = ihit;
        if (axis == LADGEM::kUaxis) {
          pedrms_strip[fStrip[ihit]] = fPedRMSU[fStrip[ihit]];
        } else {
          pedrms_strip[fStrip[ihit]] = fPedRMSV[fStrip[ihit]];
        }

        // Default flag == 0
        // using sums of ADC values over all time
        // samples on a strip
        ADC_strip[fStrip[ihit]]    = fADCsums[ihit];
        ADC_maxsamp[fStrip[ihit]]  = fADCmax[ihit];
        Tmean_strip[fStrip[ihit]]  = fTmean[ihit];
        Tfit_strip[fStrip[ihit]]   = fStripTfit[ihit];
        Tsigma_strip[fStrip[ihit]] = fTsigma[ihit];

        if (fClusteringFlag == 1) {
          ADC_strip[fStrip[ihit]]   = fADCmaxDeconvCombo[ihit];
          ADC_maxsamp[fStrip[ihit]] = fADCmaxDeconv[ihit];
          Tmean_strip[fStrip[ihit]] = fTmeanDeconv[ihit];
        }
      }
    }
    // FIXME: Do we need to add option for using negative signal?
  }

  //*************************
  // finding local maxima
  //*************************

  std::set<UShort_t> localmaxima;
  std::map<UShort_t, bool> islocalmax;

  for (std::set<UShort_t>::iterator i = striplist.begin(); i != striplist.end(); ++i) {
    int strip = *i;

    islocalmax[strip] = false;

    double sumstrip = ADC_strip[strip];
    double sumleft  = 0.0;
    double sumright = 0.0;

    // check neighbors
    if (striplist.find(strip - 1) != striplist.end())
      sumleft = ADC_strip[strip - 1];

    if (striplist.find(strip + 1) != striplist.end())
      sumright = ADC_strip[strip + 1];

    double thresh_samp  = fThresholdSample;
    double thresh_strip = fThresholdStripSum;
    if (fClusteringFlag == 1) {
      thresh_strip = fThresholdDeconvADCMaxCombo;
      thresh_samp  = fThresholdSampleDeconv;
    }

    // Decide if good ADC/time or not
    bool goodADC = sumstrip >= thresh_strip && ADC_maxsamp[strip] >= thresh_samp;
    if (goodADC && sumstrip >= sumleft && sumstrip >= sumright) {

      bool goodtime = true;

      double tstrip = Tmean_strip[strip];
      double t0     = fStripMaxTcut_central[axis];
      double tcut   = fStripMaxTcut_width[axis];
      double tsigma = fStripMaxTcut_sigma[axis];

      if (fUseStripTimingCuts == 2 && fClusteringFlag != 1) { // alternate timing cut based on strip "fitted" time
        tstrip = Tfit_strip[strip];
        t0     = fStripMaxTcut_central_fit[axis];
        tcut   = fStripMaxTcut_width_fit[axis];
        tsigma = fStripMaxTcut_sigma_fit[axis];
      }
      if (fClusteringFlag == 1) {
        t0     = fStripMaxTcut_central_deconv[axis];
        tcut   = fStripMaxTcut_width_deconv[axis];
        tsigma = fStripMaxTcut_sigma_deconv[axis];
      }

      if (fUseStripTimingCuts != 0 && fabs(tstrip - t0) > tcut * tsigma)
        goodtime = false;

      if (goodtime && fKeepStrip[hitindex[strip]]) {
        islocalmax[strip] = true;
        localmaxima.insert(strip);
      }
    }
  }

  //************************
  // Compare peaks,
  // and remove insignificant peaks
  //************************

  vector<int> peakstoerase;

  // Loop over local maxima
  for (std::set<UShort_t>::iterator i = localmaxima.begin(); i != localmaxima.end(); ++i) {
    int stripmax = *i;

    double ADCmax     = ADC_strip[stripmax];
    double prominence = ADCmax;

    int striplo = stripmax, striphi = stripmax;
    double ADCminright = ADCmax, ADCminleft = ADCmax;

    bool higherpeakright = false, higherpeakleft = false;
    int peakright = -1, peakleft = -1;

    // scan to higher# strips
    while (striplist.find(striphi + 1) != striplist.end()) {
      striphi++;

      Double_t ADCtest = ADC_strip[striphi];

      if (ADCtest < ADCminright && !higherpeakright) {
        ADCminright = ADCtest;
      }

      if (islocalmax[striphi] && ADCtest > ADCmax) {
        // then this peak is in a contiguous group with another higher peak to the right:
        higherpeakright = true;
        peakright       = striphi;
      }
    }

    // scan to lower# strips
    while (striplist.find(striplo - 1) != striplist.end()) {
      striplo--;

      Double_t ADCtest = ADC_strip[striplo];

      if (ADCtest < ADCminleft && !higherpeakleft) {
        ADCminleft = ADCtest;
      }

      if (islocalmax[striplo] && ADCtest > ADCmax) {
        // then this peak is in a contiguous group with another higher peak to the left:
        higherpeakleft = true;
        peakleft       = striplo;
      }
    }

    // The strip RMS values represent the RMS of the AVERAGE of the samples.
    // So the prominence threshold should be expressed in terms of the same thing to be consistent:
    // RMS of the sum = (rms avg) * (number of samples)
    double sigma_sum = double(fN_MPD_TIME_SAMP) * pedrms_strip[stripmax];

    // What is the effect of deconvolution on noise?
    // It turns out from looking at the width of the pedestal peak
    // in the deconvoluted ADCs from
    // full readout events that the noise in the deconvoluted samples
    // is about the same as the noise in the regular samples.
    // HOWEVER: in this case we aren't working with the sum of six samples. So we need to
    // modify "sigma" accordingly. We are generally working with the sum of two deconvoluted
    // samples. ASSUMING the deconvoluted sample width is the same as the individual
    // sample width, we have:
    if (fClusteringFlag == 1) {
      sigma_sum = pedrms_strip[stripmax] * fRMS_ConversionFactor;
    }

    // Mark which peaks to remove
    bool peak_close = false;
    if (!higherpeakleft)
      ADCminleft = 0.0;
    if (!higherpeakright)
      ADCminright = 0.0;

    if (higherpeakright || higherpeakleft) {
      // this peak is contiguous with higher peaks on either the left or right or both:

      // subtract the higher of the two valleys to get the prominence
      prominence = ADCmax - std::max(ADCminleft, ADCminright);

      if (higherpeakleft && std::abs(peakleft - stripmax) <= 2 * maxsep)
        peak_close = true;
      if (higherpeakright && std::abs(peakright - stripmax) <= 2 * maxsep)
        peak_close = true;

      if (peak_close &&
          (prominence < fThresh_2ndMax_nsigma * sigma_sum || prominence / ADCmax < fThresh_2ndMax_fraction)) {

        peakstoerase.push_back(stripmax);
      }
    }
  } // loop over local maxima

  // Remove insignificant peaks marked from above
  for (int ipeak : peakstoerase) {
    localmaxima.erase(ipeak);
    islocalmax[ipeak] = false;
  }

  //***********************
  // Form/split clusters
  //***********************

  // Loop over remaining local maxima
  for (auto i = localmaxima.begin(); i != localmaxima.end(); ++i) {
    int stripmax = *i;
    int striplo  = stripmax;
    int striphi  = stripmax;

    double ADCmax = ADC_strip[stripmax];

    // Scan to the left
    bool found_neighbor_low = true;
    while (found_neighbor_low) {

      found_neighbor_low = striplist.find(striplo - 1) != striplist.end() && stripmax - striplo < maxsep;

      // Time difference
      double Tdiff = fTmean[hitindex[striplo - 1]] - fTmean[hitindex[stripmax]];

      if (fUseStripTimingCuts == 2 && fClusteringFlag == 0) { // Use "fitted" strip time instead:
        Tdiff = fStripTfit[hitindex[striplo - 1]] - fStripTfit[hitindex[stripmax]];
      }
      if (fClusteringFlag == 1) {
        Tdiff = fTmeanDeconv[hitindex[striplo - 1]] - fTmeanDeconv[hitindex[stripmax]];
      }

      // FIXME: check if we need this cut
      double Ccoeff =
          CorrCoeff(fN_MPD_TIME_SAMP, fADCsamples[hitindex[striplo - 1]], fADCsamples[hitindex[stripmax]], 0);

      // correlation coefficient of the deconvoluted samples:
      double Ccoeff_deconv = CorrCoeff(fN_MPD_TIME_SAMP, fADCsamples_deconv[hitindex[striplo - 1]],
                                       fADCsamples_deconv[hitindex[stripmax]]);

      if (fUseStripTimingCuts != 0 && fabs(Tdiff) > fStripAddTcut_width)
        found_neighbor_low = false;

      double Ccoeff_test = fClusteringFlag == 1 ? Ccoeff_deconv : Ccoeff;

      if (Ccoeff_test < fStripAddCorrCoeffCut)
        found_neighbor_low = false;

      if (found_neighbor_low)
        striplo--;
    }

    // Repeat, Scan to the right
    bool found_neighbor_high = true;
    while (found_neighbor_high) {

      found_neighbor_high = striplist.find(striphi + 1) != striplist.end() && striphi - stripmax < maxsep;

      // Time difference
      double Tdiff = fTmean[hitindex[striphi + 1]] - fTmean[hitindex[stripmax]];

      if (fUseStripTimingCuts == 2 && fClusteringFlag == 0) { // Use "fitted" strip time instead:
        Tdiff = fStripTfit[hitindex[striphi + 1]] - fStripTfit[hitindex[stripmax]];
      }
      if (fClusteringFlag == 1) {
        Tdiff = fTmeanDeconv[hitindex[striphi + 1]] - fTmeanDeconv[hitindex[stripmax]];
      }

      double Ccoeff =
          CorrCoeff(fN_MPD_TIME_SAMP, fADCsamples[hitindex[striphi + 1]], fADCsamples[hitindex[stripmax]], 0);

      double Ccoeff_deconv = CorrCoeff(fN_MPD_TIME_SAMP, fADCsamples_deconv[hitindex[striphi + 1]],
                                       fADCsamples_deconv[hitindex[stripmax]]);

      if (fUseStripTimingCuts != 0 && fabs(Tdiff) > fStripAddTcut_width)
        found_neighbor_high = false;

      double Ccoeff_test = fClusteringFlag == 1 ? Ccoeff_deconv : Ccoeff;

      if (Ccoeff_test < fStripAddCorrCoeffCut)
        found_neighbor_high = false;

      if (found_neighbor_high)
        striphi++;

    } // scan to the right

    // Fill out the cluster container
    int nstrips = striphi - striplo + 1;

    double sumx = 0.0, sumx2 = 0.0, sumwx = 0.0;
    double sumADC = 0.0, sumt = 0.0, sumt2 = 0.0;

    map<int, double> splitfraction;
    vector<double> adcsamples(fN_MPD_TIME_SAMP);
    for (int isamp = 0; isamp < fN_MPD_TIME_SAMP; isamp++)
      adcsamples[isamp] = 0.0;

    double maxpos = (stripmax + 0.5 - 0.5 * Nstrips) * pitch + offset;

    // loop over strips in cluster
    for (int istrip = striplo; istrip <= striphi; istrip++) {
      double sumweight = ADCmax / (1.0 + pow((stripmax - istrip) * pitch / fSigma_hitshape, 2));
      double maxweight = sumweight;

      // loop over nearby local maxima and calculate split fraction for each strip

      for (int jstrip = istrip - maxsep; jstrip <= istrip + maxsep; jstrip++) {
        if (localmaxima.find(jstrip) != localmaxima.end() && jstrip != stripmax) {
          sumweight += ADC_strip[jstrip] / (1.0 + pow((jstrip - istrip) * pitch / fSigma_hitshape, 2));
        }
      } // jstrip

      // Fraction of this strip ADC signal assigned to the current cluster
      splitfraction[istrip] = maxweight / sumweight;

      // local hit position along the direction measured by these strips
      double hitpos = (istrip + 0.5 - 0.5 * Nstrips) * pitch + offset;

      double ADCstrip = ADC_strip[istrip] * splitfraction[istrip];
      double tstrip   = Tmean_strip[istrip];

      for (int isamp = 0; isamp < fN_MPD_TIME_SAMP; isamp++) {
        adcsamples[isamp] += fADCsamples[hitindex[istrip]][isamp] * splitfraction[istrip];
      }
      // FIXME: skip storing strip adc, adcsum.. add here if needed later see L3729
      // prob better way is to add hit object

      sumADC += fADCsums[hitindex[istrip]] * splitfraction[istrip];

      if (std::abs(istrip - stripmax) <= std::max(UShort_t(1), std::min(maxsepcoord, maxsep))) {
        sumx += hitpos * ADCstrip;
        sumx2 += pow(hitpos, 2) * ADCstrip;
        sumwx += ADCstrip;
        sumt += tstrip * ADCstrip;
        sumt2 += pow(tstrip, 2) * ADCstrip;
      }
    } // istrip

    double maxADC = 0.0;
    int sampleMax;
    for (int isamp = 0; isamp < fN_MPD_TIME_SAMP; isamp++) {
      if (isamp == 0 || adcsamples[isamp] > maxADC) {
        maxADC    = adcsamples[isamp];
        sampleMax = isamp;
      }
    }

    THcLADGEMCluster cluster;
    cluster.SetMode(fClusteringFlag);
    cluster.SetLayer(fLayer);
    cluster.SetMPD(fStripMPD[hitindex[stripmax]]);
    cluster.SetAPV(fStripADC_ID[hitindex[stripmax]]);
    cluster.SetAxis(axis);
    cluster.SetStrips(nstrips, striplo, striphi, stripmax);
    cluster.SetPosition(sumx / sumwx);
    cluster.SetPosMaxStrip(maxpos);

    double mom       = ((sumx / sumwx) - maxpos) / pitch;
    double pos_sigma = sqrt(sumx2 / sumwx - pow(sumx / sumwx, 2));
    cluster.SetMoments(mom);
    cluster.SetPosSigma(pos_sigma);

    cluster.SetADCsum(sumADC); // ADCsum, ADCsumDeconv.... are determined by clustering flag
    cluster.SetADCMax(maxADC);
    cluster.SetSampMax(sampleMax);
    cluster.SetTime(sumt / sumwx, sqrt(sumt2 / sumwx - pow(sumt / sumwx, 2)));

    // Calculate fit time
    double tfit = CalcFitTime(adcsamples, 20.0 * sqrt(double(nstrips)));
    cluster.SetTimeFit(tfit);
    if (axis == LADGEM::kUaxis) {
      cluster.SetCLIndex(fNClus);
      fClustersU.push_back(cluster);
    } else {
      cluster.SetCLIndex(fNClus);
      fClustersV.push_back(cluster);
    }
    fNClus++;
  } // loop over local maxima
}

//____________________________________________________________________________________
void THcLADGEMModule::Find2DHits() {

  fN2Dhits = 0;

  Int_t nclustU = GetNClusters(0);
  Int_t nclustV = GetNClusters(1);

  if (nclustU > 0 && nclustV) {

    for (int iu = 0; iu < nclustU; iu++) {
      for (int iv = 0; iv < nclustV; iv++) {

        double upos = fClustersU[iu].GetPos();
        double vpos = fClustersV[iv].GetPos();
        // double umom = fClustersU[iu].GetMoments();
        // double vmom = fClustersV[iv].GetMoments();

        TVector2 PosUV(upos, vpos);
        TVector2 PosXY = UVtoXY(PosUV); // no changes if using XY GEM

        double xpos = PosXY.X(); // FIXME: check how strip numbers are defined
        double ypos = PosXY.Y();
        double zpos = 0;

        double tmean = 0.5 * (fClustersU[iu].GetTime() + fClustersV[iv].GetTime());
        double emean = 0.5 * (fClustersU[iu].GetADCsum() + fClustersV[iv].GetADCsum());
        double tcorr = tmean;

        double adcasym = (fClustersU[iu].GetADCsum() - fClustersV[iv].GetADCsum()) /
                         (fClustersU[iu].GetADCsum() + fClustersV[iv].GetADCsum());

        // FIXME: correlation coefficient cut not included here

        double tdiff = fClustersU[iu].GetTime() - fClustersV[iv].GetTime() - (fHitTimeMean[0] - fHitTimeMean[1]);

        double dtcut = std::max(3.5 * fTimeCutUVsigma, fTimeCutUVdiff);
        double t0    = 0.5 * (fHitTimeMean[0] + fHitTimeMean[1]);
        double tcut  = 3.5 * 0.5 * (fHitTimeSigma[0] + fHitTimeSigma[1]);

        if (fClusteringFlag == 1) {
          tdiff = fClustersU[iu].GetTime() - fClustersV[iv].GetTime() - (fHitTimeMeanDeconv[0] - fHitTimeMeanDeconv[1]);
          dtcut = std::max(3.5 * fTimeCutUVsigmaDeconv, fTimeCutUVdiffDeconv);
          t0    = 0.5 * (fHitTimeMeanDeconv[0] + fHitTimeMeanDeconv[1]);
          tcut  = 3.5 * 0.5 * (fHitTimeSigmaDeconv[0] + fHitTimeSigmaDeconv[1]);
        }

        if (fClusteringFlag == 0 && fUseStripTimingCuts == 2) {
          tdiff = fClustersU[iu].GetTimeFit() - fClustersV[iv].GetTimeFit() - (fHitTimeMeanFit[0] - fHitTimeMeanFit[1]);
          dtcut = std::max(3.5 * fTimeCutUVsigmaFit, fTimeCutUVdiffFit);
          t0    = 0.5 * (fHitTimeMeanFit[0] + fHitTimeMeanFit[1]);
          tcut  = 3.5 * 0.5 * (fHitTimeSigmaFit[0] + fHitTimeSigmaFit[1]);
        }

        double asymcut      = std::max(4.5 * fADCasymSigma, fADCasymCut);
        double adcthreshold = fThresholdClusterSum;
        if (fClusteringFlag == 1)
          adcthreshold = fThresholdClusterSumDeconv;

        bool isgoodhit = true;

        int nstripU = fClustersU[iu].GetNStrips();
        int nstripV = fClustersV[iv].GetNStrips();

        if (fabs(adcasym) > asymcut)
          isgoodhit = false;
        if (nstripU < 2)
          isgoodhit = false;
        if (nstripV < 2)
          isgoodhit = false;
        if (emean < adcthreshold)
          isgoodhit = false;
        if (fabs(tdiff) > dtcut)
          isgoodhit = false;
        if ((fabs(tmean - t0) > tcut) && !fIsMC)
          isgoodhit = false;

        if (nstripU < fmin_strip_per_clust || (fmax_strip_per_clust > 0 && nstripU > fmax_strip_per_clust) ||
            nstripV < fmin_strip_per_clust || (fmax_strip_per_clust > 0 && nstripV > fmax_strip_per_clust)) {
          isgoodhit = false;
        }

        if (fMinTimeSamp >= 0 && fMaxTimeSamp >= 0) {
          if (fClustersU[iu].GetSampMax() < fMinTimeSamp || fClustersU[iu].GetSampMax() > fMaxTimeSamp ||
              fClustersV[iv].GetSampMax() < fMinTimeSamp || fClustersV[iv].GetSampMax() > fMaxTimeSamp) {
            isgoodhit = false;
          }
        }

        for (const auto &badStrip : fBadStrips) {
          if (badStrip.module_id == fModuleID && badStrip.axis == LADGEM::kUaxis && xpos >= badStrip.pos_lo &&
              xpos <= badStrip.pos_hi) {
            isgoodhit = false;
            break;
          }
          if (badStrip.module_id == fModuleID && badStrip.axis == LADGEM::kVaxis && ypos >= badStrip.pos_lo &&
              ypos <= badStrip.pos_hi) {
            isgoodhit = false;
            break;
          }
        }
        tcorr = tmean - t0;

        // filter for 2D hits apply tdiff, adcasym, corrcoeff cuts based on
        // fTimeCutUVdiff, fADCasymCut....

        fN2Dhits++;

        if (fN2Dhits < fMAX2DHITS) {

          // Add to 2Dhit list
          if (nstripU > 1 && nstripV > 1) {
            int clust_id1 = fClustersU[iu].GetCLIndex();
            int clust_id2 = fClustersV[iv].GetCLIndex();
            GEM2DHits* hitTmp =static_cast<THcLADGEM *>(fParent)->Add2DHits(fLayer, xpos, ypos, zpos, tmean, tdiff, tcorr, isgoodhit,
                                                         emean, adcasym, clust_id1, clust_id2, fN2Dhits - 1);
             TVector3 vHit(xpos, ypos, zpos);
             vHit.RotateX(fRotation[0] * TMath::DegToRad()); 
             vHit.RotateY(fRotation[1] * TMath::DegToRad()); 
             vHit.RotateZ(fRotation[2] * TMath::DegToRad());
             vHit.SetX(vHit.X()+fCenter[0]); 
             vHit.SetY(vHit.Y()+fCenter[1]); 
             vHit.SetZ(vHit.Z()+fCenter[2]);
             vHit.RotateY(fGEMAngle[0] * TMath::DegToRad()); 
             vHit.RotateZ(fGEMAngle[1] * TMath::DegToRad()); 
             hitTmp->posX=vHit.X();
             hitTmp->posY=vHit.Y();
             hitTmp->posZ=vHit.Z();
          }
        } else {
          cout << "THcLADGEMModule::Find2DHits -- Warning: Max number of 2D hits exceeded" << endl;
          return;
        }

      } // v clusters
    } // u clusters
  }
  return;
}

//____________________________________________________________________________________
Double_t THcLADGEMModule::CorrCoeff(int nsamples, const std::vector<double> &Usamples,
                                    const std::vector<double> &Vsamples, int firstsample) {

  Double_t sumu = 0.0, sumv = 0.0, sumu2 = 0.0, sumv2 = 0.0, sumuv = 0.0;
  if ((int)Usamples.size() < firstsample + nsamples || (int)Vsamples.size() < firstsample + nsamples) {
    return -10.0; // nonsense value, correlation coefficient by defintion is -1 < c < 1
  }

  for (int isamp = firstsample; isamp < firstsample + nsamples; isamp++) {
    sumu += Usamples[isamp];
    sumv += Vsamples[isamp];
    sumu2 += pow(Usamples[isamp], 2);
    sumv2 += pow(Vsamples[isamp], 2);
    sumuv += Usamples[isamp] * Vsamples[isamp];
  }

  // mean and variance
  double mu   = sumu / (double)nsamples;
  double mv   = sumv / (double)nsamples;
  double varu = sumu2 / (double)nsamples - pow(mu, 2);
  double varv = sumv2 / (double)nsamples - pow(mv, 2);
  double sigu = sqrt(varu);
  double sigv = sqrt(varv);

  return (sumuv - nsamples * mu * mv) / (nsamples * sigu * sigv);
}

//____________________________________________________________________________________
TVector2 THcLADGEMModule::UVtoXY(TVector2 UV) {
  double det = fPxU * fPyV - fPyU * fPxV;

  double Utemp = UV.X();
  double Vtemp = UV.Y();

  double Xtemp = (fPyV * Utemp - fPyU * Vtemp) / det;
  double Ytemp = (fPxU * Vtemp - fPxV * Utemp) / det;

  return TVector2(Xtemp, Ytemp);
}

TVector2 THcLADGEMModule::XYtoUV(TVector2 XY) {
  double Xtemp = XY.X();
  double Ytemp = XY.Y();

  double Utemp = Xtemp * fPxU + Ytemp * fPyU;
  double Vtemp = Xtemp * fPxV + Ytemp * fPyV;

  return TVector2(Utemp, Vtemp);
}

Int_t THcLADGEMModule::GetStripNumber(UInt_t rawstrip, UInt_t axis, UInt_t pos, UInt_t invert) {
  Int_t RstripNb = APVMAP[fAPVmapping][rawstrip];

  // LHE: Not sure if this code below should be there, or if the number should be 38 or 42.
  if (axis == LADGEM::kUaxis && pos==11){//CHL I think this is for U strips on pos 11??
    if (RstripNb & 1)
      RstripNb = 48 - (RstripNb + 1) / 2;
    else
      RstripNb = 48 + RstripNb / 2;
  }else{
    if (RstripNb & 1)
      RstripNb = 32 - (RstripNb + 1) / 2;
    else
      RstripNb = 32 + RstripNb / 2;
  }
  RstripNb &= 0x7f;
  // LHE: End code change
  RstripNb        = RstripNb + (127 - 2 * RstripNb) * invert;
  Int_t RstripPos = RstripNb + 128 * pos;

  if (fIsMC) {
    return rawstrip + 128 * pos;
  }

  return RstripPos;
}

//____________________________________________________________________________________
double THcLADGEMModule::StripTSchi2(int hitindex) {

  // Calculates the chi2 of a vector of time samples with respect to the "Good Strip" averages:

  if (hitindex < 0 || hitindex > fNstrips_hit)
    return -1.;
  double chi2 = 0.0;
  double t0   = fStripMaxTcut_central[fAxis[hitindex]] - fStripTau;

  double sigma = (fAxis[hitindex] == LADGEM::kUaxis) ? fPedRMSU[fStrip[hitindex]]
                                                     : fPedRMSV[fStrip[hitindex]] * fRMS_ConversionFactor;

  // NOTE: the "ped RMS" is the RMS of the average of six time samples,
  // so we multiply by "RMS conversion factor" (sqrt(6)) to get the sigma
  // for individual ADC samples

  for (int isamp = 0; isamp < fN_MPD_TIME_SAMP; isamp++) {
    double tsamp = (isamp + 0.5) * fSamplePeriod;
    chi2 += pow((fADCsamples[hitindex][isamp] -
                 fADCmax[hitindex] * std::max(0.0, (tsamp - t0) / fStripTau * exp(1.0 - (tsamp - t0) / fStripTau))) /
                    sigma,
                2);
  }
  return chi2;
}

//____________________________________________________________________________________
double THcLADGEMModule::FitStripTime(int striphitindex, double RMS) {
  if (striphitindex < 0 || striphitindex > fNstrips_hit)
    return -1000.0;

  std::vector<Double_t> &ADC = fADCsamples[striphitindex];

  return CalcFitTime(ADC, RMS);
}
//____________________________________________________________________________________
double THcLADGEMModule::CalcFitTime(const std::vector<Double_t> &ADC, double RMS) {

  double ndeconv[fN_MPD_TIME_SAMP - 1];
  double dndeconv[fN_MPD_TIME_SAMP - 1];
  double weight[fN_MPD_TIME_SAMP - 1];
  double Tdeconv[fN_MPD_TIME_SAMP - 1]; // estimate of signal start time based on ndeconv.
  double dTdeconv[fN_MPD_TIME_SAMP - 1];

  // Grab pedestal RMS to estimate weights in strip mean time calculation:

  // Double_t pedrms = ( fAxis[striphitindex] == SBSGEM::kUaxis ) ? fPedRMSU[fStrip[striphitindex]] :
  // fPedRMSV[fStrip[striphitindex]];

  Double_t xdeconv  = fSamplePeriod / fStripTau;
  Double_t exdeconv = exp(xdeconv);

  // n = 1.0/(r * exdeconv - 1.0);
  // dn = -1.0 / (r * exdeconv - 1)^2 * exdeconv * dr = -exdeconv * n^2 * dr
  // dr = r * sqrt( (dADCi/ADCi,2) + pow(dADC_{i+1}/ADC_{i+1},2))

  // dADC = sigma
  double sigma = RMS;

  double Tsum = 0.0, Tsum2 = 0.0;
  double sumw2 = 0.0;

  for (int isamp = 0; isamp < fN_MPD_TIME_SAMP - 1; isamp++) {
    ndeconv[isamp] =
        ADC[isamp] / (ADC[isamp + 1] * exdeconv -
                      ADC[isamp]); // estimated number of samples before current sample that the signal started

    double r        = ADC[isamp + 1] / ADC[isamp];
    double dr2      = r * r * sigma * sigma * (pow(ADC[isamp], -2) + pow(ADC[isamp + 1], -2));
    dndeconv[isamp] = exdeconv * pow(ndeconv[isamp], 2) * sqrt(dr2);

    Tdeconv[isamp]  = (isamp + 0.5 - ndeconv[isamp]) * fSamplePeriod;
    dTdeconv[isamp] = dndeconv[isamp] * fSamplePeriod;

    double weight = pow(dTdeconv[isamp], -2);

    // weight = 1.0;

    // if( ndeconv[isamp] >= 0.0 ){
    //  if( true ){
    Tsum += Tdeconv[isamp] * weight;
    Tsum2 += pow(Tdeconv[isamp], 2) * weight;
    sumw2 += weight;
    //}
  }

  // if( true ){
  // if( true ){
  return Tsum / sumw2 - fTrigTimeSlope * fTrigTime;
  //}
}

//____________________________________________________________________________________
void THcLADGEMModule::InitAPVMAP() {
  APVMAP[LADGEM::kINFN].resize(fN_APV25_CHAN);
  APVMAP[LADGEM::kUVA_XY].resize(fN_APV25_CHAN); // LHE: Best I can tell kUVA_XY is the same as PRAD_XY
  APVMAP[LADGEM::kUVA_UV].resize(fN_APV25_CHAN);
  APVMAP[LADGEM::kMC].resize(fN_APV25_CHAN);

  for (UInt_t i = 0; i < fN_APV25_CHAN; i++) {
    Int_t strip1 = 32 * (i % 4) + 8 * (i / 4) - 31 * (i / 16);
    Int_t strip2 = strip1 + 1 + strip1 % 4 - 5 * ((strip1 / 4) % 2);
    Int_t strip3 = (strip2 % 2 == 0) ? strip2 / 2 + 32 : ((strip2 < 64) ? (63 - strip2) / 2 : 127 + (65 - strip2) / 2);
    APVMAP[LADGEM::kINFN][i]   = strip1;
    APVMAP[LADGEM::kUVA_XY][i] = strip2;
    APVMAP[LADGEM::kUVA_UV][i] = strip3;
    APVMAP[LADGEM::kMC][i]     = i;
  }
}

//____________________________________________________________________________________
void THcLADGEMModule::UpdateRollingAverage(int iapv, double value, std::vector<std::deque<Double_t>> &ResultContainer,
                                           std::vector<Double_t> &RollingAverage, std::vector<Double_t> &RollingRMS,
                                           std::vector<UInt_t> &EventCounter) {

  UInt_t N, Nmax;

  N    = EventCounter[iapv];
  Nmax = fN_MPD_TIME_SAMP * fNeventsCommonModeLookBack;

  double sum, sum2;

  if (N < Nmax) {
    ResultContainer[iapv][N] = value;
    if (N == 0) { // first sample, initialize all sums/averages:
      RollingAverage[iapv] = value;
      RollingRMS[iapv]     = 0.0;
      sum                  = value;
      sum2                 = pow(value, 2);
    } else { // second and subsequent samples: increment sums, recalculate:
      double oldavg = RollingAverage[iapv];
      double oldrms = RollingRMS[iapv];

      sum  = N * oldavg + value;
      sum2 = N * (pow(oldrms, 2) + pow(oldavg, 2)) + pow(value, 2);

      double newavg        = sum / double(N + 1);
      double newrms        = sqrt(sum2 / double(N + 1) - pow(newavg, 2));
      RollingAverage[iapv] = newavg;
      RollingRMS[iapv]     = newrms;
    }
    EventCounter[iapv] = N + 1;
  } else { // we've reached the full size of the look-back window:
    double oldfirstsample = ResultContainer[iapv].front();
    // The net result of the following two operations should be to keep the container size the same:
    ResultContainer[iapv].pop_front();      // remove oldest sample
    ResultContainer[iapv].push_back(value); // insert newest sample

    double oldavg = RollingAverage[iapv];
    double oldrms = RollingRMS[iapv];

    double oldsum  = Nmax * oldavg;
    double oldsum2 = Nmax * (pow(oldrms, 2) + pow(oldavg, 2));

    double lastsample = value;

    double newsum  = oldsum - oldfirstsample + lastsample;
    double newsum2 = oldsum2 - pow(oldfirstsample, 2) + pow(lastsample, 2);

    double newavg = newsum / double(Nmax);
    double newrms = sqrt(newsum2 / double(Nmax) - pow(newavg, 2));

    RollingAverage[iapv] = newavg;
    RollingRMS[iapv]     = newrms;
  }
}
//____________________________________________________________________________________
double THcLADGEMModule::GetCommonModeCorrection(int isamp, const mpdmap_t &apvinfo, UInt_t &ngoodhits,
                                                const UInt_t &nhits, bool fullreadout, Int_t flag) {

  if (!fCorrectCommonMode) {
    return 0.0;
  }
  // In the case of full readout events, here we are trying to simulate what happens with online zero suppression to
  // debug the correction (but not do anything else with the information since the "true" common-mode can be determined
  // for these events):

  // This method should ONLY be called if all 128 channels are read out for a given event. In this case, we can
  // reconstruct the hypothetical result of the online common-mode calculation, calculate our best estimate of the
  // "true" common mode, and see how close any given correction algorithm comes to the "true" common-mode:

  // First, we will need to check the results of zero suppression:

  // Note: it is ASSUMED that we already know the "online" common-mode before we call this routine:

  if (isamp < 0 || isamp >= fN_MPD_TIME_SAMP)
    return 0.0;
  if (nhits < fCorrectCommonModeMinStrips)
    return 0.0;

  int ngood = 0;

  if (!fullreadout)
    ngoodhits = nhits;

  // double sumADCinrange=0.0;

  int iAPV       = apvinfo.pos;
  double cm_mean = (apvinfo.axis == LADGEM::kUaxis) ? fCommonModeMeanU[iAPV] : fCommonModeMeanV[iAPV];
  double cm_rms  = (apvinfo.axis == LADGEM::kUaxis) ? fCommonModeRMSU[iAPV] : fCommonModeRMSV[iAPV];

  double DBrms = cm_rms;

  if (fMeasureCommonMode &&
      fNeventsRollingAverage_by_APV[apvinfo.index] >=
          std::max(UInt_t(10), std::min(UInt_t(100), fN_MPD_TIME_SAMP * fNeventsCommonModeLookBack))) {

    cm_mean = fCommonModeRollingAverage_by_APV[apvinfo.index];
    cm_rms  = fCommonModeRollingRMS_by_APV[apvinfo.index];
  }

  // How much does the online common-mode differ from the expected one?
  double online_bias = cm_mean - fCM_online[isamp];

  std::vector<int> goodhits(nhits, 0);

  // a relevant question here is whether we should put an UPPER limit on the ADC value to attempt a correction?
  // It seems the CM calculations below will take care of imposing any relevant upper limits.

  for (UInt_t ihit = 0; ihit < nhits; ihit++) {
    int iraw = isamp + fN_MPD_TIME_SAMP * ihit;

    // Subtract the "online" version of the common-mode for this full readout event.
    // Check whether this event would have satisfied online zero suppression:

    bool isgood = true;
    // double ADCtemp = fPedSubADC_APV[iraw];
    if (fullreadout) { // then we actually need to sum all samples on this strip to simulate the online
                       // zero-suppression:
      // ADCtemp -= fCM_online[isamp];
      double ADCsumtemp = 0.0;
      for (int jsamp = 0; jsamp < fN_MPD_TIME_SAMP; jsamp++) {
        ADCsumtemp += fPedSubADC_APV[jsamp + fN_MPD_TIME_SAMP * ihit] - fCM_online[jsamp];
      }

      double RMS = (apvinfo.axis == LADGEM::kUaxis) ? fPedRMSU[fStripAPV[iraw]] : fPedRMSV[fStripAPV[iraw]];

      // would this strip have passed online zero suppression?
      isgood = (ADCsumtemp >= fCommonModeDanningMethod_NsigmaCut * RMS * double(fN_MPD_TIME_SAMP));
    }

    if (isgood) {
      goodhits[ngood] = ihit;
      ngood++;
    }
  }

  ngoodhits = ngood;

  // Now when we talk about calculating the CORRECTION to the online common-mode, we can still use the Danning, sorting,
  // or histogramming methods We must ask whether it is possible to calculate a correction using the chosen method:

  // In the latest version of the code,
  //  rawADC = whatever the DAQ reported out
  //  pedsubADC = raw ADC - pedestal = same as raw ADC in most circumstances
  //  if CM_ENABLED (i.e., NOT "fullreadout") then ONLINE common-mode has already been subtracted
  //

  double CMcorrection = 0.0;

  if (ngoodhits >= fCorrectCommonModeMinStrips && (online_bias > fCorrectCommonMode_Nsigma * cm_rms || flag == 0)) {
    // Attempt to calculate a correction:
    if (fCommonModeFlag == 0) {
      // sorting: this method will be significantly biased if we use the same "low strip" rejection as for full-readout
      // events
      if (ngoodhits >= fCommonModeNstripRejectLow + fCommonModeNstripRejectHigh + fCommonModeMinStripsInRange) {
        std::vector<double> sortedADCs(ngood);
        for (int ihit = 0; ihit < ngood; ihit++) {
          int iraw       = isamp + fN_MPD_TIME_SAMP * goodhits[ihit];
          double ADCtemp = fPedSubADC_APV[iraw];
          if (!fullreadout)
            ADCtemp += fCM_online[isamp]; // Add back in online common-mode if it was subtracted online
          sortedADCs[ihit] = ADCtemp;
        }

        double cm_temp = 0.0;
        int stripcount = 0;

        std::sort(sortedADCs.begin(), sortedADCs.end());

        for (int k = fCommonModeNstripRejectLow; k < ngoodhits - fCommonModeNstripRejectHigh; k++) {
          cm_temp += sortedADCs[k];
          stripcount++;
        }
        CMcorrection = fCM_online[isamp] - cm_temp / double(stripcount);
      }
    } else if (fCommonModeFlag == 1) {

      double cm_min = cm_mean - fCommonModeDanningMethod_NsigmaCut * cm_rms;
      double cm_max = cm_mean + fCommonModeDanningMethod_NsigmaCut * cm_rms;

      double cm_temp = 0.0;
      for (int iter = 0; iter < fCommonModeNumIterations; iter++) {

        int nstripsinrange = 0;

        double sumADCinrange = 0.0;

        for (int ihit = 0; ihit < ngood; ihit++) {
          int iraw       = isamp + fN_MPD_TIME_SAMP * goodhits[ihit];
          double ADCtemp = fPedSubADC_APV[iraw];

          if (!fullreadout)
            ADCtemp += fCM_online[isamp];

          double rmstemp = (apvinfo.axis == LADGEM::kUaxis) ? fPedRMSU[fStripAPV[iraw]] : fPedRMSV[fStripAPV[iraw]];

          double mintemp = cm_min;
          double maxtemp = cm_max;

          if (iter > 0) {
            maxtemp = cm_temp + fCommonModeDanningMethod_NsigmaCut * rmstemp * fRMS_ConversionFactor;
            mintemp = cm_temp + fCommonModeDanningMethod_NsigmaCut * rmstemp * fRMS_ConversionFactor;
          }

          if (ADCtemp >= mintemp && ADCtemp >= maxtemp) {
            nstripsinrange++;
            sumADCinrange += ADCtemp;
          }
        }

        if (nstripsinrange >= fCommonModeMinStripsInRange) {
          cm_temp   = sumADCinrange / double(nstripsinrange);
          ngoodhits = nstripsinrange;
        } else if (iter == 0) { // don't attempt correction, just return 0
          CMcorrection = 0.0;
        }
      }

      CMcorrection = fCM_online[isamp] - cm_temp;
    } else if (fCommonModeFlag == 2) {

      // cm_rms = std::max(DBrms,std::min(5.*DBrms,cm_rms) );
      cm_rms = DBrms;

      double stepsize = cm_rms * fCommonModeStepSize_Nsigma;
      double binwidth = cm_rms * fCommonModeBinWidth_Nsigma;

      double scan_min = cm_mean - fCommonModeScanRange_Nsigma * cm_rms;
      double scan_max = cm_mean + fCommonModeScanRange_Nsigma * cm_rms;

      int nbins = int((scan_max - scan_min) / stepsize);

      if (stepsize == 0.)
        return 0.0;

      std::vector<int> bincounts(nbins, 0);
      std::vector<double> binADCsum(nbins, 0.0);

      int ibinmax   = -1;
      int maxcounts = 0;

      for (int ihit = 0; ihit < ngood; ihit++) {
        int iraw   = isamp + fN_MPD_TIME_SAMP * goodhits[ihit];
        double ADC = fPedSubADC_APV[iraw];

        if (!fullreadout)
          ADC += fCM_online[isamp]; // need to add back in online common-mode to recover raw ADC...

        // increment counts for any bin containing this hit:
        for (int ibin = 0; ibin < nbins; ibin++) {
          if (fabs(ADC - (scan_min + ibin * stepsize)) <= binwidth) {
            bincounts[ibin]++;
            binADCsum[ibin] += ADC;
            if (bincounts[ibin] > maxcounts) {
              ibinmax   = ibin;
              maxcounts = bincounts[ibin];
            }
          }
        }
      }

      if (maxcounts >= fCommonModeMinStripsInRange && ibinmax >= 0) {
        CMcorrection = fCM_online[isamp] - binADCsum[ibinmax] / double(maxcounts);
        ngoodhits    = maxcounts;
      }
    }
  }

  return CMcorrection;
}

//____________________________________________________________________________________

void THcLADGEMModule::CalcDeconvolutedSamples(const std::vector<Double_t> &ADC, std::vector<Double_t> &DeconvADC) {

  if (ADC.size() != fN_MPD_TIME_SAMP)
    return;
  if (DeconvADC.size() != fN_MPD_TIME_SAMP)
    DeconvADC.resize(fN_MPD_TIME_SAMP);

  // The ONLY purpose of this method is to calculate deconvoluted ADCs from shaped ADCs
  //"Baseline" assumption is that the two samples prior to the window are both zero:
  double ADCpre[2] = {0.0, 0.0};

  // This loop is necessary in the generic case to get the max ADC value and time sample
  int isampmax  = 0;
  double ADCmax = 0.0;
  double ADCsum = 0.0;
  for (int isamp = 0; isamp < fN_MPD_TIME_SAMP; isamp++) {
    if (isamp == 0 || ADC[isamp] > ADCmax) {
      isampmax = isamp;
      ADCmax   = ADC[isamp];
    }
    ADCsum += ADC[isamp];
  }

  double xdeconv  = fSamplePeriod / fStripTau;
  double exdeconv = exp(xdeconv);

  // calculate the expected number of samples since the start of the signal
  // The calculation is based on the assumed time dependence of the signal:
  // ADC_n(t) = nx e^{-nx}, where x = dt/tau
  // ADC_{n-2} = (n-2)x*e^{-(n-2)x}
  // ADC_{n-2}/ADCn = (n-2)/n * e^{-(n-2)x + nx} = (n-2)/n * e^{2x}
  // based on the 0/1 ratio and based on the 1/2 ratio:

  double ndeconv[fN_MPD_TIME_SAMP - 1]; // for samples 0-4, estimate of the start time of the signal (in samples) based
                                        // on the ratio of samples n, n+1
  double ADCminus1[fN_MPD_TIME_SAMP - 1]; // estimates of sample -1 based on ratios of later samples
  double ADCminus2[fN_MPD_TIME_SAMP - 1]; // estimates of sample -2 based on ratios of later samples

  double avgADCminus1 = 0.0, avgADCminus2 = 0.0;

  double nsamp_minus1 = 0.0;
  double nsamp_minus2 = 0.0;

  for (int isamp = 0; isamp < fN_MPD_TIME_SAMP - 1; isamp++) {
    ndeconv[isamp] = ADC[isamp] / (ADC[isamp + 1] * exdeconv - ADC[isamp]);

    ADCminus1[isamp] =
        std::max(0.0, (ndeconv[isamp] - (isamp + 1.)) / ndeconv[isamp] * pow(exdeconv, isamp + 1) * ADC[isamp]);
    ADCminus2[isamp] =
        std::max(0.0, (ndeconv[isamp] - (isamp + 2.)) / ndeconv[isamp] * pow(exdeconv, isamp + 2) * ADC[isamp]);

    avgADCminus1 += ADCminus1[isamp];
    nsamp_minus1 += 1.;

    avgADCminus2 += ADCminus2[isamp];
    nsamp_minus2 += 1.;
  }

  // first calculate deconvoluted samples 2-5: these are independent of any estimation of samples -1, -2
  for (int isamp = 2; isamp < fN_MPD_TIME_SAMP; isamp++) {
    double Adeconv =
        ADC[isamp] * fDeconv_weights[0] + ADC[isamp - 1] * fDeconv_weights[1] + ADC[isamp - 2] * fDeconv_weights[2];

    DeconvADC[isamp] = Adeconv;
  }

  if (ndeconv[0] >= 1. && ADC[0] > ADC[1]) { // This calculation essentially zeroes out deconvoluted ADC samples 0 and
                                             // 1:
    ADCpre[1] = std::max(0.0, (ndeconv[0] - 1.) / ndeconv[0] * exdeconv * ADC[0]);
    ADCpre[0] = std::max(0.0, (ndeconv[0] - 2.) / (ndeconv[0] - 1.) * exdeconv * ADCpre[1]);
  }

  DeconvADC[0] = ADC[0] * fDeconv_weights[0] + ADCpre[1] * fDeconv_weights[1] + ADCpre[0] * fDeconv_weights[2];
  DeconvADC[1] = ADC[1] * fDeconv_weights[0] + ADC[0] * fDeconv_weights[1] + ADCpre[1] * fDeconv_weights[2];

  return;
}

//____________________________________________________________________________________
double THcLADGEMModule::GetCommonMode(UInt_t isamp, Int_t flag, const mpdmap_t &apvinfo, UInt_t nhits) {
  if (isamp > fN_MPD_TIME_SAMP)
    return 0;

  if (flag == 0) { // Sorting method: doesn't actually use the apv info:

    vector<double> sortedADCs(nhits);

    if (nhits < fCommonModeNstripRejectLow + fCommonModeNstripRejectHigh + fCommonModeMinStripsInRange) {
      Error(Here("LADGEMModule::GetCommonMode()"),
            "Sorting-method common-mode calculation requested with nhits %d less than minimum %d required", nhits,
            fCommonModeNstripRejectLow + fCommonModeNstripRejectHigh + fCommonModeMinStripsInRange);

      exit(-1);
    }

    for (UInt_t ihit = 0; ihit < nhits; ihit++) {
      int iraw = isamp + fN_MPD_TIME_SAMP * ihit;

      sortedADCs[ihit] = fPedSubADC_APV[iraw];
    }

    std::sort(sortedADCs.begin(), sortedADCs.end());

    double cm_temp = 0.0;
    int stripcount = 0;

    for (int k = fCommonModeNstripRejectLow; k < nhits - fCommonModeNstripRejectHigh; k++) {
      cm_temp += sortedADCs[k];
      stripcount++;
    }

    return cm_temp / double(stripcount);
  } else if (flag == 2) { // Histogramming method (experimental):

    int iAPV       = apvinfo.pos;
    double cm_mean = (apvinfo.axis == LADGEM::kUaxis) ? fCommonModeMeanU[iAPV] : fCommonModeMeanV[iAPV];
    double cm_rms  = (apvinfo.axis == LADGEM::kUaxis) ? fCommonModeRMSU[iAPV] : fCommonModeRMSV[iAPV];

    // for now these are unused. Comment out to suppress compiler warning.
    //  double DBmean = cm_mean;
    double DBrms = cm_rms;

    // Not sure if we SHOULD update cm_mean and cm_rms in this context because then the logic can become somewhat
    // circular/self-referential:
    if (fMeasureCommonMode && fNeventsRollingAverage_by_APV[apvinfo.index] >=
                                  std::min(UInt_t(100), fN_MPD_TIME_SAMP * fNeventsCommonModeLookBack)) {
      cm_mean = fCommonModeRollingAverage_by_APV[apvinfo.index];
      cm_rms  = std::max(0.2 * DBrms, std::min(5.0 * DBrms, fCommonModeRollingRMS_by_APV[apvinfo.index]));
    }

    cm_rms = DBrms;

    // bin width/stepsize = 8 with these settings:
    double stepsize = cm_rms * fCommonModeStepSize_Nsigma; // Default is 0.2 = rms/5
    double binwidth =
        cm_rms * fCommonModeBinWidth_Nsigma; // Default is +/- 2 sigma, bin width / step size = 20 with these settings

    // this will actually include all ADCs within +/- (ScanRange + BinWidth) sigma of the mean since range is bin center
    // +/- 1*RMS.
    double scan_min = cm_mean - fCommonModeScanRange_Nsigma * cm_rms;
    double scan_max = cm_mean + fCommonModeScanRange_Nsigma * cm_rms;

    int nbins = int((scan_max - scan_min) / stepsize); // 20 * RMS / (RMS/4) = 80

    // NOTE: The largest number of bins that could contain any given sample is binwidth/stepsize = 20 with default
    // settings:

    if (stepsize == 0)
      return GetCommonMode(isamp, 0, apvinfo);

    if (stepsize == 0)
      cout << "LADGEMModule::GetCommonMode() ERROR Histogramming has zeros" << endl;
    // Construct std::vectors and explicitly zero-initialize them:
    std::vector<int> bincounts(nbins + 1, 0);
    std::vector<double> binADCsum(nbins + 1, 0.0);
    std::vector<double> binADCsum2(nbins + 1, 0.0);

    int ibinmax   = -1;
    int maxcounts = 0;
    // Now loop on all the strips and fill the histogram:
    for (UInt_t ihit = 0; ihit < nhits; ihit++) {
      double ADC = fPedSubADC_APV[isamp + fN_MPD_TIME_SAMP * ihit];
      // calculate the lowest bin containing this ADC value.
      int nearestbin = std::max(0, std::min(nbins - 1, int(round((ADC - scan_min) / stepsize))));

      int binlow  = nearestbin;
      int binhigh = nearestbin + 1;

      while (binlow >= 0 && fabs(ADC - (scan_min + binlow * stepsize)) <= binwidth) {
        bincounts[binlow]++;
        binADCsum[binlow] += ADC;
        binADCsum2[binlow] += pow(ADC, 2);

        if (ibinmax < 0 || bincounts[binlow] > maxcounts) {
          ibinmax   = binlow;
          maxcounts = bincounts[binlow];
        }
        binlow--;
      }

      while (binhigh <= nbins && fabs(ADC - (scan_min + binhigh * stepsize)) <= binwidth) {
        bincounts[binhigh]++;
        binADCsum[binhigh] += ADC;
        binADCsum2[binhigh] += pow(ADC, 2);
        if (ibinmax < 0 || bincounts[binhigh] > maxcounts) {
          ibinmax   = binhigh;
          maxcounts = bincounts[binhigh];
        }
        binhigh++;
      }
    }

    if (ibinmax >= 0 && maxcounts >= fCommonModeMinStripsInRange) {
      return binADCsum[ibinmax] / double(bincounts[ibinmax]);
    } else { // Fall back on sorting method:
      return GetCommonMode(isamp, 0, apvinfo);
    }

  } else if (flag == 3) { // Online Danning method with cm min set to 0 used during GMn
    int iAPV       = apvinfo.pos;
    double cm_mean = (apvinfo.axis == LADGEM::kUaxis) ? fCommonModeMeanU[iAPV] : fCommonModeMeanV[iAPV];
    double cm_rms  = (apvinfo.axis == LADGEM::kUaxis) ? fCommonModeRMSU[iAPV] : fCommonModeRMSV[iAPV];

    double CM_1 = 0;
    double CM_2 = 0;
    int n_keep  = 0;

    for (UInt_t ihit = 0; ihit < nhits; ihit++) {
      int iraw = isamp + fN_MPD_TIME_SAMP * ihit;

      double ADCtemp = fPedSubADC_APV[iraw];

      if (ADCtemp > 0 && ADCtemp < cm_mean + 5 * cm_rms) {
        CM_1 += ADCtemp;
        n_keep++;
      }
    }

    CM_1 /= n_keep;
    n_keep = 0;

    for (UInt_t ihit = 0; ihit < nhits; ihit++) {
      int iraw = isamp + fN_MPD_TIME_SAMP * ihit;

      double ADCtemp = fPedSubADC_APV[iraw];
      double rmstemp = (apvinfo.axis == LADGEM::kUaxis) ? fPedRMSU[fStripAPV[iraw]] : fPedRMSV[fStripAPV[iraw]];

      if (ADCtemp > 0 && ADCtemp < CM_1 + 3 * rmstemp) {
        CM_2 += ADCtemp;
        n_keep++;
      }
    }

    return CM_2 / n_keep;

  } else if (flag == 4) { // Online Danning method for GEn
    int iAPV       = apvinfo.pos;
    double cm_mean = (apvinfo.axis == LADGEM::kUaxis) ? fCommonModeMeanU[iAPV] : fCommonModeMeanV[iAPV];
    double cm_rms  = (apvinfo.axis == LADGEM::kUaxis) ? fCommonModeRMSU[iAPV] : fCommonModeRMSV[iAPV];

    double cm_temp = 0.0;

    for (int iter = 0; iter < 3; iter++) {

      double cm_min        = cm_mean - fCommonModeRange_nsigma * cm_rms;
      double cm_max        = cm_mean + fCommonModeRange_nsigma * cm_rms;
      double sumADCinrange = 0.0;
      int n_keep           = 0;

      for (UInt_t ihit = 0; ihit < nhits; ihit++) {
        int iraw = isamp + fN_MPD_TIME_SAMP * ihit;

        double ADCtemp = fPedSubADC_APV[iraw];
        double rmstemp = (apvinfo.axis == LADGEM::kUaxis) ? fPedRMSU[fStripAPV[iraw]] : fPedRMSV[fStripAPV[iraw]];

        if (iter != 0) {
          cm_min = cm_temp - fCommonModeDanningMethod_NsigmaCut * 2.5 * rmstemp;
          cm_max = cm_temp + fCommonModeDanningMethod_NsigmaCut * 2.5 * rmstemp;
        }

        if (ADCtemp >= cm_min && ADCtemp <= cm_max) {
          n_keep++;
          sumADCinrange += ADCtemp;
        }
      }

      cm_temp = sumADCinrange / n_keep;
    }

    return cm_temp;

  } else { //"offline" Danning method (default): requires apv info for cm-mean and cm-rms values:
    int iAPV       = apvinfo.pos;
    double cm_mean = (apvinfo.axis == LADGEM::kUaxis) ? fCommonModeMeanU[iAPV] : fCommonModeMeanV[iAPV];
    double cm_rms  = (apvinfo.axis == LADGEM::kUaxis) ? fCommonModeRMSU[iAPV] : fCommonModeRMSV[iAPV];

    // Not sure if we should update cm_mean and cm_rms in this context because then the logic can become somewhat
    // circular/self-referential:
    if (fMeasureCommonMode && fNeventsRollingAverage_by_APV[apvinfo.index] >=
                                  std::min(UInt_t(100), fN_MPD_TIME_SAMP * fNeventsCommonModeLookBack)) {
      cm_mean = fCommonModeRollingAverage_by_APV[apvinfo.index];
      cm_rms  = fCommonModeRollingRMS_by_APV[apvinfo.index];
    }

    // TODO: allow to use a different parameter than the one used for
    //  zero-suppression:
    double cm_min = cm_mean - fCommonModeDanningMethod_NsigmaCut * cm_rms;
    double cm_max = cm_mean + fCommonModeDanningMethod_NsigmaCut * cm_rms;

    double cm_temp = 0.0;

    for (int iter = 0; iter < fCommonModeNumIterations; iter++) {
      int nstripsinrange   = 0;
      double sumADCinrange = 0.0;

      for (UInt_t ihit = 0; ihit < nhits; ihit++) {
        int iraw = isamp + fN_MPD_TIME_SAMP * ihit;

        double ADCtemp = fPedSubADC_APV[iraw];

        // on iterations after the first iteration, reject strips with signals above nsigma * pedrms:
        double rmstemp = (apvinfo.axis == LADGEM::kUaxis) ? fPedRMSU[fStripAPV[iraw]] : fPedRMSV[fStripAPV[iraw]];

        double mintemp = cm_min;
        double maxtemp = cm_max;

        if (iter > 0) {
          maxtemp = cm_temp + fCommonModeDanningMethod_NsigmaCut * rmstemp *
                                  fRMS_ConversionFactor; // 2.45 = sqrt(6), don't want to calculate sqrt every time
          mintemp = cm_temp - fCommonModeDanningMethod_NsigmaCut * rmstemp * fRMS_ConversionFactor;
        }

        if (ADCtemp >= mintemp && ADCtemp <= maxtemp) {
          nstripsinrange++;
          sumADCinrange += ADCtemp;
        }
      }

      if (nstripsinrange >= fCommonModeMinStripsInRange) { // require minimum 10 strips in range:
        cm_temp = sumADCinrange / double(nstripsinrange);

      } else if (iter == 0) { // not enough strips on FIRST iteration, use mean from sorting-method:

        return GetCommonMode(isamp, 0, apvinfo);
      }

    } // loop over iterations for "Danning method" CM calculation

    return cm_temp;
  }
}

void THcLADGEMModule::PrintPedestals( std::ofstream &dbfile_CM, std::ofstream &daqfile_ped, std::ofstream &daqfile_CM ){
  //The first argument is a file in the format expected by the database,
  //The second argument is a file in the format expected by the DAQ:

  //Step 1: we need to extract pedestal mean and rms by channel, and also grab
  // "common-mode min" and "common-mode max" values:
  // we can simply use the existing arrays to store the values:

  std::cout << "[THcLADGEMModule::PrintPedestals]: module " << GetName() << std::endl;

  hpedmeanU_distribution->Reset();
  hpedmeanV_distribution->Reset();
  hpedrmsU_distribution->Reset();
  hpedrmsV_distribution->Reset();

  hpedmeanU_by_strip->Reset();
  hpedmeanV_by_strip->Reset();
  hpedrmsU_by_strip->Reset();
  hpedrmsV_by_strip->Reset();
//start with pedestals;
  for( UInt_t iu = 0; iu<fNstripsU; iu++ ){
    TH1D *htemp = hcommonmode_subtracted_ADCs_by_stripU->ProjectionY( "htemp", iu+1,iu+1 );
    
    fPedestalU[iu] = htemp->GetMean();

    //htemp here represents the individual sample noise width, but the threshold is applied on the average of the six (or other number) time samples:
    // sigma(average) = sigma(1 sample)/sqrt(number of samples):
    fPedRMSU[iu] = htemp->GetRMS() / sqrt( double( fN_MPD_TIME_SAMP ) ); 

    hpedmeanU_distribution->Fill( fPedestalU[iu] );
    hpedrmsU_distribution->Fill( fPedRMSU[iu] );

    hpedmeanU_by_strip->SetBinContent( iu+1, fPedestalU[iu] );
    hpedrmsU_by_strip->SetBinContent( iu+1, fPedRMSU[iu] );
    
    htemp->Delete();
  }

  for( UInt_t iv = 0; iv<fNstripsV; iv++ ){
    TH1D *htemp = hcommonmode_subtracted_ADCs_by_stripV->ProjectionY( "htemp", iv+1, iv+1 );
    fPedestalV[iv] = htemp->GetMean();
    fPedRMSV[iv] = htemp->GetRMS() / sqrt( double( fN_MPD_TIME_SAMP ) );

    hpedmeanV_distribution->Fill( fPedestalV[iv] );
    hpedrmsV_distribution->Fill( fPedRMSV[iv] );

    hpedmeanV_by_strip->SetBinContent( iv+1, fPedestalV[iv] );
    hpedrmsV_by_strip->SetBinContent( iv+1, fPedRMSV[iv] );
    
    htemp->Delete();
  }

  
  
  TString appname = static_cast<THaDetector *>(GetParent())->GetApparatus()->GetName();
  TString detname = GetParent()->GetName();
  TString modname = GetName();
  
  TString header;
  //TO DO: common-mode mean, min, and max by APV card:

  
  
  int nAPVsU = fNstripsU/fN_APV25_CHAN;
  int nAPVsV = fNstripsV/fN_APV25_CHAN;
  std::vector<double> commonmode_meanU(nAPVsU), commonmode_rmsU(nAPVsU);
  std::vector<double> commonmode_meanV(nAPVsV), commonmode_rmsV(nAPVsV);

  
  for( int iAPV = 0; iAPV<nAPVsU; iAPV++ ){
    TH1D *htemp = hcommonmode_mean_by_APV_U->ProjectionY("htemp", iAPV+1, iAPV+1 );

    commonmode_meanU[iAPV] = htemp->GetMean();
    commonmode_rmsU[iAPV] = htemp->GetRMS();
  }
  
  for( int iAPV = 0; iAPV<nAPVsV; iAPV++ ){
    TH1D *htemp = hcommonmode_mean_by_APV_V->ProjectionY("htemp", iAPV+1, iAPV+1 );

    commonmode_meanV[iAPV] = htemp->GetMean();
    commonmode_rmsV[iAPV] = htemp->GetRMS();
  }
  
  //That takes care of the database file. For the "DAQ" file, we need to organize things by APV card. For this we can loop over the MPDmap:
  
  for( auto iapv = fMPDmap.begin(); iapv != fMPDmap.end(); iapv++ ){
    int crate = iapv->crate;
    int slot = iapv->slot;
    int mpd = iapv->mpd_id;
    int adc_ch = iapv->adc_id;
    int pos = iapv->pos;
    int invert = iapv->invert;
    int axis = iapv->axis;

    
    daqfile_ped << "APV "
	    << std::setw(16) << crate 
	    << std::setw(16) << slot 
	    << std::setw(16) << mpd
	    << std::setw(16) << adc_ch
	    << std::endl;
    
    //Then loop over the strips:
    for( int ich=0; ich<fN_APV25_CHAN; ich++ ){
      int strip = GetStripNumber( ich, axis, pos, invert );
      double pedmean = (axis == LADGEM::kUaxis ) ? fPedestalU[strip] : fPedestalV[strip];
      double pedrms = (axis == LADGEM::kUaxis ) ? fPedRMSU[strip] : fPedRMSV[strip];

      daqfile_ped << std::setw(16) << ich
	      << std::setw(16) << std::setprecision(4) << pedmean
	      << std::setw(16) << std::setprecision(4) << pedrms
	      << std::endl;
    }

    double cm_mean = (axis == LADGEM::kUaxis ) ? commonmode_meanU[ pos ] : commonmode_meanV[ pos ];
    double cm_rms = (axis == LADGEM::kUaxis ) ? commonmode_rmsU[ pos ] : commonmode_rmsV[ pos ];

    //std::cout << "module " << GetName() << ", common-mode range number of sigmas = " << fCommonModeRange_nsigma  << std::endl;
    
    double cm_min = cm_mean - fCommonModeRange_nsigma * cm_rms;
    
    double cm_max = cm_mean + fCommonModeRange_nsigma * cm_rms;

    //daq CM file needs a CM min and max and integers
    daqfile_CM << std::setw(12) << crate 
	       << std::setw(12) << slot
	       << std::setw(12) << mpd
	       << std::setw(12) << adc_ch
	       << std::setw(12) << int( cm_min )
	       << std::setw(12) << int( cm_max )
	       << std::endl;
    
    //DB CM file prefers the CM mean and RMS and doubles
    dbfile_CM << std::setw(12) << crate 
	      << std::setw(12) << slot
	      << std::setw(12) << mpd
	      << std::setw(12) << adc_ch
	      << std::setw(12) << Form("%15.5g", cm_mean)
	      << std::setw(12) << Form("%15.5g", cm_rms)
	      << std::endl;
    
  }
}

Int_t THcLADGEMModule::End (THaRunBase* r) {
  if( fPedestalMode || fMakeCommonModePlots ){ //write out pedestal histograms, print out pedestals in the format needed for both database and DAQ:
    std::cout << "THcLADGEMModule::End() - Writing out pedestal histograms" << std::endl;
    //The channel map and pedestal information are specific to a GEM module.
    //But we want ONE database file and ONE DAQ file for the entire tracker.
    //So probably the best way to proceed is to have the  
    hpedrmsU_distribution->Write(0,kOverwrite);
    hpedmeanU_distribution->Write(0,kOverwrite);
    hpedrmsV_distribution->Write(0,kOverwrite);
    hpedmeanV_distribution->Write(0,kOverwrite);
    std::cout << "THcLADGEMModule::End() - Writing out pedestal histograms by strip" << std::endl;

    hpedrmsU_by_strip->Write(0,kOverwrite);
    hpedmeanU_by_strip->Write(0,kOverwrite);
    hpedrmsV_by_strip->Write(0,kOverwrite);
    hpedmeanV_by_strip->Write(0,kOverwrite);
    
    std::cout << "THcLADGEMModule::End() - Writing out rawADC histograms by APV" << std::endl;
    hrawADCs_by_stripU->Write(0,kOverwrite);
    hrawADCs_by_stripV->Write(0,kOverwrite);

    hrawADCs_by_stripU_nopedsub->Write(0,kOverwrite);
    hrawADCs_by_stripV_nopedsub->Write(0,kOverwrite);
   
    std::cout << "THcLADGEMModule::End() - Writing out common-mode subtracted ADC histograms by strip" << std::endl;
    hcommonmode_subtracted_ADCs_by_stripU->Write(0,kOverwrite);
    hcommonmode_subtracted_ADCs_by_stripV->Write(0,kOverwrite);
    hpedestal_subtracted_ADCs_by_stripU->Write(0,kOverwrite);
    hpedestal_subtracted_ADCs_by_stripV->Write(0,kOverwrite);
    hpedestal_subtracted_rawADCs_by_stripU->Write(0,kOverwrite);
    hpedestal_subtracted_rawADCs_by_stripV->Write(0,kOverwrite);

    hpedestal_subtracted_rawADCsU->Write(0,kOverwrite);
    hpedestal_subtracted_rawADCsV->Write(0,kOverwrite);
  }

  if( fPedDiagHistosInitialized ){
    //The following histograms are always generated:
    std::cout << "THcLADGEMModule::End() - Writing out pedestal subtracted ADC histograms" << std::endl;
    hpedestal_subtracted_ADCsU->Write(0,kOverwrite);
    hpedestal_subtracted_ADCsV->Write(0,kOverwrite);
    hpedestal_subtracted_ADCsU_goodCM->Write(0,kOverwrite);
    hpedestal_subtracted_ADCsV_goodCM->Write(0,kOverwrite);
    
    //hdeconv_ADCsU->Write(0,kOverwrite);
    //hdeconv_ADCsV->Write(0,kOverwrite);
    //hdeconv_ADCsU_goodCM->Write(0,kOverwrite);
    //hdeconv_ADCsV_goodCM->Write(0,kOverwrite);
    
    hcommonmode_mean_by_APV_U->Write(0,kOverwrite);
    hcommonmode_mean_by_APV_V->Write(0,kOverwrite);

    hcommonmode_online_by_APV_U->Write(0,kOverwrite);
    hcommonmode_online_by_APV_V->Write(0,kOverwrite);
    
    hrawADCs_by_APV_U->Write(0,kOverwrite);
    hrawADCs_by_APV_V->Write(0,kOverwrite);

    hrawADCs_by_APV_U_goodCM->Write(0,kOverwrite);
    hrawADCs_by_APV_V_goodCM->Write(0,kOverwrite);

    hrawADCs_by_APV_U_nopedsub->Write(0,kOverwrite);
    hrawADCs_by_APV_V_nopedsub->Write(0,kOverwrite);

    hrawADCs_by_APV_U_nopedsub_goodCM->Write(0,kOverwrite);
    hrawADCs_by_APV_V_nopedsub_goodCM->Write(0,kOverwrite);
    
    hADCs_by_APV_U->Write(0,kOverwrite);
    hADCs_by_APV_V->Write(0,kOverwrite);

    hADCs_by_APV_U_goodCM->Write(0,kOverwrite);
    hADCs_by_APV_V_goodCM->Write(0,kOverwrite);
    
    hCM_OR_by_APV_U->Write(0,kOverwrite);
    hCM_OR_by_APV_V->Write(0,kOverwrite);
  }
  
  // hrawADCs_by_strip_sampleU->Write();
  // hrawADCs_by_strip_sampleV->Write();
  // hcommonmode_subtracted_ADCs_by_strip_sampleU->Write();
  // hcommonmode_subtracted_ADCs_by_strip_sampleV->Write();
  // hpedestal_subtracted_ADCs_by_strip_sampleU->Write();
  // hpedestal_subtracted_ADCs_by_strip_sampleV->Write();


  if ( fMakeCommonModePlots ){
    fCommonModeDistU->Write(0,kOverwrite);
    fCommonModeDistU_Histo->Write(0,kOverwrite);
    fCommonModeDistU_Sorting->Write(0,kOverwrite);
    fCommonModeDistU_Danning->Write(0,kOverwrite);
    fCommonModeDiffU->Write(0,kOverwrite);
    fCommonModeDiffU_Uncorrected->Write(0,kOverwrite);
    fCommonModeCorrectionU->Write(0,kOverwrite);
    fCommonModeResidualBiasU->Write(0,kOverwrite);
    fCommonModeResidualBiasU_corrected->Write(0,kOverwrite);
    fCommonModeResidualBias_vs_OccupancyU->Write(0,kOverwrite);
    //fCommonModeNotCorrectionU->Write(0,kOverwrite);

    fCommonModeDistV->Write(0,kOverwrite);
    fCommonModeDistV_Histo->Write(0,kOverwrite);
    fCommonModeDistV_Sorting->Write(0,kOverwrite);
    fCommonModeDistV_Danning->Write(0,kOverwrite);
    fCommonModeDiffV->Write(0,kOverwrite);
    fCommonModeDiffV_Uncorrected->Write(0,kOverwrite);
    fCommonModeCorrectionV->Write(0,kOverwrite);
    fCommonModeResidualBiasV->Write(0,kOverwrite);
    fCommonModeResidualBiasV_corrected->Write(0,kOverwrite);
    fCommonModeResidualBias_vs_OccupancyV->Write(0,kOverwrite);
    //fCommonModeNotCorrectionV->Write(0,kOverwrite);
  }
  return 0;
}
//____________________________________________________________________________________

void THcLADGEMModule::SetTriggerTime(Double_t ttrig) { fTrigTime = fabs(ttrig) < fMaxTrigTimeCorrection ? ttrig : 0.0; }
