#include "THcLADHodoPlane.h"
#include "THcGlobals.h"
#include "THcHitList.h"
#include "THcLADHodoHit.h"
#include "THcLADHodoscope.h"
#include "THcParmList.h"
#include "THcRawAdcHit.h"
#include "THcRawTdcHit.h"
#include "THcSignalHit.h"

// This is similar to THcScintillatorPlane from hcana

//_______________________________________________________________________________________
THcLADHodoPlane::THcLADHodoPlane(const char *name, const char *description, const Int_t planenum,
                                 THaDetectorBase *parent)
    : THaSubDetector(name, description, parent) {
  // constructor

  fPlaneNum      = planenum;
  fNScinHits     = 0;
  fNGoodHits     = 0;
  fNScinGoodHits = 0;

  fPosCenter = 0;
  fPosBtm    = 0;
  fPosTop    = 0;

  fHodoHits = new TClonesArray("THcLADHodoHit", MAX_HITS);

  frTopAdcErrorFlag = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcErrorFlag = new TClonesArray("THcSignalHit", MAX_HITS);

  frTopTdcTimeRaw      = new TClonesArray("THcSignalHit", MAX_HITS);
  frTopAdcPedRaw       = new TClonesArray("THcSignalHit", MAX_HITS);
  frTopAdcPulseIntRaw  = new TClonesArray("THcSignalHit", MAX_HITS);
  frTopAdcPulseAmpRaw  = new TClonesArray("THcSignalHit", MAX_HITS);
  frTopAdcPulseTimeRaw = new TClonesArray("THcSignalHit", MAX_HITS);

  frTopTdcTime      = new TClonesArray("THcSignalHit", MAX_HITS);
  frTopAdcPed       = new TClonesArray("THcSignalHit", MAX_HITS);
  frTopAdcPulseInt  = new TClonesArray("THcSignalHit", MAX_HITS);
  frTopAdcPulseAmp  = new TClonesArray("THcSignalHit", MAX_HITS);
  frTopAdcPulseTime = new TClonesArray("THcSignalHit", MAX_HITS);

  frBtmTdcTimeRaw      = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcPedRaw       = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcPulseIntRaw  = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcPulseAmpRaw  = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcPulseTimeRaw = new TClonesArray("THcSignalHit", MAX_HITS);

  frBtmTdcTime      = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcPed       = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcPulseInt  = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcPulseAmp  = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcPulseTime = new TClonesArray("THcSignalHit", MAX_HITS);

  frTopAdcSampPedRaw       = new TClonesArray("THcSignalHit", MAX_HITS);
  frTopAdcSampPulseIntRaw  = new TClonesArray("THcSignalHit", MAX_HITS);
  frTopAdcSampPulseAmpRaw  = new TClonesArray("THcSignalHit", MAX_HITS);
  frTopAdcSampPulseTimeRaw = new TClonesArray("THcSignalHit", MAX_HITS);
  frTopAdcSampPed          = new TClonesArray("THcSignalHit", MAX_HITS);
  frTopAdcSampPulseInt     = new TClonesArray("THcSignalHit", MAX_HITS);
  frTopAdcSampPulseAmp     = new TClonesArray("THcSignalHit", MAX_HITS);
  frTopAdcSampPulseTime    = new TClonesArray("THcSignalHit", MAX_HITS);

  frBtmAdcSampPedRaw       = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcSampPulseIntRaw  = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcSampPulseAmpRaw  = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcSampPulseTimeRaw = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcSampPed          = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcSampPulseInt     = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcSampPulseAmp     = new TClonesArray("THcSignalHit", MAX_HITS);
  frBtmAdcSampPulseTime    = new TClonesArray("THcSignalHit", MAX_HITS);
}

//_______________________________________________________________________________________
THcLADHodoPlane::~THcLADHodoPlane() {
  // destructor
  if (fIsSetup)
    RemoveVariables();

  delete frTopAdcErrorFlag;
  frTopAdcErrorFlag = NULL;
  delete frBtmAdcErrorFlag;
  frBtmAdcErrorFlag = NULL;

  delete fHodoHits;

  delete frTopTdcTimeRaw;
  delete frTopAdcPedRaw;
  delete frTopAdcPulseIntRaw;
  delete frTopAdcPulseAmpRaw;
  delete frTopAdcPulseTimeRaw;

  delete frTopTdcTime;
  delete frTopAdcPed;
  delete frTopAdcPulseInt;
  delete frTopAdcPulseAmp;
  delete frTopAdcPulseTime;

  delete frBtmTdcTimeRaw;
  delete frBtmAdcPedRaw;
  delete frBtmAdcPulseIntRaw;
  delete frBtmAdcPulseAmpRaw;
  delete frBtmAdcPulseTimeRaw;

  delete frBtmTdcTime;
  delete frBtmAdcPed;
  delete frBtmAdcPulseInt;
  delete frBtmAdcPulseAmp;
  delete frBtmAdcPulseTime;

  delete frTopAdcSampPedRaw;
  frTopAdcSampPedRaw = NULL;
  delete frTopAdcSampPulseIntRaw;
  frTopAdcSampPulseIntRaw = NULL;
  delete frTopAdcSampPulseAmpRaw;
  frTopAdcSampPulseAmpRaw = NULL;
  delete frTopAdcSampPulseTimeRaw;
  frTopAdcSampPulseTimeRaw = NULL;
  delete frTopAdcSampPed;
  frTopAdcSampPed = NULL;
  delete frTopAdcSampPulseInt;
  frTopAdcSampPulseInt = NULL;
  delete frTopAdcSampPulseAmp;
  frTopAdcSampPulseAmp = NULL;
  delete frTopAdcSampPulseTime;
  frTopAdcSampPulseTime = NULL;

  delete frBtmAdcSampPedRaw;
  frBtmAdcSampPedRaw = NULL;
  delete frBtmAdcSampPulseIntRaw;
  frBtmAdcSampPulseIntRaw = NULL;
  delete frBtmAdcSampPulseAmpRaw;
  frBtmAdcSampPulseAmpRaw = NULL;
  delete frBtmAdcSampPulseTimeRaw;
  frBtmAdcSampPulseTimeRaw = NULL;
  delete frBtmAdcSampPed;
  frBtmAdcSampPed = NULL;
  delete frBtmAdcSampPulseInt;
  frBtmAdcSampPulseInt = NULL;
  delete frBtmAdcSampPulseAmp;
  frBtmAdcSampPulseAmp = NULL;
  delete frBtmAdcSampPulseTime;
  frBtmAdcSampPulseTime = NULL;

  delete[] fTopPedSum;
  fTopPedSum = 0;
  delete[] fTopPedSum2;
  fTopPedSum2 = 0;
  delete[] fTopPedLimit;
  fTopPedLimit = 0;
  delete[] fTopPedCount;
  fTopPedCount = 0;
  delete[] fBtmPedSum;
  fBtmPedSum = 0;
  delete[] fBtmPedSum2;
  fBtmPedSum2 = 0;
  delete[] fBtmPedLimit;
  fBtmPedLimit = 0;
  delete[] fBtmPedCount;
  fBtmPedCount = 0;
  delete[] fTopPed;
  fTopPed = 0;
  delete[] fBtmPed;
  fBtmPed = 0;
  delete[] fTopThresh;
  fTopThresh = 0;
  delete[] fBtmThresh;
  fBtmThresh = 0;

  delete[] fPosCenter;
  fPosCenter = 0;
  delete[] fPosBtm;
  fPosBtm = 0;
  delete[] fPosTop;
  fPosTop = 0;

  // delete[] fHodoTopMinPh;
  // fHodoTopMinPh = NULL;
  // delete[] fHodoBtmMinPh;
  // fHodoBtmMinPh = NULL;
  // delete[] fHodoTopPhcCoeff;
  // fHodoTopPhcCoeff = NULL;
  // delete[] fHodoBtmPhcCoeff;
  // fHodoBtmPhcCoeff = NULL;
  // delete[] fHodoTopTimeOffset;
  // fHodoTopTimeOffset = NULL;
  // delete[] fHodoBtmTimeOffset;
  // fHodoBtmTimeOffset = NULL;
  delete[] fHodoTopAdcTimeWindowMax;
  fHodoTopAdcTimeWindowMax = NULL;
  delete[] fHodoTopAdcTimeWindowMin;
  fHodoTopAdcTimeWindowMin = NULL;
  delete[] fHodoBtmAdcTimeWindowMax;
  fHodoBtmAdcTimeWindowMax = NULL;
  delete[] fHodoBtmAdcTimeWindowMin;
  fHodoBtmAdcTimeWindowMin = NULL;
  delete[] fHodoVelFit;
  fHodoVelFit = NULL;
  delete[] fHodoCableFit;
  fHodoCableFit = NULL;
  delete[] fHodo_LCoeff;
  fHodo_LCoeff = NULL;
  delete[] fHodoVelFit_FADC;
  fHodoVelFit_FADC = NULL;
  delete[] fHodoCableFit_FADC;
  fHodoCableFit_FADC = NULL;
  delete[] fHodo_LCoeff_FADC;
  fHodo_LCoeff_FADC = NULL;
  delete[] fHodoTop_c1;
  fHodoTop_c1 = NULL;
  delete[] fHodoBtm_c1;
  fHodoBtm_c1 = NULL;
  delete[] fHodoTop_c2;
  fHodoTop_c2 = NULL;
  delete[] fHodoBtm_c2;
  fHodoBtm_c2 = NULL;

  delete[] fHodoVelLight;
  fHodoVelLight = NULL;
  // delete[] fHodoSigma;
  // fHodoSigma = NULL;

  delete[] fEdep2MeV_int;
  fEdep2MeV_int = NULL;
  delete[] fEdep2MeV_amp;
  fEdep2MeV_amp = NULL;
}

//_______________________________________________________________________________________
THaAnalysisObject::EStatus THcLADHodoPlane::Init(const TDatime &date) {

  if (IsZombie())
    return fStatus = kInitError;

  EStatus status;
  if ((status = THaSubDetector::Init(date)))
    return fStatus = status;
  return fStatus = kOK;
}

//_______________________________________________________________________________________
void THcLADHodoPlane::Clear(Option_t *opt) {

  fHodoHits->Clear();

  frTopAdcErrorFlag->Clear();
  frBtmAdcErrorFlag->Clear();

  // Ref time
  fTopTdcRefTime     = kBig;
  fTopAdcRefTime     = kBig;
  fBtmTdcRefTime     = kBig;
  fBtmAdcRefTime     = kBig;
  fTopTdcRefDiffTime = kBig;
  fTopAdcRefDiffTime = kBig;
  fBtmTdcRefDiffTime = kBig;
  fBtmAdcRefDiffTime = kBig;

  // Several arrays, vectors reset here
  frTopTdcTimeRaw->Clear();
  frTopAdcPedRaw->Clear();
  frTopAdcPulseIntRaw->Clear();
  frTopAdcPulseAmpRaw->Clear();
  frTopAdcPulseTimeRaw->Clear();

  frTopTdcTime->Clear();
  frTopAdcPed->Clear();
  frTopAdcPulseInt->Clear();
  frTopAdcPulseAmp->Clear();
  frTopAdcPulseTime->Clear();

  frTopAdcSampPedRaw->Clear();
  frTopAdcSampPulseIntRaw->Clear();
  frTopAdcSampPulseAmpRaw->Clear();
  frTopAdcSampPulseTimeRaw->Clear();
  frTopAdcSampPed->Clear();
  frTopAdcSampPulseInt->Clear();
  frTopAdcSampPulseAmp->Clear();
  frTopAdcSampPulseTime->Clear();

  frBtmTdcTimeRaw->Clear();
  frBtmAdcPedRaw->Clear();
  frBtmAdcPulseIntRaw->Clear();
  frBtmAdcPulseAmpRaw->Clear();
  frBtmAdcPulseTimeRaw->Clear();

  frBtmTdcTime->Clear();
  frBtmAdcPed->Clear();
  frBtmAdcPulseInt->Clear();
  frBtmAdcPulseAmp->Clear();
  frBtmAdcPulseTime->Clear();

  frBtmAdcSampPedRaw->Clear();
  frBtmAdcSampPulseIntRaw->Clear();
  frBtmAdcSampPulseAmpRaw->Clear();
  frBtmAdcSampPulseTimeRaw->Clear();
  frBtmAdcSampPed->Clear();
  frBtmAdcSampPulseInt->Clear();
  frBtmAdcSampPulseAmp->Clear();
  frBtmAdcSampPulseTime->Clear();

  // Waveform vectors
  fTopAdcSampWaveform.clear();
  fBtmAdcSampWaveform.clear();

  fCluster.clear();

  // Reset counters

  // Adc good hits
  fTotNumGoodTopAdcHits = 0;
  fTotNumGoodBtmAdcHits = 0;
  fTotNumGoodAdcHits    = 0;

  fTotNumTopAdcHits = 0;
  fTotNumBtmAdcHits = 0;
  fTotNumAdcHits    = 0;

  // Tdc good hits
  fTotNumGoodTopTdcHits = 0;
  fTotNumGoodBtmTdcHits = 0;
  fTotNumGoodTdcHits    = 0;

  fTotNumTopTdcHits = 0;
  fTotNumBtmTdcHits = 0;
  fTotNumTdcHits    = 0;

  fNScinHits     = 0;
  fNGoodHits     = 0; // looks like not being used.. remove it?
  fNScinGoodHits = 0; // looks like not being used.. remove it?

  // Clear occupancies
  for (UInt_t ielem = 0; ielem < fNumGoodTopAdcHits.size(); ielem++)
    fNumGoodTopAdcHits.at(ielem) = 0;
  for (UInt_t ielem = 0; ielem < fNumGoodBtmAdcHits.size(); ielem++)
    fNumGoodBtmAdcHits.at(ielem) = 0;

  for (UInt_t ielem = 0; ielem < fNumGoodTopTdcHits.size(); ielem++)
    fNumGoodTopTdcHits.at(ielem) = 0;
  for (UInt_t ielem = 0; ielem < fNumGoodBtmTdcHits.size(); ielem++)
    fNumGoodBtmTdcHits.at(ielem) = 0;

  for (UInt_t ielem = 0; ielem < fNumTopAdcHits.size(); ielem++)
    fNumTopAdcHits.at(ielem) = 0;
  for (UInt_t ielem = 0; ielem < fNumBtmAdcHits.size(); ielem++)
    fNumBtmAdcHits.at(ielem) = 0;
  for (UInt_t ielem = 0; ielem < fNumTopTdcHits.size(); ielem++)
    fNumTopTdcHits.at(ielem) = 0;
  for (UInt_t ielem = 0; ielem < fNumBtmTdcHits.size(); ielem++)
    fNumBtmTdcHits.at(ielem) = 0;

  // Clear Ped/Amps/Int/Time
  for (UInt_t ielem = 0; ielem < fGoodTopAdcPed.size(); ielem++) {
    fGoodTopAdcPed.at(ielem)         = 0.0;
    fGoodTopAdcMult.at(ielem)        = 0.0;
    fGoodTopAdcHitUsed.at(ielem)     = 0.0;
    fGoodTopAdcPulseInt.at(ielem)    = 0.0;
    fGoodTopAdcPulseAmp.at(ielem)    = 0.0;
    fGoodTopAdcPulseTime.at(ielem)   = kBig;
    fGoodTopAdcTdcDiffTime.at(ielem) = kBig;
  }
  for (UInt_t ielem = 0; ielem < fGoodBtmAdcPed.size(); ielem++) {
    fGoodBtmAdcPed.at(ielem)         = 0.0;
    fGoodBtmAdcMult.at(ielem)        = 0.0;
    fGoodBtmAdcHitUsed.at(ielem)     = 0.0;
    fGoodBtmAdcPulseInt.at(ielem)    = 0.0;
    fGoodBtmAdcPulseAmp.at(ielem)    = 0.0;
    fGoodBtmAdcPulseTime.at(ielem)   = kBig;
    fGoodBtmAdcTdcDiffTime.at(ielem) = kBig;
  }

  // Clear Good TDC Variables
  for (UInt_t ielem = 0; ielem < fGoodTopTdcTimeUnCorr.size(); ielem++) {
    fGoodTopTdcTimeUnCorr.at(ielem)   = kBig;
    fGoodTopTdcTimeCorr.at(ielem)     = kBig;
    fGoodTopTdcTimeWalkCorr.at(ielem) = kBig;
  }

  for (UInt_t ielem = 0; ielem < fGoodBtmTdcTimeUnCorr.size(); ielem++) {
    fGoodBtmTdcTimeUnCorr.at(ielem)   = kBig;
    fGoodBtmTdcTimeCorr.at(ielem)     = kBig;
    fGoodBtmTdcTimeWalkCorr.at(ielem) = kBig;
  }

  for (UInt_t ielem = 0; ielem < fGoodDiffDistTrack.size(); ielem++) {
    fGoodDiffDistTrack.at(ielem) = kBig;
    fGoodHitTimeDiff.at(ielem)   = kBig;
    fGoodHitTimeAvg.at(ielem)    = kBig;
  }

  fpTime          = -1.e4;
  fHitDistance    = kBig;
  fScinYPos       = kBig;
  fScinXPos       = kBig;
  fTrackXPosition = kBig;
  fTrackYPosition = kBig;
}

//_______________________________________________________________________________________
Int_t THcLADHodoPlane::Decode(const THaEvData &evdata) {
  // Leave it to do nothing

  return 0;
}

//_______________________________________________________________________________________

Int_t THcLADHodoPlane::ReadDatabase(const TDatime &date) {

  // Read database files as needed here
  char prefix[2];
  prefix[0] = tolower(GetParent()->GetPrefix()[0]);
  // Fix to prevent different param files for SHMS & HMS LAD hodoscopes
  prefix[0] = 'l';
  prefix[1] = '\0';

  // Get # of element for each hodo detector
  string parname     = "ladhodo_" + string(GetName()) + "_nr";
  DBRequest list_1[] = {{parname.c_str(), &fNelem, kInt}, {0}};
  gHcParms->LoadParmValues(list_1, prefix);

  delete[] fPosCenter;
  fPosCenter = new Double_t[fNelem];
  delete[] fPosBtm;
  fPosBtm = new Double_t[fNelem];
  delete[] fPosTop;
  fPosTop = new Double_t[fNelem];

  DBRequest list[] = {{Form("ladhodo_%s_zpos", GetName()), &fZpos, kDouble},
                      {Form("ladhodo_%s_dzpos", GetName()), &fDzpos, kDouble},
                      {Form("ladhodo_%s_theta", GetName()), &fTheta, kDouble},
                      {Form("ladhodo_%s_size", GetName()), &fSize, kDouble},
                      {Form("ladhodo_%s_spacing", GetName()), &fSpacing, kDouble},
                      {Form("ladhodo_%s_%s", GetName(), "btm"), fPosBtm, kDouble, fNelem},
                      {Form("ladhodo_%s_%s", GetName(), "top"), fPosTop, kDouble, fNelem},
                      {Form("ladhodo_%s_offset", GetName()), &fPosOffset, kDouble},
                      {Form("ladhodo_%s_center", GetName()), fPosCenter, kDouble, fNelem},
                      {"ladhodo_adc_mode", &fADCMode, kInt, 0, 1},
                      {"ladhodo_adc_diag_cut", &fADCDiagCut, kInt, 0, 1},
                      {"cosmicflag", &fCosmicFlag, kInt, 0, 1},
                      {"ladhodo_debug_adc", &fDebugAdc, kInt, 0, 1},
                      {"ladhodo_SampThreshold", &fSampThreshold, kDouble, 0, 1},
                      {"ladhodo_SampNSA", &fSampNSA, kInt, 0, 1},
                      {"ladhodo_SampNSAT", &fSampNSAT, kInt, 0, 1},
                      {"ladhodo_SampNSB", &fSampNSB, kInt, 0, 1},
                      {"ladhodo_OutputSampWaveform", &fOutputSampWaveform, kInt, 0, 1},
                      {"ladhodo_UseSampWaveform", &fUseSampWaveform, kInt, 0, 1},
                      {"ladhodo_TDC_match_window", &fTDC_match_window, kDouble, 0, 1},
                      {0}};

  // Set Default values

  fDebugAdc           = 0; // Set ADC debug parameter to false unless set in parameter file
  fADCMode            = kADCDynamicPedestal;
  fADCDiagCut         = 50.0;
  fCosmicFlag         = 0;
  fSampThreshold      = 5.;
  fSampNSA            = 0;  // use value stored in event 125 info
  fSampNSB            = 0;  // use value stored in event 125 info
  fSampNSAT           = 2;  // default value in THcRawHit::SetF250Params
  fOutputSampWaveform = 0;  // 0= no output , 1 = output Sample Waveform
  fUseSampWaveform    = 0;  // 0= do not use , 1 = use Sample Waveform
  fTDC_match_window   = 50; //

  gHcParms->LoadParmValues((DBRequest *)&list, prefix);

  DBRequest list5[] = {{"is_mc", &fIsMC, kInt, 0, 1}, {0}};
  fIsMC             = 0;
  gHcParms->LoadParmValues((DBRequest *)&list5, "");

  if (fCosmicFlag == 1)
    cout << " setup for cosmics in scint plane" << endl;

  // Retrieve parameters we need from parent class
  // Common for all planes

  THcLADHodoscope *parent = (THcLADHodoscope *)GetParent();

  fHodoSlop      = parent->GetHodoSlop(fPlaneNum - 1);
  fTdcOffset     = parent->GetTdcOffset(fPlaneNum - 1);
  fAdcTdcOffset  = parent->GetAdcTdcOffset(fPlaneNum - 1);
  fScinTdcMin    = parent->GetTdcMin();
  fScinTdcMax    = parent->GetTdcMax();
  fScinTdcToTime = parent->GetTdcToTime();

  // fTofTolerance=parent->GetTofTolerance();
  // fBetaNominal=parent->GetBetaNominal();
  // fStartTimeCenter=parent->GetStartTimeCenter();
  // fStartTimeSlop=parent->GetStartTimeSlop();

  // Parameters for this plane
  fHodoTopAdcTimeWindowMin = new Double_t[fNelem];
  fHodoBtmAdcTimeWindowMin = new Double_t[fNelem];
  fHodoTopAdcTimeWindowMax = new Double_t[fNelem];
  fHodoBtmAdcTimeWindowMax = new Double_t[fNelem];

  /*
  fHodoTopMinPh = new Double_t[fNelem];
  fHodoBtmMinPh = new Double_t[fNelem];
  fHodoTopPhcCoeff = new Double_t[fNelem];
  fHodoBtmPhcCoeff = new Double_t[fNelem];
  fHodoTopTimeOffset = new Double_t[fNelem];
  fHodoBtmTimeOffset = new Double_t[fNelem];
  */
  fHodoVelLight = new Double_t[fNelem];

  // fHodoSigma = new Double_t[fNelem];

  // New Time-Walk Calibration Parameters
  fHodoVelFit   = new Double_t[fNelem];
  fHodoCableFit = new Double_t[fNelem];
  fHodo_LCoeff  = new Double_t[fNelem];
  fHodoVelFit_FADC = new Double_t[fNelem];
  fHodoCableFit_FADC = new Double_t[fNelem];
  fHodo_LCoeff_FADC  = new Double_t[fNelem];
  fHodoTop_c1   = new Double_t[fNelem];
  fHodoBtm_c1   = new Double_t[fNelem];
  fHodoTop_c2   = new Double_t[fNelem];
  fHodoBtm_c2   = new Double_t[fNelem];

  fEdep2MeV_int = new Double_t[fNelem];
  fEdep2MeV_amp = new Double_t[fNelem];

  for (Int_t j = 0; j < (Int_t)fNelem; j++) {
    Int_t index                 = parent->GetScinIndex(fPlaneNum - 1, j);
    fHodoTopAdcTimeWindowMin[j] = parent->GetHodoTopAdcTimeWindowMin(index);
    fHodoTopAdcTimeWindowMax[j] = parent->GetHodoTopAdcTimeWindowMax(index);
    fHodoBtmAdcTimeWindowMin[j] = parent->GetHodoBtmAdcTimeWindowMin(index);
    fHodoBtmAdcTimeWindowMax[j] = parent->GetHodoBtmAdcTimeWindowMax(index);

    // fHodoTopMinPh[j] = parent->GetHodoTopMinPh(index);
    // fHodoBtmMinPh[j] = parent->GetHodoBtmMinPh(index);
    // fHodoTopPhcCoeff[j] = parent->GetHodoTopPhcCoeff(index);
    // fHodoBtmPhcCoeff[j] = parent->GetHodoBtmPhcCoeff(index);
    // fHodoTopTimeOffset[j] = parent->GetHodoTopTimeOffset(index);
    // fHodoBtmTimeOffset[j] = parent->GetHodoBtmTimeOffset(index);
    fHodoVelLight[j] = parent->GetHodoVelLight(index);

    // Get Time-Walk correction param
    fHodoVelFit[j]   = parent->GetHodoVelFit(index);
    fHodoCableFit[j] = parent->GetHodoCableFit(index);
    fHodo_LCoeff[j]  = parent->GetHodoLCoeff(index);
    fHodoVelFit_FADC[j]   = parent->GetHodoVelFit_FADC(index);
    fHodoCableFit_FADC[j] = parent->GetHodoCableFit_FADC(index);
    fHodo_LCoeff_FADC[j]  = parent->GetHodoLCoeff_FADC(index);
    fHodoTop_c1[j]   = parent->GetHodoTop_c1(index);
    fHodoBtm_c1[j]   = parent->GetHodoBtm_c1(index);
    fHodoTop_c2[j]   = parent->GetHodoTop_c2(index);
    fHodoBtm_c2[j]   = parent->GetHodoBtm_c2(index);
    // Double_t topsigma = parent->GetHodoTopSigma(index);
    // Double_t btmsigma = parent->GetHodoBtmSigma(index);
    // fHodoSigma[j] = TMath::Sqrt(topsigma*topsigma+btmsigma*btmsigma)/2.0;

    fEdep2MeV_int[j] = parent->GetEdep2MeV_int(index);
    fEdep2MeV_amp[j] = parent->GetEdep2MeV_amp(index);
  }

  fTdc_Thrs = parent->GetTDCThrs();

  // Create arrays to hold results here
  InitializePedestals();

  // Initialize

  fNumTopAdcHits = vector<Int_t>(fNelem, 0.0);
  fNumBtmAdcHits = vector<Int_t>(fNelem, 0.0);
  fNumTopTdcHits = vector<Int_t>(fNelem, 0.0);
  fNumBtmTdcHits = vector<Int_t>(fNelem, 0.0);

  fNumGoodTopAdcHits = vector<Int_t>(fNelem, 0.0);
  fNumGoodBtmAdcHits = vector<Int_t>(fNelem, 0.0);
  fNumGoodTopTdcHits = vector<Int_t>(fNelem, 0.0);
  fNumGoodBtmTdcHits = vector<Int_t>(fNelem, 0.0);

  fGoodTopAdcPed         = vector<Double_t>(fNelem, 0.0);
  fGoodBtmAdcPed         = vector<Double_t>(fNelem, 0.0);
  fGoodTopAdcMult        = vector<Double_t>(fNelem, 0.0);
  fGoodBtmAdcMult        = vector<Double_t>(fNelem, 0.0);
  fGoodTopAdcHitUsed     = vector<Double_t>(fNelem, 0.0);
  fGoodBtmAdcHitUsed     = vector<Double_t>(fNelem, 0.0);
  fGoodTopAdcPulseAmp    = vector<Double_t>(fNelem, 0.0);
  fGoodBtmAdcPulseAmp    = vector<Double_t>(fNelem, 0.0);
  fGoodTopAdcPulseInt    = vector<Double_t>(fNelem, 0.0);
  fGoodBtmAdcPulseInt    = vector<Double_t>(fNelem, 0.0);
  fGoodTopAdcPulseTime   = vector<Double_t>(fNelem, 0.0);
  fGoodBtmAdcPulseTime   = vector<Double_t>(fNelem, 0.0);
  fGoodTopAdcTdcDiffTime = vector<Double_t>(fNelem, 0.0);
  fGoodBtmAdcTdcDiffTime = vector<Double_t>(fNelem, 0.0);

  fGoodTopTdcTimeUnCorr   = vector<Double_t>(fNelem, 0.0);
  fGoodBtmTdcTimeUnCorr   = vector<Double_t>(fNelem, 0.0);
  fGoodTopTdcTimeCorr     = vector<Double_t>(fNelem, 0.0);
  fGoodBtmTdcTimeCorr     = vector<Double_t>(fNelem, 0.0);
  fGoodTopTdcTimeWalkCorr = vector<Double_t>(fNelem, 0.0);
  fGoodBtmTdcTimeWalkCorr = vector<Double_t>(fNelem, 0.0);
  fGoodDiffDistTrack      = vector<Double_t>(fNelem, 0.0);
  fGoodHitTimeDiff        = vector<Double_t>(fNelem, 0.0);
  fGoodHitTimeAvg         = vector<Double_t>(fNelem, 0.0);

  return 0;
}

//_______________________________________________________________________________________
Int_t THcLADHodoPlane::DefineVariables(EMode mode) {

  /*
    Initialize global variables for histograms and Root tree
  */

  if (mode == kDefine && fIsSetup)
    return kOK;
  fIsSetup = (mode == kDefine);

  // Register variables in global list

  if (fDebugAdc) {
    RVarDef vars[] = {
        {"TopAdcErrorFlag", "Error Flag for When FPGA Fails", "frTopAdcErrorFlag.THcSignalHit.GetData()"},
        {"BtmAdcErrorFlag", "Error Flag for When FPGA Fails", "frBtmAdcErrorFlag.THcSignalHit.GetData()"},

        {"TopTdcTimeRaw", "List of top raw TDC values.", "frTopTdcTimeRaw.THcSignalHit.GetData()"},
        {"TopAdcPedRaw", "List of top raw ADC pedestals", "frTopAdcPedRaw.THcSignalHit.GetData()"},
        {"TopAdcPulseIntRaw", "List of top raw ADC pulse integrals.", "frTopAdcPulseIntRaw.THcSignalHit.GetData()"},
        {"TopAdcPulseAmpRaw", "List of top raw ADC pulse amplitudes.", "frTopAdcPulseAmpRaw.THcSignalHit.GetData()"},
        {"TopAdcPulseTimeRaw", "List of top raw ADC pulse times.", "frTopAdcPulseTimeRaw.THcSignalHit.GetData()"},

        {"TopTdcTime", "List of top TDC values.", "frTopTdcTime.THcSignalHit.GetData()"},
        {"TopAdcPed", "List of top ADC pedestals", "frTopAdcPed.THcSignalHit.GetData()"},
        {"TopAdcPulseInt", "List of top ADC pulse integrals.", "frTopAdcPulseInt.THcSignalHit.GetData()"},
        {"TopAdcPulseAmp", "List of top ADC pulse amplitudes.", "frTopAdcPulseAmp.THcSignalHit.GetData()"},
        {"TopAdcPulseTime", "List of top ADC pulse times.", "frTopAdcPulseTime.THcSignalHit.GetData()"},

        {"TopAdcSampPedRaw", "Top Raw Samp ADC pedestals", "frTopAdcSampPedRaw.THcSignalHit.GetData()"},
        {"TopAdcSampPulseIntRaw", "Top Raw Samp ADC pulse integrals", "frTopAdcSampPulseIntRaw.THcSignalHit.GetData()"},
        {"TopAdcSampPulseAmpRaw", "Top Raw Samp ADC pulse amplitudes",
         "frTopAdcSampPulseAmpRaw.THcSignalHit.GetData()"},
        {"TopAdcSampPulseTimeRaw", "Top Raw Samp ADC pulse times", "frTopAdcSampPulseTimeRaw.THcSignalHit.GetData()"},
        {"TopAdcSampPed", "Top Samp ADC pedestals", "frTopAdcSampPed.THcSignalHit.GetData()"},
        {"TopAdcSampPulseInt", "Top Samp ADC pulse integrals", "frTopAdcSampPulseInt.THcSignalHit.GetData()"},
        {"TopAdcSampPulseAmp", "Top Samp ADC pulse amplitudes", "frTopAdcSampPulseAmp.THcSignalHit.GetData()"},
        {"TopAdcSampPulseTime", "Top Samp ADC pulse times", "frTopAdcSampPulseTime.THcSignalHit.GetData()"},

        {"BtmTdcTimeRaw", "List of bottom raw TDC values.", "frBtmTdcTimeRaw.THcSignalHit.GetData()"},
        {"BtmAdcPedRaw", "List of bottom raw ADC pedestals", "frBtmAdcPedRaw.THcSignalHit.GetData()"},
        {"BtmAdcPulseIntRaw", "List of bottom raw ADC pulse integrals.", "frBtmAdcPulseIntRaw.THcSignalHit.GetData()"},
        {"BtmAdcPulseAmpRaw", "List of bottom raw ADC pulse amplitudes.", "frBtmAdcPulseAmpRaw.THcSignalHit.GetData()"},
        {"BtmAdcPulseTimeRaw", "List of bottom raw ADC pulse times.", "frBtmAdcPulseTimeRaw.THcSignalHit.GetData()"},

        {"BtmTdcTime", "List of bottom TDC values.", "frBtmTdcTime.THcSignalHit.GetData()"},
        {"BtmAdcPed", "List of bottom ADC pedestals", "frBtmAdcPed.THcSignalHit.GetData()"},
        {"BtmAdcPulseInt", "List of bottom ADC pulse integrals.", "frBtmAdcPulseInt.THcSignalHit.GetData()"},
        {"BtmAdcPulseAmp", "List of bottom ADC pulse amplitudes.", "frBtmAdcPulseAmp.THcSignalHit.GetData()"},
        {"BtmAdcPulseTime", "List of bottom ADC pulse times.", "frBtmAdcPulseTime.THcSignalHit.GetData()"},

        {"BtmAdcSampPedRaw", "Bottom Raw Samp ADC pedestals", "frBtmAdcSampPedRaw.THcSignalHit.GetData()"},
        {"BtmAdcSampPulseIntRaw", "Bottom Raw Samp ADC pulse integrals",
         "frBtmAdcSampPulseIntRaw.THcSignalHit.GetData()"},
        {"BtmAdcSampPulseAmpRaw", "Bottom Raw Samp ADC pulse amplitudes",
         "frBtmAdcSampPulseAmpRaw.THcSignalHit.GetData()"},
        {"BtmAdcSampPulseTimeRaw", "Bottom Raw Samp ADC pulse times",
         "frBtmAdcSampPulseTimeRaw.THcSignalHit.GetData()"},
        {"BtmAdcSampPed", "Bottom Samp ADC pedestals", "frBtmAdcSampPed.THcSignalHit.GetData()"},
        {"BtmAdcSampPulseInt", "Bottom Samp ADC pulse integrals", "frBtmAdcSampPulseInt.THcSignalHit.GetData()"},
        {"BtmAdcSampPulseAmp", "Bottom Samp ADC pulse amplitudes", "frBtmAdcSampPulseAmp.THcSignalHit.GetData()"},
        {"BtmAdcSampPulseTime", "Bottom Samp ADC pulse times", "frBtmAdcSampPulseTime.THcSignalHit.GetData()"},

        {"numTopAdcHits", "Number of Top ADC Hits Per PMT", "fNumTopAdcHits"}, // Hodo+ ADC occupancy - vector<Int_t>
        {"numBtmAdcHits", "Number of Bottom ADC Hits Per PMT",
         "fNumBtmAdcHits"}, // Hodo- ADC occupancy - vector <Int_t>

        {"numTopTdcHits", "Number of Top TDC Hits Per PMT", "fNumTopTdcHits"}, // Hodo+ TDC occupancy - vector<Int_t>
        {"numBtmTdcHits", "Number of Bottom TDC Hits Per PMT",
         "fNumBtmTdcHits"}, // Hodo- TDC occupancy - vector <Int_t>

        {"totNumTopAdcHits", "Total Number of Top ADC Hits", "fTotNumTopAdcHits"}, // Hodo+ raw ADC multiplicity Int_t
        {"totNumBtmAdcHits", "Total Number of Bottom ADC Hits", "fTotNumBtmAdcHits"}, // Hodo- raw ADC multiplicity ""
        {"totNumAdcHits", "Total Number of PMTs Hit (as measured by ADCs)",
         "fTotNumAdcHits"}, // Hodo raw ADC multiplicity  ""

        {"totNumTopTdcHits", "Total Number of Top TDC Hits", "fTotNumTopTdcHits"},    // Hodo+ raw TDC multiplicity ""
        {"totNumBtmTdcHits", "Total Number of Bottom TDC Hits", "fTotNumBtmTdcHits"}, // Hodo- raw TDC multiplicity ""
        {"totNumTdcHits", "Total Number of PMTs Hits (as measured by TDCs)",
         "fTotNumTdcHits"}, // Hodo raw TDC multiplicity  ""
        {0}};
    DefineVarsFromList(vars, mode);
  } // end debug statement

  if (fOutputSampWaveform == 1) {
    RVarDef vars[] = {{"adcBtmSampWaveform", "FADC Btm ADCSample Waveform", "fBtmAdcSampWaveform"},
                      {"adcTopSampWaveform", "FADC Top ADCSample Waveform", "fTopAdcSampWaveform"},
                      {0}};
    DefineVarsFromList(vars, mode);
  }

  // RVarDef track_vars[] = {// Track ID
  //                         // Track Based Beta
  //                         // Delta_transverse
  //                         // Delta_longitudinal
  //                         // Matching HodoHit ID
  //                         {0}};
  // DefineVarsFromList(track_vars, mode);

  RVarDef vars[] = {
      {"nhits", "Number of paddle hits (passed TDC && ADC Min and Max cuts for either end)", "GetNScinHits() "},

      {"TopTdcCounter", "List of top TDC counter numbers.",
       "frTopTdcTimeRaw.THcSignalHit.GetPaddleNumber()"}, // Hodo+ raw TDC occupancy
      {"TopAdcCounter", "List of top ADC counter numbers.",
       "frTopAdcPulseIntRaw.THcSignalHit.GetPaddleNumber()"}, // Hodo+ raw ADC occupancy
      {"BtmTdcCounter", "List of bottom TDC counter numbers.",
       "frBtmTdcTimeRaw.THcSignalHit.GetPaddleNumber()"}, // Hodo- raw TDC occupancy
      {"BtmAdcCounter", "List of bottom ADC counter numbers.",
       "frBtmAdcPulseIntRaw.THcSignalHit.GetPaddleNumber()"}, // Hodo- raw ADC occupancy

      {"numGoodTopAdcHits", "Number of Good Top ADC Hits Per PMT",
       "fNumGoodTopAdcHits"}, // Hodo+ good ADC occupancy - vector<Int_t>
      {"numGoodBtmAdcHits", "Number of Good Bottom ADC Hits Per PMT",
       "fNumGoodBtmAdcHits"}, // Hodo- good ADC occupancy - vector <Int_t>

      {"numGoodTopTdcHits", "Number of Good Top TDC Hits Per PMT",
       "fNumGoodTopTdcHits"}, // Hodo+ good TDC occupancy - vector<Int_t>
      {"numGoodBtmTdcHits", "Number of Good Bottom TDC Hits Per PMT",
       "fNumGoodBtmTdcHits"}, // Hodo- good TDC occupancy - vector <Int_t>

      {"totNumGoodTopAdcHits", "Total Number of Good Top ADC Hits",
       "fTotNumGoodTopAdcHits"}, // Hodo+ good ADC multiplicity - Int_t
      {"totNumGoodBtmAdcHits", "Total Number of Good Bottom ADC Hits",
       "fTotNumGoodBtmAdcHits"}, // Hodo- good ADC multiplicity - Int_t
      {"totNumGoodAdcHits", "TotalNumber of Good ADC Hits Per PMT",
       "fTotNumGoodAdcHits"}, // Hodo good ADC multiplicity - Int_t

      {"totNumGoodTopTdcHits", "Total Number of Good Top TDC Hits",
       "fTotNumGoodTopTdcHits"}, // Hodo+ good TDC multiplicity - Int_t
      {"totNumGoodBtmTdcHits", "Total Number of Good Bottom TDC Hits",
       "fTotNumGoodBtmTdcHits"}, // Hodo- good TDC multiplicity - Int_t
      {"totNumGoodTdcHits", "TotalNumber of Good TDC Hits Per PMT",
       "fTotNumGoodTdcHits"}, // Hodo good TDC multiplicity - Int_t

      {"GoodTopAdcPed", "List of Top ADC pedestals (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodTopAdcPed"}, // vector<Double_t>
      {"GoodBtmAdcPed", "List of Bottom ADC pedestals (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodBtmAdcPed"}, // vector<Double_t>
      {"GoodTopAdcMult", "List of Top ADC Mult (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodTopAdcMult"}, // vector<Double_t>
      {"GoodBtmAdcMult", "List of Bottom ADC Mult (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodBtmAdcMult"}, // vector<Double_t>
      {"GoodTopAdcHitUsed", "List of Top ADC Hit Used (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodTopAdcHitUsed"}, // vector<Double_t>
      {"GoodBtmAdcHitUsed", "List of Bottom ADC Hit Used (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodBtmAdcHitUsed"}, // vector<Double_t>

      {"GoodTopAdcPulseInt", "List of top ADC values (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodTopAdcPulseInt"},
      {"GoodBtmAdcPulseInt", "List of Bottom ADC values (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodBtmAdcPulseInt"},
      {"GoodTopAdcPulseAmp", "List of top ADC peak amp (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodTopAdcPulseAmp"},
      {"GoodBtmAdcPulseAmp", "List of Bottom ADC peak amp (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodBtmAdcPulseAmp"},
      {"GoodTopAdcPulseTime", "List of top ADC time (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodTopAdcPulseTime"},
      {"GoodBtmAdcPulseTime", "List of Bottom ADC time (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodBtmAdcPulseTime"},
      {"GoodTopAdcTdcDiffTime", "List of top TDC - ADC time (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodTopAdcTdcDiffTime"},
      {"GoodBtmAdcTdcDiffTime", "List of Bottom TDC - ADC time (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodBtmAdcTdcDiffTime"},

      {"GoodTopTdcTimeUnCorr", "List of top TDC values (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodTopTdcTimeUnCorr"},
      {"GoodBtmTdcTimeUnCorr", "List of bottom TDC values (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodBtmTdcTimeUnCorr"},
      {"GoodTopTdcTimeCorr", "List of top TDC values (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodTopTdcTimeCorr"},
      {"GoodBtmTdcTimeCorr", "List of bottom TDC values (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodBtmTdcTimeCorr"},
      {"GoodTopTdcTimeWalkCorr", "List of top TDC values (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodTopTdcTimeWalkCorr"},
      {"GoodBtmTdcTimeWalkCorr", "List of bottom TDC values (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodBtmTdcTimeWalkCorr"},
      {"GoodDiffDistTrack", "List of top-bottom TDC values (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodDiffDistTrack"},
      {"GoodHitTimeDiff", "List of top-bottom TDC values (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodHitTimeDiff"},
      {"GoodHitTimeAvg", "List of top-bottom TDC values (passed TDC && ADC Min and Max cuts for either end)",
       "fGoodHitTimeAvg"},
      {"TopTdcRefTime", "Reference time of Top TDC", "fTopTdcRefTime"},
      {"BtmTdcRefTime", "Reference time of Btm TDC", "fBtmTdcRefTime"},
      {"TopAdcRefTime", "Reference time of Top ADC", "fTopAdcRefTime"},
      {"BtmAdcRefTime", "Reference time of Btm aDC", "fBtmAdcRefTime"},
      {"TopTdcRefDiffTime", "Reference Diff time of Top TDC", "fTopTdcRefDiffTime"},
      {"BtmTdcRefDiffTime", "Reference Diff time of Btm TDC", "fBtmTdcRefDiffTime"},
      {"TopAdcRefDiffTime", "Reference Diff time of Top ADC", "fTopAdcRefDiffTime"},
      {"BtmAdcRefDiffTime", "Reference Diff time of Btm aDC", "fBtmAdcRefDiffTime"},

      //{"ngoodhits", "Number of paddle hits (passed tof tolerance and used to determine the focal plane time )",
      //"GetNGoodHits() "},
      {0}};
  DefineVarsFromList(vars, mode);

  // Full Hodohit variables
  RVarDef vars2[] = {
      {"nFullHodoHits", "List of HodoHit objects", "fNScinHits"},
      {"HodoHitTime", "List of HodoHit times", "fHodoHits.THcLADHodoHit.GetScinCorrectedTime()"},
      {"HodoHitTime_FADC", "List of HodoHit times FADC", "fHodoHits.THcLADHodoHit.GetScinCorrectedTime_FADC()"},
      {"HodoHitTOF", "List of HodoHit times of flight", "fHodoHits.THcLADHodoHit.GetScinTOFCorrectedTime()"},
      {"HodoHitPaddleNum", "List of HodoHit paddle numbers", "fHodoHits.THcLADHodoHit.GetPaddleNumber()"},
      {"HodoHitPos", "List of HodoHit positions on bar", "fHodoHits.THcLADHodoHit.GetCalcPosition()"},
      {"HodoHitPos_FADC", "List of HodoHit positions on bar FADC", "fHodoHits.THcLADHodoHit.GetCalcPosition_FADC()"},
      {"HodoHitEdep", "List of HodoHit energy depositions", "fHodoHits.THcLADHodoHit.GetPaddleADC()"},
      {"HodoHitEdepAmp", "List of HodoHit ADC amplitudes", "fHodoHits.THcLADHodoHit.GetPaddleADCpeak()"},
      {"HodoHitEdep_MeV", "List of HodoHit energy depositions", "fHodoHits.THcLADHodoHit.GetPaddleADC_MeV()"},
      {"HodoHitEdepAmp_MeV", "List of HodoHit ADC amplitudes", "fHodoHits.THcLADHodoHit.GetPaddleADCpeak_MeV()"},
      {"HodoHitisFull", "List of Bool's specifing if hit has top and btm hodo hits",
       "fHodoHits.THcLADHodoHit.GetHasCorrectedTimes()"},
       {"HodoHitYPos", "List of HodoHit y positions from track intersection", "fHodoHits.THcLADHodoHit.GetCalcPosition()"},
       {"HodoHitYPos_FADC", "List of HodoHit y positions from track intersection FADC", "fHodoHits.THcLADHodoHit.GetCalcPosition_FADC()"},
      {0}};

  /*
  // cluster variables
  {"NumClus",         "Number of clusters", "fNumberClusters"},
  {"Clus.Pos",        "Position of each paddle clusters", "fCluster.THcScintPlaneCluster.GetClusterPosition()"},
  {"Clus.Size",       "Size of each paddle clusters", "fCluster.THcScintPlaneCluster.GetClusterSize()"},
  {"Clus.Flag",       "Flag of each paddle clusters", "fCluster.THcScintPlaneCluster.GetClusterFlag()"},
  {"Clus.UsedFlag",   "USed Flag of each paddle clusters", "fCluster.THcScintPlaneCluster.GetClusterUsedFlag()"},
  */

  return DefineVarsFromList(vars2, mode);
}

//_______________________________________________________________________________________
Int_t THcLADHodoPlane::ProcessHits(TClonesArray *rawhits, Int_t nexthit) {

  /*! \brief Extract scintillator paddle hits from raw data starting at "nexthit"
   * - Called by THcHodoscope::Decode
   * - Loops through "rawhits" array  starting at index of "nexthit"
   * - Assumes that the hit list is sorted by plane and looping ends when plane number of hit doesn't match fPlaneNum
   * - Fills THcSignalHit objects frTopAdcHits and frBtmAdcHit with pedestal subtracted ADC when value larger than 50
   * - For hits that have TDC value for either top or bottom PMT within  fScinTdcMin and fScinTdcMax
   *  + Creates new  fHodoHits[fNScinHits] =  THcHodoHit
   *  + Calculates pulse height correction to the top and bottom PMT times
   *  + Correct times for time traveled in paddle
   *  + Correct times for time of flight using beta from central spectrometer momentum and particle type
   *  + Calls  SetCorrectedTime method of THcHodoHit
   *  + Increments fNScinHits
   * - Returns value of nexthit + number of hits processed
   *
   */

  // Clear() is being called event by event
  // LHE: Is this true?? I don't think so. Adding Clear() here changes the output (to something that looks reasonable).
  // FixMe
  Clear();

  fTopTdcRefTime     = kBig;
  fTopAdcRefTime     = kBig;
  fBtmTdcRefTime     = kBig;
  fBtmAdcRefTime     = kBig;
  fTopTdcRefDiffTime = kBig;
  fTopAdcRefDiffTime = kBig;
  fBtmTdcRefDiffTime = kBig;
  fBtmAdcRefDiffTime = kBig;
  // counters for Tdc/Adc hits
  UInt_t nrTopTdcHits     = 0;
  UInt_t nrTopAdcHits     = 0;
  UInt_t nrBtmTdcHits     = 0;
  UInt_t nrBtmAdcHits     = 0;
  UInt_t nrSampTopAdcHits = 0;
  UInt_t nrSampBtmAdcHits = 0;
  UInt_t nrTopADCHits     = 0; // Don't really use this. Not sure how it's different from nrTopAdcHits
  UInt_t nrBtmADCHits     = 0; // Don't really use this. Not sure how it's different from nrBtmAdcHits

  Int_t nrawhits = rawhits->GetLast() + 1;
  Int_t ihit     = nexthit;

  // A THcRawHodoHit contains all the information (tdc and adc for both
  // pmts) for a single paddle for a single trigger.  The tdc information
  // might include multiple hits if it uses a multihit tdc.
  // Use "ihit" as the index over THcRawHodoHit objects.  Use
  // "thit" to index over multiple tdc hits within an "ihit".

  Bool_t problem_flag = kFALSE; // check if fTdcRefTime is filled correctly or left initialized (kBig)

  while (ihit < nrawhits) {

    // I think we use THcRawHodoHit as it is
    THcRawHodoHit *hit = (THcRawHodoHit *)rawhits->At(ihit);
    // Pos/Neg in hit class refer to Top/Btm in all LAD classes.

    if (hit->fPlane > fPlaneNum) {
      break;
    }

    Int_t padnum = hit->fCounter;
    Int_t index  = padnum - 1;

    // Top Tdc hits
    THcRawTdcHit &rawTopTdcHit = hit->GetRawTdcHitPos(); // Pos=Top
    if (rawTopTdcHit.GetNHits() > 0 && rawTopTdcHit.HasRefTime()) {
      // if (rawTopTdcHit.GetNHits() > 0) {

      // Assume fTopTdcRefTime is initialized
      if (fTopTdcRefTime == kBig) {
        fTopTdcRefTime     = rawTopTdcHit.GetRefTime();
        fTopTdcRefDiffTime = rawTopTdcHit.GetRefDiffTime();
      }

      // Set problem_flag if the it's not set correctly
      if (fTopTdcRefTime != rawTopTdcHit.GetRefTime()) {
        problem_flag = kTRUE;
      }
    }

    // Loop over multiple tdc hits within ihit
    for (UInt_t thit = 0; thit < rawTopTdcHit.GetNHits(); thit++) {
      ((THcSignalHit *)frTopTdcTimeRaw->ConstructedAt(nrTopTdcHits))->Set(padnum, rawTopTdcHit.GetTimeRaw(thit));
      ((THcSignalHit *)frTopTdcTime->ConstructedAt(nrTopTdcHits))->Set(padnum, rawTopTdcHit.GetTime(thit));

      nrTopTdcHits++; // FIXME: just use thit? or is it used somewhere else too?
      fTotNumTdcHits++;
      fTotNumTopTdcHits++;
      fNumTopTdcHits.at(padnum - 1) = padnum;
    }

    // Now, repeat for the Btm end
    THcRawTdcHit &rawBtmTdcHit = hit->GetRawTdcHitNeg(); // Neg=Btm
    if (rawBtmTdcHit.GetNHits() > 0 && rawBtmTdcHit.HasRefTime()) {
      // if (rawBtmTdcHit.GetNHits() > 0) {

      if (fBtmTdcRefTime == kBig) {
        fBtmTdcRefTime     = rawBtmTdcHit.GetRefTime();
        fBtmTdcRefDiffTime = rawBtmTdcHit.GetRefDiffTime();
      }

      if (fBtmTdcRefTime != rawBtmTdcHit.GetRefTime()) {
        problem_flag = kTRUE;
      }
    }

    for (UInt_t thit = 0; thit < rawBtmTdcHit.GetNHits(); thit++) {
      ((THcSignalHit *)frBtmTdcTimeRaw->ConstructedAt(nrBtmTdcHits))->Set(padnum, rawBtmTdcHit.GetTimeRaw(thit));
      ((THcSignalHit *)frBtmTdcTime->ConstructedAt(nrBtmTdcHits))->Set(padnum, rawBtmTdcHit.GetTime(thit));

      nrBtmTdcHits++; // FIXME: just use thit? or is it used somewhere else too?
      fTotNumTdcHits++;
      fTotNumBtmTdcHits++;
      fNumBtmTdcHits.at(padnum - 1) = padnum;
    } // thit loop

    // Top ADC hits
    THcRawAdcHit &rawTopAdcHit = hit->GetRawAdcHitPos(); // Pos=Top

    if ((rawTopAdcHit.GetNPulses() > 0 || rawTopAdcHit.GetNSamples() > 0) && rawTopAdcHit.HasRefTime()) {
      // if ((rawTopAdcHit.GetNPulses() > 0 || rawTopAdcHit.GetNSamples() > 0)) {

      if (fTopAdcRefTime == kBig) {
        fTopAdcRefTime     = rawTopAdcHit.GetRefTime();
        fTopAdcRefDiffTime = rawTopAdcHit.GetRefDiffTime();
      }

      if (fTopAdcRefTime != rawTopAdcHit.GetRefTime()) {
        problem_flag = kTRUE;
      }
    }

    if (fUseSampWaveform == 0) {

      for (UInt_t thit = 0; thit < rawTopAdcHit.GetNPulses(); thit++) {

        ((THcSignalHit *)frTopAdcPedRaw->ConstructedAt(nrTopAdcHits))->Set(padnum, rawTopAdcHit.GetPedRaw());
        ((THcSignalHit *)frTopAdcPed->ConstructedAt(nrTopAdcHits))->Set(padnum, rawTopAdcHit.GetPed());

        ((THcSignalHit *)frTopAdcPulseIntRaw->ConstructedAt(nrTopAdcHits))
            ->Set(padnum, rawTopAdcHit.GetPulseIntRaw(thit));
        ((THcSignalHit *)frTopAdcPulseInt->ConstructedAt(nrTopAdcHits))->Set(padnum, rawTopAdcHit.GetPulseInt(thit));

        ((THcSignalHit *)frTopAdcPulseAmpRaw->ConstructedAt(nrTopAdcHits))
            ->Set(padnum, rawTopAdcHit.GetPulseAmpRaw(thit));
        ((THcSignalHit *)frTopAdcPulseAmp->ConstructedAt(nrTopAdcHits))->Set(padnum, rawTopAdcHit.GetPulseAmp(thit));

        ((THcSignalHit *)frTopAdcPulseTimeRaw->ConstructedAt(nrTopAdcHits))
            ->Set(padnum, rawTopAdcHit.GetPulseTimeRaw(thit));
        ((THcSignalHit *)frTopAdcPulseTime->ConstructedAt(nrTopAdcHits))
            ->Set(padnum, rawTopAdcHit.GetPulseTime(thit) + fAdcTdcOffset);

        // Error flags 0-2
        if (rawTopAdcHit.GetPulseAmpRaw(thit) > 0)
          ((THcSignalHit *)frTopAdcErrorFlag->ConstructedAt(nrTopAdcHits))->Set(padnum, 0);
        if (rawTopAdcHit.GetPulseAmpRaw(thit) <= 0)
          ((THcSignalHit *)frTopAdcErrorFlag->ConstructedAt(nrTopAdcHits))->Set(padnum, 1);
        if (rawTopAdcHit.GetPulseAmpRaw(thit) <= 0 && rawTopAdcHit.GetNSamples() > 0)
          ((THcSignalHit *)frTopAdcErrorFlag->ConstructedAt(nrTopAdcHits))->Set(padnum, 2);

        nrTopAdcHits++;
        fTotNumAdcHits++;
        fTotNumTopAdcHits++;
        fNumTopAdcHits.at(padnum - 1) = padnum;
      }
    }

    // Use Waveform -- Do we use it? Leave it just in case
    if (rawTopAdcHit.GetNSamples() > 0) {
      rawTopAdcHit.SetSampThreshold(fSampThreshold);
      if (fSampNSA == 0)
        fSampNSA = rawTopAdcHit.GetF250_NSA();
      if (fSampNSB == 0)
        fSampNSB = rawTopAdcHit.GetF250_NSB();

      if (!fIsMC)
        rawTopAdcHit.SetF250Params(fSampNSA, fSampNSB, 4); // Set NPED =4

      if (fSampNSAT != 2)
        rawTopAdcHit.SetSampNSAT(fSampNSAT);
      rawTopAdcHit.SetSampIntTimePedestalPeak();
      fTopAdcSampWaveform.push_back(float(padnum));
      fTopAdcSampWaveform.push_back(float(rawTopAdcHit.GetNSamples()));

      for (UInt_t thit = 0; thit < rawTopAdcHit.GetNSamples(); thit++) {
        fTopAdcSampWaveform.push_back(rawTopAdcHit.GetSample(thit)); // ped subtracted sample (mV)
      }
      for (UInt_t thit = 0; thit < rawTopAdcHit.GetNSampPulses(); thit++) {
        ((THcSignalHit *)frTopAdcSampPedRaw->ConstructedAt(nrSampTopAdcHits))
            ->Set(padnum, rawTopAdcHit.GetSampPedRaw());
        ((THcSignalHit *)frTopAdcSampPed->ConstructedAt(nrSampTopAdcHits))->Set(padnum, rawTopAdcHit.GetSampPed());

        ((THcSignalHit *)frTopAdcSampPulseIntRaw->ConstructedAt(nrSampTopAdcHits))
            ->Set(padnum, rawTopAdcHit.GetSampPulseIntRaw(thit));
        ((THcSignalHit *)frTopAdcSampPulseInt->ConstructedAt(nrSampTopAdcHits))
            ->Set(padnum, rawTopAdcHit.GetSampPulseInt(thit));

        ((THcSignalHit *)frTopAdcSampPulseAmpRaw->ConstructedAt(nrSampTopAdcHits))
            ->Set(padnum, rawTopAdcHit.GetSampPulseAmpRaw(thit));
        ((THcSignalHit *)frTopAdcSampPulseAmp->ConstructedAt(nrSampTopAdcHits))
            ->Set(padnum, rawTopAdcHit.GetSampPulseAmp(thit));

        ((THcSignalHit *)frTopAdcSampPulseTimeRaw->ConstructedAt(nrSampTopAdcHits))
            ->Set(padnum, rawTopAdcHit.GetSampPulseTimeRaw(thit));
        ((THcSignalHit *)frTopAdcSampPulseTime->ConstructedAt(nrSampTopAdcHits))
            ->Set(padnum, rawTopAdcHit.GetSampPulseTime(thit) + fAdcTdcOffset);

        if (rawTopAdcHit.GetNPulses() == 0 || fUseSampWaveform == 1) {
          ((THcSignalHit *)frTopAdcPedRaw->ConstructedAt(nrTopAdcHits))->Set(padnum, rawTopAdcHit.GetSampPedRaw());
          ((THcSignalHit *)frTopAdcPed->ConstructedAt(nrTopAdcHits))->Set(padnum, rawTopAdcHit.GetSampPed());

          ((THcSignalHit *)frTopAdcPulseIntRaw->ConstructedAt(nrTopAdcHits))
              ->Set(padnum, rawTopAdcHit.GetSampPulseIntRaw(thit));
          ((THcSignalHit *)frTopAdcPulseInt->ConstructedAt(nrTopAdcHits))
              ->Set(padnum, rawTopAdcHit.GetSampPulseInt(thit));

          ((THcSignalHit *)frTopAdcPulseAmpRaw->ConstructedAt(nrTopAdcHits))
              ->Set(padnum, rawTopAdcHit.GetSampPulseAmpRaw(thit));
          ((THcSignalHit *)frTopAdcPulseAmp->ConstructedAt(nrTopAdcHits))
              ->Set(padnum, rawTopAdcHit.GetSampPulseAmp(thit));

          ((THcSignalHit *)frTopAdcPulseTimeRaw->ConstructedAt(nrTopAdcHits))
              ->Set(padnum, rawTopAdcHit.GetSampPulseTimeRaw(thit));
          ((THcSignalHit *)frTopAdcPulseTime->ConstructedAt(nrTopAdcHits))
              ->Set(padnum, rawTopAdcHit.GetSampPulseTime(thit) + fAdcTdcOffset);
          ((THcSignalHit *)frTopAdcErrorFlag->ConstructedAt(nrTopAdcHits))->Set(padnum, 3);
          if (fUseSampWaveform == 1)
            ((THcSignalHit *)frTopAdcErrorFlag->ConstructedAt(nrTopAdcHits))->Set(padnum, 0);

          ++nrTopAdcHits;
          fTotNumTopAdcHits++;
          fTotNumAdcHits++;
        }
        ++nrSampTopAdcHits;
      }
    }

    // Btm ADC hits
    THcRawAdcHit &rawBtmAdcHit = hit->GetRawAdcHitNeg(); // Neg=Btm
    if ((rawBtmAdcHit.GetNPulses() > 0 || rawBtmAdcHit.GetNSamples() > 0) && rawBtmAdcHit.HasRefTime()) {
      // if ((rawBtmAdcHit.GetNPulses() > 0 || rawBtmAdcHit.GetNSamples() > 0)) { // Remove RefTime Requirement

      if (fBtmAdcRefTime == kBig) {
        fBtmAdcRefTime     = rawBtmAdcHit.GetRefTime();
        fBtmAdcRefDiffTime = rawBtmAdcHit.GetRefDiffTime();
      }

      if (fBtmAdcRefTime != rawBtmAdcHit.GetRefTime()) {
        problem_flag = kTRUE;
      }
    }

    if (fUseSampWaveform == 0) {

      for (UInt_t thit = 0; thit < rawBtmAdcHit.GetNPulses(); thit++) {

        ((THcSignalHit *)frBtmAdcPedRaw->ConstructedAt(nrBtmAdcHits))->Set(padnum, rawBtmAdcHit.GetPedRaw());
        ((THcSignalHit *)frBtmAdcPed->ConstructedAt(nrBtmAdcHits))->Set(padnum, rawBtmAdcHit.GetPed());

        ((THcSignalHit *)frBtmAdcPulseIntRaw->ConstructedAt(nrBtmAdcHits))
            ->Set(padnum, rawBtmAdcHit.GetPulseIntRaw(thit));
        ((THcSignalHit *)frBtmAdcPulseInt->ConstructedAt(nrBtmAdcHits))->Set(padnum, rawBtmAdcHit.GetPulseInt(thit));

        ((THcSignalHit *)frBtmAdcPulseAmpRaw->ConstructedAt(nrBtmAdcHits))
            ->Set(padnum, rawBtmAdcHit.GetPulseAmpRaw(thit));
        ((THcSignalHit *)frBtmAdcPulseAmp->ConstructedAt(nrBtmAdcHits))->Set(padnum, rawBtmAdcHit.GetPulseAmp(thit));

        ((THcSignalHit *)frBtmAdcPulseTimeRaw->ConstructedAt(nrBtmAdcHits))
            ->Set(padnum, rawBtmAdcHit.GetPulseTimeRaw(thit));
        ((THcSignalHit *)frBtmAdcPulseTime->ConstructedAt(nrBtmAdcHits))
            ->Set(padnum, rawBtmAdcHit.GetPulseTime(thit) + fAdcTdcOffset);

        // Error flags 0-2
        if (rawBtmAdcHit.GetPulseAmpRaw(thit) > 0)
          ((THcSignalHit *)frBtmAdcErrorFlag->ConstructedAt(nrBtmAdcHits))->Set(padnum, 0);
        if (rawBtmAdcHit.GetPulseAmpRaw(thit) <= 0)
          ((THcSignalHit *)frBtmAdcErrorFlag->ConstructedAt(nrBtmAdcHits))->Set(padnum, 1);
        if (rawBtmAdcHit.GetPulseAmpRaw(thit) <= 0 && rawBtmAdcHit.GetNSamples() > 0)
          ((THcSignalHit *)frBtmAdcErrorFlag->ConstructedAt(nrBtmAdcHits))->Set(padnum, 2);

        nrBtmAdcHits++;
        fTotNumAdcHits++;
        fTotNumBtmAdcHits++;
        fNumBtmAdcHits.at(padnum - 1) = padnum;
      }
    }

    // Use Waveform -- Do we use it? Leave it just in case
    if (rawBtmAdcHit.GetNSamples() > 0) {
      rawBtmAdcHit.SetSampThreshold(fSampThreshold);
      if (fSampNSA == 0)
        fSampNSA = rawBtmAdcHit.GetF250_NSA();
      if (fSampNSB == 0)
        fSampNSB = rawBtmAdcHit.GetF250_NSB();

      if (!fIsMC)
        rawBtmAdcHit.SetF250Params(fSampNSA, fSampNSB, 4); // Set NPED =4

      if (fSampNSAT != 2)
        rawBtmAdcHit.SetSampNSAT(fSampNSAT);
      rawBtmAdcHit.SetSampIntTimePedestalPeak();
      fBtmAdcSampWaveform.push_back(float(padnum));
      fBtmAdcSampWaveform.push_back(float(rawBtmAdcHit.GetNSamples()));

      for (UInt_t thit = 0; thit < rawBtmAdcHit.GetNSamples(); thit++) {
        fBtmAdcSampWaveform.push_back(rawBtmAdcHit.GetSample(thit)); // ped subtracted sample (mV)
      }
      for (UInt_t thit = 0; thit < rawBtmAdcHit.GetNSampPulses(); thit++) {
        ((THcSignalHit *)frBtmAdcSampPedRaw->ConstructedAt(nrSampBtmAdcHits))
            ->Set(padnum, rawBtmAdcHit.GetSampPedRaw());
        ((THcSignalHit *)frBtmAdcSampPed->ConstructedAt(nrSampBtmAdcHits))->Set(padnum, rawBtmAdcHit.GetSampPed());

        ((THcSignalHit *)frBtmAdcSampPulseIntRaw->ConstructedAt(nrSampBtmAdcHits))
            ->Set(padnum, rawBtmAdcHit.GetSampPulseIntRaw(thit));
        ((THcSignalHit *)frBtmAdcSampPulseInt->ConstructedAt(nrSampBtmAdcHits))
            ->Set(padnum, rawBtmAdcHit.GetSampPulseInt(thit));

        ((THcSignalHit *)frBtmAdcSampPulseAmpRaw->ConstructedAt(nrSampBtmAdcHits))
            ->Set(padnum, rawBtmAdcHit.GetSampPulseAmpRaw(thit));
        ((THcSignalHit *)frBtmAdcSampPulseAmp->ConstructedAt(nrSampBtmAdcHits))
            ->Set(padnum, rawBtmAdcHit.GetSampPulseAmp(thit));

        ((THcSignalHit *)frBtmAdcSampPulseTimeRaw->ConstructedAt(nrSampBtmAdcHits))
            ->Set(padnum, rawBtmAdcHit.GetSampPulseTimeRaw(thit));
        ((THcSignalHit *)frBtmAdcSampPulseTime->ConstructedAt(nrSampBtmAdcHits))
            ->Set(padnum, rawBtmAdcHit.GetSampPulseTime(thit) + fAdcTdcOffset);

        if (rawBtmAdcHit.GetNPulses() == 0 || fUseSampWaveform == 1) {
          ((THcSignalHit *)frBtmAdcPedRaw->ConstructedAt(nrBtmAdcHits))->Set(padnum, rawBtmAdcHit.GetSampPedRaw());
          ((THcSignalHit *)frBtmAdcPed->ConstructedAt(nrBtmAdcHits))->Set(padnum, rawBtmAdcHit.GetSampPed());

          ((THcSignalHit *)frBtmAdcPulseIntRaw->ConstructedAt(nrBtmAdcHits))
              ->Set(padnum, rawBtmAdcHit.GetSampPulseIntRaw(thit));
          ((THcSignalHit *)frBtmAdcPulseInt->ConstructedAt(nrBtmAdcHits))
              ->Set(padnum, rawBtmAdcHit.GetSampPulseInt(thit));

          ((THcSignalHit *)frBtmAdcPulseAmpRaw->ConstructedAt(nrBtmAdcHits))
              ->Set(padnum, rawBtmAdcHit.GetSampPulseAmpRaw(thit));
          ((THcSignalHit *)frBtmAdcPulseAmp->ConstructedAt(nrBtmAdcHits))
              ->Set(padnum, rawBtmAdcHit.GetSampPulseAmp(thit));

          ((THcSignalHit *)frBtmAdcPulseTimeRaw->ConstructedAt(nrBtmAdcHits))
              ->Set(padnum, rawBtmAdcHit.GetSampPulseTimeRaw(thit));
          ((THcSignalHit *)frBtmAdcPulseTime->ConstructedAt(nrBtmAdcHits))
              ->Set(padnum, rawBtmAdcHit.GetSampPulseTime(thit) + fAdcTdcOffset);
          ((THcSignalHit *)frBtmAdcErrorFlag->ConstructedAt(nrBtmAdcHits))->Set(padnum, 3);
          if (fUseSampWaveform == 1)
            ((THcSignalHit *)frBtmAdcErrorFlag->ConstructedAt(nrBtmAdcHits))->Set(padnum, 0);

          ++nrBtmAdcHits;
          fTotNumBtmAdcHits++;
          fTotNumAdcHits++;
        }
        ++nrSampBtmAdcHits;
      }
    }

    // LHE. Accept more than one TDC hit per paddle
    // Determine good tdc top and btm times
    Int_t n_top_hits_max = std::max(rawTopTdcHit.GetNHits(), rawTopAdcHit.GetNSampPulses());
    Int_t n_btm_hits_max = std::max(rawBtmTdcHit.GetNHits(), rawBtmAdcHit.GetNSampPulses());

    std::vector<Bool_t> btdcraw_top(n_top_hits_max, kFALSE);
    std::vector<Bool_t> btdcraw_btm(n_btm_hits_max, kFALSE);
    std::vector<Int_t> tdc_top(n_top_hits_max, -999);
    std::vector<Int_t> tdc_btm(n_btm_hits_max, -999);
    std::vector<Double_t> good_ielem_TopTdc(n_top_hits_max, -1);
    std::vector<Double_t> good_ielem_BtmTdc(n_btm_hits_max, -1);

    // Good TDC Top times
    int n_top_hit = 0;
    int tdc_top_tmp;
    for (UInt_t thit = 0; thit < hit->GetRawTdcHitPos().GetNHits(); thit++) { // Pos=Top

      tdc_top_tmp = hit->GetRawTdcHitPos().GetTime(thit) + fTdcOffset; // Pos=Top

      if (tdc_top_tmp >= fScinTdcMin && tdc_top_tmp <= fScinTdcMax) {
        tdc_top[n_top_hit]           = tdc_top_tmp;
        btdcraw_top[n_top_hit]       = kTRUE;
        good_ielem_TopTdc[n_top_hit] = n_top_hit;
        n_top_hit++;
      }
    }

    // Good TDC Btm times
    int n_btm_hit = 0;
    int tdc_btm_tmp;
    for (UInt_t thit = 0; thit < hit->GetRawTdcHitNeg().GetNHits(); thit++) { // Neg=Btm

      tdc_btm_tmp = hit->GetRawTdcHitNeg().GetTime(thit) + fTdcOffset; // Neg=Btm

      if (tdc_btm_tmp >= fScinTdcMin && tdc_btm_tmp <= fScinTdcMax) {
        tdc_btm[n_btm_hit]           = tdc_btm_tmp;
        btdcraw_btm[n_btm_hit]       = kTRUE;
        good_ielem_BtmTdc[n_btm_hit] = n_btm_hit;
        n_btm_hit++;
      }
    }

    // Adc btm
    std::vector<Bool_t> badcraw_btm(n_btm_hit, kFALSE);
    std::vector<Double_t> adcped_btm(n_btm_hit, -999);
    std::vector<Int_t> adcmult_btm(n_btm_hit, 0);
    std::vector<Int_t> adchitused_btm(n_btm_hit, 0);
    std::vector<Double_t> adcint_btm(n_btm_hit, -999);
    std::vector<Double_t> adcamp_btm(n_btm_hit, -kBig);
    std::vector<Double_t> adctime_btm(n_btm_hit, kBig);
    std::vector<Double_t> adctdcdifftime_btm(n_btm_hit, kBig);
    std::vector<Double_t> good_ielem_BtmAdc_maxAmp(n_btm_hit, -1);
    std::vector<Double_t> good_ielem_BtmAdc_minTime(n_btm_hit, -1);

    // Adc top
    std::vector<Bool_t> badcraw_top(n_top_hit, kFALSE);
    std::vector<Double_t> adcped_top(n_top_hit, -999);
    std::vector<Int_t> adcmult_top(n_top_hit, 0);
    std::vector<Int_t> adchitused_top(n_top_hit, 0);
    std::vector<Double_t> adcint_top(n_top_hit, -999);
    std::vector<Double_t> adcamp_top(n_top_hit, -kBig);
    std::vector<Double_t> adctime_top(n_top_hit, kBig);
    std::vector<Double_t> adctdcdifftime_top(n_top_hit, kBig);
    std::vector<Double_t> good_ielem_TopAdc_maxAmp(n_top_hit, -1);
    std::vector<Double_t> good_ielem_TopAdc_minTime(n_top_hit, -1);

    if (fADCMode == kADCDynamicPedestal) {

      // Loop Here over all hits per event for btm side of plane
      // Try to find the max pulseAmp and corresponding ielem
      // within the given time window for TdcAdcTimeDiff
      for (int i_btm_hit = 0; i_btm_hit < n_btm_hit; i_btm_hit++) {

        if (good_ielem_BtmTdc[i_btm_hit] == -1)
          continue;
        Double_t max_adcamp_test     = -1000.;
        Double_t max_adctdcdiff_test = 1000.;

        for (UInt_t ielem = 0; ielem < rawBtmAdcHit.GetNPulses(); ielem++) {

          Double_t pulseAmp       = rawBtmAdcHit.GetPulseAmp(ielem);
          Double_t pulseTime      = rawBtmAdcHit.GetPulseTime(ielem) + fAdcTdcOffset;
          Double_t TdcAdcTimeDiff = tdc_btm[i_btm_hit] * fScinTdcToTime - pulseTime;

          if (rawBtmAdcHit.GetPulseAmpRaw(ielem) <= 0)
            pulseAmp = -200.; // do we want to to this? or skip simply this element?

          Bool_t pulseTimeCut =
              (TdcAdcTimeDiff > fHodoBtmAdcTimeWindowMin[index]) && (TdcAdcTimeDiff < fHodoBtmAdcTimeWindowMax[index]);
          Bool_t is_used       = false;
          Bool_t is_used_test2 = false;
          for (UInt_t i = 0; i < i_btm_hit; i++) {
            if (good_ielem_BtmAdc_maxAmp[i] == ielem) {
              is_used = true;
            }
            if (good_ielem_BtmAdc_minTime[i] == ielem) {
              is_used_test2 = true;
            }
          }
          if (!is_used && pulseTimeCut && pulseAmp > max_adcamp_test) {
            good_ielem_BtmAdc_maxAmp[i_btm_hit] = ielem;
            max_adcamp_test                     = pulseAmp;
          }
          if (!is_used_test2 && abs(TdcAdcTimeDiff) < max_adctdcdiff_test) {
            good_ielem_BtmAdc_minTime[i_btm_hit] = ielem;
            max_adctdcdiff_test                  = abs(TdcAdcTimeDiff);
          }
        }

        // good_ielem_BtmAdc_maxAmp: select the pulse that has TdcAdcTimeDiff within the given range and with maximum
        // pulseAmp good_ielem_BtmAdc_minTime: selects the pulse with minimum value of TdcAdcTimeDiff
        if (good_ielem_BtmAdc_maxAmp[i_btm_hit] == -1 && good_ielem_BtmAdc_minTime[i_btm_hit] != -1)
          good_ielem_BtmAdc_maxAmp[i_btm_hit] = good_ielem_BtmAdc_minTime[i_btm_hit];
        if (good_ielem_BtmAdc_maxAmp[i_btm_hit] == -1 && good_ielem_BtmAdc_minTime[i_btm_hit] == -1 &&
            rawBtmAdcHit.GetNPulses() > 0)
          good_ielem_BtmAdc_maxAmp[i_btm_hit] = 0;

        if (good_ielem_BtmAdc_maxAmp[i_btm_hit] != -1 &&
            good_ielem_BtmAdc_maxAmp[i_btm_hit] < rawBtmAdcHit.GetNPulses()) {
          adcped_btm[i_btm_hit]     = rawBtmAdcHit.GetPed();
          adcmult_btm[i_btm_hit]    = rawBtmAdcHit.GetNPulses();
          adchitused_btm[i_btm_hit] = good_ielem_BtmAdc_maxAmp[i_btm_hit] + 1;
          adcint_btm[i_btm_hit]     = rawBtmAdcHit.GetPulseInt(good_ielem_BtmAdc_maxAmp[i_btm_hit]);
          adcamp_btm[i_btm_hit]     = rawBtmAdcHit.GetPulseAmp(good_ielem_BtmAdc_maxAmp[i_btm_hit]);
          if (rawBtmAdcHit.GetPulseAmpRaw(good_ielem_BtmAdc_maxAmp[i_btm_hit]) <= 0)
            adcamp_btm[i_btm_hit] = -200.;
          adctime_btm[i_btm_hit] = rawBtmAdcHit.GetPulseTime(good_ielem_BtmAdc_maxAmp[i_btm_hit]) + fAdcTdcOffset;
          badcraw_btm[i_btm_hit] = kTRUE;
          adctdcdifftime_btm[i_btm_hit] = tdc_btm[i_btm_hit] * fScinTdcToTime - adctime_btm[i_btm_hit];
        }

        // Loop Here over all hits per event for top side of plane
        // Try to find the max pulseAmp and corresponding ielem
        // within the given time window for TdcAdcTimeDiff
        for (int i_top_hit = 0; i_top_hit < n_top_hit; i_top_hit++) {

          if (good_ielem_TopTdc[i_top_hit] == -1)
            continue;
          Double_t max_adcamp_test     = -1000.;
          Double_t max_adctdcdiff_test = 1000.;

          for (UInt_t ielem = 0; ielem < rawTopAdcHit.GetNPulses(); ielem++) {

            Double_t pulseAmp       = rawTopAdcHit.GetPulseAmp(ielem);
            Double_t pulseTime      = rawTopAdcHit.GetPulseTime(ielem) + fAdcTdcOffset;
            Double_t TdcAdcTimeDiff = tdc_top[i_top_hit] * fScinTdcToTime - pulseTime;

            if (rawTopAdcHit.GetPulseAmpRaw(ielem) <= 0)
              pulseAmp = -200.; // do we want to to this? or skip simply this element?

            Bool_t pulseTimeCut = (TdcAdcTimeDiff > fHodoTopAdcTimeWindowMin[index]) &&
                                  (TdcAdcTimeDiff < fHodoTopAdcTimeWindowMax[index]);
            Bool_t is_used       = false;
            Bool_t is_used_test2 = false;
            for (UInt_t i = 0; i < i_top_hit; i++) {
              if (good_ielem_TopAdc_maxAmp[i] == ielem) {
                is_used = true;
              }
              if (good_ielem_TopAdc_minTime[i] == ielem) {
                is_used_test2 = true;
              }
            }
            if (!is_used && pulseTimeCut && pulseAmp > max_adcamp_test) {
              good_ielem_TopAdc_maxAmp[i_top_hit] = ielem;
              max_adcamp_test                     = pulseAmp;
            }
            if (!is_used_test2 && abs(TdcAdcTimeDiff) < max_adctdcdiff_test) {
              good_ielem_TopAdc_minTime[i_top_hit] = ielem;
              max_adctdcdiff_test                  = abs(TdcAdcTimeDiff);
            }
          }

          // good_ielem_TopAdc_maxAmp: select the pulse that has TdcAdcTimeDiff within the given range and with maximum
          // pulseAmp good_ielem_TopAdc_minTime: selects the pulse with minimum value of TdcAdcTimeDiff
          if (good_ielem_TopAdc_maxAmp[i_top_hit] == -1 && good_ielem_TopAdc_minTime[i_top_hit] != -1)
            good_ielem_TopAdc_maxAmp[i_top_hit] = good_ielem_TopAdc_minTime[i_top_hit];
          if (good_ielem_TopAdc_maxAmp[i_top_hit] == -1 && good_ielem_TopAdc_minTime[i_top_hit] == -1 &&
              rawTopAdcHit.GetNPulses() > 0)
            good_ielem_TopAdc_maxAmp[i_top_hit] = 0;

          if (good_ielem_TopAdc_maxAmp[i_top_hit] != -1 &&
              good_ielem_TopAdc_maxAmp[i_top_hit] < rawTopAdcHit.GetNPulses()) {
            adcped_top[i_top_hit]     = rawTopAdcHit.GetPed();
            adcmult_top[i_top_hit]    = rawTopAdcHit.GetNPulses();
            adchitused_top[i_top_hit] = good_ielem_TopAdc_maxAmp[i_top_hit] + 1;
            adcint_top[i_top_hit]     = rawTopAdcHit.GetPulseInt(good_ielem_TopAdc_maxAmp[i_top_hit]);
            adcamp_top[i_top_hit]     = rawTopAdcHit.GetPulseAmp(good_ielem_TopAdc_maxAmp[i_top_hit]);
            if (rawTopAdcHit.GetPulseAmpRaw(good_ielem_TopAdc_maxAmp[i_top_hit]) <= 0)
              adcamp_top[i_top_hit] = -200.;
            adctime_top[i_top_hit] = rawTopAdcHit.GetPulseTime(good_ielem_TopAdc_maxAmp[i_top_hit]) + fAdcTdcOffset;
            badcraw_top[i_top_hit] = kTRUE;
            adctdcdifftime_top[i_top_hit] = tdc_top[i_top_hit] * fScinTdcToTime - adctime_top[i_top_hit];
          }
        }

      } // End loop over all good Btm TDC hits
    } // if kADCDynamicPedestal

    // Currently not set up to use kADCSampleIntegral or kADCSampIntDynPed
    /*
        else if (fADCMode == kADCSampleIntegral) {
          adcint_top  = hit->GetRawAdcHitPos().GetSampleIntRaw() - fTopPed[index]; // Pos=Top
          adcint_btm  = hit->GetRawAdcHitNeg().GetSampleIntRaw() - fBtmPed[index]; // Neg=Btm
          badcraw_top = badcraw_btm = kTRUE;

        } else if (fADCMode == kADCSampIntDynPed) {
          adcint_top  = hit->GetRawAdcHitPos().GetSampleInt(); // Pos=Top
          adcint_btm  = hit->GetRawAdcHitNeg().GetSampleInt(); // Neg=Btm
          badcraw_top = badcraw_btm = kTRUE;

        } else {
          // same as using kADCSampleIntegral
          adcint_top  = hit->GetRawAdcHitPos().GetPulseIntRaw() - fTopPed[index]; // Pos=Top
          adcint_btm  = hit->GetRawAdcHitNeg().GetPulseIntRaw() - fBtmPed[index]; // Neg=Btm
          badcraw_top = badcraw_btm = kTRUE;
        }
    */
    // Loop over individual TDC/ADC pairs for top and btm

    for (int i_top_hit = 0; i_top_hit < n_top_hit; i_top_hit++) {

      if (good_ielem_TopAdc_maxAmp[i_top_hit] != -1) {

        // good adc multiplicity
        fTotNumGoodTopAdcHits++;
        fTotNumGoodAdcHits++;

        // good adc occupancy
        fNumGoodTopAdcHits.at(padnum - 1)     = padnum;
        fGoodTopAdcPed.at(padnum - 1)         = adcped_top[i_top_hit];
        fGoodTopAdcMult.at(padnum - 1)        = adcmult_top[i_top_hit];
        fGoodTopAdcHitUsed.at(padnum - 1)     = adchitused_top[i_top_hit];
        fGoodTopAdcPulseInt.at(padnum - 1)    = adcint_top[i_top_hit];
        fGoodTopAdcPulseAmp.at(padnum - 1)    = adcamp_top[i_top_hit];
        fGoodTopAdcPulseTime.at(padnum - 1)   = adctime_top[i_top_hit];
        fGoodTopAdcTdcDiffTime.at(padnum - 1) = adctdcdifftime_top[i_top_hit];
      }

      // DEFINE THE "GOOD +TDC Multiplicities and Occupancies"
      if (good_ielem_TopTdc[i_top_hit] != -1) {
        fTotNumGoodTopTdcHits++;
        fTotNumGoodTdcHits++;
        // good tdc occupancy
        fNumGoodTopTdcHits.at(padnum - 1) = padnum;
      }
    }

    for (int i_btm_hit = 0; i_btm_hit < n_btm_hit; i_btm_hit++) {

      if (good_ielem_BtmAdc_maxAmp[i_btm_hit] != -1) {

        // good adc multiplicity
        fTotNumGoodBtmAdcHits++;
        fTotNumGoodAdcHits++;

        // good adc occupancy
        fNumGoodBtmAdcHits.at(padnum - 1)     = padnum;
        fGoodBtmAdcPed.at(padnum - 1)         = adcped_btm[i_btm_hit];
        fGoodBtmAdcMult.at(padnum - 1)        = adcmult_btm[i_btm_hit];
        fGoodBtmAdcHitUsed.at(padnum - 1)     = adchitused_btm[i_btm_hit];
        fGoodBtmAdcPulseInt.at(padnum - 1)    = adcint_btm[i_btm_hit];
        fGoodBtmAdcPulseAmp.at(padnum - 1)    = adcamp_btm[i_btm_hit];
        fGoodBtmAdcPulseTime.at(padnum - 1)   = adctime_btm[i_btm_hit];
        fGoodBtmAdcTdcDiffTime.at(padnum - 1) = adctdcdifftime_btm[i_btm_hit];
      }

      if (good_ielem_BtmTdc[i_btm_hit] != -1) {
        fTotNumGoodBtmTdcHits++;
        fTotNumGoodTdcHits++;
        // good tdc occupancy
        fNumGoodBtmTdcHits.at(padnum - 1) = padnum;
      }
    }

    // Match top and Bottom TDC/ADC Pairs
    vector<int> good_ielem_Top_Tdc_Full(n_top_hit, -1);
    vector<int> good_ielem_Btm_Tdc_Full(n_btm_hit, -1);
    // vector<int> good_ielem_Top_Adc_Full(n_top_hit, -1);
    // vector<int> good_ielem_Btm_Adc_Full(n_btm_hit, -1);
    int n_top_hit_full = 0;
    for (int i_top_hit = 0; i_top_hit < n_top_hit; i_top_hit++) {
      double dt = kBig;
      double good_i_elem_top_tmp, good_i_elem_btm_tmp;
      for (int i_btm_hit = 0; i_btm_hit < n_btm_hit; i_btm_hit++) {
        if (good_ielem_TopAdc_maxAmp[i_top_hit] != -1 && good_ielem_BtmAdc_maxAmp[i_btm_hit] != -1) {
          if (std::find(good_ielem_Top_Tdc_Full.begin(), good_ielem_Top_Tdc_Full.end(), good_ielem_TopTdc[i_top_hit]) !=
              good_ielem_Top_Tdc_Full.end()) {
            continue;
          }
          if (std::find(good_ielem_Btm_Tdc_Full.begin(), good_ielem_Btm_Tdc_Full.end(), good_ielem_BtmTdc[i_btm_hit]) !=
              good_ielem_Btm_Tdc_Full.end()) {
            continue;
          }
          double dt = tdc_top[i_top_hit] * fScinTdcToTime - tdc_btm[i_btm_hit] * fScinTdcToTime;
          if (abs(dt) < fTDC_match_window) { // Hard coded to 50ns
            good_ielem_Top_Tdc_Full[n_top_hit_full] = good_ielem_TopTdc[i_top_hit];
            good_ielem_Btm_Tdc_Full[n_top_hit_full] = good_ielem_BtmTdc[i_btm_hit];
            // good_ielem_Top_Adc_Full[n_top_hit_full] = good_ielem_TopAdc_maxAmp[i_top_hit];
            // good_ielem_Btm_Adc_Full[n_top_hit_full] = good_ielem_BtmAdc_maxAmp[i_btm_hit];
            n_top_hit_full++;
          }
        }
      }
      if (dt < fTDC_match_window) { // Hard coded to 50ns
        good_ielem_Top_Tdc_Full[n_top_hit_full] = good_i_elem_top_tmp;
        good_ielem_Btm_Tdc_Full[n_top_hit_full] = good_i_elem_btm_tmp;
        n_top_hit_full++;
      }
    }

    //////////////////////////////////////////////////////
    // All full hits matched by here. Now perform time walk, and other calibrations.
    //////////////////////////////////////////////////////
    for (int i_good_hit = 0; i_good_hit < n_top_hit_full; i_good_hit++) {
      if (good_ielem_Top_Tdc_Full[i_good_hit] != -1 && good_ielem_Btm_Tdc_Full[i_good_hit] != -1) {
        int i_good_top_tdc_elem = good_ielem_Top_Tdc_Full[i_good_hit];
        int i_good_btm_tdc_elem = good_ielem_Btm_Tdc_Full[i_good_hit];
        int i_good_top_adc_elem = good_ielem_Top_Tdc_Full[i_good_hit];
        int i_good_btm_adc_elem = good_ielem_Btm_Tdc_Full[i_good_hit];

        // DEFINE THE "GOOD -TDC  Multiplicities and Occupancies"

        new ((*fHodoHits)[fNScinHits])
            THcLADHodoHit(tdc_top[i_good_top_tdc_elem], tdc_btm[i_good_btm_tdc_elem], adcint_top[i_good_top_adc_elem],
                          adcint_btm[i_good_btm_adc_elem], hit->fCounter, this);

        ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetTopADCpeak(adcamp_top[i_good_top_adc_elem]);
        ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetBtmADCpeak(adcamp_btm[i_good_btm_adc_elem]);
        ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetTopADCtime(adctime_top[i_good_top_adc_elem]);
        ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetBtmADCtime(adctime_btm[i_good_btm_adc_elem]);

        // Calculate Time-Walk Correction

        // LHE: Start section 1. Need to rewrite this later
        // Define GoodTdcUnCorrTime
        if (btdcraw_top[i_good_top_tdc_elem] && badcraw_top[i_good_top_adc_elem]) {
          fGoodTopTdcTimeUnCorr.at(padnum - 1) = tdc_top[i_good_top_tdc_elem] * fScinTdcToTime;

          tw_corr_top = 1. / pow(adcamp_top[i_good_top_adc_elem] / fTdc_Thrs, fHodoTop_c2[padnum - 1]) -
                        1. / pow(200. / fTdc_Thrs, fHodoTop_c2[padnum - 1]);

          fGoodTopTdcTimeWalkCorr.at(padnum - 1) = tdc_top[i_good_top_tdc_elem] * fScinTdcToTime - tw_corr_top;
        }
        if (btdcraw_btm[i_good_btm_tdc_elem] && badcraw_btm[i_good_btm_adc_elem]) {
          fGoodBtmTdcTimeUnCorr.at(padnum - 1) = tdc_btm[i_good_btm_tdc_elem] * fScinTdcToTime;

          // tw_corr_btm = fHodoBtm_c1[padnum-1]/pow(adcamp_btm/fTdc_Thrs,fHodoBtm_c2[padnum-1]) -
          // fHodoBtm_c1[padnum-1]/pow(200./fTdc_Thrs, fHodoBtm_c2[padnum-1]);

          tw_corr_btm = 1. / pow(adcamp_btm[i_good_btm_adc_elem] / fTdc_Thrs, fHodoBtm_c2[padnum - 1]) -
                        1. / pow(200. / fTdc_Thrs, fHodoBtm_c2[padnum - 1]);

          fGoodBtmTdcTimeWalkCorr.at(padnum - 1) = tdc_btm[i_good_btm_tdc_elem] * fScinTdcToTime - tw_corr_btm;
        }
        // LHE: End section 1. Need to rewrite this later

        // Do corrections if valid TDC on both ends of bar (both ends should always be valid)
        if ((btdcraw_top[i_good_top_tdc_elem] && btdcraw_btm[i_good_btm_tdc_elem]) &&
            (badcraw_top[i_good_top_adc_elem] && badcraw_btm[i_good_btm_adc_elem])) {
          // Do the pulse height correction to the time.  (Position dependent corrections later)
          Double_t adc_timec_top = adctime_top[i_good_top_adc_elem];
          Double_t adc_timec_btm = adctime_btm[i_good_btm_adc_elem];
          Double_t timec_top, timec_btm;
          // FADC style.
          timec_top = tdc_top[i_good_top_tdc_elem] * fScinTdcToTime - tw_corr_top + fHodo_LCoeff[index];
          timec_btm = tdc_btm[i_good_btm_tdc_elem] * fScinTdcToTime - tw_corr_btm - 2 * fHodoCableFit[index] +
                      fHodo_LCoeff[index];
          adc_timec_top = adc_timec_top  + fHodo_LCoeff_FADC[index]; // Pretty sure we don't need to apply time walk corrections. Not sure why the original code did this
          adc_timec_btm = adc_timec_btm - 2 * fHodoCableFit_FADC[index] + fHodo_LCoeff_FADC[index];

          Double_t TWCorrDiff = fGoodBtmTdcTimeWalkCorr.at(padnum - 1) - 2 * fHodoCableFit[index] -
                                fGoodTopTdcTimeWalkCorr.at(padnum - 1);

          Double_t fHitDistCorr = 0.5 * TWCorrDiff * fHodoVelFit[index];

          fGoodDiffDistTrack.at(index) = fHitDistCorr;
          fGoodHitTimeDiff.at(index)   = TWCorrDiff;

          Double_t vellight = fHodoVelLight[index]; // read from hodo_cuts.param, where it is set fixed to 15.0

          Double_t dist_from_center = 0.5 * (timec_btm - timec_top) * vellight;
          Double_t scint_center     = 0.5 * (fPosBtm[padnum - 1] + fPosTop[padnum - 1]);
          Double_t hit_position     = scint_center + dist_from_center;
          hit_position              = TMath::Min(hit_position, fPosBtm[padnum - 1]);
          hit_position              = TMath::Max(hit_position, fPosTop[padnum - 1]);
          Double_t scin_corrected_time;

          scin_corrected_time = 0.5 * (timec_btm + timec_top);
          // timec_top                   = scin_corrected_time;
          // timec_btm                   = scin_corrected_time;
          Double_t adc_time_corrected = 0.5 * (adc_timec_top + adc_timec_btm);
          Double_t fHitDistCorr_FADC = 0.5 * (adc_timec_btm - adc_timec_top) * fHodoVelFit_FADC[index];
          // LHE: Need to change this to match LAD
          // if (fCosmicFlag) {
          //   toptime     = timec_top + (fZpos + (index % 2) * fDzpos) / (29.979 * fBetaNominal);
          //   btmtime     = timec_btm + (fZpos + (index % 2) * fDzpos) / (29.979 * fBetaNominal);
          //   adc_toptime = adc_time_corrected + (fZpos + (index % 2) * fDzpos) / (29.979 * fBetaNominal);
          //   adc_btmtime = adc_time_corrected + (fZpos + (index % 2) * fDzpos) / (29.979 * fBetaNominal);
          // } else {
          //   toptime     = timec_top - (fZpos + (index % 2) * fDzpos) / (29.979 * fBetaNominal);
          //   btmtime     = timec_btm - (fZpos + (index % 2) * fDzpos) / (29.979 * fBetaNominal);
          //   adc_toptime = adc_time_corrected - (fZpos + (index % 2) * fDzpos) / (29.979 * fBetaNominal);
          //   adc_btmtime = adc_time_corrected - (fZpos + (index % 2) * fDzpos) / (29.979 * fBetaNominal);
          // }

          ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetPaddleCenter(fPosCenter[index]);
          ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetCorrectedTimes(timec_top, timec_btm, scin_corrected_time);
          ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetCorrectedTimes_FADC(adc_timec_top, adc_timec_btm, adc_time_corrected);
          ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetTopADCpeak(adcamp_top[i_good_top_adc_elem]);
          ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetBtmADCpeak(adcamp_btm[i_good_btm_adc_elem]);
          ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetCalcPosition(fHitDistCorr);
          ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetCalcPosition_FADC(fHitDistCorr_FADC);

          fGoodTopTdcTimeCorr.at(padnum - 1) = timec_top;
          fGoodBtmTdcTimeCorr.at(padnum - 1) = timec_btm;
          fGoodHitTimeAvg.at(padnum - 1)     = scin_corrected_time;
        } else {
          // This condition should never happen
          cout << "THcLADHodoPlane::ProcessHits. Good hit has bad ADC or TDC on some end. Contact expert" << endl;

          /*
          Double_t timec_top, timec_btm;
          timec_top = tdc_top;
          timec_btm = tdc_btm;
          if (btdcraw_top && badcraw_top) {
            // FADC style. Removed fTofUsingInvAdc
            timec_top = tdc_top * fScinTdcToTime - tw_corr_top + fHodo_LCoeff[index];
          }
          if (btdcraw_btm && badcraw_btm) {
            // FADC style. Removed fTofUsingInvAdc
            timec_btm = tdc_btm * fScinTdcToTime - tw_corr_btm - 2 * fHodoCableFit[index] + fHodo_LCoeff[index];
          }
          Double_t adc_btm = 0., adc_top = 0.;
          if (badcraw_btm)
            adc_btm = adcamp_btm;
          if (badcraw_top)
            adc_top = adcamp_top;
          ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetPaddleCenter(fPosCenter[index]);
          ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetCorrectedTimes(timec_top, timec_btm);
          ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetBtmADCpeak(adc_btm); // needed for new TWCOrr
          ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetTopADCpeak(adc_top); // needed for new TWCOrr
          if (badcraw_btm) {
            ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetBtmADCtime(adctime_btm);
          } else {
            ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetBtmADCtime(-999.);
          }
          if (badcraw_top) {
            ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetTopADCtime(adctime_top);
          } else {
            ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetTopADCtime(-999.);
          }
          // ((THcLADHodoHit *)fHodoHits->At(fNScinHits))->SetCalcPosition(kBig); //
          fGoodTopTdcTimeCorr.at(padnum - 1)    = timec_top;
          fGoodBtmTdcTimeCorr.at(padnum - 1)    = timec_btm;
          fGoodTopTdcTimeTOFCorr.at(padnum - 1) = kBig;
          fGoodBtmTdcTimeTOFCorr.at(padnum - 1) = kBig;
          */
        }
        fNScinHits++;
      }
    }
    ihit++;
  } // while loop

  if (problem_flag) {
    cout << "THcLADHodoPlane::ProcessHits " << fPlaneNum << " " << nexthit << "/" << nrawhits << endl;
    cout << " Ref problem end *******" << endl;
  }
  return (ihit);
}

//_______________________________________________________________________________________
Int_t THcLADHodoPlane::CoarseProcess(TClonesArray &tracks) {

  // Probably we won't do anything here

  return 0;
}

//_______________________________________________________________________________________
Int_t THcLADHodoPlane::FineProcess(TClonesArray &tracks) {

  // Probably we won't do anything here

  return 0;
}

//_____________________________________________________________________________
Int_t THcLADHodoPlane::AccumulatePedestals(TClonesArray *rawhits, Int_t nexthit) {
  /*! \brief Extract the data for this plane from raw hit list THcRawHodoHit, accumulating into arrays for calculating
   * pedestals.
   *
   * - Loop through raw data for scintillator plane
   */
  Int_t nrawhits = rawhits->GetLast() + 1;
  // cout << "THcScintillatorPlane::AcculatePedestals " << fPlaneNum << " " << nexthit << "/" << nrawhits << endl;

  Int_t ihit = nexthit;
  while (ihit < nrawhits) {
    THcRawHodoHit *hit = (THcRawHodoHit *)rawhits->At(ihit);
    if (hit->fPlane > fPlaneNum) {
      break;
    }
    Int_t element = hit->fCounter - 1;                       // Should check if in range
    Int_t adctop  = hit->GetRawAdcHitPos().GetPulseIntRaw(); // Pos=Top
    Int_t adcbtm  = hit->GetRawAdcHitNeg().GetPulseIntRaw(); // Neg=Btm

    if (adctop <= fTopPedLimit[element]) {
      fTopPedSum[element] += adctop;
      fTopPedSum2[element] += adctop * adctop;
      fTopPedCount[element]++;
      if (fTopPedCount[element] == fMinPeds / 5) {
        fTopPedLimit[element] = 100 + fTopPedSum[element] / fTopPedCount[element];
      }
    }
    if (adcbtm <= fBtmPedLimit[element]) {
      fBtmPedSum[element] += adcbtm;
      fBtmPedSum2[element] += adcbtm * adcbtm;
      fBtmPedCount[element]++;
      if (fBtmPedCount[element] == fMinPeds / 5) {
        fBtmPedLimit[element] = 100 + fBtmPedSum[element] / fBtmPedCount[element];
      }
    }
    ihit++;
  }

  fNPedestalEvents++;

  return (ihit);
}

//_______________________________________________________________________________________
void THcLADHodoPlane::CalculatePedestals() {
  /*! \brief   Calculate pedestals from arrays made in THcScintillatorPlane::AccumulatePedestals
   *
   * - Calculate pedestals from arrays made in THcScintillatorPlane::AccumulatePedestals
   * - In old fortran ENGINE code, a comparison was made between calculated pedestals and the pedestals read in by the
   * FASTBUS modules for zero supression. This is not implemented.
   */
  for (Int_t i = 0; i < fNelem; i++) {

    // Top tubes
    fTopPed[i]    = ((Double_t)fTopPedSum[i]) / TMath::Max(1, fTopPedCount[i]);
    fTopThresh[i] = fTopPed[i] + 15;

    // Bottom tubes
    fBtmPed[i]    = ((Double_t)fBtmPedSum[i]) / TMath::Max(1, fBtmPedCount[i]);
    fBtmThresh[i] = fBtmPed[i] + 15;

    //    cout <<"Pedestals "<< i+1 << " " << fTopPed[i] << " " << fBtmPed[i] << endl;
  }
  //  cout << " " << endl;
}

//_____________________________________________________________________________
void THcLADHodoPlane::InitializePedestals() {
  /*! \brief   called by THcHodoPlane::ReadDatabase
   *
   * - Initialize variables used in  THcScintillatorPlane::AccumulatePedestals and
   * THcScintillatorPlane::CalculatePedestals
   * - Minimum number of pedestal events needed for calculation, fMinPeds, hadrcoded to 500
   */
  fNPedestalEvents = 0;
  fMinPeds         = 500; // In engine, this is set in parameter file
  fTopPedSum       = new Int_t[fNelem];
  fTopPedSum2      = new Int_t[fNelem];
  fTopPedLimit     = new Int_t[fNelem];
  fTopPedCount     = new Int_t[fNelem];
  fBtmPedSum       = new Int_t[fNelem];
  fBtmPedSum2      = new Int_t[fNelem];
  fBtmPedLimit     = new Int_t[fNelem];
  fBtmPedCount     = new Int_t[fNelem];

  fTopPed    = new Double_t[fNelem];
  fBtmPed    = new Double_t[fNelem];
  fTopThresh = new Double_t[fNelem];
  fBtmThresh = new Double_t[fNelem];
  for (Int_t i = 0; i < fNelem; i++) {
    fTopPedSum[i]   = 0;
    fTopPedSum2[i]  = 0;
    fTopPedLimit[i] = 1000; // In engine, this are set in parameter file
    fTopPedCount[i] = 0;
    fBtmPedSum[i]   = 0;
    fBtmPedSum2[i]  = 0;
    fBtmPedLimit[i] = 1000; // In engine, this are set in parameter file
    fBtmPedCount[i] = 0;
  }
}
//____________________________________________________________________________

ClassImp(THcLADHodoPlane)
