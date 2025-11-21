#ifndef THcLADHodoPlane_h
#define THcLADHodoPlane_h

#include "TClonesArray.h"
#include "THaSubDetector.h"
#include <vector>

using namespace std;

class THaEvData;
class THaSignalHit;

// rough class outline for LAD Hodoscope plane
// contains xx number of paddles per plane
// For each hodoscope detector, we will define it as a separate detector
// since they are not really placed together

class THcLADHodoPlane : public THaSubDetector {

public:
  THcLADHodoPlane(const char *name, const char *description, const Int_t planenum, THaDetectorBase *parent = nullptr);
  virtual ~THcLADHodoPlane();

  virtual void Clear(Option_t *opt = "");
  virtual Int_t Decode(const THaEvData &);
  virtual EStatus Init(const TDatime &date);

  virtual Int_t CoarseProcess(TClonesArray &tracks);
  virtual Int_t FineProcess(TClonesArray &tracks);
  virtual Int_t ProcessHits(TClonesArray *rawhits, Int_t nexthit);

  // these are called in THcLADHodoscope::Decode, to analyze pedestal events
  virtual Int_t AccumulatePedestals(TClonesArray *rawhits, Int_t nexthit);
  virtual void CalculatePedestals();

  // Some getter/setter functions here

  Int_t GetNelem() { return fNelem; };          // return number of paddles in this plane
  Int_t GetNScinHits() { return fNScinHits; };  // Get # hits in plane (that pass min/max TDC cuts)
  Int_t GetNGoodHits() { return fNGoodHits; };  // Get # hits in plane (used in determining focal plane time)
  Double_t GetSpacing() { return fSpacing; };   // spacing of paddles
  Double_t GetSize() { return fSize; };         // paddle size
  Double_t GetHodoSlop() { return fHodoSlop; }; // hodo slop
  Double_t GetZpos() { return fZpos; };         // return the z position
  Double_t GetDzpos() { return fDzpos; }
  Double_t GetTheta() { return fTheta; };
  Double_t GetPosBtm(Int_t PaddleNo) { return fPosBtm[PaddleNo]; };
  Double_t GetPosTop(Int_t PaddleNo) { return fPosTop[PaddleNo]; };
  Double_t GetPosOffset() { return fPosOffset; };
  Double_t GetPosCenter(Int_t PaddleNo) { return fPosCenter[PaddleNo]; }; // counting from zero!
  Double_t GetEdep2MeV_int(Int_t PaddleNo) {
    return fEdep2MeV_int[PaddleNo];
  } // conversion factor for energy deposition to MeV
  Double_t GetEdep2MeV_amp(Int_t PaddleNo) {
    return fEdep2MeV_amp[PaddleNo];
  } // conversion factor for energy deposition to MeV

  TClonesArray *GetHits() { return fHodoHits; }

protected:
  const Int_t NUM_PADDLES = 11; // number of paddles in a plane
  const Int_t MAX_HITS    = 40; // maximum number of hits in a plane (TClones array will automatically resize if needed)

  TClonesArray *fHodoHits;

  TClonesArray *frTopTdcTimeRaw;
  TClonesArray *frTopTdcTime;
  TClonesArray *frBtmTdcTimeRaw;
  TClonesArray *frBtmTdcTime;

  Double_t fTopTdcRefTime;
  Double_t fBtmTdcRefTime;
  Double_t fTopTdcRefDiffTime;
  Double_t fBtmTdcRefDiffTime;

  Double_t fTopAdcRefTime;
  Double_t fBtmAdcRefTime;
  Double_t fTopAdcRefDiffTime;
  Double_t fBtmAdcRefDiffTime;

  Double_t fTDC_match_window;
  // Counters
  Int_t fTotNumTdcHits;
  Int_t fTotNumTopTdcHits;
  Int_t fTotNumBtmTdcHits;

  Int_t fTotNumAdcHits;
  Int_t fTotNumTopAdcHits;
  Int_t fTotNumBtmAdcHits;

  vector<Int_t> fNumTopAdcHits;
  vector<Int_t> fNumBtmAdcHits;
  vector<Int_t> fNumTopTdcHits;
  vector<Int_t> fNumBtmTdcHits;

  Int_t fTotNumGoodTopAdcHits;
  Int_t fTotNumGoodBtmAdcHits;
  Int_t fTotNumGoodAdcHits;

  Int_t fTotNumGoodTopTdcHits;
  Int_t fTotNumGoodBtmTdcHits;
  Int_t fTotNumGoodTdcHits;

  vector<Int_t> fNumGoodTopAdcHits;
  vector<Int_t> fNumGoodBtmAdcHits;
  vector<Int_t> fNumGoodTopTdcHits;
  vector<Int_t> fNumGoodBtmTdcHits;

  vector<Double_t> fGoodTopAdcPed;
  vector<Double_t> fGoodTopAdcMult;
  vector<Double_t> fGoodTopAdcHitUsed;
  vector<Double_t> fGoodTopAdcPulseInt;
  vector<Double_t> fGoodTopAdcPulseAmp;
  vector<Double_t> fGoodTopAdcPulseTime;
  vector<Double_t> fGoodTopAdcTdcDiffTime;

  vector<Double_t> fGoodBtmAdcPed;
  vector<Double_t> fGoodBtmAdcMult;
  vector<Double_t> fGoodBtmAdcHitUsed;
  vector<Double_t> fGoodBtmAdcPulseInt;
  vector<Double_t> fGoodBtmAdcPulseAmp;
  vector<Double_t> fGoodBtmAdcPulseTime;
  vector<Double_t> fGoodBtmAdcTdcDiffTime;

  TClonesArray *frTopAdcPedRaw;
  TClonesArray *frTopAdcPed;
  TClonesArray *frTopAdcPulseIntRaw;
  TClonesArray *frTopAdcPulseInt;
  TClonesArray *frTopAdcPulseAmpRaw;
  TClonesArray *frTopAdcPulseAmp;
  TClonesArray *frTopAdcPulseTimeRaw;
  TClonesArray *frTopAdcPulseTime;
  TClonesArray *frTopAdcErrorFlag;

  TClonesArray *frTopAdcSampPedRaw;
  TClonesArray *frTopAdcSampPed;
  TClonesArray *frTopAdcSampPulseIntRaw;
  TClonesArray *frTopAdcSampPulseInt;
  TClonesArray *frTopAdcSampPulseAmpRaw;
  TClonesArray *frTopAdcSampPulseAmp;
  TClonesArray *frTopAdcSampPulseTimeRaw;
  TClonesArray *frTopAdcSampPulseTime;

  TClonesArray *frBtmAdcPedRaw;
  TClonesArray *frBtmAdcPed;
  TClonesArray *frBtmAdcPulseIntRaw;
  TClonesArray *frBtmAdcPulseInt;
  TClonesArray *frBtmAdcPulseAmpRaw;
  TClonesArray *frBtmAdcPulseAmp;
  TClonesArray *frBtmAdcPulseTimeRaw;
  TClonesArray *frBtmAdcPulseTime;
  TClonesArray *frBtmAdcErrorFlag;

  TClonesArray *frBtmAdcSampPedRaw;
  TClonesArray *frBtmAdcSampPed;
  TClonesArray *frBtmAdcSampPulseIntRaw;
  TClonesArray *frBtmAdcSampPulseInt;
  TClonesArray *frBtmAdcSampPulseAmpRaw;
  TClonesArray *frBtmAdcSampPulseAmp;
  TClonesArray *frBtmAdcSampPulseTimeRaw;
  TClonesArray *frBtmAdcSampPulseTime;

  vector<Double_t> fTopAdcSampWaveform;
  vector<Double_t> fBtmAdcSampWaveform;

  // Hodoscopoe "GOOD" TDC Variables
  vector<Double_t> fGoodTopTdcTimeUnCorr;
  vector<Double_t> fGoodTopTdcTimeCorr;

  vector<Double_t> fGoodBtmTdcTimeUnCorr;
  vector<Double_t> fGoodBtmTdcTimeCorr;

  // Time Walk Corrected
  vector<Double_t> fGoodTopTdcTimeWalkCorr;
  vector<Double_t> fGoodBtmTdcTimeWalkCorr;
  vector<Double_t> fGoodDiffDistTrack;
  vector<Double_t> fGoodHitTimeDiff;
  vector<Double_t> fGoodHitTimeAvg;

  TClonesArray *frTopTdcHits;
  TClonesArray *frTopAdcHits;
  TClonesArray *frTopAdcSums;
  TClonesArray *frTopAdcPeds;
  TClonesArray *frBtmTdcHits;
  TClonesArray *frBtmAdcHits;
  TClonesArray *frBtmAdcSums;
  TClonesArray *frBtmAdcPeds;

  Int_t fPlaneNum;
  Int_t fNelem;

  Int_t fNScinHits;     /* number of hits in plane (that pass min/max TDC cuts) */
  Int_t fNGoodHits;     /* number of hits in plane (used in determining focal plane time) */
  Int_t fNScinGoodHits; /* number of hits for which both ends of the paddle fired in time -- diff from fNGoodHits? */
  Double_t fpTime;

  Double_t fHitDistance;
  Double_t fScinXPos;
  Double_t fScinYPos;
  Double_t fTrackXPosition;
  Double_t fTrackYPosition;
  Double_t *fPosCenter; /* array with centers for all scintillators in the plane */

  Int_t fDebugAdc;
  Int_t fADCMode;
  Int_t fADCDiagCut;
  Int_t fCosmicFlag;
  Double_t fSampThreshold;
  Int_t fSampNSA;
  Int_t fSampNSAT;
  Int_t fSampNSB;
  Int_t fOutputSampWaveform;
  Int_t fUseSampWaveform;
  Int_t fIsMC;

  Double_t *fPosBtm; /* array with bottom position for all scintillators in the plane */
  Double_t *fPosTop; /* array with top position for all scintillators in the plane */
  Double_t fPosOffset;
  Int_t fTdcOffset;
  Double_t fAdcTdcOffset;
  Double_t fScinTdcMin;
  Double_t fScinTdcMax;
  Double_t fScinTdcToTime;
  Double_t fBetaNominal;

  Double_t *fHodoTopAdcTimeWindowMin;
  Double_t *fHodoTopAdcTimeWindowMax;
  Double_t *fHodoBtmAdcTimeWindowMin;
  Double_t *fHodoBtmAdcTimeWindowMax;

  // Hodoscope Calib Parameters
  Double_t *fHodoVelLight;

  // Time-Walk Parameters
  Double_t *fHodoVelFit;
  Double_t *fHodoCableFit;
  Double_t *fHodo_LCoeff;
  Double_t *fHodoVelFit_FADC;
  Double_t *fHodoCableFit_FADC;
  Double_t *fHodo_LCoeff_FADC;
  Double_t *fHodoTop_c1;
  Double_t *fHodoBtm_c1;
  Double_t *fHodoTop_c2;
  Double_t *fHodoBtm_c2;
  Double_t fTdc_Thrs;

  Double_t tw_corr_top;
  Double_t tw_corr_btm;

  // Edep calibrations
  Double_t *fEdep2MeV_int;
  Double_t *fEdep2MeV_amp;

  // Pedestal calculations
  Int_t fNPedestalEvents; /* Number of pedestal events */
  Int_t fMinPeds;         /* Only analyze/update if num events > */
  Int_t *fTopPedSum;      /* Accumulators for pedestals */
  Int_t *fTopPedSum2;
  Int_t *fTopPedLimit;
  Int_t *fTopPedCount;
  Int_t *fBtmPedSum;
  Int_t *fBtmPedSum2;
  Int_t *fBtmPedLimit;
  Int_t *fBtmPedCount;

  Double_t *fTopPed;
  Double_t *fTopSig;
  Double_t *fTopThresh;
  Double_t *fBtmPed;
  Double_t *fBtmSig;
  Double_t *fBtmThresh;

  Double_t afterPulse_window = 50; // ns (should make this a parameter, but hard-coded for now)

  enum { kADCStandard = 0, kADCDynamicPedestal, kADCSampleIntegral, kADCSampIntDynPed };

  // Geometry parameters
  Double_t fSpacing; /* paddle spacing */
  Double_t fSize;    /* paddle size */
  Double_t fZpos;    /* z position */
  Double_t fDzpos;
  Double_t fTheta;
  Double_t fHodoSlop;

  virtual Int_t ReadDatabase(const TDatime &date);
  virtual Int_t DefineVariables(EMode mode = kDefine);
  virtual void InitializePedestals();

  class HodoCluster {
  public:
    Int_t fClusterNumber; // cluster index
    Double_t fClusterPos;
    Double_t fClusterSize;
  };

  vector<HodoCluster> fCluster;

public:
  vector<HodoCluster> GetClusters() { return fCluster; }

  ClassDef(THcLADHodoPlane, 0);
};

#endif /* THcLADHodPlane_h */
