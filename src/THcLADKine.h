#ifndef ROOT_THcLADKine
#define ROOT_THcLADKine
//////////////////////////////////////////////////////////////////////////

#include "THaSpectrometer.h"
#include "THaTrack.h"
#include "THcGoodLADHit.h"
#include "THcLADGEM.h"
#include "THcLADHodoscope.h"
#include "THcPrimaryKine.h"
#include "THcReactionPoint.h"
#include "THcTrigDet.h"

class THcReactionPoint;
class THcLADHodoscope;
class THcLADGEM;

class THcLADKine : public THcPrimaryKine {
public:
  THcLADKine(const char *name, const char *description = "", const char *spectro = "", const char *primary_kine = "",
             const char *vertex_module = "");

  virtual ~THcLADKine();
  virtual EStatus Init(const TDatime &run_time);
  virtual Int_t ReadDatabase(const TDatime &date);

  // not sure if we ever want to override this..
  virtual Int_t Process(const THaEvData &);
  virtual void Clear(Option_t *opt = "");
  void SetApparatus(const char *name);

protected:
  TString fSpecName;
  THcLADGEM *fGEM;
  THcLADHodoscope *fHodoscope;
  TString fVertexModuleName;
  THcReactionPoint *fVertexModule;
  THaTrack *fTrack; // Golden Track
  THcTrigDet *fTrigDet;

  Int_t MAXGOODHITS = 500;
  Int_t goodhit_n;
  TClonesArray *fGoodLADHits;
  TClonesArray *fGEMTracks;
  Double_t fD0Cut_wVertex;
  Double_t fD0Cut_noVertex;
  Double_t fMax_dTrk_horiz_match;
  Double_t fMax_dTrk_vert_match;
  Double_t fTrk_dtCut;
  Int_t fNfixed_z;
  Double_t *fFixed_z;
  Double_t fglobal_time_offset;
  Double_t fTVertex;
  Double_t fRFTime;
  Double_t fTVertex_RFcorr;
  Int_t n_rf_offsets;
  Double_t *rf_offset;
  Double_t rf_period;

  Double_t fZCellMin;  // default -15.0 cm
  Double_t fZCellMax;  // default +15.0 cm
  Double_t fThetaMin;  // default 60.0 deg
  Double_t fThetaMax;  // default 170.0 deg
  Double_t fPhiMin;    // default -50.0 deg
  Double_t fPhiMax;    // default +50.0 deg

  Double_t fchisq_cut[2];//chisq difference between 1 and 2 hodo hit track fits, used to determine if we accept tracks with only 1 hodo hit (if chisq_2hit - chisq_1hit > fchisq_cut[0]), or if we have no hodo hits (if chisq_1hit - chisq_0hit < fchisq_cut[1])

  Double_t fSigma_GEM; // GEM resolution in cm, used for track fitting, should be set based on detector performance
  Double_t fSigma_Hodo; // Hodoscope resolution in cm, used for track fitting, should be set based on detector performance




  virtual Int_t DefineVariables(EMode mode = kDefine);
  void CalculateTVertex();
  Double_t CalculateToF(Double_t t_raw);
  Double_t CalculateTOFRFcorr(Double_t t_raw);
  Double_t FitTrack(TVector3 vertex, std::vector<TVector3> sp_positions, std::vector<double> sp_resolutions, double dir[3]);

  ClassDef(THcLADKine, 0)
};

#endif /* ROOT_THcLADKine */