#ifndef ROOT_THcLADKine
#define ROOT_THcLADKine
//////////////////////////////////////////////////////////////////////////

#include "THaSpectrometer.h"
#include "THaTrack.h"
#include "THcLADGEM.h"
#include "THcLADHodoscope.h"
#include "THcPrimaryKine.h"
#include "THcReactionPoint.h"
#include "THcTrigDet.h"
#include "THcGoodLADHit.h"

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
  Double_t rf_offset[2];
  virtual Int_t DefineVariables(EMode mode = kDefine);
  void CalculateTVertex();
  Double_t CalculateToF(Double_t t_raw);
  Double_t CalculateTOFRFcorr(Double_t t_raw);

  ClassDef(THcLADKine, 0)
};

#endif /* ROOT_THcLADKine */