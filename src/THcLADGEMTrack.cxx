#include "THcLADGEMTrack.h"

//_________________________________________________________
THcLADGEMTrack::THcLADGEMTrack(Int_t nlayers) {
  fProjVz    = -999.;
  fProjVy    = -999.;
  fProjVx    = -999.;
  fNSp       = 0;
  fD0        = -999;
  fhasGoodD0 = kFALSE;
  fSp        = new GEM2DHits[nlayers];
  chisq      = -999.;
  ftheta     = -999.;
  fphi       = -999.;

  // no-vertex tracking quantities
  chisq_noTrackVertex        = -999.;
  ftheta_noTrackVertex       = -999.;
  fphi_noTrackVertex         = -999.;
  fProjVx_noTrackVertex      = -999.;
  fProjVy_noTrackVertex      = -999.;
  fProjVz_noTrackVertex      = -999.;
  fD0_noTrackVertex          = -999.;
  fHasHodoHit_noTrackVertex  = 0;
  fIsGoodTrack_noTrackVertex = false;
  fBestHodoHit_noTrackVertex = nullptr;

  // x-z tracking quantities (chi-square ignores the y position)
  chisq_xz        = -999.;
  ftheta_xz       = -999.;
  fphi_xz         = -999.;
  fProjVx_xz      = -999.;
  fProjVy_xz      = -999.;
  fProjVz_xz      = -999.;
  fD0_xz          = -999.;
  fHasHodoHit_xz  = 0;
  fIsGoodTrack_xz = false;
  fBestHodoHit_xz = nullptr;

  chisq_noTrackVertex_xz        = -999.;
  ftheta_noTrackVertex_xz       = -999.;
  fphi_noTrackVertex_xz         = -999.;
  fProjVx_noTrackVertex_xz      = -999.;
  fProjVy_noTrackVertex_xz      = -999.;
  fProjVz_noTrackVertex_xz      = -999.;
  fD0_noTrackVertex_xz          = -999.;
  fHasHodoHit_noTrackVertex_xz  = 0;
  fIsGoodTrack_noTrackVertex_xz = false;
  fBestHodoHit_noTrackVertex_xz = nullptr;
}

//_________________________________________________________
THcLADGEMTrack::~THcLADGEMTrack() {
  delete[] fSp;
  fSp = nullptr;
}

//_________________________________________________________
void THcLADGEMTrack::AddSpacePoint(GEM2DHits sp) {
  fSp[fNSp] = sp;
  fNSp++;
}

//_________________________________________________________

void THcLADGEMTrack::Clear(Option_t *opt) { fNSp = 0; }

ClassImp(THcLADGEMTrack)
