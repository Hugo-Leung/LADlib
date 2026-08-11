#ifndef THcGoodLADHit_H
#define THcGoodLADHit_H
#include "TClonesArray.h"
#include "TObject.h"

class THcGoodLADHit : public TObject {
public:
  THcGoodLADHit() {
    for (int i = 0; i < 2; ++i) {
      plane[i] = paddle[i] = -1;
      hit_tof[i] = hit_tof_rfcorr[i] = hit_time[i] = hit_theta[i] = hit_phi[i] = hit_edep[i] = hit_edep_mev[i] =
          hit_edep_amp[i] = hit_edep_amp_mev[i] = hit_yPos[i] = hit_alpha[i] = hit_beta[i] = 1e30;
      is_proton[i]                                                                         = false;
    }
    trk_chiSqr = 1e30;
    track_id   = -1;
    // no-vertex tracking associations
    trk_chiSqr_noTrackVertex = 1e30;
    track_id_noTrackVertex   = -1;
    // x-z tracking associations (chi-square ignores the y position)
    trk_chiSqr_xz               = 1e30;
    track_id_xz                 = -1;
    trk_chiSqr_noTrackVertex_xz = 1e30;
    track_id_noTrackVertex_xz   = -1;
  };
  virtual ~THcGoodLADHit() = default;

  void CheckHitIndex(Int_t hit) const {
    if (hit < 0 || hit > 1) {
      throw std::out_of_range("Invalid hit index. Must be 0 or 1.");
    }
  }

  void SetPlane(Int_t hit, Int_t value) {
    CheckHitIndex(hit);
    plane[hit] = value;
  }
  void SetPaddle(Int_t hit, Int_t value) {
    CheckHitIndex(hit);
    paddle[hit] = value;
  }
  void SetTrackID(Int_t value) { track_id = value; }
  void SetTrkChiSqr(Double_t value) { trk_chiSqr = value; }
  void SetTrackID_noTrackVertex(Int_t value) { track_id_noTrackVertex = value; }
  void SetTrkChiSqr_noTrackVertex(Double_t value) { trk_chiSqr_noTrackVertex = value; }
  void SetTrackID_xz(Int_t value) { track_id_xz = value; }
  void SetTrkChiSqr_xz(Double_t value) { trk_chiSqr_xz = value; }
  void SetTrackID_noTrackVertex_xz(Int_t value) { track_id_noTrackVertex_xz = value; }
  void SetTrkChiSqr_noTrackVertex_xz(Double_t value) { trk_chiSqr_noTrackVertex_xz = value; }
  void SetIsProton(Int_t hit, Bool_t value) {
    CheckHitIndex(hit);
    is_proton[hit] = value;
  }
  void SetHitTime(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_time[hit] = value;
  }
  void SetHitTheta(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_theta[hit] = value;
  }
  void SetHitPhi(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_phi[hit] = value;
  }
  void SetHitEdep(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_edep[hit] = value;
  }
  void SetHitEdepMeV(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_edep_mev[hit] = value;
  }
  void SetHitEdepAmp(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_edep_amp[hit] = value;
  }
  void SetHitEdepAmpMeV(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_edep_amp_mev[hit] = value;
  }
  void SetHitAlpha(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_alpha[hit] = value;
  }
  void SetHitYPos(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_yPos[hit] = value;
  }
  void SetHitTOF(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_tof[hit] = value;
  }
  void SetHitTOFRFcorr(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_tof_rfcorr[hit] = value;
  }

  void SetHitBeta(Int_t hit, Double_t value) {
    CheckHitIndex(hit);
    hit_beta[hit] = value;
  }

  void CopyHit(Int_t this_plane, Int_t copy_plane, THcGoodLADHit *copyhit) {
    if (copy_plane == 0) {
      SetPlane(this_plane, copyhit->GetPlaneHit0());
      SetPaddle(this_plane, copyhit->GetPaddleHit0());
      SetTrackID(copyhit->GetTrackID());
      SetTrkChiSqr(copyhit->GetTrkChiSqr());
      SetIsProton(this_plane, copyhit->GetIsProtonHit0());
      SetHitTime(this_plane, copyhit->GetHitTimeHit0());
      SetHitTheta(this_plane, copyhit->GetHitThetaHit0());
      SetHitPhi(this_plane, copyhit->GetHitPhiHit0());
      SetHitEdep(this_plane, copyhit->GetHitEdepHit0());
      SetHitEdepMeV(this_plane, copyhit->GetHitEdepMeVHit0());
      SetHitEdepAmp(this_plane, copyhit->GetHitEdepAmpHit0());
      SetHitEdepAmpMeV(this_plane, copyhit->GetHitEdepAmpMeVHit0());
      SetHitAlpha(this_plane, copyhit->GetHitAlphaHit0());
      SetHitYPos(this_plane, copyhit->GetHitYPosHit0());
      SetHitTOF(this_plane, copyhit->GetHitTOFHit0());
      SetHitTOFRFcorr(this_plane, copyhit->GetHitTOFRFcorrHit0());
    } else if (copy_plane == 1) {
      SetPlane(this_plane, copyhit->GetPlaneHit1());
      SetPaddle(this_plane, copyhit->GetPaddleHit1());
      SetTrackID(copyhit->GetTrackID());
      SetTrkChiSqr(copyhit->GetTrkChiSqr());
      SetIsProton(this_plane, copyhit->GetIsProtonHit1());
      SetHitTime(this_plane, copyhit->GetHitTimeHit1());
      SetHitTheta(this_plane, copyhit->GetHitThetaHit1());
      SetHitPhi(this_plane, copyhit->GetHitPhiHit1());
      SetHitEdep(this_plane, copyhit->GetHitEdepHit1());
      SetHitEdepMeV(this_plane, copyhit->GetHitEdepMeVHit1());
      SetHitEdepAmp(this_plane, copyhit->GetHitEdepAmpHit1());
      SetHitEdepAmpMeV(this_plane, copyhit->GetHitEdepAmpMeVHit1());
      SetHitAlpha(this_plane, copyhit->GetHitAlphaHit1());
      SetHitYPos(this_plane, copyhit->GetHitYPosHit1());
      SetHitTOF(this_plane, copyhit->GetHitTOFHit1());
      SetHitTOFRFcorr(this_plane, copyhit->GetHitTOFRFcorrHit1());
    } else {
      throw std::out_of_range("Invalid copy_plane index. Must be 0 or 1.");
    }
  }

  // Declaring two different methods for each hit is dumb, but RDefVars can't take TClonesArray
  //  objects calling a method with a parameter (or I haven't been able to figure it out)
  //  so we have to do this
  Int_t GetPlaneHit0() const { return plane[0]; }
  Int_t GetPlaneHit1() const { return plane[1]; }

  Int_t GetPaddleHit0() const { return paddle[0]; }
  Int_t GetPaddleHit1() const { return paddle[1]; }

  Int_t GetTrackID() const { return track_id; }
  Int_t GetTrackID_noTrackVertex() const { return track_id_noTrackVertex; }
  Int_t GetTrackID_xz() const { return track_id_xz; }
  Int_t GetTrackID_noTrackVertex_xz() const { return track_id_noTrackVertex_xz; }

  Double_t GetBetaHit0() const { return hit_beta[0]; }
  Double_t GetBetaHit1() const { return hit_beta[1]; }

  Double_t GetTrkChiSqr() const { return trk_chiSqr; }
  Double_t GetTrkChiSqr_noTrackVertex() const { return trk_chiSqr_noTrackVertex; }
  Double_t GetTrkChiSqr_xz() const { return trk_chiSqr_xz; }
  Double_t GetTrkChiSqr_noTrackVertex_xz() const { return trk_chiSqr_noTrackVertex_xz; }

  Double_t GetIsProtonHit0() const { return is_proton[0]; }
  Double_t GetIsProtonHit1() const { return is_proton[1]; }

  Double_t GetHitTimeHit0() const { return hit_time[0]; }
  Double_t GetHitTimeHit1() const { return hit_time[1]; }

  Double_t GetHitThetaHit0() const { return hit_theta[0]; }
  Double_t GetHitThetaHit1() const { return hit_theta[1]; }

  Double_t GetHitPhiHit0() const { return hit_phi[0]; }
  Double_t GetHitPhiHit1() const { return hit_phi[1]; }

  Double_t GetHitEdepHit0() const { return hit_edep[0]; }
  Double_t GetHitEdepHit1() const { return hit_edep[1]; }

  Double_t GetHitEdepMeVHit0() const { return hit_edep_mev[0]; }
  Double_t GetHitEdepMeVHit1() const { return hit_edep_mev[1]; }

  Double_t GetHitEdepAmpHit0() const { return hit_edep_amp[0]; }
  Double_t GetHitEdepAmpHit1() const { return hit_edep_amp[1]; }

  Double_t GetHitEdepAmpMeVHit0() const { return hit_edep_amp_mev[0]; }
  Double_t GetHitEdepAmpMeVHit1() const { return hit_edep_amp_mev[1]; }

  Double_t GetHitAlphaHit0() const { return hit_alpha[0]; }
  Double_t GetHitAlphaHit1() const { return hit_alpha[1]; }

  Double_t GetHitYPosHit0() const { return hit_yPos[0]; }
  Double_t GetHitYPosHit1() const { return hit_yPos[1]; }

  Double_t GetHitTOFHit0() const { return hit_tof[0]; }
  Double_t GetHitTOFHit1() const { return hit_tof[1]; }

  Double_t GetHitTOFRFcorrHit0() const { return hit_tof_rfcorr[0]; }
  Double_t GetHitTOFRFcorrHit1() const { return hit_tof_rfcorr[1]; }

protected:
  Int_t plane[2];
  Int_t paddle[2];
  Int_t track_id;
  Double_t trk_chiSqr;
  Int_t track_id_noTrackVertex;         // GEM track ID from the no-vertex fit
  Double_t trk_chiSqr_noTrackVertex;    // chi-square of the associated no-vertex track
  Int_t track_id_xz;                    // GEM track ID from the vertex-constrained x-z fit
  Double_t trk_chiSqr_xz;               // chi-square of the associated x-z (no-y) track
  Int_t track_id_noTrackVertex_xz;      // GEM track ID from the no-vertex x-z fit
  Double_t trk_chiSqr_noTrackVertex_xz; // chi-square of the associated no-vertex x-z track
  Double_t is_proton[2];
  Double_t hit_time[2];
  Double_t hit_beta[2];
  Double_t hit_theta[2];
  Double_t hit_phi[2];
  Double_t hit_edep[2];
  Double_t hit_edep_mev[2];
  Double_t hit_edep_amp[2];
  Double_t hit_edep_amp_mev[2];
  Double_t hit_tof[2];
  Double_t hit_tof_rfcorr[2];
  Double_t hit_yPos[2];
  Double_t hit_alpha[2];

private:
  THcGoodLADHit(const THcGoodLADHit &);            // Prevent copy constructor
  THcGoodLADHit &operator=(const THcGoodLADHit &); // Prevent assignment operator
  ClassDef(THcGoodLADHit, 0)                       // Class for good hodo hits
};

#endif
//_____________________________________________________________________________