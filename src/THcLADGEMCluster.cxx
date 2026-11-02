#include "THcLADGEMCluster.h"
#include <iostream>

ClassImp(THcLADGEMCluster)

//__________________________________________________________________________
THcLADGEMCluster::THcLADGEMCluster()
{
  fClusteringFlag = -1;
  fNStrips = -1;
  fStripLow = -1;
  fStripHigh = -1;
  fStripMax = -1;

  fLayer = -1;
  fMPD = -1;
  fAPV = -1;
  fAxis = -1;
  fCLIndex = -1;
  fSplitSide = 0;
  fStripMaxReal = -1;

  fPos = -999.;
  fPosDiff = -999.;
  fPosMax = -999.;
  fPosSigma = -999.;
  fMoments = -999;
  fSampMax = -1;
  fADCsum = -999;
  fNHits = -1;
  fE = -999.;
  fTime = -999.;
  fTSigma = -999;
  fTimeFit = -999;
}

//__________________________________________________________________________
void THcLADGEMCluster::Clear(Option_t* opt)
{
  
  fClusteringFlag = -1;
  fNStrips = -1;
  fStripLow = -1;
  fStripHigh = -1;
  fStripMax = -1;

  fLayer = -1;
  fMPD = -1;
  fAPV = -1;
  fAxis = -1;
  fSplitSide = 0;
  fStripMaxReal = -1;

  fPos = -999.;
  fPosDiff = -999.;
  fPosMax = -999.;
  fPosSigma = -999.;
  fMoments = -999;
  fSampMax = -1;
  fADCsum = -999;
  fNHits = -1;
  fE = -999.;
  fTime = -999.;
  fTSigma = -999;
  fTimeFit = -999;

}

//__________________________________________________________________________
void THcLADGEMCluster::Print( Option_t* opt ) const
{

  std::cout << "Cluster Info: Nstrip= " << fNStrips
	    << " MaxStrip= " << fStripMax 
	    << " ADCsum= " << fADCsum
	    << " Time= " << fTime
	    << std::endl;
}
