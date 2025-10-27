// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoDreamProducerTaskForSpecificAnalysis.cxx
/// \brief Tasks that reads the track tables and creates track triplets; only three identical particles can be used
/// \author Laura Serksnyte, TU München, laura.serksnyte@tum.de

#include <vector>
#include <string>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoDream/Core/femtoDreamParticleHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamEventHisto.h"
#include "PWGCF/FemtoDream/Core/femtoDreamPairCleaner.h"
#include "PWGCF/FemtoDream/Core/femtoDreamContainerThreeBody.h"
#include "PWGCF/FemtoDream/Core/femtoDreamDetaDphiStar.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod;

struct femtoDreamProducerTaskForSpecificAnalysis {

  SliceCache cache;

  Produces<aod::StoredFDCollisions> outputCollision;
  Produces<aod::StoredFDParticles> outputParts;

  Preslice<aod::FDParticles> perCol = aod::femtodreamparticle::fdCollisionId;
  float mMassOne = -999, mMassTwo = -999, mMassThree = -999;
  int collisions = 0;

  Configurable<int> ConfNumberOfTracks{"ConfNumberOfTracks", 3, "Number of tracks"};
  Configurable<int> ConfNumberOfV0{"ConfNumberOfV0", 0, "Number of V0"};
  Configurable<int> ConfNumberOfCascades{"ConfNumberOfCascades", 0, "Number of Cascades"};

  /// Track selection
  //Configurable<float> ConfPIDthrMom{"ConfPIDthrMom", 1.f, "Momentum threshold from which TPC and TOF are required for PID"};
  //Configurable<o2::aod::femtodreamparticle::cutContainerType> ConfTPCPIDBit{"ConfTPCPIDBit", 16, "PID TPC bit from cutCulator "};
  //Configurable<o2::aod::femtodreamparticle::cutContainerType> ConfTPCTOFPIDBit{"ConfTPCTOFPIDBit", 8, "PID TPCTOF bit from cutCulator"};

  /// Partition for selected particles
  //Partition<aod::FDParticles> SelectedParts = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
  //                                            ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= ConfPIDthrMom, ncheckbit(aod::femtodreamparticle::pidcut, ConfTPCPIDBit), ncheckbit(aod::femtodreamparticle::pidcut, ConfTPCTOFPIDBit));

  /// Particle 1 (track)
  struct : ConfigurableGroup {
    std::string prefix = std::string("Track1");
    Configurable<int> PDGCode{"PDGCode", 2212, "PDG code of Particle 1 (Track)"};
    Configurable<femtodreamparticle::cutContainerType> CutBit{"CutBit", 5542474, "Particle 1 (Track) - Selection bit from cutCulator"};
    Configurable<femtodreamparticle::cutContainerType> CutBitAnti{"CutBitAnti", 5542473, "Particle 1 (Track) - Selection bit from cutCulator"};
    Configurable<femtodreamparticle::cutContainerType> TPCBit{"TPCBit", 4, "PID TPC bit from cutCulator for particle 1 (Track)"};
    Configurable<femtodreamparticle::cutContainerType> TPCBit_Reject{"TPCBit_Reject", 0, "Reject PID TPC bit from cutCulator for particle 1 (Track). Set to 0 to turn off"};
    Configurable<femtodreamparticle::cutContainerType> TPCTOFBit{"TPCTOFBit", 2, "PID TPCTOF bit from cutCulator for particle 1 (Track)"};
    Configurable<float> PIDThres{"PIDThres", 0.75, "Momentum threshold for PID selection for particle 1 (Track)"};
    Configurable<float> PtMin{"PtMin", 0., "Minimum pT of partricle 1 (Track)"};
    Configurable<float> PtMax{"PtMax", 999., "Maximum pT of partricle 1 (Track)"};
    Configurable<float> EtaMin{"EtaMin", -10., "Minimum eta of partricle 1 (Track)"};
    Configurable<float> EtaMax{"EtaMax", 10., "Maximum eta of partricle 1 (Track)"};
    Configurable<float> TempFitVarMin{"TempFitVarMin", -10., "Minimum DCAxy of partricle 1 (Track)"};
    Configurable<float> TempFitVarMax{"TempFitVarMax", 10., "Maximum DCAxy of partricle 1 (Track)"};
    Configurable<bool> DCACutPtDep{"DCACutPtDep", false, "Use pt dependent dca cut"};
  } Track1;
  Partition<FDParticles> SelectedParts = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kTrack)) &&
                                         ifnode(aod::femtodreamparticle::pt * (nexp(aod::femtodreamparticle::eta) + nexp(-1.f * aod::femtodreamparticle::eta)) / 2.f <= Track1.PIDThres, ncheckbit(aod::femtodreamparticle::pidcut, Track1.TPCBit) && ((aod::femtodreamparticle::pidcut & Track1.TPCBit_Reject) == 0u), ncheckbit(aod::femtodreamparticle::pidcut, Track1.TPCTOFBit)) &&
                                         (aod::femtodreamparticle::pt > Track1.PtMin) &&
                                         (aod::femtodreamparticle::pt < Track1.PtMax) &&
                                         (aod::femtodreamparticle::eta > Track1.EtaMin) &&
                                         (aod::femtodreamparticle::eta < Track1.EtaMax) &&
                                         ifnode(Track1.DCACutPtDep, (nabs(aod::femtodreamparticle::tempFitVar) <= 0.0105f + (0.035f / npow(aod::femtodreamparticle::pt, 1.1f))),
                                                ((aod::femtodreamparticle::tempFitVar >= Track1.TempFitVarMin) &&
                                                 (aod::femtodreamparticle::tempFitVar <= Track1.TempFitVarMax)));
  
  /// V0 selection
  Configurable<float> Conf_minInvMass_V0{"Conf_minInvMass_V0", 1.08, "Minimum invariant mass of V0 (particle)"};
  Configurable<float> Conf_maxInvMass_V0{"Conf_maxInvMass_V0", 1.15, "Maximum invariant mass of V0 (particle)"};
  Configurable<float> Conf_minInvMassAnti_V0{"Conf_minInvMassAnti_V0", 1.08, "Minimum invariant mass of V0 (antiparticle)"};
  Configurable<float> Conf_maxInvMassAnti_V0{"Conf_maxInvMassAnti_V0", 1.15, "Maximum invariant mass of V0 (antiparticle)"};
  /// Cascade selection
  //Configurable<float> Conf_minInvMass_Cascade{"Conf_minInvMass_Cascade", 1.2, "Minimum invariant mass of Cascade (particle)"};
  //Configurable<float> Conf_maxInvMass_Cascade{"Conf_maxInvMass_Cascade", 1.5, "Maximum invariant mass of Cascade (particle)"};
  
  // Partition for selected particles
  Partition<aod::FDParticles> SelectedV0s = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kV0));
  //Partition<aod::FDParticles> SelectedCascades = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kCascade));

  struct : ConfigurableGroup {
    std::string prefix = std::string("Cascade2");
    Configurable<int> PDGCode{"PDGCode", 3312, "PDG code of particle 2 (V0)"};
    Configurable<femtodreamparticle::cutContainerType> CutBit{"CutBit", 32221874, "Selection bit for particle 2 (Cascade)"};
    Configurable<femtodreamparticle::cutContainerType> ChildPos_CutBit{"ChildPos_CutBit", 278, "Selection bit for positive child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> ChildPos_TPCBit{"ChildPos_TPCBit", 1024, "PID TPC bit for positive child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> ChildNeg_CutBit{"ChildNeg_CutBit", 277, "Selection bit for negative child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> ChildNeg_TPCBit{"ChildNeg_TPCBit", 4096, "PID TPC bit for negative child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> ChildBach_CutBit{"ChildBach_CutBit", 277, "Selection bit for negative child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> ChildBach_TPCBit{"ChildBach_TPCBit", 64, "PID TPC bit for bachelor child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> CutBitAnti{"CutBitAnti", 32221874, "Selection bit for particle 2 (Cascade)"};
    Configurable<femtodreamparticle::cutContainerType> ChildPos_CutBitAnti{"ChildPos_CutBitAnti", 278, "Selection bit for positive child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> ChildPos_TPCBitAnti{"ChildPos_TPCBitAnti", 1024, "PID TPC bit for positive child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> ChildNeg_CutBitAnti{"ChildNeg_CutBitAnti", 277, "Selection bit for negative child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> ChildNeg_TPCBitAnti{"ChildNeg_TPCBitAnti", 4096, "PID TPC bit for negative child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> ChildBach_CutBitAnti{"ChildBach_CutBitAnti", 277, "Selection bit for negative child of Cascade"};
    Configurable<femtodreamparticle::cutContainerType> ChildBach_TPCBitAnti{"ChildBach_TPCBitAnti", 64, "PID TPC bit for bachelor child of Cascade"};
    Configurable<float> InvMassMin{"InvMassMin", 1.2, "Minimum invariant mass of Partricle 2 (particle) (Cascade)"};
    Configurable<float> InvMassMax{"InvMassMax", 1.4, "Maximum invariant mass of Partricle 2 (particle) (Cascade)"};
    Configurable<float> InvMassV0DaughMin{"InvMassV0DaughMin", 0., "Minimum invariant mass of the V0 Daughter"};
    Configurable<float> InvMassV0DaughMax{"InvMassV0DaughMax", 999., "Maximum invariant mass of the V0 Daughter"};
    Configurable<float> PtMin{"PtMin", 0., "Minimum pT of Partricle 2 (V0)"};
    Configurable<float> PtMax{"PtMax", 999., "Maximum pT of Partricle 2 (V0)"};
    Configurable<float> EtaMin{"EtaMin", -10., "Minimum eta of Partricle 2 (V0)"};
    Configurable<float> EtaMax{"EtaMax", 10., "Maximum eta of Partricle 2 (V0)"};
    Configurable<bool> UseChildCuts{"UseChildCuts", true, "Use cuts on the children of the Cascades additional to those of the selection of the cascade builder (for debugging purposes)"};
    Configurable<bool> UseChildPIDCuts{"UseChildPIDCuts", true, "Use PID cuts on the children of the Cascades additional to those of the selection of the cascade builder (for debugging purposes)"};
    Configurable<bool> ConfIsOmega{"ConfIsOmega", false, "Switch between Xi and Omaga Cascades: If true: Omega; else: Xi"};
    Configurable<bool> ConfRejectCompetingMass{"ConfRejectCompetingMass", false, "Reject the competing Cascade Mass (use only for debugging. More efficient to exclude it already at the producer level)"};
    Configurable<float> ConfCompetingCascadeMassLowLimit{"ConfCompetingCascadeMassLowLimit", 0., "Lower Limit of the invariant mass window within which to reject the cascade"};
    Configurable<float> ConfCompetingCascadeMassUpLimit{"ConfCompetingCascadeMassUpLimit", 0., "Upper Limit of the invariant mass window within which to reject the cascade"};
  } Cascade2;
  /// Partition for particle 2
  Partition<FDParticles> SelectedCascades = (aod::femtodreamparticle::partType == uint8_t(aod::femtodreamparticle::ParticleType::kCascade)) &&
                                             (aod::femtodreamparticle::pt > Cascade2.PtMin) &&
                                             (aod::femtodreamparticle::pt < Cascade2.PtMax) &&
                                             (aod::femtodreamparticle::eta > Cascade2.EtaMin) &&
                                             (aod::femtodreamparticle::eta < Cascade2.EtaMax) &&
                                             (aod::femtodreamparticle::mLambda > Cascade2.InvMassMin) &&
                                             (aod::femtodreamparticle::mLambda < Cascade2.InvMassMax) &&
                                             (aod::femtodreamparticle::mAntiLambda > Cascade2.InvMassV0DaughMin) &&
                                             (aod::femtodreamparticle::mAntiLambda < Cascade2.InvMassV0DaughMax);


  HistogramRegistry EventRegistry{"EventRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr uint32_t kSignPlusMask = 1 << 1;

  template <typename T>
  int getRowDaughters(int daughID, T const& vecID)
  {
    int rowInPrimaryTrackTableDaugh = -1;
    for (size_t i = 0; i < vecID.size(); i++) {
      if (vecID.at(i) == daughID) {
        rowInPrimaryTrackTableDaugh = i;
        break;
      }
    }
    return rowInPrimaryTrackTableDaugh;
  }

  void init(InitContext&)
  {
    EventRegistry.add("hStatistiscs", ";bin;Entries", kTH1F, {{3, 0, 3}});
    // get bit for the collision mask
  }
  /// This function stores accepted collisions in derived data
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param groupSelectedTracks partition for the first particle passed by the process function
  /// @param groupSelectedV0s partition for the second particle passed by the process function
  /// @param parts femtoDreamParticles table
  template <bool isMC, typename PartitionType, typename PartType>
  void createSpecifiedDerivedData(o2::aod::FDCollision& col, PartitionType groupSelectedTracks, PartitionType groupSelectedV0s, PartType parts)
  {
    /// check tracks
    int tracksCount = 0;
    int antitracksCount = 0;
    for (auto& part : groupSelectedTracks) {
      if (part.cut() & 1) {
        antitracksCount++;
      } else {
        tracksCount++;
      }
    }

    /// check V0s
    int V0Count = 0;
    int antiV0Count = 0;
    for (auto& V0 : groupSelectedV0s) {
      if ((V0.mLambda() > Conf_minInvMass_V0) && (V0.mLambda() < Conf_maxInvMass_V0)) {
        V0Count++;
      } else if ((V0.mAntiLambda() > Conf_minInvMassAnti_V0) && (V0.mAntiLambda() < Conf_maxInvMassAnti_V0)) {
        antiV0Count++;
      }
    }

    std::vector<int> tmpIDtrack;

    if ((V0Count >= ConfNumberOfV0 && tracksCount >= ConfNumberOfTracks) || (antiV0Count >= ConfNumberOfV0 && antitracksCount >= ConfNumberOfTracks)) {
      EventRegistry.fill(HIST("hStatistiscs"), 1);
      outputCollision(col.posZ(), col.multV0M(), col.multNtr(), col.sphericity(), col.magField());
      for (auto& femtoParticle : parts) {
        if (aod::femtodreamparticle::ParticleType::kTrack == femtoParticle.partType()) {
          std::vector<int> childIDs = {0, 0};
          outputParts(outputCollision.lastIndex(),
                      femtoParticle.pt(),
                      femtoParticle.eta(),
                      femtoParticle.phi(),
                      femtoParticle.partType(),
                      femtoParticle.cut(),
                      femtoParticle.pidcut(),
                      femtoParticle.tempFitVar(),
                      childIDs,
                      femtoParticle.mLambda(),
                      femtoParticle.mAntiLambda());
          tmpIDtrack.push_back(femtoParticle.index());
        }
        if (aod::femtodreamparticle::ParticleType::kV0Child == femtoParticle.partType()) {
          std::vector<int> childIDs = {0, 0};
          const auto& children = femtoParticle.childrenIds();
          int childId = (children[0] != 0) ? children[0] : children[1];
          if (childId != -1) {
            int rowInPrimaryTrackTable = getRowDaughters(childId, tmpIDtrack);
            childIDs = (children[0] != 0) ? std::vector<int>{rowInPrimaryTrackTable, 0} : std::vector<int>{0, rowInPrimaryTrackTable};
          } else {
            childIDs = (children[0] != 0) ? std::vector<int>{-1, 0} : std::vector<int>{0, -1};
          }
          outputParts(outputCollision.lastIndex(),
                      femtoParticle.pt(),
                      femtoParticle.eta(),
                      femtoParticle.phi(),
                      femtoParticle.partType(),
                      femtoParticle.cut(),
                      femtoParticle.pidcut(),
                      femtoParticle.tempFitVar(),
                      childIDs,
                      femtoParticle.mLambda(),
                      femtoParticle.mAntiLambda());
        }
        if (aod::femtodreamparticle::ParticleType::kV0 == femtoParticle.partType()) {
          // If the order in primary producer is changed of storing first pos, neg daughters and then V0 - this must be updated
          const int rowOfLastTrack = outputParts.lastIndex();
          std::vector<int> childIDs = {rowOfLastTrack - 1, rowOfLastTrack};
          outputParts(outputCollision.lastIndex(),
                      femtoParticle.pt(),
                      femtoParticle.eta(),
                      femtoParticle.phi(),
                      femtoParticle.partType(),
                      femtoParticle.cut(),
                      femtoParticle.pidcut(),
                      femtoParticle.tempFitVar(),
                      childIDs,
                      femtoParticle.mLambda(),
                      femtoParticle.mAntiLambda());
        }
      }
    } else {
      EventRegistry.fill(HIST("hStatistiscs"), 2);
    }
  }

  /// process function to create derived data with only collisions containing n tracks
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processCollisionsWithNTracksAndNV0(o2::aod::FDCollision& col,
                                          o2::aod::FDParticles& parts)
  {
    EventRegistry.fill(HIST("hStatistiscs"), 0);
    auto thegroupSelectedParts = SelectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupSelectedV0s = SelectedV0s->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

    createSpecifiedDerivedData<false>(col, thegroupSelectedParts, thegroupSelectedV0s, parts);
  }
  PROCESS_SWITCH(femtoDreamProducerTaskForSpecificAnalysis, processCollisionsWithNTracksAndNV0, "Enable producing data with ppp collisions for data", true);

  /// This function stores accepted collisions in derived data
  /// @tparam PartitionType
  /// @tparam PartType
  /// @tparam isMC: enables Monte Carlo truth specific histograms
  /// @param groupSelectedTracks partition for the first particle passed by the process function
  /// @param groupSelectedV0s partition for the second particle passed by the process function
  /// @param parts femtoDreamParticles table
  template <bool isMC, typename PartitionType, typename PartType>
  void createSpecifiedDerivedData_TrkCascade(o2::aod::FDCollision& col, PartitionType groupSelectedTracks, PartitionType groupSelectedCascades, PartType parts)
  {

    /// check tracks
    int tracksCount = 0;
    int antitracksCount = 0;
    for (auto& part : groupSelectedTracks) {
      //if (part.cut() & 1) {
      if( (part.cut() & Track1.CutBitAnti) == Track1.CutBitAnti){
        antitracksCount++;
      //} else {
      } else if( (part.cut() & Track1.CutBit) == Track1.CutBit) {
        tracksCount++;
      }
    }

    /// check Cascades
    int CascadeCount = 0;
    int antiCascadeCount = 0;
    for (auto& casc : groupSelectedCascades) {
    /* 
      if ((casc.cut() & kSignPlusMask) == kSignPlusMask) {
        antiCascadeCount++;
      } else {
        CascadeCount++;
      }
    */
      if ( (casc.cut() & Cascade2.CutBitAnti) == Cascade2.CutBitAnti){
        antiCascadeCount++;
      } else if( (casc.cut() & Cascade2.CutBit) == Cascade2.CutBit){ 
        CascadeCount++;
      }
    }
    std::vector<int> tmpIDtrack;

    if ((CascadeCount >= ConfNumberOfCascades && tracksCount >= ConfNumberOfTracks) || (antiCascadeCount >= ConfNumberOfCascades && antitracksCount >= ConfNumberOfTracks)) {
      EventRegistry.fill(HIST("hStatistiscs"), 1);
      outputCollision(col.posZ(), col.multV0M(), col.multNtr(), col.sphericity(), col.magField());

      for (auto& femtoParticle : parts) {
        if (aod::femtodreamparticle::ParticleType::kTrack == femtoParticle.partType()) {
          std::vector<int> childIDs = {0, 0};
          outputParts(outputCollision.lastIndex(),
                      femtoParticle.pt(),
                      femtoParticle.eta(),
                      femtoParticle.phi(),
                      femtoParticle.partType(),
                      femtoParticle.cut(),
                      femtoParticle.pidcut(),
                      femtoParticle.tempFitVar(),
                      childIDs,
                      femtoParticle.mLambda(),
                      femtoParticle.mAntiLambda());
          tmpIDtrack.push_back(femtoParticle.index());
        }
        
        if (aod::femtodreamparticle::ParticleType::kCascade == femtoParticle.partType()) {
          std::vector<int> childIDs = {0, 0, 0};
          const auto& posChild = parts.iteratorAt(femtoParticle.index() - 3);
          const auto& negChild = parts.iteratorAt(femtoParticle.index() - 2);
          const auto& bachChild = parts.iteratorAt(femtoParticle.index() - 1);
          
          if(femtoParticle.cut() == Cascade2.CutBit){
            
            if (Cascade2.UseChildCuts) {
              if (!(((posChild.cut() & Cascade2.ChildPos_CutBit) == Cascade2.ChildPos_CutBit) &&
                    ((negChild.cut() & Cascade2.ChildNeg_CutBit) == Cascade2.ChildNeg_CutBit) &&
                    ((bachChild.cut() & Cascade2.ChildBach_CutBit) == Cascade2.ChildBach_CutBit))) {
                continue;
              }
            }
            if (Cascade2.UseChildPIDCuts) {
              if (!(((posChild.pidcut() & Cascade2.ChildPos_TPCBit) == Cascade2.ChildPos_TPCBit) &&
                    ((negChild.pidcut() & Cascade2.ChildNeg_TPCBit) == Cascade2.ChildNeg_TPCBit) &&
                    ((bachChild.pidcut() & Cascade2.ChildBach_TPCBit) == Cascade2.ChildBach_TPCBit))) {
                continue;
              }
            }
          
          } else if (femtoParticle.cut() == Cascade2.CutBitAnti){
            
            if (Cascade2.UseChildCuts) {
              if (!(((posChild.cut() & Cascade2.ChildPos_CutBitAnti) == Cascade2.ChildPos_CutBitAnti) &&
                    ((negChild.cut() & Cascade2.ChildNeg_CutBitAnti) == Cascade2.ChildNeg_CutBitAnti) &&
                    ((bachChild.cut() & Cascade2.ChildBach_CutBitAnti) == Cascade2.ChildBach_CutBitAnti))) {
                continue;
              }
            }
            if (Cascade2.UseChildPIDCuts) {
              if (!(((posChild.pidcut() & Cascade2.ChildPos_TPCBitAnti) == Cascade2.ChildPos_TPCBitAnti) &&
                    ((negChild.pidcut() & Cascade2.ChildNeg_TPCBitAnti) == Cascade2.ChildNeg_TPCBitAnti) &&
                    ((bachChild.pidcut() & Cascade2.ChildBach_TPCBitAnti) == Cascade2.ChildBach_TPCBitAnti))) {
                continue;
              }
            }

          }
          
          //Fill positive daughter
          int rowInPrimaryTrackTablePosCasc = -1;
          rowInPrimaryTrackTablePosCasc = posChild.childrenIds()[0]; //getRowDaughters(posChild.childrenIds()[0], tmpIDtrack);
          childIDs[0] = rowInPrimaryTrackTablePosCasc;
          childIDs[1] = 0;
          childIDs[2] = 0;
          outputParts(outputCollision.lastIndex(),
                      posChild.pt(),
                      posChild.eta(),
                      posChild.phi(),
                      posChild.partType(),
                      posChild.cut(),
                      posChild.pidcut(),
                      posChild.tempFitVar(),
                      childIDs,
                      posChild.mLambda(),
                      posChild.mAntiLambda());
          const int rowOfPosCascadeTrack = outputParts.lastIndex();
                      
          // Fill negative daughter
          int rowInPrimaryTrackTableNegCasc = -1;
          rowInPrimaryTrackTableNegCasc = negChild.childrenIds()[1]; //getRowDaughters(negChild.childrenIds()[1], tmpIDtrack);
          childIDs[0] = 0;
          childIDs[1] = rowInPrimaryTrackTableNegCasc;
          childIDs[2] = 0;
          outputParts(outputCollision.lastIndex(),
                      negChild.pt(),
                      negChild.eta(),
                      negChild.phi(),
                      negChild.partType(),
                      negChild.cut(),
                      negChild.pidcut(),
                      negChild.tempFitVar(),
                      childIDs,
                      negChild.mLambda(),
                      negChild.mAntiLambda());
          const int rowOfNegCascadeTrack = outputParts.lastIndex();
          
          // Fill bachelor daughter 
          int rowInPrimaryTrackTableBachCasc = -1;
          rowInPrimaryTrackTableBachCasc = bachChild.childrenIds()[2]; //getRowDaughters(bachChild.childrenIds()[2], tmpIDtrack);
          childIDs[0] = 0;
          childIDs[1] = 0;
          childIDs[2] = rowInPrimaryTrackTableBachCasc;
          outputParts(outputCollision.lastIndex(),
                      bachChild.pt(),
                      bachChild.eta(),
                      bachChild.phi(),
                      bachChild.partType(),
                      bachChild.cut(),
                      bachChild.pidcut(),
                      bachChild.tempFitVar(),
                      childIDs,
                      bachChild.mLambda(),
                      bachChild.mAntiLambda());
          const int rowOfBachCascadeTrack = outputParts.lastIndex();


          // If the order in primary producer is changed of storing first pos, neg daughters and then V0 - this must be updated
          //const int rowOfLastTrack = outputParts.lastIndex();
          std::vector<int> CascchildIDs = {rowOfPosCascadeTrack, rowOfNegCascadeTrack, rowOfBachCascadeTrack};
          outputParts(outputCollision.lastIndex(),
                      femtoParticle.pt(),
                      femtoParticle.eta(),
                      femtoParticle.phi(),
                      femtoParticle.partType(),
                      femtoParticle.cut(),
                      femtoParticle.pidcut(),
                      femtoParticle.tempFitVar(),
                      CascchildIDs,
                      femtoParticle.mLambda(),
                      femtoParticle.mAntiLambda());
        }
      }
    } else {
      EventRegistry.fill(HIST("hStatistiscs"), 2);
    }
  }

  /// process function to create derived data with only collisions containing n tracks
  /// \param col subscribe to the collision table (Data)
  /// \param parts subscribe to the femtoDreamParticleTable
  void processCollisionsWithNTracksAndNCascades(o2::aod::FDCollision& col,
                                                o2::aod::FDParticles& parts)
  {
    EventRegistry.fill(HIST("hStatistiscs"), 0);
    auto thegroupSelectedParts = SelectedParts->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);
    auto thegroupSelectedCascades = SelectedCascades->sliceByCached(aod::femtodreamparticle::fdCollisionId, col.globalIndex(), cache);

    createSpecifiedDerivedData_TrkCascade<false>(col, thegroupSelectedParts, thegroupSelectedCascades, parts);
  }
  PROCESS_SWITCH(femtoDreamProducerTaskForSpecificAnalysis, processCollisionsWithNTracksAndNCascades, "Enable producing data with tracks and Cascades collisions for data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoDreamProducerTaskForSpecificAnalysis>(cfgc),
  };
  return workflow;
}
