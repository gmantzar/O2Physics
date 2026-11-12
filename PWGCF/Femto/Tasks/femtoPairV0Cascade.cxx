// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoPairV0Cascade.cxx
/// \brief Tasks that computes correlation between v0s and cascades
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/cascadeBuilder.h"
#include "PWGCF/Femto/Core/cascadeHistManager.h"
#include "PWGCF/Femto/Core/closePairRejection.h"
#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/pairBuilder.h"
#include "PWGCF/Femto/Core/pairCleaner.h"
#include "PWGCF/Femto/Core/pairHistManager.h"
#include "PWGCF/Femto/Core/pairProcessHelpers.h"
#include "PWGCF/Femto/Core/particleCleaner.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/v0Builder.h"
#include "PWGCF/Femto/Core/v0HistManager.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/Expressions.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <map>
#include <vector>

using namespace o2::analysis::femto;

// The femto core has no v0-cascade pair machinery, so the close pair rejection, the pair cleaner
// and the pair builder for v0-cascade pairs are defined here. They follow the v0-v0 and
// track-cascade implementations in closePairRejection.h, pairCleaner.h and pairBuilder.h.
namespace v0cascade
{

// close pair rejection
constexpr const char PrefixCprV0DaughterCascadePosDau[] = "CprV0DaughterCascadePosDau";
constexpr const char PrefixCprV0DaughterCascadeNegDau[] = "CprV0DaughterCascadeNegDau";
constexpr const char PrefixCprV0DaughterCascadeBachelor[] = "CprV0DaughterCascadeBachelor";

using ConfCprV0DaughterCascadePosDau = closepairrejection::ConfCpr<PrefixCprV0DaughterCascadePosDau>;
using ConfCprV0DaughterCascadeNegDau = closepairrejection::ConfCpr<PrefixCprV0DaughterCascadeNegDau>;
using ConfCprV0DaughterCascadeBachelor = closepairrejection::ConfCpr<PrefixCprV0DaughterCascadeBachelor>;

constexpr char PrefixCprV0CascadePosSe[] = "CPR_V0Cascade_PosDau/SE/";
constexpr char PrefixCprV0CascadeNegSe[] = "CPR_V0Cascade_NegDau/SE/";
constexpr char PrefixCprV0CascadeBachelorSe[] = "CPR_V0Cascade_Bachelor/SE/";
constexpr char PrefixCprV0CascadePosMe[] = "CPR_V0Cascade_PosDau/ME/";
constexpr char PrefixCprV0CascadeNegMe[] = "CPR_V0Cascade_NegDau/ME/";
constexpr char PrefixCprV0CascadeBachelorMe[] = "CPR_V0Cascade_Bachelor/ME/";

// pair histograms
constexpr char PrefixV0CascadeSe[] = "V0Cascade/SE/";
constexpr char PrefixV0CascadeMe[] = "V0Cascade/ME/";

// pair cleaner
constexpr char PrefixPairCleanerV0CascadeSe[] = "V0CascadeCleaner/SE/";
constexpr char PrefixPairCleanerV0CascadeMe[] = "V0CascadeCleaner/ME/";

// checks the v0 daughters against the daughters of the cascade:
//   positive v0 daughter vs positive cascade daughter
//   negative v0 daughter vs negative cascade daughter
//   v0 daughter with the same charge as the bachelor vs bachelor
template <auto& prefixPosDaus, auto& prefixNegDaus, auto& prefixBachelor>
class ClosePairRejectionV0Cascade
{
 public:
  template <typename T1, typename T2, typename T3>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> const& specsPos,
            std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> const& specsNeg,
            std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> const& specsBachelor,
            T1 const& confCprPos,
            T2 const& confCprNeg,
            T3 const& confCprBachelor)
  {
    mCtrPos.init(registry, specsPos, confCprPos, 1, 1);
    mCtrNeg.init(registry, specsNeg, confCprNeg, 1, 1);
    mCtrBachelor.init(registry, specsBachelor, confCprBachelor, 1, 1);
  }

  void setMagField(float magField)
  {
    mCtrPos.setMagField(magField);
    mCtrNeg.setMagField(magField);
    mCtrBachelor.setMagField(magField);
  }

  template <typename T1, typename T2, typename T3, typename T4>
  [[nodiscard]] bool isClosePair(T1 const& v0, T2 const& cascade, T3 const& tracks, T4 const& pairHistManager)
  {
    auto v0PosDau = tracks.rawIteratorAt(v0.posDauId() - tracks.offset());
    auto v0NegDau = tracks.rawIteratorAt(v0.negDauId() - tracks.offset());

    auto cascadePosDau = tracks.rawIteratorAt(cascade.posDauId() - tracks.offset());
    bool isClosePos = mCtrPos.isClosePair(v0PosDau, cascadePosDau, pairHistManager);

    auto cascadeNegDau = tracks.rawIteratorAt(cascade.negDauId() - tracks.offset());
    bool isCloseNeg = mCtrNeg.isClosePair(v0NegDau, cascadeNegDau, pairHistManager);

    bool isCloseBachelor = false;
    auto bachelor = tracks.rawIteratorAt(cascade.bachelorId() - tracks.offset());
    if (bachelor.signedPt() > 0.f) {
      isCloseBachelor = mCtrBachelor.isClosePair(v0PosDau, bachelor, pairHistManager);
    } else {
      isCloseBachelor = mCtrBachelor.isClosePair(v0NegDau, bachelor, pairHistManager);
    }

    return isClosePos || isCloseNeg || isCloseBachelor;
  }

 private:
  closepairrejection::CloseTrackRejection<prefixPosDaus> mCtrPos;
  closepairrejection::CloseTrackRejection<prefixNegDaus> mCtrNeg;
  closepairrejection::CloseTrackRejection<prefixBachelor> mCtrBachelor;
};

// rejects pairs where a v0 daughter is also a daughter (or the bachelor) of the cascade
template <auto& prefix>
class V0CascadePairCleaner : public paircleaner::BasePairCleaner<prefix>
{
 public:
  V0CascadePairCleaner() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  bool isCleanPair(T1 const& v0, T2 const& cascade, T3 const& trackTable, T4 const& pairHistManager) const
  {
    this->fillAll();
    auto v0PosDaughter = trackTable.rawIteratorAt(v0.posDauId() - trackTable.offset());
    auto v0NegDaughter = trackTable.rawIteratorAt(v0.negDauId() - trackTable.offset());
    auto bachelor = trackTable.rawIteratorAt(cascade.bachelorId() - trackTable.offset());
    auto cascadePosDaughter = trackTable.rawIteratorAt(cascade.posDauId() - trackTable.offset());
    auto cascadeNegDaughter = trackTable.rawIteratorAt(cascade.negDauId() - trackTable.offset());
    bool isClean = this->isCleanParticlePair(v0PosDaughter, cascadePosDaughter) && this->isCleanParticlePair(v0NegDaughter, cascadeNegDaughter) &&
                   this->isCleanParticlePair(v0PosDaughter, cascadeNegDaughter) && this->isCleanParticlePair(v0NegDaughter, cascadePosDaughter) &&
                   this->isCleanParticlePair(v0PosDaughter, bachelor) && this->isCleanParticlePair(v0NegDaughter, bachelor);
    if (!isClean) {
      this->fillBlocked(pairHistManager.getKinematic());
    }
    return isClean;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  bool isCleanPair(T1 const& v0, T2 const& cascade, T3 const& trackTable, T4 const& mcParticles, T5 const& partonicMothers, T6 const& pairHistManager) const
  {
    if (!this->isCleanPair(v0, cascade, trackTable, pairHistManager)) {
      return false;
    }
    // pair is clean
    // now check if we require common or non-common ancestry
    if (this->mMixPairsWithCommonAncestor) {
      return this->pairHasCommonAncestor(v0, cascade, mcParticles, partonicMothers);
    }
    if (this->mMixPairsWithNonCommonAncestor) {
      return this->pairHasNonCommonAncestor(v0, cascade, mcParticles, partonicMothers);
    }
    return true;
  }
};

template <auto& prefixV0,
          auto& prefixV0PosDau,
          auto& prefixV0NegDau,
          auto& prefixCascade,
          auto& prefixBachelor,
          auto& prefixCascadePosDau,
          auto& prefixCascadeNegDau,
          auto& prefixV0Me,
          auto& prefixV0PosDauMe,
          auto& prefixV0NegDauMe,
          auto& prefixCascadeMe,
          auto& prefixBachelorMe,
          auto& prefixCascadePosDauMe,
          auto& prefixCascadeNegDauMe,
          auto& prefixSe,
          auto& prefixMe,
          auto& prefixCprPosSe,
          auto& prefixCprNegSe,
          auto& prefixCprBachelorSe,
          auto& prefixCprPosMe,
          auto& prefixCprNegMe,
          auto& prefixCprBachelorMe,
          modes::V0 v0Type,
          modes::Cascade cascadeType>
class PairV0CascadeBuilder
{
 public:
  PairV0CascadeBuilder() = default;
  ~PairV0CascadeBuilder() = default;

  template <modes::Mode modeSe,
            modes::Mode modeMe,
            typename T1,
            typename T2,
            typename T3,
            typename T4,
            typename T5,
            typename T6,
            typename T7,
            typename T8,
            typename T9,
            typename T10,
            typename T11,
            typename T12,
            typename T13,
            typename T14,
            typename T15,
            typename T16,
            typename T17,
            typename T18,
            typename T19,
            typename T20,
            typename T21,
            typename T22,
            typename T23,
            typename T24>
  void init(o2::framework::HistogramRegistry* registry,
            T1 const& confCollisionBinning,
            T2 const& confV0Selection,
            T3 const& confV0Cleaner,
            T4 const& confCascadeSelection,
            T5 const& confCascadeCleaner,
            T6 const& confCprPos,
            T7 const& confCprNeg,
            T8 const& confCprBachelor,
            T9 const& confMixing,
            T10 const& confPairBinning,
            T11 const& confPairCuts,
            std::map<T12, std::vector<o2::framework::AxisSpec>> const& colHistSpec,
            std::map<T13, std::vector<o2::framework::AxisSpec>> const& v0HistSpec,
            std::map<T14, std::vector<o2::framework::AxisSpec>> const& v0PosDauHistSpec,
            std::map<T15, std::vector<o2::framework::AxisSpec>> const& v0NegDauHistSpec,
            std::map<T16, std::vector<o2::framework::AxisSpec>> const& cascadeHistSpec,
            std::map<T17, std::vector<o2::framework::AxisSpec>> const& bachelorHistSpec,
            std::map<T18, std::vector<o2::framework::AxisSpec>> const& cascadePosDauHistSpec,
            std::map<T19, std::vector<o2::framework::AxisSpec>> const& cascadeNegDauHistSpec,
            std::map<T20, std::vector<o2::framework::AxisSpec>> const& pairHistSpec,
            std::map<T21, std::vector<o2::framework::AxisSpec>> const& cprHistSpecPos,
            std::map<T22, std::vector<o2::framework::AxisSpec>> const& cprHistSpecNeg,
            std::map<T23, std::vector<o2::framework::AxisSpec>> const& cprHistSpecBachelor,
            std::map<T24, std::vector<o2::framework::AxisSpec>> const& pairCleanerHistSpec)
  {
    // a v0 and a cascade can never be identical particles
    if (confMixing.sameSpecies.value) {
      LOG(warn) << "sameSpecies=true is impossible for v0-cascade pairs. Option is ignored.";
    }

    mColHistManager.template init<modeSe>(registry, colHistSpec, confCollisionBinning);

    mV0HistManagerSe.template init<modeSe>(registry, v0HistSpec, confV0Selection, v0PosDauHistSpec, v0NegDauHistSpec);
    mCascadeHistManagerSe.template init<modeSe>(registry, cascadeHistSpec, confCascadeSelection, bachelorHistSpec, cascadePosDauHistSpec, cascadeNegDauHistSpec);
    mV0HistManagerMe.template init<modeMe>(registry, v0HistSpec, confV0Selection, v0PosDauHistSpec, v0NegDauHistSpec);
    mCascadeHistManagerMe.template init<modeMe>(registry, cascadeHistSpec, confCascadeSelection, bachelorHistSpec, cascadePosDauHistSpec, cascadeNegDauHistSpec);

    mV0Cleaner.init(confV0Cleaner);
    mCascadeCleaner.init(confCascadeCleaner);

    mPairHistManagerSe.template init<modeSe>(registry, pairHistSpec, confPairBinning, confPairCuts, confMixing);
    mPairHistManagerSe.setMass(confV0Selection.pdgCodeAbs.value, confCascadeSelection.pdgCodeAbs.value);
    mPairHistManagerSe.setCharge(1, 1);
    mCprSe.init(registry, cprHistSpecPos, cprHistSpecNeg, cprHistSpecBachelor, confCprPos, confCprNeg, confCprBachelor);

    mPairHistManagerMe.template init<modeMe>(registry, pairHistSpec, confPairBinning, confPairCuts, confMixing);
    mPairHistManagerMe.setMass(confV0Selection.pdgCodeAbs.value, confCascadeSelection.pdgCodeAbs.value);
    mPairHistManagerMe.setCharge(1, 1);
    mCprMe.init(registry, cprHistSpecPos, cprHistSpecNeg, cprHistSpecBachelor, confCprPos, confCprNeg, confCprBachelor);

    mPcSe.template init<modeSe>(registry, pairCleanerHistSpec, confPairCuts);
    mPcMe.template init<modeMe>(registry, pairCleanerHistSpec, confPairCuts);

    // setup mixing
    mMixingPolicy = static_cast<pairhistmanager::MixingPolicy>(confMixing.policy.value);
    mMixingDepth = confMixing.depth.value;
    mEventStride = confMixing.eventStride.value;
  }

  /// \return true if at least one pair passed all pair selections (usable as pair trigger)
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5>
  bool processSameEvent(T1 const& col, T2 const& trackTable, T3& v0Partition, T4& cascadePartition, T5& cache)
  {
    auto v0Slice = v0Partition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto cascadeSlice = cascadePartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (v0Slice.size() < pairbuilder::nLimitPartitionParticles || cascadeSlice.size() < pairbuilder::nLimitPartitionParticles) {
      return false;
    }
    mColHistManager.template fill<mode>(col);
    mCprSe.setMagField(col.magField());
    return pairprocesshelpers::processSameEvent<mode>(v0Slice, cascadeSlice, trackTable, col, mV0HistManagerSe, mCascadeHistManagerSe, mPairHistManagerSe, mCprSe, mPcSe);
  }

  // mc
  /// \return true if at least one pair passed all pair selections (usable as pair trigger)
  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  bool processSameEvent(T1 const& col, T2 const& mcCols, T3 const& trackTable, T4& v0Partition, T5& cascadePartition, T6 const& mcParticles, T7 const& mcMothers, T8 const& mcPartonicMothers, T9& cache)
  {
    auto v0Slice = v0Partition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto cascadeSlice = cascadePartition->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (v0Slice.size() < pairbuilder::nLimitPartitionParticles || cascadeSlice.size() < pairbuilder::nLimitPartitionParticles) {
      return false;
    }
    mColHistManager.template fill<mode>(col, mcCols);
    mCprSe.setMagField(col.magField());
    return pairprocesshelpers::processSameEvent<mode>(v0Slice, cascadeSlice, trackTable, mcParticles, mcMothers, mcPartonicMothers, col, mcCols, mV0HistManagerSe, mCascadeHistManagerSe, mPairHistManagerSe, mV0Cleaner, mCascadeCleaner, mCprSe, mPcSe);
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void processMixedEvent(T1 const& cols, T2 const& trackTable, T3& v0Partition, T4& cascadePartition, T5& cache, T6& binsVtxMult, T7& binsVtxCent, T8& binsVtxMultCent)
  {
    switch (mMixingPolicy) {
      case static_cast<int>(pairhistmanager::kVtxMult):
        pairprocesshelpers::processMixedEvent<mode>(cols, v0Partition, cascadePartition, trackTable, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mCprMe, mPcMe, mEventStride, mV0HistManagerMe, mCascadeHistManagerMe);
        break;
      case static_cast<int>(pairhistmanager::kVtxCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, v0Partition, cascadePartition, trackTable, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPcMe, mEventStride, mV0HistManagerMe, mCascadeHistManagerMe);
        break;
      case static_cast<int>(pairhistmanager::kVtxMultCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, v0Partition, cascadePartition, trackTable, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mCprMe, mPcMe, mEventStride, mV0HistManagerMe, mCascadeHistManagerMe);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12>
  void processMixedEvent(T1 const& cols, T2 const& mcCols, T3& trackTable, T4& v0Partition, T5& cascadePartition, T6 const& mcParticles, T7 const& mcMothers, T8 const& mcPartonicMothers, T9& cache, T10& binsVtxMult, T11& binsVtxCent, T12& binsVtxMultCent)
  {
    switch (mMixingPolicy) {
      case static_cast<int>(pairhistmanager::kVtxMult):
        pairprocesshelpers::processMixedEvent<mode>(cols, mcCols, v0Partition, cascadePartition, trackTable, mcParticles, mcMothers, mcPartonicMothers, cache, binsVtxMult, mMixingDepth, mPairHistManagerMe, mV0Cleaner, mCascadeCleaner, mCprMe, mPcMe, mEventStride, mV0HistManagerMe, mCascadeHistManagerMe);
        break;
      case static_cast<int>(pairhistmanager::kVtxCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, mcCols, v0Partition, cascadePartition, trackTable, mcParticles, mcMothers, mcPartonicMothers, cache, binsVtxCent, mMixingDepth, mPairHistManagerMe, mV0Cleaner, mCascadeCleaner, mCprMe, mPcMe, mEventStride, mV0HistManagerMe, mCascadeHistManagerMe);
        break;
      case static_cast<int>(pairhistmanager::kVtxMultCent):
        pairprocesshelpers::processMixedEvent<mode>(cols, mcCols, v0Partition, cascadePartition, trackTable, mcParticles, mcMothers, mcPartonicMothers, cache, binsVtxMultCent, mMixingDepth, mPairHistManagerMe, mV0Cleaner, mCascadeCleaner, mCprMe, mPcMe, mEventStride, mV0HistManagerMe, mCascadeHistManagerMe);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }

 private:
  colhistmanager::CollisionHistManager mColHistManager;
  v0histmanager::V0HistManager<prefixV0, prefixV0PosDau, prefixV0NegDau, v0Type> mV0HistManagerSe;
  cascadehistmanager::CascadeHistManager<prefixCascade, prefixBachelor, prefixCascadePosDau, prefixCascadeNegDau, cascadeType> mCascadeHistManagerSe;
  v0histmanager::V0HistManager<prefixV0Me, prefixV0PosDauMe, prefixV0NegDauMe, v0Type> mV0HistManagerMe;
  cascadehistmanager::CascadeHistManager<prefixCascadeMe, prefixBachelorMe, prefixCascadePosDauMe, prefixCascadeNegDauMe, cascadeType> mCascadeHistManagerMe;
  particlecleaner::ParticleCleaner mV0Cleaner;
  particlecleaner::ParticleCleaner mCascadeCleaner;
  pairhistmanager::PairHistManager<prefixSe, modes::Particle::kV0, modes::Particle::kCascade> mPairHistManagerSe;
  pairhistmanager::PairHistManager<prefixMe, modes::Particle::kV0, modes::Particle::kCascade> mPairHistManagerMe;
  ClosePairRejectionV0Cascade<prefixCprPosSe, prefixCprNegSe, prefixCprBachelorSe> mCprSe;
  ClosePairRejectionV0Cascade<prefixCprPosMe, prefixCprNegMe, prefixCprBachelorMe> mCprMe;
  V0CascadePairCleaner<PrefixPairCleanerV0CascadeSe> mPcSe;
  V0CascadePairCleaner<PrefixPairCleanerV0CascadeMe> mPcMe;
  pairhistmanager::MixingPolicy mMixingPolicy = pairhistmanager::MixingPolicy::kVtxMult;
  int mMixingDepth = 5;
  int mEventStride = 1;
};

} // namespace v0cascade

// The v0 always occupies slot 1 and the cascade slot 2 of the pair:
//   lambda-xi, lambda-omega, k0short-xi, k0short-omega
// Only one combination can be processed at a time.
struct FemtoPairV0Cascade {

  // setup tables
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoCollisionsWithLabel = o2::soa::Join<FemtoCollisions, o2::aod::FColLabels>;
  using FilteredFemtoCollisionsWithLabel = o2::soa::Filtered<FemtoCollisionsWithLabel>;
  using FilteredFemtoCollisionWithLabel = FilteredFemtoCollisionsWithLabel::iterator;

  using FemtoTracks = o2::aod::FTracks;
  using FemtoLambdas = o2::soa::Join<o2::aod::FLambdas, o2::aod::FLambdaMasks>;
  using FemtoK0shorts = o2::soa::Join<o2::aod::FK0shorts, o2::aod::FK0shortMasks>;
  using FemtoXis = o2::soa::Join<o2::aod::FXis, o2::aod::FXiMasks>;
  using FemtoOmegas = o2::soa::Join<o2::aod::FOmegas, o2::aod::FOmegaMasks>;

  using FemtoTracksWithLabel = o2::soa::Join<FemtoTracks, o2::aod::FTrackLabels>;
  using FemtoLambdasWithLabel = o2::soa::Join<FemtoLambdas, o2::aod::FLambdaLabels>;
  using FemtoK0shortsWithLabel = o2::soa::Join<FemtoK0shorts, o2::aod::FK0shortLabels>;
  using FemtoXisWithLabel = o2::soa::Join<FemtoXis, o2::aod::FXiLabels>;
  using FemtoOmegasWithLabel = o2::soa::Join<FemtoOmegas, o2::aod::FOmegaLabels>;

  using FemtoMcParticlesWithLabel = o2::soa::Join<o2::aod::FMcParticles, o2::aod::FMcMotherLabels>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // setup for v0 daughters
  trackhistmanager::ConfV01PosDauBinning confV0PosDauBinning;
  trackhistmanager::ConfV01NegDauBinning confV0NegDauBinning;

  // setup for cascade daughters/bachelor
  trackhistmanager::ConfCascadePosDauBinning confCascadePosDauBinning;
  trackhistmanager::ConfCascadeNegDauBinning confCascadeNegDauBinning;
  trackhistmanager::ConfCascadeBachelorBinning confBachelorBinning;

  // setup lambdas
  v0builder::ConfLambdaSelection1 confLambdaSelection;
  particlecleaner::ConfLambdaCleaner1 confLambdaCleaner;
  v0histmanager::ConfLambdaBinning1 confLambdaBinning;

  o2::framework::Partition<FemtoLambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  o2::framework::Preslice<FemtoLambdas> perColLambdas = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoLambdasWithLabel> lambdaWithLabelPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  o2::framework::Preslice<FemtoLambdasWithLabel> perColLambdasWithLabel = o2::aod::femtobase::stored::fColId;

  // setup k0shorts
  v0builder::ConfK0shortSelection1 confK0shortSelection;
  particlecleaner::ConfK0shortCleaner1 confK0shortCleaner;
  v0histmanager::ConfK0shortBinning1 confK0shortBinning;

  o2::framework::Partition<FemtoK0shorts> k0shortPartition = MAKE_K0SHORT_PARTITION(confK0shortSelection);
  o2::framework::Preslice<FemtoK0shorts> perColK0shorts = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoK0shortsWithLabel> k0shortWithLabelPartition = MAKE_K0SHORT_PARTITION(confK0shortSelection);
  o2::framework::Preslice<FemtoK0shortsWithLabel> perColK0shortsWithLabel = o2::aod::femtobase::stored::fColId;

  // setup xis
  cascadebuilder::ConfXiSelection confXiSelection;
  particlecleaner::ConfXiCleaner1 confXiCleaner;
  cascadehistmanager::ConfXiBinning confXiBinning;

  o2::framework::Partition<FemtoXis> xiPartition = MAKE_CASCADE_PARTITION(confXiSelection);
  o2::framework::Preslice<FemtoXis> perColXis = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoXisWithLabel> xiWithLabelPartition = MAKE_CASCADE_PARTITION(confXiSelection);
  o2::framework::Preslice<FemtoXisWithLabel> perColXisWithLabel = o2::aod::femtobase::stored::fColId;

  // setup omegas
  cascadebuilder::ConfOmegaSelection confOmegaSelection;
  particlecleaner::ConfOmegaCleaner1 confOmegaCleaner;
  cascadehistmanager::ConfOmegaBinning confOmegaBinning;

  o2::framework::Partition<FemtoOmegas> omegaPartition = MAKE_CASCADE_PARTITION(confOmegaSelection);
  o2::framework::Preslice<FemtoOmegas> perColOmegas = o2::aod::femtobase::stored::fColId;

  o2::framework::Partition<FemtoOmegasWithLabel> omegaWithLabelPartition = MAKE_CASCADE_PARTITION(confOmegaSelection);
  o2::framework::Preslice<FemtoOmegasWithLabel> perColOmegasWithLabel = o2::aod::femtobase::stored::fColId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::ConfPairCuts confPairCuts;
  paircleaner::ConfPairCleanerBinning confPairCleaner;

  v0cascade::PairV0CascadeBuilder<
    v0histmanager::PrefixLambda1Se,
    trackhistmanager::PrefixV01PosDaughterSe,
    trackhistmanager::PrefixV01NegDaughterSe,
    cascadehistmanager::PrefixXiSe,
    trackhistmanager::PrefixCascadeBachelorSe,
    trackhistmanager::PrefixCascadePosDaughterSe,
    trackhistmanager::PrefixCascadeNegDaughterSe,
    v0histmanager::PrefixLambda1Me,
    trackhistmanager::PrefixV01PosDaughterMe,
    trackhistmanager::PrefixV01NegDaughterMe,
    cascadehistmanager::PrefixXiMe,
    trackhistmanager::PrefixCascadeBachelorMe,
    trackhistmanager::PrefixCascadePosDaughterMe,
    trackhistmanager::PrefixCascadeNegDaughterMe,
    v0cascade::PrefixV0CascadeSe,
    v0cascade::PrefixV0CascadeMe,
    v0cascade::PrefixCprV0CascadePosSe,
    v0cascade::PrefixCprV0CascadeNegSe,
    v0cascade::PrefixCprV0CascadeBachelorSe,
    v0cascade::PrefixCprV0CascadePosMe,
    v0cascade::PrefixCprV0CascadeNegMe,
    v0cascade::PrefixCprV0CascadeBachelorMe,
    modes::V0::kLambda,
    modes::Cascade::kXi>
    pairLambdaXiBuilder;

  v0cascade::PairV0CascadeBuilder<
    v0histmanager::PrefixLambda1Se,
    trackhistmanager::PrefixV01PosDaughterSe,
    trackhistmanager::PrefixV01NegDaughterSe,
    cascadehistmanager::PrefixOmegaSe,
    trackhistmanager::PrefixCascadeBachelorSe,
    trackhistmanager::PrefixCascadePosDaughterSe,
    trackhistmanager::PrefixCascadeNegDaughterSe,
    v0histmanager::PrefixLambda1Me,
    trackhistmanager::PrefixV01PosDaughterMe,
    trackhistmanager::PrefixV01NegDaughterMe,
    cascadehistmanager::PrefixOmegaMe,
    trackhistmanager::PrefixCascadeBachelorMe,
    trackhistmanager::PrefixCascadePosDaughterMe,
    trackhistmanager::PrefixCascadeNegDaughterMe,
    v0cascade::PrefixV0CascadeSe,
    v0cascade::PrefixV0CascadeMe,
    v0cascade::PrefixCprV0CascadePosSe,
    v0cascade::PrefixCprV0CascadeNegSe,
    v0cascade::PrefixCprV0CascadeBachelorSe,
    v0cascade::PrefixCprV0CascadePosMe,
    v0cascade::PrefixCprV0CascadeNegMe,
    v0cascade::PrefixCprV0CascadeBachelorMe,
    modes::V0::kLambda,
    modes::Cascade::kOmega>
    pairLambdaOmegaBuilder;

  v0cascade::PairV0CascadeBuilder<
    v0histmanager::PrefixK0short1Se,
    trackhistmanager::PrefixV01PosDaughterSe,
    trackhistmanager::PrefixV01NegDaughterSe,
    cascadehistmanager::PrefixXiSe,
    trackhistmanager::PrefixCascadeBachelorSe,
    trackhistmanager::PrefixCascadePosDaughterSe,
    trackhistmanager::PrefixCascadeNegDaughterSe,
    v0histmanager::PrefixK0short1Me,
    trackhistmanager::PrefixV01PosDaughterMe,
    trackhistmanager::PrefixV01NegDaughterMe,
    cascadehistmanager::PrefixXiMe,
    trackhistmanager::PrefixCascadeBachelorMe,
    trackhistmanager::PrefixCascadePosDaughterMe,
    trackhistmanager::PrefixCascadeNegDaughterMe,
    v0cascade::PrefixV0CascadeSe,
    v0cascade::PrefixV0CascadeMe,
    v0cascade::PrefixCprV0CascadePosSe,
    v0cascade::PrefixCprV0CascadeNegSe,
    v0cascade::PrefixCprV0CascadeBachelorSe,
    v0cascade::PrefixCprV0CascadePosMe,
    v0cascade::PrefixCprV0CascadeNegMe,
    v0cascade::PrefixCprV0CascadeBachelorMe,
    modes::V0::kK0short,
    modes::Cascade::kXi>
    pairK0shortXiBuilder;

  v0cascade::PairV0CascadeBuilder<
    v0histmanager::PrefixK0short1Se,
    trackhistmanager::PrefixV01PosDaughterSe,
    trackhistmanager::PrefixV01NegDaughterSe,
    cascadehistmanager::PrefixOmegaSe,
    trackhistmanager::PrefixCascadeBachelorSe,
    trackhistmanager::PrefixCascadePosDaughterSe,
    trackhistmanager::PrefixCascadeNegDaughterSe,
    v0histmanager::PrefixK0short1Me,
    trackhistmanager::PrefixV01PosDaughterMe,
    trackhistmanager::PrefixV01NegDaughterMe,
    cascadehistmanager::PrefixOmegaMe,
    trackhistmanager::PrefixCascadeBachelorMe,
    trackhistmanager::PrefixCascadePosDaughterMe,
    trackhistmanager::PrefixCascadeNegDaughterMe,
    v0cascade::PrefixV0CascadeSe,
    v0cascade::PrefixV0CascadeMe,
    v0cascade::PrefixCprV0CascadePosSe,
    v0cascade::PrefixCprV0CascadeNegSe,
    v0cascade::PrefixCprV0CascadeBachelorSe,
    v0cascade::PrefixCprV0CascadePosMe,
    v0cascade::PrefixCprV0CascadeNegMe,
    v0cascade::PrefixCprV0CascadeBachelorMe,
    modes::V0::kK0short,
    modes::Cascade::kOmega>
    pairK0shortOmegaBuilder;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult> mixBinsVtxMult{{defaultVtxBins, defaultMultBins}};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Cent> mixBinsVtxCent{{defaultVtxBins, defaultCentBins}};
  o2::framework::ColumnBinningPolicy<o2::aod::femtocollisions::PosZ, o2::aod::femtocollisions::Mult, o2::aod::femtocollisions::Cent> mixBinsVtxMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}};
  pairhistmanager::ConfMixing confMixing;

  o2::framework::HistogramRegistry hRegistry{"FemtoV0Cascade", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  // setup cpr
  v0cascade::ConfCprV0DaughterCascadePosDau confCprPos;
  v0cascade::ConfCprV0DaughterCascadeNegDau confCprNeg;
  v0cascade::ConfCprV0DaughterCascadeBachelor confCprBachelor;

  void init(o2::framework::InitContext&)
  {
    bool processData = doprocessLambdaXiSameEvent || doprocessLambdaXiMixedEvent ||
                       doprocessLambdaOmegaSameEvent || doprocessLambdaOmegaMixedEvent ||
                       doprocessK0shortXiSameEvent || doprocessK0shortXiMixedEvent ||
                       doprocessK0shortOmegaSameEvent || doprocessK0shortOmegaMixedEvent;
    bool processMc = doprocessLambdaXiSameEventMc || doprocessLambdaXiMixedEventMc ||
                     doprocessLambdaOmegaSameEventMc || doprocessLambdaOmegaMixedEventMc ||
                     doprocessK0shortXiSameEventMc || doprocessK0shortXiMixedEventMc ||
                     doprocessK0shortOmegaSameEventMc || doprocessK0shortOmegaMixedEventMc;

    if (processData && processMc) {
      LOG(fatal) << "Both data and mc processing is enabled. Breaking...";
    }
    if (!processData && !processMc) {
      LOG(fatal) << "Neither data nor mc processing is enabled. Breaking...";
    }

    bool processLambdaXi = doprocessLambdaXiSameEvent || doprocessLambdaXiMixedEvent ||
                           doprocessLambdaXiSameEventMc || doprocessLambdaXiMixedEventMc;
    bool processLambdaOmega = doprocessLambdaOmegaSameEvent || doprocessLambdaOmegaMixedEvent ||
                              doprocessLambdaOmegaSameEventMc || doprocessLambdaOmegaMixedEventMc;
    bool processK0shortXi = doprocessK0shortXiSameEvent || doprocessK0shortXiMixedEvent ||
                            doprocessK0shortXiSameEventMc || doprocessK0shortXiMixedEventMc;
    bool processK0shortOmega = doprocessK0shortOmegaSameEvent || doprocessK0shortOmegaMixedEvent ||
                               doprocessK0shortOmegaSameEventMc || doprocessK0shortOmegaMixedEventMc;

    if (static_cast<int>(processLambdaXi) + static_cast<int>(processLambdaOmega) + static_cast<int>(processK0shortXi) + static_cast<int>(processK0shortOmega) > 1) {
      LOG(fatal) << "Only one of lambda-xi, lambda-omega, k0short-xi or k0short-omega processing can be enabled. Breaking...";
    }

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinsVtxMult = {{confMixing.vtxBins.value, confMixing.multBins.value}};
    mixBinsVtxCent = {{confMixing.vtxBins.value, confMixing.centBins.value}};
    mixBinsVtxMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}};

    // setup histograms
    std::map<colhistmanager::ColHist, std::vector<o2::framework::AxisSpec>> colHistSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> v0PosDauSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> v0NegDauSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> cascadePosDauSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> cascadeNegDauSpec;
    std::map<trackhistmanager::TrackHist, std::vector<o2::framework::AxisSpec>> bachelorSpec;
    std::map<v0histmanager::V0Hist, std::vector<o2::framework::AxisSpec>> lambdaHistSpec;
    std::map<v0histmanager::V0Hist, std::vector<o2::framework::AxisSpec>> k0shortHistSpec;
    std::map<cascadehistmanager::CascadeHist, std::vector<o2::framework::AxisSpec>> xiHistSpec;
    std::map<cascadehistmanager::CascadeHist, std::vector<o2::framework::AxisSpec>> omegaHistSpec;
    std::map<pairhistmanager::PairHist, std::vector<o2::framework::AxisSpec>> pairV0CascadeHistSpec;
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpecPos = closepairrejection::makeCprHistSpecMap(confCprPos);
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpecNeg = closepairrejection::makeCprHistSpecMap(confCprNeg);
    std::map<closepairrejection::CprHist, std::vector<o2::framework::AxisSpec>> cprHistSpecBachelor = closepairrejection::makeCprHistSpecMap(confCprBachelor);
    std::map<paircleaner::PairCleanerHist, std::vector<o2::framework::AxisSpec>> pairCleanerHistSpec = paircleaner::makePairCleanerHistSpecMap(confPairCleaner);

    if (processData) {
      colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
      v0PosDauSpec = trackhistmanager::makeTrackHistSpecMap(confV0PosDauBinning);
      v0NegDauSpec = trackhistmanager::makeTrackHistSpecMap(confV0NegDauBinning);
      cascadePosDauSpec = trackhistmanager::makeTrackHistSpecMap(confCascadePosDauBinning);
      cascadeNegDauSpec = trackhistmanager::makeTrackHistSpecMap(confCascadeNegDauBinning);
      bachelorSpec = trackhistmanager::makeTrackHistSpecMap(confBachelorBinning);
      pairV0CascadeHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confMixing);

      if (processLambdaXi) {
        lambdaHistSpec = v0histmanager::makeV0HistSpecMap(confLambdaBinning);
        xiHistSpec = cascadehistmanager::makeCascadeHistSpecMap(confXiBinning);
        pairLambdaXiBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confLambdaSelection, confLambdaCleaner, confXiSelection, confXiCleaner, confCprPos, confCprNeg, confCprBachelor, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec, v0PosDauSpec, v0NegDauSpec, xiHistSpec, bachelorSpec, cascadePosDauSpec, cascadeNegDauSpec, pairV0CascadeHistSpec, cprHistSpecPos, cprHistSpecNeg, cprHistSpecBachelor, pairCleanerHistSpec);
      }

      if (processLambdaOmega) {
        lambdaHistSpec = v0histmanager::makeV0HistSpecMap(confLambdaBinning);
        omegaHistSpec = cascadehistmanager::makeCascadeHistSpecMap(confOmegaBinning);
        pairLambdaOmegaBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confLambdaSelection, confLambdaCleaner, confOmegaSelection, confOmegaCleaner, confCprPos, confCprNeg, confCprBachelor, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec, v0PosDauSpec, v0NegDauSpec, omegaHistSpec, bachelorSpec, cascadePosDauSpec, cascadeNegDauSpec, pairV0CascadeHistSpec, cprHistSpecPos, cprHistSpecNeg, cprHistSpecBachelor, pairCleanerHistSpec);
      }

      if (processK0shortXi) {
        k0shortHistSpec = v0histmanager::makeV0HistSpecMap(confK0shortBinning);
        xiHistSpec = cascadehistmanager::makeCascadeHistSpecMap(confXiBinning);
        pairK0shortXiBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confK0shortSelection, confK0shortCleaner, confXiSelection, confXiCleaner, confCprPos, confCprNeg, confCprBachelor, confMixing, confPairBinning, confPairCuts, colHistSpec, k0shortHistSpec, v0PosDauSpec, v0NegDauSpec, xiHistSpec, bachelorSpec, cascadePosDauSpec, cascadeNegDauSpec, pairV0CascadeHistSpec, cprHistSpecPos, cprHistSpecNeg, cprHistSpecBachelor, pairCleanerHistSpec);
      }

      if (processK0shortOmega) {
        k0shortHistSpec = v0histmanager::makeV0HistSpecMap(confK0shortBinning);
        omegaHistSpec = cascadehistmanager::makeCascadeHistSpecMap(confOmegaBinning);
        pairK0shortOmegaBuilder.init<modes::Mode::kSe_Reco, modes::Mode::kMe_Reco>(&hRegistry, confCollisionBinning, confK0shortSelection, confK0shortCleaner, confOmegaSelection, confOmegaCleaner, confCprPos, confCprNeg, confCprBachelor, confMixing, confPairBinning, confPairCuts, colHistSpec, k0shortHistSpec, v0PosDauSpec, v0NegDauSpec, omegaHistSpec, bachelorSpec, cascadePosDauSpec, cascadeNegDauSpec, pairV0CascadeHistSpec, cprHistSpecPos, cprHistSpecNeg, cprHistSpecBachelor, pairCleanerHistSpec);
      }
    } else {
      colHistSpec = colhistmanager::makeColMcHistSpecMap(confCollisionBinning);
      v0PosDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confV0PosDauBinning);
      v0NegDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confV0NegDauBinning);
      cascadePosDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confCascadePosDauBinning);
      cascadeNegDauSpec = trackhistmanager::makeTrackMcHistSpecMap(confCascadeNegDauBinning);
      bachelorSpec = trackhistmanager::makeTrackMcHistSpecMap(confBachelorBinning);
      pairV0CascadeHistSpec = pairhistmanager::makePairMcHistSpecMap(confPairBinning, confMixing);

      if (processLambdaXi) {
        lambdaHistSpec = v0histmanager::makeV0McHistSpecMap(confLambdaBinning);
        xiHistSpec = cascadehistmanager::makeCascadeMcHistSpecMap(confXiBinning);
        pairLambdaXiBuilder.init<modes::Mode::kSe_Reco_Mc, modes::Mode::kMe_Reco_Mc>(&hRegistry, confCollisionBinning, confLambdaSelection, confLambdaCleaner, confXiSelection, confXiCleaner, confCprPos, confCprNeg, confCprBachelor, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec, v0PosDauSpec, v0NegDauSpec, xiHistSpec, bachelorSpec, cascadePosDauSpec, cascadeNegDauSpec, pairV0CascadeHistSpec, cprHistSpecPos, cprHistSpecNeg, cprHistSpecBachelor, pairCleanerHistSpec);
      }

      if (processLambdaOmega) {
        lambdaHistSpec = v0histmanager::makeV0McHistSpecMap(confLambdaBinning);
        omegaHistSpec = cascadehistmanager::makeCascadeMcHistSpecMap(confOmegaBinning);
        pairLambdaOmegaBuilder.init<modes::Mode::kSe_Reco_Mc, modes::Mode::kMe_Reco_Mc>(&hRegistry, confCollisionBinning, confLambdaSelection, confLambdaCleaner, confOmegaSelection, confOmegaCleaner, confCprPos, confCprNeg, confCprBachelor, confMixing, confPairBinning, confPairCuts, colHistSpec, lambdaHistSpec, v0PosDauSpec, v0NegDauSpec, omegaHistSpec, bachelorSpec, cascadePosDauSpec, cascadeNegDauSpec, pairV0CascadeHistSpec, cprHistSpecPos, cprHistSpecNeg, cprHistSpecBachelor, pairCleanerHistSpec);
      }

      if (processK0shortXi) {
        k0shortHistSpec = v0histmanager::makeV0McHistSpecMap(confK0shortBinning);
        xiHistSpec = cascadehistmanager::makeCascadeMcHistSpecMap(confXiBinning);
        pairK0shortXiBuilder.init<modes::Mode::kSe_Reco_Mc, modes::Mode::kMe_Reco_Mc>(&hRegistry, confCollisionBinning, confK0shortSelection, confK0shortCleaner, confXiSelection, confXiCleaner, confCprPos, confCprNeg, confCprBachelor, confMixing, confPairBinning, confPairCuts, colHistSpec, k0shortHistSpec, v0PosDauSpec, v0NegDauSpec, xiHistSpec, bachelorSpec, cascadePosDauSpec, cascadeNegDauSpec, pairV0CascadeHistSpec, cprHistSpecPos, cprHistSpecNeg, cprHistSpecBachelor, pairCleanerHistSpec);
      }

      if (processK0shortOmega) {
        k0shortHistSpec = v0histmanager::makeV0McHistSpecMap(confK0shortBinning);
        omegaHistSpec = cascadehistmanager::makeCascadeMcHistSpecMap(confOmegaBinning);
        pairK0shortOmegaBuilder.init<modes::Mode::kSe_Reco_Mc, modes::Mode::kMe_Reco_Mc>(&hRegistry, confCollisionBinning, confK0shortSelection, confK0shortCleaner, confOmegaSelection, confOmegaCleaner, confCprPos, confCprNeg, confCprBachelor, confMixing, confPairBinning, confPairCuts, colHistSpec, k0shortHistSpec, v0PosDauSpec, v0NegDauSpec, omegaHistSpec, bachelorSpec, cascadePosDauSpec, cascadeNegDauSpec, pairV0CascadeHistSpec, cprHistSpecPos, cprHistSpecNeg, cprHistSpecBachelor, pairCleanerHistSpec);
      }
    }
    hRegistry.print();
  }

  // lambda-xi
  void processLambdaXiSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoXis const& /*xis*/)
  {
    pairLambdaXiBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, lambdaPartition, xiPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processLambdaXiSameEvent, "Enable processing same event processing for lambda-xi", true);

  void processLambdaXiSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& /*lambdas*/, FemtoXisWithLabel const& /*xis*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairLambdaXiBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc>(col, mcCols, tracks, lambdaWithLabelPartition, xiWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processLambdaXiSameEventMc, "Enable processing same event processing for lambda-xi with mc information", false);

  void processLambdaXiMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoXis const& /*xis*/)
  {
    pairLambdaXiBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, lambdaPartition, xiPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processLambdaXiMixedEvent, "Enable processing mixed event processing for lambda-xi", true);

  void processLambdaXiMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& /*lambdas*/, FemtoXisWithLabel const& /*xis*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairLambdaXiBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc>(cols, mcCols, tracks, lambdaWithLabelPartition, xiWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processLambdaXiMixedEventMc, "Enable processing mixed event processing for lambda-xi with mc information", false);

  // lambda-omega
  void processLambdaOmegaSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoOmegas const& /*omegas*/)
  {
    pairLambdaOmegaBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, lambdaPartition, omegaPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processLambdaOmegaSameEvent, "Enable processing same event processing for lambda-omega", false);

  void processLambdaOmegaSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& /*lambdas*/, FemtoOmegasWithLabel const& /*omegas*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairLambdaOmegaBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc>(col, mcCols, tracks, lambdaWithLabelPartition, omegaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processLambdaOmegaSameEventMc, "Enable processing same event processing for lambda-omega with mc information", false);

  void processLambdaOmegaMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoLambdas const& /*lambdas*/, FemtoOmegas const& /*omegas*/)
  {
    pairLambdaOmegaBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, lambdaPartition, omegaPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processLambdaOmegaMixedEvent, "Enable processing mixed event processing for lambda-omega", false);

  void processLambdaOmegaMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoLambdasWithLabel const& /*lambdas*/, FemtoOmegasWithLabel const& /*omegas*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairLambdaOmegaBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc>(cols, mcCols, tracks, lambdaWithLabelPartition, omegaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processLambdaOmegaMixedEventMc, "Enable processing mixed event processing for lambda-omega with mc information", false);

  // k0short-xi
  void processK0shortXiSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/, FemtoXis const& /*xis*/)
  {
    pairK0shortXiBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, k0shortPartition, xiPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processK0shortXiSameEvent, "Enable processing same event processing for k0short-xi", false);

  void processK0shortXiSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoK0shortsWithLabel const& /*k0shorts*/, FemtoXisWithLabel const& /*xis*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairK0shortXiBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc>(col, mcCols, tracks, k0shortWithLabelPartition, xiWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processK0shortXiSameEventMc, "Enable processing same event processing for k0short-xi with mc information", false);

  void processK0shortXiMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/, FemtoXis const& /*xis*/)
  {
    pairK0shortXiBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, k0shortPartition, xiPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processK0shortXiMixedEvent, "Enable processing mixed event processing for k0short-xi", false);

  void processK0shortXiMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoK0shortsWithLabel const& /*k0shorts*/, FemtoXisWithLabel const& /*xis*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairK0shortXiBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc>(cols, mcCols, tracks, k0shortWithLabelPartition, xiWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processK0shortXiMixedEventMc, "Enable processing mixed event processing for k0short-xi with mc information", false);

  // k0short-omega
  void processK0shortOmegaSameEvent(FilteredFemtoCollision const& col, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/, FemtoOmegas const& /*omegas*/)
  {
    pairK0shortOmegaBuilder.processSameEvent<modes::Mode::kSe_Reco>(col, tracks, k0shortPartition, omegaPartition, cache);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processK0shortOmegaSameEvent, "Enable processing same event processing for k0short-omega", false);

  void processK0shortOmegaSameEventMc(FilteredFemtoCollisionWithLabel const& col, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoK0shortsWithLabel const& /*k0shorts*/, FemtoOmegasWithLabel const& /*omegas*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairK0shortOmegaBuilder.processSameEvent<modes::Mode::kSe_Reco_Mc>(col, mcCols, tracks, k0shortWithLabelPartition, omegaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processK0shortOmegaSameEventMc, "Enable processing same event processing for k0short-omega with mc information", false);

  void processK0shortOmegaMixedEvent(FilteredFemtoCollisions const& cols, FemtoTracks const& tracks, FemtoK0shorts const& /*k0shorts*/, FemtoOmegas const& /*omegas*/)
  {
    pairK0shortOmegaBuilder.processMixedEvent<modes::Mode::kMe_Reco>(cols, tracks, k0shortPartition, omegaPartition, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processK0shortOmegaMixedEvent, "Enable processing mixed event processing for k0short-omega", false);

  void processK0shortOmegaMixedEventMc(FilteredFemtoCollisionsWithLabel const& cols, o2::aod::FMcCols const& mcCols, FemtoTracksWithLabel const& tracks, FemtoK0shortsWithLabel const& /*k0shorts*/, FemtoOmegasWithLabel const& /*omegas*/, FemtoMcParticlesWithLabel const& mcParticles, o2::aod::FMcMothers const& mcMothers, o2::aod::FMcPartMoths const& mcPartonicMothers)
  {
    pairK0shortOmegaBuilder.processMixedEvent<modes::Mode::kMe_Reco_Mc>(cols, mcCols, tracks, k0shortWithLabelPartition, omegaWithLabelPartition, mcParticles, mcMothers, mcPartonicMothers, cache, mixBinsVtxMult, mixBinsVtxCent, mixBinsVtxMultCent);
  }
  PROCESS_SWITCH(FemtoPairV0Cascade, processK0shortOmegaMixedEventMc, "Enable processing mixed event processing for k0short-omega with mc information", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    adaptAnalysisTask<FemtoPairV0Cascade>(context),
  };
  return workflow;
}
