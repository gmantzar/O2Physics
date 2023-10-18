/// \brief Femtodream Tutorial 0
/// \author Luca Barioglio, Anton Riedel

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 0
// Example task illustrating how to create histograms and fill them with basic
// information. A basic event selection is applied.

struct GeorgiosTask{
  HistogramRegistry histos{
    "Histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  // Equivalent of the AliRoot task UserCreateOutputObjects
  void init(o2::framework::InitContext&)
  {
    // Define your axes
    // Constant bin width axis
    AxisSpec vtxZAxis = {100, -20, 20};
    // Variable bin width axis
    std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
                                     1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0,
                                     2.2, 2.4, 2.8, 3.2, 3.6, 4.};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    // Add histograms to histogram manager (as in the output object of in AliPhysics)
    histos.add("hP", ";#it{p} (GeV/#it{c})", kTH1F, {{35, 0.5, 4.}});
    histos.add("hPt", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
  }

  // Equivalent of the AliRoot task UserExec
  void process(aod::Track const& inputTrack)
  {

    histos.fill(HIST("hP"), inputTrack.p());
    histos.fill(HIST("hPt"), inputTrack.pt());
    
  }
    /*
    // Loop over tracks
    for (auto track : inputTracks) {
      //if (fabs(track.eta()) > 0.8) {
      // continue;
      //}
      histos.fill(HIST("hP"), track.p());
      histos.fill(HIST("hPt"), track.pt());
    }
    */

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Equivalent to the AddTask in AliPhysics
  WorkflowSpec workflow{adaptAnalysisTask<GeorgiosTask>(cfgc)};
  return workflow;
}
