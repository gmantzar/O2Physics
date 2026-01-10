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

/// \file femtoDreamCutCulator.cxx
/// \brief Executable that encodes physical selection criteria in a bit-wise
/// selection \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#include <iostream>
#include <random>
#include <string>
#include "PWGCF/FemtoDream/Utils/femtoDreamCutCulator.h"
#include "PWGCF/FemtoDream/Core/femtoDreamSelection.h"
#include "PWGCF/FemtoDream/Core/femtoDreamTrackSelection.h"
#include "PWGCF/DataModel/FemtoDerived.h"

using namespace o2::analysis::femtoDream;

/// The function takes the path to the dpl-config.json as a argument and the
/// does a Q&A session for the user to find the appropriate selection criteria
/// for the analysis task
int main(int /*argc*/, char* argv[])
{
  std::string configFileName(argv[1]);
  std::ifstream configFile(configFileName);

  if (configFile.is_open()) {
    FemtoDreamCutculator cut;
    cut.init(argv[1]);

    std::cout
      << "Do you want to work with tracks or V0s or Cascades or Resonances (T/V/C/R)? >";
    std::string choice;
    std::string flag;
    std::cin >> choice;
    flag = choice;

    bool ProducerIsReso = cut.getProducerIsReso();
    if (choice == std::string("T")) {
      if (ProducerIsReso)
      {
        /// configs in femtoDreamProducerTaskReso.cxx are required to begin with lower case
        cut.setTrackSelectionFromFile("confTrk");
        cut.setPIDSelectionFromFile("confTrk");
      } else {
        cut.setTrackSelectionFromFile("ConfTrk");
        cut.setPIDSelectionFromFile("ConfTrk");
      }
    } else if (choice == std::string("V")) {
      if(ProducerIsReso) 
      {
        std::cout << "Dow you want to select Lambda of K0short (L/K)? >";
        std::cin >> choice;

        if (choice == std::string("L")) {
          std::cout << "Do you want to select V0s or one of its children (V/T)? >";
          std::cin >> choice;

          cut.setV0SelectionFromFile("confLambda");
          cut.setTrackSelectionFromFile("confLambdaChild");
          cut.setPIDSelectionFromFile("confLambdaChild");
          
        } else if (choice == std::string("K")) {

          std::cout << "Do you want to select V0s or one of its children (V/T)? >";
          std::cin >> choice;

          cut.setV0SelectionFromFile("confK0short");
          cut.setTrackSelectionFromFile("confK0shortChild");
          cut.setPIDSelectionFromFile("confK0shortChild");
        }
      } else {

        std::cout << "Do you want to select V0s or one of its children (V/T)? >";
        std::cin >> choice;

        cut.setV0SelectionFromFile("ConfV0");
        cut.setTrackSelectionFromFile("ConfChild");
        cut.setPIDSelectionFromFile("ConfChild");
      }

    } else if (choice == std::string("R")) {

      std::cout  << "Do you want to select Phi or Kstar (P/K)? >";  
      std::cin >> choice;

      if (choice == std::string("P")) {

        std::cout << "Do you want to select Resos or one of its daughters (R/T)? >";
        std::cin >> choice;

        cut.setResoSelectionFromFile("confPhi");
        cut.setTrackSelectionFromFile("confPhiDaughter");
        cut.setPIDSelectionFromFile("confPhiDaughter");

      } else if (choice == std::string("K")) {

        std::cout << "Do you want to select Resos or one of its daughters (R/T)? >";
        std::cin >> choice;

        cut.setResoSelectionFromFile("confKstar");
        cut.setTrackSelectionFromFile("confKstarDaughter");
        cut.setPIDSelectionFromFile("confKstarDaughter");
        
      }
      
    } else if (choice == std::string("C")) {

      if (ProducerIsReso) {
        std::cout << "Do you want to select Xi or Omega (X/O)? >";
        std::cin >> choice;

        if(choice == std::string("X")) {

          std::cout << "Do you want to select cascades, V0-Daughter tracks of the cascades or the Bachelor track (C/V/B)? >";
          std::cin >> choice;
          if(choice == std::string("C")) {
            cut.setCascadeSelectionFromFile("confXi");
            choice = "C";
          } else if (choice == std::string("V")) {
            cut.setTrackSelectionFromFile("confXiV0Child");
            cut.setPIDSelectionFromFile("confXiV0Child");
            choice = "T";
          } else if (choice == std::string("B")) {
            cut.setTrackSelectionFromFile("confXiBachelor");
            cut.setPIDSelectionFromFile("confXiBachelor");
            choice = "T";
          }

        } else if (choice == std::string("O")) {

          std::cout << "Do you want to select cascades, V0-Daughter tracks of the cascades or the Bachelor track (C/V/B)? >";
          std::cin >> choice;
          if(choice == std::string("C")) {
            cut.setCascadeSelectionFromFile("confOmega");
            choice = "C";
          } else if (choice == std::string("V")) {
            cut.setTrackSelectionFromFile("confOmegaV0Child");
            cut.setPIDSelectionFromFile("confOmegaV0Child");
            choice = "T";
          } else if (choice == std::string("B")) {
            cut.setTrackSelectionFromFile("confOmegaBachelor");
            cut.setPIDSelectionFromFile("confOmegaBachelor");
            choice = "T";
          }

        }
      } else {
        std::cout << "Do you want to select cascades, V0-Daughter tracks of the cascades or the Bachelor track (C/V/B)? >";
        std::cin >> choice;
        if (choice == std::string("C")) {
          cut.setCascadeSelectionFromFile("ConfCascade");
          choice = "C";
        } else if (choice == std::string("V")) {
          cut.setTrackSelectionFromFile("ConfCascV0Child");
          cut.setPIDSelectionFromFile("ConfCascV0Child");
          choice = "T";
        } else if (choice == std::string("B")) {
          cut.setTrackSelectionFromFile("ConfCascBachelor");
          cut.setPIDSelectionFromFile("ConfCascBachelor");
          choice = "T";
        } else {
          std::cout << "Option not recognized. Break...";
          return 2;
        }
      }

    } else {
      std::cout << "Option not recognized. Break...";
      return 2;
    }
    /// \todo factor out the pid here
    /// cut.setTrackSelection(femtoDreamTrackSelection::kPIDnSigmaMax,
    /// femtoDreamSelection::kAbsUpperLimit, "ConfTrk");

    std::cout << "Do you want to manually select cuts or create systematic "
                 "variations(M/V)? >";
    std::string manual;
    std::cin >> manual;

    if (manual == std::string("M")) {
      cut.analyseCuts(choice, false, 1.0f, flag);
    } else if (manual == std::string("V")) {

      std::cout << "Do you want to set some selections (y/n)? >";
      std::string set;
      std::cin >> set;
      
      cut.Assign(set,choice); 
      
      std::ofstream out("CutCulator.txt");
      std::streambuf* coutbuf = std::cout.rdbuf(); // save old buf
      std::cout.rdbuf(out.rdbuf());                // redirect std::cout to out.txt!
      for (int i = 0; i < 1000; i++) {
        cut.analyseCuts(choice, true, 1, flag);
      }
      std::cout.rdbuf(coutbuf); // reset to standard output again
    } else {
      std::cout << "Option not recognized. Break...";
      return 2;
    }

  } else {
    std::cout << "The configuration file " << configFileName
              << " could not be found or could not be opened.";
    return 1;
  }

  return 0;
}
