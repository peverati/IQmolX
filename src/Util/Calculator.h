#pragma once
/*******************************************************************************

  This file is part of IQmol, a free molecular visualization program. See
  <http://iqmol.org> for more details.

  IQmol is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

********************************************************************************/

#include <QString>
#include <QVector>

/*!
 *  \file Calculator.h
 *
 *  \brief Minimal interface for an energy/geometry-optimization backend that
 *  works on plain arrays (atomic numbers + Cartesian coordinates), so that the
 *  molecule classes do not need to know how the calculation is performed.
 *  Implementations may run an external program or call a library.
 */

namespace IQmol {
namespace Calc {

   struct Input {
      QVector<int>    atomicNumbers;
      QVector<double> coordinates;   // x0,y0,z0,x1,y1,z1,... in Angstrom
      int             charge;
      int             multiplicity;

      Input() : charge(0), multiplicity(1) { }
   };

   struct Result {
      bool            success;       // the calculation ran and produced a result
      bool            converged;     // geometry optimization converged
      double          energy;        // Hartree
      QVector<double> coordinates;   // final coordinates, Angstrom (same layout as Input)
      QString         message;       // error description or program log tail
      QString         log;           // full program output and input (for debugging)

      // Optional details; negative / zero values mean the program did not report them
      bool            cancelled;
      double          initialEnergy; // Hartree, energy of the starting geometry
      double          gradientNorm;  // Eh/a0
      double          gapEV;         // HOMO-LUMO gap
      double          dipoleDebye;
      int             cycles;        // optimization cycles
      double          seconds;       // wall time

      Result() : success(false), converged(false), energy(0.0), cancelled(false),
         initialEnergy(0.0), gradientNorm(-1.0), gapEV(-1.0), dipoleDebye(-1.0), 
         cycles(0), seconds(0.0) { }
   };

   class Calculator {
      public:
         virtual ~Calculator() { }

         virtual QString name() const = 0;

         /// Single-point energy
         virtual Result energy(Input const&) = 0;

         /// Geometry optimization
         virtual Result optimize(Input const&) = 0;
   };

} } // end namespace IQmol::Calc
