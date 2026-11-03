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
#include <functional>

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
      QString         solvent;       // implicit solvent (ALPB name), empty for gas phase

      Input() : charge(0), multiplicity(1) { }
   };

   /// The result of a vibrational analysis
   struct VibrationalData {
      QVector<double>          frequencies;     // cm^-1, imaginary ones negative, ascending
      QVector<double>          intensities;     // IR, km/mol (empty if not available)
      QVector<double>          reducedMasses;   // amu
      QVector<QVector<double> > modes;          // per mode the Cartesian displacement of each atom
                                                // (3N numbers), normalized to a sum of squares of 1
      double                   zeroPointEnergy; // Hartree

      VibrationalData() : zeroPointEnergy(0.0) { }

      int imaginary() const {
         int n(0);
         for (double f : frequencies) if (f < 0.0) ++n;
         return n;
      }
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
      bool            hasVibrations; // the vibrational data below is available
      VibrationalData vibrations;
      bool            hasGsolv;      // a solvation free energy was reported
      double          gsolv;         // Hartree

      Result() : success(false), converged(false), energy(0.0), cancelled(false),
         initialEnergy(0.0), gradientNorm(-1.0), gapEV(-1.0), dipoleDebye(-1.0), 
         cycles(0), seconds(0.0), hasVibrations(false), hasGsolv(false), gsolv(0.0) { }
   };

   class Calculator {
      public:
         virtual ~Calculator() { }

         virtual QString name() const = 0;

         /// Single-point energy
         virtual Result energy(Input const&) = 0;

         /// Geometry optimization
         virtual Result optimize(Input const&) = 0;

         /// Vibrational frequencies, normal modes and IR intensities at the geometry of the
         /// input, which should be a stationary point.  The energy and the norm of the 
         /// gradient at that geometry are also reported.
         virtual Result frequencies(Input const&) 
         { 
            Result result;
            result.message = name() + " cannot calculate vibrational frequencies";
            return result; 
         }

         /// Text describing the progress of a calculation is passed to this function as 
         /// it becomes available.  Only used together with setInteractive(true).
         virtual void setProgressCallback(std::function<void(QString const&)>) { }

         /// When true a GUI event loop keeps running during a calculation so that the
         /// window stays responsive and cancel() can be called.  The default is to block.
         virtual void setInteractive(bool) { }

         /// Stops a running calculation (interactive mode only)
         virtual void cancel() { }
   };

} } // end namespace IQmol::Calc
