/*******************************************************************************

  This file is part of IQmol, a free molecular visualization program. See
  <http://iqmol.org> for more details.

  IQmol is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

********************************************************************************/

#include "XtbOutput.h"
#include "Vibrations.h"
#include <QRegularExpression>
#include <QStringList>
#include <cmath>


namespace IQmol {
namespace Calc {

bool parseEnergy(QString const& output, double& energy)
{
   // "| TOTAL ENERGY               -5.070544386304 Eh   |"
   QRegularExpression re("TOTAL ENERGY\\s+(-?\\d+\\.\\d+)\\s+Eh", 
      QRegularExpression::CaseInsensitiveOption);
   QRegularExpressionMatchIterator iter(re.globalMatch(output));
   bool found(false);
   while (iter.hasNext()) {   // use the last one
      energy = iter.next().captured(1).toDouble();
      found = true;
   }
   return found;
}


// Returns the last match of the first capture group, or false if there is none
static bool lastNumber(QString const& output, QString const& pattern, double& value)
{
   QRegularExpression re(pattern, QRegularExpression::MultilineOption | 
      QRegularExpression::CaseInsensitiveOption);
   QRegularExpressionMatchIterator iter(re.globalMatch(output));
   bool found(false);
   while (iter.hasNext()) {
      value = iter.next().captured(1).toDouble();
      found = true;
   }
   return found;
}


// Fills in the optional details from the xtb output
void parseDetails(QString const& output, Result& result)
{
   double x;
   if (lastNumber(output, "GRADIENT NORM\\s+(-?\\d+\\.\\d+)", x)) result.gradientNorm = x;
   if (lastNumber(output, "HOMO-LUMO GAP\\s+(-?\\d+\\.\\d+)\\s+eV", x)) result.gapEV = x;
   if (lastNumber(output, "CONVERGED AFTER\\s+(\\d+)\\s+ITERATIONS", x)) result.cycles = int(x);

   // The dipole line is the only 'full:' line with exactly four numbers, the last being the total
   if (lastNumber(output, "^\\s*full:\\s+-?\\d+\\.\\d+\\s+-?\\d+\\.\\d+\\s+-?\\d+\\.\\d+\\s+(\\d+\\.\\d+)\\s*$", x)) {
      result.dipoleDebye = x;
   }

   if (lastNumber(output, "-> Gsolv\\s+(-?\\d+\\.\\d+) Eh", x)) {
      result.hasGsolv = true;
      result.gsolv = x;
   }

   // The first summary block is the energy of the starting geometry
   QRegularExpression re(":: total energy\\s+(-?\\d+\\.\\d+) Eh");
   QRegularExpressionMatch match(re.match(output));
   if (match.hasMatch()) result.initialEnergy = match.captured(1).toDouble();
}

bool parseZeroPointEnergy(QString const& output, double& energy)
{
   return lastNumber(output, ":: zero point energy\\s+(-?\\d+\\.\\d+)\\s+Eh", energy);
}


bool parseG98(QString const& text, QVector<int> const& atomicNumbers, VibrationalData& data)
{
   data = VibrationalData();
   int const natoms(atomicNumbers.size());
   QStringList lines(text.split('\n'));
   QRegularExpression separator("\\s+");
   int const n(natoms);

   for (int i = 0; i < lines.size(); ++i) {
       if (!lines[i].trimmed().startsWith("Frequencies --")) continue;

       auto numbers = [&separator](QString const& line) {
          QVector<double> values;
          QString rest(line.mid(line.indexOf("--") + 2));
          for (QString const& token : rest.split(separator, Qt::SkipEmptyParts)) {
              bool ok(false);
              double x(token.toDouble(&ok));
              if (ok) values.append(x);
          }
          return values;
       };

       QVector<double> frequencies(numbers(lines[i]));
       QVector<double> masses, intensities;
       int k(frequencies.size());
       if (k == 0 || k > 3) return false;

       // Red. masses, Frc consts, IR Inten, Raman Activ and Depolar follow, then the atoms
       int row(i + 1);
       for ( ; row < lines.size() && !lines[row].trimmed().startsWith("Atom AN"); ++row) {
           QString line(lines[row].trimmed());
           if (line.startsWith("Red. masses --")) masses = numbers(line);
           if (line.startsWith("IR Inten    --")) intensities = numbers(line);
       }
       if (row + n >= lines.size()) return false;

       QVector<QVector<double> > block(k, QVector<double>(3 * n, 0.0));
       for (int atom = 0; atom < n; ++atom) {
           QStringList tokens(lines[row + 1 + atom].split(separator, Qt::SkipEmptyParts));
           if (tokens.size() < 2 + 3 * k) return false;
           for (int mode = 0; mode < k; ++mode) {
               for (int c = 0; c < 3; ++c) {
                   bool ok(false);
                   double x(tokens[2 + 3 * mode + c].toDouble(&ok));
                   if (!ok) return false;
                   block[mode][3 * atom + c] = x;
               }
           }
       }

       for (int mode = 0; mode < k; ++mode) {
           // mass weighted -> Cartesian displacements
           for (int atom = 0; atom < n; ++atom) {
               double mass(atomicMass(atomicNumbers[atom]));
               if (mass <= 0.0) return false;
               for (int c = 0; c < 3; ++c) block[mode][3 * atom + c] /= std::sqrt(mass);
           }
           double norm(0.0);
           for (double x : block[mode]) norm += x * x;
           norm = std::sqrt(norm);
           if (norm > 0.0) for (double& x : block[mode]) x /= norm;

           data.frequencies.append(frequencies[mode]);
           data.reducedMasses.append(mode < masses.size() ? masses[mode] : 0.0);
           if (mode < intensities.size()) data.intensities.append(intensities[mode]);
           data.modes.append(block[mode]);
       }
       i = row + n;
   }

   if (data.intensities.size() != data.frequencies.size()) data.intensities.clear();
   return !data.frequencies.isEmpty();
}

} } // end namespace IQmol::Calc
