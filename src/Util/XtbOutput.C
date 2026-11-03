/*******************************************************************************

  Copyright (C) 2026 Roberto Peverati
  Developed with the assistance of Claude (Anthropic).

  This file is part of IQmol, a free molecular visualization program. See
  <http://iqmol.org> for more details.

  IQmol is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  IQmol is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along
  with IQmol.  If not, see <http://www.gnu.org/licenses/>.

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


bool parseMolden(QString const& text, OrbitalData& data)
{
   data = OrbitalData();

   // Normalization of a primitive x^l exp(-a r^2) in the Cartesian s, p and dxy form used by
   // IQmol (the same factors as Data::Shell)
   static double const Ns(std::pow(2.0/M_PI, 0.75));
   static double const prefactor[3] = { Ns, 2.0*Ns, 4.0*Ns };
   static double const power[3]     = { 0.75, 1.25, 1.75 };
   static int const nFunctions[3]   = { 1, 3, 6 };

   enum Section { None, Atoms, Gto, Mo } section(None);
   int atom(-1);
   int primitivesLeft(0);
   bool alpha(true);
   bool haveOrbital(false);
   int nBasis(0);

   QStringList lines(text.split('\n'));
   for (QString const& raw : lines) {
       QString line(raw.trimmed());
       if (line.isEmpty()) continue;

       if (line.startsWith('[')) {
          QString key(line.toLower());
          if (key.startsWith("[atoms]"))   { section = Atoms; }
          else if (key.startsWith("[gto]")) { section = Gto; atom = -1; primitivesLeft = 0; }
          else if (key.startsWith("[mo]"))  { 
             section = Mo; 
             for (OrbitalShell const& shell : data.shells) nBasis += nFunctions[shell.angularMomentum];
          }
          else section = None;
          continue;
       }

       QStringList tokens(line.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts));

       if (section == Gto) {
          bool ok(false);
          if (tokens[0].size() == 1 && tokens[0][0].isLetter()) {
             // shell header: "s 4 1.00"
             int l(QString("spd").indexOf(tokens[0].toLower()));
             if (l < 0 || tokens.size() < 2 || atom < 0) return false;
             primitivesLeft = tokens[1].toInt(&ok);
             if (!ok) return false;
             OrbitalShell shell;
             shell.atom = atom;
             shell.angularMomentum = l;
             data.shells.append(shell);
          }else if (primitivesLeft > 0) {
             if (tokens.size() < 2 || data.shells.isEmpty()) return false;
             double exponent(tokens[0].toDouble(&ok));
             double coefficient(tokens[1].toDouble(&ok));
             if (!ok || exponent <= 0.0) return false;
             OrbitalShell& shell(data.shells.last());
             int l(shell.angularMomentum);
             shell.exponents.append(exponent);
             shell.coefficients.append(coefficient / (prefactor[l] * std::pow(exponent, power[l])));
             --primitivesLeft;
          }else {
             // atom header: "1 0"
             atom = tokens[0].toInt(&ok) - 1;
             if (!ok || atom < 0) return false;
          }
       }else if (section == Mo) {
          QString first(tokens[0].toLower());
          if (first.startsWith("sym="))   { haveOrbital = false; continue; }
          if (first.startsWith("ene=")) {
             if (tokens.size() < 2) return false;
             data.energies.append(tokens[1].toDouble());
             data.occupations.append(0.0);
             data.coefficients.append(QVector<double>());
             haveOrbital = true;
             continue;
          }
          if (first.startsWith("spin=")) {
             alpha = line.contains("alpha", Qt::CaseInsensitive);
             if (!alpha && haveOrbital) {   // not wanted, remove the orbital again
                data.energies.removeLast();
                data.occupations.removeLast();
                data.coefficients.removeLast();
                haveOrbital = false;
             }
             continue;
          }
          if (first.startsWith("occup=")) {
             if (haveOrbital && tokens.size() > 1) data.occupations.last() = tokens[1].toDouble();
             continue;
          }
          if (!haveOrbital || !alpha || tokens.size() < 2) continue;
          bool ok(false);
          double c(tokens[1].toDouble(&ok));
          if (!ok) return false;
          data.coefficients.last().append(c);
       }
   }

   if (data.shells.isEmpty() || data.coefficients.isEmpty()) return false;
   for (QVector<double> const& c : data.coefficients) {
       if (c.size() != nBasis) return false;
   }
   return true;
}

} } // end namespace IQmol::Calc
