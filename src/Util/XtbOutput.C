/*******************************************************************************

  This file is part of IQmol, a free molecular visualization program. See
  <http://iqmol.org> for more details.

  IQmol is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

********************************************************************************/

#include "XtbOutput.h"
#include <QRegularExpression>


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

} } // end namespace IQmol::Calc
