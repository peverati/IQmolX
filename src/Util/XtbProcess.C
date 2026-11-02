/*******************************************************************************

  This file is part of IQmol, a free molecular visualization program. See
  <http://iqmol.org> for more details.

  IQmol is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

********************************************************************************/

#include "XtbProcess.h"
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QProcess>
#include <QProcessEnvironment>
#include <QRegularExpression>
#include <QStandardPaths>
#include <QTemporaryDir>
#include <QTextStream>


namespace IQmol {
namespace Calc {

namespace {

// Element symbols for atomic numbers 0-118
char const* s_symbols[] = {
   "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg",
   "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn",
   "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
   "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb",
   "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
   "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir",
   "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
   "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
   "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv",
   "Ts", "Og"
};

int const s_numSymbols(sizeof(s_symbols)/sizeof(s_symbols[0]));


bool writeXyz(QString const& path, Input const& input)
{
   QFile file(path);
   if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return false;

   QTextStream out(&file);
   out << input.atomicNumbers.size() << "\n\n";
   for (int i = 0; i < input.atomicNumbers.size(); ++i) {
       int Z(input.atomicNumbers[i]);
       out << s_symbols[(Z > 0 && Z < s_numSymbols) ? Z : 0] << " "
           << QString::number(input.coordinates[3*i+0], 'f', 10) << " "
           << QString::number(input.coordinates[3*i+1], 'f', 10) << " "
           << QString::number(input.coordinates[3*i+2], 'f', 10) << "\n";
   }
   return true;
}


// Reads the Cartesian coordinates from an xyz file.  Recent versions of xtb write
// extended xyz with a long comment line, but the atom lines are standard.
bool readXyz(QString const& path, int natoms, QVector<double>& coordinates)
{
   QFile file(path);
   if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return false;

   QTextStream in(&file);
   bool ok(false);
   int n(in.readLine().trimmed().toInt(&ok));
   if (!ok || n != natoms) return false;
   in.readLine();  // comment

   coordinates.clear();
   for (int i = 0; i < natoms; ++i) {
       QStringList tokens(in.readLine().split(QRegularExpression("\\s+"),
          Qt::SkipEmptyParts));
       if (tokens.size() < 4) return false;
       for (int j = 1; j < 4; ++j) {
           double x(tokens[j].toDouble(&ok));
           if (!ok) return false;
           coordinates.append(x);
       }
   }
   return true;
}


bool parseEnergy(QString const& output, double& energy)
{
   // "| TOTAL ENERGY               -5.070544386304 Eh   |"
   QRegularExpression re("TOTAL ENERGY\\s+(-?\\d+\\.\\d+)\\s+Eh");
   QRegularExpressionMatchIterator iter(re.globalMatch(output));
   bool found(false);
   while (iter.hasNext()) {   // use the last one
      energy = iter.next().captured(1).toDouble();
      found = true;
   }
   return found;
}


QString tail(QString const& text, int lines = 15)
{
   QStringList list(text.split('\n'));
   while (!list.isEmpty() && list.last().trimmed().isEmpty()) list.removeLast();
   return list.mid(qMax(0, list.size()-lines)).join('\n');
}

} // end anonymous namespace


XtbProcess::XtbProcess(QString const& executable, Method method)
 : m_executable(executable), m_method(method), m_timeout(600)
{
}


QString XtbProcess::name() const
{
   switch (m_method) {
      case GFN2:  return "GFN2-xTB";
      case GFN1:  return "GFN1-xTB";
      case GFNFF: return "GFN-FF";
   }
   return "xTB";
}


QStringList XtbProcess::methodArguments() const
{
   switch (m_method) {
      case GFN2:  return QStringList() << "--gfn" << "2";
      case GFN1:  return QStringList() << "--gfn" << "1";
      case GFNFF: return QStringList() << "--gfnff";
   }
   return QStringList();
}


QString XtbProcess::findExecutable(QString const& hint)
{
   QStringList candidates;
   if (!hint.isEmpty()) candidates << hint;
   QString env(qEnvironmentVariable("XTB_PATH"));
   if (!env.isEmpty()) candidates << env;

   for (QString const& path : candidates) {
       QFileInfo info(path);
       if (info.isDir()) info.setFile(QDir(path).filePath("xtb"));
       if (info.isFile() && info.isExecutable()) return info.absoluteFilePath();
   }

   // A GUI application on macOS may not inherit the shell PATH, so also
   // look in the usual locations.
   QString path(QStandardPaths::findExecutable("xtb"));
   if (path.isEmpty()) {
      path = QStandardPaths::findExecutable("xtb",
         QStringList() << "/opt/homebrew/bin" << "/usr/local/bin");
   }
   return path;
}


Result XtbProcess::energy(Input const& input)
{
   return run(input, false);
}


Result XtbProcess::optimize(Input const& input)
{
   return run(input, true);
}


Result XtbProcess::run(Input const& input, bool optimize)
{
   Result result;

   int natoms(input.atomicNumbers.size());
   if (natoms == 0 || input.coordinates.size() != 3*natoms) {
      result.message = "Invalid input: inconsistent number of atoms and coordinates";
      return result;
   }
   if (input.multiplicity < 1) {
      result.message = "Invalid spin multiplicity";
      return result;
   }

   QString exe(findExecutable(m_executable));
   if (exe.isEmpty()) {
      result.message = "Unable to find the xtb executable.  Set its location in the "
                       "preferences or in the XTB_PATH environment variable.";
      return result;
   }

   QTemporaryDir dir;
   if (!dir.isValid()) {
      result.message = "Failed to create a temporary directory";
      return result;
   }

   QString inputFile(dir.filePath("input.xyz"));
   if (!writeXyz(inputFile, input)) {
      result.message = "Failed to write xtb input file";
      return result;
   }

   QStringList args;
   args << "input.xyz";
   args << methodArguments();
   args << "--chrg" << QString::number(input.charge);
   args << "--uhf"  << QString::number(input.multiplicity - 1);
   args << (optimize ? "--opt" : "--sp");

   QProcess process;
   process.setWorkingDirectory(dir.path());
   process.setProcessChannelMode(QProcess::MergedChannels);
   QProcessEnvironment env(QProcessEnvironment::systemEnvironment());
   env.insert("OMP_STACKSIZE", "1G");
   process.setProcessEnvironment(env);
   process.start(exe, args);

   if (!process.waitForStarted()) {
      result.message = "Failed to start " + exe;
      return result;
   }

   if (!process.waitForFinished(m_timeout * 1000)) {
      process.kill();
      process.waitForFinished();
      result.message = "xtb timed out after " + QString::number(m_timeout) + " seconds";
      return result;
   }

   QString output(QString::fromLocal8Bit(process.readAll()));

   if (process.exitStatus() != QProcess::NormalExit) {
      result.message = "xtb crashed\n" + tail(output);
      return result;
   }

   if (!parseEnergy(output, result.energy)) {
      result.message = "xtb failed (exit code " + QString::number(process.exitCode())
                     + ")\n" + tail(output);
      return result;
   }

   if (optimize) {
      // xtb writes the last geometry even when the optimization does not converge
      if (!readXyz(dir.filePath("xtbopt.xyz"), natoms, result.coordinates)) {
         result.message = "Failed to read the optimized geometry from xtb\n" + tail(output);
         return result;
      }
      result.converged = output.contains("GEOMETRY OPTIMIZATION CONVERGED");
      if (!result.converged) {
         result.message = "Geometry optimization did not converge\n" + tail(output);
      }
   }else {
      result.coordinates = input.coordinates;
      result.converged = true;
   }

   result.success = true;
   return result;
}

} } // end namespace IQmol::Calc
