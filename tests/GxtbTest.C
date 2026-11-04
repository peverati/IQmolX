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

#include "Util/GxtbSetup.h"
#include "Util/XtbProcess.h"
#include <QCoreApplication>
#include <QCryptographicHash>
#include <QDir>
#include <QFile>
#include <QProcess>
#include <QStandardPaths>
#include <QTemporaryDir>
#include <cstdio>

// Tests of the g-xTB backend: the download (a local archive stands in for the repository), the 
// command line and the thread limit of the Hessian (a fake xtb records them), and, if the xtb
// program of a g-xTB build is given as the argument, a real calculation of water.

using namespace IQmol::Calc;

static int s_failures(0);

static void check(bool ok, QString const& what)
{
   std::printf("%s: %s\n", ok ? "PASS" : "FAIL", qPrintable(what));
   if (!ok) ++s_failures;
}

static QString readFile(QString const& path)
{
   QFile file(path);
   return file.open(QIODevice::ReadOnly) ? QString::fromLocal8Bit(file.readAll()) : QString();
}

static void writeFile(QString const& path, QString const& text, bool executable = false)
{
   QDir().mkpath(QFileInfo(path).absolutePath());
   QFile file(path);
   file.open(QIODevice::WriteOnly);
   file.write(text.toLocal8Bit());
   file.close();
   if (executable) file.setPermissions(file.permissions() | QFileDevice::ExeOwner);
}

static QString sha256(QString const& path)
{
   QFile file(path);
   file.open(QIODevice::ReadOnly);
   QCryptographicHash hash(QCryptographicHash::Sha256);
   hash.addData(&file);
   return QString::fromLatin1(hash.result().toHex());
}


int main(int argc, char** argv)
{
   QCoreApplication app(argc, argv);
   QStandardPaths::setTestModeEnabled(true);          // never the application folder of the user
   QTemporaryDir dir;
   QString const root(dir.path());

   // --- the download ---
   QString const folder(GxtbSetup::defaultFolder());
   check(folder.contains(".qttest") && folder.endsWith("/g-xtb"), "the g-xtb folder is a test one: " + folder);
   QDir(folder).removeRecursively();
   check(GxtbSetup::defaultExecutable().isEmpty(), "nothing is downloaded at the start");
   check(GxtbSetup::findExecutable().isEmpty(), "no hint and no download: no program");
   check(GxtbSetup::findExecutable("/nonexistent/xtb").isEmpty(), "a wrong path is not replaced silently");
   check(!GxtbSetup::release().isEmpty(), "the release is named: " + GxtbSetup::release());
#ifdef Q_OS_MACOS
   check(GxtbSetup::archiveSha256().size() == 64 && GxtbSetup::archiveUrl().startsWith("https://github.com/grimme-lab/g-xtb/raw/"), 
      "the macOS download is pinned: " + GxtbSetup::archiveUrl());
#endif

   // an archive made like the real one: a folder with bin/xtb
   writeFile(root + "/src/xtb-fake/bin/xtb", "#!/bin/sh\necho fake\n", true);
   QProcess tar;
   tar.setWorkingDirectory(root + "/src");
   tar.start("tar", QStringList() << "-czf" << root + "/fake.tar.gz" << "xtb-fake");
   tar.waitForFinished(10000);
   QString const archive(root + "/fake.tar.gz");
   check(QFile::exists(archive), "the test archive was made");
   QString const url("file://" + archive);

   {
      GxtbSetup setup;
      QString log, error;
      setup.setProgressCallback([&log](QString const& text) { log += text; });
      setup.setSource(url, "0000000000000000000000000000000000000000000000000000000000000000", "fake.tar.gz");
      check(!setup.download(error) && error.contains("not the expected one"), "a wrong sum is refused: " + error.left(60));
      check(GxtbSetup::defaultExecutable().isEmpty() && !QDir(folder + ".download").exists(), 
         "nothing is installed and nothing is left behind after a wrong sum");
   }
   {
      GxtbSetup setup;
      QString log, error;
      setup.setProgressCallback([&log](QString const& text) { log += text; });
      setup.setSource(url, sha256(archive), "fake.tar.gz");
      check(setup.download(error), "a download with the right sum succeeds " + error);
      QString exe(GxtbSetup::defaultExecutable());
      check(exe.endsWith("/g-xtb/xtb-fake/bin/xtb"), "the program is found: " + exe);
      check(GxtbSetup::findExecutable() == exe, "findExecutable returns the download without a hint");
      check(GxtbSetup::findExecutable(exe) == exe, "findExecutable returns a good hint");
      check(log.contains("The SHA-256 sum is the expected one"), "the check of the sum is reported");
      check(!QDir(folder + ".download").exists(), "the work folder is removed");
   }
   {
      GxtbSetup setup;
      QString error;
      setup.setSource("file://" + root + "/nothere.tar.gz", "x", "x.tar.gz");
      check(!setup.download(error) && error.contains("download failed"), "a missing file is reported: " + error.left(60));
      check(!GxtbSetup::defaultExecutable().isEmpty(), "a failed download keeps the earlier copy");
   }

   // --- the command line, with a fake xtb that records it ---
   {
      QString fake(root + "/fakextb");
      writeFile(fake,
         "#!/bin/sh\n"
         "echo \"$@\" >> " + root + "/args.log\n"
         "echo \"threads=$OMP_NUM_THREADS\" >> " + root + "/args.log\n"
         "echo ' :: total energy             -1.000000000000 Eh    ::'\n"
         "echo 'normal termination of xtb'\n", true);
      Input water;
      water.atomicNumbers << 8 << 1 << 1;
      water.coordinates << 0.0 << 0.0 << 0.0 << 0.0 << 0.9 << 0.6 << 0.0 << -0.9 << 0.6;
      water.solvent = "water";
      water.electronicTemperature = 500.0;

      XtbProcess gxtb(fake, XtbProcess::GXTB);
      check(gxtb.name() == "g-xTB", "the method is called g-xTB");
      gxtb.energy(water);
      QString args(readFile(root + "/args.log"));
      check(args.contains("--gxtb") && !args.contains("--gfn"), "the command line has --gxtb and no --gfn: " + args.split('\n').value(0));
      check(args.contains("--gbe water") && !args.contains("--alpb"), "the solvent is a --gbe one, never --alpb");
      check(!args.contains("--etemp"), "no electronic temperature");
      check(args.contains("--chrg 0") && !args.contains("--uhf"), "the charge is passed, but no --uhf for a singlet (it would force an unrestricted wavefunction)");
      Input doublet(water);
      doublet.charge = 1;
      doublet.multiplicity = 2;
      QFile::remove(root + "/args.log");
      gxtb.energy(doublet);
      check(readFile(root + "/args.log").contains("--chrg 1 --uhf 1"), "a doublet gets --uhf 1");

      QFile::remove(root + "/args.log");
      gxtb.setSolvationModel("cosmo");
      gxtb.optimize(water);
      args = readFile(root + "/args.log");
      check(args.contains("--cosmo water") && args.contains("--opt"), "COSMO and the optimization: " + args.split('\n').value(0));
      check(!args.contains("threads=1"), "the optimization is not limited to one thread");

      QFile::remove(root + "/args.log");
      gxtb.frequencies(water);
      args = readFile(root + "/args.log");
      check(args.contains("--hess") && args.contains("threads=1"), "the Hessian runs on one thread (macOS bug of g-xTB)");

      // the plain methods do not use --gxtb and keep their threads
      QFile::remove(root + "/args.log");
      XtbProcess plain(fake, XtbProcess::GFN2);
      plain.frequencies(water);
      args = readFile(root + "/args.log");
      check(!args.contains("--gxtb") && args.contains("--alpb water") && !args.contains("threads=1"), 
         "GFN2 is unchanged: no --gxtb, --alpb, any threads");

      // the xtb of the PATH is never used for g-xTB
      qputenv("XTB_PATH", fake.toLocal8Bit());
      XtbProcess lost("", XtbProcess::GXTB);
      QDir(folder).removeRecursively();
      Result result(lost.energy(water));
      check(!result.success && result.message.contains("g-xTB program"), "without a g-xTB program there is an error: " + result.message.left(50));
      qunsetenv("XTB_PATH");
   }

   // --- a real g-xTB program, if one is given ---
   if (argc > 1 && QFile::exists(argv[1])) {
      Input water;
      water.atomicNumbers << 8 << 1 << 1;
      water.coordinates << 0.0 << 0.0 << 0.2216 << 0.0 << 1.4308 << -0.8864 << 0.0 << -1.4308 << -0.8864;   // Bohr-ish guess in Angstrom is fine: optimized below
      for (double& c : water.coordinates) c *= 0.529177;
      XtbProcess xtb(argv[1], XtbProcess::GXTB);
      Result sp(xtb.energy(water));
      check(sp.success, "real g-xTB: single point runs " + sp.message);
      check(sp.energy < -76.0 && sp.energy > -77.0, QString("real g-xTB: the energy is sensible (%1 Eh)").arg(sp.energy));
      Result opt(xtb.optimize(water));
      check(opt.success && opt.converged, "real g-xTB: the optimization converges");
      water.coordinates = opt.coordinates;
      check(opt.energy < sp.energy, "real g-xTB: optimization lowers the energy");
      water.orbitals = false;
      Result freq(xtb.frequencies(water));
      check(freq.success && freq.hasVibrations && freq.vibrations.frequencies.size() == 3, "real g-xTB: three vibrational modes");
      check(freq.hasVibrations && freq.vibrations.intensities.size() == 3, "real g-xTB: and their IR intensities");
      if (freq.hasVibrations && freq.vibrations.frequencies.size() == 3) {
         QVector<double> f(freq.vibrations.frequencies);
         check(qAbs(f[0] - 1603) < 8 && qAbs(f[1] - 3615) < 10 && qAbs(f[2] - 3672) < 10,
            QString("real g-xTB: the frequencies are those of water (%1, %2, %3)").arg(f[0]).arg(f[1]).arg(f[2]));
      }
      Input withOrbitals(water);
      withOrbitals.orbitals = true;
      Result orbitals(xtb.energy(withOrbitals));
      check(orbitals.success && orbitals.hasOrbitals, "real g-xTB: the orbitals are read from the Molden file");
      if (orbitals.hasOrbitals) {
         std::printf("   %d orbitals, %d basis functions\n", int(orbitals.orbitals.energies.size()), orbitals.orbitals.nBasis());
         check(orbitals.orbitals.energies.size() == orbitals.orbitals.nBasis(), "real g-xTB: as many orbitals as basis functions");
      }
      Input ion(water);
      ion.charge = 1;
      ion.multiplicity = 2;
      check(xtb.energy(ion).success, "real g-xTB: the cation doublet runs");
      Input bad(water);
      bad.multiplicity = 2;
      check(!xtb.energy(bad).success, "real g-xTB: an impossible spin is reported");
      Input solvated(water);
      solvated.solvent = "water";
      Result gas(xtb.energy(water)), gbe(xtb.energy(solvated));
      check(gbe.success && gbe.energy < gas.energy, "real g-xTB: water is stabilised by GBE");
   }

   // --- the real download from GitHub (about 41 MB), only when asked for: GXTB_REAL_DOWNLOAD=1 ---
   if (qEnvironmentVariableIntValue("GXTB_REAL_DOWNLOAD") == 1) {
      QDir(folder).removeRecursively();
      GxtbSetup setup;
      QString log, error;
      setup.setProgressCallback([&log](QString const& text) { log += text; });
      bool ok(setup.download(error));
      check(ok, "real download: the pinned release is fetched, checked and unpacked " + error);
      check(log.contains("The SHA-256 sum is the expected one"), "real download: the sum was checked");
      QString exe(GxtbSetup::defaultExecutable());
      check(!exe.isEmpty(), "real download: the program is found: " + exe);
      if (!exe.isEmpty()) {
         Input water;
         water.atomicNumbers << 8 << 1 << 1;
         water.coordinates << 0.0 << 0.0 << 0.2216 << 0.0 << 0.757 << -0.4692 << 0.0 << -0.757 << -0.4692;
         XtbProcess xtb(QString(), XtbProcess::GXTB);          // no path: the downloaded program
         Result r(xtb.energy(water));
         check(r.success && r.energy < -76.0 && r.energy > -77.0, QString("real download: it calculates water (%1 Eh)").arg(r.energy));
      }
   }

   QDir(folder).removeRecursively();
   std::printf("%d failure(s)\n", s_failures);
   return s_failures ? 1 : 0;
}
