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

#include "GxtbSetup.h"

#include <QCoreApplication>
#include <QCryptographicHash>
#include <QDir>
#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QProcess>
#include <QStandardPaths>
#include <QSysInfo>


namespace IQmol {
namespace Calc {

namespace {

// The release of https://github.com/grimme-lab/g-xtb that is fetched: the binaries folder at
// one commit (v2.0.1, 14 May 2026), so that the sums below stay valid when the repository changes.
char const* const s_commit  = "23f81c2db64fa4272c25b86f8776ce862aeb59e2";
char const* const s_release = "g-xTB v2.0.1 (xtb 6.7.1, build of 14 May 2026)";

struct Archive { char const* name; char const* sha256; qint64 bytes; };

Archive archiveForThisSystem()
{
#ifdef Q_OS_MACOS
   if (QSysInfo::currentCpuArchitecture() == "arm64") {
      return { "xtb-6.7.1-gxtb-140526-macos-arm64.tar.gz", 
               "035a37f241a24067eb854177ea0dd3e1d1b48aca2155cf934aae315351e098af", 41367318 };
   }
   if (QSysInfo::currentCpuArchitecture() == "x86_64") {
      return { "xtb-6.7.1-gxtb-140526-macos-x86_64.tar.gz", 
               "3f5f7ee80a66bce911cd51dcca894fa087e7039b51e8c2e036ba74208a8e3b62", 0 };
   }
#endif
   return { "", "", 0 };
}

QString executableName()
{
#ifdef Q_OS_WIN
   return "xtb.exe";
#else
   return "xtb";
#endif
}

} // end anonymous namespace


GxtbSetup::GxtbSetup() : m_interactive(false), m_cancelled(false)
{
   m_url    = archiveUrl();
   m_sha256 = archiveSha256();
   m_name   = archiveName();
}


QString GxtbSetup::defaultFolder()
{
   QString base(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation));
   if (base.isEmpty()) base = QDir::homePath() + "/.iqmolx";
   return base + "/g-xtb";
}


QString GxtbSetup::defaultExecutable()
{
   // g-xtb/<folder of the release>/bin/xtb
   QDir dir(defaultFolder());
   for (QString const& sub : dir.entryList(QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name)) {
       QFileInfo info(dir.filePath(sub + "/bin/" + executableName()));
       if (info.isFile() && info.isExecutable()) return QDir::cleanPath(info.absoluteFilePath());
   }
   return QString();
}


QString GxtbSetup::findExecutable(QString const& hint)
{
   if (!hint.isEmpty()) {
      QFileInfo info(hint);
      if (info.isFile() && info.isExecutable()) return QDir::cleanPath(info.absoluteFilePath());
      return QString();      // a path that was given and does not work is not replaced silently
   }
   return defaultExecutable();
}


QString GxtbSetup::release()        { return s_release; }
QString GxtbSetup::archiveName()    { return archiveForThisSystem().name; }
QString GxtbSetup::archiveSha256()  { return archiveForThisSystem().sha256; }
qint64  GxtbSetup::archiveBytes()   { return archiveForThisSystem().bytes; }

QString GxtbSetup::archiveUrl()
{
   QString name(archiveName());
   if (name.isEmpty()) return QString();
   return QString("https://github.com/grimme-lab/g-xtb/raw/%1/binaries/%2").arg(s_commit, name);
}


void GxtbSetup::setSource(QString const& url, QString const& sha256, QString const& name)
{
   m_url = url;
   m_sha256 = sha256;
   m_name = name;
}


/// Runs a program, passing its output on.  With a file to watch the size of the file is
/// reported as the progress of a download.
bool GxtbSetup::run(QString const& program, QStringList const& arguments, QString& error,
   QString const& watchFile, qint64 expectedBytes)
{
   auto report = [this](QString const& text) { if (m_progress) m_progress(text); };
   report("$ " + program + " " + arguments.join(' ') + "\n");

   QProcess process;
   process.setProcessChannelMode(QProcess::MergedChannels);
   process.start(program, arguments);
   if (!process.waitForStarted(10000)) {
      error = "Unable to start " + program + ": " + process.errorString();
      return false;
   }

   QString tail;
   QElapsedTimer clock;
   clock.start();
   qint64 reported(-1);
   while (process.state() != QProcess::NotRunning) {
      if (m_cancelled) {
         process.kill();
         process.waitForFinished(3000);
         error = "Cancelled";
         return false;
      }
      process.waitForReadyRead(50);
      QString text(QString::fromLocal8Bit(process.readAll()));
      if (!text.isEmpty()) { report(text); tail = (tail + text).right(1500); }
      if (!watchFile.isEmpty() && clock.elapsed() > 700) {
         clock.restart();
         qint64 megabytes(QFileInfo(watchFile).size() / 1000000);
         if (megabytes != reported) {
            reported = megabytes;
            report("   downloaded " + QString::number(megabytes) 
               + (expectedBytes > 0 ? " of " + QString::number(expectedBytes / 1000000) : QString()) + " MB\n");
         }
      }
      if (m_interactive) QCoreApplication::processEvents();
   }
   QString rest(QString::fromLocal8Bit(process.readAll()));
   if (!rest.isEmpty()) { report(rest); tail = (tail + rest).right(1500); }

   if (process.exitStatus() != QProcess::NormalExit || process.exitCode() != 0) {
      error = program + " failed (exit code " + QString::number(process.exitCode()) + ")\n" + tail.trimmed();
      return false;
   }
   return true;
}


bool GxtbSetup::download(QString& error)
{
   m_cancelled = false;
   if (m_name.isEmpty() || m_url.isEmpty()) {
      error = "There is no g-xTB download for this system.  Get a g-xTB build from "
              "https://github.com/grimme-lab/g-xtb and choose its xtb program with Browse.";
      return false;
   }

   QString folder(defaultFolder());
   QString work(folder + ".download");
   QDir(work).removeRecursively();
   if (!QDir().mkpath(work)) {
      error = "Unable to make the folder " + work;
      return false;
   }
   auto fail = [&work](QString const& why, QString& error) { 
      QDir(work).removeRecursively(); 
      error = why; 
      return false; 
   };

   auto report = [this](QString const& text) { if (m_progress) m_progress(text); };
   QString archive(work + "/" + m_name);
   QString protocol(m_url.startsWith("file:") ? "=file" : "=https");
   QStringList curl;
   curl << "--fail" << "--location" << "--silent" << "--show-error" << "--connect-timeout" << "30"
        << "--proto" << protocol << "--proto-redir" << protocol << "--output" << archive << m_url;
   QString problem;
   if (!run("curl", curl, problem, archive, archiveBytes())) {
      return fail(m_cancelled ? QString("Cancelled") : "The download failed: " + problem, error);
   }

   // The sum must be that of the release this IQmol was made for
   QFile file(archive);
   if (!file.open(QIODevice::ReadOnly)) return fail("The downloaded file cannot be read", error);
   QCryptographicHash hash(QCryptographicHash::Sha256);
   if (!hash.addData(&file)) return fail("The downloaded file cannot be read", error);
   file.close();
   QString sum(QString::fromLatin1(hash.result().toHex()));
   if (sum != m_sha256) {
      return fail("The downloaded file is not the expected one and was thrown away (SHA-256 " 
         + sum + ", expected " + m_sha256 + ").", error);
   }
   report("The SHA-256 sum is the expected one.\n");

   QString unpacked(work + "/unpacked");
   QDir().mkpath(unpacked);
   if (!run("tar", QStringList() << "-xzf" << archive << "-C" << unpacked, problem)) {
      return fail(m_cancelled ? QString("Cancelled") : "Unpacking failed: " + problem, error);
   }

   // The program must be where the release says, before the old copy is replaced
   bool found(false);
   QDir dir(unpacked);
   for (QString const& sub : dir.entryList(QDir::Dirs | QDir::NoDotAndDotDot)) {
       QFileInfo info(dir.filePath(sub + "/bin/" + executableName()));
       if (info.isFile()) found = true;
   }
   if (!found) return fail("The archive does not contain bin/" + executableName(), error);

   QDir(folder).removeRecursively();
   if (!QDir().rename(unpacked, folder)) return fail("Unable to move the program to " + folder, error);
   QDir(work).removeRecursively();
   report("\ng-xTB is ready: " + defaultExecutable() + "\n");
   return true;
}

} } // end namespace IQmol::Calc
