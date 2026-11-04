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

#include "PythonSetup.h"
#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QProcess>
#include <QProcessEnvironment>
#include <QStandardPaths>


namespace IQmol {
namespace Calc {

PythonSetup::PythonSetup() : m_interactive(false), m_cancelled(false)
{
}


static QString executableIn(QString const& path)
{
   QFileInfo info(path);
   return (info.isFile() && info.isExecutable()) ? info.absoluteFilePath() : QString();
}


QString PythonSetup::findUv(QString const& hint)
{
   QStringList candidates;
   if (!hint.isEmpty()) candidates << hint;
   QString fromEnvironment(qEnvironmentVariable("IQMOLX_UV"));
   if (!fromEnvironment.isEmpty()) candidates << fromEnvironment;

   // IQmol.app/Contents/MacOS/IQmol, with uv in Contents/Resources/bin
   QString appDir(QCoreApplication::applicationDirPath());
   if (!appDir.isEmpty()) candidates << appDir + "/../Resources/bin/uv";

   for (QString const& dir : qEnvironmentVariable("PATH").split(QDir::listSeparator(), Qt::SkipEmptyParts)) {
       candidates << dir + "/uv";
   }
   candidates << QDir::homePath() + "/.local/bin/uv" << "/opt/homebrew/bin/uv" << "/usr/local/bin/uv";

   for (QString const& path : candidates) {
       QString found(executableIn(path));
       if (!found.isEmpty()) return found;
   }
   return QString();
}


QString PythonSetup::defaultEnvironment()
{
   QString base(QStandardPaths::writableLocation(QStandardPaths::AppDataLocation));
   if (base.isEmpty()) base = QDir::homePath() + "/.iqmolx";
   return base + "/python";
}


QString PythonSetup::pythonIn(QString const& environment)
{
   return executableIn(environment + "/bin/python");
}


QString PythonSetup::wheelIn(QString const& directory)
{
   QStringList wheels(QDir(directory).entryList(QStringList("iqmolx_engines-*.whl"), 
      QDir::Files, QDir::Name | QDir::Reversed));
   return wheels.isEmpty() ? QString() : QDir::cleanPath(QDir(directory).absoluteFilePath(wheels.first()));
}


QString PythonSetup::sourceVersion(QString const& source)
{
   QFileInfo info(source);
   if (info.isFile()) {                       // iqmolx_engines-0.1.0-py3-none-any.whl
      QStringList parts(info.fileName().split('-'));
      return parts.size() > 1 ? parts[1] : QString();
   }
   QFile file(source + "/pyproject.toml");
   if (!file.open(QIODevice::ReadOnly)) return QString();
   for (QString const& line : QString::fromUtf8(file.readAll()).split('\n')) {
       QString text(line.trimmed());
       if (text.startsWith("version") && text.contains('"')) return text.section('"', 1, 1);
   }
   return QString();
}


QString PythonSetup::helperSource()
{
   // A source folder given by the developer comes first
   QString fromEnvironment(qEnvironmentVariable("IQMOLX_ENGINES_PATH"));
   if (!fromEnvironment.isEmpty()) {
      QFileInfo info(fromEnvironment + "/../pyproject.toml");
      if (info.isFile()) return QDir::cleanPath(info.absolutePath());
   }

   // IQmol.app/Contents/MacOS/IQmol, with the helper in Contents/Resources/python
   QString appDir(QCoreApplication::applicationDirPath());
   if (!appDir.isEmpty()) {
      QString resources(appDir + "/../Resources/python");
      QString wheel(wheelIn(resources));
      if (!wheel.isEmpty()) return wheel;
      QFileInfo info(resources + "/iqmolx-engines/pyproject.toml");
      if (info.isFile()) return QDir::cleanPath(info.absolutePath());
   }
   return QString();
}


/// Runs uv with the arguments, showing its output.  The command is shown first.
bool PythonSetup::run(QStringList const& arguments, QString& error)
{
   if (m_uv.isEmpty()) m_uv = findUv();
   if (m_uv.isEmpty()) {
      error = "uv was not found.  It is the tool that makes the Python environment; install "
              "it (for example with 'brew install uv') or give its location in $IQMOLX_UV.";
      return false;
   }

   auto report = [this](QString const& text) { if (m_progress) m_progress(text); };
   report("$ uv " + arguments.join(' ') + "\n");

   QProcess process;
   process.setProcessChannelMode(QProcess::MergedChannels);
   QProcessEnvironment env(QProcessEnvironment::systemEnvironment());
   env.insert("UV_NO_PROGRESS", "1");           // no animated bars in the text window
   env.insert("NO_COLOR", "1");
   process.setProcessEnvironment(env);
   process.start(m_uv, arguments);
   if (!process.waitForStarted(10000)) {
      error = "Unable to start uv: " + process.errorString();
      return false;
   }

   QString tail;
   while (process.state() != QProcess::NotRunning) {
      if (m_cancelled) {
         process.kill();
         process.waitForFinished(3000);
         error = "Cancelled";
         return false;
      }
      process.waitForReadyRead(50);
      QString text(QString::fromLocal8Bit(process.readAll()));
      if (!text.isEmpty()) {
         report(text);
         tail = (tail + text).right(1500);
      }
      if (m_interactive) QCoreApplication::processEvents();
   }
   QString rest(QString::fromLocal8Bit(process.readAll()));
   if (!rest.isEmpty()) {
      report(rest);
      tail = (tail + rest).right(1500);
   }

   if (process.exitStatus() != QProcess::NormalExit || process.exitCode() != 0) {
      error = "uv failed (exit code " + QString::number(process.exitCode()) + ")\n" + tail.trimmed();
      return false;
   }
   return true;
}


bool PythonSetup::createEnvironment(QString const& environment, QString& error)
{
   m_cancelled = false;
   QString source(helperSource());
   if (source.isEmpty()) {
      error = "The iqmolx-engines package was not found (it is a wheel in the Resources folder "
              "of the application, or set $IQMOLX_ENGINES_PATH to its src folder).";
      return false;
   }

   QDir().mkpath(QFileInfo(environment).absolutePath());
   QStringList venv;
   venv << "venv" << "--python" << "3.12" << "--python-preference" << "only-managed" 
        << "--allow-existing" << environment;
   if (!run(venv, error)) return false;

   QString python(pythonIn(environment));
   if (python.isEmpty()) {
      error = "uv made no Python in " + environment;
      return false;
   }

   QStringList pip;
   pip << "pip" << "install" << "--python" << python << source;
   if (!run(pip, error)) return false;

   if (m_progress) m_progress("\nThe Python environment is ready: " + environment + "\n");
   return true;
}


bool PythonSetup::updateHelper(QString const& python, QString& error)
{
   m_cancelled = false;
   QString source(helperSource());
   if (source.isEmpty()) {
      error = "The iqmolx-engines package of this application was not found.";
      return false;
   }
   QStringList pip;
   pip << "pip" << "install" << "--python" << python << "--reinstall-package" << "iqmolx-engines" << source;
   if (!run(pip, error)) return false;
   if (m_progress) m_progress("\nThe helper is up to date.\n");
   return true;
}


bool PythonSetup::install(QString const& python, QStringList const& requirements, QString& error)
{
   m_cancelled = false;
   if (requirements.isEmpty()) return true;
   QStringList pip;
   pip << "pip" << "install" << "--python" << python << requirements;
   if (!run(pip, error)) return false;
   if (m_progress) m_progress("\nInstalled: " + requirements.join(", ") + "\n");
   return true;
}

} } // end namespace IQmol::Calc
