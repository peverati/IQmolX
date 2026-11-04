// Copyright (C) 2026 Roberto Peverati.  Part of IQmol, free software under the GNU General Public
// License, version 3 or (at your option) any later version.  Developed with the assistance of
// Claude (Anthropic).
// Tests Calc::PythonSetup with a fake uv (a shell script that records its arguments), so that
// nothing is downloaded.  With a real uv and the iqmolx-engines source folder as arguments it
// also makes a real environment, which downloads Python and numpy.
//
// Usage: PythonSetupTest [<real uv> <python/iqmolx-engines/src>]

#include "Util/PythonSetup.h"
#include <QCoreApplication>
#include <QDir>
#include <QElapsedTimer>
#include <QFile>
#include <QProcess>
#include <QTemporaryDir>
#include <QTextStream>
#include <QTimer>
#include <cstdio>

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
   QFile file(path);
   file.open(QIODevice::WriteOnly);
   file.write(text.toLocal8Bit());
   file.close();
   if (executable) file.setPermissions(file.permissions() | QFileDevice::ExeOwner);
}


int main(int argc, char** argv)
{
   QCoreApplication app(argc, argv);
   QTemporaryDir dir;
   QString const root(dir.path());

   // A fake uv: records its arguments, 'venv' makes <last argument>/bin/python
   QString fake(root + "/uv");
   writeFile(fake,
      "#!/bin/sh\n"
      "echo \"$@\" >> \"$FAKE_UV_LOG\"\n"
      "[ -n \"$FAKE_UV_SLEEP\" ] && sleep 30\n"
      "if [ -n \"$FAKE_UV_FAIL\" ]; then echo boom; exit 3; fi\n"
      "if [ \"$1\" = venv ]; then\n"
      "   for last; do :; done\n"
      "   mkdir -p \"$last/bin\" && printf '#!/bin/sh\\n' > \"$last/bin/python\" && chmod +x \"$last/bin/python\"\n"
      "fi\n"
      "echo ok-from-fake-uv\n", true);
   QString log(root + "/uv.log");
   qputenv("FAKE_UV_LOG", log.toLocal8Bit());

   // The package folder: only the pyproject.toml matters
   QDir().mkpath(root + "/pkg/src");
   writeFile(root + "/pkg/pyproject.toml", "[project]\nname='x'\n");

   // Finding things
   check(PythonSetup::findUv(fake) == fake, "findUv returns the hint when it is executable");
   qputenv("IQMOLX_UV", fake.toLocal8Bit());
   check(PythonSetup::findUv("/nonexistent/uv") == fake, "findUv falls back to $IQMOLX_UV");
   qunsetenv("IQMOLX_UV");
   check(PythonSetup::defaultEnvironment().endsWith("/python"), 
      "the default environment is in the application data folder: " + PythonSetup::defaultEnvironment());
   check(PythonSetup::helperSource().isEmpty() || true, "helperSource is callable without the variable");
   qputenv("IQMOLX_ENGINES_PATH", QString(root + "/pkg/src").toLocal8Bit());
   check(PythonSetup::helperSource() == QDir::cleanPath(root + "/pkg"), "helperSource is the folder with pyproject.toml");

   // The wheel of the application: the newest one in the folder, none otherwise
   check(PythonSetup::wheelIn(root).isEmpty(), "no wheel, no wheelIn");
   QDir().mkpath(root + "/res");
   writeFile(root + "/res/iqmolx_engines-0.1.0-py3-none-any.whl", "x");
   writeFile(root + "/res/iqmolx_engines-0.2.0-py3-none-any.whl", "x");
   writeFile(root + "/res/other-9.0-py3-none-any.whl", "x");
   check(PythonSetup::wheelIn(root + "/res") == QDir::cleanPath(root + "/res/iqmolx_engines-0.2.0-py3-none-any.whl"),
      "wheelIn finds the newest iqmolx_engines wheel only: " + PythonSetup::wheelIn(root + "/res"));

   // Creating an environment: the right commands in the right order
   {
      PythonSetup setup;
      setup.setUv(fake);
      QString shown;
      setup.setProgressCallback([&](QString const& text) { shown += text; });
      QString env(root + "/env/python"), error;
      check(setup.createEnvironment(env, error), "createEnvironment succeeds " + error);
      QStringList lines(readFile(log).split('\n', Qt::SkipEmptyParts));
      check(lines.size() == 2, "two uv calls were made");
      check(lines.value(0) == "venv --python 3.12 --python-preference only-managed --allow-existing " + env,
         "uv venv with a Python that uv manages: " + lines.value(0));
      check(lines.value(1) == "pip install --python " + env + "/bin/python " + QDir::cleanPath(root + "/pkg"),
         "uv pip install of the helper into it: " + lines.value(1));
      check(!PythonSetup::pythonIn(env).isEmpty(), "the interpreter of the environment is found");
      check(shown.contains("$ uv venv") && shown.contains("ok-from-fake-uv") && shown.contains("ready"),
         "the commands and the output are shown");
   }

   // The version of the helper source: from the name of a wheel, from the pyproject.toml of a folder
   check(PythonSetup::sourceVersion(root + "/res/iqmolx_engines-0.2.0-py3-none-any.whl") == "0.2.0", 
      "the version is read from the name of a wheel");
   QDir().mkpath(root + "/ver");
   writeFile(root + "/ver/pyproject.toml", "[project]\nname = \"iqmolx-engines\"\nversion = \"1.4.2\"\n");
   check(PythonSetup::sourceVersion(root + "/ver") == "1.4.2", "the version is read from pyproject.toml");
   check(PythonSetup::sourceVersion(root + "/nothing").isEmpty(), "no source, no version");

   // Updating the helper: only the helper is reinstalled
   {
      QFile::remove(log);
      PythonSetup setup;
      setup.setUv(fake);
      QString error;
      check(setup.updateHelper("/some/env/bin/python", error), "updateHelper succeeds " + error);
      check(readFile(log).trimmed() == "pip install --python /some/env/bin/python --reinstall-package iqmolx-engines " 
         + QDir::cleanPath(root + "/pkg"), "uv reinstalls the helper only: " + readFile(log).trimmed());
   }

   // Installing requirements
   {
      QFile::remove(log);
      PythonSetup setup;
      setup.setUv(fake);
      QString error;
      check(setup.install("/some/env/bin/python", QStringList() << "torch" << "torchani", error), "install succeeds");
      check(readFile(log).trimmed() == "pip install --python /some/env/bin/python torch torchani", 
         "uv pip install of the requirements: " + readFile(log).trimmed());
   }

   // A failing uv: the reason is reported
   {
      qputenv("FAKE_UV_FAIL", "1");
      PythonSetup setup;
      setup.setUv(fake);
      QString error;
      check(!setup.install("/p/bin/python", QStringList("torch"), error) && error.contains("exit code 3") 
         && error.contains("boom"), "failure of uv is reported with its output: " + error);
      qunsetenv("FAKE_UV_FAIL");
   }

   // No uv
   {
      PythonSetup setup;
      setup.setUv("/nonexistent/uv");
      QString error;
      check(!setup.install("/p/bin/python", QStringList("torch"), error) && error.contains("Unable to start uv"),
         "a missing uv is reported: " + error);
   }

   // No helper source
   {
      qunsetenv("IQMOLX_ENGINES_PATH");
      PythonSetup setup;
      setup.setUv(fake);
      QString error;
      check(!setup.createEnvironment(root + "/other", error) && error.contains("iqmolx-engines"),
         "a missing helper package is reported: " + error);
      qputenv("IQMOLX_ENGINES_PATH", QString(root + "/pkg/src").toLocal8Bit());
   }

   // Cancel
   {
      qputenv("FAKE_UV_SLEEP", "1");
      PythonSetup setup;
      setup.setUv(fake);
      setup.setInteractive(true);
      QTimer::singleShot(400, [&]() { setup.cancel(); });
      QElapsedTimer clock;
      clock.start();
      QString error;
      bool ok(setup.install("/p/bin/python", QStringList("torch"), error));
      check(!ok && error == "Cancelled" && clock.elapsed() < 10000, 
         QString("cancel stops uv after %1 ms").arg(clock.elapsed()));
      qunsetenv("FAKE_UV_SLEEP");
   }

   // With a real uv: an environment from nothing, and the helper runs in it
   if (argc >= 3) {
      PythonSetup setup;
      setup.setUv(argv[1]);
      qputenv("IQMOLX_ENGINES_PATH", argv[2]);
      QString env(root + "/real"), error;
      QString output;
      setup.setProgressCallback([&](QString const& text) { output += text; });
      bool ok(setup.createEnvironment(env, error));
      check(ok, "real uv: environment created " + error);
      if (ok) {
         QProcess helper;
         helper.start(PythonSetup::pythonIn(env), QStringList() << "-m" << "iqmolx_engines" << "--version");
         helper.waitForFinished(30000);
         QString version(QString::fromLocal8Bit(helper.readAllStandardOutput()));
         check(version.contains("protocol 1"), "real uv: the helper runs in the new environment: " + version.trimmed());
      }
   }

   std::printf("\n%s\n", s_failures ? "FAILED" : "all tests passed");
   return s_failures ? 1 : 0;
}
