// Copyright (C) 2026 Roberto Peverati.  Part of IQmol, free software under the GNU General Public
// License, version 3 or (at your option) any later version.  Developed with the assistance of
// Claude (Anthropic).
// Tests the Python engines window against the demo engine of the iqmolx-engines helper, on Qt's
// offscreen platform.  The settings are stored under a different name, so the real IQmol
// preferences are not touched.
//
// Usage: PythonDialogTest <python> <python/iqmolx-engines/src>

#include "Util/Preferences.h"
#include "Util/PythonEngineDialog.h"
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QJsonDocument>
#include <QLabel>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSettings>
#include <QSpinBox>
#include <cstdio>

using namespace IQmol;

static int s_failures(0);

static void check(bool ok, QString const& what)
{
   std::printf("%s: %s\n", ok ? "PASS" : "FAIL", qPrintable(what));
   if (!ok) ++s_failures;
}

static QPushButton* button(QWidget* parent, QString const& text)
{
   for (QPushButton* b : parent->findChildren<QPushButton*>()) {
       if (b->text() == text) return b;
   }
   return 0;
}

static Calc::Input hydrogen(double r)
{
   Calc::Input input;
   input.atomicNumbers << 1 << 1;
   input.coordinates << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << r;
   return input;
}


int main(int argc, char** argv)
{
   qputenv("QT_QPA_PLATFORM", "offscreen");
   if (argc < 3) { std::printf("usage: PythonDialogTest <python> <helper src dir>\n"); return 2; }
   qputenv("IQMOLX_ENGINES_PATH", argv[2]);
   QApplication app(argc, argv);

   QString const organization("iqmol-test.invalid");
   QString const application("IQmolPythonDialogTest" + QString::number(QCoreApplication::applicationPid()));
   Preferences::UseSettingsName(organization, application);
   QSettings(QSettings::UserScope, organization, application).clear();
   Preferences::PythonEnginePath(argv[1]);

   PythonEngineDialog* d(PythonEngineDialog::instance());
   check(Preferences::PythonEngineJob() == "optimize", "the default job is a geometry optimization");
   check(button(d, "Run") && !button(d, "Run")->isEnabled(), "Run is disabled before an engine is known");

   // The window finds the engines itself the first time it is shown
   d->showWindow();
   QComboBox* engines(0);
   for (QComboBox* c : d->findChildren<QComboBox*>()) if (c->findData("demo") >= 0) engines = c;
   check(engines != 0, "the engines come from the helper");
   check(engines && engines->currentData() == "demo", "the demo engine is selected");
   check(button(d, "Run")->isEnabled(), "Run is enabled once an engine is selected");
   check(engines && engines->findData("ase") >= 0, "the ASE adapter is listed (available or not)");

   // The settings form is built from the option descriptions
   QDoubleSpinBox* depth(0);
   for (QDoubleSpinBox* s : d->findChildren<QDoubleSpinBox*>()) {
       if (qAbs(s->value() - 4.0) < 1e-9) depth = s;
   }
   check(depth != 0, "a number field for the well depth, with its default");
   check(depth && depth->minimum() == 0.1 && depth->maximum() == 20.0, "limits from the engine");
   QSpinBox* delay(0);
   for (QSpinBox* s : d->findChildren<QSpinBox*>()) if (s->maximum() == 60000) delay = s;
   check(delay != 0, "an integer field for the delay (a compute option)");

   // Changed settings reach the engine and are remembered
   if (depth) depth->setValue(2.0);
   QString problem;
   check(d->prepare(problem), "prepare " + problem);
   Calc::Result r(d->engine().energy(hydrogen(0.62)));
   check(r.success && qAbs(r.energy * 27.211386245988 + 2.0) < 1e-9, "the engine uses the window's well depth");
   QJsonDocument stored(QJsonDocument::fromJson(Preferences::PythonEngineOptions().toUtf8()));
   check(qAbs(stored.object().value("demo").toObject().value("well_depth").toDouble() - 2.0) < 1e-12,
      "the setting is stored per engine");

   // The run protocol of the window, as Molecule::runPythonJob uses it
   int runs(0);
   QObject::connect(d, &PythonEngineDialog::runRequested, [&]() {
      ++runs;
      d->started(d->engineTitle(), "single point", "water (3 atoms)", 0, 1);
      d->appendOutput("working\n");
      Calc::Result result(d->engine().energy(hydrogen(0.62)));
      d->finished(result, d->engineTitle(), "single point", false);
   });
   button(d, "Run")->click();
   check(runs == 1, "Run emits runRequested");
   check(!d->isModal() && d->isVisible(), "the window is back to normal when the job is done");
   bool summary(false);
   for (QLabel* l : d->findChildren<QLabel*>()) if (l->text().contains("Energy") && l->text().contains("Eh")) summary = true;
   check(summary, "the summary shows the energy");
   check(d->findChild<QPlainTextEdit*>()->toPlainText().contains("working"), "output is shown");

   // Charge and spin of the window replace those of the molecule when asked
   int charge(0), multiplicity(1);
   check(!d->chargeAndSpin(charge, multiplicity), "molecule values are used by default");
   for (QCheckBox* c : d->findChildren<QCheckBox*>()) if (c->text().startsWith("Use the charge")) c->setChecked(false);
   charge = 7;
   check(d->chargeAndSpin(charge, multiplicity) && charge != 7, "window values are used when asked");

   // A wrong interpreter is reported in the window
   Preferences::PythonEnginePath("/nonexistent/python");
   for (QLineEdit* e : d->findChildren<QLineEdit*>()) if (e->placeholderText().startsWith("python3")) {
      e->setText("/nonexistent/python");
      e->editingFinished();
   }
   d->refreshEngines();
   check(!button(d, "Run")->isEnabled(), "Run is disabled when the helper cannot start");
   check(!d->prepare(problem) && problem.contains("Python"), "prepare explains: " + problem);

   std::printf("\n%s\n", s_failures ? "FAILED" : "all tests passed");
   return s_failures ? 1 : 0;
}
