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
#include "Util/PythonSetup.h"
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDesktopWidget>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QTimer>
#include <QJsonDocument>
#include <QLabel>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QResizeEvent>
#include <QSettings>
#include <QSpinBox>
#include <QStandardPaths>
#include <QTemporaryDir>
#include <QDir>
#include <QFile>
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
   QStandardPaths::setTestModeEnabled(true);   // the application data folder is not the real one
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
   QLineEdit* pythonEdit(0);
   for (QLineEdit* e : d->findChildren<QLineEdit*>()) if (e->placeholderText().startsWith("python3")) pythonEdit = e;
   if (pythonEdit) {
      pythonEdit->setText("/nonexistent/python");
      pythonEdit->editingFinished();   // rebuilds the settings, so no loop over the line edits
   }
   d->refreshEngines();
   check(!button(d, "Run")->isEnabled(), "Run is disabled when the helper cannot start");
   check(!d->prepare(problem) && problem.contains("Python"), "prepare explains: " + problem);

   // --- an engine with many settings (PySCF has five) must not make the window taller than the screen:
   //     the Run and Save buttons were hidden below the Dock ---
   {
      if (pythonEdit) pythonEdit->setText(argv[1]);
      Preferences::PythonEnginePath(argv[1]);
      d->refreshEngines();
      QComboBox* combo(0);
      for (QComboBox* c : d->findChildren<QComboBox*>()) if (c->findData("pyscf") >= 0) combo = c;
      int index(combo ? combo->findData("pyscf") : -1);
      if (index < 0 || combo->itemText(index).contains("not available")) {
         std::printf("SKIP: tall engine (PySCF is not available in this Python)\n");
      }else {
         combo->setCurrentIndex(index);
         QRect screen(QApplication::desktop()->availableGeometry(d));
         d->showWindow();
         QCoreApplication::processEvents();
         check(d->minimumSizeHint().height() <= screen.height() - 40,
            QString("the window can be small enough for the screen: minimum height %1, screen %2")
               .arg(d->minimumSizeHint().height()).arg(screen.height()));
         check(screen.contains(d->frameGeometry().bottomLeft()) && d->frameGeometry().bottom() <= screen.bottom(),
            "the whole window is on the screen");
         QRect full(QApplication::desktop()->screenGeometry(d));
         check(d->frameGeometry().bottom() <= full.top() + qRound(0.80 * full.height()),
            "the window stays in the upper 80% of the screen, clear of a large Dock");
         QPushButton* run(button(d, "Run"));
         check(run && screen.contains(run->mapToGlobal(run->rect().bottomRight())), 
            "the Run button is on the screen, not behind the Dock");
         QPushButton* saveOutput(button(d, "Save Output..."));
         check(saveOutput && screen.contains(saveOutput->mapToGlobal(saveOutput->rect().bottomRight())),
            "the Save Output button is on the screen");
         combo->setCurrentIndex(combo->findData("demo"));
      }
   }

   // --- options that depend on others: the Functional field appears for "other" (PySCF) ---
   {
      QComboBox* engineCombo(0);
      for (QComboBox* c : d->findChildren<QComboBox*>()) if (c->findData("pyscf") >= 0) engineCombo = c;
      int index(engineCombo ? engineCombo->findData("pyscf") : -1);
      if (index < 0 || engineCombo->itemText(index).contains("not available")) {
         std::printf("SKIP: dependent options (PySCF is not available in this Python)\n");
      }else {
         engineCombo->setCurrentIndex(index);
         d->showWindow();
         QComboBox* method(0);
         for (QComboBox* c : d->findChildren<QComboBox*>()) if (c->findText("wb97m-v") >= 0 && c->findText("other") >= 0) method = c;
         check(method != 0, "the method menu has wb97m-v and other");
         auto field = [&](QString const& placeholder) -> QLineEdit* {
            for (QLineEdit* e : d->findChildren<QLineEdit*>()) if (e->placeholderText() == placeholder) return e;
            return 0;
         };
         auto gridSpin = [&]() -> QSpinBox* {
            for (QLabel* l : d->findChildren<QLabel*>()) if (l->text().startsWith("Grid of the VV10 part")) {
               for (QSpinBox* s : d->findChildren<QSpinBox*>()) if (s->toolTip().contains("nonlocal (VV10)")) return s;
            }
            return 0;
         };
         method->setCurrentText("b3lyp");
         QLineEdit* name(field("cam-b3lyp"));
         check(name != 0 && !name->isVisibleTo(d), "the Functional field is hidden for b3lyp");
         check(gridSpin() && !gridSpin()->isVisibleTo(d), "the grid of the VV10 part is hidden for b3lyp");
         check(button(d, "Run")->isEnabled(), "Run is on for b3lyp");

         method->setCurrentText("other");
         name = field("cam-b3lyp");
         check(name != 0 && name->isVisibleTo(d), "the Functional field appears when other is chosen");
         check(gridSpin() && gridSpin()->isVisibleTo(d), "the VV10 grid is shown for other too");
         check(!button(d, "Run")->isEnabled(), "Run waits for the name");
         QString hint;
         for (QLabel* l : d->findChildren<QLabel*>()) if (l->text().contains("Fill in")) hint = l->text();
         check(hint.contains("Functional"), "the window says what to fill in: " + hint);
         name->setText("cam-b3lyp");
         check(button(d, "Run")->isEnabled(), "Run is on once the functional is named");

         method->setCurrentText("wb97m-v");
         name = field("cam-b3lyp");
         check(name && !name->isVisibleTo(d) && button(d, "Run")->isEnabled(),
            "choosing wb97m-v hides the name again and Run stays on (its text is kept but not needed)");
         check(gridSpin() && gridSpin()->isVisibleTo(d), "the VV10 grid is shown for wb97m-v");
         method->setCurrentText("b3lyp");
         engineCombo->setCurrentIndex(engineCombo->findData("demo"));
      }
   }

   // --- Save Output...: cancelling the file dialog must leave the window open ---
   {
      QPushButton* save(button(d, "Save Output..."));
      check(save != 0 && save->isVisible(), "the Save Output button is there");
      bool dialogWasOpen(false);
      QTimer::singleShot(400, [&]() {                // press Cancel in the file dialog, once it is up
         for (QWidget* w : QApplication::topLevelWidgets()) {
            QFileDialog* f(qobject_cast<QFileDialog*>(w));
            if (f && f->isVisible()) { dialogWasOpen = true; f->reject(); }
         }
      });
      if (save) save->click();
      check(dialogWasOpen, "a file dialog opened");
      check(d->isVisible(), "the Python Engines window is still open after cancelling the save");
   }

   // --- the size: room for the output, and the size the user chose is remembered ---
   {
      QPlainTextEdit* output(d->findChild<QPlainTextEdit*>());
      check(output && output->minimumHeight() >= 120, "the output area has room to be read");
      QCoreApplication::processEvents();
      check(!Preferences::PythonEngineWindowSize().isValid(), 
         "the size the window gave itself (at the start, when shown) is not stored as the user's choice");
      d->resize(800, 700);                                  // the user dragging the window
      QCoreApplication::processEvents();
      check(Preferences::PythonEngineWindowSize() == d->size(), "the size the user chose is stored");
      check(d->width() >= 600 && d->height() >= 600, "the window is not small");
   }

   // --- a required option: Run stays off until it is filled in (needs the ASE engine) ---
   {
      if (pythonEdit) pythonEdit->setText(argv[1]);
      Preferences::PythonEnginePath(argv[1]);
      d->refreshEngines();
      QComboBox* combo(0);
      for (QComboBox* c : d->findChildren<QComboBox*>()) if (c->findData("ase") >= 0) combo = c;
      int index(combo ? combo->findData("ase") : -1);
      if (index < 0 || combo->itemText(index).contains("not available")) {
         std::printf("SKIP: required option (the ASE engine is not available in this Python)\n");
      }else {
         combo->setCurrentIndex(index);
         QLineEdit* calculator(0);
         for (QLineEdit* e : d->findChildren<QLineEdit*>()) if (e->placeholderText() == "ase.calculators.emt:EMT") calculator = e;
         check(calculator != 0, "the placeholder of the option is shown in its field");
         check(button(d, "Run") && !button(d, "Run")->isEnabled(), "Run is off while a required field is empty");
         QString hint;
         for (QLabel* l : d->findChildren<QLabel*>()) if (l->text().contains("Fill in")) hint = l->text();
         check(hint.contains("Calculator class"), "the window says what to fill in: " + hint);
         if (calculator) {
            calculator->setText("ase.calculators.emt:EMT");
            check(button(d, "Run")->isEnabled(), "Run is on once the field is filled in");
            calculator->setText("   ");
            check(!button(d, "Run")->isEnabled(), "blank is empty again");
            calculator->setText(QString());
         }
         combo->setCurrentIndex(combo->findData("demo"));      // leave things as they were
      }
   }

   // --- the Python environment, made with a fake uv that links the interpreter of this test ---
   {
      QTemporaryDir tmp;
      QString fake(tmp.path() + "/uv"), uvLog(tmp.path() + "/uv.log");
      QFile script(fake);
      script.open(QIODevice::WriteOnly);
      script.write(QString(
         "#!/bin/sh\n"
         "echo \"$@\" >> \"%1\"\n"
         "if [ \"$1\" = venv ]; then\n"
         "   for last; do :; done\n"
         "   mkdir -p \"$last/bin\" && touch \"$last/pyvenv.cfg\" && ln -sf \"%2\" \"$last/bin/python\"\n"
         "fi\n"
         "echo fake-uv-output\n").arg(uvLog, argv[1]).toLocal8Bit());
      script.close();
      script.setPermissions(script.permissions() | QFileDevice::ExeOwner);
      qputenv("IQMOLX_UV", fake.toLocal8Bit());

      QString environment(Calc::PythonSetup::defaultEnvironment());
      check(environment.contains(".qttest") && !QFile::exists(environment), 
         "the test environment folder is a test one and does not exist yet: " + environment);

      // Nothing configured: the window offers to make an environment
      Preferences::PythonEnginePath(QString());
      if (pythonEdit) pythonEdit->setText(QString());
      d->refreshEngines();
      QPushButton* setup(button(d, "Set Up Environment"));
      check(setup && !setup->isHidden() && setup->isEnabled(), "Set Up Environment is offered, as uv is there");

      if (!setup) { std::printf("FAILED: no setup button, the rest is skipped\n"); return 1; }
      setup->click();
      QFile logFile(uvLog);
      logFile.open(QIODevice::ReadOnly);
      QString uvCalls(QString::fromLocal8Bit(logFile.readAll()));
      check(uvCalls.contains("venv --python 3.12 --python-preference only-managed") 
         && uvCalls.contains("pip install --python " + environment + "/bin/python"), "uv was asked for the environment and the helper");
      check(!Calc::PythonSetup::pythonIn(environment).isEmpty(), "the environment has an interpreter");
      check(setup->isHidden(), "the button goes away once there is an environment");
      QString shown;
      for (QPlainTextEdit* t : d->findChildren<QPlainTextEdit*>()) shown += t->toPlainText();
      check(shown.contains("$ uv venv") && shown.contains("fake-uv-output"), "the commands and the output of uv are shown");
      check(Preferences::PythonEnginePath() == Calc::PythonSetup::pythonIn(environment), 
         "the new environment is the Python of the window");
      check(button(d, "Run") && button(d, "Run")->isEnabled(), "the engines are found in the new environment");
      if (!button(d, "Run")) return 1;

      // An engine whose packages are missing: selectable, with an install button
      QComboBox* combo(0);
      for (QComboBox* c : d->findChildren<QComboBox*>()) if (c->findData("ase") >= 0) combo = c;
      check(combo != 0, "the engine list has the ASE adapter");
      if (combo) combo->setCurrentIndex(combo->findData("ase"));
      check(!button(d, "Run")->isEnabled(), "Run is off for an engine that is not available");
      QPushButton* install(button(d, "Install ase"));
      check(install && !install->isHidden() && install->isEnabled(), "Install ase is offered");
      if (install) install->click();
      logFile.close();
      logFile.open(QIODevice::ReadOnly);
      uvCalls = QString::fromLocal8Bit(logFile.readAll());
      check(uvCalls.contains("pip install --python " + environment + "/bin/python ase"), 
         "uv was asked to install the requirements of the engine");
      if (combo) combo->setCurrentIndex(combo->findData("demo"));
      QPushButton* installAfter(button(d, "Install ase"));
      check(button(d, "Run")->isEnabled() && (!installAfter || !installAfter->isVisible()), 
         "a working engine has Run and no install button");

      // The folder was made by this test, in Qt's test location
      if (environment.contains(".qttest") && environment.endsWith("/python")) QDir(environment).removeRecursively();
   }

   std::printf("\n%s\n", s_failures ? "FAILED" : "all tests passed");
   return s_failures ? 1 : 0;
}
