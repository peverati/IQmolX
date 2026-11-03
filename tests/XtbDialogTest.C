// Copyright (C) 2026 Roberto Peverati.  Part of IQmol, free software under the GNU General Public
// License, version 3 or (at your option) any later version.  Developed with the assistance of
// Claude (Anthropic).
// Tests the logic of the xtb window (settings, Run/Cancel, running and finished states)
// on Qt's offscreen platform.  The settings are stored under a different name, so the
// real IQmol preferences are not touched.

#include "Util/Preferences.h"
#include "Util/XtbResultDialog.h"
#include <QApplication>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QLineEdit>
#include <QTimer>
#include <QSpinBox>
#include <QComboBox>
#include <QDir>
#include <QLabel>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSettings>
#include <QFile>
#include <cstdio>

using namespace IQmol;

static int s_failures(0);

static void check(bool ok, char const* what)
{
   std::printf("%s: %s\n", ok ? "PASS" : "FAIL", what);
   if (!ok) ++s_failures;
}

static QPushButton* button(QWidget* parent, QString const& text)
{
   for (QPushButton* b : parent->findChildren<QPushButton*>()) {
       if (b->text() == text) return b;
   }
   return 0;
}

static QComboBox* comboWith(QWidget* parent, QString const& dataValue)
{
   for (QComboBox* c : parent->findChildren<QComboBox*>()) {
       if (c->findData(dataValue) >= 0) return c;
   }
   return 0;
}


int main(int argc, char** argv)
{
   qputenv("QT_QPA_PLATFORM", "offscreen");
   QApplication app(argc, argv);

   // Settings are stored under a name of their own, never in the real IQmol preferences
   // (QSettings::setPath does not apply to the native format used on macOS)
   QString const organization("iqmol-test.invalid");
   QString const application("IQmolXtbDialogTest" + QString::number(QCoreApplication::applicationPid()));
   Preferences::UseSettingsName(organization, application);
   QSettings(QSettings::UserScope, organization, application).clear();

   // Whatever job was used last time, the window starts with a geometry optimization
   Preferences::XtbJob("optfreq");
   Preferences::XtbMethod("GFN1-xTB");
   Preferences::XtbSolvent("acetonitrile");
   Preferences::XtbBackend("program");
   Preferences::XtbElectronicTemperature(750.0);

   // With nothing stored the window starts with GFN2-xTB, a geometry optimization in the gas phase
   XtbResultDialog* d(XtbResultDialog::instance());
   check(Preferences::XtbJob() == "optimize", "a job stored by an earlier session is replaced by a geometry optimization");
   {
      QComboBox* defaultMethod(comboWith(d, "GFN2-xTB"));
      check(Preferences::XtbMethod() == "GFN2-xTB", "the default method is GFN2-xTB");
      check(defaultMethod && defaultMethod->currentData() == "GFN2-xTB", "the window starts with GFN2-xTB");
      check(defaultMethod && defaultMethod->count() == 2 && defaultMethod->findData("GFN-FF") < 0,
         "GFN-FF is not an xtb window method any more (it is a force field)");
      QList<QDoubleSpinBox*> temps(d->findChildren<QDoubleSpinBox*>());
      check(temps.size() == 1 && temps[0]->value() == 300.0 && Preferences::XtbElectronicTemperature() == 300.0,
         "the electronic temperature starts at the xtb default of 300 K");

      // A stored GFN-FF method from earlier versions and xTB default force fields are replaced
      Preferences::XtbMethod("GFN-FF");
      check(Preferences::XtbMethod() == "GFN2-xTB", "a stored GFN-FF method falls back to GFN2-xTB");
      Preferences::DefaultForceField("GFN2-xTB");
      check(Preferences::DefaultForceField() == "UFF", "a stored xTB default force field falls back to UFF");
      Preferences::DefaultForceField("GFN-FF");
      check(Preferences::DefaultForceField() == "GFN-FF", "GFN-FF is a valid default force field");
      Preferences::DefaultForceField("UFF");
      QComboBox* defaultJob(comboWith(d, "optimize"));
      check(defaultJob && defaultJob->currentData() == "optimize", "the window starts with a geometry optimization");
      QComboBox* defaultSolvent(comboWith(d, "water"));
      check(defaultSolvent && defaultSolvent->currentData() == "", "the window starts in the gas phase");
      check(Preferences::XtbMethod() == "GFN2-xTB" && Preferences::XtbSolvent().isEmpty(),
         "method and solvent of an earlier session are replaced by the defaults");
      QComboBox* defaultBackend(comboWith(d, "program"));
      check(defaultBackend && defaultBackend->currentData() == (defaultBackend->count() > 1 ? "library" : "program"),
         "the window starts with the built-in library when there is one");
      QCheckBox* molecule(0);
      for (QCheckBox* c : d->findChildren<QCheckBox*>()) if (c->text().contains("charge and spin")) molecule = c;
      check(molecule && molecule->isChecked(), "the charge and spin of the molecule are used by default");
      QDoubleSpinBox* temp(d->findChild<QDoubleSpinBox*>());
      d->show();
      check(temp && (defaultBackend->count() > 1) == !temp->isVisible(),
         "the electronic temperature is only shown for the xtb program");
   }
   Preferences::XtbBackend("program");

   // Stored settings are shown when the window opens
   Preferences::XtbMethod("GFN1-xTB");
   Preferences::XtbSolvent("water");
   Preferences::XtbJob("singlepoint");

   d->showWindow();
   check(d->isVisible(), "the window opens");

   QComboBox* method(comboWith(d, "GFN1-xTB"));
   QComboBox* solvent(comboWith(d, "water"));
   QComboBox* job(comboWith(d, "singlepoint"));
   check(method && method->currentData() == "GFN1-xTB", "method shows the stored setting");
   check(solvent && solvent->currentData() == "water", "solvent shows the stored setting");
   check(job && job->currentData() == "singlepoint", "job shows the stored setting");
   check(solvent && solvent->count() == 25, "gas phase and the 24 ALPB solvents are listed");

   QPushButton* run(button(d, "Run"));
   QPushButton* cancel(button(d, "Cancel"));
   check(run && run->isVisible(), "Run is available");
   check(cancel && !cancel->isVisible(), "Cancel is hidden until a calculation runs");

   // Changing a setting stores it
   solvent->setCurrentIndex(solvent->findData("methanol"));
   method->setCurrentIndex(method->findData("GFN2-xTB"));
   job->setCurrentIndex(job->findData("optimize"));
   check(Preferences::XtbSolvent() == "methanol", "solvent is stored when changed");
   check(Preferences::XtbMethod() == "GFN2-xTB", "method is stored when changed");

   // The electronic temperature is next to the charge and spin, and stored when changed
   QDoubleSpinBox* temperature(d->findChild<QDoubleSpinBox*>());
   check(temperature && temperature->suffix() == " K" && temperature->minimum() == 0.0, "the electronic temperature is entered in K");
   temperature->setValue(50.0);
   check(Preferences::XtbElectronicTemperature() == 50.0, "the electronic temperature is stored when changed");
   check(Preferences::XtbJob() == "optimize", "job is stored when changed");

   // Charge and spin: the molecule's values by default, an override on request
   QCheckBox* useMolecule(0);
   QCheckBox* orbitals(0);
   for (QCheckBox* c : d->findChildren<QCheckBox*>()) {
       if (c->text().contains("orbitals")) orbitals = c; else if (c->text().contains("charge and spin")) useMolecule = c;
   }
   check(orbitals && orbitals->isChecked() == Preferences::XtbOrbitals() && Preferences::XtbOrbitals(),
      "orbitals are calculated by default");
   if (orbitals) {
      orbitals->setChecked(false);
      check(!Preferences::XtbOrbitals(), "the orbitals setting is stored");
      orbitals->setChecked(true);
   }
   QList<QSpinBox*> spins(d->findChildren<QSpinBox*>());
   check(useMolecule && useMolecule->isChecked(), "the molecule's charge and spin are used by default");
   check(spins.size() == 2 && !spins[0]->isEnabled() && !spins[1]->isEnabled(), 
      "the charge and spin boxes are disabled while the molecule's values are used");
   int q(99), m(99);
   check(!d->chargeAndSpin(q, m) && q == 99 && m == 99, "no override is reported by default");
   QSpinBox* chargeSpin(spins[0]->minimum() < 0 ? spins[0] : spins[1]);
   QSpinBox* multiplicitySpin(chargeSpin == spins[0] ? spins[1] : spins[0]);
   useMolecule->setChecked(false);
   check(chargeSpin->isEnabled() && multiplicitySpin->isEnabled(), "unchecking enables the boxes");
   chargeSpin->setValue(1);
   multiplicitySpin->setValue(2);
   check(d->chargeAndSpin(q, m) && q == 1 && m == 2, "the entered charge and spin are reported");
   check(multiplicitySpin->minimum() == 1 && chargeSpin->minimum() < 0, "ranges allow anions and start at singlet");
   useMolecule->setChecked(true);

   // Run with nothing to calculate says so
   int runRequests(0);
   QObject::connect(d, &XtbResultDialog::runRequested, [&runRequests]() { ++runRequests; });
   run->click();
   check(runRequests == 1, "Run asks for the calculation");
   QLabel* summary(0);
   for (QLabel* l : d->findChildren<QLabel*>()) if (l->text().contains("Nothing was calculated")) summary = l;
   check(summary != 0, "pressing Run with no molecule explains that nothing happened");

   // A calculation: started -> running state -> finished
   int cancels(0);
   QObject::connect(d, &XtbResultDialog::cancelRequested, [&cancels]() { ++cancels; });

   d->started("GFN-FF/ALPB(methanol)", "geometry optimization", "water (3 atoms)", 0, 1);
   check(d->isVisible(), "the window is shown while running");
   check(!run->isVisible() && cancel->isVisible(), "Run is replaced by Cancel while running");
   check(!method->isEnabled() && !solvent->isEnabled() && !job->isEnabled(), "settings are locked while running");

   check(!useMolecule->isEnabled(), "the charge and spin choice is locked while running");
   check(!temperature->isEnabled(), "the electronic temperature is locked while running");
   {
      QString running;
      for (QLabel* l : d->findChildren<QLabel*>()) if (l->text().startsWith("Running on")) running = l->text();
      check(running.contains("electronic temperature 50 K"), "the temperature in use is shown while running");
   }
   d->appendOutput("cycle 1\n");
   d->appendOutput("cycle 2\n");
   QPlainTextEdit* output(d->findChild<QPlainTextEdit*>());
   check(output && output->toPlainText() == "cycle 1\ncycle 2\n", "output accumulates");

   cancel->click();
   check(cancels == 1, "Cancel button requests cancellation");
   d->close();
   check(cancels == 2 && d->isVisible(), "closing the window while running requests cancellation instead");

   IQmol::Calc::Result result;
   result.success = true;
   result.converged = true;
   result.energy = -5.070544;
   result.initialEnergy = -5.05;
   result.hasGsolv = true;
   result.gsolv = -0.022;
   result.cycles = 6;
   result.seconds = 0.1;
   d->finished(result, "GFN-FF/ALPB(methanol)", "geometry optimization", true);
   check(run->isVisible() && !cancel->isVisible(), "Run is available again when finished");
   check(method->isEnabled() && solvent->isEnabled() && job->isEnabled(), "settings are unlocked when finished");
   bool hasResult(false);
   for (QLabel* l : d->findChildren<QLabel*>()) {
       if (l->text().contains("-5.070544") && l->text().contains("Solvation free energy")) hasResult = true;
   }
   check(hasResult, "the summary shows the energy and the solvation free energy");
   bool hasChargeSpin(false);
   for (QLabel* l : d->findChildren<QLabel*>()) {
       if (l->text().contains("Charge, spin multiplicity") && l->text().contains("0, 1")) hasChargeSpin = true;
   }
   check(hasChargeSpin, "the summary shows the charge and spin that were used");
   bool hasTemperature(false);
   for (QLabel* l : d->findChildren<QLabel*>()) {
       if (l->text().contains("Electronic temperature") && l->text().contains("50 K")) hasTemperature = true;
   }
   check(hasTemperature, "the summary shows the electronic temperature that was used");

   // The values in use are shown, so that unchecking the box starts from them
   d->started("GFN2-xTB", "single point", "radical (4 atoms)", 0, 2);
   check(chargeSpin->value() == 0 && multiplicitySpin->value() == 2, "the boxes show the values in use");
   d->finished(result, "GFN2-xTB", "single point", false);
   useMolecule->setChecked(false);
   d->started("GFN2-xTB", "single point", "radical (4 atoms)", 1, 1);
   bool entered(false);
   d->finished(result, "GFN2-xTB", "single point", false);
   for (QLabel* l : d->findChildren<QLabel*>()) if (l->text().contains("entered in this window")) entered = true;
   check(entered, "the summary says when the values were entered by hand");
   useMolecule->setChecked(true);

   // A second calculation can follow
   d->started("GFN-FF", "single point", "water (3 atoms)", 0, 1);
   check(!run->isVisible(), "a second calculation can be started");
   IQmol::Calc::Result failed;
   failed.message = "xtb failed";
   d->finished(failed, "GFN-FF", "single point", false);
   check(run->isVisible(), "failure returns to the setup state");

   // Vibrational frequencies: choices, summary and warnings
   check(job->count() == 4 && job->findData("frequencies") >= 0 && job->findData("optfreq") >= 0,
      "the job list has single point, optimization, frequencies and optimization + frequencies");
   {
      d->started("GFN2-xTB", "vibrational frequencies", "water (3 atoms)", 0, 1);
      IQmol::Calc::Result freq;
      freq.success = true;
      freq.converged = true;
      freq.energy = -5.0705;
      freq.gradientNorm = 0.0004;
      freq.hasVibrations = true;
      freq.vibrations.frequencies << 1538.9 << 3643.1 << 3651.8;
      freq.vibrations.intensities << 133.3 << 6.8 << 16.6;
      freq.vibrations.zeroPointEnergy = 0.020125;
      d->finished(freq, "GFN2-xTB", "vibrational frequencies", false);
      QString text;
      for (QLabel* l : d->findChildren<QLabel*>()) if (l->text().contains("Vibrational modes")) text = l->text();
      check(text.contains("Vibrational modes") && text.contains("1538.9 to 3651.8"), "the summary shows the number and range of the modes");
      check(text.contains("Zero point energy") && text.contains("12.63 kcal/mol"), "the summary shows the zero point energy");
      check(text.contains("Strongest IR band") && text.contains("1538.9") && text.contains("133.3 km/mol"), "the summary shows the strongest IR band");
      check(text.contains("Imaginary frequencies") && !text.contains("not a minimum"), "no warning for a minimum");
      check(text.contains("Frequencies</i> entry"), "the summary says where to find the modes");

      d->started("GFN2-xTB", "vibrational frequencies", "x (3 atoms)", 0, 1);
      freq.gradientNorm = 0.02;
      freq.vibrations.frequencies.clear();
      freq.vibrations.frequencies << -85.2 << 1538.9 << 3651.8;
      freq.vibrations.intensities.clear();
      freq.vibrations.intensities << 1.0 << 133.3 << 16.6;
      d->finished(freq, "GFN2-xTB", "vibrational frequencies", false);
      QString warning;
      for (QLabel* l : d->findChildren<QLabel*>()) if (l->text().contains("Vibrational modes")) warning = l->text();
      check(warning.contains("-85.2") && warning.contains("not a minimum"), "an imaginary frequency is reported with a warning");
      check(warning.contains("not a stationary point"), "frequencies at a non-stationary geometry are flagged");

      // After an optimization the geometry is stationary by construction, so no such warning
      d->started("GFN2-xTB", "geometry optimization and frequencies", "x (3 atoms)", 0, 1);
      d->finished(freq, "GFN2-xTB", "geometry optimization and frequencies", true);
      QString afterOpt;
      for (QLabel* l : d->findChildren<QLabel*>()) if (l->text().contains("Vibrational modes")) afterOpt = l->text();
      check(!afterOpt.contains("not a stationary point"), "no stationarity warning after an optimization");
   }

   // The same window comes back with the settings after being closed
   d->close();
   check(!d->isVisible(), "the window closes when idle");
   Preferences::XtbSolvent("");
   Preferences::XtbElectronicTemperature(120.0);
   d->showWindow();
   check(d->isVisible() && solvent->currentData() == "", "reopening shows the current settings");
   check(temperature->value() == 120.0, "reopening shows the stored electronic temperature");

   // Saving the files of the xtb program
   {
      QPushButton* saveAll(button(d, "Save All Files..."));
      check(saveAll != 0 && button(d, "Save Output...") != 0, "there are buttons for the output and for all files");
      QDir tmp(QDir::tempPath());
      QString first(tmp.filePath("xtbdialogtest-first")), second(tmp.filePath("xtbdialogtest-second"));
      QDir().mkpath(first); QDir().mkpath(second);
      Calc::Result r;
      r.success = true; r.converged = true; r.energy = -1.0;
      r.directory = first; r.temporaryDirectory = true;
      d->started("GFN2-xTB", "single point", "x", 0, 1);
      d->finished(r, "GFN2-xTB", "single point", false);
      check(saveAll && saveAll->isEnabled(), "Save All Files is available after a run of the program");
      r.directory = second;
      d->started("GFN2-xTB", "single point", "x", 0, 1);
      d->finished(r, "GFN2-xTB", "single point", false);
      check(!QDir(first).exists() && QDir(second).exists(), "the temporary directory of the earlier run is removed");
      Calc::Result lib;
      lib.success = true; lib.converged = true;
      d->started("GFN2-xTB", "single point", "x", 0, 1);
      d->finished(lib, "GFN2-xTB", "single point", false);
      check(saveAll && !saveAll->isEnabled() && !QDir(second).exists(), "without files (library) there is nothing to save");
   }

   // The settings of the xtb program: path, options and working directory
   {
      QLineEdit* programEdit(0); QLineEdit* optionsEdit(0); QLineEdit* workEdit(0);
      for (QLineEdit* e : d->findChildren<QLineEdit*>()) {
          if (e->placeholderText().startsWith("path")) programEdit = e;
          else if (e->placeholderText().startsWith("e.g.")) optionsEdit = e;
          else if (e->placeholderText().startsWith("empty")) workEdit = e;
      }
      check(programEdit && optionsEdit && workEdit, "the program, options and directory fields exist");
      if (programEdit && optionsEdit && workEdit) {
         d->show();
         QComboBox* backend(comboWith(d, "program"));
         check(programEdit->isVisible() && optionsEdit->isVisible() && workEdit->isVisible(),
            "the fields are shown for the xtb program");
         optionsEdit->setText("--etemp 500 --verbose");
         optionsEdit->editingFinished();
         workEdit->setText("/tmp/xtb-work");
         workEdit->editingFinished();
         programEdit->setText("/opt/xtb/bin/xtb");
         programEdit->editingFinished();
         check(Preferences::XtbOptions() == "--etemp 500 --verbose", "the options are stored");
         check(Preferences::XtbWorkDir() == "/tmp/xtb-work", "the working directory is stored");
         check(Preferences::XtbPath() == "/opt/xtb/bin/xtb", "the program path is stored");
         if (backend && backend->count() > 1) {
            backend->setCurrentIndex(1);   // library
            check(!programEdit->isVisible() && !optionsEdit->isVisible() && !workEdit->isVisible() && !orbitals->isVisible(),
               "the fields are hidden for the library");
            backend->setCurrentIndex(0);
         }
         Preferences::XtbOptions(QString()); Preferences::XtbWorkDir(QString()); Preferences::XtbPath(QString());
      }
   }

   // Remove the settings file of this test
   {
      QSettings settings(QSettings::UserScope, organization, application);
      settings.clear();
      settings.sync();
      QFile::remove(settings.fileName());
   }

   std::printf("%d failure(s)\n", s_failures);
   return s_failures == 0 ? 0 : 1;
}
