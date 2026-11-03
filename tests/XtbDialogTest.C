// Tests the logic of the xtb window (settings, Run/Cancel, running and finished states)
// on Qt's offscreen platform.  Settings are redirected to a temporary folder so the
// real IQmol preferences are not touched.

#include "Util/Preferences.h"
#include "Util/XtbResultDialog.h"
#include <QApplication>
#include <QCheckBox>
#include <QSpinBox>
#include <QComboBox>
#include <QDir>
#include <QLabel>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSettings>
#include <QTemporaryDir>
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

   QTemporaryDir dir;
   QSettings::setDefaultFormat(QSettings::IniFormat);
   QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, dir.path());

   // Stored settings are shown when the window opens
   Preferences::XtbMethod("GFN1-xTB");
   Preferences::XtbSolvent("water");
   Preferences::XtbJob("singlepoint");

   XtbResultDialog* d(XtbResultDialog::instance());
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
   method->setCurrentIndex(method->findData("GFN-FF"));
   job->setCurrentIndex(job->findData("optimize"));
   check(Preferences::XtbSolvent() == "methanol", "solvent is stored when changed");
   check(Preferences::XtbMethod() == "GFN-FF", "method is stored when changed");
   check(Preferences::XtbJob() == "optimize", "job is stored when changed");

   // Charge and spin: the molecule's values by default, an override on request
   QCheckBox* useMolecule(0);
   for (QCheckBox* c : d->findChildren<QCheckBox*>()) useMolecule = c;
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
   d->finished(result, "GFN-FF/ALPB(methanol)", true);
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

   // The values in use are shown, so that unchecking the box starts from them
   d->started("GFN2-xTB", "single point", "radical (4 atoms)", 0, 2);
   check(chargeSpin->value() == 0 && multiplicitySpin->value() == 2, "the boxes show the values in use");
   d->finished(result, "GFN2-xTB", false);
   useMolecule->setChecked(false);
   d->started("GFN2-xTB", "single point", "radical (4 atoms)", 1, 1);
   bool entered(false);
   d->finished(result, "GFN2-xTB", false);
   for (QLabel* l : d->findChildren<QLabel*>()) if (l->text().contains("entered in this window")) entered = true;
   check(entered, "the summary says when the values were entered by hand");
   useMolecule->setChecked(true);

   // A second calculation can follow
   d->started("GFN-FF", "single point", "water (3 atoms)", 0, 1);
   check(!run->isVisible(), "a second calculation can be started");
   IQmol::Calc::Result failed;
   failed.message = "xtb failed";
   d->finished(failed, "GFN-FF", false);
   check(run->isVisible(), "failure returns to the setup state");

   // The same window comes back with the settings after being closed
   d->close();
   check(!d->isVisible(), "the window closes when idle");
   Preferences::XtbSolvent("");
   d->showWindow();
   check(d->isVisible() && solvent->currentData() == "", "reopening shows the current settings");

   std::printf("%d failure(s)\n", s_failures);
   return s_failures == 0 ? 0 : 1;
}
