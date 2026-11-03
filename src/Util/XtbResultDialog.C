/*******************************************************************************

  This file is part of IQmol, a free molecular visualization program. See
  <http://iqmol.org> for more details.

  IQmol is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

********************************************************************************/

#include "XtbResultDialog.h"
#include "FileDialog.h"
#include "Preferences.h"
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QFormLayout>
#include <QGroupBox>
#include <QDesktopWidget>
#include <QMoveEvent>
#include <QResizeEvent>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QFontDatabase>
#include <QMessageBox>
#include <QTextStream>
#include <QHBoxLayout>
#include <QLabel>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QRegularExpression>
#include <QScrollBar>
#include <QSpinBox>
#include <QVBoxLayout>


namespace IQmol {

XtbResultDialog* XtbResultDialog::instance()
{
   static XtbResultDialog* dialog(new XtbResultDialog());
   return dialog;
}


XtbResultDialog::XtbResultDialog(QWidget* parent) : QDialog(parent), m_charge(0), m_multiplicity(1), m_running(false), m_started(false), m_placed(false)
{
   // An ordinary window (not Qt::Tool): tool windows on macOS cannot be restored
   // from the Dock once they have been minimized.
   setWindowFlags(Qt::Window | Qt::WindowTitleHint | Qt::WindowSystemMenuHint | 
      Qt::WindowMinimizeButtonHint | Qt::WindowCloseButtonHint);
   setWindowTitle("xtb");
   resize(680, 600);

   // The calculation to run.  The choices are kept in the preferences.
   m_settings = new QGroupBox("Calculation", this);
   QFormLayout* form(new QFormLayout(m_settings));

   m_jobCombo = new QComboBox(m_settings);
   m_jobCombo->addItem("Single point energy", "singlepoint");
   m_jobCombo->addItem("Geometry optimization", "optimize");
   form->addRow("Job:", m_jobCombo);

   m_methodCombo = new QComboBox(m_settings);
   m_methodCombo->addItem("GFN2-xTB", "GFN2-xTB");
   m_methodCombo->addItem("GFN1-xTB", "GFN1-xTB");
   m_methodCombo->addItem("GFN-FF (force field)", "GFN-FF");
   form->addRow("Method:", m_methodCombo);

   // The solvents of the ALPB implicit solvation model (xtb name, menu text)
   m_solventCombo = new QComboBox(m_settings);
   m_solventCombo->setToolTip("Implicit solvent (ALPB model)");
   char const* solvents[][2] = {
      { "",             "Gas phase"        }, { "water",        "Water"            },
      { "acetone",      "Acetone"          }, { "acetonitrile", "Acetonitrile"     },
      { "aniline",      "Aniline"          }, { "benzaldehyde", "Benzaldehyde"     },
      { "benzene",      "Benzene"          }, { "ch2cl2",       "Dichloromethane"  },
      { "chcl3",        "Chloroform"       }, { "cs2",          "Carbon disulfide" },
      { "dioxane",      "Dioxane"          }, { "dmf",          "DMF"              },
      { "dmso",         "DMSO"             }, { "ether",        "Diethyl ether"    },
      { "ethylacetate", "Ethyl acetate"    }, { "furane",       "Furan"            },
      { "hexadecane",   "Hexadecane"       }, { "hexane",       "Hexane"           },
      { "methanol",     "Methanol"         }, { "nitromethane", "Nitromethane"     },
      { "octanol",      "Octanol"          }, { "woctanol",     "Octanol (wet)"    },
      { "phenol",       "Phenol"           }, { "thf",          "THF"              },
      { "toluene",      "Toluene"          }
   };
   for (auto const& solvent : solvents) m_solventCombo->addItem(solvent[1], solvent[0]);
   form->addRow("Solvent:", m_solventCombo);

   // By default the charge and spin of the molecule are used.  Entering other values
   // is only meant for the molecule that is open, so they are not stored.
   m_moleculeValues = new QCheckBox("Use the charge and spin of the molecule", m_settings);
   m_moleculeValues->setChecked(true);
   form->addRow(m_moleculeValues);

   m_chargeSpin = new QSpinBox(m_settings);
   m_chargeSpin->setRange(-10, 10);
   m_chargeSpin->setValue(0);
   m_chargeSpin->setEnabled(false);
   form->addRow("Charge:", m_chargeSpin);

   m_multiplicitySpin = new QSpinBox(m_settings);
   m_multiplicitySpin->setRange(1, 11);
   m_multiplicitySpin->setValue(1);
   m_multiplicitySpin->setEnabled(false);
   m_multiplicitySpin->setToolTip("Spin multiplicity 2S+1: 1 singlet, 2 doublet, 3 triplet, ...");
   form->addRow("Spin multiplicity:", m_multiplicitySpin);

   QObject::connect(m_moleculeValues, &QCheckBox::toggled, [this](bool useMolecule) {
      m_chargeSpin->setEnabled(!useMolecule);
      m_multiplicitySpin->setEnabled(!useMolecule);
   });

   m_backendCombo = new QComboBox(m_settings);
   m_backendCombo->addItem("xtb program", "program");
   m_backendCombo->setItemData(0, "Run the xtb program found at the path in the preferences", 
      Qt::ToolTipRole);
#ifdef IQMOL_WITH_LIBXTB
   m_backendCombo->addItem("Built-in library (libxtb)", "library");
   m_backendCombo->setItemData(1, "Run xtb inside IQmol, with IQmol's own geometry optimizer", 
      Qt::ToolTipRole);
   form->addRow("Backend:", m_backendCombo);
#else
   m_backendCombo->hide();   // IQmol was built without libxtb
#endif

   m_title = new QLabel(this);
   QFont font(m_title->font());
   font.setBold(true);
   font.setPointSize(font.pointSize() + 2);
   m_title->setFont(font);

   m_summary = new QLabel(this);
   m_summary->setTextFormat(Qt::RichText);
   m_summary->setTextInteractionFlags(Qt::TextSelectableByMouse);
   m_summary->setWordWrap(true);

   m_progress = new QProgressBar(this);
   m_progress->setRange(0, 0);   // busy indicator
   m_progress->setTextVisible(false);

   m_output = new QPlainTextEdit(this);
   m_output->setReadOnly(true);
   m_output->setLineWrapMode(QPlainTextEdit::NoWrap);
   m_output->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
   m_output->setMaximumBlockCount(20000);

   m_runButton    = new QPushButton("Run", this);
   m_cancelButton = new QPushButton("Cancel", this);
   m_saveButton   = new QPushButton("Save...", this);
   m_closeButton  = new QPushButton("Close", this);
   m_runButton->setDefault(true);

   QHBoxLayout* buttons(new QHBoxLayout());
   buttons->addStretch();
   buttons->addWidget(m_runButton);
   buttons->addWidget(m_cancelButton);
   buttons->addWidget(m_saveButton);
   buttons->addWidget(m_closeButton);

   QVBoxLayout* layout(new QVBoxLayout(this));
   layout->addWidget(m_settings);
   layout->addWidget(m_title);
   layout->addWidget(m_summary);
   layout->addWidget(m_progress);
   layout->addWidget(m_output, 1);
   layout->addLayout(buttons);

   m_title->setText("xtb");
   m_summary->setText("Choose the calculation and press Run to calculate the current molecule.");
   m_progress->hide();
   m_cancelButton->hide();
   loadSettings();

   connect(m_jobCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(saveSettings()));
   connect(m_methodCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(saveSettings()));
   connect(m_solventCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(saveSettings()));
   connect(m_backendCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(saveSettings()));
   connect(m_runButton, SIGNAL(clicked()), this, SLOT(run()));
   connect(m_cancelButton, SIGNAL(clicked()), this, SIGNAL(cancelRequested()));
   connect(m_saveButton, SIGNAL(clicked()), this, SLOT(save()));
   connect(m_closeButton, SIGNAL(clicked()), this, SLOT(close()));
}


void XtbResultDialog::setRunning(bool running)
{
   m_running = running;
   m_settings->setEnabled(!running);
   m_progress->setVisible(running);
   m_runButton->setVisible(!running);
   m_cancelButton->setVisible(running);
   m_saveButton->setVisible(!running);
   m_closeButton->setVisible(!running);
   if (running) m_cancelButton->setFocus(); else m_runButton->setFocus();
}


// Shows the settings stored in the preferences
void XtbResultDialog::loadSettings()
{
   auto select = [](QComboBox* combo, QString const& value) {
      int index(combo->findData(value));
      if (index >= 0) combo->setCurrentIndex(index);
   };

   // do not store what is being loaded
   QList<QComboBox*> combos;
   combos << m_jobCombo << m_methodCombo << m_solventCombo << m_backendCombo;
   for (QComboBox* combo : combos) combo->blockSignals(true);

   select(m_jobCombo, Preferences::XtbJob());
   select(m_methodCombo, Preferences::XtbMethod());
   select(m_solventCombo, Preferences::XtbSolvent());
   select(m_backendCombo, Preferences::XtbBackend());

   for (QComboBox* combo : combos) combo->blockSignals(false);
}


void XtbResultDialog::saveSettings()
{
   Preferences::XtbJob(m_jobCombo->currentData().toString());
   Preferences::XtbMethod(m_methodCombo->currentData().toString());
   Preferences::XtbSolvent(m_solventCombo->currentData().toString());
   Preferences::XtbBackend(m_backendCombo->currentData().toString());
}


void XtbResultDialog::run()
{
   saveSettings();
   m_started = false;
   runRequested();
   if (!m_started) {
      m_title->setText("xtb");
      m_summary->setText("Nothing was calculated.  Is there a molecule with atoms to run xtb on?");
   }
}


bool XtbResultDialog::chargeAndSpin(int& charge, int& multiplicity) const
{
   if (m_moleculeValues->isChecked()) return false;
   charge = m_chargeSpin->value();
   multiplicity = m_multiplicitySpin->value();
   return true;
}


void XtbResultDialog::started(QString const& method, QString const& what,
   QString const& subject, int charge, int multiplicity)
{
   // Show the values in use, so that unchecking the box starts from the molecule's values
   m_charge = charge;
   m_multiplicity = multiplicity;
   if (m_moleculeValues->isChecked()) {
      m_chargeSpin->setValue(charge);
      m_multiplicitySpin->setValue(multiplicity);
   }

   m_started = true;
   m_fileName = method + (what.contains("optim") ? "_opt" : "_sp") + ".out";
   m_fileName.replace(QRegularExpression("[^A-Za-z0-9._-]"), "_");
   m_title->setText(method + " " + what);
   m_summary->setText("Running on " + subject + ", charge " + QString::number(charge) 
      + ", spin multiplicity " + QString::number(multiplicity) + "...");
   m_output->clear();
   setRunning(true);
   // Remember the active window so that it can be reactivated when we close 
   QWidget* active(QApplication::activeWindow());
   if (active && active != this) m_previous = active;

   // Modality must not be changed while the window is visible, otherwise Qt 
   // can leave the main window disabled behind a modal state that no longer exists
   if (isVisible()) hide();
   setWindowState(windowState() & ~Qt::WindowMinimized);   // may be left over from last time
   setModal(true);   // blocks edits to the molecule while the calculation runs
   placeWindow();
   show();
   raise();
   activateWindow();
}


// The first time, the window is restored to where the user last left it, or 
// centred over the IQmol window if there is no saved position.  After that the 
// window stays wherever it was moved to.
void XtbResultDialog::placeWindow()
{
   if (m_placed || isVisible()) return;

   bool restored(false);
   QByteArray geometry(Preferences::XtbWindowGeometry());
   if (!geometry.isEmpty()) restored = restoreGeometry(geometry);

   if (!restored) {
      QWidget* anchor(QApplication::activeWindow());
      if (anchor && anchor != this) {
         QRect r(anchor->frameGeometry());
         move(r.center() - rect().center());
      }
   }
   m_placed = true;
}


void XtbResultDialog::rememberGeometry()
{
   // Ignore the moves and resizes made before the window is shown
   if (m_placed && isVisible()) Preferences::XtbWindowGeometry(saveGeometry());
}


void XtbResultDialog::moveEvent(QMoveEvent* event)
{
   QDialog::moveEvent(event);
   rememberGeometry();
}


void XtbResultDialog::resizeEvent(QResizeEvent* event)
{
   QDialog::resizeEvent(event);
   rememberGeometry();
}


void XtbResultDialog::showWindow()
{
   if (m_running) {
      raise();
      activateWindow();
      return;
   }

   loadSettings();
   setRunning(false);
   placeWindow();
   setWindowState(windowState() & ~Qt::WindowMinimized);
   showNormal();
   raise();
   activateWindow();
}


void XtbResultDialog::appendOutput(QString const& text)
{
   // appendPlainText adds a paragraph break, so insert at the end instead
   QScrollBar* bar(m_output->verticalScrollBar());
   bool atBottom(bar->value() >= bar->maximum() - 2);
   QTextCursor cursor(m_output->textCursor());
   cursor.movePosition(QTextCursor::End);
   cursor.insertText(text);
   if (atBottom) bar->setValue(bar->maximum());
}


void XtbResultDialog::finished(Calc::Result const& result, QString const& method,
   bool optimization)
{
   setRunning(false);
   hide();            // see the note in started()
   setWindowState(windowState() & ~Qt::WindowMinimized);
   setModal(false);
   showNormal();
   raise();
   activateWindow();

   QString what(optimization ? "geometry optimization" : "single point");
   QString s;

   if (result.cancelled) {
      m_title->setText(method + " " + what + " cancelled");
      m_summary->setText("The calculation was cancelled and the structure was not changed.");
      return;
   }

   if (!result.success) {
      m_title->setText(method + " " + what + " failed");
      s = result.message;
      s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\n", "<br>");
      m_summary->setText(s);
      return;
   }

   QString status(" finished");
   if (optimization) status = result.converged ? " converged" : " did not converge";
   m_title->setText(method + " " + what + status);

   double const hartreeToKcal(627.5094740631);

   auto row = [](QString const& name, QString const& value) {
      return "<tr><td>" + name + "&nbsp;&nbsp;</td><td><b>" + value + "</b></td></tr>";
   };

   s = "<table>";
   s += row("Charge, spin multiplicity", QString::number(m_charge) + ", " 
      + QString::number(m_multiplicity) + (m_moleculeValues->isChecked() ? "" 
      : " (entered in this window)"));
   s += row("Energy", QString::number(result.energy, 'f', 6) + " E<sub>h</sub>");
   if (optimization && result.initialEnergy != 0.0) {
      double dE((result.energy - result.initialEnergy) * hartreeToKcal);
      s += row("Energy change", QString::number(dE, 'f', 2) + " kcal/mol");
   }
   if (result.hasGsolv) {
      s += row("Solvation free energy", QString::number(result.gsolv, 'f', 6) + " E<sub>h</sub> ("
         + QString::number(result.gsolv * hartreeToKcal, 'f', 2) + " kcal/mol)");
   }
   if (result.gradientNorm >= 0.0) {
      s += row("Gradient norm", QString::number(result.gradientNorm, 'f', 6) + " E<sub>h</sub>/a<sub>0</sub>");
   }
   if (optimization && result.cycles > 0) {
      s += row("Optimization cycles", QString::number(result.cycles));
   }
   if (result.gapEV >= 0.0) {
      s += row("HOMO-LUMO gap", QString::number(result.gapEV, 'f', 3) + " eV");
   }
   if (result.dipoleDebye >= 0.0) {
      s += row("Dipole moment", QString::number(result.dipoleDebye, 'f', 3) + " D");
   }
   s += row("Wall time", QString::number(result.seconds, 'f', 2) + " s");
   s += "</table>";

   if (!result.converged) {
      s += "<p>The geometry has been updated to the last structure found, "
           "but it is not fully optimized.</p>";
   }
   m_summary->setText(s);
}


void XtbResultDialog::save()
{
   // Remember the folder for the rest of the session
   static QString dir(QDir::homePath());

   QString path(FileDialog::getSaveFileName(this, "Save xtb output", 
      QDir(dir).filePath(m_fileName), "xtb output (*.out *.log *.txt);;All files (*)"));
   if (path.isEmpty()) return;

   QFile file(path);
   if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
      QMessageBox::warning(this, "IQmol", "Unable to write to " + path + "\n" 
         + file.errorString());
      return;
   }

   QTextStream out(&file);
   out << m_output->toPlainText();
   if (!m_output->toPlainText().endsWith('\n')) out << "\n";
   out.flush();

   dir = QFileInfo(path).absolutePath();
}


void XtbResultDialog::done(int result)
{
   QDialog::done(result);
   if (m_previous) {
      m_previous->raise();
      m_previous->activateWindow();
   }
}


void XtbResultDialog::reject()
{
   if (m_running) {
      cancelRequested();
   }else {
      QDialog::reject();
   }
}

} // end namespace IQmol
