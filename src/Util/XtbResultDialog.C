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

#include "XtbResultDialog.h"
#include "FileDialog.h"
#include "GxtbSetup.h"
#include "Preferences.h"
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QGroupBox>
#include <QDesktopWidget>
#include <QMoveEvent>
#include <QResizeEvent>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QFileDialog>
#include <QLineEdit>
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


XtbResultDialog::XtbResultDialog(QWidget* parent) : QDialog(parent), m_charge(0), m_multiplicity(1), m_temperature(300.0), m_running(false), m_started(false), m_placed(false)
{
   // An ordinary window (not Qt::Tool): tool windows on macOS cannot be restored
   // from the Dock once they have been minimized.
   setWindowFlags(Qt::Window | Qt::WindowTitleHint | Qt::WindowSystemMenuHint | 
      Qt::WindowMinimizeButtonHint | Qt::WindowCloseButtonHint);
   setWindowTitle("xtb");
   resize(680, 600);

   // The calculation to run.  The choices are kept in the preferences.
   m_settings = new QGroupBox("Calculation", this);
   m_settings->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Maximum);   // never squeezed
   QFormLayout* form(new QFormLayout(m_settings));

   m_jobCombo = new QComboBox(m_settings);
   m_jobCombo->addItem("Single point energy", "singlepoint");
   m_jobCombo->addItem("Geometry optimization", "optimize");
   m_jobCombo->addItem("Vibrational frequencies", "frequencies");
   m_jobCombo->addItem("Geometry optimization and frequencies", "optfreq");
   form->addRow("Job:", m_jobCombo);

   m_methodCombo = new QComboBox(m_settings);
   m_methodCombo->addItem("GFN2-xTB", "GFN2-xTB");
   m_methodCombo->addItem("GFN1-xTB", "GFN1-xTB");
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
#endif
   // g-xTB is a development version of the xtb program with one more method; it is a backend of 
   // its own so that the xtb program above is never replaced by it
   m_backendCombo->addItem("g-xTB (development version)", "gxtb");
   m_backendCombo->setItemData(m_backendCombo->count() - 1, "Run the preliminary g-xTB method with the "
      "g-xTB build of the xtb program,\nwhich is downloaded from the g-xTB repository on request", 
      Qt::ToolTipRole);
   form->addRow("Backend:", m_backendCombo);

   // Settings of the xtb program, only shown when it is the selected backend
   auto fileRow = [this](QLineEdit*& edit, char const* slot, QString const& tip) {
      QWidget* row(new QWidget(m_settings));
      QHBoxLayout* rowLayout(new QHBoxLayout(row));
      rowLayout->setContentsMargins(0, 0, 0, 0);
      edit = new QLineEdit(row);
      edit->setToolTip(tip);
      QPushButton* browse(new QPushButton("Browse...", row));
      rowLayout->addWidget(edit, 1);
      rowLayout->addWidget(browse);
      connect(browse, SIGNAL(clicked()), this, slot);
      return row;
   };
   m_programRow = fileRow(m_programEdit, SLOT(browseForProgram()),
      "The xtb executable.  If empty, the copy inside the application is used if there is one,\n"
      "otherwise xtb is looked for in $XTB_PATH and the PATH.");
   m_programEdit->setPlaceholderText("path of the xtb executable");
   form->addRow("xtb program:", m_programRow);

   // The g-xTB build of the xtb program: where it is, a button to fetch it, and its solvation model
   m_gxtbRow = new QWidget(m_settings);
   QHBoxLayout* gxtbLayout(new QHBoxLayout(m_gxtbRow));
   gxtbLayout->setContentsMargins(0, 0, 0, 0);
   m_gxtbEdit = new QLineEdit(m_gxtbRow);
   m_gxtbEdit->setPlaceholderText("empty: the downloaded g-xTB");
   m_gxtbEdit->setToolTip("The xtb executable of a g-xTB build.  If empty, the copy that Download fetched is used.");
   QPushButton* gxtbBrowse(new QPushButton("Browse...", m_gxtbRow));
   m_gxtbDownloadButton = new QPushButton("Download", m_gxtbRow);
   gxtbLayout->addWidget(m_gxtbEdit, 1);
   gxtbLayout->addWidget(gxtbBrowse);
   gxtbLayout->addWidget(m_gxtbDownloadButton);
   connect(gxtbBrowse, SIGNAL(clicked()), this, SLOT(browseForGxtb()));
   connect(m_gxtbDownloadButton, SIGNAL(clicked()), this, SLOT(downloadGxtb()));
   form->addRow("g-xTB program:", m_gxtbRow);

   m_gxtbModelCombo = new QComboBox(m_settings);
   m_gxtbModelCombo->addItem("GBE (generalized Born)", "gbe");
   m_gxtbModelCombo->addItem("COSMO (ddCOSMO)", "cosmo");
   m_gxtbModelCombo->setToolTip("The implicit solvation model that is used with g-xTB when a solvent is chosen.\n"
      "Both are electrostatic only; a solvation model made for g-xTB is still being developed.\n"
      "GBE can be unstable in geometry optimizations, and the COSMO gradient is not consistent.");
   form->addRow("Solvation model:", m_gxtbModelCombo);

   m_gxtbNote = new QLabel(m_settings);
   m_gxtbNote->setWordWrap(true);
   m_gxtbNote->setTextFormat(Qt::RichText);
   m_gxtbNote->setTextInteractionFlags(Qt::TextBrowserInteraction);
   m_gxtbNote->setOpenExternalLinks(true);
   m_gxtbNote->setStyleSheet("QLabel { background: #fff3cd; color: #4d3b00; border: 1px solid #d9b94b; padding: 5px; }");
   m_gxtbNote->setText("<b>g-xTB is a preliminary development version</b> (" + Calc::GxtbSetup::release() 
      + ").  Results change between releases and differ from the preprint; not every xtb feature "
      "works with it.  Please cite T. Froitzheim, M. M&uuml;ller, A. Hansen, S. Grimme, ChemRxiv 2025, "
      "<a href=\"https://doi.org/10.26434/chemrxiv-2025-bjxvt\">doi 10.26434/chemrxiv-2025-bjxvt</a>.  "
      "Code: <a href=\"https://github.com/grimme-lab/g-xtb\">github.com/grimme-lab/g-xtb</a>.");
   form->addRow(m_gxtbNote);

   m_optionsRow = new QWidget(m_settings);
   QHBoxLayout* optionsLayout(new QHBoxLayout(m_optionsRow));
   optionsLayout->setContentsMargins(0, 0, 0, 0);
   m_optionsEdit = new QLineEdit(m_optionsRow);
   m_optionsEdit->setPlaceholderText("e.g.  --etemp 500 --verbose");
   m_optionsEdit->setToolTip("Further command line options for xtb.  Quote values that contain spaces.\n"
      "The method, charge (--chrg), spin (--uhf), solvent and electronic temperature (--etemp)\n"
      "given here replace those chosen in this window.");
   optionsLayout->addWidget(m_optionsEdit);
   form->addRow("xtb options:", m_optionsRow);

   m_workDirRow = fileRow(m_workDirEdit, SLOT(browseForWorkDir()),
      "xtb runs in this directory and the input and output files are left there.\n"
      "Leave it empty to use a temporary directory that is removed afterwards.");
   m_workDirEdit->setPlaceholderText("empty: temporary directory");
   form->addRow("Working directory:", m_workDirRow);
   // Fermi smearing of the orbital occupations.  It is a setting of the method, not of the
   // molecule, so it is stored in the preferences.
   m_temperatureSpin = new QDoubleSpinBox(m_settings);
   m_temperatureSpin->setRange(0.0, 20000.0);
   m_temperatureSpin->setDecimals(0);
   m_temperatureSpin->setSingleStep(50.0);
   m_temperatureSpin->setSuffix(" K");
   m_temperatureSpin->setValue(300.0);
   m_temperatureSpin->setToolTip("Electronic temperature of the Fermi smearing of the orbital occupations.\n"
      "xtb uses 300 K.  Below a few K, molecules with a (nearly) degenerate HOMO and LUMO cannot be calculated.");
   form->addRow("Electronic temperature:", m_temperatureSpin);

   m_orbitalsCheck = new QCheckBox("Calculate the molecular orbitals", m_settings);
   m_orbitalsCheck->setToolTip("Adds the orbitals to the molecule, where they can be plotted.\n"
      "Only with the xtb program (not the built-in library) and not for GFN-FF.\n"
      "Switch this off for very large molecules.");
   form->addRow("Orbitals:", m_orbitalsCheck);
   for (QWidget* row : QList<QWidget*>({ m_programRow, m_gxtbRow, m_optionsRow, m_workDirRow })) {
       row->setMinimumWidth(440);   // paths are long; the form would keep the fields small
   }

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
   m_saveButton   = new QPushButton("Save Output...", this);
   m_saveButton->setToolTip("Save the text of the xtb output to a file");
   m_saveAllButton = new QPushButton("Save All Files...", this);
   m_saveAllButton->setToolTip("Copy all the files that the xtb program wrote (input, output, geometry,\n"
      "molden orbitals, ...) to a folder");
   m_saveAllButton->setEnabled(false);
   m_closeButton  = new QPushButton("Close", this);
   m_runButton->setDefault(true);

   QHBoxLayout* buttons(new QHBoxLayout());
   buttons->addStretch();
   buttons->addWidget(m_runButton);
   buttons->addWidget(m_cancelButton);
   buttons->addWidget(m_saveButton);
   buttons->addWidget(m_saveAllButton);
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

   // Each session starts from the same defaults: a geometry optimization with GFN2-xTB in the
   // gas phase, using the built-in library if there is one.  The choices made during the
   // session are kept until IQmol is closed.
   Preferences::XtbJob("optimize");
   Preferences::XtbMethod("GFN2-xTB");
   Preferences::XtbSolvent(QString());
   Preferences::XtbElectronicTemperature(300.0);
#ifdef IQMOL_WITH_LIBXTB
   Preferences::XtbBackend("library");
#else
   Preferences::XtbBackend("program");
#endif
   loadSettings();

   connect(m_jobCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(saveSettings()));
   connect(m_methodCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(saveSettings()));
   connect(m_solventCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(saveSettings()));
   connect(m_backendCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(saveSettings()));
   connect(m_temperatureSpin, SIGNAL(valueChanged(double)), this, SLOT(saveSettings()));
   connect(m_orbitalsCheck, SIGNAL(toggled(bool)), this, SLOT(saveSettings()));
   connect(m_programEdit, SIGNAL(editingFinished()), this, SLOT(saveSettings()));
   connect(m_gxtbEdit, SIGNAL(editingFinished()), this, SLOT(saveSettings()));
   connect(m_gxtbModelCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(saveSettings()));
   connect(m_optionsEdit, SIGNAL(editingFinished()), this, SLOT(saveSettings()));
   connect(m_workDirEdit, SIGNAL(editingFinished()), this, SLOT(saveSettings()));
   connect(m_backendCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updateProgramFields()));
   updateProgramFields();
   connect(m_runButton, SIGNAL(clicked()), this, SLOT(run()));
   connect(m_cancelButton, SIGNAL(clicked()), this, SIGNAL(cancelRequested()));
   connect(m_saveButton, SIGNAL(clicked()), this, SLOT(save()));
   connect(m_saveAllButton, SIGNAL(clicked()), this, SLOT(saveAll()));
   // temporary directories of the xtb program are not left behind
   connect(qApp, &QCoreApplication::aboutToQuit, this, &XtbResultDialog::forgetTemporaryDirectory);
   // the electronic temperature goes back to the xtb default when IQmol is closed
   connect(qApp, &QCoreApplication::aboutToQuit, []() { Preferences::XtbElectronicTemperature(300.0); });
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
   m_saveAllButton->setVisible(!running);
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
   combos << m_jobCombo << m_methodCombo << m_solventCombo << m_backendCombo << m_gxtbModelCombo;
   for (QComboBox* combo : combos) combo->blockSignals(true);
   m_temperatureSpin->blockSignals(true);
   m_orbitalsCheck->blockSignals(true);

   select(m_jobCombo, Preferences::XtbJob());
   select(m_methodCombo, Preferences::XtbMethod());
   select(m_solventCombo, Preferences::XtbSolvent());
   select(m_backendCombo, Preferences::XtbBackend());
   m_programEdit->setText(Preferences::XtbPath());
   select(m_gxtbModelCombo, Preferences::XtbGxtbSolvation());
   m_gxtbEdit->setText(Preferences::XtbGxtbPath());
   m_optionsEdit->setText(Preferences::XtbOptions());
   m_workDirEdit->setText(Preferences::XtbWorkDir());
   m_temperatureSpin->setValue(Preferences::XtbElectronicTemperature());
   m_orbitalsCheck->setChecked(Preferences::XtbOrbitals());

   for (QComboBox* combo : combos) combo->blockSignals(false);
   m_temperatureSpin->blockSignals(false);
   m_orbitalsCheck->blockSignals(false);
   updateProgramFields();   // the backend may have changed without a signal
}


void XtbResultDialog::saveSettings()
{
   Preferences::XtbJob(m_jobCombo->currentData().toString());
   Preferences::XtbMethod(m_methodCombo->currentData().toString());
   Preferences::XtbSolvent(m_solventCombo->currentData().toString());
   Preferences::XtbBackend(m_backendCombo->currentData().toString());
   Preferences::XtbElectronicTemperature(m_temperatureSpin->value());
   Preferences::XtbOrbitals(m_orbitalsCheck->isChecked());
   Preferences::XtbPath(m_programEdit->text().trimmed());
   Preferences::XtbGxtbPath(m_gxtbEdit->text().trimmed());
   Preferences::XtbGxtbSolvation(m_gxtbModelCombo->currentData().toString());
   Preferences::XtbOptions(m_optionsEdit->text().trimmed());
   Preferences::XtbWorkDir(m_workDirEdit->text().trimmed());
}


// The program, its options and its directory only matter when the program is the backend
void XtbResultDialog::updateProgramFields()
{
   QString backend(m_backendCombo->currentData().toString());
   bool const gxtb(backend == "gxtb");
   bool const program(backend == "program");
   QFormLayout* form(qobject_cast<QFormLayout*>(m_settings->layout()));
   if (!form) return;
   auto show = [form](QWidget* row, bool visible) {
      row->setVisible(visible);
      if (QWidget* label = form->labelForField(row)) label->setVisible(visible);
   };
   // the xtb program and g-xTB both run a program, with its options, directory and orbitals
   show(m_programRow, program);
   show(m_temperatureSpin, program);
   show(m_optionsRow, program || gxtb);
   show(m_workDirRow, program || gxtb);
   show(m_orbitalsCheck, program || gxtb);
   // g-xTB is its own method, with its own program and solvation models
   show(m_gxtbRow, gxtb);
   show(m_gxtbModelCombo, gxtb);
   m_gxtbNote->setVisible(gxtb);
   show(m_methodCombo, !gxtb);
   m_solventCombo->setToolTip(gxtb ? "Implicit solvent (GBE or COSMO model of g-xTB)" : "Implicit solvent (ALPB model)");
   m_gxtbDownloadButton->setToolTip("Downloads " + Calc::GxtbSetup::release() + " from github.com/grimme-lab/g-xtb ("
      + (Calc::GxtbSetup::archiveBytes() > 0 ? "about " + QString::number(Calc::GxtbSetup::archiveBytes() / 1000000) + " MB" : QString("a few tens of MB"))
      + ") and checks its SHA-256 sum.\nIt is kept in " + Calc::GxtbSetup::defaultFolder() + ".  Nothing is downloaded unless you press this.");
   m_gxtbDownloadButton->setEnabled(!Calc::GxtbSetup::archiveName().isEmpty());
   m_gxtbDownloadButton->setText(Calc::GxtbSetup::defaultExecutable().isEmpty() ? "Download" : "Download again");
}


void XtbResultDialog::browseForProgram()
{
   QString start(m_programEdit->text().isEmpty() ? QDir::homePath() : m_programEdit->text());
   QString path(FileDialog::getOpenFileName(this, "Select the xtb executable", start));
   if (!path.isEmpty()) {
      m_programEdit->setText(path);
      saveSettings();
   }
}


void XtbResultDialog::browseForGxtb()
{
   QString start(m_gxtbEdit->text().isEmpty() ? QDir::homePath() : m_gxtbEdit->text());
   QString path(FileDialog::getOpenFileName(this, "Select the xtb executable of a g-xTB build", start));
   if (!path.isEmpty()) {
      m_gxtbEdit->setText(path);
      saveSettings();
   }
}


/// Fetches the g-xTB build of the xtb program.  The output of the transfer is shown in the output
/// area and Cancel works, like for a calculation.
void XtbResultDialog::downloadGxtb()
{
   m_title->setText("Downloading g-xTB");
   m_summary->setText("From github.com/grimme-lab/g-xtb, " + Calc::GxtbSetup::release() + ".  This can take a minute.");
   m_output->clear();
   setRunning(true);
   QApplication::processEvents();

   Calc::GxtbSetup setup;
   setup.setInteractive(true);
   setup.setProgressCallback([this](QString const& text) { appendOutput(text); });
   QMetaObject::Connection cancelConnection(connect(this, &XtbResultDialog::cancelRequested, 
      [&setup]() { setup.cancel(); }));
   QString error;
   bool ok(setup.download(error));
   disconnect(cancelConnection);
   setRunning(false);

   if (ok) {
      m_title->setText("g-xTB downloaded");
      m_summary->setText("The g-xTB program is in " + Calc::GxtbSetup::defaultFolder().toHtmlEscaped() + ".");
      m_gxtbEdit->clear();      // use the downloaded copy
      saveSettings();
   }else if (setup.cancelled()) {
      m_title->setText("Download cancelled");
      m_summary->setText("Nothing was installed.");
   }else {
      m_title->setText("Download failed");
      m_summary->setText(error.toHtmlEscaped().replace("\n", "<br>"));
   }
   updateProgramFields();
}


void XtbResultDialog::browseForWorkDir()
{
   QString start(m_workDirEdit->text().isEmpty() ? QDir::homePath() : m_workDirEdit->text());
   QString path(FileDialog::getExistingDirectory(this, "Select the working directory", start));
   if (!path.isEmpty()) {
      m_workDirEdit->setText(path);
      saveSettings();
   }
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
   // the electronic temperature is only used by the xtb program
   QString backend(m_backendCombo->currentData().toString());
   m_temperature = (backend == "library" || backend == "gxtb") ? -1.0 : m_temperatureSpin->value();
   if (m_moleculeValues->isChecked()) {
      m_chargeSpin->setValue(charge);
      m_multiplicitySpin->setValue(multiplicity);
   }

   m_started = true;
   m_fileName = method + (what.contains("optim") ? "_opt" : "_sp") + ".out";
   m_fileName.replace(QRegularExpression("[^A-Za-z0-9._-]"), "_");
   m_title->setText(method + " " + what);
   m_summary->setText("Running on " + subject + ", charge " + QString::number(charge) 
      + ", spin multiplicity " + QString::number(multiplicity)
      + (m_temperature >= 0.0 ? ", electronic temperature " + QString::number(m_temperature, 'f', 0) + " K" : QString())
      + "...");
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


// Removes the temporary directory of an earlier run, which is kept for saving
void XtbResultDialog::forgetTemporaryDirectory()
{
   if (!m_temporaryDirectory.isEmpty()) QDir(m_temporaryDirectory).removeRecursively();
   m_temporaryDirectory.clear();
   m_directory.clear();
}


void XtbResultDialog::saveAll()
{
   static QString dir(QDir::homePath());

   if (m_directory.isEmpty() || !QDir(m_directory).exists()) {
      QMessageBox::information(this, "IQmol", "There are no xtb files to save.");
      return;
   }

   QString target(FileDialog::getExistingDirectory(this, 
      "Save all the xtb files to this folder (it can be created here)", dir));
   if (target.isEmpty()) return;
   dir = target;

   QDir source(m_directory);
   if (QFileInfo(target).canonicalFilePath() == QFileInfo(m_directory).canonicalFilePath()) {
      QMessageBox::information(this, "IQmol", "That is the folder the files are in already.");
      return;
   }

   QStringList names(source.entryList(QDir::Files | QDir::NoDotAndDotDot));
   QStringList existing;
   for (QString const& name : names) if (QFileInfo::exists(QDir(target).filePath(name))) existing << name;
   if (!existing.isEmpty()) {
      QString msg("These files exist in " + target + " and will be replaced:\n\n" 
         + existing.mid(0, 8).join("\n") + (existing.size() > 8 ? "\n..." : "") + "\n\nContinue?");
      if (QMessageBox::question(this, "IQmol", msg) != QMessageBox::Yes) return;
   }

   int copied(0);
   QStringList failed;
   for (QString const& name : names) {
       QString to(QDir(target).filePath(name));
       QFile::remove(to);
       if (QFile::copy(source.filePath(name), to)) ++copied; else failed << name;
   }

   QString msg("Saved " + QString::number(copied) + " files to " + target + ".");
   if (!failed.isEmpty()) msg += "\n\nCould not copy: " + failed.join(", ");
   QMessageBox::information(this, "IQmol", msg);
}


void XtbResultDialog::finished(Calc::Result const& result, QString const& method,
   QString const& what, bool optimization)
{
   // The files of the program, if it left any, can be saved; those of an earlier run go
   if (m_temporaryDirectory != result.directory) forgetTemporaryDirectory();
   m_directory = result.directory;
   if (result.temporaryDirectory) m_temporaryDirectory = result.directory;
   m_saveAllButton->setEnabled(!m_directory.isEmpty());

   setRunning(false);
   hide();            // see the note in started()
   setWindowState(windowState() & ~Qt::WindowMinimized);
   setModal(false);
   showNormal();
   raise();
   activateWindow();

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
   if (optimization && result.converged && result.hasVibrations) status = " finished";
   m_title->setText(method + " " + what + status);

   double const hartreeToKcal(627.5094740631);

   auto row = [](QString const& name, QString const& value) {
      return "<tr><td>" + name + "&nbsp;&nbsp;</td><td><b>" + value + "</b></td></tr>";
   };

   s = "<table>";
   s += row("Charge, spin multiplicity", QString::number(m_charge) + ", " 
      + QString::number(m_multiplicity) + (m_moleculeValues->isChecked() ? "" 
      : " (entered in this window)"));
   if (m_temperature >= 0.0) s += row("Electronic temperature", QString::number(m_temperature, 'f', 0) + " K");
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
   if (result.hasVibrations) {
      Calc::VibrationalData const& v(result.vibrations);
      int imaginary(v.imaginary());
      int real(v.frequencies.size() - imaginary);
      s += row("Vibrational modes", QString::number(v.frequencies.size()));
      if (real > 0) {
         s += row("Frequencies", QString::number(v.frequencies[imaginary], 'f', 1) + " to " 
            + QString::number(v.frequencies.last(), 'f', 1) + " cm<sup>-1</sup>");
      }
      s += row("Zero point energy", QString::number(v.zeroPointEnergy, 'f', 6) + " E<sub>h</sub> ("
         + QString::number(v.zeroPointEnergy * hartreeToKcal, 'f', 2) + " kcal/mol)");
      if (!v.intensities.isEmpty()) {
         int strongest(0);
         for (int i = 1; i < v.intensities.size(); ++i) {
             if (v.intensities[i] > v.intensities[strongest]) strongest = i;
         }
         s += row("Strongest IR band", QString::number(v.frequencies[strongest], 'f', 1) + " cm<sup>-1</sup> ("
            + QString::number(v.intensities[strongest], 'f', 1) + " km/mol)");
      }
      QString imaginaryText(QString::number(imaginary));
      if (imaginary > 0) {
         imaginaryText += " (";
         for (int i = 0; i < imaginary && i < 3; ++i) {
             imaginaryText += (i ? ", " : "") + QString::number(v.frequencies[i], 'f', 1);
         }
         imaginaryText += (imaginary > 3 ? ", ..." : "") + QString(" cm<sup>-1</sup>)");
      }
      s += row("Imaginary frequencies", imaginary > 0 ? "<span style=\"color:#c03030\">" + imaginaryText + "</span>" : imaginaryText);
   }
   s += row("Wall time", QString::number(result.seconds, 'f', 2) + " s");
   s += "</table>";

   if (result.hasVibrations) {
      int imaginary(result.vibrations.imaginary());
      if (imaginary > 0) {
         s += "<p>The structure is not a minimum: " + QString::number(imaginary) 
            + (imaginary == 1 ? " imaginary frequency indicates" : " imaginary frequencies indicate")
            + " a saddle point (or an incompletely optimized structure).</p>";
      }
      if (!optimization && result.gradientNorm > 0.002) {
         s += "<p><b>The geometry is not a stationary point</b> (gradient norm " 
            + QString::number(result.gradientNorm, 'f', 4) + " E<sub>h</sub>/a<sub>0</sub>), "
              "so the frequencies are not reliable.  Use <i>Geometry optimization and frequencies</i>.</p>";
      }
      s += "<p>The modes are in the <i>Frequencies</i> entry of the molecule: select a mode there to "
           "animate it, or open it to see the IR spectrum.</p>";
   }
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
