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

#include "PythonEngineDialog.h"
#include "Preferences.h"
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFile>
#include <QFileInfo>
#include <QFontDatabase>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QJsonDocument>
#include <QLabel>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QDesktopWidget>
#include <QResizeEvent>
#include <QPushButton>
#include <QRegularExpression>
#include <QScrollBar>
#include <QStandardItemModel>
#include <QSpinBox>
#include <QTextBrowser>
#include <QTextStream>
#include <QVBoxLayout>
#include <limits>


namespace IQmol {

namespace {

QString const s_jobs[][2] = {
   { "singlepoint", "Single point energy" },
   { "optimize",    "Geometry optimization" },
   { "frequencies", "Vibrational frequencies" },
   { "optfreq",     "Optimization and frequencies" },
};

// Text for the settings area.  A word-wrapped label under-reports its height in a form, so the
// text is cut off; a borderless text area that is as high as its text, for the width it will
// have, does not have that problem.
QTextBrowser* infoText(QWidget* parent)
{
   QTextBrowser* text(new QTextBrowser(parent));
   text->setFrameShape(QFrame::NoFrame);
   text->setReadOnly(true);
   text->setOpenLinks(false);
   text->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
   text->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
   text->viewport()->setAutoFillBackground(false);
   return text;
}

void setInfo(QTextBrowser* text, QString const& html, int width)
{
   text->setHtml(html);
   text->document()->setTextWidth(width);
   text->setFixedHeight(static_cast<int>(text->document()->size().height()) + 6);
}


QString plain(QString text)
{
   text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\n", "<br>");
   return text;
}

} // end anonymous namespace


PythonEngineDialog* PythonEngineDialog::instance()
{
   static PythonEngineDialog* dialog(new PythonEngineDialog());
   return dialog;
}


PythonEngineDialog::PythonEngineDialog(QWidget* parent) : QDialog(parent), 
   m_running(false), m_started(false), m_listed(false)
{
   setWindowTitle("Python Engines");
   setMinimumWidth(560);

   m_settings = new QGroupBox("Calculation", this);
   QFormLayout* form(new QFormLayout(m_settings));
   form->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);

   m_pythonEdit = new QLineEdit(this);
   m_pythonEdit->setPlaceholderText("python3 (the first one on the PATH)");
   m_pythonEdit->setToolTip("The Python interpreter, or environment, in which iqmolx-engines is "
      "installed.\nThe engines and their models live there, not in IQmol.");
   QPushButton* browse(new QPushButton("Browse...", this));
   QPushButton* refresh(new QPushButton("Find Engines", this));
   refresh->setToolTip("Start the helper and list the engines it offers");
   QHBoxLayout* pythonRow(new QHBoxLayout());
   pythonRow->setContentsMargins(0, 0, 0, 0);
   pythonRow->addWidget(m_pythonEdit, 1);
   pythonRow->addWidget(browse);
   pythonRow->addWidget(refresh);
   form->addRow("Python:", pythonRow);

   m_environmentInfo = infoText(this);
   m_setupButton = new QPushButton("Set Up Environment", this);
   m_setupButton->setToolTip("Make a Python environment for the engines with uv.  uv downloads a "
      "Python of its own,\nso none has to be installed.  Nothing is installed anywhere else.");
   QHBoxLayout* environmentRow(new QHBoxLayout());
   environmentRow->setContentsMargins(0, 0, 0, 0);
   environmentRow->addWidget(m_environmentInfo, 1);
   environmentRow->addWidget(m_setupButton);
   form->addRow(QString(), environmentRow);

   m_engineCombo = new QComboBox(this);
   m_installButton = new QPushButton("Install", this);
   m_installButton->hide();
   QHBoxLayout* engineRow(new QHBoxLayout());
   engineRow->setContentsMargins(0, 0, 0, 0);
   engineRow->addWidget(m_engineCombo, 1);
   engineRow->addWidget(m_installButton);
   form->addRow("Engine:", engineRow);

   m_engineInfo = infoText(this);
   form->addRow(QString(), m_engineInfo);

   m_optionsBox = new QGroupBox("Engine settings", this);
   m_optionsForm = new QFormLayout(m_optionsBox);
   form->addRow(m_optionsBox);

   m_jobCombo = new QComboBox(this);
   for (auto const& job : s_jobs) m_jobCombo->addItem(job[1], job[0]);
   form->addRow("Job:", m_jobCombo);

   m_moleculeValues = new QCheckBox("Use the charge and spin of the molecule", this);
   m_moleculeValues->setChecked(true);
   form->addRow(QString(), m_moleculeValues);
   m_chargeSpin = new QSpinBox(this);
   m_chargeSpin->setRange(-10, 10);
   m_multiplicitySpin = new QSpinBox(this);
   m_multiplicitySpin->setRange(1, 10);
   QHBoxLayout* spinRow(new QHBoxLayout());
   spinRow->setContentsMargins(0, 0, 0, 0);
   spinRow->addWidget(new QLabel("Charge", this));
   spinRow->addWidget(m_chargeSpin);
   spinRow->addWidget(new QLabel("Spin multiplicity", this));
   spinRow->addWidget(m_multiplicitySpin);
   spinRow->addStretch();
   form->addRow(QString(), spinRow);
   m_chargeSpin->setEnabled(false);
   m_multiplicitySpin->setEnabled(false);
   connect(m_moleculeValues, &QCheckBox::toggled, [this](bool on) {
      m_chargeSpin->setEnabled(!on);
      m_multiplicitySpin->setEnabled(!on);
   });

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
   m_output->setMinimumHeight(130);        // room to read what the engine and uv write; more when the window is larger

   m_runButton    = new QPushButton("Run", this);
   m_cancelButton = new QPushButton("Cancel", this);
   m_saveButton   = new QPushButton("Save Output...", this);
   m_closeButton  = new QPushButton("Close", this);
   m_runButton->setDefault(true);

   m_runHint = new QLabel(this);
   m_runHint->setWordWrap(true);

   QHBoxLayout* buttons(new QHBoxLayout());
   buttons->addWidget(m_runHint, 1);
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

   m_title->setText("Python engines");
   m_summary->setText("Choose an engine and press Run to calculate the current molecule.");
   m_progress->hide();
   m_cancelButton->hide();

   m_pythonEdit->setText(Preferences::PythonEnginePath());
   int job(m_jobCombo->findData(Preferences::PythonEngineJob()));
   if (job >= 0) m_jobCombo->setCurrentIndex(job);

   connect(browse, SIGNAL(clicked()), this, SLOT(browseForPython()));
   connect(refresh, SIGNAL(clicked()), this, SLOT(refreshEngines()));
   connect(m_setupButton, SIGNAL(clicked()), this, SLOT(setUpEnvironment()));
   connect(m_installButton, SIGNAL(clicked()), this, SLOT(installPackages()));
   connect(m_pythonEdit, SIGNAL(editingFinished()), this, SLOT(saveSettings()));
   connect(m_pythonEdit, &QLineEdit::editingFinished, [this]() { updateEnvironmentRow(); engineChanged(); });
   connect(m_engineCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(engineChanged()));
   connect(m_jobCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(saveSettings()));
   connect(m_runButton, SIGNAL(clicked()), this, SLOT(run()));
   connect(m_cancelButton, SIGNAL(clicked()), this, SIGNAL(cancelRequested()));
   connect(m_saveButton, SIGNAL(clicked()), this, SLOT(save()));
   connect(m_closeButton, SIGNAL(clicked()), this, SLOT(close()));
   // the helper and its model go when IQmol is closed
   connect(qApp, &QCoreApplication::aboutToQuit, [this]() { m_engine.shutdown(); });

   updateEnvironmentRow();
   showEngine();

   QSize size(Preferences::PythonEngineWindowSize());
   QRect screen(QApplication::desktop()->availableGeometry(this));
   resize(qMax(size.width(), 600), qMin(qMax(size.height(), 600), screen.height() - 40));
}


void PythonEngineDialog::resizeEvent(QResizeEvent* event)
{
   QDialog::resizeEvent(event);
   if (isVisible()) Preferences::PythonEngineWindowSize(size());   // kept for the next session
}


QString PythonEngineDialog::job() const
{
   return m_jobCombo->currentData().toString();
}


QString PythonEngineDialog::engineTitle() const
{
   QString title(m_descriptor.value("title").toString());
   return title.isEmpty() ? m_engineCombo->currentData().toString() : title;
}


bool PythonEngineDialog::chargeAndSpin(int& charge, int& multiplicity) const
{
   if (m_moleculeValues->isChecked()) return false;
   charge = m_chargeSpin->value();
   multiplicity = m_multiplicitySpin->value();
   return true;
}


// ----- engines -----

void PythonEngineDialog::browseForPython()
{
   QString start(m_pythonEdit->text().isEmpty() ? QDir::homePath() : m_pythonEdit->text());
   QString path(QFileDialog::getOpenFileName(this, "Python interpreter", start));
   if (path.isEmpty()) return;
   m_pythonEdit->setText(path);
   saveSettings();
   refreshEngines();
}


/// Starts the helper (again, if the interpreter changed) and lists its engines
void PythonEngineDialog::refreshEngines()
{
   m_listed = true;
   saveSettings();
   updateEnvironmentRow();
   QString wanted(Preferences::PythonEngineName());

   m_engine.setPython(effectivePython());
   m_engine.setPythonPath(qEnvironmentVariable("IQMOLX_ENGINES_PATH"));
   m_engine.setEngine(wanted);

   QApplication::setOverrideCursor(Qt::WaitCursor);
   QString error;
   m_descriptors = m_engine.listEngines(&error);
   QApplication::restoreOverrideCursor();

   m_engineCombo->blockSignals(true);
   m_engineCombo->clear();
   int select(-1);
   for (QJsonValue const& value : m_descriptors) {
       QJsonObject d(value.toObject());
       QString title(d.value("title").toString(d.value("name").toString()));
       bool available(d.value("available").toBool(true));
       m_engineCombo->addItem(available ? title : title + " (not available)", 
          d.value("name").toString());
       int index(m_engineCombo->count() - 1);
       // An engine that is not available can be selected: the window then says what is
       // missing and offers to install it.  An available one is preferred as the default.
       if (!available) m_engineCombo->setItemData(index, d.value("unavailable_reason").toString(), 
          Qt::ToolTipRole);
       if (d.value("name").toString() == wanted || (select < 0 && available)) select = index;
   }
   if (select >= 0) m_engineCombo->setCurrentIndex(select);
   m_engineCombo->blockSignals(false);

   if (m_descriptors.isEmpty()) {
      m_summary->setText(plain(error.isEmpty() ? QString("The Python helper offers no engines.") : error));
      m_title->setText("The Python helper is not available");
   }
   engineChanged();
}


void PythonEngineDialog::engineChanged()
{
   QString name(m_engineCombo->currentData().toString());
   m_descriptor = QJsonObject();
   for (QJsonValue const& value : m_descriptors) {
       if (value.toObject().value("name").toString() == name) m_descriptor = value.toObject();
   }
   if (!name.isEmpty()) Preferences::PythonEngineName(name);
   showEngine();
}


QString PythonEngineDialog::describeElements(QJsonValue const& value)
{
   if (!value.isArray()) return "all elements";
   QStringList numbers;
   for (QJsonValue const& z : value.toArray()) numbers << QString::number(z.toInt());
   return "atomic numbers " + numbers.join(", ");
}


/// Shows what is known of the selected engine: its description, license and settings
void PythonEngineDialog::showEngine()
{
   bool const known(!m_descriptor.isEmpty());
   bool const available(known && m_descriptor.value("available").toBool(true));
   m_engineInfo->setVisible(known);
   buildOptions(m_descriptor);
   updateRunButton();

   QStringList requirements;
   for (QJsonValue const& r : m_descriptor.value("requirements").toArray()) requirements << r.toString();
   m_installButton->setVisible(known && !available && !requirements.isEmpty());
   m_installButton->setText("Install " + requirements.join(", "));
   QString python(effectivePython());
   bool const canInstall(!python.isEmpty() && QFile::exists(QFileInfo(python).absolutePath() + "/../pyvenv.cfg")
      && !Calc::PythonSetup::findUv().isEmpty());
   m_installButton->setEnabled(canInstall);
   m_installButton->setToolTip(canInstall 
      ? "Install " + requirements.join(", ") + " into the Python environment with uv"
      : "Packages can be installed into an environment made with uv (Set Up Environment).\n"
        "Otherwise install " + requirements.join(", ") + " in the Python you chose, or install uv.");
   if (!known) return;

   QString s("<small>" + plain(m_descriptor.value("description").toString()));
   if (!available) {
      s = "<small><b>Not available:</b> " + plain(m_descriptor.value("unavailable_reason").toString()) 
        + "<br>" + plain(m_descriptor.value("description").toString());
   }
   QStringList properties;
   for (QJsonValue const& p : m_descriptor.value("properties").toArray()) properties << p.toString();
   s += "<br>Elements: " + describeElements(m_descriptor.value("elements"));
   s += ". Provides: " + properties.join(", ") + ".";
   QString license(m_descriptor.value("license").toString());
   if (!license.isEmpty()) s += "<br>License: " + plain(license) + ".";
   if (m_descriptor.value("requires_download").toBool()) {
      s += "<br>The model is downloaded the first time it is used.";
   }
   if (m_descriptor.value("requires_login").toBool()) {
      s += "<br>The model needs a login (for example to Hugging Face), which you do in "
           "a terminal, not here.";
   }
   QString citation(m_descriptor.value("citation").toString());
   if (!citation.isEmpty()) s += "<br>Please cite: " + plain(citation);
   setInfo(m_engineInfo, s + "</small>", 415);

   // what the engine cannot handle is not offered
   bool charge(m_descriptor.value("charge").toBool(true));
   bool spin(m_descriptor.value("multiplicity").toBool(true));
   m_moleculeValues->setToolTip(QString(charge ? "" : "This engine is for neutral molecules.  ") 
      + (spin ? "" : "This engine is for closed-shell molecules (multiplicity 1)."));
}


QJsonObject PythonEngineDialog::storedOptions(QString const& engine) const
{
   QJsonDocument doc(QJsonDocument::fromJson(Preferences::PythonEngineOptions().toUtf8()));
   return doc.object().value(engine).toObject();
}


/// Makes one input widget for each option of the engine, as its descriptor says
void PythonEngineDialog::buildOptions(QJsonObject const& descriptor)
{
   while (m_optionsForm->rowCount() > 0) m_optionsForm->removeRow(0);
   m_options.clear();

   QJsonObject stored(storedOptions(descriptor.value("name").toString()));
   double const big(1e9);

   for (QJsonValue const& value : descriptor.value("options").toArray()) {
       QJsonObject spec(value.toObject());
       QString name(spec.value("name").toString()), type(spec.value("type").toString());
       QJsonValue current(stored.contains(name) ? stored.value(name) : spec.value("default"));
       QWidget* widget(0);

       if (type == "bool") {
          QCheckBox* box(new QCheckBox(this));
          box->setChecked(current.toBool());
          connect(box, SIGNAL(toggled(bool)), this, SLOT(saveSettings()));
          widget = box;
       }else if (type == "int") {
          QSpinBox* box(new QSpinBox(this));
          box->setRange(static_cast<int>(spec.value("min").toDouble(-big)),
                        static_cast<int>(spec.value("max").toDouble(big)));
          box->setValue(current.toInt());
          connect(box, SIGNAL(valueChanged(int)), this, SLOT(saveSettings()));
          widget = box;
       }else if (type == "float") {
          QDoubleSpinBox* box(new QDoubleSpinBox(this));
          box->setDecimals(4);
          box->setRange(spec.value("min").toDouble(-big), spec.value("max").toDouble(big));
          box->setSingleStep(0.1);
          box->setValue(current.toDouble());
          connect(box, SIGNAL(valueChanged(double)), this, SLOT(saveSettings()));
          widget = box;
       }else if (type == "choice") {
          QComboBox* box(new QComboBox(this));
          for (QJsonValue const& c : spec.value("choices").toArray()) box->addItem(c.toString());
          box->setCurrentText(current.toString());
          connect(box, SIGNAL(currentIndexChanged(int)), this, SLOT(saveSettings()));
          widget = box;
       }else {
          QLineEdit* edit(new QLineEdit(this));
          edit->setText(current.toString());
          edit->setMinimumWidth(300);
          edit->setPlaceholderText(spec.value("placeholder").toString());
          connect(edit, SIGNAL(editingFinished()), this, SLOT(saveSettings()));
          connect(edit, SIGNAL(textChanged(QString const&)), this, SLOT(updateRunButton()));
          widget = edit;
       }

       widget->setToolTip(spec.value("description").toString());
       m_optionsForm->addRow(spec.value("title").toString(name) + ":", widget);
       OptionWidget entry;
       entry.spec = spec;
       entry.widget = widget;
       m_options.append(entry);
   }
   m_optionsBox->setVisible(!m_options.isEmpty());
}


/// The values in the option widgets that belong to 'scope' ("load" or "compute")
QJsonObject PythonEngineDialog::optionValues(QString const& scope) const
{
   QJsonObject values;
   for (OptionWidget const& o : m_options) {
       if (o.spec.value("scope").toString("load") != scope) continue;
       QString name(o.spec.value("name").toString()), type(o.spec.value("type").toString());
       if (type == "bool") {
          values[name] = static_cast<QCheckBox*>(o.widget)->isChecked();
       }else if (type == "int") {
          values[name] = static_cast<QSpinBox*>(o.widget)->value();
       }else if (type == "float") {
          values[name] = static_cast<QDoubleSpinBox*>(o.widget)->value();
       }else if (type == "choice") {
          values[name] = static_cast<QComboBox*>(o.widget)->currentText();
       }else {
          values[name] = static_cast<QLineEdit*>(o.widget)->text();
       }
   }
   return values;
}


void PythonEngineDialog::saveSettings()
{
   Preferences::PythonEnginePath(m_pythonEdit->text().trimmed());
   Preferences::PythonEngineJob(job());

   // the settings of the engines, per engine
   QString name(m_descriptor.value("name").toString());
   if (name.isEmpty() || m_options.isEmpty()) return;
   QJsonDocument doc(QJsonDocument::fromJson(Preferences::PythonEngineOptions().toUtf8()));
   QJsonObject all(doc.object()), mine(optionValues("load"));
   QJsonObject compute(optionValues("compute"));
   for (QString const& key : compute.keys()) mine[key] = compute.value(key);
   all[name] = mine;
   Preferences::PythonEngineOptions(QString::fromUtf8(QJsonDocument(all).toJson(QJsonDocument::Compact)));
}


bool PythonEngineDialog::prepare(QString& problem)
{
   if (!m_listed || m_descriptor.isEmpty()) refreshEngines();
   if (m_descriptor.isEmpty()) {
      problem = m_descriptors.isEmpty() ? "The Python helper could not be started.  Choose the "
         "Python interpreter in the Python Engines window."
         : "No Python engine is selected, or the selected one is not available.";
      return false;
   }
   saveSettings();
   m_engine.setPython(effectivePython());
   m_engine.setEngine(m_descriptor.value("name").toString());
   m_engine.setLoadOptions(optionValues("load"));
   m_engine.setComputeOptions(optionValues("compute"));
   return true;
}


/// Run is on when the engine can run and every required text field is filled in; when it is off
/// the window says why, next to the button
void PythonEngineDialog::updateRunButton()
{
   bool const known(!m_descriptor.isEmpty());
   bool const available(known && m_descriptor.value("available").toBool(true));

   QStringList missing;
   for (OptionWidget const& o : m_options) {
       if (!o.spec.value("required").toBool() || o.spec.value("type").toString() != "string") continue;
       if (static_cast<QLineEdit*>(o.widget)->text().trimmed().isEmpty()) {
          missing << o.spec.value("title").toString(o.spec.value("name").toString());
       }
   }

   m_runButton->setEnabled(available && missing.isEmpty());
   if (!known) {
      m_runHint->clear();
   }else if (!available) {
      m_runHint->setText("<small>This engine is not available: see above.</small>");
   }else if (!missing.isEmpty()) {
      m_runHint->setText("<small>Fill in: <b>" + plain(missing.join(", ")) + "</b></small>");
   }else {
      m_runHint->clear();
   }
}


// ----- the Python environment -----

/// The interpreter that runs the helper: the one in the field, else the IQMolX environment if 
/// there is one, else (empty) python3 from the PATH
QString PythonEngineDialog::effectivePython() const
{
   QString python(m_pythonEdit->text().trimmed());
   if (!python.isEmpty()) return python;
   return Calc::PythonSetup::pythonIn(Calc::PythonSetup::defaultEnvironment());
}


void PythonEngineDialog::updateEnvironmentRow()
{
   QString environment(Calc::PythonSetup::defaultEnvironment());
   bool const have(!Calc::PythonSetup::pythonIn(environment).isEmpty());
   bool const uv(!Calc::PythonSetup::findUv().isEmpty());

   m_pythonEdit->setPlaceholderText(have ? "the IQMolX environment" 
      : "python3 (the first one on the PATH)");
   m_setupButton->setVisible(!have);
   m_setupButton->setEnabled(uv);
   if (have) {
      setInfo(m_environmentInfo, "<small>IQMolX Python environment: " + plain(environment) + "</small>", 415);
   }else if (uv) {
      setInfo(m_environmentInfo, "<small>There is no IQMolX Python environment yet.  Set Up Environment "
         "makes one with uv: it downloads a Python of its own (about 25 MB) and installs the "
         "helper.  Or choose a Python that has iqmolx-engines.</small>", 255);
   }else {
      setInfo(m_environmentInfo, "<small>There is no IQMolX Python environment, and uv, which makes "
         "it, was not found.  Install uv (for example <tt>brew install uv</tt>) or choose a Python "
         "that has iqmolx-engines.</small>", 255);
   }
}


/// Runs a job of PythonSetup with its output in the window and Cancel working.  The helper is
/// stopped first, as the job changes what it can import.
bool PythonEngineDialog::runSetup(QString const& title, 
   std::function<bool(Calc::PythonSetup&, QString&)> job)
{
   m_engine.shutdown();
   m_title->setText(title);
   m_summary->setText("This can take a minute.  The output of uv is shown below.");
   m_output->clear();
   setRunning(true);
   QApplication::processEvents();

   Calc::PythonSetup setup;
   setup.setInteractive(true);
   setup.setProgressCallback([this](QString const& text) { appendOutput(text); });
   QMetaObject::Connection cancelConnection(connect(this, &PythonEngineDialog::cancelRequested, 
      [&setup]() { setup.cancel(); }));

   QString error;
   bool ok(job(setup, error));
   disconnect(cancelConnection);
   setRunning(false);

   if (ok) {
      m_title->setText(title + " finished");
      m_summary->setText("Done.");
   }else if (setup.cancelled()) {
      m_title->setText(title + " cancelled");
      m_summary->setText("Cancelled.  What was already downloaded is kept.");
   }else {
      m_title->setText(title + " failed");
      m_summary->setText(plain(error));
   }
   return ok;
}


void PythonEngineDialog::setUpEnvironment()
{
   QString environment(Calc::PythonSetup::defaultEnvironment());
   bool ok(runSetup("Python environment", [environment](Calc::PythonSetup& setup, QString& error) {
      return setup.createEnvironment(environment, error);
   }));
   updateEnvironmentRow();
   if (!ok) return;

   // The environment that was just made is what the window uses from now on
   m_pythonEdit->setText(Calc::PythonSetup::pythonIn(environment));
   saveSettings();
   updateEnvironmentRow();
   refreshEngines();
}


void PythonEngineDialog::installPackages()
{
   QStringList requirements;
   for (QJsonValue const& r : m_descriptor.value("requirements").toArray()) requirements << r.toString();
   QString python(effectivePython());
   if (requirements.isEmpty() || python.isEmpty()) return;

   QString title("Installing " + requirements.join(", "));
   bool ok(runSetup(title, [python, requirements](Calc::PythonSetup& setup, QString& error) {
      return setup.install(python, requirements, error);
   }));
   if (ok) refreshEngines();
}


// ----- running -----

void PythonEngineDialog::setRunning(bool running)
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


void PythonEngineDialog::run()
{
   saveSettings();
   m_started = false;
   runRequested();
   if (!m_started) {
      m_title->setText("Python engines");
      m_summary->setText("Nothing was calculated.  Is there a molecule with atoms to run on?");
   }
}


void PythonEngineDialog::started(QString const& method, QString const& what,
   QString const& subject, int charge, int multiplicity)
{
   if (m_moleculeValues->isChecked()) {
      m_chargeSpin->setValue(charge);
      m_multiplicitySpin->setValue(multiplicity);
   }

   m_started = true;
   m_fileName = method + (what.contains("optim") ? "_opt" : "_sp") + ".out";
   m_fileName.replace(QRegularExpression("[^A-Za-z0-9._-]"), "_");
   m_title->setText(method + " " + what);
   m_summary->setText("Running on " + plain(subject) + ", charge " + QString::number(charge) 
      + ", spin multiplicity " + QString::number(multiplicity) + "...");
   m_output->clear();
   setRunning(true);
   QWidget* active(QApplication::activeWindow());
   if (active && active != this) m_previous = active;

   // Modality must not be changed while the window is visible (see the xtb window)
   if (isVisible()) hide();
   setWindowState(windowState() & ~Qt::WindowMinimized);
   setModal(true);   // blocks edits to the molecule while the calculation runs
   if (active && active != this) move(active->frameGeometry().center() - rect().center());
   show();
   raise();
   activateWindow();
}


void PythonEngineDialog::showWindow()
{
   if (m_running) {
      raise();
      activateWindow();
      return;
   }
   setRunning(false);
   updateEnvironmentRow();
   setWindowState(windowState() & ~Qt::WindowMinimized);
   showNormal();
   raise();
   activateWindow();
   if (!m_listed) refreshEngines();
}


void PythonEngineDialog::appendOutput(QString const& text)
{
   // appendPlainText adds a paragraph break, so insert at the end instead
   QScrollBar* bar(m_output->verticalScrollBar());
   bool atBottom(bar->value() >= bar->maximum() - 2);
   QTextCursor cursor(m_output->textCursor());
   cursor.movePosition(QTextCursor::End);
   cursor.insertText(text);
   if (atBottom) bar->setValue(bar->maximum());
}


void PythonEngineDialog::finished(Calc::Result const& result, QString const& method,
   QString const& what, bool optimization)
{
   setRunning(false);
   hide();            // see the note in started()
   setWindowState(windowState() & ~Qt::WindowMinimized);
   setModal(false);
   showNormal();
   raise();
   activateWindow();

   if (result.cancelled) {
      m_title->setText(method + " " + what + " cancelled");
      m_summary->setText("The calculation was cancelled and the structure was not changed.");
      return;
   }

   if (!result.success) {
      m_title->setText(method + " " + what + " failed");
      m_summary->setText(plain(result.message));
      appendOutput("\nFailed: " + result.message + "\n");   // not only above, where it is easy to miss
      return;
   }

   QString status(" finished");
   if (optimization) status = result.converged ? " converged" : " did not converge";
   if (optimization && result.converged && result.hasVibrations) status = " finished";
   m_title->setText(method + " " + what + status);

   double const eVPerHartree(27.211386245988);
   double const kcalPerHartree(627.5094740631);

   auto row = [](QString const& name, QString const& value) {
      return "<tr><td>" + name + "&nbsp;&nbsp;</td><td><b>" + value + "</b></td></tr>";
   };

   QString s("<table>");
   s += row("Energy", QString::number(result.energy, 'f', 8) + " Eh &nbsp;("
      + QString::number(result.energy * eVPerHartree, 'f', 4) + " eV)");
   if (optimization) {
      s += row("Energy change", QString::number((result.energy - result.initialEnergy) * kcalPerHartree, 'f', 2)
         + " kcal/mol in " + QString::number(result.cycles) + " cycles");
   }
   if (result.gradientNorm >= 0.0) s += row("Gradient norm", QString::number(result.gradientNorm, 'f', 6) + " Eh/Bohr");
   if (result.dipoleDebye >= 0.0) s += row("Dipole moment", QString::number(result.dipoleDebye, 'f', 3) + " D");
   if (result.hasVibrations) {
      int imaginary(result.vibrations.imaginary());
      s += row("Frequencies", QString::number(result.vibrations.frequencies.size()) + " modes, "
         + (imaginary == 0 ? QString("no imaginary frequencies") 
                           : QString::number(imaginary) + " imaginary"));
      s += row("Zero point energy", QString::number(result.vibrations.zeroPointEnergy, 'f', 6) + " Eh");
   }
   s += row("Time", QString::number(result.seconds, 'f', 1) + " s");
   s += "</table>";
   if (!result.message.isEmpty()) s += "<br>" + plain(result.message);
   m_summary->setText(s);
}


void PythonEngineDialog::save()
{
   QString path(QFileDialog::getSaveFileName(this, "Save output", m_fileName));
   if (path.isEmpty()) return;
   QFile file(path);
   if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
      QTextStream(&file) << m_output->toPlainText();
   }
}


void PythonEngineDialog::done(int result)
{
   QDialog::done(result);
   if (m_previous) {
      m_previous->raise();
      m_previous->activateWindow();
   }
}


void PythonEngineDialog::reject()
{
   if (m_running) {
      cancelRequested();
   }else {
      QDialog::reject();
   }
}

} // end namespace IQmol
