#pragma once

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

#include "PythonEngine.h"
#include "PythonSetup.h"
#include <QDialog>
#include <QJsonObject>
#include <QPointer>
#include <QVector>
#include <functional>

class QCheckBox;
class QComboBox;
class QFormLayout;
class QGroupBox;
class QLabel;
class QLineEdit;
class QPlainTextEdit;
class QProgressBar;
class QPushButton;
class QTextBrowser;
class QSpinBox;

/*!
 *  \file PythonEngineDialog.h
 *
 *  \brief The window of the Python engines (iqmolx-engines helper).  It is where the
 *  Python interpreter and the engine are chosen, where the settings of the engine are made,
 *  and where a calculation is started and followed.  The list of engines and the settings 
 *  are not known to IQmol: they are built from the capabilities the helper reports, so a 
 *  new engine needs no change here.
 */

namespace IQmol {

   class PythonEngineDialog : public QDialog {

      Q_OBJECT

      public:
         /// A single window is shared by all calculations
         static PythonEngineDialog* instance();

         /// The calculator behind the window.  It keeps the helper process (and with it the
         /// loaded model) alive between calculations.  Call prepare() before using it.
         Calc::PythonEngine& engine() { return m_engine; }

         /// What the window runs: "singlepoint", "optimize", "frequencies" or "optfreq"
         QString job() const;

         /// Applies the settings of the window to engine().  Returns false, and says why in 
         /// 'problem', if no engine can be used.
         bool prepare(QString& problem);

         /// Name of the selected engine as shown in titles, e.g. "Demo (Morse toy potential)"
         QString engineTitle() const;

         /// True if the charge and spin entered in this window are to be used instead of the
         /// values of the molecule
         bool chargeAndSpin(int& charge, int& multiplicity) const;

         /// Clears the window and shows it in its 'running' state (modal, so the molecule
         /// cannot be edited meanwhile)
         void started(QString const& method, QString const& what, QString const& subject,
            int charge, int multiplicity);

         void appendOutput(QString const& text);

         /// Switches to the 'finished' state and shows the summary of the result
         void finished(Calc::Result const&, QString const& method, QString const& what,
            bool optimization);

      public Q_SLOTS:
         /// Opens the window (or brings it to the front) without starting anything
         void showWindow();

         /// Asks the helper for its engines and shows them
         void refreshEngines();

      Q_SIGNALS:
         /// The Run button was pressed, the settings are in the preferences
         void runRequested();
         void cancelRequested();

      private Q_SLOTS:
         void run();
         void saveSettings();
         void engineChanged();
         void browseForPython();
         void setUpEnvironment();
         void installPackages();
         void updateRunButton();
         void save();

      protected:
         void reject();
         void resizeEvent(QResizeEvent*);
         void done(int result);

      private:
         struct OptionWidget {
            QJsonObject spec;
            QWidget*    widget;
         };

         explicit PythonEngineDialog(QWidget* parent = 0);

         void setRunning(bool running);
         void fitToScreen();
         void sizeTo(int width, int height);
         QString effectivePython() const;
         void updateEnvironmentRow();
         bool runSetup(QString const& title, 
            std::function<bool(Calc::PythonSetup&, QString&)> job);
         void showEngine();
         void buildOptions(QJsonObject const& descriptor);
         QJsonObject optionValues(QString const& scope) const;
         QJsonObject storedOptions(QString const& engine) const;
         static QString describeElements(QJsonValue const&);

         Calc::PythonEngine m_engine;
         QJsonArray     m_descriptors;     // from the helper
         QJsonObject    m_descriptor;      // of the selected engine
         QVector<OptionWidget> m_options;

         QGroupBox*     m_settings;
         QLineEdit*     m_pythonEdit;
         QTextBrowser*  m_environmentInfo;
         QPushButton*   m_setupButton;
         QPushButton*   m_installButton;
         QComboBox*     m_engineCombo;
         QTextBrowser*  m_engineInfo;
         QGroupBox*     m_optionsBox;
         QFormLayout*   m_optionsForm;
         QComboBox*     m_jobCombo;
         QCheckBox*     m_moleculeValues;
         QSpinBox*      m_chargeSpin;
         QSpinBox*      m_multiplicitySpin;
         QLabel*        m_title;
         QLabel*        m_summary;
         QProgressBar*  m_progress;
         QPlainTextEdit* m_output;
         QLabel*        m_runHint;       // why Run is off, when it is
         QPushButton*   m_runButton;
         QPushButton*   m_cancelButton;
         QPushButton*   m_saveButton;
         QPushButton*   m_closeButton;
         QString        m_fileName;
         bool           m_running;
         bool           m_started;       // set by started(), to detect that nothing was run
         bool           m_listed;        // the engines were requested at least once
         bool           m_userSized;     // the user dragged the window to a size of their own
         QSize          m_ownSize;       // the size the window gave itself last
         QPointer<QWidget> m_previous;   // window that was active before we opened
   };

} // end namespace IQmol
