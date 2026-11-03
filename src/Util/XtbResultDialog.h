#pragma once
/*******************************************************************************

  This file is part of IQmol, a free molecular visualization program. See
  <http://iqmol.org> for more details.

  IQmol is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

********************************************************************************/

#include "Calculator.h"
#include <QDialog>
#include <QPointer>

class QComboBox;
class QGroupBox;
class QLabel;
class QPlainTextEdit;
class QProgressBar;
class QPushButton;

/*!
 *  \file XtbResultDialog.h
 *
 *  \brief The xtb window.  It is where an xtb calculation is set up (job, method,
 *  solvent and backend, stored in the preferences) and started, and where the progress
 *  (live program output) and a summary of the result are shown.
 */

namespace IQmol {

   class XtbResultDialog : public QDialog {

      Q_OBJECT

      public:
         /// A single window is shared by all calculations
         static XtbResultDialog* instance();

         /// Clears the window and shows it in its 'running' state.  The dialog is
         /// modal while the calculation runs so the molecule cannot be edited.
         void started(QString const& method, QString const& what, QString const& subject);

         /// Brings the window to the front without starting anything, so the output of 
         /// the last calculation can be viewed again
         void showWindow();

         void appendOutput(QString const& text);

         /// Switches to the 'finished' state and shows the summary of the result
         void finished(Calc::Result const&, QString const& method, bool optimization);

      Q_SIGNALS:
         /// The Run button was pressed, the settings are in the preferences
         void runRequested();
         void cancelRequested();

      private Q_SLOTS:
         void save();
         void run();
         void saveSettings();

      protected:
         void moveEvent(QMoveEvent*);
         void resizeEvent(QResizeEvent*);
         void reject();
         void done(int result);   // Escape or window close while running means cancel

      private:
         explicit XtbResultDialog(QWidget* parent = 0);

         void setRunning(bool running);
         void loadSettings();
         void placeWindow();
         void rememberGeometry();

         QGroupBox*      m_settings;
         QComboBox*      m_jobCombo;
         QComboBox*      m_methodCombo;
         QComboBox*      m_solventCombo;
         QComboBox*      m_backendCombo;
         QPushButton*    m_runButton;
         QLabel*         m_title;
         QLabel*         m_summary;
         QPlainTextEdit* m_output;
         QProgressBar*   m_progress;
         QPushButton*    m_cancelButton;
         QPushButton*    m_saveButton;
         QPushButton*    m_closeButton;
         QString         m_fileName;   // suggested name for the saved output
         bool            m_running;
         bool            m_started;    // set by started(), to detect that nothing was run
         QPointer<QWidget> m_previous; // window that was active before we opened
         bool            m_placed;     // true once the window has been positioned
   };

} // end namespace IQmol
