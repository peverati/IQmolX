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

class QLabel;
class QPlainTextEdit;
class QProgressBar;
class QPushButton;

/*!
 *  \file XtbResultDialog.h
 *
 *  \brief Window that shows the progress of an xtb calculation (live program
 *  output) and a summary of the result once it has finished.
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

         void appendOutput(QString const& text);

         /// Switches to the 'finished' state and shows the summary of the result
         void finished(Calc::Result const&, QString const& method, bool optimization);

      Q_SIGNALS:
         void cancelRequested();

      protected:
         void reject();   // Escape or window close while running means cancel

      private:
         explicit XtbResultDialog(QWidget* parent = 0);

         void setRunning(bool running);

         QLabel*         m_title;
         QLabel*         m_summary;
         QPlainTextEdit* m_output;
         QProgressBar*   m_progress;
         QPushButton*    m_cancelButton;
         QPushButton*    m_closeButton;
         bool            m_running;
   };

} // end namespace IQmol
