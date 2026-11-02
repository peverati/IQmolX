/*******************************************************************************

  This file is part of IQmol, a free molecular visualization program. See
  <http://iqmol.org> for more details.

  IQmol is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

********************************************************************************/

#include "XtbResultDialog.h"
#include <QFontDatabase>
#include <QHBoxLayout>
#include <QLabel>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QScrollBar>
#include <QVBoxLayout>


namespace IQmol {

XtbResultDialog* XtbResultDialog::instance()
{
   static XtbResultDialog* dialog(new XtbResultDialog());
   return dialog;
}


XtbResultDialog::XtbResultDialog(QWidget* parent) : QDialog(parent), m_running(false)
{
   // Qt::Tool keeps the window above the main window on macOS
   setWindowFlags(Qt::Tool | Qt::WindowTitleHint | Qt::WindowCloseButtonHint);
   setWindowTitle("xtb");
   resize(640, 480);

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

   m_cancelButton = new QPushButton("Cancel", this);
   m_closeButton  = new QPushButton("Close", this);
   m_closeButton->setDefault(true);

   QHBoxLayout* buttons(new QHBoxLayout());
   buttons->addStretch();
   buttons->addWidget(m_cancelButton);
   buttons->addWidget(m_closeButton);

   QVBoxLayout* layout(new QVBoxLayout(this));
   layout->addWidget(m_title);
   layout->addWidget(m_summary);
   layout->addWidget(m_progress);
   layout->addWidget(m_output, 1);
   layout->addLayout(buttons);

   connect(m_cancelButton, SIGNAL(clicked()), this, SIGNAL(cancelRequested()));
   connect(m_closeButton, SIGNAL(clicked()), this, SLOT(close()));
}


void XtbResultDialog::setRunning(bool running)
{
   m_running = running;
   m_progress->setVisible(running);
   m_cancelButton->setVisible(running);
   m_closeButton->setVisible(!running);
   if (running) m_cancelButton->setFocus(); else m_closeButton->setFocus();
}


void XtbResultDialog::started(QString const& method, QString const& what,
   QString const& subject)
{
   m_title->setText(method + " " + what);
   m_summary->setText("Running on " + subject + "...");
   m_output->clear();
   setRunning(true);
   setModal(true);   // blocks edits to the molecule while the calculation runs
   show();
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
   setModal(false);
   show();      // setModal hides the window
   raise();

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
   s += row("Energy", QString::number(result.energy, 'f', 6) + " E<sub>h</sub>");
   if (optimization && result.initialEnergy != 0.0) {
      double dE((result.energy - result.initialEnergy) * hartreeToKcal);
      s += row("Energy change", QString::number(dE, 'f', 2) + " kcal/mol");
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


void XtbResultDialog::reject()
{
   if (m_running) {
      cancelRequested();
   }else {
      QDialog::reject();
   }
}

} // end namespace IQmol
