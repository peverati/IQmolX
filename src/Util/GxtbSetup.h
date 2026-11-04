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

#include <QString>
#include <QStringList>
#include <functional>

/*!
 *  \file GxtbSetup.h
 *
 *  \brief Fetches the development version of g-xTB (https://github.com/grimme-lab/g-xtb), a
 *  build of the xtb program that has the g-xTB method (xtb --gxtb), into the application data
 *  folder.  The file is pinned to one release of the repository, its SHA-256 sum is checked
 *  before it is unpacked, and nothing is downloaded unless the user asks for it.  g-xTB is not
 *  part of IQmol: its parameters have no license of their own in the repository, and the
 *  program is a preliminary version that changes between releases.
 *
 *  The transfer uses curl and the unpacking tar, which every supported system has.
 */

namespace IQmol {
namespace Calc {

   class GxtbSetup {
      public:
         GxtbSetup();

         /// Where the program is kept: in the application data folder of the user
         static QString defaultFolder();

         /// The xtb executable of the downloaded release, empty if there is none
         static QString defaultExecutable();

         /// The given path if it is an executable, otherwise the downloaded program.  Empty if
         /// there is neither.  (The plain xtb of $XTB_PATH or the PATH is never used for
         /// g-xTB: it does not have the method.)
         static QString findExecutable(QString const& hint = QString());

         /// Which release this is and what it is called on this system.  Empty file name if
         /// there is no download for this system
         static QString release();
         static QString archiveName();
         static QString archiveUrl();
         static QString archiveSha256();
         static qint64  archiveBytes();

         /// For tests: another file to fetch (a file: URL works), with its sum and size
         void setSource(QString const& url, QString const& sha256, QString const& name);

         void setInteractive(bool interactive) { m_interactive = interactive; }
         void setProgressCallback(std::function<void(QString const&)> callback) { m_progress = callback; }
         void cancel() { m_cancelled = true; }
         bool cancelled() const { return m_cancelled; }

         /// Downloads, checks and unpacks the program into defaultFolder(), replacing an 
         /// earlier copy.  False and the reason in 'error' on failure
         bool download(QString& error);

      private:
         bool run(QString const& program, QStringList const& arguments, QString& error, 
            QString const& watchFile = QString(), qint64 expectedBytes = 0);

         QString m_url;
         QString m_sha256;
         QString m_name;
         bool    m_interactive;
         bool    m_cancelled;
         std::function<void(QString const&)> m_progress;
   };

} } // end namespace IQmol::Calc
