/*******************************************************************************

  Copyright (C) 2011-2015 Andrew Gilbert

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

#include "ViewerModelView.h"
#include "QVariantPtr.h"
#include "Layer/MoleculeLayer.h"
#include "Layer/OrbitalsLayer.h"
#include "Layer/SurfaceLayer.h"

#include <QContextMenuEvent>
#include <QModelIndex>
#include <QMenu>
#include <QMouseEvent>
#include <QStyleOptionViewItem>

#include <QDebug>


namespace IQmol {

ViewerModelView::ViewerModelView(QWidget* parent) : QTreeView(parent)
{ 
   setDragDropMode(QAbstractItemView::DropOnly); 
   setDropIndicatorShown(true);
   setIndentation(15);
   setSelectionMode(QAbstractItemView::ExtendedSelection);
   setEditTriggers( QAbstractItemView::EditKeyPressed);
   setExpandsOnDoubleClick(false);

   // The native check boxes of the items are not drawn in the dark tree on macOS, so there
   // would be no way to see, or change, which layers are shown.
   setStyleSheet(
      "QTreeView::indicator { width: 14px; height: 14px; }"
      "QTreeView::indicator:unchecked { image: url(:/icons/CheckUnchecked.png); }"
      "QTreeView::indicator:checked { image: url(:/icons/CheckChecked.png); }"
      "QTreeView::indicator:indeterminate { image: url(:/icons/CheckPartial.png); }");
}


// Clicking one of the surfaces of an orbital layer (not its check box) shows that surface on
// its own: with several orbitals calculated they would otherwise all be drawn on top of
// each other.  The other surfaces can be shown again with their check boxes.
void ViewerModelView::mousePressEvent(QMouseEvent* event)
{
   QTreeView::mousePressEvent(event);
   if (event->button() != Qt::LeftButton || event->modifiers() != Qt::NoModifier) return;

   QModelIndex index(indexAt(event->pos()));
   if (!index.isValid()) return;

   QStyleOptionViewItem option;
   option.initFrom(this);
   option.rect = visualRect(index);
   option.features |= QStyleOptionViewItem::HasCheckIndicator;
   if (style()->subElementRect(QStyle::SE_ItemViewItemCheckIndicator, &option, this)
      .adjusted(-2, -2, 2, 2).contains(event->pos())) return;

   Layer::Base* base(QVariantPtr<Layer::Base>::toPointer(index.data(Qt::UserRole+1)));
   Layer::Surface* surface(qobject_cast<Layer::Surface*>(base));
   if (!surface) return;

   QStandardItem* parent(surface->QStandardItem::parent());
   if (!parent || !qobject_cast<Layer::Orbitals*>(
      QVariantPtr<Layer::Base>::toPointer(parent->data()))) return;

   for (int row = 0; row < parent->rowCount(); ++row) {
       QStandardItem* sibling(parent->child(row));
       if (sibling != surface && qobject_cast<Layer::Surface*>(
          QVariantPtr<Layer::Base>::toPointer(sibling->data(Qt::UserRole+1)))) {
          sibling->setCheckState(Qt::Unchecked);
       }
   }
   surface->setCheckState(Qt::Checked);
}


void ViewerModelView::contextMenuEvent(QContextMenuEvent*)
{
   QItemSelectionModel* selectModel(selectionModel());
   QModelIndexList selection(selectModel->selectedRows());

   Layer::Base* base;
   QMenu menu(this);

   if (selection.size() == 1) {
      QModelIndex target(selection.first());
      base = QVariantPtr<Layer::Base>::toPointer(target.data(Qt::UserRole+1));
	  // QVariantPtr is not guaranteed to return a valid Layer::Base pointer.
      if ( (base = qobject_cast<Layer::Base*>(base)) ) {
         QList<QAction*> actions(base->getActions());
         for (auto action : actions) menu.addAction(action);
      }
   }else {
      // Custom actions when more than one Layer is selected
      bool allMolecules(true);
      for (auto item : selection) {
          base = QVariantPtr<Layer::Base>::toPointer(item.data(Qt::UserRole+1));
          if (//!qobject_cast<Layer::Group*>(base) && 
              !qobject_cast<Layer::Molecule*>(base) ) {
              allMolecules = false;
              break;
          }
      }
      if (allMolecules) {
         QAction* action(new QAction("Merge Molecules")); 
         connect(action, SIGNAL(triggered()), this, SLOT(mergeSelection()));
         menu.addAction(action);

         action = new QAction("Show Molecules"); 
         connect(action, SIGNAL(triggered()), this, SLOT(showMolecules()));
         menu.addAction(action);

         action = new QAction("Hide Molecules"); 
         connect(action, SIGNAL(triggered()), this, SLOT(hideMolecules()));
         menu.addAction(action);
      }
   }

   bool canPromote(true);
   for (auto item : selection) {
       base = QVariantPtr<Layer::Base>::toPointer(item.data(Qt::UserRole+1));
       if (!qobject_cast<Layer::Group*>(base) && 
           !qobject_cast<Layer::Atom*>(base)  && 
           !qobject_cast<Layer::Bond*>(base)  &&
           !qobject_cast<Layer::Molecule*>(base) ) {
           canPromote = false;
           break;
       }
   }

   if (canPromote) {
      QAction* action(new QAction("New Molecule From Selection")); 
      connect(action, SIGNAL(triggered()), this, SLOT(newMoleculeFromSelection()));
      menu.addAction(action);
   }

   if (!menu.isEmpty()) menu.exec(QCursor::pos());
}


void ViewerModelView::newMoleculeFromSelection()
{
   QItemSelectionModel* selectModel(selectionModel());
   QModelIndexList selection(selectModel->selectedRows());

   newMoleculeFromSelection(selection);
}


void ViewerModelView::mergeSelection()
{
   QItemSelectionModel* selectModel(selectionModel());
   QModelIndexList selection(selectModel->selectedRows());

   mergeSelection(selection);
}


void ViewerModelView::deleteSelection()
{
   QItemSelectionModel* selectModel(selectionModel());
   QModelIndexList selection(selectModel->selectedRows());

   deleteSelection(selection);
}


void ViewerModelView::showMolecules()
{
   QItemSelectionModel* selectModel(selectionModel());
   QModelIndexList selection(selectModel->selectedRows());

   showMolecules(selection);
}


void ViewerModelView::hideMolecules()
{
   QItemSelectionModel* selectModel(selectionModel());
   QModelIndexList selection(selectModel->selectedRows());

   hideMolecules(selection);
}

} // end namespace IQmol
