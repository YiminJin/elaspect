/*
  Copyright (C) 2023 by Yimin Jin.

  This file is part of elASPECT.

  elASPECT is modified from the free software ASPECT; you can 
  redistribute it and/or modify it under the terms of the GNU 
  General Public License as published by the Free Software 
  Foundation; either version 2, or (at your option) any later 
  version.

  elASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with elASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <elaspect/global.h>
#include <elaspect/quadrature_point_data.h>

namespace elaspect
{
  template <int dim>
  QPDHandler<dim>::QPDHandler(const Triangulation<dim> &tria_,
                              const unsigned int n_quadrature_points,
                              const unsigned int n_components)
    : tria(&tria_, typeid(*this).name())
    , n_q_points(n_quadrature_points)
    , n_comp(n_components)
  {}


  template <int dim>
  void QPDHandler<dim>::reinit()
  {
    // resize the data vector
    unsigned int n_locally_owned_active_cells = 0;
    for (const auto &cell : tria->active_cell_iterators())
      if (cell->is_locally_owned())
        ++n_locally_owned_active_cells;

    data.resize(n_locally_owned_active_cells * n_q_points * n_comp);
    std::fill(data.begin(), data.end(), 0);

    // map the active cell index to data address
    index_data_map.resize(tria->n_active_cells());
    std::fill(index_data_map.begin(), index_data_map.end(), nullptr);
    auto p = data.begin();
    for (const auto &cell : tria->active_cell_iterators())
      if (cell->is_locally_owned())
      {
        index_data_map[cell->active_cell_index()] = &(*p);
        std::advance(p, n_q_points * n_comp);
      }

    Assert(p == data.end(), ExcInternalError());
  }


  template <int dim>
  typename QPDHandler<dim>::cell_iterator
  QPDHandler<dim>::begin(const unsigned int level) const
  {
    typename Triangulation<dim>::cell_iterator cell =
      this->get_triangulation().begin(level);
    if (cell == this->get_triangulation().end(level))
      return end(level);
    return cell_iterator(*cell, this);
  }


  template <int dim>
  typename QPDHandler<dim>::active_cell_iterator
  QPDHandler<dim>::begin_active(const unsigned int level) const
  {
    cell_iterator i = begin(level);
    if (i.state() != IteratorState::valid)
      return i;
    while (i->has_children())
      if ((++i).state() != IteratorState::valid)
        return i;
    return i;
  }


  template <int dim>
  typename QPDHandler<dim>::cell_iterator
  QPDHandler<dim>::end() const
  {
    return cell_iterator(&this->get_triangulation(), -1, -1, this);
  }


  template <int dim>
  typename QPDHandler<dim>::cell_iterator
  QPDHandler<dim>::end(const unsigned int level) const
  {
    typename Triangulation<dim>::cell_iterator cell =
      this->get_triangulation().end(level);
    if (cell.state() != IteratorState::valid)
      return end();
    return cell_iterator(*cell, this);
  }


  template <int dim>
  typename QPDHandler<dim>::active_cell_iterator
  QPDHandler<dim>::end_active(const unsigned int level) const
  {
    typename Triangulation<dim>::cell_iterator cell =
      this->get_triangulation().end_active(level);
    if (cell.state() != IteratorState::valid)
      return active_cell_iterator(end());
    return active_cell_iterator(*cell, this);
  }


  template <int dim>
  IteratorRange<typename QPDHandler<dim>::cell_iterator>
  QPDHandler<dim>::cell_iterators() const
  {
    return IteratorRange<typename QPDHandler<dim>::cell_iterator>(
      begin(), end());
  }


  template <int dim>
  IteratorRange<typename QPDHandler<dim>::active_cell_iterator>
  QPDHandler<dim>::active_cell_iterators() const
  {
    return IteratorRange<typename QPDHandler<dim>::active_cell_iterator>(
      begin_active(), end());
  }
}


// explicit instantiations
namespace elaspect
{
#define INSTANTIATE(dim) \
  template class QPDHandler<dim>; \
  template class QPDAccessor<dim>; 

  ELASPECT_INSTANTIATE(INSTANTIATE)
}
