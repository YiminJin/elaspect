#ifndef _elaspect_quadrature_point_data_h
#define _elaspect_quadrature_point_data_h

#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/iterator_range.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria.h>

namespace elaspect
{
  using namespace dealii;

  template <int dim> class QPDAccessor;

template <int dim>
  class QPDHandler
  {
    public:
      using cell_iterator        = TriaIterator<QPDAccessor<dim>>;
      using active_cell_iterator = TriaActiveIterator<QPDAccessor<dim>>;

      QPDHandler(const Triangulation<dim> &triangulation,
                 const unsigned int n_quadrature_points,
                 const unsigned int n_components);

      void reinit();

      cell_iterator 
      begin(const unsigned int level = 0) const;

      active_cell_iterator 
      begin_active(const unsigned int level = 0) const;

      cell_iterator 
      end() const;

      cell_iterator 
      end(const unsigned int level) const;

      active_cell_iterator
      end_active(const unsigned int level) const;

      IteratorRange<typename QPDHandler<dim>::cell_iterator>
      cell_iterators() const;

      IteratorRange<typename QPDHandler<dim>::active_cell_iterator>
      active_cell_iterators() const;

      unsigned int n_quadrature_points() const;

      unsigned int n_components() const;

      const Triangulation<dim> &
      get_triangulation() const;

    private:
      SmartPointer<const Triangulation<dim>, QPDHandler<dim>> tria;

      const unsigned int n_q_points;

      const unsigned int n_comp;

      std::vector<double *> index_data_map;

      std::vector<double> data;

      friend class QPDAccessor<dim>;
  };


  template <int dim>
  class QPDAccessor : public CellAccessor<dim, dim>
  {
    public:
      using BaseClass    = TriaAccessor<dim, dim, dim>;
      using AccessorData = QPDHandler<dim>;

      QPDAccessor();

      QPDAccessor(const Triangulation<dim> *tria,
                  const int                 level,
                  const int                 index,
                  const QPDHandler<dim> *   qpd_handler);

      QPDAccessor(const QPDAccessor<dim> &) = default;

void copy_from(const TriaAccessorBase<dim, dim, dim> &a);

      void copy_from(const QPDAccessor<dim> &a);

      bool operator==(const QPDAccessor<dim> &a) const;

      bool operator!=(const QPDAccessor<dim> &a) const;

      QPDAccessor<dim> &operator=(const QPDAccessor<dim> &a) = delete;

      QPDAccessor<dim> &operator=(QPDAccessor<dim> &&) = default;

      double get_scalar(const unsigned int quadrature_point,
                        const unsigned int first_component) const;

      Tensor<1,dim>
      get_vector(const unsigned int quadrature_point,
                 const unsigned int first_component) const;

      SymmetricTensor<2,dim>
      get_symmetric_tensor(const unsigned int quadrature_point,
                           const unsigned int first_component) const;

      void get(const unsigned int component,
               std::vector<double> &values) const;

      void set(const unsigned int quadrature_point,
               const unsigned int first_component,
               const double       value) const;

      void set(const unsigned int quadrature_point,
               const unsigned int first_component,
               const Tensor<1,dim> &value) const;

      void set(const unsigned int quadrature_point,
               const unsigned int first_component,
               const SymmetricTensor<2,dim> &value) const;

      void set(const unsigned int component,
               const std::vector<double> &values) const;

      DeclExceptionMsg(ExcInvalidObject,
                       "This accessor object has not been "
                       "associated with any QPDHandler object.");

      DeclException0(ExcCantCompareIterators);

    private:
      double *handle(const unsigned int quadrature_point,
                     const unsigned int component) const;

      QPDHandler<dim> *qpd_handler;
  };


  /*------------------------ inline functions: QPDHandler --------------------------*/

  template <int dim>
  inline unsigned int
  QPDHandler<dim>::n_quadrature_points() const
  {
    return n_q_points;
  }


  template <int dim>
  inline unsigned int
  QPDHandler<dim>::n_components() const
  {
    return n_comp;
  }


  template <int dim>
  inline const Triangulation<dim> &
  QPDHandler<dim>::get_triangulation() const
  {
    Assert(tria != nullptr,
           ExcMessage("This QPDHandler object has not been associated "
                      "with a triangulation."));
    return *tria;
  }


  /*------------------------ inline functions: QPDAccessor --------------------------*/

  template <int dim>
  inline QPDAccessor<dim>::QPDAccessor()
  {
    Assert(false, ExcInvalidObject());
  }


  template <int dim>
  inline
  QPDAccessor<dim>::QPDAccessor(const Triangulation<dim> *tria,
                                const int                 level,
                                const int                 index,
                                const QPDHandler<dim> *   qpd_handler_)
    : CellAccessor<dim, dim>(tria, level, index)
    , qpd_handler(const_cast<QPDHandler<dim>*>(qpd_handler_))
  {
    Assert(tria == nullptr || &qpd_handler->get_triangulation() == tria,
           ExcMessage("You can't create a QPD accessor in which the QPDHandler object "
                      "uses a different triangulation than the one you pass as argument."));
  }


  template <int dim>
  inline void
  QPDAccessor<dim>::copy_from(const TriaAccessorBase<dim, dim, dim> &a)
  {
    Assert(qpd_handler != nullptr, ExcInvalidObject());
    BaseClass::copy_from(a);
  }

  
  template <int dim>
  inline void
  QPDAccessor<dim>::copy_from(const QPDAccessor<dim> &a)
  {
    BaseClass::copy_from(a);
    qpd_handler = a.qpd_handler;
  }


  template <int dim>
  inline bool
  QPDAccessor<dim>::operator==(const QPDAccessor<dim> &a) const
  {
    Assert(qpd_handler == a.qpd_handler, ExcCantCompareIterators());
    return BaseClass::operator==(a);
  }


  template <int dim>
  inline bool
  QPDAccessor<dim>::operator!=(const QPDAccessor<dim> &a) const
  {
    Assert(qpd_handler == a.qpd_handler, ExcCantCompareIterators());
    return BaseClass::operator!=(a);
  }


  template <int dim>
  inline double *
  QPDAccessor<dim>::handle(const unsigned int quadrature_point,
                           const unsigned int component) const
  {
    Assert(this->active_cell_index() < qpd_handler->index_data_map.size(),
           ExcIndexRange(this->active_cell_index(), 0, qpd_handler->index_data_map.size()));
    Assert(quadrature_point < qpd_handler->n_q_points,
           ExcIndexRange(quadrature_point, 0, qpd_handler->n_q_points));
    Assert(component < qpd_handler->n_comp,
           ExcIndexRange(component, 0, qpd_handler->n_comp));

    return qpd_handler->index_data_map[this->active_cell_index()]
           + qpd_handler->n_comp * quadrature_point + component;
  }


  template <int dim>
  inline double
  QPDAccessor<dim>::get_scalar(const unsigned int quadrature_point,
                               const unsigned int first_component) const
  {
    double *p = handle(quadrature_point, first_component);
    Assert(p != nullptr, ExcInternalError());
    return *p;
  }


  template <int dim>
  inline Tensor<1,dim>
  QPDAccessor<dim>::get_vector(const unsigned int quadrature_point,
                               const unsigned int first_component) const
  {
    double *p = handle(quadrature_point, first_component);
    Assert(p != nullptr, ExcInternalError());
    return Tensor<1,dim>(ArrayView<const double>(p, dim));
  }


  template <int dim>
  inline SymmetricTensor<2,dim>
  QPDAccessor<dim>::get_symmetric_tensor(const unsigned int quadrature_point,
                                         const unsigned int first_component) const
  {
    SymmetricTensor<2,dim> result;
    double *p = handle(quadrature_point, first_component);
    Assert(p != nullptr, ExcInternalError());
    for (unsigned int c = 0; c < SymmetricTensor<2,dim>::n_independent_components; ++c, ++p)
      result.access_raw_entry(c) = *p;

    return result;
  }


  template <int dim>
  inline void
  QPDAccessor<dim>::get(const unsigned int component,
                        std::vector<double> &values) const
  {
    Assert(this->active_cell_index() < qpd_handler->index_data_map.size(),
           ExcIndexRange(this->active_cell_index(), 0, qpd_handler->index_data_map.size()));
    Assert(component < qpd_handler->n_comp,
           ExcIndexRange(component, 0, qpd_handler->n_comp));
    Assert(values.size() == qpd_handler->n_q_points,
           ExcDimensionMismatch(values.size(), qpd_handler->n_q_points));

    double *p = handle(0, component);
    Assert(p != nullptr, ExcInternalError());
    for (auto val = values.begin(); val != values.end(); ++val)
    {
      *val = *p;
      p += qpd_handler->n_comp;
    }
  }


  template <int dim>
  inline void
  QPDAccessor<dim>::set(const unsigned int quadrature_point,
                        const unsigned int first_component,
                        const double       value) const
  {
    double *p = handle(quadrature_point, first_component);
    Assert(p != nullptr, ExcInternalError());
    *p = value;
  }


  template <int dim>
  inline void
  QPDAccessor<dim>::set(const unsigned int quadrature_point,
                        const unsigned int first_component,
                        const Tensor<1,dim> &value) const
  {
    double *p = handle(quadrature_point, first_component);
    Assert(p != nullptr, ExcInternalError());
    for (unsigned int d = 0; d < dim; ++d, ++p)
      *p = value[d];
  }


  template <int dim>
  inline void
  QPDAccessor<dim>::set(const unsigned int quadrature_point,
                        const unsigned int first_component,
                        const SymmetricTensor<2,dim> &value) const
  {
    ArrayView<const double> view = make_array_view(value);
    double *p = handle(quadrature_point, first_component);
    Assert(p != nullptr, ExcInternalError());
    for (auto val = view.begin(); val != view.end(); ++val, ++p)
      *p = *val;
  }


  template <int dim>
  inline void
  QPDAccessor<dim>::set(const unsigned int component,
                        const std::vector<double> &values) const
  {
    Assert(this->active_cell_index() < qpd_handler->index_data_map.size(),
           ExcIndexRange(this->active_cell_index(), 0, qpd_handler->index_data_map.size()));
    Assert(component < qpd_handler->n_comp,
           ExcIndexRange(component, 0, qpd_handler->n_comp));
    Assert(values.size() == qpd_handler->n_q_points,
           ExcDimensionMismatch(values.size(), qpd_handler->n_q_points));

    double *p = handle(0, component);
    Assert(p != nullptr, ExcInternalError());
    for (auto val = values.begin(); val != values.end(); ++val)
    {
      *p = *val;
      p += qpd_handler->n_comp;
    }
  }
}

#endif
