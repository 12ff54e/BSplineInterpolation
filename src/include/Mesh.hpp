#pragma once

#include <array>
#include <vector>

#include "util.hpp"

/**
 * @brief A multi dimension mesh storing data on each mesh point
 *
 * @tparam T Type of data stored
 * @tparam D Dimension
 * @tparam Alloc Allocator type, defaulted to std::allocator<T>
 */
template <typename T, unsigned D, typename Alloc = std::allocator<T>>
class Mesh {
   public:
    using size_type = unsigned int;
    using val_type = T;
    const static size_type dim = D;
    using Indices = std::array<size_type, dim>;
    using allocator_type = Alloc;

   private:
    using iterator = typename std::vector<val_type>::const_iterator;

    template <typename U>
    class skip_iterator {
       public:
        using value_type = U;
        using difference_type = std::ptrdiff_t;
        using pointer = U*;
        using reference = U&;
        using iterator_category = std::random_access_iterator_tag;

       private:
        pointer ptr;
        size_type step_length;

       public:
        skip_iterator(value_type* ptr, size_type step_length)
            : ptr(ptr), step_length(step_length) {}

        // allow iterator to const_iterator conversion
        template <typename V>
        skip_iterator(
            skip_iterator<typename std::enable_if<
                std::is_same<V, typename std::remove_const<value_type>::type>::
                    value,
                V>::type> other)
            : ptr(other.ptr), step_length(other.step_length) {}

        // forward iterator requirement

        reference operator*() { return *ptr; }
        reference operator->() { return *ptr; }

        bool operator==(skip_iterator other) {
            return this->ptr == other.ptr &&
                   this->step_length == other.step_length;
        }
        bool operator!=(skip_iterator other) { return !operator==(other); }

        skip_iterator& operator++() {
            ptr += step_length;
            return *this;
        }
        skip_iterator operator++(int) {
            skip_iterator tmp(*this);
            operator++();
            return tmp;
        }

        // bidirectional iterator requirement

        skip_iterator& operator--() {
            ptr -= step_length;
            return *this;
        }
        skip_iterator operator--(int) {
            skip_iterator tmp(*this);
            operator--();
            return tmp;
        }

        // random access iterator requirement

        skip_iterator& operator+=(size_type n) {
            ptr += n * step_length;
            return *this;
        }
        skip_iterator operator+(size_type n) {
            skip_iterator tmp(*this);
            return tmp += n;
        }
        friend skip_iterator operator+(size_type n, skip_iterator it) {
            return it += n;
        }

        skip_iterator& operator-=(size_type n) {
            ptr -= n * step_length;
            return *this;
        }
        skip_iterator operator-(size_type n) {
            skip_iterator tmp(*this);
            return tmp -= n;
        }

        difference_type operator-(skip_iterator other) {
            return (ptr - other.ptr) / step_length;
        }

        reference operator[](size_type n) { return *(ptr + n * step_length); }

        bool operator<(skip_iterator other) { return other - *this > 0; }
        bool operator>(skip_iterator other) { return other < *this; }
        bool operator<=(skip_iterator other) { return !(*this > other); }
        bool operator>=(skip_iterator other) { return !(*this < other); }
    };

    /**
     * @brief Stores the mesh content in row-major format.
     */
    std::vector<val_type, allocator_type> storage;
    /**
     * @brief The i-th element stores the sub-mesh size when the first (dim-i)
     * coordinates are specified.
     */
    std::array<size_type, dim + 1> __dim_acc_size;
    std::array<size_type, dim> __dim_size;

    // auxilary methods

    /**
     * @brief Set the __dim_acc_size object. Check description of __dim_acc_size
     * for details.
     */
    void set_dim_acc_size() {
        for (size_type i = 0; i <= dim; ++i) {
            __dim_acc_size[i] = 1;
            for (size_type j = 0; j < i; ++j) {
                __dim_acc_size[i] *= __dim_size[dim - j - 1];
            }
        }
        return;
    }

    /**
     * @brief Convert multi-dimesion index to one dimension index in storage
     * vector.
     *
     * @param ind
     * @param indices
     * @return size_type
     */
    template <typename... _Indices>
    size_type indexing(size_type ind, _Indices... indices) const {
        return ind * __dim_acc_size[sizeof...(indices)] + indexing(indices...);
    }

    size_type indexing(Indices ind_arr) const {
        size_type ind{};
        for (size_type i = 0; i < dim; ++i) {
            ind += ind_arr[i] * __dim_acc_size[dim - i - 1];
        }
        return ind;
    }

    /**
     * @brief Convert one dimension index in storage vector to multi-dimension
     * indices
     *
     * @param total_ind
     * @return std::array<size_type, dim>
     */
    Indices dimwise_indices(size_type total_ind) {
        Indices indices;

        for (unsigned i = 0; i < dim; ++i) {
            indices[i] = total_ind / __dim_acc_size[dim - i - 1];
            total_ind %= __dim_acc_size[dim - i - 1];
        }

        return indices;
    }

    constexpr size_type indexing() const { return 0; }

    template <typename U, unsigned DD>
    friend class InterpolationFunction;  // friend class forward declaration,
                                         // InterpolationFunction need access
                                         // to indexing function during
                                         // interpolation procedure

   public:
    explicit Mesh(std::initializer_list<size_type> il,
                  const allocator_type& alloc = allocator_type()) {
        std::copy(il.begin(), il.end(), __dim_size.begin());
        set_dim_acc_size();
        storage.resize(__dim_acc_size.back(), val_type{});
    }

    explicit Mesh(size_type n, const allocator_type& alloc = allocator_type()) {
        std::fill(__dim_size.begin(), __dim_size.end(), n);
        set_dim_acc_size();
        storage.resize(__dim_acc_size.back(), val_type{});
    }

    template <typename InputIter,
              typename = typename std::enable_if<
                  dim == 1u &&
                  std::is_convertible<typename std::iterator_traits<
                                          InputIter>::iterator_category,
                                      std::input_iterator_tag>::value>::type>
    explicit Mesh(InputIter iter_begin,
                  InputIter iter_end,
                  const allocator_type& alloc = allocator_type())
        : storage(iter_begin, iter_end, alloc),
          __dim_size{(size_type)storage.size()},
          __dim_acc_size{1u, (size_type)storage.size()} {}

    template <typename Array,
              typename = typename std::enable_if<
                  dim == 1u && util::is_iteratable<Array>::value>::type>
    explicit Mesh(const Array& array,
                  const allocator_type& alloc = allocator_type())
        : Mesh(array.begin(), array.end(), alloc) {}

   public:
    // properties

    size_type size() const { return storage.size(); }

    size_type dim_size(size_type dim_ind) const { return __dim_size[dim_ind]; }

    // modifiers

    void resize(Indices sizes) {
        __dim_size = sizes;
        set_dim_acc_size();
        storage.resize(__dim_acc_size.back());
    }

    // element access

    template <typename... _Indices>
    val_type& operator()(_Indices... indices) {
        return storage[indexing(indices...)];
    }

    template <typename... _Indices>
    val_type operator()(_Indices... indices) const {
        return storage[indexing(indices...)];
    }

    const val_type* data() const { return storage.data(); }

    // iterator

    /**
     * @brief Begin const_iterator to underlying container.
     *
     * @return iterator
     */
    iterator begin() const { return storage.cbegin(); }
    /**
     * @brief End const_iterator to underlying container.
     *
     * @return iterator
     */
    iterator end() const { return storage.cend(); }

    skip_iterator<val_type> begin(size_type dim_ind, Indices indices) {
        indices[dim_ind] = 0;
        return skip_iterator<val_type>(storage.data() + indexing(indices),
                                       __dim_acc_size[dim - dim_ind - 1]);
    }
    skip_iterator<val_type> end(size_type dim_ind, Indices indices) {
        indices[dim_ind] = __dim_size[dim_ind];
        return skip_iterator<val_type>(storage.data() + indexing(indices),
                                       __dim_acc_size[dim - dim_ind - 1]);
    }
    skip_iterator<const val_type> begin(size_type dim_ind,
                                        Indices indices) const {
        indices[dim_ind] = 0;
        return skip_iterator<const val_type>(storage.data() + indexing(indices),
                                             __dim_acc_size[dim - dim_ind - 1]);
    }
    skip_iterator<const val_type> end(size_type dim_ind,
                                      Indices indices) const {
        indices[dim_ind] = __dim_size[dim_ind];
        return skip_iterator<const val_type>(storage.data() + indexing(indices),
                                             __dim_acc_size[dim - dim_ind - 1]);
    }

    Indices iter_indices(iterator iter) {
        return dimwise_indices(std::distance(begin(), iter));
    }

    // debug

    const std::array<unsigned, dim + 1>& dim_acc_size() const {
        return __dim_acc_size;
    }
};
