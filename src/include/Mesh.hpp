#ifndef INTP_MESH
#define INTP_MESH

#include <array>
#include <vector>

#include "util.hpp"

namespace intp {

template <size_t D>
class MeshDimension {
   public:
    using size_type = size_t;
    constexpr static size_type dim = D;
    using index_type = std::array<size_type, dim>;

   private:
    index_type dim_size_;

   public:
    MeshDimension() = default;

    MeshDimension(index_type dim_size) : dim_size_(dim_size) {}

    template <typename... Args,
              typename = typename std::enable_if<sizeof...(Args) == dim>::type>
    MeshDimension(Args... args) : dim_size_{static_cast<size_type>(args)...} {}

    MeshDimension(size_type n) {
        std::fill(dim_size_.begin(), dim_size_.end(), n);
    }

    // properties

    size_type size() const {
        size_type s = 1;
        for (auto&& d : dim_size_) { s *= d; }
        return s;
    }

    size_type dim_size(size_type dim_ind) const { return dim_size_[dim_ind]; }
    size_type& dim_size(size_type dim_ind) { return dim_size_[dim_ind]; }

    size_type dim_acc_size(size_type dim_ind) const {
        size_type das = 1;
        for (size_type d = 0; d < dim_ind; ++d) {
            das *= dim_size_[dim - d - 1];
        }
        return das;
    }

    /**
     * @brief Convert to underlying dimensions array.
     */
    operator index_type() const { return dim_size_; }

    /**
     * @brief Convert multi-dimension index to one dimension index in storage
     * vector.
     *
     * @return size_type
     */
    template <typename... Indices>
    size_type indexing_safe(Indices... indices) const {
        return indexing_safe(index_type{static_cast<size_type>(indices)...});
    }

    size_type indexing_safe(const index_type& ind_arr) const {
        size_type ind{};
        size_type sub_mesh_size = 1;
        for (size_type d = 0; d < dim; ++d) {
            const size_type d_r = dim - d - 1;
            INTP_ASSERT(ind_arr[d_r] >= 0 && ind_arr[d_r] < dim_size_[d_r],
                        std::string("Mesh access out of range at dim ") +
                            std::to_string(d_r) + std::string(", index ") +
                            std::to_string(ind_arr[d_r]) +
                            std::string(" is out of [0, ") +
                            std::to_string(dim_size_[d_r] - 1) +
                            std::string("]."));
            ind += ind_arr[d_r] * sub_mesh_size;
            sub_mesh_size *= dim_size_[d_r];
        }
        return ind;
    }

    size_type indexing(const index_type& ind_arr) const {
        size_type ind{};
        size_type sub_mesh_size = 1;
        for (size_type d = 0; d < dim; ++d) {
            ind += ind_arr[dim - d - 1] * sub_mesh_size;
            sub_mesh_size *= dim_size_[dim - d - 1];
        }
        return ind;
    }

    /**
     * @brief Convert one dimension index in storage vector to multi-dimension
     * indices
     *
     */
    index_type dimwise_indices(size_type total_ind) const {
        index_type indices;

        for (size_type d = 0; d < dim; ++d) {
            indices[dim - d - 1] = total_ind % dim_size_[dim - d - 1];
            total_ind /= dim_size_[dim - d - 1];
        }

        return indices;
    }

    // modifiers

    void resize(index_type sizes) { dim_size_ = sizes; }
};

/**
 * @brief A multi dimension mesh storing data on each mesh point
 *
 * @tparam T Type of data stored
 * @tparam D Dimension
 * @tparam Alloc Allocator type, defaulted to std::allocator<T>
 */
template <typename T, size_t D, typename Alloc = std::allocator<T>>
class Mesh {
   public:
    using size_type = size_t;
    using val_type = T;
    const static size_type dim = D;
    using index_type = typename MeshDimension<dim>::index_type;
    using allocator_type = Alloc;

   private:
    using container_type = std::vector<val_type, allocator_type>;
    using const_iterator = typename container_type::const_iterator;

    template <typename U>
    class skip_iterator {
       public:
        using value_type = U;
        using difference_type = std::ptrdiff_t;
        using pointer = U*;
        using const_pointer = const U*;
        using reference = U&;
        using iterator_category = std::random_access_iterator_tag;

       private:
        pointer ptr_;
        difference_type step_length_;

       public:
        skip_iterator(value_type* ptr, difference_type step_length)
            : ptr_(ptr), step_length_(step_length) {}

        // allow iterator to const_iterator conversion
        operator skip_iterator<const T>() { return {ptr_, step_length_}; }
        // cast to (const) underlying pointer type
        explicit operator const_pointer() const { return ptr_; }

        // forward iterator requirement

        reference operator*() { return *ptr_; }
        reference operator->() { return ptr_; }

        bool operator==(const skip_iterator& other) const {
            return this->ptr_ == other.ptr_ &&
                   this->step_length_ == other.step_length_;
        }
        bool operator!=(const skip_iterator& other) const {
            return !operator==(other);
        }

        skip_iterator& operator++() {
            ptr_ += step_length_;
            return *this;
        }
        skip_iterator operator++(int) {
            skip_iterator tmp(*this);
            operator++();
            return tmp;
        }

        // bidirectional iterator requirement

        skip_iterator& operator--() {
            ptr_ -= step_length_;
            return *this;
        }
        skip_iterator operator--(int) {
            skip_iterator tmp(*this);
            operator--();
            return tmp;
        }

        // random access iterator requirement

        skip_iterator& operator+=(difference_type n) {
            ptr_ += n * step_length_;
            return *this;
        }
        skip_iterator operator+(difference_type n) {
            skip_iterator tmp(*this);
            return tmp += n;
        }
        friend skip_iterator operator+(difference_type n, skip_iterator it) {
            return it += n;
        }

        skip_iterator& operator-=(difference_type n) {
            ptr_ -= n * step_length_;
            return *this;
        }
        skip_iterator operator-(difference_type n) {
            skip_iterator tmp(*this);
            return tmp -= n;
        }

        difference_type operator-(skip_iterator other) {
            return (ptr_ - other.ptr_) / step_length_;
        }

        reference operator[](difference_type n) {
            return *(ptr_ + n * step_length_);
        }

        bool operator<(const skip_iterator& other) const {
            return other - *this > 0;
        }
        bool operator>(const skip_iterator& other) const {
            return other < *this;
        }
        bool operator<=(const skip_iterator& other) const {
            return !(*this > other);
        }
        bool operator>=(const skip_iterator& other) const {
            return !(*this < other);
        }
    };

    /**
     * @brief Stores the mesh content in row-major format.
     */
    container_type storage_;

    MeshDimension<dim> dimension_;

   public:
    explicit Mesh(const MeshDimension<dim>& mesh_dimension,
                  const allocator_type& alloc = allocator_type())
        : storage_(alloc), dimension_(mesh_dimension) {
        storage_.resize(dimension_.size(), val_type{});
    }

    explicit Mesh(size_type n, const allocator_type& alloc = allocator_type())
        : Mesh(MeshDimension<dim>(n), alloc) {}

    template <typename... Args,
              typename = typename std::enable_if<sizeof...(Args) == dim>::type,
              typename = typename std::enable_if<std::is_integral<
                  typename std::common_type<Args...>::type>::value>::type>
    explicit Mesh(Args... args)
        : Mesh(MeshDimension<dim>(static_cast<size_type>(args)...)) {}

    template <typename InputIter,
              typename = typename std::enable_if<
                  dim == 1u &&
                  std::is_convertible<typename std::iterator_traits<
                                          InputIter>::iterator_category,
                                      std::input_iterator_tag>::value>::type>
    explicit Mesh(std::pair<InputIter, InputIter> range,
                  const allocator_type& alloc = allocator_type())
        : storage_(range.first, range.second, alloc),
          dimension_{static_cast<size_type>(storage_.size())} {}

    // convert constructor from mesh using different allocator
    template <typename Allocator>
    Mesh(const Mesh<val_type, dim, Allocator>& mesh,
         const allocator_type& alloc = allocator_type())
        : Mesh(mesh.dimension(), alloc) {
        std::copy(mesh.begin(), mesh.end(), storage_.begin());
    }

    // properties

    size_type size() const { return storage_.size(); }

    size_type dim_size(size_type dim_ind) const {
        return dimension_.dim_size(dim_ind);
    }

    /**
     * @brief Get the underlying mesh dimension object
     *
     */
    const MeshDimension<dim>& dimension() const { return dimension_; }

    // modifiers

    void resize(index_type sizes) {
        dimension_.resize(sizes);
        storage_.resize(dimension_.size());
    }

    // element access

    template <typename... Indices>
    val_type& operator()(Indices... indices) {
        return storage_[dimension_.indexing_safe(indices...)];
    }

    template <typename... Indices>
    const val_type& operator()(Indices... indices) const {
        return storage_[dimension_.indexing_safe(indices...)];
    }

    val_type& operator()(index_type indices) {
        return storage_[dimension_.indexing(indices)];
    }

    const val_type& operator()(index_type indices) const {
        return storage_[dimension_.indexing(indices)];
    }

    const val_type* data() const { return storage_.data(); }

    // iterator

    /**
     * @brief Begin const_iterator to underlying container.
     *
     * @return iterator
     */
    const_iterator begin() const { return storage_.cbegin(); }
    /**
     * @brief End const_iterator to underlying container.
     *
     * @return iterator
     */
    const_iterator end() const { return storage_.cend(); }

    skip_iterator<val_type> begin(size_type dim_ind, index_type indices) {
        indices[dim_ind] = 0;
        return skip_iterator<val_type>(
            storage_.data() + dimension_.indexing_safe(indices),
            static_cast<typename skip_iterator<val_type>::difference_type>(
                dimension_.dim_acc_size(dim - dim_ind - 1)));
    }
    skip_iterator<val_type> end(size_type dim_ind, index_type indices) {
        indices[dim_ind] = dimension_.dim_size(dim_ind);
        return skip_iterator<val_type>(
            storage_.data() + dimension_.indexing(indices),
            static_cast<typename skip_iterator<val_type>::difference_type>(
                dimension_.dim_acc_size(dim - dim_ind - 1)));
    }
    skip_iterator<const val_type> begin(size_type dim_ind,
                                        index_type indices) const {
        indices[dim_ind] = 0;
        return skip_iterator<const val_type>(
            storage_.data() + dimension_.indexing_safe(indices),
            static_cast<typename skip_iterator<val_type>::difference_type>(
                dimension_.dim_acc_size(dim - dim_ind - 1)));
    }
    skip_iterator<const val_type> end(size_type dim_ind,
                                      index_type indices) const {
        indices[dim_ind] = dimension_.dim_size(dim_ind);
        return skip_iterator<const val_type>(
            storage_.data() + indexing(indices),
            dimension_.dim_acc_size(dim - dim_ind - 1));
    }

    index_type iter_indices(const_iterator iter) const {
        return dimension_.dimwise_indices(
            static_cast<size_type>(std::distance(begin(), iter)));
    }

    index_type iter_indices(skip_iterator<const val_type> skip_iter) const {
        return dimension_.dimwise_indices(static_cast<size_type>(
            static_cast<const typename skip_iterator<const val_type>::pointer>(
                skip_iter) -
            data()));
    }
};

}  // namespace intp

#endif
