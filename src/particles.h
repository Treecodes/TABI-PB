#ifndef H_TABIPB_PARTICLES_STRUCT_H
#define H_TABIPB_PARTICLES_STRUCT_H

#include <array>
#include <vector>
#include <cstdlib>

#include "params.h"


template <typename order_iterator, typename value_iterator>
static void apply_order(order_iterator order_begin, order_iterator order_end, value_iterator v_begin);

template <typename order_iterator, typename value_iterator>
static void apply_unorder(order_iterator order_begin, order_iterator order_end, value_iterator v_begin);


class Particles
{
protected:
    const struct Params& params_;
    
    std::size_t num_;
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> z_;
    std::vector<std::size_t> order_;
    
public:
    Particles(const struct Params& params) : params_(params) {};
    ~Particles() = default;
   
    std::size_t num() const { return num_; };
    const double* x_ptr() const { return x_.data(); };
    const double* y_ptr() const { return y_.data(); };
    const double* z_ptr() const { return z_.data(); };

    int partition_8(std::size_t, std::size_t, std::array<std::size_t, 16>&);
    const std::array<double, 6> bounds(std::size_t begin, std::size_t end) const;
    
    virtual void reorder() = 0;
    virtual void unorder() = 0;
    
    virtual void copyin_to_device() const = 0;
    virtual void delete_from_device() const = 0;
};


template <typename order_iterator, typename value_iterator>
static void apply_order(order_iterator order_begin, order_iterator order_end, value_iterator v_begin)
{
    using value_t = typename std::iterator_traits< value_iterator >::value_type;
    using index_t = typename std::iterator_traits< order_iterator >::value_type;
    
    auto v_end = v_begin + std::distance(order_begin, order_end) + 1;
    std::vector<value_t> tmp(v_begin, v_end);

    std::for_each(order_begin, order_end,
                  [&tmp, &v_begin](index_t idx){ *v_begin = tmp[idx]; std::advance(v_begin, 1); });
}


template <typename order_iterator, typename value_iterator>
static void apply_unorder(order_iterator order_begin, order_iterator order_end, value_iterator v_begin)
{
    using value_t = typename std::iterator_traits< value_iterator >::value_type;
    using index_t = typename std::iterator_traits< order_iterator >::value_type;
    
    auto v_end = v_begin + std::distance(order_begin, order_end);
    std::vector<value_t> tmp(v_begin, v_end);
    
    auto v_iter = v_begin;
    std::for_each(order_begin, order_end,
                  [&tmp, &v_iter](index_t idx){ tmp[idx] = *v_iter; std::advance(v_iter, 1); });
    
    std::copy(tmp.begin(), tmp.end(), v_begin);
}

#endif /* H_TABIPB_PARTICLES_STRUCT_H */
