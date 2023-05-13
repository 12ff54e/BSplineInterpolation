/* This is a C++ MEX file for MATLAB.
C++ B-spline interpolation interface
https://github.com/12ff54e/BSplineInterpolation
*/

#include <iostream>
#include <array>
#include <Interpolation.hpp>
#include "mex.hpp"
#include "mexAdapter.hpp"

template<typename T> struct make_inverse_index_sequence_impl;
template<std::size_t... Is> struct make_inverse_index_sequence_impl<std::index_sequence<Is...>>
{
    using type = std::index_sequence<sizeof...(Is) - Is - 1 ...>;
};
template<std::size_t N> using make_inverse_index_sequence = typename make_inverse_index_sequence_impl<std::make_index_sequence<N>>::type;


using namespace matlab::engine;
using namespace  matlab::data;

decltype(auto) get_matlab_arr(auto&& x, auto i, auto... is)
{
    if constexpr (sizeof...(is) == 0)
        return x[i];
    else
        return get_matlab_arr(x[i], is...);
}

template<std::size_t I = 0, typename F, std::size_t N>
void range_invoke_row_major(F&& f, const std::array<std::size_t, N>& range, auto... is)
{
    if constexpr (I == N)
        f(is...);
    else
    {
        for (std::size_t i = 0; i < range[I]; i++)
            range_invoke_row_major<I + 1>(std::forward<F>(f), range, is..., i);
    }
}

template<std::size_t I = 0, typename F, std::size_t N>
void range_invoke_col_major(F&& f, const std::array<std::size_t, N>& range, auto... is)
{
    if constexpr (I == N)
        f(is...);
    else
    {
        //std::cout << "range " << I << " dim " << range[N - 1 - I] << std::endl;
        for (std::size_t i = 0; i < range[N - 1 - I]; i++)
            range_invoke_col_major<I + 1>(std::forward<F>(f), range, i, is...);
    }
}

decltype(auto) reverse_invoke(auto&& f, auto&&... is)
{
    const auto tup = std::tuple{ is... };
    return [&] <std::size_t... Is>(std::index_sequence<Is...>)->decltype(auto) {
        return f(std::get<Is>(tup)...);
    }(make_inverse_index_sequence<sizeof...(is)>{});
}


//constexpr std::size_t max_dim = 3;
//using interp_variant = decltype([]<std::size_t... Is>(std::index_sequence<Is...>) {
//    return  std::variant<std::integral_constant<std::size_t, Is>...>{};
//}(std::make_index_sequence<max_dim>{}));
//
//template<class... Ts> struct Overloaded : Ts... { using Ts::operator()...; };
//template<class... Ts> Overloaded(Ts...) -> Overloaded<Ts...>;


class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;
    std::shared_ptr<MATLABEngine> matlabPtr = getEngine();

public:
    template<std::size_t Dim>
    void interpolation_nd(const uint64_t order, const TypedArray<bool>& is_periodic, const TypedArray<double>& range, const TypedArray<double>& mesh_in, const TypedArray<double>& coor_in, const TypedArray<std::size_t>& derivative_in, TypedArray<double>& result)
    {
        //// debug print info
        //std::cout << "dimensions = " << Dim << " interp order " << order << std::endl;
        //std::cout << "range dimensions " << range.getDimensions()[0] << ", " << range.getDimensions()[1] << std::endl;
        //for (std::size_t i = 0; i < Dim; i++)
        //{
        //    std::cout << "dim " << i  << " size " << mesh_in.getDimensions()[i] << " is periodic :" << is_periodic[i] << ", range (" << range[i][0] << "," << range[i][1] << ")" << std::endl;
        //}
        ////std::cout << "interp number = " << interp_num << std::endl;

        // dim
        const std::array<std::size_t, Dim> Dims = [&]<std::size_t... Is>(auto & grid, std::index_sequence<Is...>) {
            return std::array<std::size_t, Dim>{grid.getDimensions()[Is]...};
        }(mesh_in, std::make_index_sequence<Dim>{});

        // mesh
        auto f_nd = [&]<std::size_t... Is>(std::index_sequence<Is...>) {
            //((std::cout << "Mesh size = " << Dims[Is] + std::size_t{ is_periodic[Is] } << std::endl), ...);
            return  intp::Mesh<double, Dim>{ Dims[Is] + std::size_t{ is_periodic[Is] }... };
        }(make_inverse_index_sequence<Dim>{});

        range_invoke_col_major([&](auto... is)
            {
                reverse_invoke(f_nd, is...) = get_matlab_arr(mesh_in, is...);
            }, Dims);

        //if constexpr (Dim==2)
        //{
        //    for (std::size_t i = 0; i < Dims[0]; i++)
        //        for (std::size_t j = 0; j < Dims[1]; j++)
        //            std::cout << f_nd(j, i) << std::endl;
        //}

		//if constexpr (Dim == 1)
		//{
  //          std::cout << "mesh value " << std::endl;
		//	for (std::size_t i = 0; i < Dims[0]; i++)
		//			std::cout << i << " : " << f_nd(i) << std::endl;
  //          f_nd(12) = 0;
		//}


        // interpolation initial
        auto interpolation_function = [&]<std::size_t... Is>(std::index_sequence<Is...>) {
            return  intp::InterpolationFunction<double, Dim>{
                order, { is_periodic[Is]... }, f_nd,
                    std::make_pair(double{ range[Is][0] }, double{ range[Is][1] })...};
        }(make_inverse_index_sequence<Dim>{});

        // interpolation
        for (std::size_t i = 0; i < coor_in.getDimensions()[0]; i++)
        {
            result[i] = [&]<std::size_t... Is>(std::index_sequence<Is...>) {
                //return  interpolation_function(double{ coor_in[i][Is] }...);
                return  interpolation_function.derivative({ double{ coor_in[i][Is] }... },{ derivative_in[Is]...});
            }(make_inverse_index_sequence<Dim>{});
        }

        //if constexpr (Dim == 3) {
        //    std::cout << interpolation_function({ 4,3,2 }) << std::endl;
        //    std::cout << interpolation_function.derivative({ 4,3,2 }, { 0,0,0 }) << std::endl;
        //}
    }


    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) override {
        // inputs : ( order, is periodic, range, array, array interpolation, derivative )
        // outputs : ( result )
        checkArguments(outputs, inputs);

        const uint64_t order = inputs[0][0];
        const TypedArray<bool> is_periodic = std::move(inputs[1]);
        const TypedArray<double> range = std::move(inputs[2]);
        const TypedArray<double> mesh_in = std::move(inputs[3]);
        const TypedArray<double> coor_in = std::move(inputs[4]);
        const TypedArray<std::size_t> derivative_in = std::move(inputs[5]);
        TypedArray<double> result = factory.createArray<double>({ coor_in.getDimensions()[0] });

        std::size_t dim = 0;
        for (std::size_t i = 0; i < mesh_in.getDimensions().size(); i++)
            if (mesh_in.getDimensions()[i] > 1) dim++;

        //switch (const auto dim = mesh_in.getDimensions().size())
        switch (dim)
        {
        case 1:
            interpolation_nd<1>(order, is_periodic, range, mesh_in, coor_in, derivative_in, result);
            break;
        case 2:
            interpolation_nd<2>(order, is_periodic, range, mesh_in, coor_in, derivative_in, result);
            break;
        case 3:
            interpolation_nd<3>(order, is_periodic, range, mesh_in, coor_in, derivative_in, result);
            break;
        default:
            std::cout << "unsupport dim, you need add dim " << dim << " in bspline.cpp." << std::endl;
            throw;
        }

        TypedArray<double> outputArray = factory.createArray({ result.getDimensions()[0] }, result.cbegin(), result.cend());
		outputs[0] = std::move(outputArray);

    }


    void checkArguments(matlab::mex::ArgumentList& outputs, matlab::mex::ArgumentList& inputs) {

        if (inputs.size() != 6) {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("Six input required") }));
        }

        if(inputs[0].getType()!=ArrayType::UINT64 || inputs[1].getType()!=ArrayType::LOGICAL || inputs[2].getType()!=ArrayType::DOUBLE || inputs[3].getType() != ArrayType::DOUBLE || inputs[4].getType() != ArrayType::DOUBLE || inputs[5].getType() != ArrayType::UINT64) {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("Input type error") }));
        }

        const uint64_t order = inputs[0][0];
        const TypedArray<bool> is_periodic = inputs[1];
        const TypedArray<double> range = inputs[2];
        const TypedArray<double> mesh_in = inputs[3];
        const TypedArray<double> coor_in = inputs[4];
        const TypedArray<std::size_t> derivative_in = inputs[5];
        std::size_t dim = 0;
        for (std::size_t i = 0; i < mesh_in.getDimensions().size(); i++)
            if (mesh_in.getDimensions()[i] > 1) dim++;

        if (range.getDimensions().size() != 2 || range.getDimensions()[0]!=dim || range.getDimensions()[1] != 2) {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("Range require Dim * 2 array") }));
        }

        if(mesh_in.getDimensions()[0]==1)
        {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("Mesh should be col vector for 1 dimension") }));
        }

        if(is_periodic.getDimensions()[0]!= dim || derivative_in.getDimensions()[0] != dim)
        {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("is_periodic, derivative_in require length dim array") }));
        }

        if (coor_in.getDimensions().size() != 2 || coor_in.getDimensions()[1] != dim)
        {
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar("Interpolate coordinate require N * dim array") }));
        }

    }

};


//// Pass data array to MATLAB sqrt function 
//// And return results.
//auto result = matlabPtr->feval(u"sqrt", inputArray);

//// const auto dim = 3;
       // const auto dim = mesh_in.getDimensions().size();

       //// std::array<bool, 3> is_periodic2{};

       // const std::size_t interp_num = coor_in.getDimensions()[0];
       // TypedArray<double> result = factory.createArray<double>({ interp_num });

       // //TypedArray<double> Yp = factory.createArray<double>({ numRows, numColumns });
       //// yPrime(Yp, Y);
       //// outputs[0] = std::move(Yp);
       // //std::cout << numRows << " " << numColumns << std::endl;
       // //for (auto& i: Y)
       // //    std::cout << i << " ";

       // //const std::array<std::size_t,Dim> Dims=[&]<std::size_t... Is>(std::index_sequence<Is...>){
       // //    ((std::cout << inputs[3].getDimensions()[Is] << " "), ...); std::cout << std::endl;
       // //    return std::array<std::size_t, Dim>{inputs[3].getDimensions()[Is]...};
       // //}(std::make_index_sequence<Dim>{});


       // // debug print info
       // std::cout << "dimensions = " << dim << " interp order " << order << std::endl;
       // std::cout << "range dimensions " << range.getDimensions()[0] << ", " << range.getDimensions()[1] << std::endl;
       // for(std::size_t i=0; i<dim; i++)
       // {
       //     std::cout << "dim " << i << " is periodic :" << is_periodic[i] << ", range (" << range[i][0] << "," << range[i][1] << ")" << std::endl;
       // }
       // std::cout << "interp number = " << interp_num << std::endl;


       //// const TypedArray<double> mesh = std::move(inputs[3]);


       // //const TypedArray<double> Y = inputs[0];

       // const std::array<std::size_t, Dim> Dims = [&]<std::size_t... Is>(auto & grid, std::index_sequence<Is...>) {
       //     ((std::cout << grid.getDimensions()[Is] << " "), ...); std::cout << std::endl;
       //     return std::array<std::size_t, Dim>{grid.getDimensions()[Is]...};
       // }(mesh_in, std::make_index_sequence<Dim>{});
       // // reverse index sequence!!!


       // std::cout << get_matlab_arr(mesh_in, 3, 2, 1) << std::endl; // reverse_invoke
       // std::cout << reverse_invoke([&](auto... is)
       // {
       //       return  get_matlab_arr(mesh_in, is...);
       // },1,2,3) << std::endl;

       // //// print mesh
       // //for (std::size_t i = 0; i < Dims[2]; i++)
       // //    for (std::size_t j = 0; j < Dims[1]; j++)
       // //        for (std::size_t k = 0; k < Dims[0]; k++)
       // //            std::cout  << mesh_in[k][j][i] << " "; // << "(" << k << "," << j << "," << i << ")"
       // //std::cout << std::endl;

       // //range_invoke_col_major([&](auto... is)
       // //    {
       // //        //std::cout << "(";
       // //        //((std::cout << is << " "),...);
       // //        //std::cout << ")";
       // //        std::cout << get_matlab_arr(mesh_in, is...) << " ";
       // //    }, Dims);
       // //std::cout << std::endl;

       // // mesh
       // intp::Mesh<double, 3> f3d{ Dims[2], Dims[1], Dims[0]};
       // for (size_t i = 0; i < f3d.dim_size(0); ++i) {
       //     for (size_t j = 0; j < f3d.dim_size(1); ++j)
       //         for (size_t k = 0; k < f3d.dim_size(2); ++k)
       //     { f3d(i, j, k) = mesh_in[k][j][i]; }
       // }
       // std::cout << f3d( 3, 2, 1) << std::endl;


       // //intp::Mesh<double, Dim> f_nd{ Dims };

       // auto f_nd = [&]<std::size_t... Is>(std::index_sequence<Is...>) {
       //     return  intp::Mesh<double, Dim>{ Dims[Is] + std::size_t{ is_periodic[Is] }... };
       // }(make_inverse_index_sequence<Dim>{});

       // range_invoke_col_major([&](auto... is)
       //     {
       //         reverse_invoke(f_nd, is...) = get_matlab_arr(mesh_in, is...);
       //     }, Dims);

       // std::cout << f_nd(3, 2, 1) << std::endl;


       ///* [&]<std::size_t I>(this auto self, std::integral_constant<std::size_t,I>, auto... is)
       // {
       //     for (size_t i = 0; i < f_nd.dim_size(I); ++i)
       //         if constexpr (I == 0)
       //         {
          //          [&] <std::size_t... Is>(std::index_sequence<Is...>, auto tup)
          //          {
       //                 f_nd(is...) = 
       //             }(std::make_index_sequence<Dim>{}, std::tuple{ is... });
       //         }
       //         else
       //             self(std::integral_constant<std::size_t, I - 1>{}, i, is...);
       // }(std::integral_constant<std::size_t, Dim>{});*/



       // // interp initial
       // auto interp_function = [&]<std::size_t... Is>(std::index_sequence<Is...>) {
       //     return  intp::InterpolationFunction<double, Dim>{
       //         order, { is_periodic[Is]... }, f_nd,
       //             std::make_pair(double(range[Is][0]), double(range[Is][1]))...};
       // }(make_inverse_index_sequence<Dim>{});

       // // interp
       // for (std::size_t i = 0; i < interp_num; i++)
       // {
       //    // std::array<std::size_t, 3> indexs{ coor_in[i][0],coor_in[i][1],coor_in[i][2] };
       //     double x{ coor_in[i][2] }, y{ coor_in[i][1] }, z{ coor_in[i][0] };
       //     /*result[i] = interp_function(x,y,z);*/
       //     //result[i] = interp_function(z, y, x);
       //     result[i] = [&]<std::size_t... Is>(std::index_sequence<Is...>) {
       //         return  interp_function(double{ coor_in[i][Is] }...);
       //     }(make_inverse_index_sequence<Dim>{});

       //     std::cout << "index " << x << " " << y << " " << z << " " << result[i] << std::endl;// <<"index " << x << " " << y << " " << z << " "
       // }
       // std::cout << std::endl;
