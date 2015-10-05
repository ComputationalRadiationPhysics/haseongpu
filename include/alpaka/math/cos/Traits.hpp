/**
* \file
* Copyright 2014-2015 Benjamin Worpitz
*
* This file is part of alpaka.
*
* alpaka is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* alpaka is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with alpaka.
* If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <alpaka/core/Common.hpp>   // ALPAKA_FN_HOST_ACC

#include <type_traits>              // std::enable_if, std::is_base_of, std::is_same, std::decay

namespace alpaka
{
    namespace math
    {
        namespace traits
        {
            //#############################################################################
            //! The cos trait.
            //#############################################################################
            template<
                typename T,
                typename TArg,
                typename TSfinae = void>
            struct Cos;
        }

        //-----------------------------------------------------------------------------
        //! Computes the cosine (measured in radians).
        //!
        //! \tparam T The type of the object specializing Cos.
        //! \tparam TArg The arg type.
        //! \param cos The object specializing Cos.
        //! \param arg The arg.
        //-----------------------------------------------------------------------------
        ALPAKA_NO_HOST_ACC_WARNING
        template<
            typename T,
            typename TArg>
        ALPAKA_FN_HOST_ACC auto cos(
            T const & cos,
            TArg const & arg)
        -> decltype(
            traits::Cos<
                T,
                TArg>
            ::cos(
                cos,
                arg))
        {
            return
                traits::Cos<
                    T,
                    TArg>
                ::cos(
                    cos,
                    arg);
        }

        namespace traits
        {
            //#############################################################################
            //! The Cos specialization for classes with CosBase member type.
            //#############################################################################
            template<
                typename T,
                typename TArg>
            struct Cos<
                T,
                TArg,
                typename std::enable_if<
                    std::is_base_of<typename T::CosBase, typename std::decay<T>::type>::value
                    && (!std::is_same<typename T::CosBase, typename std::decay<T>::type>::value)>::type>
            {
                //-----------------------------------------------------------------------------
                //
                //-----------------------------------------------------------------------------
                ALPAKA_NO_HOST_ACC_WARNING
                ALPAKA_FN_HOST_ACC static auto cos(
                    T const & cos,
                    TArg const & arg)
                -> decltype(
                    math::cos(
                        static_cast<typename T::CosBase const &>(cos),
                        arg))
                {
                    // Delegate the call to the base class.
                    return
                        math::cos(
                            static_cast<typename T::CosBase const &>(cos),
                            arg);
                }
            };
        }
    }
}
