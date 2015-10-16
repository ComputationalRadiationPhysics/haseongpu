#pragma once

// ALPAKA
#include <alpaka/alpaka.hpp>

template<typename T_Acc>
struct RandomGenerator {

    using Gen = decltype(alpaka::rand::generator::createDefault(std::declval<T_Acc const &>(),
								std::declval<uint32_t &>(),
								std::declval<uint32_t &>()));
    using Dist = decltype(alpaka::rand::distribution::createUniformReal<float>(std::declval<T_Acc const &>()));
    
    RandomGenerator(T_Acc const &acc, unsigned const seed, unsigned const subsequence) : gen(alpaka::rand::generator::createDefault(acc, seed, subsequence)),
											 dist(alpaka::rand::distribution::createUniformReal<float>(acc)) {

    }
    
    double operator()(){
	return dist(gen);
    }

private:
    Gen  gen;
    Dist dist;    
};
