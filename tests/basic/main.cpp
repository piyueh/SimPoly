/**
 * \file tests/basic/main.cpp
 * \brief Unit tests for functions in the namespace op.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-01-29
 */


# include <random>
# include <chrono>

# include <gtest/gtest.h>


std::default_random_engine generator;
    

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator = std::default_random_engine(seed);
    
    return RUN_ALL_TESTS();
}
