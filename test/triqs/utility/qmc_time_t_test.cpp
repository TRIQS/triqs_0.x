#include <iostream>

#include <triqs/utility/time_pt.hpp>

using triqs::utility::time_pt;

void print_me (double x) { std::cout  << " print " << x<< std::endl ;}

int main(int argc, char **argv)
{
    time_pt zero(0.0, 40.0);
    time_pt beta(40.0, 40.0);

    time_pt t1(10.0, beta), t2(30.0, beta);

    std::cout << "zero: " << zero << std::endl;
    std::cout << "beta: " << beta << std::endl;
    std::cout << "t1: " << t1 << "" << std::endl;
    std::cout << "t2: " << t2 << "" << std::endl;

    std::cout << "t1 + t2: " << t1+t2 << "" << std::endl;
    std::cout << "t1 + t2 - (as double) : " << t1+t2 - (double(t1) + double(t2)) << "" << std::endl;

    std::cout << "t1 - t2: " << t1-t2 << "" << std::endl;
    std::cout << "t1 - t2 - (as double) : " << t1-t2 - (double(t1) - double(t2)) << "" << std::endl;

    std::cout << "t1 * t2: " << t1*t2 << "" << std::endl;
    std::cout << "t1 * t2 - (as double) : " << t1*t2 - (double(t1) * double(t2)) << "" << std::endl;

    std::cout << "t1 / t2: " << t1/t2 << "" << std::endl;
    std::cout << "t1 / t2 - (as double) : " << t1/t2 - (double(t1) / double(t2)) << "" << std::endl;

    std::cout << "t1 + beta: " << t1 + time_pt(beta,beta) << "" << std::endl;
    std::cout << "t1 - beta: " << t1 - time_pt(beta,beta) << "" << std::endl;

    print_me(t1);// to check cast

    std::cout << " -t1: " << - t1 << "" << std::endl;
    
    std::cout << "epsilon: " << time_pt::epsilon(beta) << "" << std::endl;
    std::cout << "- epsilon: " << - time_pt::epsilon(beta) << "" << std::endl;
    std::cout << "t1 - epsilon :" << t1- time_pt::epsilon(beta) << "" << std::endl;

    std::cout  << " num limit max " << std::numeric_limits<uint64_t>::max() <<std::endl ;
    uint64_t n = std::numeric_limits<uint64_t>::max()-2;
    n++;
    std::cout  << " num limit-1  " << n <<std::endl ;
    n++;
    std::cout  << " num limit+0  " << n <<std::endl ;
    n++;
    std::cout  << " num limit+1  " << n <<std::endl ;
    n++;
    std::cout  << " num limit+2  " << n <<std::endl ;



    return 0;
}
