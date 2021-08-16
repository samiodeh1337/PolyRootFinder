#pragma once
// ConsoleApplication6.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
using namespace std;


class Poly {
private:
    vector<long double> coff;
    int right_limit_sign;
    int left_limit_sign;

    void find_limits();
public:

    Poly(vector<long double> x);
    Poly(const Poly& p);
    Poly(Poly&& p) noexcept;
    ~Poly();
    long double& operator [](int n);
    /*bool operator == (Poly& other);
    bool operator != (Poly& other);*/
    Poly operator / (long double x);
    int get_right_sign();
    int get_left_sign();
    int Degree();
    long double eval(long double x);
    vector<long double>& getdata();
    void print();
};


