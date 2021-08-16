// ConsoleApplication6.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include<iterator> 
#include <algorithm>
#include "Poly.h"
#include <math.h> 
#include <limits.h>
#include "FindRoots.h"

using namespace std;

Poly::Poly(vector<long double> x): right_limit_sign(0), left_limit_sign(0) {
    //cout << "Poly::Poly(vector<long double> x)" << endl;
    //reserve the size of vector to not copy data every allocation - optimization
    coff.reserve(x.size());
    vector<long double>::iterator it;
    for (it = x.begin(); it < x.end(); it++) 
        coff.emplace_back(*it);
    find_limits();
}

Poly::Poly(const Poly& p){
    //cout << "Poly::Poly(const Poly& p)" << endl;
    coff.reserve(p.coff.size());
    for (int i = 0; i < p.coff.size(); i++) {
        coff.emplace_back(p.coff[i]);
    }
    find_limits();
}
Poly::Poly(Poly&& p) noexcept {
    //cout << "Poly::Poly(const Poly&& p)" << endl;
    coff.reserve(p.coff.size());
    coff = p.coff;
    vector<long double>().swap(p.coff);
    /*for (int i = 0; i < p.coff.size(); i++) {
        coff.push_back(p.coff[i]);
    }*/
    find_limits();
}
void Poly::find_limits() {
    if (coff.front() >= 0) {
        if (Degree() % 2 == 0) {
            right_limit_sign = 1;
            left_limit_sign = 1;
        }
        else {
            right_limit_sign = 1;
            left_limit_sign = -1;
        }
    }
    else {
        if (Degree() % 2 == 0) {
            right_limit_sign = -1;
            left_limit_sign = -1;
        }
        else {
            right_limit_sign = -1;
            left_limit_sign = 1;
        }
    }
}
Poly::~Poly() {
    vector<long double>().swap(coff);
}
int Poly::get_right_sign() {
    return right_limit_sign;
}
int Poly::get_left_sign() {
    return left_limit_sign;
}

long double Poly::eval(long double x) {
    /*long double val = 0;
    for (int k = 0; k <= Degree(); k++) {
        val += coff[k] * pow(x, Degree() - k);
    }
    return val;*/
    // Evaluate value of polynomial using Horner's method - optimization
    long double val = coff.front(); // Initialize result
    for (int i = 1; i <= Degree(); i++) 
        val = val * x + coff[i];
    return val;
}
long double& Poly::operator [](int n){
    return coff[n];
}
/*
bool Poly::operator == (Poly& other) {
    if (coff.size() != other.coff.size()) {
        return false;
    }
    for (int i = 0; i != coff.size(); ++i) {
        if (coff[i] != other[i]) {
            return false;
        }
    }
    return true;
}

bool Poly::operator != (Poly& other) {
    return !(*this == other);
}
   */

Poly Poly::operator / (long double x) {
    //copy constructor for vector class
    vector<long double> a = coff;
    for (long double& mem : a) {
        mem /= x;
    }
    Poly b(a);
    return b;
}

int Poly::Degree() {
    if (coff.empty()) {
        return -1;
    }
    else {
        return static_cast<int>(coff.size()) - 1;
    }
}

vector<long double>& Poly::getdata() {
    return coff;
}


void Poly::print() {
    cout << "print()" << endl;
    int deg = Degree();
    vector<long double>::iterator it;
    for (it = coff.begin(); it < coff.end() -1 ; it++)
        cout << *it << "*x^" << deg-- << " + ";
    cout << *it << endl;
}


