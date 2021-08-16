#include <iostream>
#include <vector>
#include <chrono>
#include "Poly.h"
#include "FindRoots.h"
#include "inputs.h"
#include <algorithm>

int main()
{
  
    int index = 0;
    Poly B(m);
    auto start = chrono::high_resolution_clock::now();

        vector<long double> solutions = Solve(B, false, 1.0e-37);
        cout << "Solution: " << endl;
        vector<long double>::iterator it;
        for (it = solutions.begin(); it < solutions.end(); it++)
            cout << "Root[" << index++ << "]" << " = " << *it << endl;
        cout << endl;

    

    auto elapsed = chrono::high_resolution_clock::now() - start;
    long double mseconds = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
    cout << "time elapsed: " << mseconds << " ms" << endl;
    return 0;
}