#pragma once
#include "Poly.h"
#include <functional>
vector<long double> Solve(Poly Pol_solve,bool show_progress, long double epsilon);
long double bisection_n_newtoon(Poly& pol, Poly& derv, long double& a, long double& b, function<int(long double)>& value_sign,long double eps);
vector<long double> external_root(Poly& pol, Poly& diff, long double extreme_p, long double limit_sign, long double init_step, function<int(long double)>& value_sign, long double eps);
vector<long double> internal_roots(Poly& pol, Poly& diff, vector<long double>& diff_sols, function<int(long double)>& value_sign, long double eps);