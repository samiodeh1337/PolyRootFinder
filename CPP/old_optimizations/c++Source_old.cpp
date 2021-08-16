#include "FindRoots.h";
#include <algorithm>

#define epsilon 1.0e-120

vector<long double> Solve(Poly Pol_solve) {
    /* auto max = max_element(begin(coff), end(coff));
    long double maxcoff = *max;*/
    vector<long double> coff = Pol_solve.getdata();
    /*long double maxco = 0;
    for (auto it = coff.begin(); it < coff.end(); it++)
        if (fabs(*it) > maxco)
            maxco = fabs(*it);
    */
    auto find_max_fabs = [](vector<long double> m) {
        return *max_element(m.begin(), m.end(), [](const long double& a, const long double& b)
            {
                return fabs(a) < fabs(b);
            });
    };

    /*vector<long double> normalized(coff.begin(), coff.end());

    for (int i = 0; i < normalized.size(); i++)
    {
        normalized[i] = normalized[i] / maxco;
    }*/
    //noneed//transform(normalized.begin(), normalized.end(), normalized.begin(), [&maxcoff](auto& c) {return  c / fabs(maxcoff); });
    //Poly pol(normalized,normalized.size() -1);
    Poly pol = Pol_solve / find_max_fabs(coff);

    vector<Poly> derv;
    derv.reserve(Pol_solve.Degree()+1);
    derv.push_back(pol);
    int degree = Pol_solve.Degree();
    int index = 0;
    while (degree > 1) {
        vector<long double> dervcoff(derv.at(index).getdata().begin(), derv.at(index).getdata().end() - 1);
        int sizeofup = derv.at(index).Degree();
        transform(dervcoff.begin(), dervcoff.end(), dervcoff.begin(), [&sizeofup](auto& c) {return  c * (sizeofup--); });
        /*vector<long double>::iterator it;
        long double maxco = 0;
        for (it = dervcoff.begin(); it < dervcoff.end(); it++)
            if (fabs(*it) > maxco)
                maxco = fabs(*it);*/
        long double max = find_max_fabs(dervcoff);
        for (long double& mem : dervcoff)
        {
            mem /= max;
        }


        //Poly pol(dervcoff, dervcoff.size()-1);
        //derv.push_back(pol);

        derv.emplace_back(dervcoff, dervcoff.size() - 1);
        index++;
        degree--;
    }

    int stack_size = derv.size();
    vector<long double> sols_derv = { -derv[stack_size - 1][1] / derv[stack_size - 1][0] };
    vector<long double> sols_poly;
    cout << sols_derv.front() << endl;
    cout << stack_size << endl;
    int number_of_sols = sols_derv.size();

    long double left_ex, right_ex, right_limit_sign, left_limit_sign;
    long double value_ex, sign_ex, init_step, a, b;
    long double x, sign_x, value_x;
    long double val_a, val_b, sign_a, sign_b, sol;


    for (int i = stack_size - 1; i > 0; i--) {
        vector<long double>().swap(sols_poly);
        if (number_of_sols == 0) {
            left_ex = right_ex = 0;
        }
        else {
            left_ex = sols_derv[0];
            right_ex = sols_derv[number_of_sols - 1];
        }

        right_limit_sign = derv[i - 1].get_right_sign();
        left_limit_sign = derv[i - 1].get_left_sign();



       /* auto Poly_at_x = [&derv, i](long double x) {
            long double val = 0;
            for (int k = 0; k <= derv[i - 1].Degree(); k++) {
                val += derv[i - 1][k] * pow(x, derv[i - 1].Degree() - k);
            }
            return val;
        };*/

        auto value_sign = [](long double x) {
            if (x > 0) {
                return 1;
            }
            return -1;
        };


        value_ex = derv[i - 1].eval(left_ex);
        init_step = -0.1;
        if (fabs(value_ex) < epsilon) {
            //this is a solution!
            sols_poly.push_back(left_ex);
        }
        else {
            sign_ex = value_sign(value_ex);
            long double x = left_ex;
            if (sign_ex != left_limit_sign) {
                do {
                    x += init_step;
                    value_x = derv[i - 1].eval(x);
                    sign_x = value_sign(value_x);
                } while (sign_x == sign_ex);
                a = x;
                b = left_ex;
                sol = bisection_n_newtoon(derv[i - 1], derv[i], a, b);
                sols_poly.push_back(sol);
            }

        }



        //middle
        for (int temp = 0; temp < number_of_sols - 1; temp++) {
            //i,i+1
            a = sols_derv[temp];
            b = sols_derv[temp + 1];
            val_a = derv[i - 1].eval(a);
            val_b = derv[i - 1].eval(b);
            if (fabs(val_a) >= epsilon) {
                if (fabs(val_b) <= epsilon) {
                    sols_poly.push_back(b);
                }
                else {
                    sign_a = value_sign(val_a);
                    sign_b = value_sign(val_b);
                    if (sign_a != sign_b) {
                        sol = bisection_n_newtoon(derv[i - 1], derv[i], a, b);
                        sols_poly.push_back(sol);

                    }
                }
            }
        }

        //last

        value_ex = derv[i - 1].eval(right_ex);
        init_step = 0.1;
        if (fabs(value_ex) < epsilon) {
            //this is a solution!
            sols_poly.push_back(right_ex);
        }
        else {
            sign_ex = value_sign(value_ex);
            long double x = right_ex;
            if (sign_ex != right_limit_sign) {
                do {
                    x += init_step;
                    value_x = derv[i - 1].eval(x);
                    sign_x = value_sign(value_x);
                } while (sign_x == sign_ex);
                a = right_ex;
                b = x;
                sol = bisection_n_newtoon(derv[i - 1], derv[i], a, b);
                sols_poly.push_back(sol);
            }

        }
        vector<long double>().swap(sols_derv);
        number_of_sols = sols_poly.size();
        for (int j = 0; j < sols_poly.size(); j++) {
            sols_derv.push_back(sols_poly[j]);
        }

    }

    return sols_poly;

}


long double bisection_n_newtoon(Poly& pol, Poly& derv, long double a, long double b) {
    int i;
    long double x, f, fa, fb, m, fm, fdiff_m, old_a, old_b, next_x_newton;

    int n_bsic = 0, n_newt = 0, iter = 0;

    /*auto Poly_at_x = [&pol](long double x) {
        long double val = 0;
        for (int k = 0; k <= pol.Degree(); k++) {
            val += pol[k] * pow(x, pol.Degree() - k);
        }
        return val;
    };

    auto derv_at_x = [&derv](long double x) {
        long double val = 0;
        for (int k = 0; k <= derv.Degree(); k++) {
            val += derv[k] * pow(x, derv.Degree() - k);
        }
        return val;
    };*/

    fa = pol.eval(a);
    fb = pol.eval(b);
    m = (a + b) / 2.0;
    fm = pol.eval(m);
    while (fabs(fm) > epsilon) {
        old_a = a;
        old_b = b;
        iter++;
        fm = pol.eval(m);
        fdiff_m = derv.eval(m);
        next_x_newton = m - fm / fdiff_m;
        if (next_x_newton >= a && next_x_newton <= b) {
            m = next_x_newton;
            n_newt++;
        }
        else {
            m = (a + b) / 2.0;
            n_bsic++;
        }
        fm = pol.eval(m);;
        if (fa * fm < 0) {
            b = m;
            fb = fm;
        }
        if (fb * fm < 0) {
            a = m;
            fa = fm;
        }
        if (a == old_a && b == old_b) {
            break;
        }
    }
    return m;
}