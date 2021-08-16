#include "FindRoots.h";
#include <algorithm>
#include <thread>
#define epsilon 1.0e-140
//#define epsilon 0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001

void external_root(Poly pol,Poly diff,long double extreme_p, long double limit_sign,long double init_step, vector<long double>& sols_poly) {
    long double value_ex, sign_ex, value_x, sign_x, a, b, sol;
    //vector<long double> sols_poly;
    auto value_sign = [](long double x) {
        if (x > 0) {
            return 1;
        }
        return -1;
    };

    value_ex = pol.eval(extreme_p);
    if (fabs(value_ex) <= epsilon) {
        sols_poly.emplace_back(extreme_p);
    }
    else {
        sign_ex = value_sign(value_ex);
        long double x = extreme_p;
        if (sign_ex != limit_sign) {
            do {
                x += init_step;
                //init_step = init_step * 2.0;
                value_x = pol.eval(x);
                sign_x = value_sign(value_x);
            } while (sign_x == sign_ex);
            if (init_step > 0) {
                a = extreme_p;
                b = x;
            }
            else {
                a = x;
                b = extreme_p;
            }
            sols_poly.emplace_back(bisection_n_newtoon(pol, diff, a, b));
        }
    }
}

void internal_roots(Poly pol, Poly diff, vector<long double> diff_sols,vector<long double>& sols_poly) {
    long double a, b, val_a, val_b, sign_a, sign_b, sol;
    //vector<long double> sols_poly;

    auto value_sign = [](long double x) {
        if (x > 0) {
            return 1;
        }
        return -1;
    };

    for (int temp = 0; temp < diff_sols.size() - 1; temp++) {
        //i,i+1
        a = diff_sols[temp];
        b = diff_sols[temp + 1];
        val_a = pol.eval(a);
        val_b = pol.eval(b);
        if (fabs(val_a) >= epsilon) {
            if (fabs(val_b) <= epsilon) {
                sols_poly.emplace_back(b);
            }
            else {
                sign_a = value_sign(val_a);
                sign_b = value_sign(val_b);
                if (sign_a != sign_b) {
                    sol = bisection_n_newtoon(pol, diff, a, b);
                    sols_poly.emplace_back(sol);

                }
            }
        }
    }
}

void loading_bar(float& progress,int barWidth,int n) {
    cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) cout << (char)254u;
        else cout << " ";
    }
    cout << "] " << int(progress * 100.0) << " %\r";
    cout.flush();
    progress += (1 / static_cast<float>(n));
}

vector<long double> Solve(Poly Pol_solve,bool show_progress) {

    vector<long double> coff = Pol_solve.getdata();
    auto find_max_fabs = [](vector<long double> m) {
        return *max_element(m.begin(), m.end(), [](long double& a, long double& b){
                return fabs(a) < fabs(b); });
    };
    auto find_max = [](vector<long double> m) {
        return *max_element(m.begin(), m.end());
    };

    Poly pol = Pol_solve / (find_max_fabs(coff)* 9999999999999999999.0);
    cout << "Calculating Poly diff's..." << endl;
    vector<Poly> derv;
    derv.reserve(Pol_solve.Degree()+1);
    derv.push_back(pol);
    int degree = Pol_solve.Degree();
    int index = 0;
    float progress = 0.0;
    while (degree > 1) {
        vector<long double> dervcoff(derv.at(index).getdata().begin(), derv.at(index).getdata().end() - 1);
        int sizeofup = derv.at(index).Degree();
        transform(dervcoff.begin(), dervcoff.end(), dervcoff.begin(), [&sizeofup](auto& c) {return  c * (sizeofup--); });
        long double max = find_max_fabs(dervcoff)* 9999999999999999999.0;
        for (long double& mem : dervcoff)
            mem /= max;
            //mem /= static_cast<long double>(fabs(max+1));
        //Poly pol(dervcoff, dervcoff.size()-1);
        //derv.push_back(pol);
        derv.emplace_back(dervcoff, dervcoff.size() - 1);
        if(show_progress)
            loading_bar(progress, 70, Pol_solve.Degree()-1);
        index++;
        degree--;
    }
    if (show_progress) {
        loading_bar(progress, 70, Pol_solve.Degree() - 1);
        cout << endl;
    }
    int stack_size = derv.size();

    vector<long double> sols_derv = { -derv[stack_size - 1][1] / derv[stack_size - 1][0] };
    vector<long double> sols_poly;
    sols_poly.reserve(10);
    int number_of_sols = sols_derv.size();

    long double left_ex, right_ex, right_limit_sign, left_limit_sign;
    auto value_sign = [](long double x) {
        if (x > 0) {
            return 1;
        }
        return -1;
    };
    vector<long double> external_left,internal, external_right;
    vector<long double> thread_sols_mid, thread_left_sols, thread_right_sols;
    thread thread_1, thread_2, thread_3;
    cout << "Finding Roots of the Poly..." << endl;
    progress = 0.0;
    for (int i = stack_size - 1; i > 0; i--) {
        if (number_of_sols == 0) {
            left_ex = right_ex = 0;
        }
        else {
            left_ex = sols_derv[0];
            right_ex = sols_derv[number_of_sols - 1];
        }

        right_limit_sign = derv[i - 1].get_right_sign();
        left_limit_sign = derv[i - 1].get_left_sign();

        if(show_progress)
            loading_bar(progress, 70, stack_size-1);
   
        vector<long double>().swap(sols_poly);
        sols_poly.reserve(20);
        long double step = 0.001;
   

        vector<long double>().swap(thread_left_sols);
        vector<long double>().swap(thread_sols_mid);
        vector<long double>().swap(thread_right_sols);

        thread_1 = thread(external_root, derv[i - 1], derv[i], left_ex, left_limit_sign, -step, std::ref(thread_left_sols));
        
        

        //left external root
        /*external_left = external_root(derv[i - 1], derv[i], left_ex, left_limit_sign, -step);
        sols_poly.insert(sols_poly.end(), external_left.begin(), external_left.end());*/

        
        thread_2 = thread(internal_roots, derv[i - 1], derv[i], sols_derv, std::ref(thread_sols_mid));
        
        

        //internal roots
        //internal = internal_roots(derv[i - 1], derv[i], sols_derv);
        //sols_poly.insert(sols_poly.end(), internal.begin(), internal.end());

        thread_3 = thread(external_root, derv[i - 1], derv[i], right_ex, right_limit_sign, step, std::ref(thread_right_sols));
        
       

        //right external root
       /* external_right = external_root(derv[i - 1], derv[i], right_ex, right_limit_sign, step);
        sols_poly.insert(sols_poly.end(), external_right.begin(), external_right.end());*/


        thread_1.join();
        thread_2.join();
        thread_3.join();
        sols_poly.insert(sols_poly.end(), thread_left_sols.begin(), thread_left_sols.end());
        sols_poly.insert(sols_poly.end(), thread_sols_mid.begin(), thread_sols_mid.end());
        sols_poly.insert(sols_poly.end(), thread_right_sols.begin(), thread_right_sols.end());
        

        vector<long double>().swap(sols_derv);
        number_of_sols = sols_poly.size();
        for (int j = 0; j < sols_poly.size(); j++) {
            sols_derv.emplace_back(sols_poly[j]);
        }
        

    }
    if (show_progress)
        loading_bar(progress, 70, stack_size - 1);
    cout << endl;
    return sols_poly;

}


long double bisection_n_newtoon(Poly pol, Poly derv, long double a, long double b) {
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
    while (fabs(fm) > epsilon && iter < 100) {
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
        //cout << " a = " << a << " , b= " << b << ", fa = " << fa << ", fb = " << fb << ", m = " << m << ", fm = " << fm << endl;
    }
    return m;
}