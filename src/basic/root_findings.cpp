/**
 * \file root_findings.cpp
 * \brief Root-finding functions of polynomials defined by std::valarray.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-01-29
 */


# include <algorithm>
# include <numeric>
# include <cmath>

# include "basic.h"
# include "exceptions.h"


namespace simpoly
{
namespace basic
{

CArry low_degree_roots0(const CArry &P) { return CArry(0); }

CArry low_degree_roots1(const CArry &P) { return {-P[0]/P[1]}; }

CArry low_degree_roots2(const CArry &P)
{
    Cmplx twoA = 2.0 * P[2];
    Cmplx sqFourAC = std::sqrt(4.0 * P[2] * P[0]);
    return {(-P[1]+sqFourAC)/twoA, (-P[1]-sqFourAC)/twoA};
}

CArry use_low_degree_formula(const CArry &P)
{
    switch (P.size())
    {
        case 1: // degree 0 (i.e., constant)
            return low_degree_roots0(P);
            break;
        case 2: // degree 1 (i.e., linear)
            return low_degree_roots1(P);
            break;
        case 3: // degree 2 (i.e., quadratic)
            return low_degree_roots2(P);
            break;
        case 4: // degree 2 (i.e., cubic)
        case 5: // degree 2 (i.e., quartic)
       default:
           throw exceptions::PolynomialErrorGeneral(
                   __FL__, "Not implemented yet.");
   }
}

template <typename T>
T newton_raphson(const Arry<T> &coeffs, const T guess, const double tol)
{
    CHECK_COEFS(coeffs, 1e-12);

    T ans = guess;
    Arry<T> d = derivative(coeffs);

    for(unsigned iter=0; iter<10000; ++iter)
    {
        T value = evaluate(coeffs, ans),
          d_value = evaluate(d, ans),
          diff;

        if (value == 0.0) break; // found root; note we use exact 0.0 here

        if (d_value == 0.0) // to avoid devide by zero
        {
            ans *= 1.0001;
            value = evaluate(coeffs, ans);
            d_value = evaluate(d, ans);
        }
        diff = value / d_value;
        ans -= diff;
        if (std::abs(diff)/std::abs(((ans==0.0)?1.0:ans)) < tol) break;
        if (iter > 10000) throw exceptions::InfLoop(__FILE__, __LINE__);
    }

    return ans;
}


CArry aberth(const CArry &coeffs, const CArry &guess, const double tol)
{
    CHECK_COEFS(coeffs, 1e-12);

# ifdef NDEBUG
    if (guess.size() == 0) throw exceptions::PolynomialErrorGeneral(
            __FL__, "The length of initial guess can not be zero.");
# endif

    // alias to the length of provided coefficient array
    const auto &len = coeffs.size();

    // use exact solution for low-degree polynomials
    if (len < 4) return use_low_degree_formula(coeffs);

    // initialize initial guess through copying
    CArry rts(guess);

    // create a vector indicating if each root has been found
    std::vector<bool> stop(guess.size(), false);

    // derivative
    CArry d = derivative(coeffs);

    // an index to record the number of while iteration
    long iter = 0;

    while (std::any_of(stop.begin(), stop.end(), [](bool b){return !b;}))
    {
        for(unsigned int i=0; i<guess.size(); ++i)
        {
            const Cmplx &zi = rts[i]; // alias
            Cmplx temp, delta;

            temp = evaluate(coeffs, zi) / evaluate(d, zi);

            delta = std::accumulate(
                rts.begin(), rts.begin()+i, Cmplx(0.0, 0.0),
                [&zi] (const Cmplx &x, const Cmplx &y) -> Cmplx
                    { return x + 1.0 / (zi - y); });

            delta += std::accumulate(
                rts.begin()+i+1, rts.end(), Cmplx(0.0, 0.0),
                [&zi] (const Cmplx &x, const Cmplx &y) -> Cmplx
                    { return x + 1.0 / (zi - y); });

            delta *= temp;
            delta = 1.0 - delta;
            delta = temp / delta;

            if ((std::abs(delta)/std::abs(rts[i])) < tol) stop[i] = true;

            rts[i] -= delta;
            if (evaluate(coeffs, rts[i]) == 0.0) stop[i] = true; // exact zero
        }

        // check the number of iterations
        iter += 1;
        if (iter > 10000) throw exceptions::InfLoop(__FILE__, __LINE__);
    }

    return rts;
}


CArry aberth(const CArry &coeffs, const DArry &guess, const double tol)
{
    CHECK_COEFS(coeffs, 1e-12);

    // use exact solution for low-degree polynomials
    if (coeffs.size() < 4) return use_low_degree_formula(coeffs);

    // make a copy of guess with complex type, and perturbation in imag
    CArry G = to_CArry(guess);

    // add perturbation to imaginary parts
    std::for_each(std::begin(G), std::end(G), [&tol](Cmplx &i){i.imag(tol);});

    return aberth(coeffs, G, tol);
}


CArry aberth(const DArry &coeffs, const CArry &guess, const double tol)
{
    CHECK_COEFS(coeffs, 1e-12);

    // use exact solution for low-degree polynomials
    if (coeffs.size() < 4) return use_low_degree_formula(to_CArry(coeffs));

    return aberth(to_CArry(coeffs), guess, tol);
}


CArry aberth(const DArry &coeffs, const DArry &guess, const double tol)
{
    CHECK_COEFS(coeffs, 1e-12);

    // use exact solution for low-degree polynomials
    if (coeffs.size() < 4) return use_low_degree_formula(to_CArry(coeffs));

    // make a copy of guess with complex type
    CArry G = to_CArry(guess);

    // add perturbation to imaginary parts
    std::for_each(std::begin(G), std::end(G), [&tol](Cmplx &i){i.imag(tol);});

    return aberth(coeffs, G, tol);
}


CArry aberth(const CArry &coeffs, const double tol)
{
    CHECK_COEFS(coeffs, 1e-12);

    // use exact solution for low-degree polynomials
    if (coeffs.size() < 4) return use_low_degree_formula(coeffs);

    // initialize initial guess through copying
    CArry guess(coeffs.size()-1);

    const auto &bg = guess.data();
    std::for_each(guess.begin(), guess.end(),
        [&bg](Cmplx &x){x = std::pow(Cmplx(0.5, 0.5), double(&x-bg));});

    return aberth(coeffs, guess, tol);
}


CArry aberth(const DArry &coeffs, const double tol)
{
    CHECK_COEFS(coeffs, 1e-12);

    // use exact solution for low-degree polynomials
    if (coeffs.size() < 4) return use_low_degree_formula(to_CArry(coeffs));

    // initialize initial guess through copying
    CArry guess(coeffs.size()-1);

    const auto &bg = guess.data();
    std::for_each(guess.begin(), guess.end(),
        [&bg](Cmplx &x){x = std::pow(Cmplx(0.5, 0.5), double(&x-bg));});

    return aberth(coeffs, guess, tol);
}


CArry yan_and_chieng_2006(const CArry &coeffs, const double tol)
{
    CHECK_COEFS(coeffs, 1e-12);

    // use exact solution for low-degree polynomials
    if (coeffs.size() < 4) return use_low_degree_formula(coeffs);

    CArry     drv = derivative(coeffs); // direvative
    CArry     agcd = GCD(coeffs, drv); // approximated GCD
    CArry     q_coeffs = divide(coeffs, agcd); // f(x) / GCD
    CArry     q_drv = divide(drv, agcd); // direvative / GCD
    CArry     drv_q_coeffs = derivative(q_coeffs);

    // simple roots from q_coeffs
    CArry       simples1 = aberth(q_coeffs, tol);

    // refine simple roots with original polynomial `coeffs`
    CArry       simples2 = aberth(coeffs, simples1, tol);

    // an array for final results
    CArry       result(coeffs.size() - 1);

    unsigned k = 0; // index for result
    for(unsigned i=0; i<simples1.size(); ++i)
    {
        // claculate multiplicity
        unsigned m = int((evaluate(q_drv, simples1[i]) /
                evaluate(drv_q_coeffs, simples1[i])).real() + 0.5);

        // choose the best one
        if (m > 1)
        {
            double e1 = std::norm(evaluate(drv, simples1[i])) +
                        std::norm(evaluate(coeffs, simples1[i])),
                   e2 = std::norm(evaluate(drv, simples2[i])) +
                        std::norm(evaluate(coeffs, simples2[i]));

            if (e1 < e2) result[k] = simples1[i];
            else result[k] = simples2[i];
            k += 1;

            // duplicate multiple roots
            for(unsigned mi=1; mi<m; ++mi)
            {
                result[k] = result[k-1];
                k += 1;
            }
        }
        else
        {
            double e1 = std::abs(evaluate(coeffs, simples1[i])),
                   e2 = std::abs(evaluate(coeffs, simples2[i]));

            if (e1 < e2) result[k] = simples1[i];
            else result[k] = simples2[i];
            k += 1;
        }
    }

    return result;
}


CArry yan_and_chieng_2006(const DArry &coeffs, const double tol)
{
    CHECK_COEFS(coeffs, 1e-12);

    // use exact solution for low-degree polynomials
    if (coeffs.size() < 4) return use_low_degree_formula(to_CArry(coeffs));

    CArry C = to_CArry(coeffs);
    return yan_and_chieng_2006(C, tol);
}


// explicit instantiation
template double newton_raphson(const DArry &coeffs, const double guess, const double tol);
template Cmplx newton_raphson(const CArry &coeffs, const Cmplx guess, const double tol);

} // end of namespace basic
} // end of namespace simpoly
