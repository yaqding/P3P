#pragma once

#include <new/solve_cubic_new.h>
#include <new/solve_eig0_new.h>
#include <new/refine_lambda_new.h>
#include <new/p3p_timers_new.h>

using std::cout;
using std::endl;

namespace cvl
{

    template <class T, int refinement_iterations = 5>
    int p3p_new(Vector3<T> y1,
                Vector3<T> y2,
                Vector3<T> y3,
                Vector3<T> x1,
                Vector3<T> x2,
                Vector3<T> x3,
                Vector<cvl::Matrix<T, 3, 3>, 4> &Rs,
                Vector<Vector3<T>, 4> &Ts)
    {

        // normalize the length of ys, we could expect it, but lets not...
        TIC(lt1);
        y1.normalize();
        y2.normalize();
        y3.normalize();

        T b12 = -2.0 * (y1.dot(y2));
        T b13 = -2.0 * (y1.dot(y3));
        T b23 = -2.0 * (y2.dot(y3));

        // implicit creation of Vector3<T> can be removed
        Vector3<T> d12 = x1 - x2;
        Vector3<T> d13 = x1 - x3;
        Vector3<T> d23 = x2 - x3;
        Vector3<T> d12xd13(d12.cross(d13));

        T a12 = d12.squaredLength();
        T a13 = d13.squaredLength();
        T a23 = d23.squaredLength();

        // if(abs(D1.determinant())<1e-5 || fabs(D2.determinant())<1e-5)        cout<<"det(D): "<<D1.determinant()<<" "<<D2.determinant()<<endl;
        TOC(lt1);
        TIC(lt2);

        // a*g^3 + b*g^2 + c*g + d = 0
        T c31 = -0.5 * b13;
        T c23 = -0.5 * b23;
        T c12 = -0.5 * b12;
        T blob = (c12 * c23 * c31 - 1.0);

        T s31_squared = 1.0 - c31 * c31;
        T s23_squared = 1.0 - c23 * c23;
        T s12_squared = 1.0 - c12 * c12;

        T p3 = (a13 * (a23 * s31_squared - a13 * s23_squared));

        T p2 = 2.0 * blob * a23 * a13 + a13 * (2.0 * a12 + a13) * s23_squared + a23 * (a23 - a12) * s31_squared;

        T p1 = a23 * (a13 - a23) * s12_squared - a12 * a12 * s23_squared - 2.0 * a12 * (blob * a23 + a13 * s23_squared);

        T p0 = a12 * (a12 * s23_squared - a23 * s12_squared);

        TOC(lt2);
        TIC(lt3);
        T g = 0;

        p3 = 1.0 / p3;
        p2 *= p3;
        p1 *= p3;
        p0 *= p3;

        T a = p1 - p2 * p2 / 3.0;
        T b = (2.0 * p2 * p2 * p2 - 9.0 * p2 * p1) / 27.0 + p0;
        T dc = b * b / 4.0 + a * a * a / 27.0;
        T threshold = 1.0e-30;
        if (dc > threshold)
        {
            T cc = std::sqrt(dc);
            b *= -0.5;
            g = std::cbrt(b + cc) + std::cbrt(b - cc) - p2 / 3.0;
        }
        else if (dc < -threshold)
        {
            T cc = 3.0 * b / (2.0 * a) * std::sqrt(-3.0 / a);
            g = 2.0 * std::sqrt(-a / 3.0) * std::cos(std::acos(cc) / 3.0) - p2 / 3.0;
        }
        else if (dc <= threshold && dc >= -threshold && a!=0)
        {
            g = 3.0 * b / a - p2 / 3.0;
        }
        else
        {
            g = - p2 / 3.0;
        }

        // g = cubick_new(p2, p1, p0);

        TOC(lt3);
        // we can swap D1,D2 and the coeffs!
        // oki, Ds are:
        // D1=M12*XtX(2,2) - M23*XtX(1,1);
        // D2=M23*XtX(3,3) - M13*XtX(2,2);

        //[    a23 - a23*g,                 (a23*b12)/2,              -(a23*b13*g)/2]
        //[    (a23*b12)/2,           a23 - a12 + a13*g, (a13*b23*g)/2 - (a12*b23)/2]
        //[ -(a23*b13*g)/2, (a13*b23*g)/2 - (a12*b23)/2,         g*(a13 - a23) - a12]
        TIC(lt4);

        // gain 13 ns...
        T A00 = a23 * (1.0 - g);
        T A01 = (a23 * b12) * 0.5;
        T A02 = (a23 * b13 * g) * (-0.5);
        T A11 = a23 - a12 + a13 * g;
        T A12 = b23 * (a13 * g - a12) * 0.5;
        T A22 = g * (a13 - a23) - a12;


        // Matrix<T, 2, 2> Aiv(A01, A02,
        //                   A11, A12);
        // Vector2d Ab(-A00,-A01);
        // Vector2d yz = Aiv.inverse()*Ab;
        // Vector3<T> Bp;
        // Bp(0) = 1; Bp(1) = yz(0); Bp(2) = yz(1);
        // Bp.normalize();
        // Bp = Bp*std::sqrt(A01*A01+A02*A02+A12*A12-A00*A11-A00*A22-A22*A11);

        Matrix<T, 3, 2> lm;

        if (std::abs(A00) >= std::abs(A11) && std::abs(A00) > std::abs(A22))
        {
            T co0 = 1.0/A00;
            T co1 = A01*co0;
            T co2 = A02*co0;
            T co3 = A11*co0;
            T co4 = A12*co0;
            T del = std::sqrt(co1*co1-co3);

            lm(0, 0) = 1.0;
            lm(1, 0) = co1+del;
            lm(2, 0) = (lm(1, 0)*co2-co4)/del;
            lm(0, 1) = 1.0;
            lm(1, 1) = co1-del;
            lm(2, 1) = 2*co2-lm(2, 0);

        }
        else if (std::abs(A11) >= std::abs(A00) && std::abs(A11) > std::abs(A22))
        {
            T co0 = 1.0/A11;
            T co1 = A01*co0;
            T co2 = A12*co0;
            T co3 = A00*co0;
            T co4 = A02*co0;
            T del = std::sqrt(co1*co1-co3);

            lm(0, 0) = co1+del;
            lm(1, 0) = 1.0;
            lm(2, 0) = (lm(0, 0)*co2-co4)/del;
            lm(0, 1) = co1-del;
            lm(1, 1) = 1.0;
            lm(2, 1) = 2*co2-lm(2, 0);

        }
        else
        {
            T co0 = 1.0/A22;
            T co1 = A02*co0;
            T co2 = A12*co0;
            T co3 = A00*co0;
            T co4 = A01*co0;
            T del = std::sqrt(co1*co1-co3);

            lm(0, 0) = co1+del;
            lm(1, 0) = (lm(0, 0)*co2-co4)/del;
            lm(2, 0) = 1.0;
            lm(0, 1) = co1-del;
            lm(1, 1) = 2*co2-lm(1, 0);
            lm(2, 1) = 1.0;
        }


        // Matrix<T, 3, 3> A(A00, A01, A02,
        //                   A01, A11, A12,
        //                   A02, A12, A22);

        // T B00 = -A11 * A22 + A12 * A12;
        // T B01 = -A12 * A02 + A01 * A22;
        // T B02 = -A01 * A12 + A11 * A02;
        // T B11 = -A00 * A22 + A02 * A02;
        // T B12 = -A01 * A02 + A00 * A12;
        // T B22 = -A00 * A11 + A01 * A01;

        // Vector3<T> Bp;
        // T Bm;
        // if (B00 >= B11 && B00 >= B22)
        // {
        //     Bm = 1.0 / std::sqrt(B00);
        //     Bp(0) = B00 * Bm;
        //     Bp(1) = B01 * Bm;
        //     Bp(2) = B02 * Bm;
        // }
        // else if (B11 > B00 && B11 > B22)
        // {
        //     Bm = 1.0 / std::sqrt(B11);
        //     Bp(0) = B01 * Bm;
        //     Bp(1) = B11 * Bm;
        //     Bp(2) = B12 * Bm;
        // }
        // else
        // {
        //     Bm = 1.0 / std::sqrt(B22);
        //     Bp(0) = B02 * Bm;
        //     Bp(1) = B12 * Bm;
        //     Bp(2) = B22 * Bm;
        // }

        // Matrix<T, 3, 3> Mp;
        // Mp = Bp.crossMatrix();

        // Matrix<T, 3, 3> Cc;
        // Cc = A + Mp;

        // Matrix<T, 3, 2> lm;

        // lm(0, 0) = Cc(0, 0);
        // lm(1, 0) = Cc(1, 0);
        // lm(2, 0) = Cc(2, 0);
        // lm(0, 1) = Cc(0, 0);
        // lm(1, 1) = Cc(0, 1);
        // lm(2, 1) = Cc(0, 2);


        // TOC(lt4);
        // TIC(lt5);

        int valid = 0;
        Vector<Vector<T, 3>, 4> Ls;

        // use the t=Vl with t2,st2,t3 and solve for t3 in t2

        // Vector2<T> ss(v,-v);
        for (int i = 0; i < 2; ++i)
        {

            // u = V(:, 1) - sV(:,2)

            

            T u1 = lm(0, i);
            T u2 = lm(1, i);
            T u3 = lm(2, i);

            // we are computing lambda using a linear relation
            // u1*l1 + u2*l2 + u3*l3=0
            // li>0, implies all three ui=0 is degenerate...
            // if two are zero the third must be
            // hence at most one can be zero.
            // divide by the largest for best numerics,
            // simple version, use the bigger of u1, u2, one will always be non-zero
            if (std::abs(u1) < std::abs(u2))
            {
                // solve for l2
                // T a= (a23 - a12)*u3*u3 - a12*u2*u2 + a12*b23*u2*u3;
                // T b= (T(2)*a23*u1*u3 - T(2)*a12*u1*u3 + a12*b23*u1*u2 - a23*b12*u2*u3)/a;
                // T c= (a23*u1*u1 - a12*u1*u1 + a23*u2*u2 - a23*b12*u1*u2)/a;

                T u33 = u3 * u3,
                  u22 = u2 * u2,
                  u23 = u2 * u3,
                  u13 = u1 * u3,
                  u12 = u1 * u2,
                  u11 = u1 * u1,
                  a12b23 = a12 * b23,
                  a23b12 = a23 * b12;
                T a = 1.0 / ((a23 - a12) * u33 - a12 * u22 + a12b23 * u23);
                T b = ((T(2) * a23 - T(2) * a12) * u13 + a12b23 * u12 - a23b12 * u23) * a;
                T c = ((a23 - a12) * u11 + a23 * u22 - a23b12 * u12) * a;

                Vector2<T> taus;
                if (!root2real(b, c, taus[0], taus[1]))
                    continue;
                for (T tau : taus)
                {
                    if (tau <= 0)
                        continue;
                    //(tau^2 + b13*tau + 1)*l1^2 = a13
                    // positive only
                    T l1 = std::sqrt(a13 / (tau * (tau + b13) + T(1.0)));
                    T l3 = tau * l1;
                    T l2 = -(u1 * l1 + u3 * l3) / u2;
                    if (l2 <= 0)
                        continue;
                    Ls[valid] = {l1, l2, l3};
                    ++valid;
                }
            }
            else
            { // solve for l1

                T w2 = T(1.0) / (-u1);
                T w0 = u2 * w2;
                T w1 = u3 * w2;

                T a = T(1.0) / ((a13 - a12) * w1 * w1 - a12 * b13 * w1 - a12);
                T b = (a13 * b12 * w1 - a12 * b13 * w0 - T(2.0) * w0 * w1 * (a12 - a13)) * a;
                T c = ((a13 - a12) * w0 * w0 + a13 * b12 * w0 + a13) * a;

                Vector2<T> taus;
                if (!root2real(b, c, taus[0], taus[1]))
                    continue;
                for (T tau : taus)
                {
                    if (tau <= 0)
                        continue;
                    T d = a23 / (tau * (b23 + tau) + T(1.0));
                    T l2 = std::sqrt(d);
                    T l3 = tau * l2;
                    T l1 = w0 * l2 + w1 * l3;
                    if (l1 <= 0)
                        continue;

                    Ls[valid] = {l1, l2, l3};
                    ++valid;
                }
            }

            // if (c < -1.0e-26 && )
            if (valid>0 && dc > 0.000000001) 
                // std::cout << dc <<" :" << valid <<std::endl;
                break;

            // if(valid>0 && sigma_p<0) break;
        }







        TOC(lt5);
        TIC(lt6);

        for (int i = 0; i < valid; ++i)
        {
            gauss_newton_refineL_new<T, refinement_iterations>(Ls[i], a12, a13, a23, b12, b13, b23);
        }
        TOC(lt6);

        TIC(lt8);
        Vector3<T> ry1, ry2, ry3;
        Vector3<T> yd1;
        Vector3<T> yd2;
        Vector3<T> yd1xd2;
        Matrix<T, 3, 3> X(d12(0), d13(0), d12xd13(0),
                          d12(1), d13(1), d12xd13(1),
                          d12(2), d13(2), d12xd13(2));
        // guaranteed to exist due to cross, and non degen
        X = X.inverse();

        for (int i = 0; i < valid; ++i)
        {
            // cout<<"Li="<<Ls(i)<<endl;

            // compute the rotation:
            ry1 = y1 * Ls(i)(0);
            ry2 = y2 * Ls(i)(1);
            ry3 = y3 * Ls(i)(2);

            yd1 = ry1 - ry2;
            yd2 = ry1 - ry3;
            yd1xd2 = yd1.cross(yd2);

            Matrix<T, 3, 3> Y(yd1(0), yd2(0), yd1xd2(0),
                              yd1(1), yd2(1), yd1xd2(1),
                              yd1(2), yd2(2), yd1xd2(2));

            Rs[i] = Y * X;

            Ts[i] = (ry1 - Rs[i] * x1);
        }

        TOC(lt8);
        return valid;
    }
}
