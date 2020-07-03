#pragma once

#include <linear/vec.hpp>
#include <lie/exp_utils.hpp>

#define _USE_MATH_DEFINES
#include <math.h>

namespace sight
{
    // Represents a group element of SO(3)
    // which is a rotation matrix of det +1.
    template <typename S>
    class SO3
    {
    public:

        // Default constructor leaves R invalid
        SO3() {}
        
        SO3(std::initializer_list<S> l)
        {
            std::copy(l.begin(), l.end(), &R[0]);
        }

        SO3 operator*(const SO3& s) const
        {
            SO3 res;
            res.R[0] = R[0] * s.R[0] + R[1] * s.R[3] + R[2] * s.R[6];
            res.R[1] = R[0] * s.R[1] + R[1] * s.R[4] + R[2] * s.R[7];
            res.R[2] = R[0] * s.R[2] + R[1] * s.R[5] + R[2] * s.R[8];

            res.R[3] = R[3] * s.R[0] + R[4] * s.R[3] + R[5] * s.R[6];
            res.R[4] = R[3] * s.R[1] + R[4] * s.R[4] + R[5] * s.R[7];
            res.R[5] = R[3] * s.R[2] + R[4] * s.R[5] + R[5] * s.R[8];

            res.R[6] = R[6] * s.R[0] + R[7] * s.R[3] + R[8] * s.R[6];
            res.R[7] = R[6] * s.R[1] + R[7] * s.R[4] + R[8] * s.R[7];
            res.R[8] = R[6] * s.R[2] + R[7] * s.R[5] + R[8] * s.R[8];

            return res;
        }

        Vec3<S> operator*(const Vec3<S>& v) const
        {
            Vec3<S> Rv;
            Rv[0] = R[0] * v[0] + R[1] * v[1] + R[2] * v[2];
            Rv[1] = R[3] * v[0] + R[4] * v[1] + R[5] * v[2];
            Rv[2] = R[6] * v[0] + R[7] * v[1] + R[8] * v[2];

            return Rv;
        }

        SO3 Inverse() const
        {
            SO3 inv;
            inv.R[0] = R[0];
            inv.R[1] = R[3];
            inv.R[2] = R[6];

            inv.R[3] = R[1];
            inv.R[4] = R[4];
            inv.R[5] = R[7];

            inv.R[6] = R[2];
            inv.R[7] = R[5];
            inv.R[8] = R[8];

            return inv;
        }

        void SetIdentity()
        {
            R[0] = S(1);
            R[1] = S(0);
            R[2] = S(0);

            R[3] = S(0);
            R[4] = S(1);
            R[5] = S(0);

            R[6] = S(0);
            R[7] = S(0);
            R[8] = S(1);
        }

        Vec3<S> Log() const
        {
            // Since R = I + a[w] + b[w]^2
            //
            // We have:
            //       R - R^T = 2a[w]
            //
            // Since I and b[w]^2 are symmetric, but [w] is skew-symmetric.
            //
            // With a = sin(t),
            //
            //      [w] = (R - R^T) / (2sin(t))
            //
            // The trace of a rotation matrix is equal to the sum of its
            // eigenvalues. With this we can show that:
            //
            //     Trace(R) = 1 + 2 * cos(t)
            //
            // Thus,
            //
            //     cos(t) = (Trace(R) - 1) / 2

            // Initialize v to be aw, where a = sin(t)
            Vec3<S> v = {
                (R[7] - R[5]) / S(2),
                (R[2] - R[6]) / S(2),
                (R[3] - R[1]) / S(2)
            };

            const S trace = R[0] + R[4] + R[8];
            const S cost = (trace - 1) / S(2);
            const S sint_sq = v.SquaredNorm();

            // Small values of theta.
            //
            // We could use cos(x) where x is the constant
            // threshold used for exponentiation. However,
            // it comes very close to 1 for both floats and
            // doubles. So we just choose something close
            // to 1 instead.
            if (cost > S(.999))
            {
                // Currently, v = sin(t) * w
                //
                // We desire, v = tw, so we have to find a constant,
                // c, such that c * sin(t) = t;
                //
                // So, c = t / sin(t)
                //
                // We can calculate c = arcsin(sin(x)) / sin(x). So we
                // calculate the Taylor Series approximation at x = 0 of:
                //
                //     f(x) = arcsin(x) / x
                const S t_sint = Arcsin_Taylor(sint_sq);
                
                v[0] *= t_sint;
                v[1] *= t_sint;
                v[2] *= t_sint;

                return v;
            }

            // Normal values of theta.
            if (cost > S(-.999))
            {
                const S t = acos(cost);
                const S sint = sqrt(sint_sq);
                const S scale = t / sint;

                v[0] *= scale;
                v[1] *= scale;
                v[2] *= scale;

                return v;
            }

            // Thetas close to pi.
            //
            // As theta approaches pi, sin(t)/t approaches 0, which
            // can make our approach of leveraging R - R^T erroneous.
            //
            // To avoid the cancellation issues, we look at R + R^T:
            //
            // R + R^T = 2(I + b[w]^2)
            //
            // b = (1 - cos(t)) / t^2
            //
            // Now we must operate with values of [w]^2, but we do
            // this for accuracy's sake.
            //
            // So, from the previous equation we have:
            //
            // WW = [w]^2 = (.5 * (R + R^T) - I) / b
            //
            // WW = b[w]^2 = b[{x, y, z}]^2
            //
            //     [-y^2 - z^2     xy          xz   ]
            // b * [    xy     -x^2 - z^2      yz   ]
            //     [    xz         yz     -x^2 - y^2]
            //
            // If we consider WW{00} - WW{11} - WW{22}, we receive:
            //
            //   (y^2 + z^2) - (-x^2 - z^2) - (-x^2 - y^2) = 2x^2
            //
            // Thankfully, this expression surprisingly simplifies down
            // to subtracting cos(t) - thanks to the trace of the matrix!
            
            // We now calculate theta via arccos(cos(t)), but, we
            // opt for pi - arcsin(sin(t)) = arccos(cost(t)) as it
            // has more accuracy towards pi.
            const S sint = sqrt(sint_sq);
            const S t = S(M_PI) - std::asin(sint); // == acos(cos(t))
            const S binv = (t * t) / (S(1) - cost);

            const S x2 = (R[0] - cost) * binv;
            const S y2 = (R[4] - cost) * binv;
            const S z2 = (R[8] - cost) * binv;

            const S xy = ((R[1] + R[3]) / S(2)) * binv;
            const S xz = ((R[2] + R[6]) / S(2)) * binv;
            const S yz = ((R[5] + R[7]) / S(2)) * binv;

            // If we have either of x, y, or z, we can now
            // calculate {x, y, z} via division. We choose
            // the element with greatest magnitude so
            // our sqrt operation is accurate.
            //
            // Since the squaring operation removes the sign
            // of the operation, we can recover it by looking
            // at the original estimate of v.xyz = sin(t) * w
            if (x2 > y2)
            {
                if (x2 > z2)
                {
                    // x is largest
                    const S sign = (v(0) < 0 ? S(-1) : S(1));
                    const S x = sign * sqrt(x2);

                    v(0) = x;
                    v(1) = xy / x;
                    v(2) = xz / x;
                }
                else
                {
                    // z is largest
                    const S sign = (v(2) < 0 ? S(-1) : S(1));
                    const S z = sign * sqrt(z2);

                    v(0) = xz / z;
                    v(1) = yz / z;
                    v(2) = z;
                }
            }
            else
            {
                if (y2 > z2)
                {
                    // y is largest
                    const S sign = (v(1) < 0 ? S(-1) : S(1));
                    const S y = sign * sqrt(y2);

                    v(0) = xy / y;
                    v(1) = y;
                    v(2) = yz / y;
                }
                else
                {
                    // z is largest
                    const S sign = (v(2) < 0 ? S(-1) : S(1));
                    const S z = sign * sqrt(z2);

                    v(0) = xz / z;
                    v(1) = yz / z;
                    v(2) = z;
                }
            }
            
            return v;
        }

        static SO3 Exp(
            const Vec3<S>& v,
            const S eps = SO3ExpCoeffs<S>::THETA_SQ_EPS)
        {
            // t = v.norm()
            // w = v / t
            //
            // R = I + a[w]+ b[w]^2
            //
            // a = sin(t)
            // b = (1 - cos(t))
            const S thetaSq = v.SquaredNorm();
            const SO3ExpCoeffs<S> c(thetaSq, eps);

            SO3 res;
            ComputeExpMatrix3x3(&res.R[0], c.sint_t, c.omcost_t2, v(0), v(1), v(2));
            return res;
        }

        static SO3 Identity()
        {
            SO3 s;
            s.SetIdentity();
            return s;
        }

        inline S& operator[](int i) { return R[i]; }
        inline const S& operator[](int i) const { return R[i]; }

        inline S& operator()(int i, int j) { return R[i * 3 + j]; }
        inline const S& operator()(int i, int j) const { return R[i * 3 + j]; }

        // Row-order 3x3 rotation matrix
        S R[9];
    };

} // namespace sight