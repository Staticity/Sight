#pragma once

#include <lie/so3.hpp>
#include <linear/vec.hpp>

namespace sight
{
    template <typename S>
    class SE3
    {
    public:
        SE3()
            : R()
            , t()
        {}

        SE3(const SO3<S>& R)
            : R(R)
            , t(S())
        {}

        SE3(const SO3<S>& R, const Vec3<S>& t)
            : R(R)
            , t(t)
        {}

        inline S& operator()(int i, int j)
        {
            return j < 3 ? R(i, j) : t(i);
        }

        inline const S& operator()(int i, int j) const
        {
            return j < 3 ? R(i, j) : t(i);
        }

        SE3 operator*(const SE3& s) const
        {
            SE3 res;
            res.R = R * s.R;
            res.t = R * s.t + t;
            return res;
        }

        SE3 operator*(const SO3<S>& s) const
        {
            SE3 res;
            res.R = R * s.R;
            res.t = t;
            return res;
        }

        Vec3<S> operator*(const Vec3<S>& v) const
        {
            return R * v + t;
        }

        Vec3<S> Position() const
        {
            return -(R.Inverse() * t);
        }

        SE3 Inverse() const
        {
            SE3 inv;
            inv.R = R.Inverse();
            inv.t = -(inv.R * t);
            return inv;
        }

        void SetIdentity()
        {
            R.SetIdentity();
            t.Fill(S());
        }

        static SE3 Identity()
        {
            SE3 T;
            T.SetIdentity();
            return T;
        }

        Vec<S, 6> Log() const
        {
            Vec<S, 6> log;
            
            // Re-use SO(3)'s log operation
            const Vec3<S> logR = R.Log();
            log[0] = logR[0];
            log[1] = logR[1];
            log[2] = logR[2];

            // We know that our translational component is of the form:
            //
            //     t = Vu
            //
            // And we'd like to solve for u. Thus, if we can find V^-1,
            // we can get u like so:
            //
            //    u = V.inv() * t
            //
            // If V is of the form:
            //
            //    V = I + a[w] + b[w]^2
            //
            //    a = sin(t) / t
            //    b = (1 - cos(t)) / t^2
            //    c = (t - sin(t)) / t^3 (unused above)
            //
            // We can find its inverse in a similar form:
            //
            //    V.inv() = I + d[w] + e[w]^2
            //
            //    d = -.5
            //    e = (b/2 - c) / a
            //    e = (b - a/2) / (b * t^2)
            //
            // We can choose between either expression for e,
            // whichever provides better accuracy for values of
            // theta.
            //
            // TODO: Prove the equation of V^-1

            const S t2 = logR.SquaredNorm();
            const SE3ExpCoeffs<S> coeffs(t2);

            // Match the comment names
            const S a = coeffs.sint_t;
            const S b = coeffs.omcost_t2;
            const S c = coeffs.tmsint_t3;
            
            const S d = -.5;
            S e;

            // Values of theta close to 0.
            if (t2 < SE3ExpCoeffs<S>::THETA_SQ_EPS)
            {
                // Taylor series approximation for e at x = 0
                e = S(1.0 / 12) + t2 * (S(1.0 / 720) + t2 * S(1.0 / 30240));
            }
            else if (t2 > S(9))
            {
                // Values of theta somewhat around pi (not crucial)
                //
                // We use the second expression for e.
                e = (b - a/2) / (b * t2);
            }
            else
            {
                // Regular values of theta.
                //
                // We use the first expression for e
                e = (b / S(2) - c) / a;
            }

            S Vinv[3 * 3];
            ComputeExpMatrix3x3(&Vinv[0], d, e, logR[0], logR[1], logR[2]);
            
            // u = V.inv() * t
            log[3] = Vinv[0] * t(0) + Vinv[1] * t(1) + Vinv[2] * t(2);
            log[4] = Vinv[3] * t(0) + Vinv[4] * t(1) + Vinv[5] * t(2);
            log[5] = Vinv[6] * t(0) + Vinv[7] * t(1) + Vinv[8] * t(2);

            return log;
        }

        // `v` is the set of coefficients for the
        // generators of SE(3). Stored as follows:
        //
        // rx, ry, rz, ux, uy, uz
        static SE3 Exp(
            const Vec<S, 6>& v,
            const S eps = SO3ExpCoeffs<S>::THETA_SQ_EPS)
        {
            const S rx = v(0);
            const S ry = v(1);
            const S rz = v(2);
            const S ux = v(3);
            const S uy = v(4);
            const S uz = v(5);

            const S thetaSq = rx * rx + ry * ry + rz * rz;
            const SE3ExpCoeffs<S> c(thetaSq, eps);
            
            SE3 T;
            // R = I + a[v] + b[v]^2
            // a = sin(t) / t
            // b = (1 - cos(t)) / t^2
            ComputeExpMatrix3x3(&T.R[0], c.sint_t, c.omcost_t2, rx, ry, rz);

            // t = Vu
            //
            //
            // where V = I + a[v] + b[v]^2
            //       a = (1 - cos(t)) / t
            //       b = (t - sin(t)) / t^3
            S V[3 * 3];
            const S theta = sqrt(thetaSq);
            ComputeExpMatrix3x3(&V[0], c.omcost_t2, c.tmsint_t3, rx, ry, rz);

            // t = Vu
            T.t(0) = V[0] * ux + V[1] * uy + V[2] * uz;
            T.t(1) = V[3] * ux + V[4] * uy + V[5] * uz;
            T.t(2) = V[6] * ux + V[7] * uy + V[8] * uz;

            return T;
        }
        
        SO3<S> R;
        Vec3<S> t;
    };
}