#pragma once

#include <limits>

namespace sight
{
    // Contains relevant constants related
    // to the exponential mapping of SO3.
    template <typename S>
    struct SO3ExpCoeffs
    {
        static constexpr S THETA_SQ_EPS = std::numeric_limits<S>::epsilon() * 25;

        SO3ExpCoeffs(
            const S thetaSq,
            const S eps = THETA_SQ_EPS)
        {
            const S t2 = thetaSq;
            
            // For small angles
            if (thetaSq <= eps)
            {
                // Use Taylor Series expansion at theta=0. Note that we nest
                // multiplication by -tt, which will alternate the signs.
                sint_t = SinThetaOverTheta_Taylor(thetaSq);
                omcost_t2 = OneMinusCosThetaOverThetaSq_Taylor(thetaSq);
            }
            else
            {
                // Should we waste this sqrt operation - not save it?
                const S t = sqrt(thetaSq);

                // For normal angles
                sint_t = sin(t) / t;
                omcost_t2 = (1 - cos(t)) / thetaSq;
            }
        }

        // sin(t) / t
        S sint_t;

        // (1 - cos(t)) / t^2
        S omcost_t2;

        // (1 - sin(t) / t) / t^2
        S omsint_t_t2;
    };

    template <typename S>
    struct SE3ExpCoeffs
    {
        static constexpr S THETA_SQ_EPS = std::numeric_limits<S>::epsilon() * 25;

        SE3ExpCoeffs(
            const S thetaSq,
            const S eps = THETA_SQ_EPS)
        {
            const S t2 = thetaSq;
            
            // For small angles
            if (thetaSq <= eps)
            {
                // Use Taylor Series expansion at theta=0. Note that we nest
                // multiplication by -tt, which will alternate the signs.
                sint_t = SinThetaOverTheta_Taylor(thetaSq);
                omcost_t2 = OneMinusCosThetaOverThetaSq_Taylor(thetaSq);
                tmsint_t3 = ThetaMinusSinThetaOverThetaCubed_Taylor(thetaSq);
            }
            else
            {
                // Should we waste this sqrt operation - not save it?
                const S t = sqrt(thetaSq);

                // For normal angles
                const S sint = sin(t);
                sint_t = sint / t;
                omcost_t2 = (1 - cos(t)) / thetaSq;
                tmsint_t3 = (t - sint) / (t * thetaSq);
            }
        }

        // sin(t) / t
        S sint_t;

        // (1 - cos(t)) / t
        S omcost_t2;

        // (t - sin(t)) / t^3
        S tmsint_t3;

    };

    // Computes sin(t) / t
    template <typename S>
    inline S SinThetaOverTheta_Taylor(const S thetaSq)
    {
        const S t2 = thetaSq;
        return S(1) - t2 * (S(1.0 / 6) - t2 * (S(1.0 / 120) - t2 * S(1.0 / 5040)));
    }

    // Computes (1 - cos(t)) / t^2
    template <typename S>
    inline S OneMinusCosThetaOverThetaSq_Taylor(const S thetaSq)
    {
        const S t2 = thetaSq;
        return S(.5) - t2 * (S(1.0 / 24) - t2 * (S(1.0 / 720) - t2 * S(1.0 / 40320)));
    }

    // Computes (t - sin(t)/t) / t^3
    template <typename S>
    inline S ThetaMinusSinThetaOverThetaCubed_Taylor(const S thetaSq)
    {
        const S t2 = thetaSq;
        return S(1.0 / 6) - t2 * (S(1.0 / 120) - t2 * (S(1.0 / 5040) - t2 * S(1.0 / 362880)));
    }

    // Colmputes t / sin(t)
    template <typename S>
    inline S ThetaOverSinTheta_Taylor(const S thetaSq)
    {
        const S t2 = thetaSq;
        return S(1) + t2 * (S(1.0 / 6) + t2 * (S(7.0 / 360) + t2 * S(31.0 / 15120)));
    }

    // Computes arcsin(x)
    template <typename S>
    inline S Arcsin_Taylor(const S xSq)
    {
        return S(1) + xSq * (S(1.0 / 6) + xSq * (S(3.0 / 40) + xSq * S(5.0 / 112)));
    }

#if 0
    template <typename S>
    inline S SinThetaOverTheta(const S theta)
    {
        return sin(theta) / theta;
    }

    template <typename S>
    inline S OneMinusCosThetaOverThetaSq(const S theta)
    {
        return (1 - cos(theta)) / (theta * theta);
    }

    template <typename S>
    inline S OneMinusCosThetaOverThetaSq(const S theta, const S thetaSq)
    {
        return (1 - cos(theta)) / (thetaSq);
    }

    inline S ThetaMinusSinThetaOverThetaCubed(const S theta)
    {
        const S thetaCubed = theta * theta * theta;
        return (theta - sin(theta)) / thetaCubed;
    }

    inline S ThetaMinusSinThetaOverThetaCubed(const S theta, const S thetaSq)
    {
        const S thetaCubed = theta * thetaSq;
        return (theta - sin(theta)) / thetaCubed;
    }
#endif

    // Computes the matrix:
    //
    //    M = 1 + a[w] + b[w]^2
    //
    // where [w] is a skew symmetric matrix.
    //
    // and w = {x, y, z}
    template <typename S>
    void ComputeExpMatrix3x3(S M[3 * 3], S a, S b, S x, S y, S z)
    {
        const S ax = a * x;
        const S ay = a * y;
        const S az = a * z;

        const S bx = b * x;
        const S by = b * y;
        const S bz = b * z;

        const S bx2 = bx * x;
        const S by2 = by * y;
        const S bz2 = bz * z;

        const S bxy = bx * y;
        const S bxz = bx * z;
        const S byz = by * z;

        // K = [v] =
        //
        // [ 0 -z  y]
        // [ z  0 -x]
        // [-y  x  0]
        //
        // K^2 =
        //
        // [-y^2 - z^2      xy           xz    ]
        // [    xy      -x^2 - z^2       yz    ]
        // [    xz          yz      -x^2 - y^2 ]
        M[0] = 1 - by2 - bz2;
        M[4] = 1 - bx2 - bz2;
        M[8] = 1 - bx2 - by2;

        M[1] = bxy - az;
        M[2] = bxz + ay;

        M[3] = bxy + az;
        M[5] = byz - ax;

        M[6] = bxz - ay;
        M[7] = byz + ax;
    }

}