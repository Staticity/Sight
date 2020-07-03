import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

def OneMinusCosOverThetaSq(theta, nTerms=3):
    p = 1.0
    d = 2.0

    c = 2

    res = 0.0
    for i in range(nTerms):
        term = p / d if i % 2 == 0 else -p / d
        res += term
        p *= theta * theta
        d *= (c + 1) * (c + 2)
        c += 2
    # print()
    return res

def SinThetaOverTheta(theta, nTerms=3):
    p = 1.0
    d = 1.0
    c = 1
    res = 0.0

    for i in range(nTerms):
        sign = 1 if i % 2 == 0 else -1
        term = sign * p / d
        res += term
        p *= theta * theta
        d *= (c + 1) * (c + 2)
        c += 2
    return res

if __name__ == '__main__':
    # find the theta which minimizes the sum residuals over
    # all of the inputs.
    hi_value = .1
    test_values = np.linspace(0.0, hi_value, 100)
    test_thresh = np.linspace(0.0, .5, 1000)

    errors = []
    for thresh in test_thresh:
        errSq = 0.0
        sin_err = []
        cos_err = []
        for value in test_values:
            exp_sin = SinThetaOverTheta(value, 100)
            exp_cos = OneMinusCosOverThetaSq(value, 100)
            if value < thresh:
                test_sin = SinThetaOverTheta(value, 4)
                test_cos = OneMinusCosOverThetaSq(value, 4)
            else:
                test_sin = math.sin(value) / value
                test_cos = (1 - math.cos(value)) / (value ** 2)

            sin_err.append(test_sin - exp_sin)
            cos_err.append(test_cos - exp_cos)
            errSq += (test_sin - exp_sin) ** 2
            errSq += (test_cos - exp_cos) ** 2
        
        # fig, ax = plt.subplots(1, 1)
        # ax.set_title(f'Thresh: {thresh}')
        # ax.plot(test_values, sin_err, c='r')
        # ax.plot(test_values, cos_err, c='g')
        # plt.show()
        errors.append(errSq)
    
    plt.plot(test_thresh, errors)
    plt.show()


    # theta = .54321
    # print(SinThetaOverTheta(theta, 100))
    # print(OneMinusCosOverThetaSq(theta, 100))
