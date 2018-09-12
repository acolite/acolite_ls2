# -*- coding: utf-8 -*-
from __future__ import absolute_import, division

import numpy as np

# lsqfit_regression
# http://www.mbari.org/staff/etp3/regress/rules.htm

"""
For Model II regressions, neither X nor Y is an INDEPENDENT variable but both
are assumed to be DEPENDENT on some other parameter which is often unknown.
Neither are "controlled", both are measured, and both include some error. We do
not seek an equation of how Y varies in response to a change in X, but rather
we look for how they both co-vary in time or space in response to some other
variable or process. There are several possible Model II regressions. Which one
is used depends upon the specifics of the case. See Ricker (1973) or Sokal and
Rohlf (1995, pp. 541-549) for a discussion of which may apply. For convenience,
I have also compiled some rules of thumb.
"""


def lsqfity(X, Y):
    """
    Calculate a "MODEL-1" least squares fit.

    The line is fit by MINIMIZING the residuals in Y only.

    The equation of the line is:     Y = my * X + by.

    Equations are from Bevington & Robinson (1992)
    Data Reduction and Error Analysis for the Physical Sciences, 2nd Ed."
    pp: 104, 108-109, 199.

    Data are input and output as follows:

    my, by, ry, smy, sby = lsqfity(X,Y)
    X     =    x data (vector)
    Y     =    y data (vector)
    my    =    slope
    by    =    y-intercept
    ry    =    correlation coefficient
    smy   =    standard deviation of the slope
    sby   =    standard deviation of the y-intercept

    """

    X, Y = list(map(np.asanyarray, (X, Y)))

    # Determine the size of the vector.
    n = len(X)

    # Calculate the sums.

    Sx = np.sum(X)
    Sy = np.sum(Y)
    Sx2 = np.sum(X ** 2)
    Sxy = np.sum(X * Y)
    Sy2 = np.sum(Y ** 2)

    # Calculate re-used expressions.
    num = n * Sxy - Sx * Sy
    den = n * Sx2 - Sx ** 2

    # Calculate my, by, ry, s2, smy and sby.
    my = num / den
    by = (Sx2 * Sy - Sx * Sxy) / den
    ry = num / (np.sqrt(den) * np.sqrt(n * Sy2 - Sy ** 2))

    diff = Y - by - my * X

    s2 = np.sum(diff * diff) / (n - 2)
    smy = np.sqrt(n * s2 / den)
    sby = np.sqrt(Sx2 * s2 / den)

    return my, by, ry, smy, sby
