import numpy as np


SECONDS_PER_DAY = 60. * 60. * 24.
SECONDS_PER_MICROSECOND = .000001


def convert_times_to_seconds_since_start( times ):
    """
    Given an array of datetime objects, return a NumPy array of seconds since
    the earliest time. This assumes that the given times are in chronological
    order.
    """
    seconds = np.empty( times.shape, np.double )

    for index in xrange( times.shape[ 0 ] ):
        delta = times[ index ] - times[ 0 ]
        seconds[ index ] = \
            delta.days * SECONDS_PER_DAY + \
            delta.seconds + \
            delta.microseconds * SECONDS_PER_MICROSECOND

    return seconds
