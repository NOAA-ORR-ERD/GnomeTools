#!/usr/bin/env python

"""
test of the BuildStartTimes code

not much here, but enough to test the Gaps syntax...

"""

import nose

import sys

import datetime 

sys.path.append("../Scripts")

import BuildStartTimes

def test_simple_gaps():
    run_time = datetime.timedelta(5) # days
    gaps = (( datetime.datetime(1999, 1, 3), datetime.datetime(1999, 2, 2)  ),
            ( datetime.datetime(2001, 4, 5, 12), datetime.datetime(2001, 4, 20) ),
            )
    gapset = BuildStartTimes.SimpleGapSet(gaps)

    # these should be in the gaps
    assert gapset.TimeInGap(datetime.datetime(1999, 1, 3), run_time)
    assert gapset.TimeInGap(datetime.datetime(1999, 2, 1), run_time)
    assert gapset.TimeInGap(datetime.datetime(2001, 4, 10), run_time)
      # on the edge
    assert gapset.TimeInGap(datetime.datetime(2001, 4, 20), run_time)
    assert gapset.TimeInGap(datetime.datetime(2001, 4, 5, 12), run_time)
      # outside, but within run_time
    assert gapset.TimeInGap(datetime.datetime(2001, 4, 1, 12), run_time)

    # these should not be in the gaps
    assert not gapset.TimeInGap(datetime.datetime(2001, 4, 21, 12), run_time)
    assert not gapset.TimeInGap(datetime.datetime(2000, 2, 1), run_time)
    assert not gapset.TimeInGap(datetime.datetime(1998, 12, 28, 23), run_time)


def test_empty_gaps():
    run_time = datetime.timedelta(hours=72) # hrs 
    gapset = BuildStartTimes.EmptyGapSet()

    # a few arbitrary times -- they should all return False
    assert not gapset.TimeInGap(datetime.datetime(1999, 1, 3), run_time)
    assert not gapset.TimeInGap(datetime.datetime(2099, 1, 3), run_time)
    assert not gapset.TimeInGap(datetime.datetime(1099, 12, 31), run_time)
    
    