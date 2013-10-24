#!/usr/bin/env python

"""
tests for reading moss files

designed for pytest
"""

import os

from post_gnome.readers import read_moss

sample_moss_file_root = "sample_data/tests/tests"
sample_moss_file_root2 = "sample_data/rings/rings"

def test_ReadSpillInfo():
    ## filename without the extensions.
    filename = sample_moss_file_root

    spill_info = read_moss.ReadSpillInfo(filename)

    assert spill_info['ISSUED'] == '13:20, 7/24/09'
    assert spill_info['VALIDFOR'] == '06:00, 7/17/09'

    # there really should be more...

def test_ReadMossLeFiles():
    ## fixme -- lot's more to test!
    ## filename without the extensions.
    filename = sample_moss_file_root

    blackLEs, redLEs = read_moss.ReadMossLeFiles(filename)

    assert len(blackLEs) == 1000
    assert len(redLEs) == 1000

    assert blackLEs[0].latitudeStr == '24.85016'
    assert blackLEs[0].longitudeStr == '-96.84725'

    assert redLEs[0].latitudeStr == '24.85206'
    assert redLEs[0].longitudeStr == '-96.86061'

def test_ReadMossPolygons():
    ## fixme -- lot's more to test!
    ## filename without the extensions.
    filename = sample_moss_file_root2

    (landPolygons,
     heavyPolygons,
     mediumPolygons,
     lightPolygons,
     uncertaintyPolygons) = read_moss.ReadMossPolygons(filename, printDiagnostic=False)

    # print  (landPolygons,
    #  heavyPolygons,
    #  mediumPolygons,
    #  lightPolygons,
    #  uncertaintyPolygons)


