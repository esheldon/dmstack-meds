"""
TODO: 

    - include coadd image
    - seg maps

    - deal with fullStamp going off image, stamps going off edge
        - should just fill in part of the stamp
        - currently not writing anything
            - set weight==0 and some bitmask for missing pixels off edge

    - stamps too small, resize to "5-sigma"
    - make sure stamps are the same size in every band

    - set weight map zero for bad bits in bitmask
    - scale the images to common zero point
    - write the images to MEDS (easy)

    - make sure we agree on coordinate conventions


BUGS found in dmstack

    - trying to import test.py, which means we can't have a test.py in our cwd!
    - tries to connect to display, causing crash
    - installation instructions wrong: location of eups-setups.sh

"""

from __future__ import print_function
import numpy

import lsst.daf.persistence as dafPersist
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage


class LSSTProducer(object):
    """
    Class to help make MEDS files from LSST DM stack outputs.

    Usage:

     - Construct an instance of the class, corresponding to a single coadd patch.
     - Call getCatalog to get information about each object
     - Call getStamps on each source returned by getCatalog
       to get list of postage stamp Exposure objects that contain all the information
       we need about that postage stamp.
    """

    def __init__(self, butler, tract, patch, filter, limit=None):
        self.butler = butler
        self.meas = self.butler.get("deepCoadd_meas", tract=tract, patch=patch, filter=filter)
        self.coadds = self.butler.get("deepCoadd_calexp", tract=tract, patch=patch, filter=filter)
        self.ccds = self.coadds.getInfo().getCoaddInputs().ccds

        self.limit=limit

        self.loadImages()

    def getOverlappingEpochs(self, source):
        """
        determine epochs that overlap this source

        This is an iterator that yields (ccdRecord, Footprint)
        """
        for ccdRecord in self.ccds:
            if ccdRecord.contains(source.getCoord()):
                # get object footprint on the calexp
                fp = source.getFootprint().transform(
                    self.coadds.getWcs(),
                    ccdRecord.getWcs(),
                    ccdRecord.getBBox()
                )
                if fp.getArea() != 0:
                    yield ccdRecord, fp

    def getSourceBBox(self, source):
        """
        run through all overalapping footprints and get max size

        return the number of overlaps and max size
        """
        ncutout = 0
        maxWidth = 0
        for ccdRecord, footprint in self.getOverlappingEpochs(source):
            ncutout += 1
            maxWidth = max(maxWidth, footprint.getBBox().getWidth())
            maxWidth = max(maxWidth, footprint.getBBox().getHeight())
        return ncutout, maxWidth

    def _get_struct(self, n):
        dt = [
            ('id','i8'),
            ('box_size','i4'),
            ('ra','f8'),
            ('dec','f8'),
            ('ncutout','i4'),
        ]

        return numpy.zeros(n, dtype=dt)

    def getCatalog(self):
        if not hasattr(self,'catalog'):
            self.makeCatalog()

        return self.catalog

    def makeCatalog(self):
        """
        Make the catalog for all objects in the coadd patch.

        If `limit` was set in construction, only that many objects will be
        used.  The objects are selected from near the middle of the catalog to
        avoid just returning garbage on the edge.
        """

        limit=self.limit
        if limit is None:
            meas = self.meas
        else:
            start = len(self.meas)//2 - limit//2
            stop = start + limit
            meas = self.meas[start:stop]

        result = []
        for source in meas:
            ncutout, width = self.getSourceBBox(source)
            result.append((source.getId(), ncutout, width, source.getCoord()))

        n = len(result)
        data = self._get_struct(n)
        for i,objdata in enumerate(result):
            coord=objdata[3]

            data['id'][i]  = objdata[0]
            data['box_size'][i] = objdata[2]
            data['ra'][i] = coord.getRa().asDegrees()
            data['dec'][i] = coord.getDec().asDegrees()
            data['ncutout'][i]  = objdata[1]

        self.catalog=data

    def getDataId(self, ccdRecord):
        """Make a calexp data ID from a CCD ExposureRecord.

        Must be overridden for cameras whose data IDs have something other
        than ("visit", "ccd") as keys.
        """
        return dict(visit=ccdRecord["visit"], ccd=ccdRecord["ccd"])

    def loadImages(self):
        self.calexps = {}
        for ccdRecord in self.ccds:
            self.calexps[ccdRecord.getId()] = self.butler.get("calexp", self.getDataId(ccdRecord))

    def getStamps(self, obj_data):
        """
        """
        source = self.meas.find(obj_data['id'])  # find src record by ID
        stamps = []
        for ccdRecord, footprint in self.getOverlappingEpochs(source):
            calexp = self.calexps[ccdRecord.getId()]

            box_size = obj_data['box_size']

            extent = afwGeom.Extent2I(box_size, box_size)
            fullBBox = afwGeom.Box2I(
                footprint.getBBox().getMin(), 
                extent,
            )
            if calexp.getBBox().contains(fullBBox):
                fullStamp = calexp.Factory(calexp, fullBBox, afwImage.PARENT, True)
            else:
                fullStamp = None
            position = ccdRecord.getWcs().skyToPixel(source.getCoord())
            stamps.append((fullStamp, position))
        return stamps


def test_make_producer(limit=10):
    butler = dafPersist.Butler("/u/ki/boutigny/ki19/MACSJ2243/output/coadd_dir_cc/")

    tract = 0
    patch = "1,4"
    filter = "r"

    producer = LSSTProducer(
        butler,
        tract,
        patch,
        filter,
        limit=limit,
    )

    return producer

def test():
    """
    test making a producer
    """

    producer = test_make_producer(limit=10)

    cat = producer.getCatalog()


    stamps = producer.getStamps(cat[1])
    stamp, orig_pos = stamps[1]
    ncutout = len(stamps)
    
    numbers = numpy.arange(cat.size)

    # image
    arr=stamp.getMaskedImage().getImage().getArray()
    # variance
    var=stamp.getMaskedImage().getVariance().getArray()
    weight=1.0/var

    # bit mask
    mask=stamp.getMaskedImage().getMask().getArray()

    # no seg map yet

    # need to add to the object data; means writing object data last
    # instead of first
    wcs = stamp.getWcs().linearizePixelToSky(orig_pos, afwGeom.arcseconds)
    jacobian = wcs.getLinear().getMatrix()
    print("jacobian:",jacobian)

    # psf image
    psfobj=stamp.getPsf()
    psfim = psfobj.computeKernelImage(orig_pos).getArray()

    # unpack position in original image
    orig_row = orig_pos.getY()
    orig_col = orig_pos.getX()

    orig_start = stamp.getXY0()
    orig_start_row = orig_start.getY()
    orig_start_col = orig_start.getX()

    # location in the cutout
    cutout_cen = orig_pos - orig_start


    cutout_row = cutout_cen.getY()
    cutout_col = cutout_cen.getX()

    file_id=-9999

    return producer
 
