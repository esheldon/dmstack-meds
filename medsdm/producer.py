"""
TODO: 

    - seg maps

    - deal with fullStamp going off image, stamps going off edge
        - should just fill in part of the stamp
        - currently not writing anything
            - set weight==0 and some bitmask for missing pixels off edge

    - round stamps to 2^N or 3*2^N

    - provide ability to get sky variance only
    - write psf image

    - make sure we agree on coordinate conventions

    - copy flags_to_check into meta data
    - flags_to_check also for coadd
    - implement coadd only mode
    - implement subtracted nbrs for coadd
    - include coadd catalog  flags in object_data extension


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
import lsst.afw.table as afwTable
import lsst.afw.geom.ellipses as afwEllipses


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
        self.ref = self.butler.get("deepCoadd_ref", tract=tract, patch=patch,
                                   flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)
        self.coadd = self.butler.get("deepCoadd_calexp", tract=tract, patch=patch, filter=filter)
        self.ccds = self.coadd.getInfo().getCoaddInputs().ccds

        self.limit=limit

        self.loadImages()

    @staticmethod
    def projectBox(source, wcs, radius):
        pixel = afwGeom.Point2I(wcs.skyToPixel(source.getCoord()))
        box = afwGeom.Box2I()
        box.include(pixel)
        box.grow(radius)
        return box

    @staticmethod
    def getBoxRadiusFromWidth(width):
        return (width - 1)//2

    @staticmethod
    def getBoxWidthFromRadius(radius):
        return radius*2 + 1

    @staticmethod
    def computeBoxRadius(source):
        """
        Calculate the postage stamp "radius" for a source.

        TODO: make RADIUS_FACTOR and MIN_RADIUS configurable.
        """
        RADIUS_FACTOR = 5.0
        MIN_RADIUS = 16.0
        sigma = afwEllipses.Axes(source.getShape()).getA()
        if not (sigma >= MIN_RADIUS):  # handles small objects and NaNs
            return int(numpy.ceil(MIN_RADIUS))
        else:
            return int(numpy.ceil(RADIUS_FACTOR*sigma))

    def findOverlappingEpochs(self, source, radius=None, include_coadd=True):
        """Determine the epoch images that overlap a coadd source.

        Returns a list of tuples of `(box, ccd)`, where `box` is
        the postage-stamp bounding box in pixel coordinates and `ccd`
        is a an `ExposureRecord` containing CCD metadata.

        TODO: at present this skips images for which the postage stamp
        lies on the image boundary, because we don't yet have code to
        pad them.
        """
        result = []
        if include_coadd:
            sourceBox = self.projectBox(source, self.coadd.getWcs(), radius)
            imageBox = self.coadd.getBBox()
            if not imageBox.contains(sourceBox):
                raise NotImplementedError("Cannot process sources on the edge of the coadd")
            result.append((None, sourceBox))
        for ccd in self.ccds:
            sourceBox = self.projectBox(source, ccd.getWcs(), radius)
            imageBox = ccd.getBBox()
            if not imageBox.contains(sourceBox):
                continue
            result.append((ccd, sourceBox))
        return result

    def _get_struct(self, n):
        dt = [
            ('id','i8'),
            ('number','i8'),
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
            ref = self.ref
        else:
            start = len(self.ref)//2 - limit//2
            stop = start + limit
            ref = self.ref[start:stop]

        result = []
        nChildKey = ref.schema.find("deblend_nChild").key
        for source in ref:
            if source.get(nChildKey) != 0:
                # Skip parent objects, since we'll also process their children.
                continue
            radius = self.computeBoxRadius(source)
            epochs = self.findOverlappingEpochs(source, radius=radius)
            result.append((source.getId(),
                           len(epochs),
                           self.getBoxWidthFromRadius(radius),
                           source.getCoord()))

        n = len(result)
        data = self._get_struct(n)
        for i, (objId, nEpochs, width, coord) in enumerate(result):
            data['id'][i]  = objId
            data['number'][i]  = objId
            data['box_size'][i] = width
            data['ra'][i] = coord.getRa().asDegrees()
            data['dec'][i] = coord.getDec().asDegrees()
            data['ncutout'][i]  = nEpochs

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
        TODO

        Currently calexp.getBBox().contains(fullBBox) is checked which returns
        False if the full stamp is not contained
        """
        source = self.ref.find(obj_data['id'])  # find src record by ID
        stamps = []
        width = obj_data['box_size']
        radius = self.getBoxRadiusFromWidth(width)
        coaddFluxMag0 = self.coadd.getCalib().getFluxMag0()[0]
        for ccdRecord, bbox in self.findOverlappingEpochs(source, radius=radius):
            if ccdRecord is None:
                # this is a coadd stamp
                fullStamp = self.coadd.Factory(self.coadd, bbox=bbox, origin=afwImage.PARENT,
                                               deep=True)
                position = self.coadd.getWcs().skyToPixel(source.getCoord())
            else:
                calexp = self.calexps[ccdRecord.getId()]
                calexpFluxMag0 = calexp.getCalib().getFluxMag0()[0]
                fluxScaling = coaddFluxMag0/calexpFluxMag0
                assert bbox.getWidth() == width and bbox.getHeight() == width
                assert calexp.getBBox().contains(bbox)
                fullStamp = calexp.Factory(calexp, bbox=bbox, origin=afwImage.PARENT, deep=True)
                mi = fullStamp.getMaskedImage()
                mi *= fluxScaling
                position = ccdRecord.getWcs().skyToPixel(source.getCoord())
            stamps.append((fullStamp, position))
        return stamps


def test_make_producer(limit=10):
    butler = dafPersist.Butler("/datasets/hsc/repo/rerun/private/hchiang2/RC/DM-10129")

    tract = 8766
    patch = "4,4"
    filter = "HSC-I"

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


    index=1
    stamps = producer.getStamps(cat[index])
    stamp, orig_pos = stamps[1]


    flag_dict = stamp.mask.getMaskPlaneDict()
    flags_to_ignore = ['DETECTED','DETECTED_NEGATIVE']
    flags_to_check = 0
    for key in flag_dict:
        if key not in flags_to_ignore:
            flags_to_check |= stamp.mask.getPlaneBitMask(key)

    print("flags to check:",flags_to_check)
    ncutout = len(stamps)
    
    numbers = numpy.arange(cat.size)

    # image
    arr=stamp.image.array
    # variance
    var=stamp.variance.array
    weight=1.0/var

    # bit mask
    mask=stamp.mask.array

    # no seg map yet

    # need to add to the object data; means writing object data last
    # instead of first
    wcs = stamp.getWcs().linearizePixelToSky(orig_pos, afwGeom.arcseconds)
    jacobian = wcs.getLinear().getMatrix()
    print("jacobian:",jacobian)

    # psf image of coadd, with a little padding this should be
    # bigger than any SE psf
    coadd_psfobj = producer.coadd.getPsf()
    coadd_pos = producer.ref.find(cat['id'][index]).getCentroid()
    coadd_psfim = coadd_psfobj.computeKernelImage(coadd_pos).array

    # psf image
    psfobj=stamp.getPsf()
    psfim = psfobj.computeKernelImage(orig_pos).array

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
 
