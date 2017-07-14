"""
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

from .defaults import DEFAULT_PRODUCER_CONFIG 

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

    def __init__(self, butler, tract, patch, filter, limit=None, config=None):

        self.setConfig(config)
        self.butler = butler
        self.ref = self.butler.get("deepCoadd_ref", tract=tract, patch=patch,
                                   flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)
        self.coadd = self.butler.get("deepCoadd_calexp", tract=tract, patch=patch, filter=filter)
        self.ccds = self.coadd.getInfo().getCoaddInputs().ccds
        self.coaddSegMap = None

        self.limit=limit

        self.loadImages()

    def setConfig(self, input_config):

        config={}
        config.update(DEFAULT_PRODUCER_CONFIG)
        if input_config is not None:
            config.update(input_config)

        self.config=config

    def makeCoaddSegMap(self, radius=5):
        """Build a *really* naive segmentation map by inserting small circular regions
        for each source, going from faintest to brightest in blends.
        """
        result = afwImage.ImageI(self.coadd.getBBox())
        obj_data = self.getCatalog()
        id_to_number = {obj["id"]: obj["number"] for obj in obj_data}
        for parent in self.ref.getChildren(0):
            children = self.ref.getChildren(parent.getId())
            if len(children) == 0:
                data = [(parent.getPsfFlux(), parent.getCentroid(), parent.getId())]
            else:
                data = [(child.getPsfFlux(), child.getCentroid(), child.getId())
                        for child in children]
                data.sort()
            for flux, centroid, objId in data:
                number = id_to_number.get(objId, None)
                if number is None:
                    # This should only happen when we're limited to a small number of objects
                    # for testing
                    continue
                stencil = afwGeom.SpanSet.fromShape(radius, afwGeom.Stencil.CIRCLE,
                                                    afwGeom.Point2I(centroid))
                stencil.setImage(result, id_to_number[objId], region=result.getBBox(), doClip=True)
        return result

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

    def computeBoxRadius(self, source):
        """
        Calculate the postage stamp "radius" for a source.

        TODO: make RADIUS_FACTOR and MIN_RADIUS configurable.
        """
        conf = self.config

        min_radius = conf['min_box_size']/2
        max_radius = conf['max_box_size']/2

        sigma = afwEllipses.Axes(source.getShape()).getA()

        if numpy.isnan(sigma):
            sigma = 1.0

        rad = conf['radius_factor']*sigma
        if rad < min_radius:
            rad = min_radius
        elif rad > max_radius:
            rad = max_radius

        return int(numpy.ceil(rad))

        #if not (sigma >= min_radius):  # handles small objects and NaNs
        #    return int(numpy.ceil(min_radius))
        #else:
        #    return int(numpy.ceil(conf['radius_factor']*sigma))

    def findOverlappingEpochs(self, source, radius=None):
        """Determine the epoch images that overlap a coadd source.

        Returns a list of tuples of `(box, ccd)`, where `box` is
        the postage-stamp bounding box in pixel coordinates and `ccd`
        is a an `ExposureRecord` containing CCD metadata.

        TODO: at present this skips images for which the postage stamp
        lies on the image boundary, because we don't yet have code to
        pad them.
        """
        result = []
        if self.config['include_coadd']:
            sourceBox = self.projectBox(source, self.coadd.getWcs(), radius)
            imageBox = self.coadd.getBBox()
            result.append((None, sourceBox))
        for ccd in self.ccds:
            sourceBox = self.projectBox(source, ccd.getWcs(), radius)
            imageBox = ccd.getBBox()
            if not imageBox.overlaps(sourceBox):
                continue
            result.append((ccd, sourceBox))
        return result

    def _get_struct(self, n):
        dt = [
            ('id','i8'),
            ('number','i4'),
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
            data['number'][i]  = i + 1
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

    @staticmethod
    def getPaddedSubImage(original, bbox):
        region = original.getBBox()
        if region.contains(bbox):
            return original.Factory(original, bbox, afwImage.PARENT, True)
        result = original.Factory(bbox)
        bbox2 = afwGeom.Box2I(bbox)
        bbox2.clip(region)
        if isinstance(original, afwImage.Exposure):
            result.setPsf(original.getPsf())
            result.setWcs(original.getWcs())
            result.setCalib(original.getCalib())
            result.image.array[:, :] = float("nan")
            result.variance.array[:, :] = float("inf")
            result.mask.array[:, :] = numpy.uint16(result.mask.getPlaneBitMask("NO_DATA"))
            subIn = afwImage.MaskedImageF(original.maskedImage, bbox=bbox2,
                                          origin=afwImage.PARENT, deep=False)
            result.maskedImage.assign(subIn, bbox=bbox2, origin=afwImage.PARENT)
        elif isinstance(original, afwImage.ImageI):
            result.array[:, :] = 0
            subIn = afwImage.ImageI(original, bbox=bbox2,
                                   origin=afwImage.PARENT, deep=False)
            result.assign(subIn, bbox=bbox2, origin=afwImage.PARENT)
        else:
            raise ValueError("Image type not supported")
        return result

    def getStamps(self, obj_data):
        """
        TODO

        Currently calexp.getBBox().contains(fullBBox) is checked which returns
        False if the full stamp is not contained
        """

        conf=self.config
        if self.coaddSegMap is None:
            self.coaddSegMap = self.makeCoaddSegMap(radius=conf['fake_seg_radius'])

        source = self.ref.find(obj_data['id'])  # find src record by ID
        stamps = []
        width = obj_data['box_size']
        radius = self.getBoxRadiusFromWidth(width)
        coaddFluxMag0 = self.coadd.getCalib().getFluxMag0()[0]
        for ccdRecord, bbox in self.findOverlappingEpochs(source, radius=radius):
            if ccdRecord is None:
                # this is a coadd stamp
                fullStamp = self.getPaddedSubImage(self.coadd, bbox=bbox)
                position = self.coadd.getWcs().skyToPixel(source.getCoord())
                segmap = self.getPaddedSubImage(self.coaddSegMap, bbox=bbox)
            else:
                calexp = self.calexps[ccdRecord.getId()]
                calexpFluxMag0 = calexp.getCalib().getFluxMag0()[0]
                fluxScaling = coaddFluxMag0/calexpFluxMag0
                assert bbox.getWidth() == width and bbox.getHeight() == width
                fullStamp = self.getPaddedSubImage(calexp, bbox=bbox)

                # scales both the image and variance image
                fullStamp.maskedImage *= fluxScaling

                position = ccdRecord.getWcs().skyToPixel(source.getCoord())
                segmap = None
            stamps.append((fullStamp, position, segmap))
        return stamps


def test_make_producer(limit=10, config=None):
    #butler = dafPersist.Butler("/datasets/hsc/repo/rerun/private/hchiang2/RC/DM-10129")
    butler = dafPersist.Butler("/project/hsc_rc/w_2017_26/DM-11165")

    tract = 8766
    patch = "4,4"
    filter = "HSC-I"

    producer = LSSTProducer(
        butler,
        tract,
        patch,
        filter,
        limit=limit,
        config=config,
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
    stamp, orig_pos, seg_map = stamps[1]


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
 
