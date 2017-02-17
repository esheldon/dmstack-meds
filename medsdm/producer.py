import lsst.daf.persistence as dafPersist
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage


class LSSTProducer(object):
    """Class to help make MEDS files from LSST DM stack outputs.

    Usage:

     - Construct an instance of the class, corresponding to a single coadd patch.
     - Call makeCatalog to get a list of tuples of (source_id, num_epochs, width).
     - Call loadImages to load all of the CCD-level images into memory.
     - Call makeStamps on each (source_id, num_epochs, width) returned by makeCatalog
       to get list of postage stamp Exposure objects that contain all the information
       we need about that postage stamp.
    """

    def __init__(self, butler, tract, patch, filter):
        self.butler = butler
        self.meas = self.butler.get("deepCoadd_meas", tract=tract, patch=patch, filter=filter)
        self.coadds = self.butler.get("deepCoadd_calexp", tract=tract, patch=patch, filter=filter)
        self.ccds = self.coadds.getInfo().getCoaddInputs().ccds

    def getOverlappingEpochs(self, source):
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
        count = 0
        maxWidth = 0
        for ccdRecord, footprint in self.getOverlappingEpochs(source):
            count += 1
            maxWidth = max(maxWidth, footprint.getBBox().getWidth())
            maxWidth = max(maxWidth, footprint.getBBox().getHeight())
        return count, maxWidth

    def makeCatalog(self, limit=None):
        """Return a list of tuples of (id, num_epochs, stamp_width) for all
        objects in the coadd patch.

        If `limit` is not None, only tuples for that many objects will be
        returned.  The objects are selected from near the middle of the
        catalog to avoid just returning garbage on the edge.
        """
        if limit is None:
            meas = self.meas
        else:
            start = len(self.meas)//2 - limit//2
            stop = start + limit
            meas = self.meas[start:stop]
        result = []
        for source in meas:
            count, width = self.getSourceBBox(source)
            result.append((source.getId(), count, width, source.getCoord()))

        return result

    def makeDataId(self, ccdRecord):
        """Make a calexp data ID from a CCD ExposureRecord.

        Must be overridden for cameras whose data IDs have something other
        than ("visit", "ccd") as keys.
        """
        return dict(visit=ccdRecord["visit"], ccd=ccdRecord["ccd"])

    def loadImages(self):
        self.calexps = {}
        for ccdRecord in self.ccds:
            self.calexps[ccdRecord.getId()] = self.butler.get("calexp", self.makeDataId(ccdRecord))

    def makeStamps(self, sourceId, count, width):
        """
        TODO: remove count
              deal with fullStamp going off image
        """
        source = self.meas.find(sourceId)  # find src record by ID
        stamps = []
        for ccdRecord, footprint in self.getOverlappingEpochs(source):
            calexp = self.calexps[ccdRecord.getId()]

            extent = afwGeom.Extent2I(width, width)
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


def test_make_producer():
    butler = dafPersist.Butler("/u/ki/boutigny/ki19/MACSJ2243/output/coadd_dir_cc/")

    tract = 0
    patch = "1,4"
    filter = "r"

    producer = LSSTProducer(
        butler,
        tract,
        patch,
        filter,
    )

    return producer

def test():
    """
    test making a producer
    """

    producer = test_make_producer()

    cat = producer.makeCatalog(limit=10)
    producer.loadImages()

    objdata = cat[1]
    skycoord = objdata[3]
    ra = skycoord.getRa().asDegrees()
    dec = skycoord.getDec().asDegrees()

    stamps = producer.makeStamps(objdata[0], objdata[1], objdata[2])
    stamp, pos = stamps[1]

    arr=stamp.getMaskedImage().getImage().getArray()
    var=stamp.getMaskedImage().getVariance().getArray()
    mask=stamp.getMaskedImage().getMask().getArray()

    # this will be in arcsec
    wcs = stamp.getWcs().linearizePixelToSky(pos, afwGeom.arcseconds)
    jacobian = wcs.getLinear().getMatrix()

    psfobj=stamp.getPsf()

    psfim = psfobj.computeKernelImage(pos).getArray()

    cutout_cen = pos - stamp.getXY0()

    row = cutout_cen.getY()
    col = cutout_cen.getX()

    return producer
 
