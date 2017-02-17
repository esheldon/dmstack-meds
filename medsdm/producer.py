import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage


class LSSTProducer(object):

    def __init__(self, butler, tract, patch, filter):
        self.butler = butler
        self.meas = self.butler.get("deepCoadd_meas", tract=tract, patch=patch, filter=filter)
        self.coadd = self.butler.get("deepCoadd_calexp", tract=tract, patch=patch, filter=filter)
        self.ccds = self.coadds.getInfo().getCoaddInputs().ccds

    def getOverlappingEpochs(self, source):
        for ccdRecord in self.ccds:
            if ccdRecord.contains(source.getCoord()):
                # get object footprint on the calexp
                fp = source.getFootprint().transform(
                    self.coadd.getWcs(),
                    ccdRecord.getWcs(),
                    ccdRecord.getBBox()
                )
                if fp.getArea() != 0:
                    yield ccdRecord, fp

    def getSourceBBox(self, source):
        count = 0
        maxWidth = 0
        maxHeight = 0
        for ccdRecord, footprint in self.getOverlappingEpochs(source):
            count += 1
            maxWidth = max(maxWidth, fp.getBBox().getWidth())
            maxHeight = max(maxHeight, fp.getBBox().getHeight())
        return count, maxWidth, maxHeight

    def makeCatalog(self):
        result = []
        for src in self.meas:
            count, width, height = self.getSourceBBox(source)
            result.append((source.getId(), count, width, height))
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
            self.calexp[ccdRecord.getId()] = self.butler.get("calexp", self.makeDataId(ccdRecord))

    def makeStamps(self, sourceId, count, width, height):
        source = self.meas.find(sourceId)  # find src record by ID
        stamps = []
        for ccdRecord, footprint in self.getOverlappingEpochs(source):
            calexp = self.calexps(ccdRecord.getId())
            fullBBox = afwGeom.Box2I(footprint.getBBox().getMin(), width, height)
            if calexp.getBBox().contains(fullBBox):
                fullStamp = calexp.Factory(calexp, fullBBox, afwImage.PARENT, True)
            else:
                fullStamp = None
            stamps.append(fullStamp)
        return stamps
