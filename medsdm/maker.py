from __future__ import print_function
import numpy
import meds
import fitsio

class DMMedsMaker(meds.MEDSMaker):
    """

    - we don't have access to files
    - so the image_info will currently be blank
    - Need ncutout, box size for every object to reserve the image on disk
        - We might not know the stamp sizes, so maybe for now fix them to 48x48
        - won't know number of epochs ahead of time
            - two pass for now?
    - need scale to put all on same zero point, if not already done

    - will probably have cutout_row
    - make sure to get orig_start_row
    """

    def _build_meds_layout(self):
        """
        This fudges some things for now
        """

        # box sizes are even
        obj_data=self.obj_data

        # assuming we knew ncutout from the beginning
        #self.obj_data = self._make_resized_data(obj_data)
        self._set_start_rows_and_pixel_count()

    def write(self, filename, producer):
        """
        build the meds layout and write images from
        the producer
        """
        self._build_meds_layout()
        self._write_data(filename, producer)

    def _write_data(self, filename, producer):
        """
        run through and write cutouts from each SE file
        for each image
        """

        print("opening output MEDS file: '%s'" % filename)
        with fitsio.FITS(filename,'rw',clobber=True) as fits:
            self.fits=fits

            self._write_object_data()
            self._write_image_info()
            self._write_metadata()

            self._reserve_mosaic_images()

            for iobj,obj in enumerate(producer):
                self._write_object_cutouts(iobj,obj)

        print('output is in:',filename)

    def _write_object_cutouts(self, iobj, obj):
        """
        write the cutouts for the specified type
        """

        obj_data=self.obj_data
        nobj=obj_data.size
        assert iobj < nobj

        for cutout_type in self['cutout_types']:
            print('%d: writing %s cutouts' % (iobj,cutout_type))

            cutout_hdu = self._get_cutout_hdu(cutout_type)

            obs_list = self._extract_images(obj, cutout_type)

            ncut =len(obs_list)
            nexp = obj_data['ncutout'][iobj]
            if ncut != nexp:
                raise ValueError("expected %d cutouts, got %d" % (nexp,ncut))

            for icut, im_data in obs_list:
                image = imdata # some unpacking needs to be done here
                self._write_cutout(
                    iobj,
                    icut,
                    cutout_hdu,
                    im_data,
                    cutout_type,
                )

def make_test_image_info():
    return numpy.zeros(10, dtype=[('dummy','i4')])

def make_test_obj_data():
    n=10
    ncutout=10
    box_size=48
    pixscale=0.2

    dtype=[
        ('id','i8'),
        ('number','i8'),
        ('ra','f8'),
        ('dec','f8'),
        ('box_size','i4'),
        ('ncutout','i4'),

        ('file_id','i4',ncutout),

        # will be filled in
        ('start_row','i8',ncutout),

        ('orig_row','f8',ncutout),
        ('orig_col','f8',ncutout),
        ('orig_start_row','i8',ncutout),
        ('orig_start_col','i8',ncutout),

        ('cutout_row','f8',ncutout),
        ('cutout_col','f8',ncutout),

        ('dudcol','f8',ncutout),
        ('dudrow','f8',ncutout),
        ('dvdcol','f8',ncutout),
        ('dvdrow','f8',ncutout),
    ]

    data = numpy.zeros(
        n,
        dtype=dtype,
    )

    data['id'] = numpy.arange(n)
    data['number'] = data['id']
    data['ra'] = 200.0
    data['dec'] = -15.0
    data['box_size'] = box_size

    for n in ['file_id','start_row',
              'orig_row','orig_col',
              'orig_start_row','orig_start_col']:
        data[n] = -9999

    data['ncutout'] = 10
    data['dudcol'] = pixscale
    data['dvdrow'] = pixscale

    for i in xrange(data.size):
        data['cutout_row'][i,:] = (box_size-1.0)/2.0
        data['cutout_col'][i,:] = (box_size-1.0)/2.0

    return data

def test():

    obj_data = make_test_obj_data()
    ii = make_test_image_info()
    
    maker = DMMedsMaker(
        obj_data,
        ii,
    )
    # dummy variable
    producer={}
    maker.write("test.fits", producer)

def test_producer():
    """
    """
    from .producer import LSSTProducer
    import lsst.daf.persistence as dafPersist

    bpath="/u/ki/boutigny/ki19/MACSJ2243/output/u/ki/boutigny/ki19/MACSJ2243/output/coadd_dir_cc"

    butler = dafPersist.Butler(bpath)
    tract = 0
    patch = "1,4"
    filter = "r"

    producer = LSSTProducer(
        butler,
        tract,
        patch,
        filter,
    )
    
    arr=stamp.getMaskedImage().getImage().getArray()
    var=stamp.getMaskedImage().getVariance().getArray()
    mask=stamp.getMaskedImage().getMask().getArray()
    psfobj=stamp.getPsf()
    psfim = psfobj.computeKernelImage(pos1, pos2)

    pass
    
    # create the producer

    # make a catalog

    # copy info into numpy array

    # learn how to extract image info


if __name__=="__main__":
    test()
