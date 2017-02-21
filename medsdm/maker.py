from __future__ import print_function
import numpy
import meds
import fitsio

import lsst.afw.geom as afwGeom

import esutil as eu
from meds import util

from .defaults import default_config

class DMMedsMaker(meds.MEDSMaker):
    """
    - check box sizes are 2^N 3*2^N ?
    - we don't have access to files
    - so the image_info will currently be blank
    - need scale to put all on same zero point, if not already done
    """

    def __init__(self,
                 producer,
                 config=None,
                 meta_data=None):

        self._load_config(config)
        self._set_extra_config()

        # make copies since we may alter some things
        self._set_meta_data(meta_data)

        self.producer = producer

        print("    getting the catalog")
        obj_data = self.producer.getCatalog()
        self._set_obj_data(obj_data)

        # fake the image info for now
        self.image_info=util.get_image_info_struct(
            10,
            10,
        )

    def _build_meds_layout(self):
        """
        This fudges some things for now
        """

        # box sizes are even
        obj_data=self.obj_data

        # assuming we knew ncutout from the beginning
        #self.obj_data = self._make_resized_data(obj_data)
        self._set_start_rows_and_pixel_count()

    def write(self, filename):
        """
        build the meds layout and write images from
        the producer
        """
        self._build_meds_layout()
        self._write_data(filename)

    def _write_data(self, filename):
        """
        run through and write cutouts from each SE file
        for each image
        """

        print("    opening output MEDS file: '%s'" % filename)

        with fitsio.FITS(filename,'rw',clobber=True) as fits:
            self.fits=fits

            self._reserve_mosaic_images()

            # write images and also fill in object data structure
            nobj=self.obj_data.size
            for iobj in xrange(nobj):
                self._write_object_cutouts(iobj)

            # We filled this on the fly, write last
            self._write_object_data()
            self._write_image_info()
            self._write_metadata()

        print('    output is in:',filename)

    def _write_object_cutouts(self, iobj):
        """
        write the cutouts for the specified type
        """

        obj_data=self.obj_data
        nobj=obj_data.size
        assert iobj < nobj

        image_data = self.producer.getStamps(obj_data[iobj])
        ncut =len(image_data)
        nexp = obj_data['ncutout'][iobj]
        if ncut != nexp:
            raise ValueError("expected %d cutouts, got %d" % (nexp,ncut))


        # fill in obj_data for the stamps
        self._fill_obj_data(iobj, image_data)

        # write image data
        for cutout_type in self['cutout_types']:
            print('    %d: writing %s cutouts' % (iobj,cutout_type))

            cutout_hdu = self._get_cutout_hdu(cutout_type)

            for icut, idata in enumerate(image_data):
                print("    icut:",icut)
                stamp, orig_pos = idata

                if stamp is None:
                    print("    warning: stamp is None")
                    continue

                im_data = self._extract_image(stamp, cutout_type, obj_data['box_size'][iobj])

                self._write_cutout(
                    iobj,
                    icut,
                    cutout_hdu,
                    im_data,
                    cutout_type,
                )

    def _write_cutout(self,
                      iobj,
                      icut,
                      cutout_hdu,
                      im_data,
                      cutout_type):
        """
        extract a cutout and write it to the mosaic image
        """
        start_row = self.obj_data['start_row'][iobj,icut]

        cutout_hdu.write(im_data, start=start_row)


    def _extract_image(self, stamp, cutout_type, dim):
        mim  = stamp.getMaskedImage()
        if cutout_type == 'image':
            data = mim.getImage().getArray()
        elif cutout_type == 'bmask':
            data = mim.getMask().getArray()
        elif cutout_type=='seg':
            raise NotImplementedError("implement seg map")
        elif cutout_type=='weight':
            var   = mim.getVariance().getArray()
            data = var.copy()
            data[:,:]=0
            w=numpy.where(var > 0)
            if data[0].size > 0:
                data[w] = 1.0/var[w]

        else:
            raise NotImplementedError("bad image cutout_type: '%s'" % cutout_type)


        eshape=(dim,dim)
        if data.shape != eshape:
            raise ValueError("expected dims %s, got %s" % (eshape,data.shape))

        tn='%s_dtype' % cutout_type
        dtype=self[tn]
        data = numpy.ascontiguousarray(
            data,
            dtype=dtype,
        )

        return data

    def _fill_obj_data(self, iobj, image_data):

        obj_data=self.obj_data
        for icut,idata in enumerate(image_data):
            stamp, orig_pos = idata
            if stamp is None:
                continue

            wcs = stamp.getWcs().linearizePixelToSky(
                orig_pos,
                afwGeom.arcseconds,
            )
            jacobian = wcs.getLinear().getMatrix()

            obj_data['orig_row'][iobj,icut] = orig_pos.getY()
            obj_data['orig_col'][iobj,icut] = orig_pos.getX()


            orig_start = stamp.getXY0()
            obj_data['orig_start_row'][iobj,icut] = orig_start.getY()
            obj_data['orig_start_col'][iobj,icut] = orig_start.getX()

            # location in the cutout
            cutout_cen = orig_pos - orig_start

            obj_data['cutout_row'][iobj,icut] = cutout_cen.getY()
            obj_data['cutout_col'][iobj,icut] = cutout_cen.getX()


            # TODO determine the actual convention
            obj_data['dudrow'][iobj,icut] = jacobian[0,0]
            obj_data['dudcol'][iobj,icut] = jacobian[0,1]
            obj_data['dvdrow'][iobj,icut] = jacobian[1,0]
            obj_data['dvdcol'][iobj,icut] = jacobian[1,1]


    def _get_full_obj_data(self, obj_data):
        """
        For the dmstack maker, we expect ncutout to be
        in the input data
        """

        nmax = obj_data['ncutout'].max()
        if nmax < 2:
            nmax = 2

        self._set_extra_fields(obj_data, nmax)

        nobj = obj_data.size
        new_obj_data = util.get_meds_output_struct(
            nobj,
            nmax,
            extra_fields=self['extra_fields'],
        )
        eu.numpy_util.copy_fields(obj_data, new_obj_data)

        return new_obj_data

    def _get_minimal_meds_input(self):
        extra_fields=[('ncutout','i4')]
        return get_meds_input_struct(1, extra_fields=extra_fields)

    def _load_config(self, config):
        """
        load the default config, then load the input config
        """


        # now config for this class
        # first the defaults
        this_config = {}
        this_config.update(default_config)

        # now override
        if config is not None:
            this_config.update(config)

        # first load the defaults from the parent
        super(DMMedsMaker,self)._load_config(this_config)


def make_test_image_info():
    # this will set ncutout_max
    return numpy.zeros(10, dtype=[('dummy','i4')])

def make_test_obj_data(image_info):
    n=10
    ncutout=10
    ncutout_max=image_info.size

    box_size=48
    pixscale=0.2

    extra_fields=[('number','i8')]
    data = meds.util.get_meds_output_struct(
        n,
        ncutout_max,
        extra_fields=extra_fields,
    )

    data['id'] = numpy.arange(n)
    data['number'] = data['id']
    data['ra'] = 200.0
    data['dec'] = -15.0
    data['box_size'] = -9999 # to be filled i

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

def test(limit=10):
    from .producer import test_make_producer
    producer = test_make_producer(limit=limit)

    config={
        'cutout_types':['image','weight','bmask'],
    }
    maker = DMMedsMaker(producer, config=config)
    maker.write("test.fits")


if __name__=="__main__":
    test()
