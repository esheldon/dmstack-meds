from __future__ import print_function
try:
    xrange
except:
    xrange=range

import os
import numpy
import esutil as eu

import meds
import fitsio

import lsst.afw.geom as afwGeom

from .defaults import DEFAULT_MAKER_CONFIG
from .pbar import prange

class DMMedsMaker(meds.MEDSMaker):
    """
    Wrapper class for the meds.MEDSMaker to adapt
    outputs from the LSST dmstack
    """

    def __init__(self,
                 producer,
                 config=None,
                 meta_data=None):

        # fake this
        self.psf_data={}

        self._load_config(config)
        self._set_extra_config()

        # make copies since we may alter some things
        self._set_meta_data(meta_data)

        self.producer = producer

        print("    getting the catalog")
        obj_data = self.producer.getCatalog()
        self._set_obj_data(obj_data)

        self._begin_image_info()
        # fake the image info for now
        #self.image_info=meds.util.get_image_info_struct(
        #    10,
        #    10,
        #)

    def write(self, filename):
        """
        build the meds layout and write images from
        the producer
        """
        self._build_meds_layout()
        self._write_data(filename)

    def _build_meds_layout(self):
        """
        build the layout for the FITS file on disk
        """

        # box sizes are even
        obj_data=self.obj_data

        # assuming we knew ncutout from the beginning
        #self.obj_data = self._make_resized_data(obj_data)
        self._set_start_rows_and_pixel_count()

        self._set_psf_layout()


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
            self.current_psf_position=0
            print("writing cutouts")
            for iobj in prange(nobj):
                try:
                    self._write_object_cutouts(iobj)
                except:
                    print("Failed generating cutouts for object %d with ID %d" % (iobj, self.obj_data[iobj]["id"]))
                    raise

            # We filled this on the fly, write last
            self._write_object_data()
            self._write_image_info()
            self._write_metadata()

        print('    output is in:',filename)

    def _write_image_info(self):
        """
        combine the info structs and write to disk
        """

        self.image_info = numpy.hstack(self._image_info_list)
        super(DMMedsMaker,self)._write_image_info()

    def _write_object_cutouts(self, iobj):
        """
        write the cutouts for the specified type
        """

        obj_data=self.obj_data
        nobj=obj_data.size
        assert iobj < nobj

        image_data = self.producer.getStamps(obj_data[iobj])

        ncut = len(image_data)
        nexp = obj_data['ncutout'][iobj]
        if ncut != nexp:
            raise ValueError("expected %d cutouts, got %d" % (nexp,ncut))


        # fill in obj_data for the stamps
        self._fill_obj_data(iobj, image_data)

        box_size = obj_data['box_size'][iobj]
        # write image data
        for cutout_type in self['cutout_types'] + ['psf']:
            #print('    %d: writing %s cutouts' % (iobj,cutout_type))

            cutout_hdu = self._get_cutout_hdu(cutout_type)

            for icut, idata in enumerate(image_data):
                stamp = idata['stamp']
                seg_map = idata['seg_map']

                if stamp is None:
                    print("    stamp",icut,"is None")
                    continue

                if idata['image_id'] not in self._idlist:
                    self._append_image_info(idata)

                obj_data['file_id'][iobj,icut] = \
                        self._idlist.index(idata['image_id'])

                if cutout_type == 'seg':
                    if seg_map is None and icut != 0:
                        assert self['fake_se_seg']
                        # grab the image and make a fake
                        # seg map like that

                        im_data = numpy.zeros( [box_size]*2, dtype='i4')
                    else:
                        im_data = numpy.array(seg_map.array, dtype='i4', copy=False)

                elif cutout_type=='psf':
                    # psfs are variable in size
                    im_data = self._extract_psf_image(stamp, idata['image_pos'])

                    obj_data['psf_box_size'][iobj,icut] = im_data.shape[0]
                    obj_data['psf_start_row'][iobj,icut] = self.current_psf_position

                    # increment for the next write
                    self.current_psf_position += im_data.size

                else:
                    im_data = self._extract_image(
                        stamp,
                        cutout_type,
                        box_size,
                    )

                self._write_cutout(
                    iobj,
                    icut,
                    cutout_hdu,
                    im_data,
                    cutout_type,
                )

    def _extract_psf_image(self, stamp, orig_pos):
        """
        get the psf associated with this stamp
        """
        psfobj=stamp.getPsf()
        psfim = psfobj.computeKernelImage(orig_pos).array.astype('f4')
        psfim = numpy.array(psfim, dtype='f4', copy=False)

        d=psfim.shape
        if d[0] != d[1]:
            #print("trimming psf",d)
            if d[0] > d[1]:
                bigger=d[0]
                smaller=d[1]
            else:
                bigger=d[1]
                smaller=d[0]

            diff = bigger-smaller
            assert (diff % 2) == 0

            beg = diff//2
            end = bigger - diff//2

            if d[0] > d[1]:
                psfim = psfim[beg:end,:]
            else:
                psfim = psfim[:,beg:end]

        return psfim

    def _write_cutout(self,
                      iobj,
                      icut,
                      cutout_hdu,
                      im_data,
                      cutout_type):
        """
        extract a cutout and write it to the mosaic image
        """

        if cutout_type=='psf':
            start_row = self.obj_data['psf_start_row'][iobj,icut]
        else:
            start_row = self.obj_data['start_row'][iobj,icut]

        cutout_hdu.write(im_data, start=start_row)


    def _extract_image(self, stamp, cutout_type, dim):
        """
        wrapper to extract image,bmask or weight images
        """
        if cutout_type == 'image':
            data = stamp.image.array
        elif cutout_type == 'bmask':
            data = stamp.mask.array
        elif cutout_type=='weight':
            data = self._extract_weight(stamp)
        else:
            raise NotImplementedError("bad image cutout_type "
                                      "for _extract_image: '%s'" % cutout_type)


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

    def _extract_weight(self, stamp):
        """
        extract a weight map

        Areas with NO_DATA will get zero weight.

        Because the variance map includes the full poisson variance, which
        depends on the signal, we instead extract the median of the parts of
        the image without NO_DATA set
        """
        var_image  = stamp.variance.array
        weight = var_image.copy()

        weight[:,:]=0

        zlogic = var_image > 0

        no_data_logic = numpy.logical_not(
            stamp.mask.array & stamp.mask.getPlaneBitMask("NO_DATA")
        )
        w=numpy.where(zlogic & no_data_logic)

        if w[0].size > 0:
            medvar = numpy.median(var_image[w])
            weight[w] = 1.0/medvar

        return weight

    def _fill_obj_data(self, iobj, image_data):
        """
        fill in the data  for each object, such as
        wcs jacobians and info about location of the postage
        stamps, etc.
        """
        obj_data=self.obj_data
        for icut,idata in enumerate(image_data):
            stamp = idata['stamp']
            if stamp is None:
                continue

            orig_pos = idata['image_pos']

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
        new_obj_data = meds.util.get_meds_output_struct(
            nobj,
            nmax,
            extra_fields=self['extra_fields'],
        )
        eu.numpy_util.copy_fields(obj_data, new_obj_data)

        return new_obj_data

    def _set_extra_fields(self, obj_data, nmax):
        self['extra_fields'] = [
            ('number','i4'),
            ('psf_box_size','i8',nmax),
            ('psf_start_row','i8',nmax),
        ]

    def _get_minimal_meds_input(self):
        extra_fields=[('ncutout','i4')]
        return meds.util.get_meds_input_struct(1, extra_fields=extra_fields)

    def _append_image_info(self, idata):
        info = self._get_empty_image_info()
        info['image_id'] = idata['image_id']
        if 'visit' in idata:
            info['visit'] = idata['visit']
            info['ccd']   = idata['ccd']
        else:
            info['visit'] = -9999
            info['ccd']   = -9999

        self._image_info_list.append(info)
        self._idlist.append(idata['image_id'])


    def _begin_image_info(self):
        self._idlist=[]
        self._image_info_list=[]

        cat = self.producer.getCatalog()
        stamps = self.producer.getStamps(cat[0])
        coadd_idata = stamps[0]

        self._append_image_info(coadd_idata)


    def _get_empty_image_info(self):
        return meds.util.get_image_info_struct(
            1,
            10,
            extra_dtype=[('visit','i8'),('ccd','i2')],
        )



    def _set_psf_layout(self):
        """
        set the box sizes and start row for each psf image

        we don't yet actually know the box sizes, so we will assume
        the coadd psf is maximal, with some padding
        """

        obj_data=self.obj_data

        producer = self.producer

        cat = producer.getCatalog()
        stamps = producer.getStamps(cat[0])

        sdata = stamps[0]
        psfobj=sdata['stamp'].getPsf()
        psfim = psfobj.computeKernelImage(sdata['image_pos']).array

        # some padding
        psf_size = max(psfim.shape)+2

        # now assume all are the same size for reserving the
        # data on disk. Not all pixels will be used

        #obj_data['psf_box_size'] = psf_size
        total_psf_pixels = 0
        psf_npix = psf_size*psf_size

        #psf_start_row = 0
        for i in xrange(obj_data.size):
            for j in xrange(obj_data['ncutout'][i]):
                #obj_data['psf_start_row'][i,j] = psf_start_row

                #psf_start_row += psf_npix
                total_psf_pixels += psf_npix


        self.total_psf_pixels = total_psf_pixels

    '''
    def _get_psf_dtype(self, nmax):
        """
        LSST psfs are variable in size
        """
        return [
            ('psf_box_size','i8',nmax),
            ('psf_start_row','i8',nmax),
        ]
    '''

    def _load_config(self, config):
        """
        load the default config, then load the input config
        """


        # now config for this class
        # first the defaults
        this_config = {}
        this_config.update(DEFAULT_MAKER_CONFIG)

        # now override
        if config is not None:
            this_config.update(config)

        # first load the defaults from the parent
        super(DMMedsMaker,self)._load_config(this_config)

def test(filter="HSC-I",tract=8766, patch="4,4", limit=10):
    from .producer import test_make_producer

    pstr = patch.replace(',','')
    outdir='output'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    fstr=filter[4:4+1].lower()
    fname='test-medsdm-%s-tract%06d-patch%s' % (fstr,tract, pstr)

    if limit is not None:
        fname = '%s-limit%d' % (fname,limit)

    fname='%s.fits' % fname
    fname=os.path.join(outdir, fname)
    print("writing:",fname)

    producer_config={
        'fake_seg_radius':5,
        #'min_box_size':48,
    }

    producer = test_make_producer(
        filter,
        tract=tract,
        patch=patch,
        limit=limit,
        config=producer_config,
    )

    maker = DMMedsMaker(producer)
    maker.write(fname)

if __name__=="__main__":
    import sys
    import yaml
    from argparse import ArgumentParser, RawTextHelpFormatter
    import lsst.daf.persistence as dafPersist
    from .producer import LSSTProducer

    usage = """
    Extracts MEDS file cutouts from the DM coadds and individual visits for a
    given (tract, patch).
    """
    parser = ArgumentParser(description=usage,
                         formatter_class=RawTextHelpFormatter)
    parser.add_argument('--prefix', type=str, default='./medsdm',
                        help='Output prefix to use for the meds files.')
    parser.add_argument('--config', type=str, default=None,
                        help='Yaml configuration file for medsdm')
    parser.add_argument('repo', type=str,
                     help='Filepath to LSST DM Stack Butler repository.')
    parser.add_argument('tract', type=int,
                     help='Skymap tract to process.')
    parser.add_argument('patch', type=str,
                     help='patch to process (e.g. 0,1).')
    parser.add_argument('filter', type=str,
                     help='Filter to process (e.g. i).')
    args = parser.parse_args(sys.argv[1:])

    # Extracts output filename
    fname=args.prefix+'-%s-tract%06d-patch%s' % (args.filter, args.tract, args.patch)

    # Extracts configuration if different from default
    if args.config is not None:
        with open(args.config, 'r') as ymlfile:
            config = yaml.load(ymlfile)
        config_producer = config['producer']
        config_maker = config['maker']
    else:
        config_producer = None
        config_maker =  None

    # Create producer using default configuration
    butler =  dafPersist.Butler(args.repo)
    producer = LSSTProducer(butler,
                            tract=args.tract,
                            patch=args.patch,
                            filter=args.filter,
                            config= args.config['producer'] if args.config is not None else None)

    # Create maker instance, with default configuration for now
    maker = DMMedsMaker(producer,
                        config= args.config['maker'] if args.config is not None else None)
    maker.write(fname)
