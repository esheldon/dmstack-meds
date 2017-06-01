"""
default configuration

These get overridden by user configuration on constructing a DESMEDSMaker
object

"""

# sextractor names
DEFAULT_MAKER_CONFIG = {
    # need this for lsst
    #'refband':'i',
    
    # types of cutout images to make. Note psf is
    # automatically made
    'cutout_types': ['image','weight','seg','bmask'],

    # currently must be True, since there
    # are no SE seg maps
    'fake_se_seg':True,

    # for fpacking the file.  Put into header
    'fpack_dims': [10240,1],

    # need this for lsst
    #'magzp_ref':30.0,

    'variable_psf_box_size': True,
}

DEFAULT_PRODUCER_CONFIG = {
    'fake_seg_radius':5,
    'min_box_size':32,
    'max_box_size':256,
    'radius_factor':5.0,
    'include_coadd':True,
}


