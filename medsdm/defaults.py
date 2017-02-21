"""
default configuration

These get overridden by user configuration on constructing a DESMEDSMaker
object

"""

# sextractor names
default_config = {
    # need this for lsst
    'refband':'i',
    
    # types of cutout images to make
    'cutout_types': ['image','weight','seg','bmask'],

    # for fpacking the file.  Put into header
    'fpack_dims': [10240,1],

    # need this for lsst
    'magzp_ref':30.0,

    # 2**N or 3*2**N for fast FFTS
    'allowed_box_sizes': [
        2,3,4,6,8,12,16,24,32,48,
        64,96,128,192,256,
        384,512,768,1024,1536,
        2048,3072,4096,6144
    ],
    'min_box_size': 32,
    'max_box_size': 256,

}
