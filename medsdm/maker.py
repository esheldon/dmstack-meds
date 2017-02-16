import meds

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
    def blah(self):
        pass
