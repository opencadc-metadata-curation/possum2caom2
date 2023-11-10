"""

POSSUM_archive_thumbnails.py:

Generates thumbnails for the CADC archive. Designed primarily for frequency
cubes, but with some thought given to making it work for FDF cubes as well.


Created on Fri Nov  3 09:15:05 2023
@author: cvaneck
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from astropy.nddata import block_reduce
import argparse
import re
import os
import glob

from caom2 import ProductType, ReleaseType
from caom2pipe.manage_composable import PreviewVisitor

# def command_line():
#     """Function for calling from the command line. Takes input directory, optional
#     output filename or overwrite order.
#     """

#     descStr = """Create thumbnail and preview images for a specified FITS file
#     or directory.
#     """

#     parser = argparse.ArgumentParser(description=descStr,
#                                  formatter_class=argparse.RawTextHelpFormatter)
#     parser.add_argument("filename",metavar="FITS_file",
#                         help="Input FITS file or directory.")
#     args = parser.parse_args()

#     if os.path.isfile(args.filename):
#         create_thumbnail_and_preview(args.filename)
#     elif os.path.isdir(args.filename):
#         files=glob.glob(os.path.join(args.filename,'*.fits'))
#         for file in files:
#             create_thumbnail_and_preview(file)


def create_thumbnail_and_preview(filename):
    """Create a 1024x1024 preview and 256x256 thumbnail image of the supplied
    FITS file. Assumes input file is a 2048x2048 POSSUM tile.
    If the file is part of a Stokes Q/U or real/imaginary pair,
    the previous image is of polarized intensity.
    If the file is a cube, it is collapsed along the third axis.
    """
    # Get data and header:
    data = pf.getdata(filename)

    # Check if image is part of a Q,U pair:
    m1 = re.search('_[qu]_', filename)
    m2 = re.search('_real_|_im_', filename)
    # If so, grab other file and generate polarized intensity thumbnail
    if m1 is not None:
        if m1[0] == '_q_':
            data_alt = pf.getdata(filename.replace('_q_', '_u_'))
        elif m1[0] == '_u_':
            data_alt = pf.getdata(filename.replace('_u_', '_q_'))

        data = collapse_qu_cubes(data, data_alt, mode='mean')

    if m2 is not None:
        if m2[0] == '_real_':
            data_alt = pf.getdata(filename.replace('_real_', '_im_'))
        elif m2[0] == '_im_':
            data_alt = pf.getdata(filename.replace('_im_', '_real_'))

        data = collapse_qu_cubes(data, data_alt, mode='maximum')

    # If data is still 3D (i.e., a cube that isn't part of a Q,U pair), collapse
    # the 3rd dimension.
    data = np.squeeze(data)
    if data.ndim == 3:
        data = np.nanmean(data, axis=0)

    outfile = filename.replace('.fits', '_preview.jpg')
    array_to_jpeg(data, outfile, downscale_factor=2)

    outfile = filename.replace('.fits', '_thumbnail.jpg')
    array_to_jpeg(data, outfile, downscale_factor=8)


# def collapse_qu_cubes(Qdata,Udata,mode='mean'):
#         """Converts Q and U into polarized intensity, then collapses along the
#         third axis. Two collapse modes are supported: 'mean' and 'maximum'.
#         The mean mode does no correction for polarization bias, so it'll end
#         up with a positive bias (which should hopefully go away when converting
#         data values to colors).
#         Assumes input data is in FITS axis ordering (frequency/Faraday depth first).
#         """
#         Pdata = np.squeeze(np.sqrt(Qdata**2+Udata**2))
#         if mode == 'mean':
#             data=np.mean(Pdata,axis=0)
#         elif mode == 'maximum':
#             data=np.max(Pdata,axis=0)
#         else:
#             raise Exception('Invalid collapse mode specified. Only "mean" and "maximum" supported.')

#         return data


class PossumPreview(PreviewVisitor):
    def __init__(self, **kwargs):
        super().__init__(ReleaseType.META, **kwargs)

    def generate_plots(self, obs_id):
        """Create a 1024x1024 preview and 256x256 thumbnail image of the supplied
        FITS file. Assumes input file is a 2048x2048 POSSUM tile.
        If the file is part of a Stokes Q/U or real/imaginary pair,
        the previous image is of polarized intensity.
        If the file is a cube, it is collapsed along the third axis.
        """
        self._logger.debug(f'Begin generate_plots for {obs_id}')
        count = 0
        # Get data and header:
        data = pf.getdata(self._science_fqn)

        # Check if image is part of a Q,U pair:
        m1 = re.search('_[qu]_', self._science_fqn)
        m2 = re.search('_real_|_im_', self._science_fqn)
        # If so, grab other file and generate polarized intensity thumbnail
        if m1 is not None:
            if m1[0] == '_q_':
                data_alt = pf.getdata(self._science_fqn.replace('_q_', '_u_'))
            elif m1[0] == '_u_':
                data_alt = pf.getdata(self._science_fqn.replace('_u_', '_q_'))

            data = self._collapse_qu_cubes(data, data_alt, mode='mean')

        if m2 is not None:
            if m2[0] == '_real_':
                data_alt = pf.getdata(self._science_fqn.replace('_real_', '_im_'))
            elif m2[0] == '_im_':
                data_alt = pf.getdata(self._science_fqn.replace('_im_', '_real_'))

            data = self._collapse_qu_cubes(data, data_alt, mode='maximum')

        # If data is still 3D (i.e., a cube that isn't part of a Q,U pair), collapse
        # the 3rd dimension.
        data = np.squeeze(data)
        if data.ndim == 3:
            data = np.nanmean(data, axis=0)

        self._array_to_jpeg(data, self._preview_fqn, downscale_factor=2)
        count += 1
        self.add_preview(self._storage_name.prev_uri, self._storage_name.prev, ProductType.PREVIEW, ReleaseType.DATA)
        self.add_to_delete(self._preview_fqn)

        self._array_to_jpeg(data, self._thumb_fqn, downscale_factor=8)
        count += 1
        self.add_preview(
            self._storage_name.thumb_uri, self._storage_name.thumb, ProductType.THUMBNAIL, ReleaseType.META
        )
        self.add_to_delete(self._thumb_fqn)
        self._logger.debug('End generate_plots')
        return count

    def _array_to_jpeg(self, data, jpeg_filename, downscale_factor=8):
        """
        Convert data array to greyscale JPEG file. Must already be collapsed to 2D

        Parameters:
        - data: 2D array of image values
        - jpeg_filename: Path where the output JPEG will be saved.
        - downscale_factor: Factor by which the image should be downscaled.
        """

        # Check if there are any singleton dimensions and squeeze them out
        data = np.squeeze(data)
        # Replace NaNs because matplotlib can't handle them:
        data_cleaned = np.where(np.isfinite(data), data, 0.0)

        # Downsample the 2D image
        data_downsampled = block_reduce(data_cleaned, downscale_factor)
        # Save the downsampled image as a JPEG
        plt.imsave(
            jpeg_filename,
            data_downsampled,
            vmin=np.percentile(data_downsampled, 1),
            vmax=np.percentile(data_downsampled, 99),
            cmap='gray',
            origin='lower',
        )

    def _collapse_qu_cubes(self, Qdata, Udata, mode='mean'):
        """Converts Q and U into polarized intensity, then collapses along the
        third axis. Two collapse modes are supported: 'mean' and 'maximum'.
        The mean mode does no correction for polarization bias, so it'll end
        up with a positive bias (which should hopefully go away when converting
        data values to colors).
        Assumes input data is in FITS axis ordering (frequency/Faraday depth first).
        """
        Pdata = np.squeeze(np.sqrt(Qdata**2 + Udata**2))
        if mode == 'mean':
            data = np.nanmean(Pdata, axis=0)
        elif mode == 'maximum':
            data = np.nanmax(Pdata, axis=0)
        else:
            raise Exception('Invalid collapse mode specified. Only "mean" and "maximum" supported.')

        return data


def visit(observation, **kwargs):
    previewer = PossumPreview(**kwargs)
    return previewer.visit(observation)
