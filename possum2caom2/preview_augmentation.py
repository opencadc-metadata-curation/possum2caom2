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
import re
import os

from caom2 import ProductType, ReleaseType
from caom2pipe.manage_composable import PreviewVisitor


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
                pair_file = self._science_fqn.replace('_q_', '_u_')
                self._get_pair_file(pair_file, self.storage_name.file_uri.replace('_q_', '_u_'))
            elif m1[0] == '_u_':
                pair_file = self._science_fqn.replace('_u_', '_q_')
                self._get_pair_file(pair_file, self.storage_name.file_uri.replace('_u_', '_q_'))
            data_alt = pf.getdata(pair_file)

            data = self._collapse_qu_cubes(data, data_alt, mode='mean')

        if m2 is not None:
            if m2[0] == '_real_':
                pair_file = self._science_fqn.replace('_real_', '_im_')
                self._get_pair_file(pair_file, self.storage_name.file_uri.replace('_real_', '_im_'))
            elif m2[0] == '_im_':
                pair_file = self._science_fqn.replace('_im_', '_real_')
                self._get_pair_file(pair_file, self.storage_name.file_uri.replace('_im_', '_real_'))
            data_alt = pf.getdata(pair_file)
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

    def _get_pair_file(self, pair_fqn, pair_uri):
        if not os.path.exists(pair_fqn):
            self._clients.data_client.get(os.path.dirname(pair_fqn), pair_uri)


def visit(observation, **kwargs):
    previewer = PossumPreview(**kwargs)
    return previewer.visit(observation)
