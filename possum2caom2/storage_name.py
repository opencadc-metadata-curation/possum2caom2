# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2024.                            (c) 2024.
#  Government of Canada                 Gouvernement du Canada
#  National Research Council            Conseil national de recherches
#  Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#  All rights reserved                  Tous droits réservés
#
#  NRC disclaims any warranties,        Le CNRC dénie toute garantie
#  expressed, implied, or               énoncée, implicite ou légale,
#  statutory, of any kind with          de quelque nature que ce
#  respect to the software,             soit, concernant le logiciel,
#  including without limitation         y compris sans restriction
#  any warranty of merchantability      toute garantie de valeur
#  or fitness for a particular          marchande ou de pertinence
#  purpose. NRC shall not be            pour un usage particulier.
#  liable in any event for any          Le CNRC ne pourra en aucun cas
#  damages, whether direct or           être tenu responsable de tout
#  indirect, special or general,        dommage, direct ou indirect,
#  consequential or incidental,         particulier ou général,
#  arising from the use of the          accessoire ou fortuit, résultant
#  software.  Neither the name          de l'utilisation du logiciel. Ni
#  of the National Research             le nom du Conseil National de
#  Council of Canada nor the            Recherches du Canada ni les noms
#  names of its contributors may        de ses  participants ne peuvent
#  be used to endorse or promote        être utilisés pour approuver ou
#  products derived from this           promouvoir les produits dérivés
#  software without specific prior      de ce logiciel sans autorisation
#  written permission.                  préalable et particulière
#                                       par écrit.
#
#  This file is part of the             Ce fichier fait partie du projet
#  OpenCADC project.                    OpenCADC.
#
#  OpenCADC is free software:           OpenCADC est un logiciel libre ;
#  you can redistribute it and/or       vous pouvez le redistribuer ou le
#  modify it under the terms of         modifier suivant les termes de
#  the GNU Affero General Public        la “GNU Affero General Public
#  License as published by the          License” telle que publiée
#  Free Software Foundation,            par la Free Software Foundation
#  either version 3 of the              : soit la version 3 de cette
#  License, or (at your option)         licence, soit (à votre gré)
#  any later version.                   toute version ultérieure.
#
#  OpenCADC is distributed in the       OpenCADC est distribué
#  hope that it will be useful,         dans l’espoir qu’il vous
#  but WITHOUT ANY WARRANTY;            sera utile, mais SANS AUCUNE
#  without even the implied             GARANTIE : sans même la garantie
#  warranty of MERCHANTABILITY          implicite de COMMERCIALISABILITÉ
#  or FITNESS FOR A PARTICULAR          ni d’ADÉQUATION À UN OBJECTIF
#  PURPOSE.  See the GNU Affero         PARTICULIER. Consultez la Licence
#  General Public License for           Générale Publique GNU Affero
#  more details.                        pour plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est
#  <http://www.gnu.org/licenses/>.      pas le cas, consultez :
#                                       <http://www.gnu.org/licenses/>.
#
#  $Revision: 4 $
#
# ***********************************************************************
#

from os.path import basename
from urllib import parse as parse
from caom2pipe.manage_composable import CadcException, StorageName


class PossumName(StorageName):
    """
    From AusSRC:
    PSM_pilot1_944MHz_18asec_2226-5552_11268_i_v1.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_q_v1.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_u_v1.fits

    From POSSUM Group:
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_ampPeakPIfitEff.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_coeff0err.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_coeff0.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_coeff1err.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_coeff1.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_dAmpPeakPIfit.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_dFDFcorMAD.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_dPhiPeakPIfit_rm2.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_dPolAngle0Fit_deg.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_FDF_imag_dirty.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_FDF_real_dirty.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_FDF_tot_dirty.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_lam0Sq_m2.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_mad_chanwidth.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_max_freq.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_min_freq.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_peakFDFimagFit.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_peakFDFrealFit.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_phiPeakPIfit_rm2.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_polAngle0Fit_deg.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_RMSF1D.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_RMSF_FWHM.fits
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_snrPIfit.fits

    Guessing on the rename for this:
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_RMSF_imag.fits.tar.gz
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_RMSF_real.fits.tar.gz
    to:
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_RMSF_imag.fits.fz
    PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_RMSF_real.fits.fz

    """

    POSSUM_NAME_PATTERN = '*'

    def __init__(self, entry):
        self._healpix_index = None
        self._spatial_resolution = None
        # the renamed files in a list
        self._stage_names = []
        super().__init__(file_name=basename(entry.replace('.header', '')), source_names=[entry])

    def __str__(self):
        temp = super().__str__()
        return f'{temp}\n          rename: {self.rename("")}\n   healpix_index: {self._healpix_index}\n      resolution: {self._spatial_resolution}'

    @property
    def file_uri(self):
        """The CADC Storage URI for the file."""
        current_scheme = self._get_scheme()
        return self._get_uri(self._file_name.replace('.gz', '').replace('.header', ''), current_scheme)

    @property
    def healpix_index(self):
        return self._healpix_index

    @property
    def is_bintable(self):
        return '_FDFs.fits' in self._file_name or '_spectra.fits' in self._file_name

    @property
    def prev(self):
        """The preview file name for the file."""
        return f'{self._obs_id}_{self._product_id}_prev.jpg'

    @property
    def spatial_resolution(self):
        return self._spatial_resolution

    @property
    def stage_names(self):
        return self._stage_names

    @property
    def thumb(self):
        """The thumbnail file name for the file."""
        return f'{self._obs_id}_{self._product_id}_prev_256.jpg'

    def rename(self, band, version='v1', resolution='20asec'):
        if 'asec' in self._file_id:
            self._logger.error('second time through')
            return self._file_name
        self._logger.error(self._file_name)
        # Jennifer West, 08-04-24
        # band1 == 944MHz, band2 == 1296MHz
        # is version == v1
        # resolution == 20asec
        x = self._file_name.split('.fits')
        bits = x[0].split('.')
        temp1 = '_'.join(ii for ii in bits)
        if 'band1' in self._file_name:
            y = temp1.replace('band1', f'944MHz_{resolution}')
        elif 'band2' in self._file_name:
            y = temp1.replace('band2', f'1296MHz_{resolution}')
        else:
            y = temp1
        temp2 = y.replace('POSSUM', 'PSM')
        if len(x) == 1:
            temp = f'{temp2}_{version}.fits'
        else:
            temp = f'{temp2}_{version}.fits{x[1]}'
        # insert = f'PSM_{band}_{resolution}'
        # temp2 = temp1.replace('PSM', insert)
        # if len(x) == 1:
        #     temp = f'{temp2}_{version}.fits'
        # else:
        #     temp = f'{temp2}_{version}.fits{x[1]}'
        self._stage_names.append(temp)
        return temp

    def set_destination_uris(self):
        for entry in self._source_names:
            temp = parse.urlparse(entry)
            base_name = basename(temp.path)
            current_scheme = StorageName.scheme
            if '_p3d_' in base_name or '_p1d_' in base_name:
                current_scheme = StorageName.preview_scheme
            if '.fits' in entry:
                self._destination_uris.append(
                    self._get_uri(base_name.replace('.gz', '').replace('.header', ''), current_scheme)
                )
            else:
                self._destination_uris.append(self._get_uri(base_name, current_scheme))

    def set_obs_id(self):
        # picking the common prefix, e.g. 944MHz_pilot1_18asec_2226-5552_11268, and then re-organize it a bit
        # leave off the "PSM" because collection is POSSUM
        self._logger.error(self._file_id)
        bits = self._file_id.split('_')
        if len(bits) > 2:
            if 'pilot' in self._file_id:
                self._obs_id = f'{bits[2]}_{bits[3]}_{bits[4]}_{bits[5]}_{bits[1]}'
                self._healpix_index = bits[5]
                self._spatial_resolution = float(bits[3].replace('asec', ''))
            else:
                if len(bits) == 8:
                    self._obs_id = '_'.join(ii for ii in bits[1:-1])
                    self._healpix_index = bits[-4]
                    self._spatial_resolution = float(bits[2].replace('asec', ''))
                else:
                    self._obs_id = '_'.join(ii for ii in bits[1:])
                    self._healpix_index = bits[-2]
                    self._spatial_resolution = float(bits[2].replace('asec', ''))
        else:
            bits = self._file_id.split('.')
            self._obs_id = f'{bits[1]}_{bits[2]}'
            self._healpix_index = bits[2]

    def set_product_id(self):
        bits = self._file_id.split('_')
        self._logger.error(self._file_id)
        if len(bits) > 2:
            index = 5
            if 'pilot' in self._file_id:
                index = 6
            elif len(bits) == 7 and bits[-1] == 'v1':
                self._logger.error(f'len bits {len(bits)} {self._file_id} {bits}')
                index = 5
            elif len(bits) == 8 and bits[-1] == 'v1':
                index = 6
            self._logger.error(f'index {index} {bits[index]}')
            if '_p3d_' in self._file_id:
                self._product_id = '3d_pipeline'
            elif '_p1d_' in self._file_id:
                # catalog in csv, spectra, FDF in BINTABLE
                # self._product_id = '1d_pipeline'
                self._product_id = self._file_id
            elif bits[index] == 'i':
                self._product_id = 'raw_i'
            elif bits[index] == 'q' or bits[index] == 'u':
                self._product_id = 'raw_qu'
            elif bits[index] == 't0' or bits[index] == 't1':
                # Cameron Van Eck - 23-10-23
                # “mfs_i_t0" or “multifrequencysynthesis_i_t0” for the image product ProductID
                self._product_id = f'multifrequencysynthesis_{bits[7]}_{bits[6]}'
            else:
                raise CadcException(f'Unexpected file naming pattern {self._file_id}')
        else:
            bits = self._file_id.split('.')
            index = 3
            if self._file_id.startswith('POSSUM'):
                index = 4
            if bits[index] == 'i':
                self._product_id = 'raw_i'
            elif bits[index] == 'q' or bits[index] == 'u':
                self._product_id = 'raw_qu'
            else:
                raise CadcException(f'Unexpected raw file naming pattern {self._file_id}')

    def un_name(self):
        """Undo the renaming from Pawsey Acacia, to be able to work backwards for sc2 Observation metadata."""
        return self._obs_id.split('asec_')[1].split('_v1')[0]

    def _get_scheme(self):
        if self._product_id in ['raw_i', 'raw_qu'] or self._product_id.startswith('multifrequencysynthesis_'):
            result = StorageName.scheme
        else:
            result = StorageName.preview_scheme
        return result
