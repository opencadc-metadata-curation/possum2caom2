# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2020.                            (c) 2020.
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

"""
This module implements the ObsBlueprint mapping.

Temporal WCS:

Cameron Van Eck - 17-10-23 - At present, most of the tiles don't have any observation date in the headers. A few do,
but they are likely not accurate (they somehow survived mosaicking?). Strongly prefer not to have observation dates
attached to these data -- since the relationship between obs. date and files is complicated, users should consult
our metadata database to determine which files are relevant.
"""

import logging
import traceback

from datetime import datetime, timedelta
from os.path import basename
from urllib import parse as parse

from astropy.io import fits

from caom2 import CalibrationLevel, DataProductType, ProductType, ReleaseType
from caom2utils.blueprints import _to_float
from caom2utils.wcs_parsers import FitsWcsParser
from caom2pipe.astro_composable import get_datetime_mjd
from caom2pipe import caom_composable as cc
from caom2pipe.manage_composable import CadcException, make_datetime, StorageName, ValueRepairCache


__all__ = ['PossumName', 'mapping_factory']


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
        super(PossumName, self).__init__(file_name=basename(entry.replace('.header', '')), source_names=[entry])

    @property
    def file_uri(self):
        """The CADC Storage URI for the file."""
        current_scheme = self._get_scheme()
        return self._get_uri(self._file_name.replace('.gz', '').replace('.header', ''), current_scheme)

    @property
    def prev(self):
        """The preview file name for the file."""
        return f'{self._obs_id}_{self._product_id}_prev.jpg'

    @property
    def thumb(self):
        """The thumbnail file name for the file."""
        return f'{self._obs_id}_{self._product_id}_prev_256.jpg'

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
        bits = self._file_id.split('_')
        self._obs_id = f'{bits[2]}_{bits[3]}_{bits[4]}_{bits[5]}_{bits[1]}'

    def set_product_id(self):
        bits = self._file_id.split('_')
        if '_p3d_' in self._file_id:
            self._product_id = '3d_pipeline'
        elif '_p1d_' in self._file_id:
            self._product_id = '1d_pipeline'
        elif bits[6] == 'i':
            self._product_id = 'raw_i'
        elif bits[6] == 'q' or bits[6] == 'u':
            self._product_id = 'raw_qu'
        elif bits[6] == 't0' or bits[6] == 't1':
            # Cameron Van Eck - 23-10-23
            # “mfs_i_t0" or “multifrequencysynthesis_i_t0” for the image product ProductID
            self._product_id = f'multifrequencysynthesis_{bits[7]}_{bits[6]}'
        else:
            raise CadcException(f'Unexcepted file naming pattern {self._file_id}')

    def _get_scheme(self):
        if self._product_id in ['raw_i', 'raw_qu'] or self._product_id.startswith('multifrequencysynthesis_'):
            result = StorageName.scheme
        else:
            result = StorageName.preview_scheme
        return result


class PossumValueRepair(ValueRepairCache):
    VALUE_REPAIR = {
        'chunk.custom.axis.axis.cunit': {
            'rad / m2': 'rad/m**2',
        }
    }

    def __init__(self):
        self._value_repair = PossumValueRepair.VALUE_REPAIR
        self._key = None
        self._values = None
        self._logger = logging.getLogger(self.__class__.__name__)


class Possum1DMapping(cc.TelescopeMapping):
    value_repair = PossumValueRepair()

    def __init__(self, storage_name, headers, clients, observable, observation, config):
        super().__init__(storage_name, headers, clients, observable, observation, config)
        # Cameron Van Eck - 23-10-23
        # Set release date to be 12 months after ingest. That’s the current POSSUM policy: data goes public 12
        # months after being generated. It doesn't have to be particularly precise: date of ingest + increment year by 1
        self._1_year_after = datetime.now() + timedelta(days=365)

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        self._logger.debug('Begin accumulate_bp.')
        super().accumulate_blueprint(bp)
        # JW - 17-10-23 - Use ASKAP
        bp.set('Observation.instrument.name', 'ASKAP')
        bp.set('Observation.metaRelease', self._1_year_after)
        bp.set('Observation.proposal.id', '_get_proposal_id()')
        bp.set_default('Observation.telescope.name', 'ASKAP')
        bp.set_default('Observation.telescope.geoLocationX', -2558266.717765)
        bp.set_default('Observation.telescope.geoLocationY', 5095672.176508)
        bp.set_default('Observation.telescope.geoLocationZ', -2849020.838078)
        bp.set('Plane.calibrationLevel', CalibrationLevel.CALIBRATED)
        bp.set('Plane.dataProductType', '_get_data_product_type()')
        bp.set('Plane.metaRelease', self._1_year_after)
        bp.set('Plane.dataRelease', self._1_year_after)
        bp.set('Artifact.productType', ProductType.SCIENCE)
        bp.set('Artifact.releaseType', ReleaseType.DATA)
        self._logger.debug('Done accumulate_bp.')

    def update(self, file_info):
        """Called to fill multiple CAOM model elements and/or attributes
        (an n:n relationship between TDM attributes and CAOM attributes).
        """
        self._logger.debug(f'Begin update for {self._observation.observation_id}.')
        try:
            super().update(file_info)
            Possum1DMapping.value_repair.repair(self._observation)
            self._logger.debug('Done update.')
            return self._observation
        except CadcException as e:
            tb = traceback.format_exc()
            self._logger.debug(tb)
            self._logger.error(e)
            self._logger.error(f'Terminating ingestion for {self._observation.observation_id}')
            return None

    def _get_data_product_type(self, ext):
        naxis = self._headers[ext].get('NAXIS')
        result = DataProductType.CUBE
        if naxis == 0:
            result = DataProductType.MEASUREMENTS
        elif naxis == 2:
            result = DataProductType.IMAGE
        elif naxis == 4:
            naxis3 = self._headers[ext].get('NAXIS3')
            naxis4 = self._headers[ext].get('NAXIS4')
            if naxis3 == 1 and naxis4 == 1:
                result = DataProductType.IMAGE
        return result

    def _get_position_resolution(self, ext):
        result = None
        # JW - 17-10-23 - Use either BMAJ or BMIN
        # Cameron Van Eck - 19-10-23 - Prefer BMAJ
        bmaj = self._headers[ext].get('BMAJ')
        if bmaj:
            # Cameron Van Eck - 23-10-23
            # FITS header value is in degrees, convert to arcseconds
            result = bmaj * 3600.0
        return result

    def _get_proposal_id(self, ext):
        # Cameron Van Eck - 23-10-23
        # For proposalID: All pilot data can have value “AS103". All full-survey data will have value “AS203”.
        result = 'AS203'
        if '_pilot' in self._storage_name.file_name:
            result = 'AS103'
        return result

    def _update_artifact(self, artifact):
        delete_these = []
        for part in artifact.parts.values():
            if len(part.chunks) == 0:
                delete_these.append(part.name)
            else:
                for chunk in part.chunks:
                    if (
                        chunk.custom is None
                        and chunk.energy is None
                        and chunk.observable is None
                        and chunk.polarization is None
                        and chunk.position is None
                        and chunk.time is None
                    ) or (  # handle the Taylor BINTABLE extension case
                        chunk.custom is None
                        and chunk.energy is None
                        and chunk.observable is None
                        and chunk.polarization is None
                        and chunk.position is None
                        and chunk.time is not None
                    ):
                        delete_these.append(part.name)
                        break

        for entry in delete_these:
            artifact.parts.pop(entry)
            self._logger.info(f'Deleting part {entry} from artifact {artifact.uri}')

    @staticmethod
    def _from_pc_to_cd(from_header, to_header):
        cd1_1 = from_header.get('CD1_1')
        # caom2IngestSitelle.py, l745
        # CW
        # Be able to handle any of the 3 wcs systems used
        if cd1_1 is None:
            pc1_1 = from_header.get('PC1_1')
            if pc1_1 is not None:
                cdelt1 = _to_float(from_header.get('CDELT1'))
                if cdelt1 is None:
                    cd1_1 = _to_float(from_header.get('PC1_1'))
                    cd1_2 = _to_float(from_header.get('PC1_2'))
                    cd2_1 = _to_float(from_header.get('PC2_1'))
                    cd2_2 = _to_float(from_header.get('PC2_2'))
                else:
                    cdelt2 = _to_float(from_header.get('CDELT2'))
                    cd1_1 = cdelt1 * _to_float(from_header.get('PC1_1'))
                    cd1_2 = cdelt1 * _to_float(from_header.get('PC1_2'))
                    cd2_1 = cdelt2 * _to_float(from_header.get('PC2_1'))
                    cd2_2 = cdelt2 * _to_float(from_header.get('PC2_2'))
                to_header['CD1_1'] = cd1_1
                to_header['CD1_2'] = cd1_2
                to_header['CD2_1'] = cd2_1
                to_header['CD2_2'] = cd2_2


class InputTileMapping(Possum1DMapping):
    def __init__(self, storage_name, headers, clients, observable, observation, config):
        super().__init__(storage_name, headers, clients, observable, observation, config)

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        self._logger.debug('Begin accumulate_bp.')
        super().accumulate_blueprint(bp)
        bp.set('Plane.calibrationLevel', CalibrationLevel.CALIBRATED)
        bp.clear('Plane.provenance.name')
        bp.add_attribute('Plane.provenance.name', 'ORIGIN')
        # JW - 17-10-23 - Use AusSRC for producer
        bp.set('Plane.provenance.producer', 'AusSRC')
        bp.set_default('Plane.provenance.reference', 'https://possum-survey.org/')
        bp.add_attribute('Plane.provenance.lastExecuted', 'DATE')
        bp.set_default('Plane.provenance.project', 'POSSUM')
        bp.configure_position_axes((1, 2))
        bp.set('Chunk.position.resolution', '_get_position_resolution()')

        bp.configure_energy_axis(3)
        bp.set_default('Chunk.energy.specsys', 'TOPOCENT')

        bp.configure_polarization_axis(4)
        self._logger.debug('Done accumulate_bp.')

    def _update_artifact(self, artifact):
        self._logger.debug(f'Begin _update_artifact for {artifact.uri}')
        super()._update_artifact(artifact)
        for part in artifact.parts.values():
            for chunk in part.chunks:
                if chunk.energy is not None:
                    # JW - 17-10-23 - remove restfrq
                    chunk.energy.restfrq = None
                self._update_chunk_position(chunk)
        self._logger.debug('End _update_artifact')

    def _update_chunk_position(self, chunk):
        self._logger.debug(f'Begin update_position_function for {self._storage_name.obs_id}')
        if chunk.position is not None:
            header = self._headers[0]
            cd1_1 = header.get('CD1_1')
            if cd1_1 is None:
                hdr = fits.Header()
                Possum1DMapping._from_pc_to_cd(header, hdr)
                for kw in [
                    'CDELT1',
                    'CDELT2',
                    'CRPIX1',
                    'CRPIX2',
                    'CRVAL1',
                    'CRVAL2',
                    'CTYPE1',
                    'CTYPE2',
                    'CUNIT1',
                    'CUNIT2',
                    'NAXIS',
                    'NAXIS1',
                    'NAXIS2',
                    'DATE-OBS',
                    'EQUINOX',
                ]:
                    hdr[kw] = header.get(kw)
                wcs_parser = FitsWcsParser(hdr, self._storage_name.obs_id, 0)
                wcs_parser.augment_position(chunk)
        self._logger.debug(f'End update_function_position for {self._storage_name.obs_id}')


class TaylorMapping(InputTileMapping):
    def __init__(self, storage_name, headers, clients, observable, observation, config):
        super().__init__(storage_name, headers, clients, observable, observation, config)

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        self._logger.debug('Begin accumulate_bp.')
        super().accumulate_blueprint(bp)
        bp.set('Plane.provenance.name', '_get_plane_provenance_name()')
        bp.set('Plane.provenance.version', '_get_plane_provenance_version()')

        bp.configure_time_axis(5)
        bp.set('Chunk.time.axis.axis.ctype', 'TIME')
        bp.set('Chunk.time.axis.axis.cunit', 'd')
        bp.set('Chunk.time.axis.function.naxis', 1)
        bp.set('Chunk.time.axis.function.refCoord.pix', 0.5)
        bp.set('Chunk.time.axis.function.refCoord.val', '_get_time_function_refcoord_val()')

        self._logger.debug('Done accumulate_bp.')

    def _get_plane_provenance_name(self, ext):
        origin = self._headers[ext].get('ORIGIN')
        result = None
        if origin:
            result = origin
            bits = origin.split(' ')
            if len(bits) >= 3:
                other_bits = bits[2].split(':')
                result = f'{bits[0]} {other_bits[0]}'
        return result

    def _get_plane_provenance_version(self, ext):
        origin = self._headers[ext].get('ORIGIN')
        result = None
        if origin:
            bits = origin.split(' ')
            if len(bits) >= 3:
                other_bits = bits[2].split(':')
                result = f'{bits[1]} {other_bits[1]}'
        return result

    def _get_time_function_refcoord_val(self, ext):
        date_obs = self._headers[ext].get('DATE-OBS')
        if date_obs is not None:
            result = get_datetime_mjd(date_obs)
        return result

    def _update_artifact(self, artifact):
        super()._update_artifact(artifact)
        for part in artifact.parts.values():
            for chunk in part.chunks:
                if chunk.time_axis is not None:
                    chunk.time_axis = None


class OutputSpatialTemporal(Possum1DMapping):
    def __init__(self, storage_name, headers, clients, observable, observation, config):
        super().__init__(storage_name, headers, clients, observable, observation, config)

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        self._logger.debug('Begin accumulate_bp.')
        super().accumulate_blueprint(bp)

        bp.set('Plane.provenance.name', 'POSSUM')
        bp.clear('Plane.provenance.lastExecuted')
        bp.add_attribute('Plane.provenance.lastExecuted', 'DATE')

        bp.configure_position_axes((1, 2))
        bp.clear('Chunk.position.axis.function.cd11')
        bp.clear('Chunk.position.axis.function.cd22')
        bp.add_attribute('Chunk.position.axis.function.cd11', 'CDELT1')
        bp.set('Chunk.position.axis.function.cd12', 0.0)
        bp.set('Chunk.position.axis.function.cd21', 0.0)
        bp.add_attribute('Chunk.position.axis.function.cd22', 'CDELT2')
        bp.set('Chunk.position.resolution', '_get_position_resolution()')

        bp.configure_time_axis(5)
        bp.set('Chunk.time.axis.axis.ctype', 'TIME')
        bp.set('Chunk.time.axis.axis.cunit', 'd')
        bp.set('Chunk.time.axis.function.naxis', 1)
        bp.set('Chunk.time.axis.function.refCoord.pix', 0.5)

        self._logger.debug('Done accumulate_bp.')

    def update(self, file_info):
        """Called to fill multiple CAOM model elements and/or attributes
        (an n:n relationship between TDM attributes and CAOM attributes).
        """
        super().update(file_info)
        for plane in self._observation.planes.values():
            for artifact in plane.artifacts.values():
                if artifact.uri != self._storage_name.file_uri:
                    continue
                for part in artifact.parts.values():
                    for chunk in part.chunks:
                        if chunk.time is not None:
                            chunk.time_axis = None
                            if (
                                chunk.time.axis is not None
                                and chunk.time.axis.function is not None
                                and chunk.time.axis.function.ref_coord is not None
                            ):
                                # because the CD matrix is present, but not correctly for TIME (index 5)
                                chunk.time.axis.function.ref_coord.val = make_datetime(
                                    self._headers[0].get('DATE-OBS')
                                ).timestamp()
                        if chunk.energy is not None:
                            chunk.energy_axis = None
        return self._observation


class Output3DMapping(OutputSpatialTemporal):
    def __init__(self, storage_name, headers, clients, observable, observation, config):
        super().__init__(storage_name, headers, clients, observable, observation, config)

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)
        bp.configure_polarization_axis(3)
        bp.configure_custom_axis(4)
        self._logger.debug('Done accumulate_bp.')

class OutputFWHM(OutputSpatialTemporal):
    def __init__(self, storage_name, headers, clients, observable, observation, config):
        super().__init__(storage_name, headers, clients, observable, observation, config)

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        super().accumulate_blueprint(bp)
        bp.configure_custom_axis(3)
        bp.configure_polarization_axis(4)
        self._logger.debug('Done accumulate_bp.')


def mapping_factory(storage_name, headers, clients, observable, observation, config):
    if storage_name.product_id == '1d_pipeline':
        result = Possum1DMapping(storage_name, headers, clients, observable, observation, config)
    elif storage_name.product_id == '3d_pipeline':
        if '_FWHM' in storage_name.file_name:
            result = OutputFWHM(storage_name, headers, clients, observable, observation, config)
        else:
            result = Output3DMapping(storage_name, headers, clients, observable, observation, config)
    elif storage_name.product_id.startswith('multifrequencysynthesis_'):
        result = TaylorMapping(storage_name, headers, clients, observable, observation, config)
    else:
        result = InputTileMapping(storage_name, headers, clients, observable, observation, config)
    logging.error(f'Constructed {result.__class__.__name__} for mapping.')
    return result
