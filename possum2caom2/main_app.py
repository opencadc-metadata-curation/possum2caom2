# -*- coding: utf-8 -*-
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
This module implements the ObsBlueprint mapping, as well as the workflow 
entry point that executes the workflow.
"""

import logging
import traceback

from math import sqrt
from os.path import basename

from caom2 import CalibrationLevel, DataProductType, ProductType, ReleaseType
from caom2pipe import caom_composable as cc
from caom2pipe.manage_composable import CadcException, make_time, StorageName, ValueRepairCache


__all__ = ['PossumName', 'APPLICATION', 'mapping_factory']


APPLICATION = 'possum2caom2'


class PossumName(StorageName):
    """
    Inputs:
    PSM_944MHz_18asec_2226-5552_11268_i_v1.fits
    PSM_944MHz_18asec_2226-5552_11268_q_v1.fits
    PSM_944MHz_18asec_2226-5552_11268_u_v1.fits

    Outputs:
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
    def is_1d_output(self):
        return 'pilot' in self._file_id and 'RMSF1D' in self._file_id

    @property
    def is_output(self):
        return 'pilot' in self._file_id

    def is_valid(self):
        return True

    def set_obs_id(self):
        # picking the common prefix to start with, e.g. 944MHz_18asec_2226-5552_11268
        # leave off the "PSM" and "PSM_pilot1" because it's inconsistent between the inputs and the outputs
        bits = self._file_id.split('_')

        if self.is_output:
            self._obs_id = f'{bits[2]}_{bits[3]}_{bits[4]}_{bits[5]}'
        else:
            self._obs_id = f'{bits[1]}_{bits[2]}_{bits[3]}_{bits[4]}'

    def set_product_id(self):
        self._product_id = self._file_id.split(self._obs_id)[-1].lstrip('_')
        if '_p3d_' in self._file_id:
            self._product_id = '3d_pipeline'
        elif '_p1d_' in self._file_id:
            self._product_id = '1d_pipeline'


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

    def __init__(self, storage_name, headers, clients):
        super().__init__(storage_name, headers, clients)

    def accumulate_blueprint(self, bp, application=None):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        self._logger.debug('Begin accumulate_bp.')
        super().accumulate_blueprint(bp, APPLICATION)
        bp.set('Observation.metaRelease', '2025-01-01T00:00:00.000')
        bp.set('Plane.calibrationLevel', CalibrationLevel.CALIBRATED)
        bp.set('Plane.dataProductType', '_get_data_product_type()')
        bp.set('Plane.metaRelease', '2025-01-01T00:00:00.000')
        bp.set('Plane.dataRelease', '2025-01-01T00:00:00.000')
        bp.set('Artifact.productType', ProductType.SCIENCE)
        bp.set('Artifact.releaseType', ReleaseType.DATA)
        self._logger.debug('Done accumulate_bp.')

    def update(self, observation, file_info):
        """Called to fill multiple CAOM model elements and/or attributes
        (an n:n relationship between TDM attributes and CAOM attributes).
        """
        self._logger.debug(f'Begin update for {observation.observation_id}.')
        try:
            super().update(observation, file_info)
            Possum1DMapping.value_repair.repair(observation)
            self._logger.debug('Done update.')
            return observation
        except CadcException as e:
            tb = traceback.format_exc()
            self._logger.debug(tb)
            self._logger.error(e)
            self._logger.error(
                f'Terminating ingestion for {observation.observation_id}'
            )
            return None

    def _get_data_product_type(self, ext):
        naxis = self._headers[ext].get('NAXIS')
        result = DataProductType.CUBE
        if naxis == 0:
            result = DataProductType.MEASUREMENTS
        elif naxis == 2:
            result = DataProductType.IMAGE
        return result


class PossumInputMapping(Possum1DMapping):
    def __init__(self, storage_name, headers, clients):
        super().__init__(storage_name, headers, clients)

    def accumulate_blueprint(self, bp, application=None):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        self._logger.debug('Begin accumulate_bp.')
        super().accumulate_blueprint(bp, APPLICATION)
        bp.configure_position_axes((1, 2))
        bp.configure_time_axis(3)
        bp.configure_energy_axis(4)
        bp.configure_polarization_axis(5)

        bp.set('Plane.calibrationLevel', CalibrationLevel.RAW_STANDARD)
        self._logger.debug('Done accumulate_bp.')


class PossumOutputMapping(Possum1DMapping):
    def __init__(self, storage_name, headers, clients):
        super().__init__(storage_name, headers, clients)

    def accumulate_blueprint(self, bp, application=None):
        """Configure the telescope-specific ObsBlueprint at the CAOM model
        Observation level."""
        self._logger.debug('Begin accumulate_bp.')
        super().accumulate_blueprint(bp, APPLICATION)

        bp.configure_position_axes((1, 2))
        bp.configure_polarization_axis(3)
        bp.configure_custom_axis(4)

        bp.configure_time_axis(5)

        bp.set('Plane.provenance.name', 'POSSUM')
        bp.clear('Plane.provenance.lastExecuted')
        bp.add_attribute('Plane.provenance.lastExecuted', 'DATE')

        bp.clear('Chunk.position.axis.function.cd11')
        bp.clear('Chunk.position.axis.function.cd22')
        bp.add_attribute('Chunk.position.axis.function.cd11', 'CDELT1')
        bp.set('Chunk.position.axis.function.cd12', 0.0)
        bp.set('Chunk.position.axis.function.cd21', 0.0)
        bp.add_attribute('Chunk.position.axis.function.cd22', 'CDELT2')
        bp.set('Chunk.position.resolution', '_get_position_resolution()')

        bp.set('Chunk.time.axis.axis.ctype', 'TIME')
        bp.set('Chunk.time.axis.axis.cunit', 'd')
        bp.set('Chunk.time.axis.function.naxis', 1)
        bp.set('Chunk.time.axis.function.refCoord.pix', 0.5)

        self._logger.debug('Done accumulate_bp.')

    def update(self, observation, file_info):
        """Called to fill multiple CAOM model elements and/or attributes
        (an n:n relationship between TDM attributes and CAOM attributes).
        """
        super().update(observation, file_info)
        for plane in observation.planes.values():
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
                                chunk.time.axis.function.ref_coord.val = make_time(self._headers[0].get('DATE-OBS')).timestamp()
                        if chunk.energy is not None:
                            chunk.energy_axis = None
        return observation

    def _get_position_resolution(self, ext):
        bmaj = self._headers[ext]['BMAJ']
        bmin = self._headers[ext]['BMIN']
        # From
        # https://open-confluence.nrao.edu/pages/viewpage.action?pageId=13697486
        # Clare Chandler via JJK - 21-08-18
        result = None
        if bmaj is not None and bmaj != 'INF' and bmin is not None and bmin != 'INF':
            result = 3600.0 * sqrt(bmaj * bmin)
        return result


def mapping_factory(storage_name, headers, clients):
    if storage_name.is_1d_output:
        result = Possum1DMapping(storage_name, headers, clients)
    elif storage_name.is_output:
        result = PossumOutputMapping(storage_name, headers, clients)
    else:
        result = PossumInputMapping(storage_name, headers, clients)
    return result
