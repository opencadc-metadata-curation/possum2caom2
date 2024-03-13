# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2019.                            (c) 2019.
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

from glob import glob
from os.path import basename

from possum2caom2 import PossumName


def pytest_generate_tests(metafunc):
    test_data_dir = f'{metafunc.config.invocation_dir}/data'
    obs_id_list = glob(f'{test_data_dir}/**/*.fits.header')
    metafunc.parametrize('test_name', obs_id_list)


def test_storage_name(test_config, test_name):
    test_obs_ids = [
        '1368MHz_18asec_2031-5249_11073_pilot1',
        '1367MHz_18asec_2013-5553_11261_pilot1',
        '944MHz_18asec_2226-5552_11268_pilot1',
        '1367MHz_18asec_2039-5115_10973_pilot1',
        '1368MHz_18asec_2013-5552_11261_pilot1',
    ]
    test_f_name = basename(test_name)
    test_uri = f'{test_config.scheme}:{test_config.collection}/{test_f_name.replace(".header", "")}'
    for entry in [
        test_f_name,
        test_uri,
        f'https://localhost:8020/{test_f_name}',
        f'vos:goliaths/test/{test_f_name}',
        f'/tmp/{test_f_name}',
    ]:
        test_subject = PossumName(entry)
        assert test_subject.obs_id in test_obs_ids, f'wrong obs id {test_f_name} {test_subject}'
        assert test_subject.source_names == [entry], f'wrong source names {test_f_name}'
        if 'p3d' in entry:
            assert test_subject.product_id == '3d_pipeline', f'wrong product id {test_subject.product_id}'
        elif 'p1d' in entry:
            assert test_subject.product_id == '1d_pipeline', f'wrong product id {test_subject.product_id}'
        else:
            if '_t0_' in entry:
                assert (
                    test_subject.product_id == 'multifrequencysynthesis_i_t0'
                ), f'wrong product id {test_subject.product_id}'
            elif '_t1_' in entry:
                assert (
                    test_subject.product_id == 'multifrequencysynthesis_i_t1'
                ), f'wrong product id {test_subject.product_id}'
            elif '_i_' in entry:
                assert test_subject.product_id == 'raw_i', f'wrong product id {test_subject.product_id}'
            else:
                assert test_subject.product_id == 'raw_qu', f'wrong product id {test_subject.product_id}'
        test_check_uri = f'{test_config.scheme}:{test_config.collection}/{test_f_name.replace(".header", "")}'
        assert test_subject.file_uri == test_check_uri, f'wrong file uri {test_subject}'
        assert test_subject.destination_uris == [test_check_uri], f'wrong uris {test_subject}'
        assert test_subject.prev == f'{test_subject.obs_id}_{test_subject.product_id}_prev.jpg', 'preview uri'
        assert test_subject.thumb == f'{test_subject.obs_id}_{test_subject.product_id}_prev_256.jpg', 'thumbnail uri'
        assert (
            test_subject.prev_uri == f'{test_config.preview_scheme}:{test_config.collection}/{test_subject.prev}'
        ), 'preview uri'
        assert (
            test_subject.thumb_uri == f'{test_config.preview_scheme}:{test_config.collection}/{test_subject.thumb}'
        ), 'thumbnail uri'
