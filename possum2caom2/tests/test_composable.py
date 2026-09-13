# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2026.                            (c) 2026.
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

import os

from datetime import datetime
from mock import patch, PropertyMock

from caom2pipe.astro_composable import make_headers_from_file
from caom2pipe.data_source_composable import StateRunnerMeta
from caom2pipe.manage_composable import Config, State, StorageName, TaskType
from possum2caom2.storage_name import PossumName
from possum2caom2 import composable


@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
@patch(
    'caom2pipe.data_source_composable.ListDirTimeBoxDataSource.' 'get_time_box_work',
    autospec=True,
)
@patch(
    'caom2pipe.data_source_composable.ListDirTimeBoxDataSource.end_dt',
    new_callable=PropertyMock(return_value=datetime(year=2019, month=3, day=7, hour=19, minute=5)),
)
@patch('caom2pipe.execute_composable.OrganizeExecutes.do_one')
@patch('possum2caom2.possum_execute.RCloneClients')
def test_run_by_state(
    clients_mock, do_one_mock, end_time_mock, get_work_mock, access_mock, test_config, tmp_path, change_test_dir
):
    access_mock.return_value = 'https://localhost'
    test_config.interval = 3600
    test_config.logging_level = 'DEBUG'
    test_config.change_working_directory(tmp_path.as_posix())
    test_config.proxy_file_name = 'test_proxy.pem'
    test_config.task_types = [TaskType.INGEST]
    start_time = datetime(year=2019, month=3, day=3, hour=19, minute=5)
    State.write_bookmark(test_config.state_fqn, test_config.bookmark, start_time)
    Config.write_to_file(test_config)
    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')

    test_f_name = 'PSM_pilot1_944MHz_18asec_2226-5552_11268_i_v1.fits'
    test_obs_id = '944MHz_18asec_2226-5552_11268_pilot1_v1'
    do_one_mock.return_value = (0, None)
    get_work_mock.side_effect = lambda x, y, z: [
        StateRunnerMeta(
            os.path.join(os.path.join(tmp_path, 'test_files'), test_f_name),
            datetime(year=2019, month=10, day=23, hour=16, minute=27, second=19),
        ),
    ]

    # execution
    test_result = composable._run_incremental()
    assert test_result == 0, 'mocking correct execution'
    assert do_one_mock.called, 'should have been called'
    args, kwargs = do_one_mock.call_args
    test_storage = args[0]
    assert isinstance(test_storage, PossumName), type(test_storage)
    assert test_storage.obs_id == test_obs_id, f'wrong obs id {test_storage.obs_id}'
    assert test_storage.file_name == test_f_name, 'wrong file name'
    assert test_storage.file_uri == f'{test_config.scheme}:{test_config.collection}/{test_f_name}', 'wrong uri'


@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
@patch('caom2pipe.execute_composable.OrganizeExecutes.do_one')
def test_run(run_mock, access_mock, test_config, tmp_path):
    run_mock.return_value = (0, None)
    access_mock.return_value = 'https://localhost'
    test_f_id = 'PSM_pilot1_944MHz_18asec_2226-5552_11268_p3d_v1_snrPIfit'
    test_f_name = f'{test_f_id}.fits'
    orig_cwd = os.getcwd()
    try:
        os.chdir(tmp_path.as_posix())
        test_config.change_working_directory(tmp_path.as_posix())
        test_config.proxy_file_name = 'test_proxy.fqn'
        test_config.task_types = [TaskType.INGEST]
        test_config.write_to_file(test_config)

        with open(test_config.proxy_fqn, 'w') as f:
            f.write('test content')
        with open(test_config.work_fqn, 'w') as f:
            f.write(test_f_name)

        try:
            # execution
            test_result = composable._run()
        except Exception as e:
            assert False, e

        assert test_result == 0, 'wrong return value'
        assert run_mock.called, 'should have been called'
        args, kwargs = run_mock.call_args
        test_storage = args[0]
        assert isinstance(test_storage, StorageName), type(test_storage)
        assert test_storage.file_name == test_f_name, 'wrong file name'
        assert test_storage.source_names[0] == test_f_name, 'wrong fname on disk'
    finally:
        os.chdir(orig_cwd)


@patch('possum2caom2.possum_execute.get_local_headers_from_fits')
@patch('possum2caom2.possum_execute.RCloneTransfer.get')
@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
@patch('caom2pipe.execute_composable.OrganizeExecutes.do_one')
def test_run_remote_todo(run_mock, access_mock, get_mock, headers_mock, test_config, tmp_path, change_test_dir):
    run_mock.return_value = (0, None)
    access_mock.return_value = 'https://localhost'
    test_pawsey_path = 'components/0049-51/survey/i/PSM.0049-51.11400.q.fits'
    test_f_id = 'PSM_943MHz_20asec_1026-5855_11400_q_v1'
    test_f_name = f'{test_f_id}.fits'
    test_config.change_working_directory(tmp_path.as_posix())
    test_config.proxy_file_name = 'test_proxy.fqn'
    test_config.task_types = [TaskType.STORE, TaskType.INGEST, TaskType.MODIFY]
    test_config.data_sources = ['pawsey_test/acacia/possum1234']
    test_config.data_source_extensions = ['.fits', '.fits.header']
    test_config.logging_level = 'DEBUG'
    test_config.write_to_file(test_config)

    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')
    with open(test_config.work_fqn, 'w') as f:
        f.write(test_pawsey_path)

    get_mock.side_effect = _copy_pawsey_file
    headers_mock.side_effect = _pawsey_headers_mock

    try:
        # execution
        test_result = composable._run_remote_todo()
    except Exception as e:
        assert False, e

    assert test_result == 0, 'wrong return value'
    assert run_mock.called, 'should have been called'
    args, _ = run_mock.call_args
    test_storage = args[0]
    assert isinstance(test_storage, StorageName), type(test_storage)
    assert test_storage.file_name == test_f_name, 'wrong file name'
    assert os.path.basename(test_storage.source_names[0]) == test_f_name, \
        f'wrong fname on disk got {test_storage.source_names[0]} expected {test_f_name}'
    assert get_mock.call_count == 1, 'number of rclone copy calls'
    get_mock.assert_called_with(
        test_pawsey_path, f'{tmp_path.as_posix()}/{os.path.basename(test_pawsey_path)}'
    ), 'wrong rclone call'


def _copy_pawsey_file(ignore_source, destination_fqn):
    with open(destination_fqn, 'w') as f:
        f.write('test content')


def _pawsey_headers_mock(_):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    return make_headers_from_file(
        os.path.join(os.path.join(this_dir, 'data/possum'), 'PSM_943MHz_20asec_1026-5855_11400_q_v1.fits.header'),
    )

