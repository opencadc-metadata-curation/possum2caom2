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

import logging

from datetime import datetime
from os import listdir
from os.path import exists
from shutil import copyfile

from caom2utils.data_util import get_local_file_headers
from possum2caom2 import possum_execute
from caom2pipe.manage_composable import ExecutionReporter, make_datetime, Observable, State, TaskType
from mock import ANY, call, Mock, patch, PropertyMock


# def test_choose_store(test_config):
#     test_config.data_sources = ['rclone_remote']
#     test_config.task_types = [TaskType.STORE]
#     test_subject = RCloneOrganizeExecutes(
#         test_config,
#         meta_visitors=[],
#         data_visitors=[],
#         chooser=None,
#         store_transfer=Mock(),
#         modify_transfer=None,
#         metadata_reader=None,
#         clients=Mock(),
#         observable=Mock(),
#         reporter=Mock(),
#     )
#     test_subject.choose()
#     assert len(test_subject._executors) == 1, 'length'
#     assert isinstance(test_subject._executors[0], RCloneStore), 'type'


# def test_choose_visit(test_config):
#     test_config.data_sources = ['rclone_remote']
#     test_config.task_types = [TaskType.STORE, TaskType.INGEST]
#     test_subject = RCloneOrganizeExecutes(
#         test_config,
#         meta_visitors=[],
#         data_visitors=[],
#         chooser=None,
#         store_transfer=Mock(),
#         modify_transfer=None,
#         metadata_reader=None,
#         clients=Mock(),
#         observable=Mock(),
#         reporter=Mock(),
#     )
#     test_subject.choose()
#     assert len(test_subject._executors) == 1, 'length'
#     assert isinstance(test_subject._executors[0], RCloneNoFheadStoreVisit), 'type'


# def test_rclone_store(test_config, tmp_path, change_test_dir):
#     test_config.change_working_directory(tmp_path)
#     test_config.data_sources = ['rclone_test']
#     working_path = tmp_path.joinpath(test_config.data_sources[0])
#     working_path.mkdir()
#     test_observable = Mock()
#     test_clients = Mock()
#     test_metadata_reader = Mock()
#     test_transferrer = Mock()

#     test_subject = RCloneStore(test_config, test_observable, test_clients, test_metadata_reader, test_transferrer)
#     test_subject.execute(None)
#     assert test_observable.mock_calls == [], 'observable'
#     assert test_clients.mock_calls == [], 'clients'
#     assert test_metadata_reader.mock_calls == [], 'metadata_reader'
#     assert test_transferrer.mock_calls == [call.get(None, f'{tmp_path}/{test_config.data_sources[0]}')], 'transferrer'


# def test_rclone_visit(test_config):
#     test_observable = Mock()
#     test_clients = Mock()
#     test_metadata_reader = Mock()
#     test_transferrer = Mock()
#     test_meta_visitors = []
#     test_data_visitors = []
#     test_subject = RCloneNoFheadStoreVisit(
#         test_config,
#         test_clients,
#         test_transferrer,
#         test_meta_visitors,
#         test_data_visitors,
#         test_metadata_reader,
#         test_observable,
#     )
#     test_subject.execute(None)
#     assert test_observable.mock_calls == [], 'observable'
#     assert test_clients.mock_calls == [], 'clients'
#     assert test_metadata_reader.mock_calls == [], 'metadata_reader'
#     assert test_transferrer.mock_calls == [], 'transferrer'


def test_execution_unit_start_stop(test_config, tmp_path):
    kwargs = {
        'clients': None,
        'data_source': None,
        'metadata_reader': None,
        'observable': None,
        'reporter': None,
        'prev_exec_dt': make_datetime('2023-10-28T20:47:49.000000000Z'),
        'exec_dt': make_datetime('2023-11-28T20:47:49.000000000Z'),
    }
    test_config.change_working_directory(tmp_path)
    test_config.cleanup_files_when_storing = True
    test_subject = possum_execute.ExecutionUnit(test_config, **kwargs)

    # preconditions
    test_files = listdir(tmp_path)
    assert len(test_files) == 0, 'directory should be empty'

    test_subject.start()

    test_files = listdir(tmp_path)
    assert len(test_files) == 1, 'directory should have a workspace directory'
    assert '2023-10-28T20_47_49_2023-11-28T20_47_49' in test_files, 'wrong working directory'

    test_subject.stop()

    test_files = listdir(tmp_path)
    assert len(test_files) == 0, 'post-condition directory should be cleaned up and empty'

    test_config.cleanup_files_when_storing = False
    test_subject = possum_execute.ExecutionUnit(test_config, **kwargs)

    # preconditions
    test_files = listdir(tmp_path)
    assert len(test_files) == 0, 'directory should be empty'

    test_subject.start()

    test_files = listdir(tmp_path)
    assert len(test_files) == 1, 'directory should have a workspace directory'
    assert '2023-10-28T20_47_49_2023-11-28T20_47_49' in test_files, 'wrong working directory'

    test_subject.stop()

    test_files = listdir(tmp_path)
    assert len(test_files) == 1, 'post-condition directory should exist'


# need test_config parameter so StorageName.collection is set
def test_remote_metadata_reader_file_info(test_config, test_data_dir):
    input_file = f'{test_data_dir}/storage_mock/rclone_lsjson.json'
    test_file_uri = 'cadc:POSSUM/PSM.0049-51.10887.i.fits'
    test_subject = possum_execute.RemoteMetadataReader()

    with open(input_file) as f:
        test_subject.set_file_info(f.read())

    assert len(test_subject.file_info) == 4, 'wrong number of results'
    logging.error(test_subject.file_info.keys())
    test_result = test_subject.file_info.get(test_file_uri)
    assert test_result is not None, 'expect a result'
    assert test_result.size == 4831848000, 'wrong size'
    assert test_result.file_type == 'application/fits', 'wrong file type'
    assert test_result.lastmod == datetime(2023, 11, 18, 20, 47, 50), 'wrong modification time'

    # TODO - are these StorageName instances usable?
    test_storage_name = test_subject.storage_names.get(test_file_uri)
    assert test_storage_name.file_uri == test_file_uri, 'wrong file uri'


@patch('possum2caom2.possum_execute.exec_cmd')
@patch('possum2caom2.possum_execute.exec_cmd_info')
def test_remote_data_source(exec_cmd_info_mock, exec_cmd_mock, test_data_dir, test_config, tmp_path):
    with open(f'{test_data_dir}/storage_mock/rclone_lsjson.json') as f:
        exec_cmd_info_mock.return_value = f.read()

    exec_cmd_mock.side_effect = Mock()

    test_config.change_working_directory(tmp_path)
    test_start_key = 'test/acacia_possum/pawsey0980'
    test_start_time = make_datetime('2023-10-28T20:47:49.000000000Z')
    test_end_time = make_datetime('2023-11-28T20:47:49.000000000Z')
    State.write_bookmark(test_config.state_fqn, test_start_key, test_start_time)
    test_metadata_reader = possum_execute.RemoteMetadataReader()
    mock_1 = Mock()
    mock_2 = Mock()
    mock_3 = Mock()
    kwargs = {
        'clients': mock_1,
        'observable': mock_2,
        'reporter': mock_3,
    }
    test_subject = possum_execute.RemoteIncrementalDataSource(
        test_config,
        test_start_key,
        test_metadata_reader,
        **kwargs,
    )
    test_subject.initialize_start_dt()
    assert test_subject.start_dt == datetime(2023, 10, 28, 20, 47, 49), 'start_dt'
    test_subject.initialize_end_dt()
    assert test_subject.end_dt == datetime(2023, 11, 18, 20, 47, 50), 'end_dt'
    test_result = test_subject.get_time_box_work(test_start_time, test_end_time)
    assert test_result is not None, 'expect a result'
    assert test_result._clients == mock_1, 'clients'
    assert test_result._observable == mock_2, 'observable'
    assert test_result._reporter == mock_3, 'reporter'
    assert test_result._metadata_reader == test_metadata_reader, 'reader'


def test_state_runner_reporter(test_config, tmp_path, change_test_dir):
    # make sure the StateRunner goes through at least one time box check, and creates the
    # right log locations
    test_config.change_working_directory(tmp_path)
    test_config.task_types = [TaskType.STORE, TaskType.INGEST, TaskType.MODIFY]
    test_config.data_sources = ['test/acacia:possum1234']
    test_config.interval = 60 * 48  # work in time-boxes of 2 days => 60m * 48h
    test_organizer = Mock()
    # the time-box times, or, this is "when" the code looks
    test_start_time = make_datetime('2023-10-28T20:47:49.000000000Z')
    test_end_time = make_datetime('2023-11-28T20:47:49.000000000Z')
    test_data_source = Mock()
    end_time_mock = PropertyMock(return_value=test_end_time)
    start_time_mock = PropertyMock(return_value=test_start_time)
    type(test_data_source).end_dt = end_time_mock
    type(test_data_source).start_dt = start_time_mock
    # the execution times, or, this is "what" the code finds
    test_entry_time = make_datetime('2023-11-28T08:47:49.000000000Z')
    execution_unit_mock = Mock()
    type(execution_unit_mock).entry_dt = test_entry_time
    type(execution_unit_mock).num_entries = 1
    test_data_source.get_time_box_work.return_value = execution_unit_mock
    test_data_sources = [test_data_source]
    test_observable = Mock()
    test_reporter = ExecutionReporter(test_config, test_observable)
    test_subject = possum_execute.ExecutionUnitStateRunner(
        test_config,
        test_organizer,
        test_data_sources,
        test_observable,
        test_reporter,
    )
    test_result = test_subject.run()
    assert test_result is not None, 'expect a result'
    assert test_result == 0, 'happy path'
    assert test_organizer.mock_calls == [call.choose], 'organizer'
    assert test_data_source.initialize_start_dt.called, 'initialize_start_dt'
    assert test_data_source.initialize_start_dt.call_count == 1, 'initialize_start_dt count'
    assert test_data_source.initialize_end_dt.called, 'initialize_end_dt'
    assert test_data_source.initialize_end_dt.call_count == 1, 'initialize_end_dt count'
    assert test_data_source.get_time_box_work.called, 'get_time_box_work'
    # 16 == number of days / 2 between the start and end times
    assert test_data_source.get_time_box_work.call_count == 16, 'get_time_box_work count'
    test_observable.assert_has_calls([]), 'observable calls'
    assert exists(test_config.failure_fqn), 'failure'
    assert exists(test_config.progress_fqn), 'progress'
    assert exists(test_config.success_fqn), 'success'
    assert exists(test_config.retry_fqn), 'retries'
    assert exists(test_config.total_retry_fqn), 'total_retries'
    assert exists(test_config.report_fqn), 'report'
    assert test_reporter.success == 0, 'reporter success'
    assert test_reporter.all == 0, 'reporter all'


@patch('caom2utils.data_util.get_local_headers_from_fits')
@patch('possum2caom2.possum_execute.exec_cmd')
@patch('possum2caom2.possum_execute.exec_cmd_info')
def test_state_runner_nominal_multiple_files(
    exec_cmd_info_mock, exec_cmd_mock, header_mock, test_config, test_data_dir, tmp_path, change_test_dir
):
    # test that three file get processed properly, and get left behind
    with open(f'{test_data_dir}/storage_mock/rclone_lsjson.json') as f:
        exec_cmd_info_mock.return_value = f.read()

    w_test_file = 'PSM.0049-51.11092.w.fits.header'
    q_test_file = 'PSM.1136-64.11836.q.fits.header'
    u_test_file = 'PSM.1136-64.11485.u.fits.header'
    time_box_dir_name = '2023-10-28T20_47_49_2023-10-30T20_47_49'
    def _exec_cmd_mock(arg1):
        if arg1 == (
            'rclone copy --min 2023-10-28T20:47:49 --max 2023-10-30T20:47:49 --includes *.fits *.fits.header'
        ):
            copyfile(
                f'{test_data_dir}/casda/PSM_pilot1_1367MHz_18asec_2013-5553_11261_t0_i_v1.fits.header',
                f'{tmp_path}/{time_box_dir_name}/{w_test_file}',
            )
            copyfile(
                f'{test_data_dir}/casda/PSM_pilot1_1368MHz_18asec_2031-5249_11073_i_v1.fits.header',
                f'{tmp_path}/{time_box_dir_name}/{q_test_file}',
            )
            copyfile(
                f'{test_data_dir}/casda/PSM_pilot1_1368MHz_18asec_2031-5249_11073_q_v1.fits.header',
                f'{tmp_path}/{time_box_dir_name}/{u_test_file}',
            )
    exec_cmd_mock.side_effect = _exec_cmd_mock
    header_mock.side_effect = get_local_file_headers
    test_config.change_working_directory(tmp_path)
    test_config.cleanup_files_when_storing = False
    test_config.task_types = [TaskType.STORE, TaskType.INGEST, TaskType.MODIFY]
    test_config.data_sources = ['test/acacia/possum1234']
    test_config.data_source_extensions = ['.fits', '.fits.header']
    test_config.logging_level = 'INFO'
    test_config.interval = 60 * 48  # work in time-boxes of 2 days => 60m * 48h
    test_config.observe_execution = True
    test_organizer = Mock()
    # the time-box times, or, this is "when" the code looks
    test_start_time = make_datetime('2023-10-28T20:47:49.000000000Z')
    test_end_time = make_datetime('2023-11-18T20:47:50.000000000Z')
    State.write_bookmark(test_config.state_fqn, test_config.data_sources[0], test_start_time)
    test_metadata_reader = possum_execute.RemoteMetadataReader()
    test_observable = Observable(test_config)
    test_reporter = ExecutionReporter(test_config, test_observable)
    test_clients = Mock()
    test_clients.metadata_client.read.return_value = None
    kwargs = {
        'clients': test_clients,
        'observable': test_observable,
        'reporter': test_reporter
    }
    test_data_source = possum_execute.RemoteIncrementalDataSource(
        test_config,
        test_config.data_sources[0],
        test_metadata_reader,
        **kwargs,
    )
    test_data_source.reporter = test_reporter
    test_data_sources = [test_data_source]
    test_subject = possum_execute.ExecutionUnitStateRunner(
        test_config,
        test_organizer,
        test_data_sources,
        test_observable,
        test_reporter,
    )
    test_result = test_subject.run()
    assert test_result is not None, 'expect a result'
    assert test_result == 0, 'happy path'
    assert test_organizer.mock_calls == [call.choose], 'organizer'
    assert test_clients.mock_calls == [
        call.data_client.put(f'{tmp_path}/{time_box_dir_name}', f'cadc:POSSUM/{u_test_file.replace(".header", "")}'),
        call.metadata_client.read('POSSUM', '1136-64_11485'),
        call.metadata_client.create(ANY),
        call.data_client.put(f'{tmp_path}/{time_box_dir_name}', f'cadc:POSSUM/{q_test_file.replace(".header", "")}'),
        call.metadata_client.read('POSSUM', '1136-64_11836'),
        call.metadata_client.create(ANY),
        call.data_client.put(f'{tmp_path}/{time_box_dir_name}', f'cadc:POSSUM/{w_test_file.replace(".header", "")}'),
        call.metadata_client.read('POSSUM', '0049-51_11092'),
        call.metadata_client.create(ANY),
    ], f'clients {test_clients.mock_calls}'
    assert exists(test_config.rejected_fqn), f'rejected {test_config.rejected_fqn}'
    assert exists(test_config.observable_directory), f'metrics {test_config.observable_directory}'
    assert test_data_source.end_dt == test_end_time, 'end_dt'
    assert test_clients.data_client.put.called, 'data client put'
    assert test_reporter.all == 3, 'wrong file count'
    assert test_reporter.success == 3, 'wrong file count'
    left_behind = listdir(f'{tmp_path}/{time_box_dir_name}')
    assert len(left_behind) == 3, 'no files cleaned up'


def test_state_runner_clean_up_when_storing(test_config):
    # one file gets cleaned up
    test_config.cleanup_files_when_storing = True
    with open(f'{test_data_dir}/storage_mock/rclone_lsjson.json') as f:
        exec_cmd_info_mock.return_value = f.read()

    w_test_file = 'PSM.0049-51.11092.w.fits.header'
    q_test_file = 'PSM.1136-64.11836.q.fits.header'
    u_test_file = 'PSM.1136-64.11485.u.fits.header'
    time_box_dir_name = '2023-10-28T20_47_49_2023-10-30T20_47_49'
    def _exec_cmd_mock(arg1):
        if arg1 == (
            'rclone copy --min 2023-10-28T20:47:49 --max 2023-10-30T20:47:49 --includes *.fits *.fits.header'
        ):
            copyfile(
                f'{test_data_dir}/casda/PSM_pilot1_1367MHz_18asec_2013-5553_11261_t0_i_v1.fits.header',
                f'{tmp_path}/{time_box_dir_name}/{w_test_file}',
            )
            copyfile(
                f'{test_data_dir}/casda/PSM_pilot1_1368MHz_18asec_2031-5249_11073_i_v1.fits.header',
                f'{tmp_path}/{time_box_dir_name}/{q_test_file}',
            )
            copyfile(
                f'{test_data_dir}/casda/PSM_pilot1_1368MHz_18asec_2031-5249_11073_q_v1.fits.header',
                f'{tmp_path}/{time_box_dir_name}/{u_test_file}',
            )
    exec_cmd_mock.side_effect = _exec_cmd_mock
    header_mock.side_effect = get_local_file_headers
    test_config.change_working_directory(tmp_path)
    test_config.cleanup_files_when_storing = False
    test_config.task_types = [TaskType.STORE, TaskType.INGEST, TaskType.MODIFY]
    test_config.data_sources = ['test/acacia/possum1234']
    test_config.data_source_extensions = ['.fits', '.fits.header']
    test_config.logging_level = 'INFO'
    test_config.interval = 60 * 48  # work in time-boxes of 2 days => 60m * 48h
    test_config.observe_execution = True
    test_organizer = Mock()
    # the time-box times, or, this is "when" the code looks
    test_start_time = make_datetime('2023-10-28T20:47:49.000000000Z')
    test_end_time = make_datetime('2023-11-18T20:47:50.000000000Z')
    State.write_bookmark(test_config.state_fqn, test_config.data_sources[0], test_start_time)
    test_metadata_reader = possum_execute.RemoteMetadataReader()
    test_observable = Observable(test_config)
    test_reporter = ExecutionReporter(test_config, test_observable)
    test_clients = Mock()
    test_clients.metadata_client.read.return_value = None
    kwargs = {
        'clients': test_clients,
        'observable': test_observable,
        'reporter': test_reporter
    }
    test_data_source = possum_execute.RemoteIncrementalDataSource(
        test_config,
        test_config.data_sources[0],
        test_metadata_reader,
        **kwargs,
    )
    test_data_source.reporter = test_reporter
    test_data_sources = [test_data_source]
    test_subject = possum_execute.ExecutionUnitStateRunner(
        test_config,
        test_organizer,
        test_data_sources,
        test_observable,
        test_reporter,
    )
    test_result = test_subject.run()
    assert test_result is not None, 'expect a result'
    assert test_result == 0, 'happy path'
    assert test_organizer.mock_calls == [call.choose], 'organizer'
    assert test_clients.mock_calls == [
        call.data_client.put(f'{tmp_path}/{time_box_dir_name}', f'cadc:POSSUM/{u_test_file.replace(".header", "")}'),
        call.metadata_client.read('POSSUM', '1136-64_11485'),
        call.metadata_client.create(ANY),
        call.data_client.put(f'{tmp_path}/{time_box_dir_name}', f'cadc:POSSUM/{q_test_file.replace(".header", "")}'),
        call.metadata_client.read('POSSUM', '1136-64_11836'),
        call.metadata_client.create(ANY),
        call.data_client.put(f'{tmp_path}/{time_box_dir_name}', f'cadc:POSSUM/{w_test_file.replace(".header", "")}'),
        call.metadata_client.read('POSSUM', '0049-51_11092'),
        call.metadata_client.create(ANY),
    ], f'clients {test_clients.mock_calls}'
    assert exists(test_config.rejected_fqn), f'rejected {test_config.rejected_fqn}'
    assert exists(test_config.observable_directory), f'metrics {test_config.observable_directory}'
    assert test_data_source.end_dt == test_end_time, 'end_dt'
    assert test_clients.data_client.put.called, 'data client put'
    assert test_reporter.all == 3, 'wrong file count'
    assert test_reporter.success == 3, 'wrong file count'
    left_behind = listdir(f'{tmp_path}/{time_box_dir_name}')
    assert len(left_behind) == 3, 'no files cleaned up'
    assert False
