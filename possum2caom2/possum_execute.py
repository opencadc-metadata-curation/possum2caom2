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

import json
import logging
import os

from cadcdata.storageinv import FileInfo
from caom2pipe.client_composable import ClientCollection
from caom2pipe.data_source_composable import IncrementalDataSource, ListDirSeparateDataSource
from caom2pipe.execute_composable import OrganizeExecutes
from caom2pipe.manage_composable import Config, create_dir, exec_cmd, exec_cmd_info, ExecutionReporter, increment_time
from caom2pipe.manage_composable import make_datetime, Observable, StorageName, TaskType
from caom2pipe.name_builder_composable import EntryBuilder
from caom2pipe.reader_composable import FileMetadataReader
from caom2pipe.run_composable import set_logging, StateRunner, TodoRunner
from caom2pipe.transfer_composable import CadcTransfer, Transfer
from caom2utils.data_util import get_file_type
from possum2caom2 import fits2caom2_augmentation, preview_augmentation
from possum2caom2.storage_name import PossumName


__all__ = ['DATA_VISITORS', 'META_VISITORS', 'remote_execution']
META_VISITORS = [fits2caom2_augmentation]
DATA_VISITORS = [preview_augmentation]

"""
TODO List:
1. config switch to keep the files, or not
2. redo a single file?
3. which files of the 8510 do you actually want? all of them copied to CANFAR? which ones archived?

"""
class RCloneClients(ClientCollection):

    def __init__(self, config):
        super().__init__(config)
        # TODO credentials
        self._rclone_client = None

    @property
    def rclone_client(self):
        return self._rclone_client


class RemoteMetadataReader(FileMetadataReader):
    def __init__(self):
        super().__init__()
        self._storage_names = {}
        self._max_dt = None

    @property
    def max_dt(self):
        return self._max_dt

    @property
    def storage_names(self):
        return self._storage_names

    def _retrieve_file_info(self, key, source_name):
        raise NotImplementedError

    def get_time_box_work_parameters(self, prev_exec_time, exec_time):
        self._logger.debug(f'Begin get_time_box_work_parameters from {prev_exec_time} to {exec_time}')
        count = 0
        max_time_box = prev_exec_time
        for entry in self._file_info.values():
            if prev_exec_time < entry.lastmod <= exec_time:
                count += 1
                max_time_box = max(prev_exec_time, entry.lastmod)
        self._logger.debug(f'End get_time_box_work_parameters with count {count}')
        return count, max_time_box

    def set(self, storage_name):
        self._logger.debug(f'Begin set for {storage_name.file_name}')
        if isinstance(storage_name, StorageName):
            self.set_headers(storage_name)
        else:
            raise NotImplementedError
        self._logger.debug('End set')

    def set_file_info(self, storage_name):
        """
        Path elements from the JSON listing:
        components
        components mfs
        components mfs i
        components mfs mfs
        components mfs w
        components survey
        components survey i
        components survey q
        components survey u
        components survey w
        :param storage_name str JSON
        """
        self._logger.debug('Begin set_file_info with rclone lsjson output')
        content = json.loads(storage_name)
        for entry in content:
            name = entry.get('Name')
            if name.startswith('PSM') and '.fits' in name:
                # keys are destination URIs
                storage_name = PossumName(name)
                if storage_name.file_uri not in self._file_info:
                    self._logger.debug(f'Retrieve FileInfo for {storage_name.file_uri}')
                    self._file_info[storage_name.file_uri] = FileInfo(
                        id=storage_name.file_uri,
                        file_type=get_file_type(name),
                        size=entry.get('Size'),
                        lastmod=make_datetime(entry.get('ModTime')),
                    )
                    if self._max_dt:
                        self._max_dt = max(self._file_info[storage_name.file_uri].lastmod, self._max_dt)
                    else:
                        self._max_dt = self._file_info[storage_name.file_uri].lastmod
                    self._storage_names[storage_name.file_uri] = storage_name
        self._logger.debug('End set_file_info')


class RemoteIncrementalDataSource(IncrementalDataSource):

    def __init__(self, config, start_key, metadata_reader, **kwargs):
        super().__init__(config, start_key)
        self._data_source_extensions = ' '.join(f'*{ii}' for ii in config.data_source_extensions)
        self._metadata_reader = metadata_reader
        self._kwargs = kwargs

    def _initialize_end_dt(self):
        self._logger.debug('Begin _initialize_end_dt')
        output = exec_cmd_info(f'rclone lsjson {self._start_key} --min {self._start_dt} --includes {self._data_source_extensions}')
        self._metadata_reader.set_file_info(output)
        self._end_dt = self._metadata_reader.max_dt
        self._logger.debug('End _initialize_end_dt')

    def get_time_box_work(self, prev_exec_dt, exec_dt):
        self._logger.debug('Begin get_time_box_work')
        self._kwargs['prev_exec_dt'] = prev_exec_dt
        self._kwargs['exec_dt'] = exec_dt
        self._kwargs['metadata_reader'] = self._metadata_reader
        execution_unit = ExecutionUnit(self._config, **self._kwargs)
        execution_unit.start()
        # get the files from the DataSource to the staging space
        exec_cmd(f'rclone copy --min {prev_exec_dt.isoformat()} --max {exec_dt.isoformat()} --includes {self._data_source_extensions}')
        execution_unit.num_entries, execution_unit.entry_dt = self._metadata_reader.get_time_box_work_parameters(prev_exec_dt, exec_dt)
        if execution_unit.num_entries == 0:
            execution_unit.stop()
        self._logger.debug('End get_time_box_work')
        return execution_unit


class ExecutionUnitStateRunner(StateRunner):
    """
    This class brings together the mechanisms for identifying the time-boxed lists of work to be done (DataSource
    specializations), and the mechanisms for translating a unit of work into something that CaomExecute can work with.

    For retries, accumulate the retry-able entries in a single file for each time-box interval, for each data source.
    After all the incremental execution, attempt the retries.
    """

    def __init__(
        self,
        config,
        organizer,
        data_sources,
        observable,
        reporter,
    ):
        super().__init__(
            config=config,
            organizer=organizer,
            builder=None,
            data_sources=data_sources,
            metadata_reader=None,
            observable=observable,
            reporter=reporter,
        )

    def _process_data_source(self, data_source):
        """
        Uses an iterable with an instance of StateRunnerMeta.

        :return: 0 for success, -1 for failure
        """
        data_source.initialize_start_dt()
        data_source.initialize_end_dt()
        prev_exec_time = data_source.start_dt
        incremented = increment_time(prev_exec_time, self._config.interval)
        exec_time = min(incremented, data_source.end_dt)

        self._logger.info(f'Starting at {prev_exec_time}, ending at {data_source.end_dt}')
        result = 0
        if prev_exec_time == data_source.end_dt:
            self._logger.info(f'Start time is the same as end time {prev_exec_time}, stopping.')
            exec_time = prev_exec_time
        else:
            cumulative = 0
            result = 0
            while exec_time <= data_source.end_dt:
                self._logger.info(f'Processing {data_source.start_key} from {prev_exec_time} to {exec_time}')
                save_time = exec_time
                self._reporter.set_log_location(self._config)
                work = data_source.get_time_box_work(prev_exec_time, exec_time)
                if work.num_entries > 0:
                    try:
                        # work.start()
                        self._logger.info(f'Processing {work.num_entries} entries.')
                        work.do()
                    finally:
                        work.stop()
                    save_time = min(work.entry_dt, exec_time)
                    self._record_retries()
                self._record_progress(work.num_entries, cumulative, prev_exec_time, save_time)
                data_source.save_start_dt(save_time)
                if exec_time == data_source.end_dt:
                    # the last interval will always have the exec time equal to the end time, which will fail the
                    # while check so leave after the last interval has been processed
                    #
                    # but the while <= check is required so that an interval smaller than exec_time -> end_time will
                    # get executed, so don't get rid of the '=' in the while loop comparison, just because this one
                    # exists
                    break
                prev_exec_time = exec_time
                new_time = increment_time(prev_exec_time, self._config.interval)
                exec_time = min(new_time, data_source.end_dt)
                cumulative += work.num_entries

        data_source.save_start_dt(exec_time)
        msg = f'Done for {data_source.start_key}, saved state is {exec_time}'
        self._logger.info('=' * len(msg))
        self._logger.info(msg)
        self._logger.info(f'{self._reporter.success} of {self._reporter.all} records processed correctly.')
        self._logger.info('=' * len(msg))
        self._logger.debug(f'End _process_data_source with result {result}')
        return result


class ExecutionUnit:
    """
    Could be:
    - 1 file
    - 1 rclone timebox
    - 1 group of files for horizontal scaling deployment

    Temporal Cohesion between logging setup/teardown and workspace setup/teardown.
    """

    def __init__(self, config, **kwargs):
        """
        :param root_directory str staging space location
        :param label str name of the execution unit. Should be unique and conform to posix directory naming standards.
        """
        self._log_fqn = None
        self._logging_level = None
        self._log_handler = None
        self._task_types = config.task_types
        self._config = config
        self._entry_dt = None
        self._clients = kwargs.get('clients')
        # self._data_source = kwargs.get('data_source')
        self._metadata_reader = kwargs.get('metadata_reader')
        self._observable = kwargs.get('observable')
        self._reporter = kwargs.get('reporter')
        self._prev_exec_dt = kwargs.get('prev_exec_dt')
        self._exec_dt = kwargs.get('exec_dt')
        self._label = (
            f'{self._prev_exec_dt.isoformat().replace(":", "_").replace(".", "_")}_'
            f'{self._exec_dt.isoformat().replace(":", "_").replace(".", "_")}'
        )
        self._working_directory = os.path.join(config.working_directory, self._label)
        if config.log_to_file:
            if config.log_file_directory:
                self._log_fqn = os.path.join(config.log_file_directory, self._label)
            else:
                self._log_fqn = os.path.join(config.working_directory, self._label)
            self._logging_level = config.logging_level
        self._num_entries = None
        self._logger = logging.getLogger(self.__class__.__name__)

    @property
    def entry_dt(self):
        return self._entry_dt

    @entry_dt.setter
    def entry_dt(self, value):
        self._entry_dt = value

    @property
    def label(self):
        return self._label

    @property
    def num_entries(self):
        return self._num_entries

    @num_entries.setter
    def num_entries(self, value):
        self._num_entries = value

    @property
    def working_directory(self):
        return self._working_directory

    def do(self):
        """Make the execution unit one time-boxed copy from the DataSource to staging space, followed by a TodoRunner
        pointed to the staging space, and using that staging space with use_local_files: True. """
        self._logger.debug(f'Begin do for {self._num_entries} entries')
        result = None
        # set a Config instance to use the staging space with 'use_local_files: True'
        todo_config = Config()
        for attr in dir(self._config):
            if attr.startswith('__') or attr in [
                'bookmark', 'is_connected', 'report_fqn', 'time_zone', 'total_retry_fqn', 'use_vos'
            ]:
                continue
            value = getattr(self._config, attr)
            setattr(todo_config, attr, value)
        todo_config.use_local_files = True
        todo_config.data_sources = [self._working_directory]
        self._logger.debug(f'do config for TodoRunner: {todo_config}')
        organizer = OrganizeExecutes(
            todo_config,
            META_VISITORS,
            DATA_VISITORS,
            None,
            Transfer(),
            CadcTransfer(self._clients.cadc_client),
            self._metadata_reader,
            self._clients,
            self._observable,
            self._reporter,
        )
        local_data_source = ListDirSeparateDataSource(todo_config)
        local_data_source.reporter = self._reporter
        builder = EntryBuilder(PossumName)
        # start a TodoRunner with the new Config instance and the new data_source
        todo_runner = TodoRunner(
            todo_config,
            organizer,
            builder=builder,
            data_sources=[local_data_source],
            metadata_reader=self._metadata_reader,
            observable=self._observable,
            reporter=self._reporter,
        )
        result = todo_runner.run()
        if todo_config.cleanup_files_when_storing:
            result |= todo_runner.run_retry()
            todo_runner.report()
        self._logger.debug(f'End do with result {result}')
        return result

    def start(self):
        self._set_up_file_logging()
        self._create_workspace()

    def stop(self):
        self._clean_up_workspace()
        self._unset_file_logging()

    def _create_workspace(self):
        """Create the working area if it does not already exist."""
        self._logger.debug(f'Create working directory {self._working_directory}')
        create_dir(self._working_directory)

    def _clean_up_workspace(self):
        """Remove a directory and all its contents. Only do this if there is not a 'SCRAPE' task type, since the
        point of scraping is to be able to look at the pipeline execution artefacts once the processing is done.
        """
        if os.path.exists(self._working_directory) and TaskType.SCRAPE not in self._task_types and self._config.cleanup_files_when_storing:
            for ii in os.listdir(self._working_directory):
                os.remove(os.path.join(self._working_directory, ii))
            os.rmdir(self._working_directory)
            self._logger.debug(f'Removed working directory {self._working_directory} and contents.')
        self._logger.debug('End _clean_up_workspace')

    def _set_up_file_logging(self):
        """Configure logging to a separate file for each execution unit.

        If log_to_file is set to False, don't create a separate log file for each entry, because the application
        should leave as small a logging trace as possible.
        """
        if self._log_fqn and self._logging_level:
            self._log_handler = logging.FileHandler(self._log_fqn)
            formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(name)-12s:%(lineno)d:%(message)s')
            self._log_handler.setLevel(self._logging_level)
            self._log_handler.setFormatter(formatter)
            logging.getLogger().addHandler(self._log_handler)

    def _unset_file_logging(self):
        """Turn off the logging to the separate file for each entry being
        processed."""
        if self._log_handler:
            logging.getLogger().removeHandler(self._log_handler)
            self._log_handler.flush()
            self._log_handler.close()


class ExecutionUnitOrganizeExecutes(OrganizeExecutes):
    """A class to do nothing except be "not None" when called."""

    def __init__(self):
        pass

    def choose(self):
        # do nothing for the over-arching StateRunner
        pass

    def do_one(self, _):
        raise NotImplementedError


def remote_execution():
    """When running remotely, do a time-boxed 2-stage execution:
    1. stage 1 - the work is to use rclone to retrieve files to a staging area
    2. stage 2 - with the files in the staging area, use the pipeline as usual to store the files and create and
                 store the CAOM2 records, thumbnails, and previews

    Stage 1 is controlled with the ExecutionUnitStateRunner.
    Stage 2 is controlled with a TodoRunner, that is created for every ExecutionUnitStateRunner time-box that brings
    over files.
    """
    config = Config()
    config.get_executors()
    set_logging(config)
    observable = Observable(config)
    reporter = ExecutionReporter(config, observable)
    reporter.set_log_location(config)
    metadata_reader = RemoteMetadataReader()
    organizer = ExecutionUnitOrganizeExecutes()
    clients = RCloneClients(config)
    kwargs = {
        'clients': clients,
        'observable': observable,
        'reporter': reporter,
    }
    data_sources = []
    for entry in config.data_sources:
        data_source = RemoteIncrementalDataSource(
            config,
            entry,  # should look like "acacia_possum:pawsey0980" acacia_possum => rclone named config, pawsey0980 => root bucket
            metadata_reader,
            **kwargs,
        )
        data_source.reporter = reporter
        data_sources.append(data_source)
    runner = ExecutionUnitStateRunner(
        config,
        organizer,
        sources=data_sources,
        observable=observable,
        reporter=reporter,
    )
    result = runner.run()
    result |= runner.run_retry()
    runner.report()
    return result
