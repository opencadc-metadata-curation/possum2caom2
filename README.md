# possum2caom2
Starting point to build an application to generate CAOM2 Observations from FITS files.

# How To Run POSSUM Testing

In an empty directory:

1. This is the working directory, so it should probably have some space.

1. In the `main` branch of this repository, find the file `Dockerfile`. In the
`scripts` directory, find the files
`docker-entrypoint.sh`, and `config.yml`. Copy these files to the working directory.

   ```
   wget https://raw.github.com/opencadc-metadata-curation/possum2caom2/main/Dockerfile
   wget https://raw.github.com/opencadc-metadata-curation/possum2caom2/main/scripts/docker-entrypoint.sh
   wget https://raw.github.com/opencadc-metadata-curation/possum2caom2/main/scripts/config.yml
   ```

1. Make `docker-entrypoint.sh` executable.

1. `config.yml` is configuration information for the ingestion. It will work with
the files named and described here. For a complete description of its
content, see
https://github.com/opencadc-metadata-curation/collection2caom2/wiki/config.yml.

1. The ways to tell this tool the work to be done:

   1. provide a file containing the list of file ids to process, one file id
   per line, and the config.yml file containing the entries 'use_local_files'
   set to False, and 'task_types' set to -ingest -modify. The 'todo'
   file may provided in one of two ways:
      1. named 'todo.txt' in this directory, as specified in config.yml, or
      1. as the fully-qualified name with the --todo parameter

   1. provide the files to be processed in the working directory, and the
   config.yml file containing the entries 'use_local_files' set to True,
   and 'task_types' set to -store -ingest -modify.
      1. The store task does not have to be present, unless the files on disk
      are newer than the same files at CADC.

   1. provide the files to be processed in a Pawsey acacia remote, with the
   config.yml file containing the entries 'use_local_files' set to False,
   the 'task_types' set to -store -ingest -modify, and the `data_sources` set to
   `<acacia remote>/possum/tiles`.
      1. The `data_sources` entry requires the `rclone` configuration to be set up.
      Use `pawsey` in the `acacia remote` name to get the correct `rclone` syntax
      for the commands.

   1. provide the files to be processed in the working directory, and the
   config.yml file containing the entries 'use_local_files' set to True,
   and 'task_types' set to -scrape.
      1. This configuration will not attempt to write files or CAOM2 records
      to CADC. It is a good way to craft the content of the CAOM2 record without
      continually updating database content.

1. To build the container image, run this:

   ```
   docker build -f Dockerfile -t possum_run_cli ./
   ```

1. In the working directory, place a CADC proxy certificate. The Docker image can be used to create a
proxy certificate as follows. You will be prompted for the password for your CADC user:

   ```
   user@dockerhost:<cwd># docker run --rm -ti -v ${PWD}:/usr/src/app  -v <fully-qualified path to staging directory>:/data --user $(id -u):$(id -g) -e HOME=/usr/src/app --name possum_run_cli possum_run cadc-get-cert --days-valid 10 --cert-filename /usr/src/app/cadcproxy.pem -u <your CADC username>
   ```

1. To set up the `rclone` configuration for Pawsey s3 acacia storage, run the image, and then run `rclone config` from within the image. Follow `rclone config` steps as described here: https://www.youtube.com/watch?v=mOp7NJpwzac&t=1507s. Note that this will leave the `.config/rclone/rclone.conf` file on disk, which is why the last step is to set permissions on the `rclone` configuration file:

   ```
   user@dockerhost:<cwd># docker run --rm -ti -v <cwd>:/usr/src/app --user $(id -u):$(id -g) -e HOME=/usr/src/app --name possum_run_cli possum_run_cli /bin/bash
   cadcops@d51a02720ea6:~$ rclone config
   cadcops@d51a02720ea6:~$ <follow the rclone config steps>
   cadcops@d51a02720ea6:~$ exit
   user@dockerhost:<cwd># chmod 600 .config/rclone/rclone.conf
   ```

1. To run the application where it will retrieve files from the remote Pawsey s3 acacia storage:

   ```
   user@dockerhost:<cwd># docker run --rm -ti -v <cwd>:/usr/src/app --user $(id -u):$(id -g) -e HOME=/usr/src/app --name possum_run_cli possum_run_cli possum_run_remote
   ```

1. To edit and test the application from inside a container:

   ```
   user@dockerhost:<cwd># git clone https://github.com/opencadc-metadata-curation/possum2caom2.git
   user@dockerhost:<cwd># docker run --rm -ti -v <cwd>:/usr/src/app -v <fully-qualified path to staging directory>:/data --user $(id -u):$(id -g) -e HOME=/usr/src/app --name possum_run_cli possum_run_cli /bin/bash
   root@53bef30d8af3:/usr/src/app# pip install -e ./possum2caom2
   root@53bef30d8af3:/usr/src/app# pip install mock pytest
   root@53bef30d8af3:/usr/src/app# cd possum2caom2/possum2caom2/tests
   root@53bef30d8af3:/usr/src/app# pytest
   ```

1. For some instructions that might be helpful on using containers, see:
https://github.com/opencadc-metadata-curation/collection2caom2/wiki/Docker-and-Collections

