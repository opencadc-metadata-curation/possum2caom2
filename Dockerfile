ARG OPENCADC_PYTHON_VERSION=3.12
FROM opencadc/matplotlib:${OPENCADC_PYTHON_VERSION}-slim AS builder
ARG OPENCADC_PYTHON_VERSION

RUN apt-get update --no-install-recommends && \
    apt-get install -y build-essential busybox curl git libcfitsio-bin && \
    rm -rf /var/lib/apt/lists/ /tmp/* /var/tmp/*

WORKDIR /usr/src/app
ADD https://rclone.org/install.sh /usr/src/app
RUN chmod 755 /usr/src/app/install.sh && /bin/bash -c "/usr/src/app/install.sh"

ARG OPENCADC_BRANCH=main
ARG OPENCADC_REPO=opencadc

RUN git clone https://github.com/${OPENCADC_REPO}/caom2tools.git && \
    cd caom2tools && \
    git checkout ${OPENCADC_BRANCH} && \
    pip install ./caom2utils && \
    cd ..

RUN pip install git+https://github.com/${OPENCADC_REPO}/caom2pipe@${OPENCADC_BRANCH}#egg=caom2pipe

RUN pip install git+https://github.com/${OPENCADC_REPO}/possum2caom2@${OPENCADC_BRANCH}#egg=possum2caom2

RUN useradd --create-home --shell /bin/bash cadcops
USER cadcops

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

