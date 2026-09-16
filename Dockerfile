# Use full python for building
FROM python:3.11-bullseye as build

# Create virtualenv for deployment (python3 executable will be in here and will give /opt/venv/lib precedence over PYTHONPATH)
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Set working directory
WORKDIR /app_repository

# Setup
COPY README.md ./
COPY setup.py ./
COPY epiread_tools/ ./epiread_tools/
RUN python3 setup.py install
# Dev only
COPY tests/ ./tests/

# Deconvolution models
COPY deconvolution_models/ ./deconvolution_models/
RUN pip install --upgrade pip
RUN pip install ./deconvolution_models

#COPY deconvolution_models/deconvolution_models/ ./deconvolution_models/
# COPY deconvolution_models/setup.py ./
# RUN python3 setup.py install


# Wanted to use slim for deployment, but did not have basic
# things needed by nextflow (like "ps")
#FROM python:3.11-slim-bullseye
FROM python:3.11-bullseye

# Install required system dependencies
RUN apt-get update && \
    apt-get -y install \
    build-essential \
    unzip

# Create user and switch to it
RUN groupadd -g 1002 appuser && useradd -u 1002 -g appuser -m appuser
WORKDIR /home/appuser

# Use the previously built virtualenv
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
COPY --from=build /opt/venv /opt/venv

# Debugging only
#COPY --from=build /app_repository /app_repository 

# Dev only. REMOVE FOR PRODUCTION
# NOTE: IF YOU RUN /opt/venv/bin/deconvolution, it will still get the /opt/venv libraries
ENV PYTHONPATH="/nf-core/epiread-tools/deconvolution_models/deconvolution_models:/nf-core/epiread-tools/epiread-tools:$PYTHONPATH"
COPY deconvolution_models/demo ./demo/

USER appuser

# Nextflow did not like setting an entrypoint.
#ENTRYPOINT /bin/bash
