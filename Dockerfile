# Use full python for building
FROM python:3.11-bullseye as build

# Create virtualenv for deployment
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

# Use slim for deployment
FROM python:3.11-slim-bullseye

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
# Dev only
COPY --from=build /app_repository /app_repository 


USER appuser

# Set the entry point command
ENTRYPOINT /bin/bash
#CMD [ "bash" ]
