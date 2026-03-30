FROM condaforge/mambaforge:latest

USER root

# project dir
ENV PROJECT_DIR=/app
RUN mkdir -p $PROJECT_DIR
COPY . $PROJECT_DIR
WORKDIR $PROJECT_DIR

# build conda env
ENV ENV_PREFIX=$PROJECT_DIR/env
COPY --chown=$user:$user docker/environment_docker.yml /tmp/environment_docker.yml

COPY --chown=$user:$user docker/entrypoint.sh /usr/local/bin/
RUN chmod u+x /usr/local/bin/entrypoint.sh

# install environment
RUN mamba env update -n base -f /tmp/environment_docker.yml && \
    mamba clean --all --yes

# Pre-create Snakemake conda environments so they don't need to be built at runtime.
# Uses a dummy config - no rules are executed, only environments are resolved and installed.
COPY --chown=root:root docker/config_docker.yaml /tmp/config_docker.yaml
RUN touch dummy.txt && mkdir -p /dummy/db && \
    snakemake --conda-create-envs-only --use-conda --conda-frontend mamba \
    --cores 1 --configfile /tmp/config_docker.yaml && \
    rm dummy.txt && rm -rf /dummy/db && \
    mamba clean --all --yes

# download baktadb
RUN mkdir -p bakta_db && \
    bakta_db download --output bakta_db --type light
    
ENV PATH=$PATH:/opt/conda/bin

ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]
