# micromamba base image
FROM mambaorg/micromamba:2.4.0

# setting to root for installing add. pip packages
USER root

# Editor, curl, tmux, git etc. + cmake and build-essential for building packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    nano git curl \
    cmake build-essential && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Work directory inside the container
WORKDIR /workspace

# COPY env files into /tmp
COPY envs/ /tmp/env_specs/

# Copying the entire project from prj root
COPY . /workspace

# permissions to allow writing
RUN chown -R $MAMBA_USER:$MAMBA_USER /workspace

# Create all conda environments, install base utils, install cli, and clean up in one layer
RUN micromamba create -y -n scanpy_legacy -f /tmp/env_specs/scanpy_env.yml && \
    micromamba create -y -n cplex_aman_new -f /tmp/env_specs/cplex_env.yml && \
    micromamba create -y -n gecko_aman    -f /tmp/env_specs/gecko_env.yml    && \
    micromamba run -n gecko_aman pip install --force-reinstall --no-cache-dir \
        git+https://github.com/ginkgobioworks/geckopy.git && \
    micromamba run -n cplex_aman_new pip install "numpy<2" && \
    micromamba install -y python=3.10 pyyaml rich typer prefect && \
    micromamba run -n base python -m pip install -e . && \
    micromamba clean -afy && \
    rm -rf /tmp/env_specs

# Suppress Prefect telemetry and noise
ENV DO_NOT_TRACK=1
ENV PREFECT_LOGGING_LEVEL=ERROR

# Activate; ## To have an environment active during a RUN command
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# revert
USER $MAMBA_USER

# No ENTRYPOINT since multiple envs
CMD ["bash"]
