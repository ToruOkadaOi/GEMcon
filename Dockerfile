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

# COPY env file into /tmp
COPY envs/ /tmp/env_specs/

# Copying the entire project from prj root
COPY . /workspace

# permissions to allow writing
RUN chown -R $MAMBA_USER:$MAMBA_USER /workspace

# Create the gemcon env, pin troppo/cobamp without dep resolution
# (cobamp 0.2.1 has a stale numpy pin), install the CLI, clean up.
RUN micromamba create -y -n gemcon -f /tmp/env_specs/environment.yml && \
    micromamba run -n gemcon pip install --no-deps --force-reinstall \
        troppo==0.0.7 cobamp==0.2.1 && \
    micromamba run -n gemcon python -m pip install -e . && \
    micromamba clean -afy && \
    rm -rf /tmp/env_specs

# Suppress Prefect telemetry and noise
ENV DO_NOT_TRACK=1
ENV PREFECT_LOGGING_LEVEL=ERROR

# Activate gemcon during RUN commands (and at container start)
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV ENV_NAME=gemcon

# revert
USER $MAMBA_USER

CMD ["bash"]