# micromamba base image. #TODO: Should I use mamba instead??
FROM mambaorg/micromamba:2.4.0

# setting to root for installing add. pip packages
USER root

# Work directory inside the container
WORKDIR /workspace

# Copying the entire project from prj root
COPY . /workspace

# COPY env files into /tmp
COPY envs/ /tmp/env_specs/

# permissions to allow writing
RUN chown -R $MAMBA_USER:$MAMBA_USER /workspace

# Create all conda environments
RUN micromamba create -y -n scanpy_legacy -f /tmp/env_specs/scanpy_env.yml && \
    micromamba create -y -n cplex_aman    -f /tmp/env_specs/cplex_env.yml    && \
    micromamba create -y -n gecko_aman    -f /tmp/env_specs/gecko_env.yml

# Activate; ## To have an environment active during a RUN command; ### https://micromamba-docker.readthedocs.io/en/latest/quick_start.html
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# revert
USER $MAMBA_USER

# No ENTRYPOINT since multiple envs
CMD ["bash"]