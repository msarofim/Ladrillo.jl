FROM ubuntu:noble AS facts-core

# set env variables for RCT and apt
ENV DEBIAN_FRONTEND="noninteractive"

# set shell for RUN instruction
SHELL [ "/bin/bash", "-c" ]

# installs dependences
RUN apt-get update &&\
    apt-get install -y \
        python3-full \
        python3-pip \
        python3-venv \
        git \
        sudo \
        libhdf5-dev \
        libnetcdf-dev \
        gnupg \
        curl \
        apt-transport-https \
        dnsutils \
        iputils-ping

# Creates and activates python3 virtual environment
RUN python3 -m venv factsVe
RUN source /factsVe/bin/activate

# installs required FACTS packages
RUN /factsVe/bin/pip install --no-cache-dir --upgrade \
        setuptools==75.8.2 pip==25.0.1 wheel==0.45.1

RUN /factsVe/bin/pip install --no-cache-dir \
        radical.entk==1.92.0 \
        radical.gtod==1.90.0 \
        radical.pilot==1.92.0 \
        radical.saga==1.90.0 \
        radical.utils==1.91.1 \
        pyyaml==6.0.2 \
        xarray==2025.1.2 \
        numpy==2.2.3 \
        netcdf4==1.7.2 \
        h5py==3.13.0 \
        dask==2025.2.0

# Pre-install module python_dependencies (scipy/matplotlib/dask[array]) so the per-task
# pip installs FACTS runs at pre_exec become no-ops -> no concurrent-pip race on the shared
# /factsVe venv when many modules launch simultaneously (e.g. with an external/dummy climate step).
RUN /factsVe/bin/pip install --no-cache-dir scipy matplotlib 'dask[array]'

# Neuter pip. FACTS' per-task pre_exec runs "pip install --upgrade pip; pip install <deps>" with BARE
# pip, concurrently across tasks: the "--upgrade pip" rewrites pip itself -> races/corruption, and any
# task that doesn't activate the venv hits system pip -> "externally-managed-environment". ALL module
# python deps are baked above, so make every pip a no-op (exit 0): never installs, upgrades, or touches
# the venv. (If internal-FaIR is ever needed, bake fair==1.6.4 above instead of relying on runtime pip.)
RUN for p in /factsVe/bin/pip /factsVe/bin/pip3 /usr/bin/pip /usr/bin/pip3 /usr/local/bin/pip /usr/local/bin/pip3; do \
        if [ -e "$p" ] || [ -L "$p" ]; then printf '#!/bin/sh\nexit 0\n' > "$p" && chmod 0555 "$p"; fi; \
    done

# installs R for ubuntu 20.04 Focal
RUN apt-get update &&\
    apt-get install -y --no-install-recommends r-base \
        cmake \
        libopenblas-dev \
        gfortran


ARG NB_USER=jovyan
ARG NB_UID=1001
ENV USER=${NB_USER}
ENV NB_UID=${NB_UID}
ENV HOME=/home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

# Make sure the contents of our repo are in ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
RUN chown -R ${NB_UID} /factsVe
RUN chown -R ${NB_UID} /usr/local/lib/R/site-library
USER ${NB_USER}

# make directory for radical pilot sandbox
RUN mkdir -p ~/radical.pilot.sandbox

# Activate venv
RUN echo -e ". /factsVe/bin/activate" > ~/.bashrc

# create resource_local.json with the specified contents as well as
# create resource_local.json with the specified contents
RUN mkdir -p ~/.radical/pilot/configs
RUN printf '{\n\
    "local": {\n\
        "rp_version"    : "installed",\n\
        "virtenv_mode"  : "use"\n\
    }\n\
}' > ~/.radical/pilot/configs/resource_local.json

# return a shell to the user
CMD /bin/bash

FROM facts-core AS facts-jupyter

USER root
RUN /factsVe/bin/pip install --no-cache-dir \
	matplotlib==3.10.1 \
	notebook==7.3.2 \
	jupyterlab==4.3.5 \
	jupyter==1.1.1
RUN chown -R ${NB_UID} /factsVe
EXPOSE 8888
USER ${NB_USER}

#ENTRYPOINT ["jupyter", "lab","--ip=0.0.0.0","--allow-root"]
