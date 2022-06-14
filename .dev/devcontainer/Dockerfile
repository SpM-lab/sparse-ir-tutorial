FROM julia:1.7.2

# create user with a home directory
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    python3 \
    python3-dev \
    python3-distutils \
    curl \
    ca-certificates \
    git \
    wget \
    zip \
    libjpeg-dev \
    vim \
    openssh-server \
    tree \
    sudo \
    && \
    apt-get clean && rm -rf /var/cache/apt/archives/* /var/lib/apt/lists/* # clean up
RUN curl -kL https://bootstrap.pypa.io/get-pip.py | python3

# sudo
RUN echo "%${USER}    ALL=(ALL)   NOPASSWD:    ALL" >> /etc/sudoers.d/${USER} && \
    chmod 0440 /etc/sudoers.d/${USER}

WORKDIR ${HOME}
USER ${NB_USER}

RUN chown -R ${NB_UID} ${HOME}

# Install default Python libraries
RUN pip3 install \
    jupyter-book \
    jupytext \
    sparse-ir \
    xprec \
    matplotlib \
    ghp-import \
    --user

# Install default Julia packages
RUN mkdir -p ${HOME}/.julia/config && \
    echo '\
    # set environment variables\n\
    ENV["PYTHON"]=Sys.which("python3")\n\
    ENV["JUPYTER"]=Sys.which("jupyter")\n\
    ' >> ${HOME}/.julia/config/startup.jl && cat ${HOME}/.julia/config/startup.jl
#RUN julia -e 'using Pkg; Pkg.add(["IJulia", "Plots", "SparseIR", "Revise", "FFTW", "Roots", "OMEinsum", "GR", "FastGaussQuadrature", "LaTeXStrings"])'
RUN julia -e 'using Pkg; Pkg.add(["IJulia", "Plots", "Revise", "FFTW", "Roots", "OMEinsum", "GR", "FastGaussQuadrature", "LaTeXStrings"])'

RUN echo "export PATH=$HOME/.local/bin:\$PATH" >> ${HOME}/.bashrc

# Workdir for VS code
USER root
WORKDIR /workspace/sparse_ir_tutorial
RUN chown -R ${NB_UID} /workspace/sparse_ir_tutorial


WORKDIR ${HOME}
USER ${USER}
EXPOSE 8000

#RUN echo "alias python=python3" >> ${HOME}/.bashrc
COPY docker-entrypoint.sh /tmp
ENTRYPOINT ["/tmp/docker-entrypoint.sh"]

CMD ["julia"]