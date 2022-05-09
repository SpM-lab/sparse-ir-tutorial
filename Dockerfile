FROM julia:1.7.0

# create user with a home directory
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

USER root

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
    && \
    apt-get clean && rm -rf /var/cache/apt/archives/* /var/lib/apt/lists/* # clean up


RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    vim \
    openssh-server \
    tree \
    && \
    apt-get clean && rm -rf /var/cache/apt/archives/* /var/lib/apt/lists/* # clean up

RUN curl -kL https://bootstrap.pypa.io/get-pip.py | python3 && \
    pip3 install \
    jupyter-book \
    sparse-ir \
    xprec \
    matplotlib \
    && \
    echo Done

WORKDIR ${HOME}
USER ${USER}

RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

RUN mkdir -p ${HOME}/.julia/config && \
    echo '\
    # set environment variables\n\
    ENV["PYTHON"]=Sys.which("python3")\n\
    ENV["JUPYTER"]=Sys.which("jupyter")\n\
    ' >> ${HOME}/.julia/config/startup.jl && cat ${HOME}/.julia/config/startup.jl

RUN julia -e 'using Pkg; Pkg.add(["IJulia", "PyPlot", "Plots", "SparseIR", "Revise", "FFTW", "Roots", "OMEinsum", "GR"])'


USER ${USER}

WORKDIR /workspace/sparse_ir_tutorial

COPY ./src /workspace/sparse_ir_tutorial/src

USER root
RUN chown -R ${NB_UID} /workspace/sparse_ir_tutorial
USER ${USER}


WORKDIR ${HOME}
USER ${USER}
EXPOSE 8000

CMD ["julia"]
