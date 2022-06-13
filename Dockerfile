FROM shinaoka/sparse-ir-tutorial-binder:latest

ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

#COPY . ${HOME}
COPY ./README.md ${HOME}/README.md
COPY ./src ${HOME}/src

USER root
RUN chown -R ${NB_UID} ${HOME}
