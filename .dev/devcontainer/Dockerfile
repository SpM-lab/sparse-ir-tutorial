FROM shinaoka/sparse-ir

# User
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

# Workdir for VS code
USER root
WORKDIR /workspace/sparse_ir_tutorial
RUN chown -R ${NB_UID} /workspace/sparse_ir_tutorial

USER ${NB_USER}

EXPOSE 8888

CMD ["julia"]
