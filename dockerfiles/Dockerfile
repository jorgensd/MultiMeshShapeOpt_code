# This is a Dockerfile used to install the necessary dependencies for
# running the examples in the MultiMesh Shape Optimization article
# Ipopt installation based on pyadjoints installation scrips of Ipopt,
# made by Simon W. Funke.
#
# Authors:
# Jørgen S. Dokken <dokken92@gmail.com>

FROM quay.io/fenicsproject/stable:2018.1.0.r3

USER fenics
ENV IPOPT_VER=3.12.9
ENV GMSH_VER=3.0.6
ENV FEMORPH_BRANCH=dokken/restructuring


COPY ipopt.conf $FENICS_HOME/ipopt.conf
COPY gmsh.conf $FENICS_HOME/gmsh.conf
COPY femorph.conf $FENICS_HOME/femorph.conf
COPY mmshapeopt.conf $FENICS_HOME/mmshapeopt.conf

USER root
RUN apt-get update && \
apt-get install -y libgl1-mesa-glx libxcursor1 libxft2 libxinerama1 libglu1-mesa && \
apt-get clean

RUN /bin/bash -l -c "source $FENICS_HOME/mmshapeopt.conf"

RUN /bin/bash -l -c "source $FENICS_HOME/mmshapeopt.conf &&\
                     source $FENICS_HOME/femorph.conf &&\
                     install_femorph "

RUN /bin/bash -l -c "source $FENICS_HOME/gmsh.conf &&\
                     get_gmsh "

RUN /bin/bash -l -c "source $FENICS_HOME/mmshapeopt.conf && \
                     source $FENICS_HOME/ipopt.conf && \
                     update_ipopt && \
                     update_pyipopt"

USER fenics
RUN echo "source $FENICS_HOME/mmshapeopt.conf" >> $FENICS_HOME/.bash_profile
RUN /bin/bash -l -c "python3 -c \"import dolfin; import pyipopt\""
RUN /bin/bash -l -c "python3 -c \"import femorph\""
RUN /bin/bash -l -c "gmsh -version"

USER root
