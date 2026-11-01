# OPAL (classic, with OPAL-cycl) and PyOpal, built from source.
# OPAL-2024.1.0 is the latest tag of github.com/OPALX-project/OPAL. OPALX runs from its own image (DESIGN.md section 25).
FROM ubuntu:24.04@sha256:fe81560990fae86dc4743202b2cb4d3df98d0bcf949173cecef15ce1127ef852

# Sources pinned to commits: H5hut-2.0.0rc7 and OPAL-2024.1.0.
ARG H5HUT_SHA=c41bd6ec549f1c4e7c90f8d52655054ab56bca57
ARG OPAL_SHA=fb50729ede0a4895ef0be8a7b1a1a26ebd035beb
ARG JOBS=4

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential cmake git ca-certificates autoconf automake libtool pkg-config \
        libopenmpi-dev openmpi-bin libhdf5-openmpi-dev libgsl-dev zlib1g-dev \
        libboost-filesystem-dev libboost-iostreams-dev libboost-regex-dev \
        libboost-serialization-dev libboost-system-dev libboost-timer-dev libboost-python-dev \
        python3-dev python3-numpy \
    && rm -rf /var/lib/apt/lists/*

ENV OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Debian splits parallel HDF5 into /usr/include/hdf5/openmpi and /usr/lib/<triplet>/hdf5/openmpi.
RUN git clone --filter=blob:none https://github.com/H5hut/H5hut.git /src/h5hut \
    && git -C /src/h5hut checkout --detach ${H5HUT_SHA} \
    && cd /src/h5hut && ./autogen.sh \
    && HDF5_INC=$(dirname $(find /usr/include -path '*openmpi*' -name hdf5.h | head -1)) \
    && HDF5_LIB=$(dirname $(find /usr/lib -path '*openmpi*' -name libhdf5.so | head -1)) \
    && CC=mpicc CXX=mpicxx CPPFLAGS="-I$HDF5_INC" LDFLAGS="-L$HDF5_LIB" LIBS="-lm" \
        ./configure --prefix=/opt/h5hut --enable-parallel --enable-shared \
    && make -j${JOBS} && make install

ENV H5HUT_PREFIX=/opt/h5hut

# OPAL-2024.1.0 sets -Werror; GCC 13 raises new warnings (OpenMPI C++ header, IPPL LField),
# so -Werror is dropped. Warnings stay visible in the build log.
# The single-particle orbit file is written with 9 digits; 17 are needed for the benchmarks.
RUN git clone --filter=blob:none https://github.com/OPALX-project/OPAL.git /src/opal \
    && git -C /src/opal checkout --detach ${OPAL_SHA} \
    && sed -i '/add_compile_options (-Werror)/d' /src/opal/CMakeLists.txt \
    && sed -i 's/outfTrackOrbit_m.precision(8);/outfTrackOrbit_m.precision(17);/' \
        /src/opal/src/Algorithms/ParallelCyclotronTracker.cpp \
    && cmake -S /src/opal -B /src/opal/build \
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/opal \
        -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
        -DWILL_BUILD_SHARED_LIBRARY=ON -DBUILD_OPAL_PYTHON=ON \
        -DPYOPAL_INSTALL_DIR=$(python3 -c 'import sysconfig; print(sysconfig.get_paths()["platlib"])') \
    && cmake --build /src/opal/build -j${JOBS} \
    && cmake --install /src/opal/build

ENV PATH=/opt/opal/bin:$PATH LD_LIBRARY_PATH=/opt/opal/lib:/opt/h5hut/lib
WORKDIR /data
