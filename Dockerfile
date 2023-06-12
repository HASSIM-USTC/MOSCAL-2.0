FROM ubuntu:latest
COPY gperftools /home/gperf
COPY json11 /home/json11
COPY jemalloc /home/jemalloc
COPY folly /home/folly
COPY googletest /home/googletest
COPY gflags /home/gflags
COPY fmt /home/fmt
COPY eigen /home/eigen
RUN sed -i 's@//.*archive.ubuntu.com@//mirrors.ustc.edu.cn@g' /etc/apt/sources.list\
    && sed -i 's/security.ubuntu.com/mirrors.ustc.edu.cn/g' /etc/apt/sources.list\
    && apt-get update\
    && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata\
    && apt-get install -y libssl-dev make cmake gdb\
    curl libunwind-dev libblas-dev liblapack-dev liblapack3 libopenblas-base libopenblas-dev python3 python3-matplotlib python3-numpy python3-sympy graphviz ghostscript nano autoconf\    
    libunwind8-dev libboost-all-dev libevent-dev libdouble-conversion-dev libgoogle-glog-dev \
    libgflags-dev libiberty-dev liblz4-dev liblzma-dev libsnappy-dev zlib1g-dev binutils-dev \
    libjemalloc-dev pkg-config libunwind-dev libelf-dev libdwarf-dev
RUN cd /home/gperf \
    && ./configure \
    && make -j$(nproc)\
    && make install 
RUN cd /home/json11 \
    && cmake . \
    && make -j$(nproc)\
    && make install
RUN cd /home/googletest\
    && cmake .\
    && make -j6\
    && make install
RUN cd /home/fmt\
    && mkdir _build\
    && cd _build\
    && cmake ..\
    && make -j$(nproc)\
    && make install
RUN cd /home/gflags\
    && mkdir _build\
    && cd _build\
    && cmake ..\
    && make -j$(nproc)\
    && make install
RUN cd /home/folly\
    && mkdir _build && cd _build\
    && cmake configure ..\
    && make -j $(nproc)\
    && make install
RUN cd /home/eigen\
    && mkdir _build && cd _build\
    && cmake ..\
    && make install\
    && /sbin/ldconfig\
    && rm -rf /home/folly /home/googletest /home/fmt /home/eigen /home/gflags /home/gperf /home/json11 /home/jemalloc
COPY jemalloc /home/jemalloc
RUN apt-get install -y python3-pandas python3-scipy python3-matplotlib python3-numpy python3-sympy python3-pip \
    && python3 -m pip install --upgrade termcolor tqdm\
    && cd /home/jemalloc \
    && . ./autogen.sh\
    && make -j$(nproc)\
    && make install\
    && rm -rf /home/jemalloc