# based on the dockerfile from https://github.com/statgen/popscle

# Source Image
FROM ubuntu:18.04

# Set the working directory to /popscle_out
# WORKDIR /popscle_out


# supposed run in a separate folder docker_conda
# Copy the current directory contents into the container 
# COPY . .

# Set noninterative mode
ENV DEBIAN_FRONTEND noninteractive

# apt update and install global requirements
RUN apt-get clean all && \
     apt-get update && \
     apt-get upgrade -y && \
     apt-get install -y  \
         autoconf \
         build-essential \
         cmake \
         git \
         wget \
         libbz2-dev \
         libcurl4-openssl-dev \
         libssl-dev \
         zlib1g-dev \
         liblzma-dev

# apt clean and remove cached source lists
RUN apt-get clean && \
     rm -rf /var/lib/apt/lists/*

 RUN git clone https://github.com/samtools/htslib.git
 RUN cd htslib && \
     git submodule update --init --recursive && \
     autoheader && \
     autoconf && \
     ./configure --prefix=/usr/local/ && \
     make && \
     make install
  
# Install popscle
RUN git clone https://github.com/statgen/popscle.git
RUN cd popscle && \
     mkdir build && \
     cd build && \
     cmake .. && \
     make
RUN cp /popscle/bin/popscle /usr/local/bin

ENV PATH=${PATH}:/usr/local/bin
