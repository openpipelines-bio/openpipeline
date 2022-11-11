FROM python:3.9

RUN apt-get update && apt-get --no-install-recommends install -y wget gcc make libbz2-dev \
        zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev python3-pip && \
	apt-get clean && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /usr/src

RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
	tar jxf samtools-1.16.1.tar.bz2 && \
	rm samtools-1.16.1.tar.bz2 && \
	cd samtools-1.16.1 && \
	./configure --prefix $(pwd) && \
	make

RUN pip3 install umi_tools

ENV PATH=${PATH}:/usr/src/samtools-1.16.1 

