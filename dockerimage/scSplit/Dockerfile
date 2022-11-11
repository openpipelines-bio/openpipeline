FROM python:3.6

RUN apt-get update && apt-get --no-install-recommends install -y python3-pip && \
	apt-get clean && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN python -m pip install --upgrade pip
RUN pip3 install numpy pandas pysam PyVCF scikit-learn scipy statistics scSplit
                    
