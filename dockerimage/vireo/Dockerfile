FROM python:3.9

RUN apt-get update && apt-get --no-install-recommends install -y python3-pip && \
	apt-get clean && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN pip3 install threadpoolctl vireoSNP
                    
