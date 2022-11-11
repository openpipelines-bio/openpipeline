FROM biocontainers/freebayes:v1.2.0-2-deb_cv1
USER root
RUN apt-get update && apt-get install -y python3.6
RUN apt-get install -y libvcflib-tools libvcflib-dev
RUN apt-get update && apt-get --no-install-recommends install -y parallel && \
	apt-get clean && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

