FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive 
RUN apt-get update && apt-get install -y --no-install-recommends r-base build-essential
RUN apt install software-properties-common -y 
RUN add-apt-repository "ppa:deadsnakes/ppa"
RUN apt install -y python3.8