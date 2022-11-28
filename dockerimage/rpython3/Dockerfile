FROM python:3.8

RUN apt-get update && apt-get install -y --no-install-recommends r-base build-essential libstdc++6 
RUN apt-get install software-properties-common -y
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 1E9377A2BA9EF27F
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get update 
RUN apt-get upgrade -y
RUN apt-get dist-upgrade