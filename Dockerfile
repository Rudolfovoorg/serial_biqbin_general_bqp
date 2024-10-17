FROM ubuntu:22.04

RUN apt-get update
RUN apt install build-essential -y
RUN apt-get install libopenblas-dev -y


WORKDIR /solver
COPY . .

RUN chmod -R 777 .

# Run it as some non root user 
USER 1001:1001

RUN make clean
RUN make
RUN make test



