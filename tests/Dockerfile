# Dockerfile for running `workflow-pathoscope` tests.

FROM workflow-pathoscope

RUN apt-get update -y
RUN apt-get install -y gcc python3-dev

WORKDIR test

COPY . .
RUN pip install tox

ENTRYPOINT ["tox"]

