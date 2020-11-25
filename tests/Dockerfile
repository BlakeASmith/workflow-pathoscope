# Dockerfile for running `workflow-pathoscope` tests.

FROM python:3.9-buster

WORKDIR test

COPY . .
RUN pip install tox

ENTRYPOINT ["tox"]

