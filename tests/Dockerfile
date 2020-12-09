# Dockerfile for running `workflow-pathoscope` tests.

FROM workflow-pathoscope

RUN apt-get update -y
RUN apt-get install -y gcc python3-dev
RUN pip install pytest pytest-asyncio pytest-aiohttp pytest-cov pytest-xdist

WORKDIR test

COPY . .


ENTRYPOINT ["pytest"]