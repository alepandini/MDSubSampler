# to rebuild the image, run:
# docker build -t subsampler .
FROM python:3.9-buster

RUN apt-get update && \
  apt-get install -y --no-install-recommends \
  gcc g++ graphviz

WORKDIR /app
COPY pyproject.toml pyproject.toml
RUN curl -sSL https://install.python-poetry.org | python3 -
ENV PATH="${PATH}:/root/.local/bin"
RUN poetry config virtualenvs.create false
RUN poetry install