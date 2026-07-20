FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        gfortran \
        ninja-build \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .

RUN python -m pip install --upgrade pip \
    && python -m pip install --no-cache-dir ".[server]"

CMD ["sh", "-c", "exec gunicorn PyMieSim.application.wsgi:server --bind 0.0.0.0:${PORT:-8050} --workers 1 --threads 2 --timeout 3600"]

