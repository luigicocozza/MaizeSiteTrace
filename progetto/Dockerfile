FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV STREAMLIT_SERVER_HEADLESS=true
ENV STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

COPY requirements_docker_runtime.txt /app/requirements_docker_runtime.txt

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r /app/requirements_docker_runtime.txt

COPY app /app/app
COPY api /app/api
COPY docker /app/docker

RUN mkdir -p /app/data/processed /app/outputs /app/models /app/logs /app/docs && \
    chmod +x /app/docker/start_services.sh

EXPOSE 8000
EXPOSE 8501

CMD ["/app/docker/start_services.sh"]
