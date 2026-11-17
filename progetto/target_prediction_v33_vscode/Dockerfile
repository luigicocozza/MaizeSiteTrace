FROM python:3.11-slim

WORKDIR /app

RUN pip install --no-cache-dir streamlit==1.50.0

COPY app_v33.py /app/app_v33.py

EXPOSE 8502

CMD ["streamlit", "run", "/app/app_v33.py", "--server.address=0.0.0.0", "--server.port=8502", "--server.headless=true", "--server.maxUploadSize=250"]
