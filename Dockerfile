FROM python:3.11.11-slim-bookworm AS builder

# Prevent Python from writing .pyc files and enable unbuffered logging
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Install build-time dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    build-essential \
    cmake \
    zlib1g-dev \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# Create a virtual environment to isolate the installation
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install uncalled4 directly from the GitHub repository into the venv
RUN pip install --no-cache-dir git+https://github.com/skovaka/uncalled4.git

FROM python:3.11.11-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PATH="/opt/venv/bin:$PATH"

# Install runtime dependencies only
RUN apt-get update && apt-get install -y --no-install-recommends \
    samtools \
    libhdf5-103-1 \
    && rm -rf /var/lib/apt/lists/*

# Copy the pre-built virtual environment from the builder
COPY --from=builder /opt/venv /opt/venv

# Create a working directory for data analysis
WORKDIR /data

# Create a non-root user and grant access to the working directory
RUN useradd -m -u 1000 appuser && chown appuser:appuser /data
USER appuser

# Set the default entrypoint to the uncalled4 executable
ENTRYPOINT ["uncalled4"]
CMD ["--help"]