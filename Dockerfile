# Multi-platform Docker image for tskibd
FROM ubuntu:22.04 AS builder

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git \
    ca-certificates \
    build-essential \
    g++ \
    pkg-config \
    meson \
    ninja-build \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /build

# Clone the repository with submodules
RUN git clone --recursive https://github.com/bguo068/tskibd.git . && \
    if [ -f tskit/c/VERSION ]; then mv tskit/c/VERSION tskit/c/VERSION.txt; fi

# Build tskibd
RUN meson setup build && \
    ninja -C build tskibd

# Runtime stage - smaller image
FROM ubuntu:22.04

# Install runtime dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libstdc++6 \
    && rm -rf /var/lib/apt/lists/*

# Copy the compiled binary
COPY --from=builder /build/build/tskibd /usr/local/bin/tskibd

# Create working directory for data
WORKDIR /data

# Set entrypoint
ENTRYPOINT ["tskibd"]
CMD ["--help"]
