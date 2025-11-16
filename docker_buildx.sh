#! /usr/bin/env bash
docker buildx build \
  --platform linux/amd64,linux/arm64 \
  -t bguo068/tskibd:v0.0.2 \
  -t bguo068/tskibd:latest \
  --push .
