#!/bin/bash

# =============================================================================
# VIRAL SCREENING PIPELINE DOCKER BUILD SCRIPT
# =============================================================================
# 
# DESCRIPTION:
# Builds the viral screening pipeline Docker container with proper optimization
# and tagging for deployment to DNAnexus.
#
# USAGE:
# ./build.sh [tag]
# Example: ./build.sh v1.0.0
#
# AUTHOR: Viral-GWAS Project
# VERSION: 1.0
# DATE: January 2025
# =============================================================================

set -euo pipefail

# Default tag
DEFAULT_TAG="viral-screen:latest"
TAG="${1:-$DEFAULT_TAG}"

# Build context
BUILD_CONTEXT="."

echo "=========================================="
echo "BUILDING VIRAL SCREENING PIPELINE DOCKER"
echo "=========================================="
echo "Tag: $TAG"
echo "Build context: $BUILD_CONTEXT"
echo "Start time: $(date)"
echo "=========================================="

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
    echo "Error: Docker is not running or not accessible"
    exit 1
fi

# Build the Docker image with optimization
echo "Building Docker image..."
docker build \
    --tag "$TAG" \
    --build-arg BUILDKIT_INLINE_CACHE=1 \
    --progress=plain \
    --no-cache \
    "$BUILD_CONTEXT"

# Check if build was successful
if [ $? -eq 0 ]; then
    echo "=========================================="
    echo "BUILD SUCCESSFUL"
    echo "=========================================="
    echo "Image: $TAG"
    echo "End time: $(date)"
    
    # Show image info
    echo ""
    echo "Image information:"
    docker images "$TAG"
    
    # Show image size
    IMAGE_SIZE=$(docker images --format "table {{.Repository}}:{{.Tag}}\t{{.Size}}" | grep "$TAG" | awk '{print $2}')
    echo "Image size: $IMAGE_SIZE"
    
    echo ""
    echo "To test the container:"
    echo "docker run --rm -v /path/to/data:/data -v /path/to/results:/results $TAG"
    echo ""
    echo "To push to registry:"
    echo "docker tag $TAG your-registry/$TAG"
    echo "docker push your-registry/$TAG"
    
else
    echo "=========================================="
    echo "BUILD FAILED"
    echo "=========================================="
    exit 1
fi 