#!/bin/bash

# =============================================================================
# CLOUD BUILD SCRIPT FOR VIRAL SCREENING PIPELINE
# =============================================================================
# 
# DESCRIPTION:
# Simple script to build the viral screening Docker container on Linux
# for deployment to DNAnexus and other cloud platforms.
#
# USAGE:
# ./build_cloud.sh [tag]
# Example: ./build_cloud.sh viral-screen:v1.0.0
#
# AUTHOR: Viral-GWAS Project
# VERSION: 1.0
# DATE: January 2025
# =============================================================================

set -euo pipefail

# Default tag
DEFAULT_TAG="viral-screen:latest"
TAG="${1:-$DEFAULT_TAG}"

echo "=========================================="
echo "BUILDING VIRAL SCREENING PIPELINE"
echo "=========================================="
echo "Tag: $TAG"
echo "Start time: $(date)"
echo "=========================================="

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
    echo "Error: Docker is not running"
    exit 1
fi

# Build the Docker image
echo "Building Docker image..."
docker build --tag "$TAG" .

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