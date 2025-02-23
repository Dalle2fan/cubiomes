#!/bin/bash

# Exit on any error
set -e

echo "Checking cubiomes library..."

# Check if git is installed
if ! command -v git &> /dev/null; then
    echo "Error: git is not installed"
    exit 1
fi

# Check if gcc is installed
if ! command -v gcc &> /dev/null; then
    echo "Error: gcc is not installed"
    exit 1
fi

# Check if make is installed
if ! command -v make &> /dev/null; then
    echo "Error: make is not installed"
    exit 1
fi

# Clone the repository if it doesn't exist
if [ ! -d "cubiomes" ]; then
    echo "Cloning cubiomes repository..."
    git clone https://github.com/Cubitect/cubiomes.git
fi

# Change to the cubiomes directory
cd cubiomes

# Only build if library doesn't exist
if [ ! -f "libcubiomes.a" ]; then
    echo "Building library..."
    make clean
    make

    # Verify the build artifacts exist
    if [ ! -f "libcubiomes.a" ]; then
        echo "Error: Static library build failed"
        exit 1
    fi

    echo "Build completed successfully!"
    echo "Build artifacts:"
    ls -l libcubiomes.a
else
    echo "Library already built, skipping compilation"
fi

cd ..