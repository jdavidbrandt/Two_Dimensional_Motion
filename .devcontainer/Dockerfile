FROM mcr.microsoft.com/devcontainers/cpp:1-ubuntu-24.04
LABEL maintainer="Joshua Brandt <jdbrandt92@gmail.com>"
ARG KEY_NAME="id_ed25519"

ARG REINSTALL_CMAKE_VERSION_FROM_SOURCE="3.22.2"

COPY ./reinstall-cmake.sh /tmp/

RUN if [ "${REINSTALL_CMAKE_VERSION_FROM_SOURCE}" != "none" ]; then \
        chmod +x /tmp/reinstall-cmake.sh && /tmp/reinstall-cmake.sh ${REINSTALL_CMAKE_VERSION_FROM_SOURCE}; \
    fi \
    && rm -f /tmp/reinstall-cmake.sh



#Barebones install
RUN apt-get update && \
    apt-get install -y sudo curl gnupg gnuplot \
    build-essential git locales nodejs wget neovim npm fonts-powerline && \
    locale-gen en_US.UTF-8 && \
    adduser --quiet --disabled-password \
    --shell /bin/sh --home /home/devuser/ \
    --gecos "User" devuser
    
#Build and install raylib with dependencies
RUN apt-get update && \
    apt-get install -y libasound2-dev libx11-dev libxrandr-dev libxi-dev libgl1-mesa-dev libglu1-mesa-dev \
    libxcursor-dev libxinerama-dev libwayland-dev libxkbcommon-dev && \
    git clone https://github.com/raylib-extras/raylib-quickstart.git two_dimensional_motion

