#!/bin/bash

# for building FTorch
# SPDX-FileCopyrightText: Copyright (c) 2023 - 2024 NVIDIA CORPORATION & AFFILIATES.
# SPDX-FileCopyrightText: All rights reserved.
# SPDX-License-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# export CUDA_PATH=/usr/local/cuda
export PY_SITE_PATH=/usr/local/lib/python3.10/dist-packages
# export PY_SITE_PATH=/media/peter/samlinux/miniconda3/envs/torchcu124/lib/python3.11/site-packages
export CMAKE_PREFIX_PATH="${PY_SITE_PATH}/torch/"
echo "Ftorch install. Trying to use Python libraries from $PY_SITE_PATH"

CONFIG=Release
DEVICETYPE=CUDA

INSTALL_PATH=/opt/ftorch
# export Torch_DIR="${PY_SITE_PATH}/torch/share/cmake/Torch/"
# cd FTorch/

mkdir build 
cd build 

cmake .. \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH \
    -DCMAKE_INSTALL_PREFIX=$INSTALL_PATH \
    -DCMAKE_BUILD_TYPE=$CONFIG \
    -DGPU_DEVICE=$DEVICETYPE \
    -DMULTI_GPU=OFF

cmake --build . --target install
