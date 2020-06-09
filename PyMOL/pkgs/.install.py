# (c) 2012-2016 Anaconda, Inc. / https://anaconda.com
# All Rights Reserved
#
# conda is distributed under the terms of the BSD 3-clause license.
# Consult LICENSE.txt or http://opensource.org/licenses/BSD-3-Clause.
'''
We use the following conventions in this module:

    dist:        canonical package name, e.g. 'numpy-1.6.2-py26_0'

    ROOT_PREFIX: the prefix to the root environment, e.g. /opt/anaconda

    PKGS_DIR:    the "package cache directory", e.g. '/opt/anaconda/pkgs'
                 this is always equal to ROOT_PREFIX/pkgs

    prefix:      the prefix of a particular environment, which may also
                 be the root environment

Also, this module is directly invoked by the (self extracting) tarball
installer to create the initial environment, therefore it needs to be
standalone, i.e. not import any other parts of `conda` (only depend on
the standard library).
'''
import os
import re
import sys
import json
import shutil
import stat
from os.path import abspath, dirname, exists, isdir, isfile, islink, join
from optparse import OptionParser


on_win = bool(sys.platform == 'win32')
try:
    FORCE = bool(int(os.getenv('FORCE', 0)))
except ValueError:
    FORCE = False

LINK_HARD = 1
LINK_SOFT = 2  # never used during the install process
LINK_COPY = 3
link_name_map = {
    LINK_HARD: 'hard-link',
    LINK_SOFT: 'soft-link',
    LINK_COPY: 'copy',
}
SPECIAL_ASCII = '$!&\%^|{}[]<>~`"\':;?@*#'

# these may be changed in main()
ROOT_PREFIX = sys.prefix
PKGS_DIR = join(ROOT_PREFIX, 'pkgs')
SKIP_SCRIPTS = False
IDISTS = {
  "apbs-1.5-1": {
    "md5": "d96bf8ff23f4d4b59cb9383a498f1afe",
    "url": "https://conda.anaconda.org/schrodinger/win-64/apbs-1.5-1.tar.bz2"
  },
  "asn1crypto-0.24.0-py37_0": {
    "md5": "3a4916054f7d730b17e5d8a4343d9b59",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/asn1crypto-0.24.0-py37_0.tar.bz2"
  },
  "biopython-1.73-py37he774522_0": {
    "md5": "1748806f1ff3d65e5c51b396faf5e038",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/biopython-1.73-py37he774522_0.tar.bz2"
  },
  "blas-1.0-mkl": {
    "md5": "e8aa6b7daaf0925245c148aaeaa0722e",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/blas-1.0-mkl.tar.bz2"
  },
  "ca-certificates-2019.5.15-0": {
    "md5": "6d581fc902e0485e9c3259773c6e43a9",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/ca-certificates-2019.5.15-0.tar.bz2"
  },
  "certifi-2019.3.9-py37_0": {
    "md5": "5807c9c1b64980de9b0d1e9807d1a4af",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/certifi-2019.3.9-py37_0.tar.bz2"
  },
  "cffi-1.12.3-py37h7a1dbc1_0": {
    "md5": "2d6e69aff537a2d370ef906ab19bda48",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/cffi-1.12.3-py37h7a1dbc1_0.tar.bz2"
  },
  "chardet-3.0.4-py37_1": {
    "md5": "fa614ba863adbffcc5d5d449a8a19add",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/chardet-3.0.4-py37_1.tar.bz2"
  },
  "conda-4.6.14-py37_0": {
    "md5": "5c4d29f4ba9e9ca552a77dfa9c130920",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/conda-4.6.14-py37_0.tar.bz2"
  },
  "console_shortcut-0.1.1-3": {
    "md5": "aaad91749d7bb128bd51c0df273285e6",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/console_shortcut-0.1.1-3.tar.bz2"
  },
  "cryptography-2.7-py37h7a1dbc1_0": {
    "md5": "9e4c5385c50cf3e3c1657c85fac64217",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/cryptography-2.7-py37h7a1dbc1_0.tar.bz2"
  },
  "freemol-1.158-py37_1": {
    "md5": "1ee42a40adc98d49925c595982f440e9",
    "url": "https://conda.anaconda.org/schrodinger/win-64/freemol-1.158-py37_1.tar.bz2"
  },
  "freetype-2.9.1-ha9979f8_1": {
    "md5": "f3263f3249686ede2d2bea458e30e61f",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/freetype-2.9.1-ha9979f8_1.tar.bz2"
  },
  "glew-2.0.0-0": {
    "md5": "ffa7b14eaaf3f7bef79498558495f1d0",
    "url": "https://conda.anaconda.org/schrodinger/win-64/glew-2.0.0-0.tar.bz2"
  },
  "h5py-2.9.0-py37h5e291fa_0": {
    "md5": "60cd72ab81b1eb0842df40a93ab7ce5b",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/h5py-2.9.0-py37h5e291fa_0.tar.bz2"
  },
  "hdf5-1.10.4-h7ebc959_0": {
    "md5": "de2d3c434711af9c473efc9741189816",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/hdf5-1.10.4-h7ebc959_0.tar.bz2"
  },
  "icc_rt-2019.0.0-h0cc432a_1": {
    "md5": "ce2949c239b5cd45848a0a7865d30520",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/icc_rt-2019.0.0-h0cc432a_1.tar.bz2"
  },
  "icu-58.2-ha66f8fd_1": {
    "md5": "8b9a078b693623d58afd56f4f2a162c3",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/icu-58.2-ha66f8fd_1.tar.bz2"
  },
  "idna-2.8-py37_0": {
    "md5": "6afef2ec3bc9568cb43c92080008a802",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/idna-2.8-py37_0.tar.bz2"
  },
  "intel-openmp-2019.4-245": {
    "md5": "cd40f1846d3ada9ad3cfa61ef5fa7502",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/intel-openmp-2019.4-245.tar.bz2"
  },
  "jpeg-9b-hb83a4c4_2": {
    "md5": "7abafa9f9a2c609e1b77424f5ffa6a8a",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/jpeg-9b-hb83a4c4_2.tar.bz2"
  },
  "libiconv-1.15-h1df5818_7": {
    "md5": "cf3a121e3e744f0520d8fee6acdd2869",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/libiconv-1.15-h1df5818_7.tar.bz2"
  },
  "libpng-1.6.37-h2a8f88b_0": {
    "md5": "351328bf2bee7fd7ee7e9c158aeb1be2",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/libpng-1.6.37-h2a8f88b_0.tar.bz2"
  },
  "libtiff-4.0.10-hb898794_2": {
    "md5": "181742712383c04626ec76353dee85ba",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/libtiff-4.0.10-hb898794_2.tar.bz2"
  },
  "libxml2-2.9.9-h464c3ec_0": {
    "md5": "6ed19ec01f8b50c99b09e9134742a478",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/libxml2-2.9.9-h464c3ec_0.tar.bz2"
  },
  "mengine-1-0": {
    "md5": "530c0e49f31ad8cb643cee5d8a73ae01",
    "url": "https://conda.anaconda.org/schrodinger/win-64/mengine-1-0.tar.bz2"
  },
  "menuinst-1.4.16-py37he774522_0": {
    "md5": "90c96852dd3e41ee18c7c7def9c5f26e",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/menuinst-1.4.16-py37he774522_0.tar.bz2"
  },
  "mkl-2019.4-245": {
    "md5": "0f897fde2379b1adf2f6df4e79f22276",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/mkl-2019.4-245.tar.bz2"
  },
  "mkl_fft-1.0.12-py37h14836fe_0": {
    "md5": "879f59b0e6a920a546834d815aa21eda",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/mkl_fft-1.0.12-py37h14836fe_0.tar.bz2"
  },
  "mkl_random-1.0.2-py37h343c172_0": {
    "md5": "d45af154246685f482763506afd9cea0",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/mkl_random-1.0.2-py37h343c172_0.tar.bz2"
  },
  "mpeg_encode-1-0": {
    "md5": "c53591aeecb04fe948b7146652f3afe9",
    "url": "https://conda.anaconda.org/schrodinger/win-64/mpeg_encode-1-0.tar.bz2"
  },
  "mtz2ccp4_px-1.0-1": {
    "md5": "a4383b328d1c53a7ac7cab4305bb2e70",
    "url": "https://conda.anaconda.org/schrodinger/win-64/mtz2ccp4_px-1.0-1.tar.bz2"
  },
  "numpy-1.16.4-py37h19fb1c0_0": {
    "md5": "298f0baa621f915920d4d1fe92d61590",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/numpy-1.16.4-py37h19fb1c0_0.tar.bz2"
  },
  "numpy-base-1.16.4-py37hc3f5095_0": {
    "md5": "8869840ec723ace7061dcd53ff816eed",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/numpy-base-1.16.4-py37hc3f5095_0.tar.bz2"
  },
  "olefile-0.46-py37_0": {
    "md5": "1b98cdb7836ace32e351cbb9e497934d",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/olefile-0.46-py37_0.tar.bz2"
  },
  "openssl-1.1.1c-he774522_1": {
    "md5": "3cfeb411b7add6bc2e5c83d79f80deea",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/openssl-1.1.1c-he774522_1.tar.bz2"
  },
  "pdb2pqr-2.1.1-py37_1": {
    "md5": "410aecee085e29ebc92bfb3789649aa1",
    "url": "https://conda.anaconda.org/schrodinger/win-64/pdb2pqr-2.1.1-py37_1.tar.bz2"
  },
  "pillow-6.0.0-py37hdc69c19_0": {
    "md5": "247aa045caf28aecd6b383aaeeb6081c",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/pillow-6.0.0-py37hdc69c19_0.tar.bz2"
  },
  "pip-19.1.1-py37_0": {
    "md5": "90183042f1835c7930edc2779749650e",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/pip-19.1.1-py37_0.tar.bz2"
  },
  "pmw-2.0.1-py37_2": {
    "md5": "18f79820d9e2472f04a1a636fc37542e",
    "url": "https://conda.anaconda.org/schrodinger/win-64/pmw-2.0.1-py37_2.tar.bz2"
  },
  "pycosat-0.6.3-py37hfa6e2cd_0": {
    "md5": "0286cada8369b363e53969d840b58f27",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/pycosat-0.6.3-py37hfa6e2cd_0.tar.bz2"
  },
  "pycparser-2.19-py37_0": {
    "md5": "0cd552bb636e7933b51d0acef7a03304",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/pycparser-2.19-py37_0.tar.bz2"
  },
  "pymol-2.3.2-py37_1": {
    "md5": "d3383274b8cb9fa0daeee94233db0cee",
    "url": "https://conda.anaconda.org/schrodinger/win-64/pymol-2.3.2-py37_1.tar.bz2"
  },
  "pymol-upgrade-helper-1.0-1": {
    "md5": "9e394c07c1a86a9ca6b15dad72edcdfd",
    "url": "https://conda.anaconda.org/schrodinger/win-64/pymol-upgrade-helper-1.0-1.tar.bz2"
  },
  "pymol-web-examples-2.0-2": {
    "md5": "4267e27a1c109008878a6d5481bbce8b",
    "url": "https://conda.anaconda.org/schrodinger/noarch/pymol-web-examples-2.0-2.tar.bz2"
  },
  "pymol-web-legacy-2.0-py_1": {
    "md5": "b73e3cd531bd91866f1297788a8665f1",
    "url": "https://conda.anaconda.org/schrodinger/noarch/pymol-web-legacy-2.0-py_1.tar.bz2"
  },
  "pyopenssl-19.0.0-py37_0": {
    "md5": "c5d8afb5521012da459f88b62e0dfcc7",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/pyopenssl-19.0.0-py37_0.tar.bz2"
  },
  "pyqt-5.9.2-py37h6538335_2": {
    "md5": "e6802caccb172f1d7e5a9ba855050085",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/pyqt-5.9.2-py37h6538335_2.tar.bz2"
  },
  "pyreadline-2.1-py37_1": {
    "md5": "ad6e2a3e2e7eebb38cf66aeb27d73afd",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/pyreadline-2.1-py37_1.tar.bz2"
  },
  "pysocks-1.7.0-py37_0": {
    "md5": "e4e08777e59494ec127e19c9b2dd3d3d",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/pysocks-1.7.0-py37_0.tar.bz2"
  },
  "python-3.7.3-h8c8aaf0_1": {
    "md5": "da72a954aec12511e9a7d2474b853116",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/python-3.7.3-h8c8aaf0_1.tar.bz2"
  },
  "pywin32-223-py37hfa6e2cd_1": {
    "md5": "f109e28513833e67c4e0ef4a5cdb9cae",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/pywin32-223-py37hfa6e2cd_1.tar.bz2"
  },
  "qt-5.9.7-vc14h73c81de_0": {
    "md5": "2c80f3edd65223da2062f84292c9eac9",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/qt-5.9.7-vc14h73c81de_0.tar.bz2"
  },
  "requests-2.22.0-py37_0": {
    "md5": "0c6895204f6b9a733a0145b5dfd60890",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/requests-2.22.0-py37_0.tar.bz2"
  },
  "rigimol-1.3-2": {
    "md5": "35e50e07a352ff860c3d2f580e0b0178",
    "url": "https://conda.anaconda.org/schrodinger/win-64/rigimol-1.3-2.tar.bz2"
  },
  "ruamel_yaml-0.15.46-py37hfa6e2cd_0": {
    "md5": "163e6e1230812be0d46c0b8a0e1734ff",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/ruamel_yaml-0.15.46-py37hfa6e2cd_0.tar.bz2"
  },
  "setuptools-41.0.1-py37_0": {
    "md5": "af7f411babddc47623c43e39250e8328",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/setuptools-41.0.1-py37_0.tar.bz2"
  },
  "sip-4.19.8-py37h6538335_0": {
    "md5": "bc320ab7c6c8c8562913a06c5646aaa5",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/sip-4.19.8-py37h6538335_0.tar.bz2"
  },
  "six-1.12.0-py37_0": {
    "md5": "f5721c7381878ad9882026f356f5233a",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/six-1.12.0-py37_0.tar.bz2"
  },
  "sqlite-3.28.0-he774522_0": {
    "md5": "08bd6dc305316bf448e33004421313af",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/sqlite-3.28.0-he774522_0.tar.bz2"
  },
  "tk-8.6.8-hfa6e2cd_0": {
    "md5": "13cd026af1b90010901b7339bf19c35a",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/tk-8.6.8-hfa6e2cd_0.tar.bz2"
  },
  "urllib3-1.24.2-py37_0": {
    "md5": "e7b54e3583353a12c06a0a898f6357eb",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/urllib3-1.24.2-py37_0.tar.bz2"
  },
  "vc-14.1-h0510ff6_4": {
    "md5": "00bfe39f46f67409a376939cd2ab5039",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/vc-14.1-h0510ff6_4.tar.bz2"
  },
  "vs2013_runtime-12.0.21005-1": {
    "md5": "fc71f16dbab3acbc27ec911e506dcbd3",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/vs2013_runtime-12.0.21005-1.tar.bz2"
  },
  "vs2015_runtime-14.15.26706-h3a45250_4": {
    "md5": "850af9220d561dd550056957c981e62e",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/vs2015_runtime-14.15.26706-h3a45250_4.tar.bz2"
  },
  "wheel-0.33.4-py37_0": {
    "md5": "3af12b81202a646b264110460fd03fab",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/wheel-0.33.4-py37_0.tar.bz2"
  },
  "win_inet_pton-1.1.0-py37_0": {
    "md5": "4ec69295c06d6bb5d7735b98ecb393e6",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/win_inet_pton-1.1.0-py37_0.tar.bz2"
  },
  "wincertstore-0.2-py37_0": {
    "md5": "0418cf653523e6ef0174c259efaa78b4",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/wincertstore-0.2-py37_0.tar.bz2"
  },
  "xz-5.2.4-h2fa13f4_4": {
    "md5": "4d85f36274b4d66be67bb1fffd81613f",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/xz-5.2.4-h2fa13f4_4.tar.bz2"
  },
  "yaml-0.1.7-hc54c509_2": {
    "md5": "eaeee418b8ac33d392ea03bb9b6f389d",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/yaml-0.1.7-hc54c509_2.tar.bz2"
  },
  "zlib-1.2.11-h62dcd97_3": {
    "md5": "6f96fd91475cc78aabf76d2e1a9ed91f",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/zlib-1.2.11-h62dcd97_3.tar.bz2"
  },
  "zstd-1.3.7-h508b16e_0": {
    "md5": "251803bb127a0b26ee8c92106d3ab1ae",
    "url": "https://repo.anaconda.com/pkgs/main/win-64/zstd-1.3.7-h508b16e_0.tar.bz2"
  }
}
C_ENVS = {
  "root": [
    "python-3.7.3-h8c8aaf0_1",
    "apbs-1.5-1",
    "blas-1.0-mkl",
    "ca-certificates-2019.5.15-0",
    "glew-2.0.0-0",
    "icc_rt-2019.0.0-h0cc432a_1",
    "intel-openmp-2019.4-245",
    "mengine-1-0",
    "mpeg_encode-1-0",
    "mtz2ccp4_px-1.0-1",
    "pymol-upgrade-helper-1.0-1",
    "pymol-web-examples-2.0-2",
    "vs2013_runtime-12.0.21005-1",
    "vs2015_runtime-14.15.26706-h3a45250_4",
    "mkl-2019.4-245",
    "rigimol-1.3-2",
    "vc-14.1-h0510ff6_4",
    "icu-58.2-ha66f8fd_1",
    "jpeg-9b-hb83a4c4_2",
    "libiconv-1.15-h1df5818_7",
    "openssl-1.1.1c-he774522_1",
    "sqlite-3.28.0-he774522_0",
    "tk-8.6.8-hfa6e2cd_0",
    "xz-5.2.4-h2fa13f4_4",
    "yaml-0.1.7-hc54c509_2",
    "zlib-1.2.11-h62dcd97_3",
    "hdf5-1.10.4-h7ebc959_0",
    "libpng-1.6.37-h2a8f88b_0",
    "libxml2-2.9.9-h464c3ec_0",
    "zstd-1.3.7-h508b16e_0",
    "asn1crypto-0.24.0-py37_0",
    "certifi-2019.3.9-py37_0",
    "chardet-3.0.4-py37_1",
    "console_shortcut-0.1.1-3",
    "freetype-2.9.1-ha9979f8_1",
    "idna-2.8-py37_0",
    "libtiff-4.0.10-hb898794_2",
    "olefile-0.46-py37_0",
    "pdb2pqr-2.1.1-py37_1",
    "pmw-2.0.1-py37_2",
    "pycosat-0.6.3-py37hfa6e2cd_0",
    "pycparser-2.19-py37_0",
    "pyreadline-2.1-py37_1",
    "pywin32-223-py37hfa6e2cd_1",
    "qt-5.9.7-vc14h73c81de_0",
    "ruamel_yaml-0.15.46-py37hfa6e2cd_0",
    "sip-4.19.8-py37h6538335_0",
    "six-1.12.0-py37_0",
    "win_inet_pton-1.1.0-py37_0",
    "wincertstore-0.2-py37_0",
    "cffi-1.12.3-py37h7a1dbc1_0",
    "freemol-1.158-py37_1",
    "menuinst-1.4.16-py37he774522_0",
    "pillow-6.0.0-py37hdc69c19_0",
    "pyqt-5.9.2-py37h6538335_2",
    "pysocks-1.7.0-py37_0",
    "setuptools-41.0.1-py37_0",
    "cryptography-2.7-py37h7a1dbc1_0",
    "wheel-0.33.4-py37_0",
    "pip-19.1.1-py37_0",
    "pyopenssl-19.0.0-py37_0",
    "urllib3-1.24.2-py37_0",
    "requests-2.22.0-py37_0",
    "conda-4.6.14-py37_0",
    "pymol-web-legacy-2.0-py_1",
    "biopython-1.73-py37he774522_0",
    "h5py-2.9.0-py37h5e291fa_0",
    "mkl_fft-1.0.12-py37h14836fe_0",
    "mkl_random-1.0.2-py37h343c172_0",
    "numpy-1.16.4-py37h19fb1c0_0",
    "numpy-base-1.16.4-py37hc3f5095_0",
    "pymol-2.3.2-py37_1"
  ]
}



def _link(src, dst, linktype=LINK_HARD):
    if linktype == LINK_HARD:
        if on_win:
            from ctypes import windll, wintypes
            CreateHardLink = windll.kernel32.CreateHardLinkW
            CreateHardLink.restype = wintypes.BOOL
            CreateHardLink.argtypes = [wintypes.LPCWSTR, wintypes.LPCWSTR,
                                       wintypes.LPVOID]
            if not CreateHardLink(dst, src, None):
                raise OSError('win32 hard link failed')
        else:
            os.link(src, dst)
    elif linktype == LINK_COPY:
        # copy relative symlinks as symlinks
        if islink(src) and not os.readlink(src).startswith(os.path.sep):
            os.symlink(os.readlink(src), dst)
        else:
            shutil.copy2(src, dst)
    else:
        raise Exception("Did not expect linktype=%r" % linktype)


def rm_rf(path):
    """
    try to delete path, but never fail
    """
    try:
        if islink(path) or isfile(path):
            # Note that we have to check if the destination is a link because
            # exists('/path/to/dead-link') will return False, although
            # islink('/path/to/dead-link') is True.
            os.unlink(path)
        elif isdir(path):
            shutil.rmtree(path)
    except (OSError, IOError):
        pass


def yield_lines(path):
    for line in open(path):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        yield line


prefix_placeholder = ('/opt/anaconda1anaconda2'
                      # this is intentionally split into parts,
                      # such that running this program on itself
                      # will leave it unchanged
                      'anaconda3')

def read_has_prefix(path):
    """
    reads `has_prefix` file and return dict mapping filenames to
    tuples(placeholder, mode)
    """
    import shlex

    res = {}
    try:
        for line in yield_lines(path):
            try:
                parts = [x.strip('"\'') for x in shlex.split(line, posix=False)]
                # assumption: placeholder and mode will never have a space
                placeholder, mode, f = parts[0], parts[1], ' '.join(parts[2:])
                res[f] = (placeholder, mode)
            except (ValueError, IndexError):
                res[line] = (prefix_placeholder, 'text')
    except IOError:
        pass
    return res


def exp_backoff_fn(fn, *args):
    """
    for retrying file operations that fail on Windows due to virus scanners
    """
    if not on_win:
        return fn(*args)

    import time
    import errno
    max_tries = 6  # max total time = 6.4 sec
    for n in range(max_tries):
        try:
            result = fn(*args)
        except (OSError, IOError) as e:
            if e.errno in (errno.EPERM, errno.EACCES):
                if n == max_tries - 1:
                    raise Exception("max_tries=%d reached" % max_tries)
                time.sleep(0.1 * (2 ** n))
            else:
                raise e
        else:
            return result


class PaddingError(Exception):
    pass


def binary_replace(data, a, b):
    """
    Perform a binary replacement of `data`, where the placeholder `a` is
    replaced with `b` and the remaining string is padded with null characters.
    All input arguments are expected to be bytes objects.
    """
    def replace(match):
        occurances = match.group().count(a)
        padding = (len(a) - len(b)) * occurances
        if padding < 0:
            raise PaddingError(a, b, padding)
        return match.group().replace(a, b) + b'\0' * padding

    pat = re.compile(re.escape(a) + b'([^\0]*?)\0')
    res = pat.sub(replace, data)
    assert len(res) == len(data)
    return res


def update_prefix(path, new_prefix, placeholder, mode):
    if on_win:
        # force all prefix replacements to forward slashes to simplify need
        # to escape backslashes - replace with unix-style path separators
        new_prefix = new_prefix.replace('\\', '/')

    path = os.path.realpath(path)
    with open(path, 'rb') as fi:
        data = fi.read()
    if mode == 'text':
        new_data = data.replace(placeholder.encode('utf-8'),
                                new_prefix.encode('utf-8'))
    elif mode == 'binary':
        if on_win:
            # anaconda-verify will not allow binary placeholder on Windows.
            # However, since some packages might be created wrong (and a
            # binary placeholder would break the package, we just skip here.
            return
        new_data = binary_replace(data, placeholder.encode('utf-8'),
                                  new_prefix.encode('utf-8'))
    else:
        sys.exit("Invalid mode:" % mode)

    if new_data == data:
        return
    st = os.lstat(path)
    # unlink in case the file is memory mapped
    exp_backoff_fn(os.unlink, path)
    with open(path, 'wb') as fo:
        fo.write(new_data)
    os.chmod(path, stat.S_IMODE(st.st_mode))


def name_dist(dist):
    if hasattr(dist, 'name'):
        return dist.name
    else:
        return dist.rsplit('-', 2)[0]


def create_meta(prefix, dist, info_dir, extra_info):
    """
    Create the conda metadata, in a given prefix, for a given package.
    """
    # read info/index.json first
    with open(join(info_dir, 'index.json')) as fi:
        meta = json.load(fi)
    # add extra info
    meta.update(extra_info)
    # write into <prefix>/conda-meta/<dist>.json
    meta_dir = join(prefix, 'conda-meta')
    if not isdir(meta_dir):
        os.makedirs(meta_dir)
    with open(join(meta_dir, dist + '.json'), 'w') as fo:
        json.dump(meta, fo, indent=2, sort_keys=True)


def run_script(prefix, dist, action='post-link'):
    """
    call the post-link (or pre-unlink) script, and return True on success,
    False on failure
    """
    path = join(prefix, 'Scripts' if on_win else 'bin', '.%s-%s.%s' % (
            name_dist(dist),
            action,
            'bat' if on_win else 'sh'))
    if not isfile(path):
        return True
    if SKIP_SCRIPTS:
        print("WARNING: skipping %s script by user request" % action)
        return True

    if on_win:
        try:
            args = [os.environ['COMSPEC'], '/c', path]
        except KeyError:
            return False
    else:
        shell_path = '/bin/sh' if 'bsd' in sys.platform else '/bin/bash'
        args = [shell_path, path]

    env = os.environ
    env['PREFIX'] = prefix

    import subprocess
    try:
        subprocess.check_call(args, env=env)
    except subprocess.CalledProcessError:
        return False
    return True


url_pat = re.compile(r'''
(?P<baseurl>\S+/)                 # base URL
(?P<fn>[^\s#/]+)                  # filename
([#](?P<md5>[0-9a-f]{32}))?       # optional MD5
$                                 # EOL
''', re.VERBOSE)

def read_urls(dist):
    try:
        data = open(join(PKGS_DIR, 'urls')).read()
        for line in data.split()[::-1]:
            m = url_pat.match(line)
            if m is None:
                continue
            if m.group('fn') == '%s.tar.bz2' % dist:
                return {'url': m.group('baseurl') + m.group('fn'),
                        'md5': m.group('md5')}
    except IOError:
        pass
    return {}


def read_no_link(info_dir):
    res = set()
    for fn in 'no_link', 'no_softlink':
        try:
            res.update(set(yield_lines(join(info_dir, fn))))
        except IOError:
            pass
    return res


def linked(prefix):
    """
    Return the (set of canonical names) of linked packages in prefix.
    """
    meta_dir = join(prefix, 'conda-meta')
    if not isdir(meta_dir):
        return set()
    return set(fn[:-5] for fn in os.listdir(meta_dir) if fn.endswith('.json'))


def link(prefix, dist, linktype=LINK_HARD, info_dir=None):
    '''
    Link a package in a specified prefix.  We assume that the packacge has
    been extra_info in either
      - <PKGS_DIR>/dist
      - <ROOT_PREFIX>/ (when the linktype is None)
    '''
    if linktype:
        source_dir = join(PKGS_DIR, dist)
        info_dir = join(source_dir, 'info')
        no_link = read_no_link(info_dir)
    else:
        info_dir = info_dir or join(prefix, 'info')

    files = list(yield_lines(join(info_dir, 'files')))
    # TODO: Use paths.json, if available or fall back to this method
    has_prefix_files = read_has_prefix(join(info_dir, 'has_prefix'))

    if linktype:
        for f in files:
            src = join(source_dir, f)
            dst = join(prefix, f)
            dst_dir = dirname(dst)
            if not isdir(dst_dir):
                os.makedirs(dst_dir)
            if exists(dst):
                if FORCE:
                    rm_rf(dst)
                else:
                    raise Exception("dst exists: %r" % dst)
            lt = linktype
            if f in has_prefix_files or f in no_link or islink(src):
                lt = LINK_COPY
            try:
                _link(src, dst, lt)
            except OSError:
                pass

    for f in sorted(has_prefix_files):
        placeholder, mode = has_prefix_files[f]
        try:
            update_prefix(join(prefix, f), prefix, placeholder, mode)
        except PaddingError:
            sys.exit("ERROR: placeholder '%s' too short in: %s\n" %
                     (placeholder, dist))

    if not run_script(prefix, dist, 'post-link'):
        sys.exit("Error: post-link failed for: %s" % dist)

    meta = {
        'files': files,
        'link': ({'source': source_dir,
                  'type': link_name_map.get(linktype)}
                 if linktype else None),
    }
    try:    # add URL and MD5
        meta.update(IDISTS[dist])
    except KeyError:
        meta.update(read_urls(dist))
    meta['installed_by'] = 'PyMOL-2.3.2-Windows-x86_64.exe'
    create_meta(prefix, dist, info_dir, meta)


def duplicates_to_remove(linked_dists, keep_dists):
    """
    Returns the (sorted) list of distributions to be removed, such that
    only one distribution (for each name) remains.  `keep_dists` is an
    interable of distributions (which are not allowed to be removed).
    """
    from collections import defaultdict

    keep_dists = set(keep_dists)
    ldists = defaultdict(set) # map names to set of distributions
    for dist in linked_dists:
        name = name_dist(dist)
        ldists[name].add(dist)

    res = set()
    for dists in ldists.values():
        # `dists` is the group of packages with the same name
        if len(dists) == 1:
            # if there is only one package, nothing has to be removed
            continue
        if dists & keep_dists:
            # if the group has packages which are have to be kept, we just
            # take the set of packages which are in group but not in the
            # ones which have to be kept
            res.update(dists - keep_dists)
        else:
            # otherwise, we take lowest (n-1) (sorted) packages
            res.update(sorted(dists)[:-1])
    return sorted(res)


def yield_idists():
    for line in open(join(PKGS_DIR, 'urls')):
        m = url_pat.match(line)
        if m:
            fn = m.group('fn')
            yield fn[:-8]


def remove_duplicates():
    idists = list(yield_idists())
    keep_files = set()
    for dist in idists:
        with open(join(ROOT_PREFIX, 'conda-meta', dist + '.json')) as fi:
            meta = json.load(fi)
        keep_files.update(meta['files'])

    for dist in duplicates_to_remove(linked(ROOT_PREFIX), idists):
        print("unlinking: %s" % dist)
        meta_path = join(ROOT_PREFIX, 'conda-meta', dist + '.json')
        with open(meta_path) as fi:
            meta = json.load(fi)
        for f in meta['files']:
            if f not in keep_files:
                rm_rf(join(ROOT_PREFIX, f))
        rm_rf(meta_path)


def determine_link_type_capability():
    src = join(PKGS_DIR, 'urls')
    dst = join(ROOT_PREFIX, '.hard-link')
    assert isfile(src), src
    assert not isfile(dst), dst
    try:
        _link(src, dst, LINK_HARD)
        linktype = LINK_HARD
    except OSError:
        linktype = LINK_COPY
    finally:
        rm_rf(dst)
    return linktype


def link_dist(dist, linktype=None):
    if not linktype:
        linktype = determine_link_type_capability()
    prefix = prefix_env('root')
    link(prefix, dist, linktype)


def link_idists():
    linktype = determine_link_type_capability()
    for env_name in sorted(C_ENVS):
        dists = C_ENVS[env_name]
        assert isinstance(dists, list)
        if len(dists) == 0:
            continue

        prefix = prefix_env(env_name)
        for dist in dists:
            assert dist in IDISTS
            link_dist(dist, linktype)

        for dist in duplicates_to_remove(linked(prefix), dists):
            meta_path = join(prefix, 'conda-meta', dist + '.json')
            print("WARNING: unlinking: %s" % meta_path)
            try:
                os.rename(meta_path, meta_path + '.bak')
            except OSError:
                rm_rf(meta_path)


def prefix_env(env_name):
    if env_name == 'root':
        return ROOT_PREFIX
    else:
        return join(ROOT_PREFIX, 'envs', env_name)


def post_extract(env_name='root'):
    """
    assuming that the package is extracted in the environment `env_name`,
    this function does everything link() does except the actual linking,
    i.e. update prefix files, run 'post-link', creates the conda metadata,
    and removed the info/ directory afterwards.
    """
    prefix = prefix_env(env_name)
    info_dir = join(prefix, 'info')
    with open(join(info_dir, 'index.json')) as fi:
        meta = json.load(fi)
    dist = '%(name)s-%(version)s-%(build)s' % meta
    if FORCE:
        run_script(prefix, dist, 'pre-unlink')
    link(prefix, dist, linktype=None)
    shutil.rmtree(info_dir)


def multi_post_extract():
    # This function is called when using the --multi option, when building
    # .pkg packages on OSX.  I tried avoiding this extra option by running
    # the post extract step on each individual package (like it is done for
    # the .sh and .exe installers), by adding a postinstall script to each
    # conda .pkg file, but this did not work as expected.  Running all the
    # post extracts at end is also faster and could be considered for the
    # other installer types as well.
    for dist in yield_idists():
        info_dir = join(ROOT_PREFIX, 'info', dist)
        with open(join(info_dir, 'index.json')) as fi:
            meta = json.load(fi)
        dist = '%(name)s-%(version)s-%(build)s' % meta
        link(ROOT_PREFIX, dist, linktype=None, info_dir=info_dir)


def main():
    global SKIP_SCRIPTS, ROOT_PREFIX, PKGS_DIR

    p = OptionParser(description="conda post extract tool used by installers")

    p.add_option('--skip-scripts',
                 action="store_true",
                 help="skip running pre/post-link scripts")

    p.add_option('--rm-dup',
                 action="store_true",
                 help="remove duplicates")

    p.add_option('--multi',
                 action="store_true",
                 help="multi post extract usecase")

    p.add_option('--link-dist',
                 action="store",
                 default=None,
                 help="link dist")

    p.add_option('--root-prefix',
                 action="store",
                 default=abspath(join(__file__, '..', '..')),
                 help="root prefix (defaults to %default)")

    p.add_option('--post',
                 action="store",
                 help="perform post extract (on a single package), "
                      "in environment NAME",
                 metavar='NAME')

    opts, args = p.parse_args()
    ROOT_PREFIX = opts.root_prefix.replace('//', '/')
    PKGS_DIR = join(ROOT_PREFIX, 'pkgs')

    if args:
        p.error('no arguments expected')

    if FORCE:
        print("using -f (force) option")

    if opts.post:
        post_extract(opts.post)
        return

    if opts.skip_scripts:
        SKIP_SCRIPTS = True

    if opts.rm_dup:
        remove_duplicates()
        return

    if opts.multi:
        multi_post_extract()
        return

    if opts.link_dist:
        link_dist(opts.link_dist)
        return

    if IDISTS:
        link_idists()
    else:
        post_extract()


def warn_on_special_chrs():
    if on_win:
        return
    for c in SPECIAL_ASCII:
        if c in ROOT_PREFIX:
            print("WARNING: found '%s' in install prefix." % c)


if __name__ == '__main__':
    main()
    warn_on_special_chrs()
