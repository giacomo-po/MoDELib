testmg.py : generate a configuration and populate evl/ F/ folders
testdc.py : read a configuration and run a few glide steps

To run testmg.py or testdc.py, this directory must contain:
 - a link or copy of ddpy.<python-arch>.so, generated from tools/DDpy
 - modelib input and output directories:
   * evl, F, inputFiles, where inputFiles may contain files with paths
     that need to be adjusted:
      - MaterialsLibrary, MeshLibrary, MicrostructureLibrary, NoiseLibrary
   * the contents of evl/ and F/ may be initiated using either modelib's
       microstructureGenerator program or testmg.py

