# Manifold-Model-Mesh-Shape
An implementation of the Manifold Model for Shape Deformation outlined in a paper by Freifeld and Black

# Copyright Information
The Geomproc library was authored by Dr. Oliver Van Kiack and provided under an MIT license  

The datasets used are from the [TOSCA](http://tosca.cs.technion.ac.il/data/toscahires.txt) and [DFAUST](https://dfaust.is.tuebingen.mpg.de/license.html) for academic purposes (course project) and are copyright their original publishers. Clicking either dataset will bring up the relevant license information for said dataset


## Prerequisites 
The model requires python3 and numpy to run, disturbing the file structure for the datasets will almost certainly cause issues as well.

## Running information
Assuming the prerequisites are met, the model can be run from the command line using 
```python ManifoldMeshModel.py```   
for a Windows system, substitute python3 for unix systems

By default, the model will load the 11 cat meshes from the TOSCA dataset and provide some analytic information from them. 